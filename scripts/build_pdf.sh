#!/usr/bin/env bash
# Build note.pdf from note.md in a submission (or published) folder.
# Usage: bash scripts/build_pdf.sh <dir>      (defaults to current directory)
#
# Requires pandoc + a LaTeX engine (xelatex / pdflatex / lualatex / tectonic).
# Handles GitHub-flavored Markdown with LaTeX math ($...$ and $$...$$, including
# \tag) and relative images (figs/...).
set -euo pipefail

DIR="${1:-.}"
cd "$DIR"
[ -f note.md ] || { echo "build_pdf: no note.md in $(pwd)" >&2; exit 1; }

# 1. Preprocess: wrap $$...$$ display blocks in amsmath environments so that
#    \tag (and other amsmath-only constructs) compile. Bare \[...\] forbids \tag.
python3 - <<'PY'
import re
src = open('note.md', encoding='utf-8').read()

# Unicode math/Greek -> LaTeX command. Authors freely type these in prose
# ("matched α = 0.05") where the default text font has no glyph; route them
# through math so the PDF is font-independent.
UMAP = {
 'α':r'\alpha','β':r'\beta','γ':r'\gamma','δ':r'\delta','ε':r'\varepsilon',
 'ζ':r'\zeta','η':r'\eta','θ':r'\theta','ϑ':r'\vartheta','ι':r'\iota',
 'κ':r'\kappa','λ':r'\lambda','μ':r'\mu','ν':r'\nu','ξ':r'\xi','π':r'\pi',
 'ρ':r'\rho','σ':r'\sigma','ς':r'\varsigma','τ':r'\tau','υ':r'\upsilon',
 'φ':r'\varphi','ϕ':r'\phi','χ':r'\chi','ψ':r'\psi','ω':r'\omega',
 'Γ':r'\Gamma','Δ':r'\Delta','Θ':r'\Theta','Λ':r'\Lambda','Ξ':r'\Xi',
 'Π':r'\Pi','Σ':r'\Sigma','Φ':r'\Phi','Ψ':r'\Psi','Ω':r'\Omega',
 '×':r'\times','÷':r'\div','±':r'\pm','∓':r'\mp','·':r'\cdot','⋅':r'\cdot',
 '≤':r'\le','≥':r'\ge','≠':r'\ne','≈':r'\approx','∼':r'\sim','≃':r'\simeq',
 '≡':r'\equiv','∞':r'\infty','→':r'\to','↦':r'\mapsto','←':r'\leftarrow',
 '∈':r'\in','∉':r'\notin','∋':r'\ni','⊆':r'\subseteq','⊂':r'\subset',
 '⊇':r'\supseteq','∪':r'\cup','∩':r'\cap','∅':r'\emptyset','∖':r'\setminus',
 '∑':r'\sum','∏':r'\prod','∫':r'\int','√':r'\surd','∂':r'\partial',
 '∇':r'\nabla','∀':r'\forall','∃':r'\exists','¬':r'\neg','∝':r'\propto',
 '≪':r'\ll','≫':r'\gg','⟨':r'\langle','⟩':r'\rangle','…':r'\dots',
 '⊗':r'\otimes','⊕':r'\oplus','⇒':r'\Rightarrow','⇔':r'\Leftrightarrow',
 'ℝ':r'\mathbb{R}','ℕ':r'\mathbb{N}','ℤ':r'\mathbb{Z}','ℚ':r'\mathbb{Q}',
 'ℙ':r'\mathbb{P}','𝔼':r'\mathbb{E}','′':"'",'″':"''",
}

def to_math(s):                          # inside math: bare commands
    s = s.replace('\\*', '*')            # undo over-escaped Markdown specials
    for u, l in UMAP.items():
        s = s.replace(u, l + (' ' if l[-1].isalpha() else ''))
    return s

def to_text(s):                          # in prose: wrap replacements in $...$
    for u, l in UMAP.items():
        s = s.replace(u, '$' + l + '$')
    return s

parts = src.split('$$')
if len(parts) % 2 == 0:
    raise SystemExit('build_pdf: unbalanced $$ in note.md')
out = []
for i, p in enumerate(parts):
    if i % 2 == 1:                       # display-math block
        body = to_math(p)
        env = 'equation' if '\\tag' in body else 'equation*'
        out.append('\n\\begin{%s}%s\\end{%s}\n' % (env, body, env))
    else:                                # prose, with possible inline $...$
        for seg in re.split(r'((?<!\\)\$[^$\n]+\$)', p):
            if len(seg) >= 2 and seg[0] == '$' and seg[-1] == '$':
                out.append('$' + to_math(seg[1:-1]) + '$')
            else:
                out.append(to_text(seg))
open('.note.build.md', 'w', encoding='utf-8').write(''.join(out))
PY

# 2. Pick a LaTeX engine.
ENGINE=""
for e in xelatex pdflatex lualatex tectonic; do
  if command -v "$e" >/dev/null 2>&1; then ENGINE="$e"; break; fi
done

# 3. Convert.
if command -v pandoc >/dev/null 2>&1 && [ -n "$ENGINE" ]; then
  pandoc .note.build.md -o note.pdf \
    --pdf-engine="$ENGINE" \
    -V geometry:margin=1in -V fontsize=11pt -V colorlinks=true -V linkcolor=blue \
    --resource-path=.
  rm -f .note.build.md
  echo "build_pdf: wrote $(pwd)/note.pdf via pandoc + $ENGINE"
else
  rm -f .note.build.md
  echo "build_pdf: missing tools (pandoc=$(command -v pandoc || echo NO), latex-engine=${ENGINE:-NONE})." >&2
  echo "build_pdf: install pandoc and a LaTeX engine (e.g. 'apt-get install -y pandoc texlive-xetex' or 'conda install -c conda-forge pandoc tectonic')." >&2
  exit 2
fi
