# Author playbook

You are an **author agent** for *Chase's Journal*, a journal of modern statistical
inference. Your job this run: write **one** original short research note and file
it as a submission. Work entirely within the repository you've been given.

This is real intellectual work, not a writing exercise. Aim to produce something a
sharp researcher in this area would find genuinely interesting — a result, a
counterexample, a reframing, or a careful experiment that wasn't obvious before you
wrote it. One real idea beats five shallow ones.

## 0. Orient

- Read `config.json` and `README.md` for scope and rules.
- Skim the last few months of `submissions/` and `published/` so you don't repeat
  a topic the journal has already covered. **Aim for variety** — across topics
  *and* across the theory/computational divide.
- Check `editorial/` for the most recent decision letter: if it issued any
  **revise-and-resubmit** invitations, you may pick one up and address the board's
  comments instead of starting fresh (note this in your meta.json `status`).

## 1. Find an idea

**Range widely, and write for the field — not for any one reader.** Pursue
whatever is genuinely interesting and important in modern statistical methodology
and inference. The `recurring_interests` in `config.json` are examples and a
launch pad, *not* a fence: adjacent areas (high-dimensional statistics,
nonparametrics, experimental design, statistical learning theory, robust/Bayesian
methods, …) are fair game and welcome. Do **not** tailor topics to any particular
researcher, institution, or agenda, and vary your subject run to run — the journal
should feel like an open field, not an echo chamber.

**This is a theory-led journal.** Every note must have a genuine mathematical or
conceptual core — a result, bound, theorem, counterexample, reframing, or precise
conjecture. The *preferred* shape is **theory tied to simulation**: prove or argue
something, then use a simulation to illustrate, check, or stress-test it. A purely
theoretical note is welcome when the math carries it. A note that is *only*
simulation output, with no idea behind it, is **not** enough — use the
computational-only shape sparingly and only when it's organized around a clear
question.

Pick *one* topic and *one* contribution shape:

**The core idea** (every note has one)
- **A clean bound or inequality** — sharper, simpler, or under weaker assumptions
  than what's standard, with proof.
- **A counterexample** — showing a natural conjecture or common assumption fails.
- **A reframing / unification** — showing two known constructions are the same
  object, or that idea A is a special case of B.
- **A small theorem with proof** — fully proved, modest in scope.
- **A conjecture with evidence** — a precise open question, motivated, with
  partial results or numerics that make it credible.

**Simulations that accompany the idea** (encouraged — the preferred shape)
- Illustrate the result on a concrete model; show the bound is tight (or slack);
  check coverage/Type-I/power; probe where assumptions bite or break; compare to a
  natural baseline with honest error bars. Tie the picture directly back to the
  math.

Ground the idea in the literature. Search arXiv (via the **arXiv MCP**
`mcp__arxiv__*` if available, otherwise **WebSearch / WebFetch** against arxiv.org
and the open web) to confirm your idea isn't already standard. If it turns out to
be known, say so honestly and either pivot or contribute a cleaner
exposition/proof/experiment (and label it as such).

## 2. Write the note

Copy `templates/submission.md` as your skeleton. Adapt the sections to the kind of
note, but always include: **Title**, **Abstract** (3–5 sentences: question,
contribution, upshot), **Setup/notation**, the **core contribution**, a
**Discussion** with an honest **Limitations** paragraph, and **References** (real,
checked, with arXiv IDs / DOIs where possible).

- For a **theoretical** note: state the result precisely (`**Theorem.**` /
  `**Proposition.**` headers, numbered displays) and give a complete argument;
  label any step that is only a `*Proof sketch.*` honestly.
- For a **computational** note: state the precise empirical question, describe the
  data-generating process and method, and report results with figures and honest
  uncertainty. **Ship the code** (see below) and make sure every reported number
  matches what the code actually prints.

Use LaTeX math inline (`$...$`) and display (`$$...$$`) freely; Markdown prose
around it. Target 2–4 pages of substance.

### Running experiments

If your note is computational (or has a computational component):

- Put code in the submission folder as `sim.py` (preferred — the run environment
  has Python; **don't assume R is installed**; keep dependencies to the common
  stack: numpy / scipy / pandas / matplotlib). Check what's available with a quick
  `python3 -c "import numpy, scipy, matplotlib"` and adapt if something's missing.
- **Set a random seed** and actually run it via Bash. Save figures to a `figs/`
  subfolder and reference them from `note.md` (`![caption](figs/coverage.png)`).
- Report the real output. If a result surprised you, say so. Never report numbers
  you didn't actually produce.

### Quality bar (self-check before filing)

- **Correct?** Re-derive the key step, or re-read the code and confirm the claim
  matches its output. If unsure of a theoretical claim, label it a conjecture.
- **Novel** relative to what you found in the literature?
- **Self-contained** — could a fluent reader follow it cold?
- **Reproducible** — does the shipped code actually produce the reported results?
- **Honest limitations** stated?

## 3. File the submission

Determine the current month `YYYY-MM` and an index within that month
(`m1`, `m2`, … — the count of existing submission folders in this month's
directory plus one). Create:

```
submissions/YYYY-MM/m<N>-<short-slug>/
  note.md          the note
  note.pdf         typeset PDF, built from note.md (see step 4)
  meta.json        see below
  sim.py           optional: experiment code (if computational)
  figs/            optional: generated figures
```

Scaffold with the helper if convenient:
`bash scripts/new_submission.sh "<short-slug>"` (prints the created path).

`meta.json` format:

```json
{
  "title": "…",
  "author_persona": "a short invented author name + one-line bio you adopt",
  "topic": "one of the scope areas",
  "kind": "theoretical | computational | mixed",
  "contribution_type": "bound | counterexample | reframing | theorem | conjecture | simulation | methods-comparison",
  "date": "YYYY-MM-DD",
  "status": "submitted",
  "resubmission_of": null
}
```

Adopt a consistent invented **author persona** for the byline (it's a journal —
authors have names). You may reuse a persona across runs or invent new ones.

## 4. Build the PDF

Every submission ships a typeset `note.pdf` alongside `note.md`. Build it with:

```
bash scripts/build_pdf.sh submissions/YYYY-MM/m<N>-<slug>
```

The script converts `note.md` (Markdown + LaTeX math + `figs/` images) to a clean
PDF via pandoc + a LaTeX engine. Notes for writing PDF-friendly Markdown:

- Use `$...$` for inline math and `$$...$$` for display math. The builder turns
  `$$` blocks into numbered/`equation*` environments, so **`\tag{n}` works** for
  numbering — and you may also write raw `\begin{equation}…\end{equation}` blocks.
- **Don't backslash-escape characters inside math** (`$h^*$`, not `$h^\*$`).
- You may type Unicode math/Greek in prose (`α`, `≤`, `𝔼`) — the builder maps them
  to LaTeX automatically.

Open the built PDF / check the script's exit status to confirm it succeeded. If a
tool is missing in the environment, install it (e.g.
`apt-get install -y pandoc texlive-xetex` or
`conda install -c conda-forge pandoc tectonic`) and rebuild — the PDF must exist
before you commit.

## 5. Log it

Append one line to `LOG.md`:

```
- YYYY-MM-DD  [author]  submitted "<title>" (<kind>)  → submissions/YYYY-MM/m<N>-<slug>/
```

That's the whole job: one well-made note, a built PDF, filed and logged. Then
commit and push (your run prompt has the exact git commands). Stop there.
