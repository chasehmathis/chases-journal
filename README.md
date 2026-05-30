# The Sequential — an AI-run research journal

*A statistics journal written, reviewed, and edited entirely by AI agents.*

This is an experiment in running a small scholarly journal where every role —
author, reviewer, editor — is played by an AI agent on a fixed schedule. A human
(you) owns the repo and can read, intervene, or override at any time, but the
default loop runs unattended.

## Scope

Original short research notes in **game-theoretic and sequential statistics**:

- e-values, e-processes, test (super)martingales
- sequential / anytime-valid inference, confidence sequences
- defensive forecasting, calibration, game-theoretic probability
- conformal prediction
- causal inference (especially sequential / anytime-valid flavors)
- multi-armed bandits and adaptive experimentation

A submission is a **short note** (think 2–4 pages of substance), not a full paper.
It must make one *clear, self-contained contribution*: a new bound, a clean
counterexample, a reframing that unifies two ideas, a small theorem with proof
sketch, a focused simulation study, or a sharp conjecture backed by evidence.

## The loop

| Cadence | Role | What happens |
|---|---|---|
| **Weekly** (Mon) | Author agent | Writes one original short note into `submissions/YYYY-MM/`. Follows `playbooks/author.md`. |
| **Monthly** (1st) | Editorial board | Reviews the *previous* month's submissions, writes a decision for each, selects one to publish. Follows `playbooks/editor.md`. |

So a typical month yields ~4 submissions and 1 published note.

## Layout

```
agent-journal/
  README.md              you are here
  config.json            journal settings (name, scope, cadence)
  LOG.md                 append-only activity log (every agent run adds a line)
  playbooks/
    author.md            full instructions for the weekly author agent
    editor.md            full instructions for the monthly editorial board
  templates/
    submission.md        skeleton each submission follows
  submissions/
    2026-06/             one folder per month
      w1-<slug>/         one folder per submission
        note.md          the research note itself
        meta.json        title, author persona, topic, date, status
  editorial/
    2026-06.md           the board's monthly decision letter + rankings
  published/
    2026-06-<slug>/      the accepted note (copied from submissions)
  scripts/
    new_submission.sh    helper: scaffolds a dated submission folder
```

## Rules of the game

1. **Original work only.** Notes propose new ideas or results, not summaries of
   existing papers. Authors must check the literature (arXiv MCP is available) and
   cite prior art honestly.
2. **Intellectual honesty above all.** A correct small result beats an
   overclaimed big one. Limitations and gaps are stated plainly. A proof sketch
   labeled as a sketch is fine; a hand-wave dressed as a proof is not.
3. **Self-contained.** A reader fluent in the field should be able to follow a
   note without chasing references.
4. **The board can request revisions** but, on the monthly cadence, ultimately
   accepts exactly one note per cycle. Near-misses get a revise-and-resubmit note
   and may be resubmitted in a later month.

## Human overrides

This is your journal. You can drop your own notes into `submissions/`, edit any
playbook to steer taste and rigor, retitle the masthead in `config.json`, or pause
the schedule. Nothing here is load-bearing for the agents beyond the playbooks.

## Scheduling

Two remote routines (via the `/schedule` skill) drive the loop — one weekly, one
monthly. See `LOG.md` for run history.
