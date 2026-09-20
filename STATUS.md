# asinh_rust_operator — status, morning of 2026-09-20

Built overnight against the goal in `~/tercen/goals/2026-09-19-asinh-flowvs.md`. **Local git only:
no remote, no image published, not installed anywhere.** One commit.

## Where it got to

The first crosstab operator of the port, following the `create-rust-operator` skill.

| check | result |
|---|---|
| parity with the R operator's own golden | 1,000 cells, worst relative difference **2.2e-16** |
| `cargo test` | 13 tests green (unit, golden, spec drift) |
| `cargo clippy -D warnings`, `cargo fmt --check` | clean |
| image | 6.4 MB compressed, 24.9 MB on disk, under the 20 MB static-tier gate |
| runs as `--user 1000:1000` | yes, exits 1 with "TERCEN_TASK_ID is not set" |
| `operator.json` spec | mirrored from the R operator, with a test that it matches what the writer emits |

Parity uses the R operator's published fixtures and projection, not goldens invented here, and
needs no Tercen instance: `.ri` indexes the row facet table in file order and `.ci` the column one.

## Reviewed and restructured, 2026-09-20

A structural review found seven things and all are fixed: five undeclared properties, a spec that
described one output while `auto` emits two, estimation memory governed by a property rather than
a budget, three passes over the crosstab where one does, `auto` reading the input twice in collect
mode, a hard-coded sample factor, and shared modules that had already diverged from read_fcs.
`CLAUDE.md` has the reasoning. Both cargo configurations build, test and lint clean.

## First runs against Tercen, 2026-09-20

Studio 1.1.8, via `dev/setup_crabs.py`.

| run | result |
|---|---|
| `fixed` on the R operator's crabs fixture, 1,000 cells | server ingested and linked it; the exported table matches the R golden to **2.2e-16** |
| `auto` on a synthetic 8-channel set, 160,000 cells | ran in 1.0 s, peak RSS 28.5 MB; the Cofactors table landed with all six columns |

Two bugs came out of the first run, both about talking to the server rather than arithmetic: an
unbounded chunk loop (fixed by counting rows against the schema) and the wrong TSON layout for
streamed tables (both layouts are accepted now). Neither was reachable from a unit test.

The `auto` run also exposed something about flowVS rather than the operator: the estimate can
swing between 2 and 99 for the same channel depending on the subsample size, with a healthy
Bartlett statistic at both. See CLAUDE.md. **Do not freeze a cofactor table from this operator
until the stability check exists.**

## What is not done

- **Memory is still unmeasured at scale.** The runs above peaked at 11.6 MB and 28.5 MB, which
  says nothing about the 20 M-cell cap the 600 MB booking is designed around, and the streaming
  path has never run against Tercen at all.
- **`tests/test.json` points at the R operator's goldens.** Correct for parity, but the platform
  unit test should be regenerated from a Studio run of *this* operator once it installs.
- **`auto` has no *stability* check**: see above. The path itself now works end to end.
- **`manual` mode has no end-to-end test**, only a unit test, because the R operator ships no
  fixture for it. It is the mode CYTOSHRINK will actually use, with a cofactor table as the
  second row factor.
- **No `auto`.** Deliberate: that is flowVS, now implemented in `~/tercen/flowvs-rs`, which could
  be wired in here later — though `flowvs-rust-plan.md` §10 Q2 recommends against it for
  Tercen-old, since the R operator already covers that case.

## Suggested next steps

1. A dev run on Studio against a real crosstab, then refit `memory_model.json` from the measured peak.
2. A container run under `--memory 600M` on the largest crosstab you care about.
3. `manual` mode with a cofactor table, which is the pipeline's real use.
4. Only then a repository, a tag and an install — all of which are yours to decide.
