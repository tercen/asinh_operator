## 0.1.8 (2026-09-28) — typed sidecars for the auto test's tables

0.1.7 got past the row order and failed on `bad.column.type -- cells_used -- expected double found
int32`: with no `.schema` sidecar the platform types a fixture's columns from the CSV, and the
cofactor table's `cells_used` / `seed` are int32. `tests/auto_table1..4.csv.schema` now declare
every column's type, in the CSV's column order (the platform matches sidecar to CSV by position).

## 0.1.7 (2026-09-28) — the auto test's fixtures, in the order the platform stores rows

Installing 0.1.6 on Studio failed the `asinh_auto_shape` test with
`bad.value -- At column .ri -- at line 2 -- val = 1 refVal = 0`. The values were right; the
**row order** was not. The platform compares each stored table line by line, and a `(.ri, .ci)`
result is stored in the order the operator wrote it, which is the order the crosstab streamed
in: 2 x 2 blocks of channel by cell (`(0,0) (0,1) (1,0) (1,1) (0,2) (0,3) ...`), not row-major.
The fixed test never saw this because its expected table came from a real run; the auto test's
had been written by hand, row-major.

`tests/auto_table1..4.csv` are now the tables of a Studio run of this operator on the test's own
projection (`crabs-long.csv`, `observation` on columns, `variable` on rows, `method = auto`),
exported with `tercenctl data export-csv`. Every value agrees with the previous fixtures to
9e-16; only the order changed. Expected tables for a `.ri/.ci` result must always be taken from
a run on the platform — never written or sorted by hand.

## 0.1.6 (2026-09-22) — memory model reshaped from three production runs

Measured on tercen.com (`stats_d_actual_ram_peak` / `_anon`), fixed cofactors:

| values | CPUs | wall | booked (0.1.5) | peak | anon |
|---|---|---|---|---|---|
| 9.0 M | 1 | 0.9 min | 425 MB | 319 MB | 178 MB |
| 19.1 M | 2 | 1.7 min | 689 MB | **631 MB (92 %)** | 346 MB |
| 45.1 M | 4 | 4.2 min | 1371 MB | 766 MB | 131 MB |

Anonymous memory does not grow with the crosstab — the operator streams it in 1 M-value chunks — and
the page cache of the result file saturates because the writer releases every 256 MB. The 0.1.5
model (25 B/value + 180 MB) had the wrong shape: 92 % used at 19 M values, and 8 GB demanded for the
312 M-value cohort. Now **2 B/value + 900 MB**: 19 M → 938 MB (67 %), 45 M → 990 MB, 312 M → 1.5 GB.

Throughput is ~180,000 values/s regardless of CPUs (0.9 / 1.7 / 4.2 min): the time is the crosstab
paging out of the server and the result going back, not the arcsinh. Extra booked cores are unused
until the reader fetches page ranges concurrently — a design change, not done here.

## 0.1.5: `auto` no longer takes a sample per cell (2026-09-21)

Faris ran `auto` on tercen.com and the task was killed for memory. Not a booking problem: with no
`sample_factor` set, the operator took the **first column factor** as the sample — and on a
cytometry projection that is the event id, so every cell became its own sample. Then the
estimation floor of 500 cells per sample beat the 128 MB budget, and each reservoir
preallocated 4 KB. On a 100,000-cell crosstab in Studio that was 15,000 "samples" (one server
page of the column table — the column reader was also single-page) and **1.33 GB** against a
600 MB booking.

Three fixes, each measured: with no `sample_factor`, a column factor that has one value per
column is treated as *no* sample factor and the crosstab is one sample, with a warning that
says what to set; the budget now refuses rather than allocating when the floor cannot be met
(`… cannot be estimated … Set 'sample_factor' …`); reservoirs grow on demand. The column reader
walks every page, as the cell reader already did.

| projection | before | after |
|---|---|---|
| 100,000 cells x 20 channels, no sample factor | 1.33 GB | **60 MB** |
| 500,000 cells x 20 channels, 12 samples | — | **211 MB** |

That is 18.8 bytes per value in collect mode (the three columns, 16 B, plus the result buffer),
and the estimation is bounded whatever the sample count. The memory model was a constant
600 MB — designed, never measured, and noted as such since the first release. It is now
`0.000025 x n_main + 180 MB`: 430 MB booked where 211 was used at 10 M values. Above
`collect_max_cells` (20 M) the operator streams and the booking over-estimates; lower the
property if a cohort-scale step is refused for memory.

## 0.1.4: the Cofactors table could not be built on (2026-09-21)

`auto` emitted its Cofactors table as a second relation with an explicit `JoinOperator` on an
empty `ColumnPair` — a shape copied from `read_fcs`, where it is right because an import operator
has no crosstab on the left. Here it produced a composite the query engine refuses:
`bad relation -- !relation.hasAnyAttributes(attrs)` on any downstream step. Every `auto` result
was a dead end, and no `cargo test` could see it, because the join happens inside Tercen after
the operator's bytes.

The server's rule, found by trying the alternatives on Studio: every table in `tables` must
carry `.ri` or `.ci` (`OperatorResult -- .ci or .ri attribute is required` otherwise), and the
server joins it by those. The Cofactors table now carries `.ri`, one row per channel, and no
join is declared. Downstream projection: `DoneState`; the tree shows the server keying it on
`.ri` and attaching it through the row factor.

`tests/asinh_auto.json` is the platform's unit test for that shape. Its numbers are
self-recorded on `crabs-long.csv` — a dataset with one sample and 200 cells per channel, on
which flowVS finds no populations and every channel comes back `unstable`. It guards the join,
not the estimate; the estimate is guarded by flowvs-rs's parity tests against R. Known and
left alone: an `unstable` channel still uses flowVS's number rather than `scale`, because the
number is finite; whether it should fall back is a separate decision.

# asinh_operator (Rust, 2.x; developed as asinh_rust_operator) — status, morning of 2026-09-20

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

The `auto` run also raised a question about flowVS rather than the operator, and measuring it
answered it: on **real** data the cofactor moves by 5–20 % under resampling, which is the order
the plan already tolerates, so a frozen table is sound as long as the subsample size and seed
travel with it. On synthetic data with a very tight negative population the estimate collapses
to a degenerate minimum instead (2 versus 99), which is the dim-channel failure the cofactor
check saw on the study panel. See CLAUDE.md for the numbers and the fix.

## What is not done

- ~~Memory unmeasured at scale~~ — done. The streaming path ran on a real 19.53 M-cell crosstab:
  17.5 s end to end, peak RSS 35.5 MB, spill files 8 B/cell. The 600 MB booking is generous; the
  collect path at its 20 M cap is what it is really sized for.
- **Superseded:** The runs above peaked at 11.6 MB and 28.5 MB, which
  says nothing about the 20 M-cell cap the 600 MB booking is designed around, and the streaming
  path has never run against Tercen at all.
- ~~`tests/test.json` points at the R operator's goldens~~ — the auto test's tables now come from a
  Studio run of this operator (0.1.7); the fixed test's came from one already.
- ~~`auto` reports no fragility signal~~ — done. The cofactor table carries a status
  (resolved / fragile / floored / unstable), what flowVS said, the σ_neg cofactor and the
  runner-up. The floor is opt-in via `cofactor_floor`; 2.5 is the cofactor check's convention.
- **Thresholds are uncalibrated.** `fragile` fires when a runner-up is within 10 % on objective
  and more than 1.5× away in cofactor. Those are starting values; the 32-channel cofactor check
  is the set to calibrate them on.
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
