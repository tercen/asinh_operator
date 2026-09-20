# asinh_rust_operator — maintenance notes

Rust port of `tercen/asinh_operator` (`fixed` and `manual`). Built 2026-09-19/20 following the
`create-rust-operator` skill, as the first **crosstab** operator of the CYTOSHRINK port —
`read_fcs_rust_operator` proved document input, not this.

## What it is

Per-cell transform: `asinh(y / cofactor)`. `fixed` takes the cofactor from the `scale` property,
`manual` from the second row factor (R: `ctx$rnames[[2]]`). Output is one column,
`<namespace>.asinh`, on `.ri` / `.ci`.

`auto` (flowVS estimation) is deliberately absent. `props::settings_from_ctx` rejects it with a
message pointing at the R operator; the estimator is being ported as `flowvs-rs`, per
`flowvs-rust-plan.md` §10 Q2 ("a Rust operator on the old runtime would only buy speed on a
once-per-study step").

## Parity

`tests/r_goldens.rs` reproduces the **R operator's own** golden — its `crabs-long.csv` input, its
`table1.csv` output, its projection from `test_workflow.json` — to a worst relative difference of
**2.2e-16** over 1,000 cells. `.ri` indexes `table3.csv` (the row facet table) in file order and
`.ci` indexes `table2.csv`, which is what makes the golden reproducible with no Tercen instance.

`manual` has no R fixture (the R operator ships none), so it is checked against the definition.

## Structure, after the 2026-09-20 review

The operator is one crate with three methods and two size paths, and the review that produced
this section found seven things wrong with the first cut. What the layout means now:

| module | holds |
|---|---|
| `props` | the settings, and every numeric property parsed as f64 then cast |
| `input` | the crosstab: cell count, chunk decode, row labels, column groups, manual cofactors |
| `algorithm` | the transform itself, pure functions over slices |
| `cofactors` | `auto`: the subsample plan, the reservoirs, and the call into `flowvs` |
| `output` | the result envelope, the value/index columns and the cofactor table |
| `lib` | orchestration: which cofactors, which write path, upload |

`tson`, `upload`, `pagecache` and `progress` are **copies** of read_fcs's modules. They are kept
byte-identical on purpose — they diverged within a day of being copied — and they should become
a shared crate when a third operator needs them. That is the one structural debt here.

## Memory

Two paths, chosen by cell count (`output::COLLECT_MAX_CELLS`, 20 M):

- **collect**: one pass, buffering `.ri`, `.ci` and the values — 16 B/cell, so 320 MB at the cap.
- **stream**: also **one** pass. TSON is column-major, so the obvious way to avoid buffering is a
  pass per column, and that was the first implementation. It is the wrong trade: transfer beats
  the arithmetic by about twenty to one, so three passes is the most expensive thing the operator
  could do. Instead the value column is declared first and streamed straight through, while
  `.ri` and `.ci` go to spill files in the exact little-endian layout TSON wants and are poured
  in afterwards. One pass, one chunk of memory, 8 B/cell of scratch disk released as it goes.

Estimation for `auto` never holds more than `cofactors::ESTIMATE_BUDGET_BYTES` (128 MB): the
per-sample cell count is reduced to fit the projection, so a user raising `estimate_max_cells` to
200,000 cannot push the operator past its booking. In collect mode the estimate is taken from the
values already in memory, so only the streaming path spends a pass on it.

Peak therefore does not grow with the projection, so `memory_model.json` books a **constant**
600 MB (`intercept` 420 + 1.5 × `offset` 120) with a zero-exponent feature, because the install
rejects an empty `features` array. **Refit once a real run reports `stats_d_actual_ram_peak`.**

The streamed path declares `nRows` before it has the rows, so `check_count` fails the run if the
table returns a different number of values than its schema promised, rather than writing a result
that decodes into nonsense.

## Why the result is written by hand

`save_table` encodes the whole `OperatorResult` in memory, and this result has exactly one row per
cell. The envelope written here is the same one `tercen-rs` builds for a per-cell result —
`{kind: OperatorResult, tables: [table], joinOperators: []}` — streamed to a file and uploaded in
chunks. `tson.rs`, `upload.rs`, `pagecache.rs` and `progress.rs` are the modules from
`read_fcs_rust_operator`; `tson.rs` gained `i32_list_header`/`i32_chunk` so `.ri`/`.ci` can be
streamed too.

## `auto`, and why the cofactors are an output

`auto` estimates one cofactor per channel with `flowvs` and applies it in the same run, which is
what the R asinh operator does and what the R logicle operator does with `estimateLogicle`. The
difference is that those keep the estimate to themselves. Here it lands in a `Cofactors` table
with Bartlett's statistic, an unstable flag, and the subsample size and seed that produced it.

That table is the point. An automatic estimate silently changes whenever the data flowing through
the step changes, so two timepoints of one study stop being comparable. Reviewing the table and
feeding it back as the cofactor row factor with `method = manual` is how a study gets frozen, and
the plan (`flowvs-rust-plan.md` §11) says to freeze once per study from the batch controls.

Samples come from the column factor named by `sample_factor`, or the first one. It matters
because flowVS pools populations across samples, and in a cytometry projection the columns are
events while the sample is something like `filename`.

`auto` is behind a cargo feature because `flowvs` is a path dependency and the image build copies
this repository alone. Drop the feature once that crate has a remote: shipping an image where a
declared method sometimes exists is worse than either alternative.

## What the first auto run showed about flowVS

On a synthetic eight-channel set (four samples, 5,000 events, a 70/30 mixture per channel), the
same channel estimates a cofactor of **2.0 from 3,000 cells per sample and 99 from 5,000**. The
objective is bimodal, both minima are real, and Bartlett's statistic is small and healthy-looking
at both, so nothing in the output says the answer was a coin toss. CD8 did the opposite: 90 at
both sizes.

Two consequences.

`estimate_max_cells` is a **scientific** parameter, not a performance knob. It belongs in the
Cofactors table, which is why it is written there with the seed.

The `unstable` flag as implemented is too weak. It only fires when flowVS finds fewer than two
usable populations, which is total failure; it does not fire when the search lands in a different
local minimum. The honest diagnostic is to estimate twice on disjoint subsamples and report the
spread, which costs one more pass over the same cells. That is worth doing before anyone freezes
a table from this operator, and it is the guardrail `flowvs-rust-plan.md` §5 gestures at.

This is not a Rust artefact: the crate reproduces the R implementation's published cofactors to
1e-15, and R would swing the same way on the same subsamples.

## The operator spec

`operator.json`'s `operatorSpec` is mirrored from the R operator, including the `axis` entry that
names the y factor, and extended with a **conditional** output: `auto` declares the cofactor
table beside the values, the other methods declare only the values. The condition strings are
written the way the platform tests them (`DataStep._matchesCondition` looks for the property name
and then its value inside the same string).

`tests/r_goldens.rs` checks all of it: every property the code reads is declared, the method enum
matches, both alternatives list exactly the columns the writers produce, and the conditions match
the right methods. Nothing here is data-dependent, so unlike read_fcs no `allowAdditionalAttributes`
is needed.

## Docker

`scratch` + one musl binary. The builder normalises the package version before `cargo chef
prepare`, because the recipe carries that version and a release commit always bumps it: without
this the ~10 minute dependency layer is never reused between releases (measured on read_fcs:
12m50s with a warm cache that the version bump invalidated).

## Next

1. Run it in Studio against a real crosstab and measure peak RSS, then refit the memory model.
2. Container run as uid 1000 under `--memory 600M` on that crosstab.
3. Regenerate `tests/test.json`'s goldens from a Studio run of **this** operator once it installs,
   so the Tercen unit test covers the Rust output rather than the R one.
4. `manual` mode end to end with a cofactor table, which is how the CYTOSHRINK pipeline will use it.
