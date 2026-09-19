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

## Memory

Two paths, chosen by cell count (`output::COLLECT_MAX_CELLS`, 20 M):

- **collect**: one pass, buffering `.ri`, `.ci` and the values — 16 B/cell, so 320 MB at the cap.
- **stream**: three passes, one per output column, holding one chunk (1 M cells). The `.ri`/`.ci`
  passes fetch a single column each; the value pass fetches `.ri` and `.y` because a per-channel
  cofactor needs the row index.

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

## The operator spec

`operator.json`'s `operatorSpec` is mirrored from the R operator so both declare the same shape,
including the `axis` entry that names the y factor. `tests/r_goldens.rs` asserts the declared
output attribute and the column the writer produces have not drifted apart. Nothing here is
data-dependent, so — unlike read_fcs — no `allowAdditionalAttributes` is needed.

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
