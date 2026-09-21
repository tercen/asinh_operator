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

## Loading the task (0.1.1)

`src/context.rs` builds the context from the **task** and never fetches the workflow.

`ProductionContext::from_task_id` does fetch it, to pull colour and palette settings off the
step, and it fails the whole run when the step is not in the saved workflow document:
`Step '…' not found in workflow`. That happens intermittently in normal use — change a property,
run, and the task can reference a step the saved workflow has not caught up with; a rerun
usually succeeds. Faris hit it on tercen.com on 2026-09-20, at cofactor 10 and again at 7, with
successful runs in between.

Everything a transform needs is on the task: the `CubeQuery` carries the table hashes and the
operator settings, and the schema ids are on the task or on the `CubeQueryTask` that produced
it. Colours, palettes, chart kind and axis tables are for plot operators.

`read_fcs_rust_operator` still calls `from_task_id` and has the same intermittent failure
waiting for it; the same module should move there. The general fix belongs upstream in
`tercen-rs`, where colour extraction should warn rather than fail an operator that never asked
for colours.

## Reading a crosstab: one response is many documents (0.1.2)

`TableStreamer::stream_tson` concatenates every gRPC message of a response, and the server pages
its answer — about 15,000 rows a page on Studio. So the buffer holds **one complete TSON document
per page, back to back**, and a decoder that reads the first one silently discards the rest.

That was pathological rather than merely wrong. Asking for a million rows transferred a million
rows, used fifteen thousand, then asked again from a slightly later offset, so the same data
crossed the wire dozens of times: 800 bytes on the wire per cell of three columns that need
sixteen.

`input::decode_chunk` now reads every document in the buffer. On a 19.53 M-cell crosstab
(93 files × 5,000 events × 43 channels, the read_fcs output):

| | before | after |
|---|---|---|
| read, transform and write | 646.9 s | **8.2 s** |
| end to end | 653.3 s | **17.5 s** |
| peak RSS | 61.5 MB | 35.5 MB |
| wire cost | ~800 B/cell | 16 B/cell, exactly the columns |

`CHUNK` is 200,000 because that measured fastest (3.0 M cells/s, against 0.8 M at 15,000 and
2.4 M at 1,000,000). The reader is bounded by the schema's row count either way.

Decoding advances the cursor by re-encoding each document, because `rustson::decode` takes the
cursor by value and reports nothing. A `decode_from(&mut Cursor)` upstream would remove that.

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

`auto` is a cargo feature, **on by default**, pinning `tercen/flowvs-rs` at a tag. It was off
while that crate was local-only, because cargo resolves a `path` dependency even when its feature
is off and the image build copies this repository alone — which is why the Dockerfile and CI used
to strip the line. Both strips are gone. Keep the pin on a tag rather than a branch: the cofactors
the crate returns are this operator's output, and output must not move under a rebuild.

`cofactor_floor` defaults to **2.5** here while the crate defaults it off. The crate's job is to be
flowVS; the operator's is to hand a biologist a number they will not check, and flowVS can prefer a
degenerate cofactor near 2 on an ordinary channel and call it `resolved`. `README.md` has the
measurement. Setting the property to 0 reproduces flowVS exactly, failure included.

## How stable is a flowVS cofactor? (measured 2026-09-20)

Worth knowing before anyone freezes a table, and the answer differs sharply between a channel
with a clear positive population and one without.

**Real data** (`flowvs_input.csv`, 3 channels, 12 samples, ~5,000 cells each) resampled:

| | CD4 | CD8 | CD3 |
|---|---|---|---|
| 1,000 → 5,000 cells per sample | 6723 → 6389 | 4937 → 4625 | 5361 → 6147 |
| spread across sizes | 6% | 20% | 15% |
| spread across three seeds at 3,000 | 4% | 12% | 10% |

So on a well-behaved channel the estimate wobbles by roughly **5 to 20 %**, whichever way it is
resampled. That is the same order as the tolerance `flowvs-rust-plan.md` §3 already accepts, and
as the 12.7 % it records between the R port and the C original. It is invisible on a histogram.
Freezing a table is fine, provided the subsample size and seed are recorded — they are, in the
Cofactors table — and kept fixed for the study.

**Synthetic data with a very tight negative population** behaves completely differently: the same
channel gave a cofactor of **2.0 from 3,000 cells and 99 from 5,000**, both genuine minima of the
objective, both with a small Bartlett statistic. The small answer is the degenerate one, where
asinh is effectively a log of the noise and the variances match for the wrong reason.

That is not only a synthetic curiosity: it is the dim-channel failure the cofactor check already
saw on the study panel, where flowVS collapsed on TCRγδ, CD56 and CD45. The principled fix is the
negative-spread floor in `flowvs-rust-plan.md` §5 — a cofactor below the spread of the negative
population cannot be stabilising anything — and it belongs in the crate, not here.

Both guards are now implemented in the crate and surfaced here.

**The floor** (`cofactor_floor`, a multiple of σ_neg, **0 = off**) refuses a cofactor below the
negative population's spread and uses `factor · σ_neg` instead. Checked twice: on the real
three-channel case the σ_neg values come out at 3908, 4019 and 2917 against flowVS's 6450, 4767
and 6317, so nothing is floored and the published answers stand; on the synthetic set every
channel was floored from about 1–4 up to about 75, which is 2.5 × the 30-unit negative spread
that data was generated with.

**The runner-up** is free — the search computes a best cofactor per interval and kept only the
winner — but weaker than it sounds. On the degenerate synthetic case it reported *resolved*,
because the second-best interval was genuinely worse. It catches a different failure: two
comparable minima far apart.

A channel now carries a `status` of `resolved`, `fragile`, `floored` or `unstable`, and the
Cofactors table carries what flowVS said, the σ_neg cofactor and the runner-up beside the number
that was used. Anything but `resolved` is also a warning in the task log.

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

## An extra table needs `.ri` or `.ci`; never declare a join on an empty key

Every table in `tables` must carry `.ri` and/or `.ci`, and the server joins it by those — a
`.ri`-only table is a per-row annotation, which is what the Cofactors table is. Declaring a
`JoinOperator` with an empty `ColumnPair` instead (copied from an import operator, where there
is no crosstab to join against) gave a composite that no downstream step could query. The only
test that sees this is the platform's `OperatorUnitTest`, because it diffs the assembled
relations; `tests/asinh_auto.json` exists for that reason. Measured on Studio, 2026-09-21.
