# asinh_operator

Arcsinh transform for Tercen, with optional cofactor estimation. **Version 2 is a Rust
implementation** that replaces the R one (1.x, kept on the
[`r-legacy`](https://github.com/tercen/asinh_operator/tree/r-legacy) branch and the
`r-legacy-1.2.0` tag). It was developed as `tercen/asinh_rust_operator` and merged here with its
history.

## Changes from 1.x

- **`fixed` and `manual` give the same values as 1.x.** `tests/r_goldens.rs` runs 1.x's own crabs
  test against its expected output (gate 1e-12, measured 2.2e-16). `manual` applies the same rule
  as 1.x, `asinh(y / cofactor)` with the cofactor from the second row factor; 2.0 refuses a zero
  or non-finite cofactor instead of writing Inf/NaN.
- **`auto` can differ from 1.x.** It estimates on a seeded subsample (`estimate_max_cells`, 3,000
  cells per sample, per `sample_factor`) instead of every cell, and by default refuses a
  degenerate cofactor (`cofactor_floor` = 2.5; set 0 for plain flowVS). The flowVS arithmetic
  itself matches R to 2e-16 (`tercen/flowvs-rs`).
- **`auto` adds a second table**, one row per channel: the cofactor used, its Bartlett statistic,
  status, and the alternatives considered, so the estimate can be reviewed and frozen (feed it
  back as the `manual` row factor).
- **New properties:** `sample_factor`, `estimate_max_cells`, `seed`, `threads`, `cofactor_floor`.
- **gRPC operator**, static image; needs a Tercen server with gRPC operator support.

## Input

A crosstab: **rows** are channels, **columns** are events, **y** is the value to transform. The
`manual` method needs a **second row factor** holding each channel's cofactor, exactly as the R
operator does.

## Properties

| Name | Type | Default | Description |
|---|---|---:|---|
| `method` | Enumerated | `fixed` | `fixed` uses `scale` for every channel; `manual` takes a cofactor per channel from the second row factor; `auto` estimates them with flowVS. |
| `scale` | Double | 5 | Cofactor for `fixed`, and the fallback for a channel `auto` cannot resolve. |
| `sample_factor` | String | *(none)* | `auto`: the column factor that identifies a sample, for example the file name. Unset, the whole crosstab is one sample — the operator says so in its log — and a first column factor that is one value per column (the event id) is never used as a sample. |
| `estimate_max_cells` | Double | 3000 | `auto`: cells per sample used for estimation. |
| `seed` | Double | 1 | `auto`: seed for that subsample, so the estimate repeats. |
| `threads` | Double | 4 | `auto`: channels estimated at once. |
| `cofactor_floor` | Double | 2.5 | `auto`: refuse a cofactor below this multiple of the negative population's spread, and use the floor instead. 0 reproduces flowVS, including the failure below. |
| `signifLevel` | Double | 0.05 | `auto`: flowVS peak-detection significance. |
| `bwCorr` | Double | 1.0 | `auto`: flowVS bandwidth correction. |

## Output

One column, `<namespace>.asinh`, one value per cell, carried on `.ri` / `.ci` — the same attribute
name the R operator declares, so a workflow can swap one for the other without re-projecting.

`auto` emits a second table, `Cofactors`: the channel, the cofactor used, Bartlett's statistic, a
status of `resolved`, `fragile`, `floored` or `unstable`, what flowVS itself returned, the
cofactor implied by the negative population's spread, the runner-up from the search, and the
subsample size and seed. **Review it, then freeze it**: feed it back as the cofactor row factor with
`method = manual`. An automatic estimate changes with whatever data flows through the step, so
timepoints estimated separately are not comparable, and the point of a study-wide cofactor table
is that they are.

## Parity

`cargo test` reproduces the R operator's own published golden (`tests/table1.csv`, 1,000 cells
from `crabs-long.csv`) to a worst relative difference of **2.2e-16**. The fixtures and the
projection come from that operator's `test_workflow.json`, so the check needs no Tercen instance.

Cofactor estimation lives in the [`flowvs`](https://github.com/tercen/flowvs-rs) crate, which
reproduces the reference implementation's published cofactors to 1e-15.

### Why the floor is on

flowVS chooses the cofactor that equalises population variances. On an ordinary two-population
channel that search can find a second minimum near a cofactor of 2, where `asinh` has become a
logarithm of the negative population and the equal variances are a coincidence — and it can beat
the real minimum. It is reported `resolved`, because the two objectives are only a factor of 2.5
apart, well outside what the fragile test looks for.

Measured on a synthetic mixture (`the_floor_rescues_a_channel_flowvs_gets_wrong`): flowVS returns
**2.39**, the floor gives **76.0**, and R's published answer on data of the same shape is
**79.87**. So the operator defaults `cofactor_floor` to 2.5 while the crate leaves it off: the
crate's job is to be flowVS, the operator's is to hand a biologist a number they will not check.
The Cofactors table marks such a channel `floored` and keeps what flowVS returned beside it.

## Size and speed

The crosstab is never held whole. Under 20 M cells the operator buffers one pass; above that it
streams, writing the value column as it reads and spilling the two index columns to disk, so it
still makes **one** pass over the input. Estimation is bounded to 128 MB of subsample whatever
the properties ask for, and runs across channels in parallel.

`auto` needs the `auto` cargo feature, which is **on** by default and pins
`tercen/flowvs-rs` at a tag — the cofactors it returns are output, so they must not move under a
rebuild. Build with `--no-default-features` to skip the estimator entirely.
