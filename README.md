# asinh_rust_operator

Arcsinh transform for Tercen, with optional cofactor estimation. Rust port of
[`tercen/asinh_operator`](https://github.com/tercen/asinh_operator).

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
