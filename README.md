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
| `sample_factor` | String | *(first)* | `auto`: which column factor identifies the sample. |
| `estimate_max_cells` | Double | 3000 | `auto`: cells per sample used for estimation. |
| `seed` | Double | 1 | `auto`: seed for that subsample, so the estimate repeats. |
| `threads` | Double | 4 | `auto`: channels estimated at once. |
| `signifLevel` | Double | 0.05 | `auto`: flowVS peak-detection significance. |
| `bwCorr` | Double | 1.0 | `auto`: flowVS bandwidth correction. |

## Output

One column, `<namespace>.asinh`, one value per cell, carried on `.ri` / `.ci` — the same attribute
name the R operator declares, so a workflow can swap one for the other without re-projecting.

`auto` emits a second table, `Cofactors`: the channel, the cofactor used, Bartlett's statistic, an
unstable flag for channels where flowVS found nothing to stabilise, and the subsample size and
seed. **Review it, then freeze it**: feed it back as the cofactor row factor with
`method = manual`. An automatic estimate changes with whatever data flows through the step, so
timepoints estimated separately are not comparable, and the point of a study-wide cofactor table
is that they are.

## Parity

`cargo test` reproduces the R operator's own published golden (`tests/table1.csv`, 1,000 cells
from `crabs-long.csv`) to a worst relative difference of **2.2e-16**. The fixtures and the
projection come from that operator's `test_workflow.json`, so the check needs no Tercen instance.

Cofactor estimation lives in the `flowvs` crate, which reproduces the reference implementation's
published cofactors to 1e-15.

## Size and speed

The crosstab is never held whole. Under 20 M cells the operator buffers one pass; above that it
streams, writing the value column as it reads and spilling the two index columns to disk, so it
still makes **one** pass over the input. Estimation is bounded to 128 MB of subsample whatever
the properties ask for, and runs across channels in parallel.

`auto` needs the `auto` cargo feature, which is off by default until the `flowvs` crate has a
git remote.
