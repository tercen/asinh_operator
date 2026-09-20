//! `method = auto`: estimate one cofactor per channel with flowVS, then transform with it.
//!
//! Two decisions shape this.
//!
//! **It estimates from a subsample.** The search evaluates its objective around a hundred times
//! per channel, each pass running a kernel density estimate over every value. On a cohort-scale
//! crosstab that is not something to do by accident, and flowVS practice is a few thousand cells
//! per sample anyway. The subsample is seeded, and both the count and the seed are written into
//! the output so the estimate can be repeated exactly.
//!
//! **The estimate is an output, not a secret.** It lands in a `Cofactors` table beside the
//! transformed values, with Bartlett's statistic and an unstable flag, so it can be reviewed and
//! then frozen: feed that table back as the cofactor row factor with `method = manual`, and the
//! run stops depending on which cells happened to flow through it.
use std::collections::HashMap;

use anyhow::{Result, bail};
use tercen_rs::context::ContextBase;

use crate::input;
use crate::output::CofactorRow;
use crate::props::Settings;

/// Estimation never holds more than this, whatever the properties say. 128 MB against a 600 MB
/// booking leaves room for the transform.
pub const ESTIMATE_BUDGET_BYTES: usize = 128 * 1000 * 1000;
/// Below this a sample has noise, not populations.
pub const MIN_CELLS_PER_SAMPLE: usize = 500;

/// Cells kept per (channel, sample) so the whole subsample fits the budget.
///
/// It depends on the projection, not on the machine, so two runs of the same step keep the same
/// cells and produce the same cofactors.
pub fn cells_per_sample(asked: usize, n_channels: usize, n_samples: usize) -> usize {
    let groups = n_channels.max(1) * n_samples.max(1);
    let affordable = ESTIMATE_BUDGET_BYTES / (groups * std::mem::size_of::<f64>());
    asked.min(affordable).max(MIN_CELLS_PER_SAMPLE)
}

/// Deterministic per-group stream: the same (seed, channel, sample) always keeps the same cells.
struct Xorshift(u64);

impl Xorshift {
    fn new(seed: u64, a: usize, b: usize) -> Self {
        // Mix the three so neighbouring groups do not share a stream.
        let mut s = seed
            .wrapping_mul(0x9E37_79B9_7F4A_7C15)
            .wrapping_add((a as u64).wrapping_mul(0xBF58_476D_1CE4_E5B9))
            .wrapping_add((b as u64).wrapping_mul(0x94D0_49BB_1331_11EB));
        if s == 0 {
            s = 0x2545_F491_4F6C_DD1D;
        }
        Self(s)
    }
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
}

/// A reservoir per (channel, sample), so one pass over the crosstab is enough.
struct Reservoir {
    keep: Vec<f64>,
    seen: usize,
    cap: usize,
    rng: Xorshift,
}

impl Reservoir {
    fn new(cap: usize, seed: u64, a: usize, b: usize) -> Self {
        Self {
            keep: Vec::with_capacity(cap.min(4096)),
            seen: 0,
            cap,
            rng: Xorshift::new(seed, a, b),
        }
    }
    fn offer(&mut self, v: f64) {
        if !v.is_finite() {
            return;
        }
        self.seen += 1;
        if self.keep.len() < self.cap {
            self.keep.push(v);
        } else {
            let j = (self.rng.next() % self.seen as u64) as usize;
            if j < self.cap {
                self.keep[j] = v;
            }
        }
    }
}

/// What the estimation pass produced.
pub struct Estimated {
    /// One cofactor per row index, ready for the transform.
    pub per_row: Vec<f64>,
    /// The table that goes into the result.
    pub table: Vec<CofactorRow>,
}

/// Estimate cofactors for every channel in the projection.
///
/// Channels come from the first row factor; samples from the column factor named by the
/// `sample_factor` property, or the first one. R uses each
/// column as its own sample when no column factor is present; with event-level columns that
/// leaves one value per group and nothing to estimate, so here the whole crosstab is treated as
/// one sample instead, which is the useful reading of the same situation.
pub async fn estimate(ctx: &ContextBase, s: &Settings) -> Result<Estimated> {
    let plan = Plan::new(ctx, s).await?;
    let n_cells = input::cell_count(ctx).await?;
    let mut res: HashMap<(usize, usize), Reservoir> = HashMap::new();
    input::for_each_chunk(ctx, &[".ri", ".ci", ".y"], n_cells, 1_000_000, |c| {
        plan.offer(&mut res, &c.ri, &c.ci, &c.y, s.seed);
        Ok(())
    })
    .await?;
    Ok(plan.finish(res, s))
}

/// The same estimate from values already in memory — the collect path has them, so reading the
/// crosstab twice would be waste, and on a small projection that second pass is most of the run.
pub async fn estimate_from_cells(
    ctx: &ContextBase,
    s: &Settings,
    ri: &[i32],
    ci: &[i32],
    y: &[f64],
) -> Result<Estimated> {
    let plan = Plan::new(ctx, s).await?;
    let mut res: HashMap<(usize, usize), Reservoir> = HashMap::new();
    plan.offer(&mut res, ri, ci, y, s.seed);
    Ok(plan.finish(res, s))
}

/// What estimation needs to know before it sees a single value.
struct Plan {
    channel_names: Vec<String>,
    sample_of_col: Vec<usize>,
    n_channels: usize,
    n_samples: usize,
    cells: usize,
}

impl Plan {
    async fn new(ctx: &ContextBase, s: &Settings) -> Result<Self> {
        if s.estimate_max_cells < MIN_CELLS_PER_SAMPLE {
            bail!(
                "property 'estimate_max_cells' is {}, and flowVS needs at least {} cells per \
                 sample to find populations",
                s.estimate_max_cells,
                MIN_CELLS_PER_SAMPLE
            );
        }
        let channel_names = input::row_labels(ctx).await?;
        let sample_of_col = input::column_groups(ctx, &s.sample_factor).await?;
        let n_channels = channel_names.len();
        let n_samples = sample_of_col
            .iter()
            .copied()
            .max()
            .map(|m| m + 1)
            .unwrap_or(1);
        // The subsample is what estimation costs in memory, and `estimate_max_cells` is a
        // property, so an operator that trusted it could be asked for gigabytes and be killed.
        // Bound the total rather than the property: the booking holds whatever a user types.
        let cells = cells_per_sample(s.estimate_max_cells, n_channels, n_samples);
        if cells < s.estimate_max_cells {
            tracing::warn!(
                asked = s.estimate_max_cells,
                using = cells,
                n_channels,
                n_samples,
                budget_mb = ESTIMATE_BUDGET_BYTES / 1_000_000,
                "estimate_max_cells reduced to stay inside the estimation memory budget"
            );
        }
        tracing::info!(
            n_channels,
            n_samples,
            cells_per_sample = cells,
            "estimating cofactors (flowVS)"
        );
        Ok(Self {
            channel_names,
            sample_of_col,
            n_channels,
            n_samples,
            cells,
        })
    }

    fn offer(
        &self,
        res: &mut HashMap<(usize, usize), Reservoir>,
        ri: &[i32],
        ci: &[i32],
        y: &[f64],
        seed: u64,
    ) {
        for k in 0..y.len() {
            let r = ri[k].max(0) as usize;
            let c = ci[k].max(0) as usize;
            let sample = self.sample_of_col.get(c).copied().unwrap_or(0);
            res.entry((r, sample))
                .or_insert_with(|| Reservoir::new(self.cells, seed, r, sample))
                .offer(y[k]);
        }
    }

    fn finish(&self, res: HashMap<(usize, usize), Reservoir>, s: &Settings) -> Estimated {
        let channels: Vec<Vec<Vec<f64>>> = (0..self.n_channels)
            .map(|ri| {
                (0..self.n_samples)
                    .filter_map(|sa| res.get(&(ri, sa)).map(|r| r.keep.clone()))
                    .filter(|v| v.len() >= 10)
                    .collect()
            })
            .collect();

        let opts = flowvs::estimate::Options {
            signif_level: s.signif_level,
            bw_corr: s.bw_corr,
            threads: s.threads,
        };
        let t = std::time::Instant::now();
        let est = flowvs::estimate::estimate_cofactors(&channels, opts);
        tracing::info!(
            secs = format!("{:.1}", t.elapsed().as_secs_f64()),
            threads = s.threads,
            "cofactors estimated"
        );

        let mut per_row = Vec::with_capacity(self.n_channels);
        let mut table = Vec::with_capacity(self.n_channels);
        for (ri, e) in est.iter().enumerate() {
            let unstable = e.objective >= flowvs::estimate::MAX_BT || !e.cofactor.is_finite();
            // An unstable channel gets the fixed scale rather than a number flowVS did not
            // really find; the flag and the statistic say so in the output.
            let cofactor = if unstable || e.cofactor <= 0.0 {
                s.scale
            } else {
                e.cofactor
            };
            if unstable {
                tracing::warn!(
                    channel = %self.channel_names[ri],
                    "flowVS found no usable populations; falling back to scale = {}",
                    s.scale
                );
            }
            per_row.push(cofactor);
            table.push(CofactorRow {
                channel: self.channel_names[ri].clone(),
                cofactor,
                objective: e.objective,
                unstable,
                cells_used: self.cells.min(i32::MAX as usize) as i32,
                seed: s.seed.min(i32::MAX as u64) as i32,
            });
        }
        Estimated { per_row, table }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_subsample_stays_inside_its_budget() {
        // a cohort-scale projection: the ask is cut down
        let cells = cells_per_sample(200_000, 43, 93);
        assert!(cells < 200_000);
        assert!(43 * 93 * cells * 8 <= ESTIMATE_BUDGET_BYTES);
        // a small one: the ask is honoured
        assert_eq!(cells_per_sample(3000, 20, 10), 3000);
    }

    #[test]
    fn a_tiny_budget_never_goes_below_the_floor() {
        assert_eq!(cells_per_sample(3000, 10_000, 10_000), MIN_CELLS_PER_SAMPLE);
    }

    #[test]
    fn the_subsample_is_the_same_on_every_run() {
        let take = |seed: u64| {
            let mut r = Reservoir::new(50, seed, 3, 7);
            for i in 0..5000 {
                r.offer(i as f64);
            }
            r.keep
        };
        assert_eq!(take(1), take(1));
        assert_ne!(take(1), take(2));
    }
}
