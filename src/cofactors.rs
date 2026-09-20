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

use anyhow::{Result, anyhow, bail};
use tercen_rs::context::ContextBase;

use crate::input;
use crate::output::CofactorRow;
use crate::props::Settings;

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
/// Channels come from the first row factor and samples from the first column factor. R uses each
/// column as its own sample when no column factor is present; with event-level columns that
/// leaves one value per group and nothing to estimate, so here the whole crosstab is treated as
/// one sample instead, which is the useful reading of the same situation.
pub async fn estimate(ctx: &ContextBase, s: &Settings) -> Result<Estimated> {
    if s.estimate_max_cells < 100 {
        bail!(
            "property 'estimate_max_cells' is {} — flowVS needs at least a few hundred cells per \
             sample to find populations",
            s.estimate_max_cells
        );
    }
    let channel_names = input::row_labels(ctx).await?;
    let sample_of_col = input::column_groups(ctx).await?;
    let n_channels = channel_names.len();
    let n_samples = sample_of_col
        .iter()
        .copied()
        .max()
        .map(|m| m + 1)
        .unwrap_or(1);
    tracing::info!(
        n_channels,
        n_samples,
        cells_per_sample = s.estimate_max_cells,
        "estimating cofactors (flowVS)"
    );

    let mut res: HashMap<(usize, usize), Reservoir> = HashMap::new();
    let mut err: Option<anyhow::Error> = None;
    ctx.streamer()
        .stream_table_chunked(
            ctx.qt_hash(),
            Some(vec![".ri".into(), ".ci".into(), ".y".into()]),
            1_000_000,
            |bytes| {
                match input::decode_chunk(&bytes, true, true, true) {
                    Ok(c) => {
                        for k in 0..c.y.len() {
                            let ri = c.ri[k].max(0) as usize;
                            let ci = c.ci[k].max(0) as usize;
                            let sample = sample_of_col.get(ci).copied().unwrap_or(0);
                            res.entry((ri, sample))
                                .or_insert_with(|| {
                                    Reservoir::new(s.estimate_max_cells, s.seed, ri, sample)
                                })
                                .offer(c.y[k]);
                        }
                    }
                    Err(e) => err = Some(e),
                }
                Ok(())
            },
        )
        .await
        .map_err(|e| anyhow!("stream the crosstab for estimation: {e}"))?;
    if let Some(e) = err {
        return Err(e);
    }

    let channels: Vec<Vec<Vec<f64>>> = (0..n_channels)
        .map(|ri| {
            (0..n_samples)
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

    let mut per_row = Vec::with_capacity(n_channels);
    let mut table = Vec::with_capacity(n_channels);
    for (ri, e) in est.iter().enumerate() {
        let unstable = e.objective >= flowvs::estimate::MAX_BT || !e.cofactor.is_finite();
        // An unstable channel gets the fixed scale rather than a number flowVS did not really
        // find; the flag and the statistic say so in the output.
        let cofactor = if unstable || e.cofactor <= 0.0 {
            s.scale
        } else {
            e.cofactor
        };
        if unstable {
            tracing::warn!(
                channel = %channel_names[ri],
                "flowVS found no usable populations; falling back to scale = {}",
                s.scale
            );
        }
        per_row.push(cofactor);
        table.push(CofactorRow {
            channel: channel_names[ri].clone(),
            cofactor,
            objective: e.objective,
            unstable,
            cells_used: s.estimate_max_cells.min(i32::MAX as usize) as i32,
            seed: s.seed.min(i32::MAX as u64) as i32,
        });
    }
    Ok(Estimated { per_row, table })
}
