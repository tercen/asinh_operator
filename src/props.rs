//! Operator properties, with the R operator's defaults (`asinh_operator/main.R`).
//!
//! `ctx$op.value("method", default = "fixed")`, `scale` (integer, 5), and — for the `auto`
//! method only — `signifLevel` and `bwCorr`, which this port does not implement (see `Method`).
use anyhow::{Result, bail};
use tercen_rs::PropertyReader;
use tercen_rs::context::ContextBase;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    /// One cofactor for every cell, from the `scale` property.
    Fixed,
    /// A cofactor per channel, taken from the **second** row factor.
    Manual,
    /// Estimate a cofactor per channel with flowVS, then apply it. The estimate is emitted as a
    /// second output table so it can be reviewed, frozen and fed back in as `manual`.
    ///
    /// Only present with the `auto` feature, which pulls in the `flowvs` crate.
    #[cfg(feature = "auto")]
    Auto,
}

#[derive(Debug, Clone)]
pub struct Settings {
    pub method: Method,
    pub scale: f64,
    /// Cells per sample used for estimation (`auto`). flowVS practice is a few thousand; the
    /// whole crosstab would mean a hundred objective evaluations per channel over every event.
    pub estimate_max_cells: usize,
    /// Seed for that subsample, so an estimate is repeatable.
    pub seed: u64,
    /// Channels estimated at once. Peak memory grows with it, so it is explicit, not "all cores".
    pub threads: usize,
    /// flowVS significance level for peak detection.
    pub signif_level: f64,
    /// flowVS bandwidth correction (1.0 = the R port's default).
    pub bw_corr: f64,
    /// Which column factor identifies the sample for estimation. Empty means the first one.
    ///
    /// It matters: flowVS pools populations across samples, and in a cytometry projection the
    /// columns are events while the sample is a factor like `filename`, which need not be first.
    pub sample_factor: String,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            method: Method::Fixed,
            scale: 5.0,
            estimate_max_cells: 3000,
            seed: 1,
            threads: 4,
            signif_level: 0.05,
            bw_corr: 1.0,
            sample_factor: String::new(),
        }
    }
}

/// Read the properties off the task's `CubeQueryTask` snapshot.
pub fn settings_from_ctx(ctx: &ContextBase) -> Result<Settings> {
    let pr = PropertyReader::from_operator_settings(ctx.operator_settings());
    // Tercen serialises a numeric property as e.g. "5.0" even where the operator wants an
    // integer, so parse every number as f64 and cast. Parsing as i32 makes a legitimate "5.0"
    // fall back to the default *silently*, which is worse than failing (create-rust-operator §2).
    // Every numeric property is parsed as f64 and cast (see below).
    let num = |name: &str, default: f64| -> Result<f64> {
        let raw = pr.get_string(name, &default.to_string());
        raw.trim()
            .parse::<f64>()
            .map_err(|_| anyhow::anyhow!("property '{name}' is not a number: '{raw}'"))
    };
    let raw_scale = pr.get_string("scale", "5");
    let scale: f64 = raw_scale
        .trim()
        .parse()
        .map_err(|_| anyhow::anyhow!("property 'scale' is not a number: '{raw_scale}'"))?;
    if !(scale.is_finite() && scale != 0.0) {
        bail!("property 'scale' must be a non-zero finite number, got {scale}");
    }
    let method = match pr
        .get_string("method", "fixed")
        .trim()
        .to_ascii_lowercase()
        .as_str()
    {
        "fixed" => Method::Fixed,
        "manual" => Method::Manual,
        "auto" => {
            #[cfg(feature = "auto")]
            {
                Method::Auto
            }
            #[cfg(not(feature = "auto"))]
            {
                bail!(
                    "method 'auto' estimates cofactors with flowVS; this build was compiled \
                     without the `auto` feature. Use 'manual' with a cofactor row factor, or \
                     the R `asinh_operator`."
                )
            }
        }
        other => bail!("property 'method' must be 'fixed', 'manual' or 'auto', got '{other}'"),
    };
    let estimate_max_cells = num("estimate_max_cells", 3000.0)?;
    let seed = num("seed", 1.0)?;
    let threads = num("threads", 4.0)?;
    let signif_level = num("signifLevel", 0.05)?;
    let bw_corr = num("bwCorr", 1.0)?;
    if !(0.0..1.0).contains(&signif_level) {
        bail!("property 'signifLevel' must be in [0, 1), got {signif_level}");
    }
    if bw_corr <= 0.0 {
        bail!("property 'bwCorr' must be > 0, got {bw_corr}");
    }
    Ok(Settings {
        method,
        scale,
        estimate_max_cells: estimate_max_cells.max(0.0) as usize,
        seed: seed.max(0.0) as u64,
        threads: threads.max(0.0) as usize,
        signif_level,
        bw_corr,
        sample_factor: pr.get_string("sample_factor", "").trim().to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_match_the_r_operator() {
        let s = Settings::default();
        assert_eq!(s.method, Method::Fixed);
        assert_eq!(s.scale, 5.0);
    }
}
