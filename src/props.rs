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
}

#[derive(Debug, Clone)]
pub struct Settings {
    pub method: Method,
    pub scale: f64,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            method: Method::Fixed,
            scale: 5.0,
        }
    }
}

/// Read the properties off the task's `CubeQueryTask` snapshot.
pub fn settings_from_ctx(ctx: &ContextBase) -> Result<Settings> {
    let pr = PropertyReader::from_operator_settings(ctx.operator_settings());
    // Tercen serialises a numeric property as e.g. "5.0" even where the operator wants an
    // integer, so parse every number as f64 and cast. Parsing as i32 makes a legitimate "5.0"
    // fall back to the default *silently*, which is worse than failing (create-rust-operator §2).
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
        "auto" => bail!(
            "method 'auto' estimates cofactors with flowVS, which this Rust port does not \
             implement yet — use the R `asinh_operator` for 'auto', or 'manual' here with a \
             cofactor row factor (see flowvs-rust-plan.md)"
        ),
        other => bail!("property 'method' must be 'fixed', 'manual' or 'auto', got '{other}'"),
    };
    Ok(Settings { method, scale })
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
