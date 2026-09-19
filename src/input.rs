//! Reading the crosstab.
//!
//! The projection (R `main.R`): rows are channels, columns are events/observations, y is the
//! value. In `manual` mode a **second** row factor holds the cofactor for each channel
//! (`ctx$rnames[[2]]`), so the row-facet table is fetched once and indexed by `.ri`.
//!
//! The cell table itself is never held whole: `stream_table_chunked` hands back TSON chunks and
//! the caller folds them. A chunk is decoded with `rustson` rather than polars, so a chunk of
//! 1 M cells costs the three column vectors and nothing else.
use anyhow::{Result, anyhow, bail};
use tercen_rs::context::ContextBase;

/// One decoded chunk of the crosstab: the columns asked for, in row order.
#[derive(Debug, Default)]
pub struct Chunk {
    pub ri: Vec<i32>,
    pub ci: Vec<i32>,
    pub y: Vec<f64>,
}

/// Number of cells in the projected crosstab, from the table schema (`Schema.nRows`).
pub async fn cell_count(ctx: &ContextBase) -> Result<usize> {
    let schema = ctx
        .streamer()
        .get_schema(ctx.qt_hash())
        .await
        .map_err(|e| anyhow!("get schema {}: {e}", ctx.qt_hash()))?;
    use tercen_rs::client::proto::e_schema;
    let n = match schema.object {
        Some(e_schema::Object::Schema(s)) => s.n_rows,
        Some(e_schema::Object::Tableschema(s)) => s.n_rows,
        Some(e_schema::Object::Computedtableschema(s)) => s.n_rows,
        Some(e_schema::Object::Cubequerytableschema(s)) => s.n_rows,
        None => bail!("schema {} has no object", ctx.qt_hash()),
    };
    usize::try_from(n).map_err(|_| anyhow!("schema reports {n} rows"))
}

/// Cofactors per row index, for `manual`. `None` when the projection has no second row factor.
pub async fn row_cofactors(ctx: &ContextBase) -> Result<Vec<f64>> {
    let rnames = ctx
        .rnames()
        .await
        .map_err(|e| anyhow!("read row factor names: {e}"))?;
    if rnames.len() < 2 {
        bail!(
            "method 'manual' needs a cofactor row factor after the channel name — the row \
             projection has {:?}",
            rnames
        );
    }
    let name = rnames[1].clone();
    // The row table is one row per channel, so it is small; fetch it whole, and decode it with
    // the same TSON reader as the cell chunks rather than pulling in polars for one column.
    let bytes = ctx
        .streamer()
        .stream_tson(ctx.row_hash(), Some(vec![name.clone()]), 0, -1)
        .await
        .map_err(|e| anyhow!("read row factor '{name}': {e}"))?;
    let out = column_as_f64(&bytes, &name)?;
    for (i, v) in out.iter().enumerate() {
        if !v.is_finite() || *v == 0.0 {
            bail!("cofactor at row {i} is {v}; it must be a non-zero finite number");
        }
    }
    if out.is_empty() {
        bail!("row factor '{name}' has no values");
    }
    Ok(out)
}

/// Read one named column of a TSON table as f64.
pub fn column_as_f64(bytes: &[u8], name: &str) -> Result<Vec<f64>> {
    let v = rustson::decode_bytes(bytes).map_err(|e| anyhow!("decode table: {e:?}"))?;
    let rustson::Value::MAP(m) = v else {
        bail!("table is not a map")
    };
    let Some(rustson::Value::LST(cols)) = m.get("columns") else {
        bail!("table has no columns")
    };
    for c in cols {
        let rustson::Value::MAP(c) = c else { continue };
        let Some(rustson::Value::STR(n)) = c.get("name") else {
            continue;
        };
        if n == name {
            let values = c
                .get("values")
                .ok_or_else(|| anyhow!("column '{name}' has no values"))?;
            return as_f64(values, name);
        }
    }
    bail!("column '{name}' not found in the table")
}

/// Decode one TSON chunk of the cell table into the requested columns.
pub fn decode_chunk(bytes: &[u8], want_ri: bool, want_ci: bool, want_y: bool) -> Result<Chunk> {
    let v = rustson::decode_bytes(bytes).map_err(|e| anyhow!("decode chunk: {e:?}"))?;
    let rustson::Value::MAP(m) = v else {
        bail!("chunk is not a map")
    };
    let rustson::Value::LST(cols) = m
        .get("columns")
        .ok_or_else(|| anyhow!("chunk has no columns"))?
    else {
        bail!("chunk columns is not a list")
    };
    let mut out = Chunk::default();
    for c in cols {
        let rustson::Value::MAP(c) = c else {
            bail!("column is not a map")
        };
        let Some(rustson::Value::STR(name)) = c.get("name") else {
            bail!("column has no name")
        };
        let values = c
            .get("values")
            .ok_or_else(|| anyhow!("column '{name}' has no values"))?;
        match name.as_str() {
            ".ri" if want_ri => out.ri = as_i32(values, name)?,
            ".ci" if want_ci => out.ci = as_i32(values, name)?,
            ".y" if want_y => out.y = as_f64(values, name)?,
            _ => {}
        }
    }
    Ok(out)
}

fn as_i32(v: &rustson::Value, name: &str) -> Result<Vec<i32>> {
    Ok(match v {
        rustson::Value::LSTI32(v) => v.clone(),
        rustson::Value::LSTF64(v) => v.iter().map(|x| *x as i32).collect(),
        rustson::Value::LSTU8(v) => v.iter().map(|x| *x as i32).collect(),
        other => bail!("column '{name}' is not an integer list ({other:?})"),
    })
}

fn as_f64(v: &rustson::Value, name: &str) -> Result<Vec<f64>> {
    Ok(match v {
        rustson::Value::LSTF64(v) => v.clone(),
        rustson::Value::LSTI32(v) => v.iter().map(|x| *x as f64).collect(),
        rustson::Value::LSTU8(v) => v.iter().map(|x| *x as f64).collect(),
        other => bail!("column '{name}' is not a numeric list ({other:?})"),
    })
}
