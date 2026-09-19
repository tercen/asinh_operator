//! The result: one table, `.ri` / `.ci` / `<namespace>.asinh`, one row per cell.
//!
//! Shape is what `tercen-rs`'s `save_table` builds for a per-cell result — `{kind:
//! OperatorResult, tables: [table], joinOperators: []}` — but written as a **stream**, because
//! the result has exactly as many rows as the projected crosstab has cells and `save_table`
//! encodes the whole thing in memory first (create-rust-operator §4).
//!
//! TSON is column-major, so a column has to be written from start to end before the next one
//! begins. Two ways to do that:
//!
//! * **collect** (default under [`COLLECT_MAX_CELLS`]): one pass over the crosstab, keeping
//!   `.ri`, `.ci` and the transformed values (16 B/cell), then write.
//! * **stream**: three passes over the crosstab, one per column, keeping only the current chunk.
//!   Three times the server round trips, constant memory. Used above the threshold.
use std::io::Write;

use anyhow::{Result, anyhow};

use crate::tson::TsonWriter;

/// 20 M cells ≈ 320 MB of buffers, which fits comfortably under a 1 GB booking.
pub const COLLECT_MAX_CELLS: usize = 20_000_000;

pub struct ColSpec<'a> {
    pub name: &'a str,
    pub ty: &'a str, // "double" | "int32"
}

/// The three columns of the result, in the order R saves them.
pub fn result_columns(namespace: &str) -> [ColSpec<'_>; 3] {
    [
        ColSpec {
            name: ".ri",
            ty: "int32",
        },
        ColSpec {
            name: ".ci",
            ty: "int32",
        },
        ColSpec {
            name: namespace,
            ty: "double",
        },
    ]
}

/// The result column's name: `<namespace>.asinh` (R `ctx$addNamespace()`).
pub fn value_column(namespace: &str) -> String {
    format!("{namespace}.asinh")
}

pub fn write_header<W: Write>(
    w: &mut TsonWriter<W>,
    table_name: &str,
    n_rows: usize,
    cols: &[ColSpec],
) -> Result<()> {
    w.map(3)?;
    w.key("kind")?;
    w.str("OperatorResult")?;
    w.key("tables")?;
    w.list(1)?;

    w.map(4)?;
    w.key("kind")?;
    w.str("Table")?;
    w.key("nRows")?;
    w.i32(i32::try_from(n_rows).map_err(|_| {
        anyhow!(
            "the result would have {n_rows} rows, more than a Tercen table can hold (i32::MAX). \
             Project fewer cells, or split the step."
        )
    })?)?;
    w.key("properties")?;
    w.map(4)?;
    w.key("kind")?;
    w.str("TableProperties")?;
    w.key("name")?;
    w.str(table_name)?;
    w.key("sortOrder")?;
    w.list(0)?;
    w.key("ascending")?;
    w.bool(false)?;
    w.key("columns")?;
    w.list(cols.len())?;
    Ok(())
}

pub fn write_column_header<W: Write>(
    w: &mut TsonWriter<W>,
    c: &ColSpec,
    n_rows: usize,
) -> Result<()> {
    w.map(6)?;
    w.key("kind")?;
    w.str("Column")?;
    w.key("name")?;
    w.str(c.name)?;
    w.key("type")?;
    w.str(c.ty)?;
    w.key("nRows")?;
    w.i32(n_rows as i32)?;
    w.key("size")?;
    w.i32(n_rows as i32)?;
    w.key("values")?;
    Ok(())
}

/// Close the result after the last column (the empty join list).
pub fn write_footer<W: Write>(w: &mut TsonWriter<W>) -> Result<()> {
    w.key("joinOperators")?;
    w.list(0)?;
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_value_column_carries_the_namespace() {
        assert_eq!(value_column("ds0"), "ds0.asinh");
    }

    #[test]
    fn a_result_bigger_than_an_i32_is_an_error_not_a_panic() {
        let mut w = TsonWriter::new(Vec::new()).unwrap();
        let cols = result_columns("ds0.asinh");
        let e = write_header(&mut w, "t", i32::MAX as usize + 1, &cols).unwrap_err();
        assert!(e.to_string().contains("more than a Tercen table can hold"));
    }
}
