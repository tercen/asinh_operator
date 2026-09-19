//! Parity against the **R operator's own goldens**, not goldens invented here.
//!
//! `tercen/asinh_operator/tests` ships `crabs-long.csv` (the input), `table2.csv` and
//! `table3.csv` (the column and row facet tables the platform derives from it) and
//! `table1.csv` (the result it published). `test_workflow.json` pins the projection: columns =
//! `observation`, rows = `variable`, y = `measurement`, no properties, so `method` is `fixed`
//! and `scale` is 5.
//!
//! `.ri` indexes the row facet table in file order, `.ci` the column facet table, which is what
//! makes the golden reproducible without a Tercen instance: build the same crosstab, transform
//! it, and compare cell by cell.
use asinh_operator::algorithm::asinh_scaled;
use std::collections::HashMap;

/// Minimal CSV reader for these fixtures: one header line, comma separated, optional quotes,
/// no embedded commas.
fn read_csv(path: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    let mut lines = text.lines();
    let unquote = |s: &str| s.trim().trim_matches('"').to_string();
    let header: Vec<String> = lines.next().unwrap().split(',').map(unquote).collect();
    let rows = lines
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split(',').map(unquote).collect())
        .collect();
    (header, rows)
}

fn column(path: &str, name: &str) -> Vec<String> {
    let (header, rows) = read_csv(path);
    let i = header
        .iter()
        .position(|h| h == name)
        .unwrap_or_else(|| panic!("{path} has no column {name}, only {header:?}"));
    rows.into_iter().map(|r| r[i].clone()).collect()
}

#[test]
fn fixed_method_reproduces_the_r_operator_golden() {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests");
    // the crosstab, keyed the way Tercen indexes it
    let (header, rows) = read_csv(&format!("{dir}/crabs-long.csv"));
    let col = |n: &str| header.iter().position(|h| h == n).unwrap();
    let (i_obs, i_var, i_val) = (col("observation"), col("variable"), col("measurement"));
    let mut cells: HashMap<(String, String), f64> = HashMap::new();
    for r in &rows {
        cells.insert(
            (r[i_var].clone(), r[i_obs].clone()),
            r[i_val].parse().unwrap(),
        );
    }
    let variables = column(&format!("{dir}/table3.csv"), "variable"); // .ri order
    let observations = column(&format!("{dir}/table2.csv"), "observation"); // .ci order
    assert_eq!(variables.len(), 5);
    assert_eq!(observations.len(), 200);

    let (gh, golden) = read_csv(&format!("{dir}/table1.csv"));
    let (g_ri, g_ci, g_v) = (
        gh.iter().position(|h| h == ".ri").unwrap(),
        gh.iter().position(|h| h == ".ci").unwrap(),
        gh.iter().position(|h| h == "ds0.asinh").unwrap(),
    );
    assert_eq!(golden.len(), 1000, "the golden should cover every cell");

    let mut worst = 0.0f64;
    for g in &golden {
        let ri: usize = g[g_ri].parse().unwrap();
        let ci: usize = g[g_ci].parse().unwrap();
        let want: f64 = g[g_v].parse().unwrap();
        let y = cells[&(variables[ri].clone(), observations[ci].clone())];
        let got = asinh_scaled(y, 5.0); // scale property default, as the test json leaves it unset
        let rel = (got - want).abs() / want.abs().max(1.0);
        worst = worst.max(rel);
    }
    assert!(
        worst <= 1e-12,
        "worst relative difference from the R operator's golden is {worst:e}"
    );
    println!("R golden parity: 1000 cells, worst relative difference {worst:e}");
}

#[test]
fn manual_method_applies_one_cofactor_per_row() {
    // The same crosstab, but each variable gets its own cofactor, as a second row factor would
    // supply. Checked against the definition rather than a golden, since the R operator ships no
    // fixture for `manual`.
    let cofactors = [1.0, 2.0, 5.0, 150.0, 1000.0];
    let mut values = vec![100.0; 5];
    let ri: Vec<i32> = (0..5).collect();
    asinh_operator::algorithm::asinh_per_row(&mut values, &ri, &cofactors).unwrap();
    for (k, c) in cofactors.iter().enumerate() {
        assert!((values[k] - (100.0f64 / c).asinh()).abs() < 1e-15);
    }
    // and a larger cofactor must compress more
    assert!(values.windows(2).all(|w| w[0] > w[1]));
}

/// `operator.json`'s spec says what a caller will get before the step runs
/// (`DataStep.getPredictedAttributes`). It is hand-written, so compare it with what the code
/// actually produces: the R operator names the output attribute `asinh`, and this port writes
/// `<namespace>.asinh`, which is the same thing once the platform applies the namespace.
#[test]
fn operator_spec_matches_the_result_column() {
    let manifest = concat!(env!("CARGO_MANIFEST_DIR"), "/operator.json");
    let spec: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(manifest).unwrap()).unwrap();
    let attrs = &spec["operatorSpec"]["outputSpecsV2"][0]["joinOperators"][0]["rightRelation"]
        ["attributes"];
    let declared: Vec<&str> = attrs
        .as_array()
        .expect("outputSpecsV2[0].joinOperators[0].rightRelation.attributes")
        .iter()
        .map(|a| a["name"].as_str().unwrap())
        .collect();
    assert_eq!(
        declared,
        ["asinh"],
        "the spec must declare exactly the column the operator writes"
    );
    let ns = "ds0";
    assert_eq!(
        asinh_operator::output::value_column(ns),
        format!("{ns}.{}", declared[0]),
        "the written column and the declared attribute have drifted apart"
    );
    // The input spec has to name a y axis, or the platform cannot tell a user what to project.
    let axis = &spec["operatorSpec"]["inputSpecs"][0]["axis"][0]["metaFactors"][0];
    assert_eq!(axis["crosstabMapping"], "y");
    assert_eq!(axis["cardinality"], "1");
}
