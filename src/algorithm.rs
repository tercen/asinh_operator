//! The transform. Pure functions over slices: no Tercen types, so it is testable on its own.
//!
//! R: `asinh(.y / cofactor)`. `f64::asinh` and R's `asinh` are both the platform libm, so this
//! reproduces the R operator's numbers to the last bit on the same hardware; the parity test
//! allows 1e-12 relative to leave room for a different libm.

/// `asinh(y / cofactor)` for one value.
#[inline]
pub fn asinh_scaled(y: f64, cofactor: f64) -> f64 {
    (y / cofactor).asinh()
}

/// In place, one cofactor for every value.
pub fn asinh_fixed(values: &mut [f64], cofactor: f64) {
    for v in values.iter_mut() {
        *v = asinh_scaled(*v, cofactor);
    }
}

/// In place, a cofactor per row index (`.ri` → `cofactors[ri]`).
pub fn asinh_per_row(values: &mut [f64], ri: &[i32], cofactors: &[f64]) -> Result<(), usize> {
    for (v, r) in values.iter_mut().zip(ri) {
        let i = usize::try_from(*r).map_err(|_| 0usize)?;
        let c = *cofactors.get(i).ok_or(i)?;
        *v = asinh_scaled(*v, c);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_the_definition() {
        // asinh(x) = ln(x + sqrt(x^2+1)); check a few points against that closed form
        for (y, c) in [(0.0f64, 5.0f64), (1.0, 5.0), (-3.5, 2.0), (1e6, 150.0)] {
            let x: f64 = y / c;
            let want = (x + (x * x + 1.0).sqrt()).ln();
            assert!((asinh_scaled(y, c) - want).abs() <= 1e-12 * want.abs().max(1.0));
        }
    }

    #[test]
    fn zero_maps_to_zero_and_sign_is_kept() {
        assert_eq!(asinh_scaled(0.0, 5.0), 0.0);
        assert!(asinh_scaled(-1.0, 5.0) < 0.0);
    }

    #[test]
    fn per_row_uses_the_row_index() {
        let mut v = vec![10.0, 10.0];
        asinh_per_row(&mut v, &[0, 1], &[5.0, 100.0]).unwrap();
        assert!(v[0] > v[1], "a smaller cofactor must give a larger value");
    }

    #[test]
    fn a_row_index_outside_the_cofactor_table_is_an_error_not_a_panic() {
        let mut v = vec![1.0];
        assert_eq!(asinh_per_row(&mut v, &[7], &[5.0]), Err(7));
    }
}
