//! Length-dependent loop penalties for DNA bulges and internal loops.
//!
//! Source: SantaLucia & Hicks (2004), *Annu. Rev. Biophys. Biomol. Struct.*
//! 33:415-440, Table 4 (bulge column). Convention: ΔH° = 0 for loop terms; ΔS°
//! is back-derived from ΔG°₃₇ via ΔG = ΔH − T·ΔS, i.e. ΔS° = −ΔG°₃₇·1000/T_REF.
//!
//! Lengths outside the tabulated set use the Jacobson-Stockmayer entropy
//! extrapolation (Equation 7 in the same paper):
//!
//! ```text
//!   ΔG°₃₇(loop-n) = ΔG°₃₇(loop-x) + 2.44 · R · T_REF · ln(n/x)
//! ```
//!
//! where `x` is the largest tabulated length ≤ `n`.

use crate::thermo::ThermoError;

/// SantaLucia & Hicks 2004 Table 4 — bulge ΔG°₃₇ (kcal/mol) at 1 M NaCl.
const BULGE_G37: &[(usize, f64)] = &[
    (1, 4.0),
    (2, 2.9),
    (3, 3.1),
    (4, 3.2),
    (5, 3.3),
    (6, 3.5),
    (7, 3.7),
    (8, 3.9),
    (9, 4.1),
    (10, 4.3),
    (12, 4.5),
    (14, 4.8),
    (16, 5.0),
    (18, 5.2),
    (20, 5.3),
    (25, 5.6),
    (30, 5.9),
];

const R: f64 = 1.987;
const T_REF: f64 = 310.15;

/// Returns `(ΔH°, ΔS°)` for a bulge of `n` bases, in `(kcal/mol, cal/(mol·K))`.
///
/// `ΔH°` is zero by the SantaLucia & Hicks 2004 convention. `ΔS°` is derived
/// from the tabulated `ΔG°₃₇` so that `ΔG = ΔH − T·ΔS` evaluates back to the
/// table value at `T = 310.15 K`.
///
/// For lengths between tabulated entries (e.g. n = 11, 13, 15, ...) and for
/// `n > 30`, the Jacobson-Stockmayer extrapolation is applied from the largest
/// tabulated length ≤ `n`.
pub fn bulge_length(n: usize) -> Result<(f64, f64), ThermoError> {
    if n == 0 {
        return Err(ThermoError::UnknownLoopLength(0));
    }
    let g37 = lookup_or_extrapolate(n);
    let dh = 0.0;
    let ds = -g37 * 1000.0 / T_REF;
    Ok((dh, ds))
}

/// v1 fallback: internal loops use the bulge length column until the proper
/// internal-loop column from SantaLucia & Hicks 2004 Table 4 is extracted.
/// Underestimates loop destabilization by ~0.5 kcal/mol on average.
pub fn internal_loop_length(n: usize) -> Result<(f64, f64), ThermoError> {
    bulge_length(n)
}

fn lookup_or_extrapolate(n: usize) -> f64 {
    if let Some(&(_, g)) = BULGE_G37.iter().find(|&&(k, _)| k == n) {
        return g;
    }
    let (x, g_x) = BULGE_G37
        .iter()
        .rev()
        .find(|&&(k, _)| k <= n)
        .copied()
        .unwrap_or((1, 4.0));
    // 2.44 · R · T_REF · ln(n/x) is in cal/mol; scale to kcal/mol.
    g_x + 2.44 * R * T_REF * ((n as f64 / x as f64).ln()) / 1000.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn zero_length_is_error() {
        assert!(matches!(bulge_length(0), Err(ThermoError::UnknownLoopLength(0))));
    }

    #[test]
    fn length_1_is_4_kcal() {
        let (dh, ds) = bulge_length(1).unwrap();
        assert_eq!(dh, 0.0);
        // ΔS = -4.0 * 1000 / 310.15 = -12.897 cal/(mol·K)
        assert!(close(ds, -12.897, 0.01), "ds = {ds}");
        // Round-trip: ΔG37 = ΔH - T·ΔS/1000 = 0 - 310.15 · (-12.897)/1000 ≈ 4.0
        let g37 = dh - T_REF * ds / 1000.0;
        assert!(close(g37, 4.0, 1e-6), "ΔG37 = {g37}");
    }

    #[test]
    fn length_2_uses_table_value() {
        let (_, ds) = bulge_length(2).unwrap();
        // ΔG37(2) = 2.9 → ΔS = -2.9·1000/310.15 ≈ -9.350
        assert!(close(ds, -9.350, 0.01), "ds = {ds}");
    }

    #[test]
    fn nontabulated_length_extrapolates_from_below() {
        // n = 11 not tabulated; nearest x ≤ 11 is 10 with g_x = 4.3.
        // ΔG37(11) = 4.3 + 2.44·1.987·310.15·ln(11/10)/1000
        //         = 4.3 + 1503.04 / 1000 · 0.09531
        //         ≈ 4.3 + 0.1432
        //         ≈ 4.443
        let (_, ds) = bulge_length(11).unwrap();
        let g37 = -T_REF * ds / 1000.0;
        let expected = 4.3 + 2.44 * R * T_REF * (11.0_f64 / 10.0).ln() / 1000.0;
        assert!(close(g37, expected, 1e-9), "ΔG37(11) = {g37} expected {expected}");
    }

    #[test]
    fn long_loop_extrapolates_past_30() {
        // n = 50: extrapolate from x = 30, g_x = 5.9.
        let (_, ds) = bulge_length(50).unwrap();
        let g37 = -T_REF * ds / 1000.0;
        let expected = 5.9 + 2.44 * R * T_REF * (50.0_f64 / 30.0).ln() / 1000.0;
        assert!(close(g37, expected, 1e-9), "ΔG37(50) = {g37} expected {expected}");
        // Sanity: should be larger than 5.9 (longer = more destabilizing).
        assert!(g37 > 5.9);
    }

    #[test]
    fn monotone_nondecreasing_above_length_2() {
        // The literature table is non-monotonic at length 1 (4.0 then 2.9 at 2),
        // but from length 2 onward it increases. Verify our lookup preserves this.
        let mut last = f64::NEG_INFINITY;
        for n in 2..=30 {
            let (_, ds) = bulge_length(n).unwrap();
            let g37 = -T_REF * ds / 1000.0;
            assert!(g37 >= last - 1e-9, "non-monotone at n={n}: {g37} < {last}");
            last = g37;
        }
    }

    #[test]
    fn internal_loop_matches_bulge_fallback() {
        // v1: same table. Once we extract the proper internal-loop column,
        // this test will need to be updated.
        for n in [1, 3, 5, 10] {
            assert_eq!(bulge_length(n).unwrap(), internal_loop_length(n).unwrap());
        }
    }
}
