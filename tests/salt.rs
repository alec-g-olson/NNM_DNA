//! Salt correction tests for [`Thermo::tm_with_salt`].
//!
//! Three monovalent-cation correction methods (matching primer3-py's
//! `salt_corrections_method` parameter):
//!
//! - `SantaLucia` (default): SantaLucia, J. (1998). PNAS 95:1460.
//!   ΔS°(salt) = ΔS° + 0.368·(n_bp − 1)·ln([Na⁺]).
//! - `Schildkraut`: Schildkraut & Lifson (1965). *Biopolymers* 3:195.
//!   Tm(salt) = Tm(1 M) + 16.6·log₁₀([Na⁺]).
//! - `Owczarzy`: Owczarzy et al. (2004). *Biochemistry* 43:3537.
//!   1/Tm(salt) = 1/Tm(1 M) + (4.29·fGC − 3.95)·10⁻⁵·ln[Na⁺] + 9.40·10⁻⁶·(ln[Na⁺])².

use nnm_dna::{duplex_thermo, Duplex, SaltCorrection, Sequence};

const R: f64 = 1.987;

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

fn thermo_of(seq: &str) -> nnm_dna::Thermo {
    duplex_thermo(&Duplex::perfect(Sequence::parse(seq).unwrap())).unwrap()
}

#[test]
fn at_1m_all_methods_match_uncorrected_tm() {
    // At [Na⁺] = 1 M, every correction reduces to the uncorrected Tm.
    //   SantaLucia: ln(1) = 0  → ΔS unchanged
    //   Schildkraut: log₁₀(1) = 0  → no offset
    //   Owczarzy: ln(1) = 0  → no inverse-Tm term
    let t = thermo_of("CGTTGA");
    let tm_ref = t.tm(1e-4);
    for method in [
        SaltCorrection::SantaLucia,
        SaltCorrection::Schildkraut,
        SaltCorrection::Owczarzy,
    ] {
        let tm_salt = t.tm_with_salt(1e-4, 1.0, method);
        assert!(
            close(tm_salt, tm_ref, 1e-9),
            "{method:?} at 1 M differs from uncorrected: {tm_salt} vs {tm_ref}"
        );
    }
}

#[test]
fn santalucia_formula_at_50mm() {
    // Hand-calculated: for CGTTGA (n_bp=6) at C_T = 100 μM, 50 mM Na⁺:
    //   ΔS_corrected = ΔS + 0.368·(6-1)·ln(0.050)
    //                = -115.4 + 0.368·5·(-2.9957)
    //                = -115.4 - 5.512
    //                = -120.912 cal/(mol·K)
    //   Tm = -41200 / (-120.912 + 1.987·ln(2.5e-5)) - 273.15
    //      = -41200 / (-141.967) - 273.15
    //      = 290.21 - 273.15
    //      ≈ 17.06 °C
    let t = thermo_of("CGTTGA");
    let tm = t.tm_with_salt(1e-4, 0.050, SaltCorrection::SantaLucia);

    // Direct computation
    let n = (t.n_bp as f64 - 1.0).max(0.0);
    let ds_corrected = t.ds + 0.368 * n * 0.050_f64.ln();
    let dh_cal = t.dh * 1000.0;
    let x = if t.self_comp { 1.0 } else { 4.0 };
    let expected = dh_cal / (ds_corrected + R * (1e-4_f64 / x).ln()) - 273.15;

    assert!(close(tm, expected, 1e-9), "tm = {tm}, expected = {expected}");
    // And confirm the qualitative direction: lower salt → lower Tm
    assert!(tm < t.tm(1e-4), "expected lower Tm at 50 mM");
}

#[test]
fn schildkraut_offset_is_16_6_log10() {
    // For CGTTGA at 50 mM Na⁺:
    //   offset = 16.6 · log₁₀(0.050) = 16.6 · (-1.30103) ≈ -21.597 °C
    let t = thermo_of("CGTTGA");
    let tm_1m = t.tm(1e-4);
    let tm_50mm = t.tm_with_salt(1e-4, 0.050, SaltCorrection::Schildkraut);
    let expected_offset = 16.6 * 0.050_f64.log10();
    assert!(
        close(tm_50mm - tm_1m, expected_offset, 1e-9),
        "offset = {} vs expected {}",
        tm_50mm - tm_1m,
        expected_offset
    );
    assert!(close(expected_offset, -21.597, 0.005));
}

#[test]
fn owczarzy_inverse_tm_formula() {
    // CGTTGA has GC fraction 3/6 = 0.5.
    // At 50 mM Na⁺:
    //   1/Tm_salt - 1/Tm_1M = (4.29·0.333 - 3.95)·1e-5·ln(0.050) + 9.4e-6·ln(0.050)²
    //                       = (1.4286 - 3.95)·1e-5·(-2.9957) + 9.4e-6·8.974
    //                       = (-2.5214)·1e-5·(-2.9957) + 8.4357e-5
    //                       = 7.5527e-5 + 8.4357e-5
    //                       = 1.5988e-4
    let t = thermo_of("CGTTGA");
    let tm_1m_k = t.tm(1e-4) + 273.15;
    let tm_50mm_k = t.tm_with_salt(1e-4, 0.050, SaltCorrection::Owczarzy) + 273.15;

    let f_gc = t.gc_fraction;
    let ln_na = 0.050_f64.ln();
    let expected_delta_inv_tm = (4.29 * f_gc - 3.95) * 1e-5 * ln_na + 9.40e-6 * ln_na * ln_na;
    let actual_delta_inv_tm = 1.0 / tm_50mm_k - 1.0 / tm_1m_k;
    assert!(
        close(actual_delta_inv_tm, expected_delta_inv_tm, 1e-12),
        "Δ(1/Tm) = {actual_delta_inv_tm} vs expected {expected_delta_inv_tm}"
    );
}

#[test]
fn lower_salt_means_lower_tm_for_all_methods() {
    let t = thermo_of("ACGTACGTACGT");
    let tm_1m = t.tm(1e-4);
    for method in [
        SaltCorrection::SantaLucia,
        SaltCorrection::Schildkraut,
        SaltCorrection::Owczarzy,
    ] {
        let tm_50mm = t.tm_with_salt(1e-4, 0.050, method);
        let tm_10mm = t.tm_with_salt(1e-4, 0.010, method);
        assert!(tm_50mm < tm_1m, "{method:?}: 50 mM Tm ({tm_50mm}) not less than 1 M ({tm_1m})");
        assert!(tm_10mm < tm_50mm, "{method:?}: 10 mM Tm ({tm_10mm}) not less than 50 mM ({tm_50mm})");
    }
}

#[test]
fn gc_fraction_populated_correctly() {
    // CGTTGA: C, G, G are G/C → 3/6 = 0.5
    let t = thermo_of("CGTTGA");
    assert!(close(t.gc_fraction, 0.5, 1e-12), "got {}", t.gc_fraction);
    assert_eq!(t.n_bp, 6);

    // GCGCGC: all GC
    let t = thermo_of("GCGCGC");
    assert!(close(t.gc_fraction, 1.0, 1e-12));
    assert_eq!(t.n_bp, 6);

    // AAAA: no GC
    let t = thermo_of("AAAA");
    assert!(close(t.gc_fraction, 0.0, 1e-12));
    assert_eq!(t.n_bp, 4);
}

#[test]
fn default_method_is_santalucia() {
    // Sanity: SaltCorrection::default() == SantaLucia, matching primer3-py's
    // default of "santalucia".
    let default_method: SaltCorrection = Default::default();
    assert_eq!(default_method, SaltCorrection::SantaLucia);
}
