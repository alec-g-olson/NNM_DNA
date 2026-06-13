//! Integration tests against worked examples from SantaLucia 1998 and standard
//! NN-model walkthroughs. Perfect Watson-Crick duplexes only.

use n3dna::{duplex_thermo, Duplex, Sequence};

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

fn perfect_thermo(seq: &str) -> n3dna::Thermo {
    let dx = Duplex::perfect(Sequence::parse(seq).unwrap());
    duplex_thermo(&dx).unwrap()
}

/// 5'-CGTTGA-3' / 3'-GCAACT-5' — non-self-complementary.
///
/// Hand-calc breakdown:
///   NN stacks  CG, GT, TT, TG, GA
///     ΔH:  -10.6 + -8.4 + -7.9 + -8.5 + -8.2  =  -43.6 kcal/mol
///     ΔS:  -27.2 + -22.4 + -22.2 + -22.7 + -22.2 = -116.7 cal/(mol·K)
///   Init: one G·C end (C·G pair) + one A·T end (A·T pair)
///     ΔH:  +0.1 + +2.3 = +2.4
///     ΔS:  -2.8 + +4.1 = +1.3
///   Total:  ΔH = -41.2 kcal/mol, ΔS = -115.4 cal/(mol·K)
///   ΔG°37 = -41.2 - 310.15 * (-115.4/1000) = -5.41 kcal/mol
#[test]
fn cgttga_thermodynamics() {
    let t = perfect_thermo("CGTTGA");
    assert!(!t.self_comp);
    assert!(close(t.dh, -41.2, 0.05), "dh = {}", t.dh);
    assert!(close(t.ds, -115.4, 0.05), "ds = {}", t.ds);
    assert!(close(t.delta_g_37(), -5.41, 0.05), "ΔG37 = {}", t.delta_g_37());
}

#[test]
fn cgttga_tm() {
    let t = perfect_thermo("CGTTGA");
    // At C_T = 100 μM (non-self-comp, x=4):
    //   Tm = -41200 / (-115.4 + 1.987*ln(2.5e-5)) - 273.15  ≈  28.78 °C
    let tm = t.tm(1.0e-4);
    assert!(close(tm, 28.78, 0.1), "Tm = {tm}");
}

/// 5'-CGCGCG-3' — self-complementary.
///
/// NN stacks CG, GC, CG, GC, CG → 3·(-10.6, -27.2) + 2·(-9.8, -24.4)
///     ΔH = -51.4, ΔS = -130.4
/// Init: both ends G·C → 2·(0.1, -2.8) = (0.2, -5.6)
///     ΔH = -51.2, ΔS = -136.0
/// Symmetry correction: (0.0, -1.4)
///     ΔH = -51.2, ΔS = -137.4
/// ΔG°37 = -51.2 - 310.15 * (-137.4/1000) ≈ -8.59 kcal/mol
#[test]
fn cgcgcg_self_complementary() {
    let t = perfect_thermo("CGCGCG");
    assert!(t.self_comp);
    assert!(close(t.dh, -51.2, 0.05), "dh = {}", t.dh);
    assert!(close(t.ds, -137.4, 0.05), "ds = {}", t.ds);
    assert!(close(t.delta_g_37(), -8.59, 0.05), "ΔG37 = {}", t.delta_g_37());
}

#[test]
fn self_complementary_tm_uses_ct_directly() {
    let t = perfect_thermo("CGCGCG");
    // At C_T = 100 μM:
    //   Tm = -51200 / (-137.4 + 1.987*ln(1e-4)) - 273.15
    //      = -51200 / (-137.4 - 18.299) - 273.15
    //      = -51200 / -155.699 - 273.15
    //      ≈ 55.78 °C
    let tm = t.tm(1.0e-4);
    assert!(close(tm, 55.78, 0.1), "Tm = {tm}");
}

#[test]
fn at_only_duplex_uses_at_initiation_only() {
    let t = perfect_thermo("AAAA");
    assert!(!t.self_comp);
    let expected_dh = 3.0 * -7.9 + 2.0 * 2.3;
    let expected_ds = 3.0 * -22.2 + 2.0 * 4.1;
    assert!(close(t.dh, expected_dh, 1e-9), "dh = {} expected {}", t.dh, expected_dh);
    assert!(close(t.ds, expected_ds, 1e-9), "ds = {} expected {}", t.ds, expected_ds);
}
