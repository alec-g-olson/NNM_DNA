//! End-to-end thermo tests for duplexes containing internal bulges.
//!
//! Worked-example values are hand-calculated from:
//! - SantaLucia 1998 unified NN ΔH°/ΔS° (`src/params.rs`)
//! - SantaLucia & Hicks 2004 Table 4 bulge length penalties (`src/params_loop.rs`)
//! - Equation 13: bulge ΔG = length-penalty + intervening-NN (length-1 only)
//!   + 0.5·(# AT flanks) (length-1 only)
//!
//! All values at 1 M NaCl, 37 °C.

use nnm_dna::{aligned_duplex_thermo, AlignedDuplex, Duplex, Sequence};

const T_REF: f64 = 310.15;

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

fn thermo_aligned(top_str: &str, bot_str: &str) -> nnm_dna::Thermo {
    let top = Sequence::parse(top_str).unwrap();
    let bot = Sequence::parse(bot_str).unwrap();
    let aligned = AlignedDuplex::align(top, bot).unwrap();
    aligned_duplex_thermo(&aligned).unwrap()
}

/// 5'-A C G T G-3' / 3'-T G C   C-5'  (bottom 5'→3' = "CCGT")
///
/// Single T bulge on top at position 3. Flanking pairs: G·C at left close
/// (top G at index 2), G·C at right close (top G at index 4). Both flanks
/// are G·C → no AT bonus.
///
/// Stacks: WcStem [A, C, G] → AC (-8.4, -22.4) + CG (-10.6, -27.2)
/// Bulge length-1 penalty: ΔH=0, ΔS = -4.0·1000/310.15 = -12.897
/// Intervening NN: nn_params(G, G) → GG/CC = (-8.0, -19.9)
/// No AT bonus.
/// WcStem [G] alone: no stacks.
///
/// Init: A·T (5' terminus) + G·C (3' terminus) = (+2.3, +4.1) + (+0.1, -2.8)
///                                              = (+2.4, +1.3)
///
/// Total ΔH = +2.4 + (-8.4) + (-10.6) + 0.0 + (-8.0) = -24.6 kcal/mol
/// Total ΔS = +1.3 + (-22.4) + (-27.2) + (-12.897) + (-19.9) = -81.097 cal/(mol·K)
/// ΔG°₃₇   = -24.6 - 310.15·(-81.097)/1000 = +0.55 kcal/mol  (barely unstable)
#[test]
fn single_top_bulge_with_gc_flanks() {
    let t = thermo_aligned("ACGTG", "CCGT");
    assert!(close(t.dh, -24.6, 0.05), "dh = {}", t.dh);
    assert!(close(t.ds, -81.097, 0.05), "ds = {}", t.ds);
    assert!(close(t.delta_g_37(), 0.55, 0.05), "ΔG37 = {}", t.delta_g_37());
}

/// Single G bulge between two A·T pairs — both flanks are AT, so the AT
/// bonus fires twice (+1.0 kcal/mol total).
///
/// 5'-A G T-3' / 3'-T   A-5'  (bottom 5'→3' = "AT")
///
/// WcStem [A]: no stacks. WcStem [T]: no stacks.
/// Bulge length-1: ΔH=0, ΔS=-12.897.
/// Intervening NN: nn_params(A, T) → AT/AT = (-7.2, -20.4)
/// AT bonus: 2·0.5 = 1.0 kcal/mol ΔG37 → ΔS contribution = -1.0·1000/310.15 = -3.224
/// Init: 2× A·T = 2·(+2.3, +4.1) = (+4.6, +8.2)
///
/// Total ΔH = +4.6 + 0.0 + (-7.2) = -2.6
/// Total ΔS = +8.2 + (-12.897) + (-20.4) + (-3.224) = -28.321
#[test]
fn single_top_bulge_with_at_flanks() {
    let t = thermo_aligned("AGT", "AT");
    assert!(close(t.dh, -2.6, 0.05), "dh = {}", t.dh);
    assert!(close(t.ds, -28.321, 0.05), "ds = {}", t.ds);
    // ΔG37 = -2.6 - 310.15·(-28.321)/1000 = +6.18 kcal/mol (very unstable, only 2 bp)
    assert!(close(t.delta_g_37(), 6.18, 0.05), "ΔG37 = {}", t.delta_g_37());
}

/// Three-base bulge (AAA on top) — length ≥ 2 means **no** intervening NN
/// and **no** AT bonus per Eq. 13. Just the length penalty.
///
/// 5'-A C G A A A T G-3' / 3'-T G C       A C-5'  (bottom 5'→3' = "CACGT")
///
/// WcStem [A, C, G]: AC + CG = (-19.0, -49.6)
/// Bulge length-3: ΔH=0, ΔS = -3.1·1000/310.15 = -9.995
/// WcStem [T, G]: TG (= CA/GT canonical) = (-8.5, -22.7)
/// Init: A·T + G·C = (+2.4, +1.3)
///
/// Total ΔH = +2.4 - 19.0 + 0.0 - 8.5 = -25.1
/// Total ΔS = +1.3 - 49.6 - 9.995 - 22.7 = -80.995
#[test]
fn three_base_bulge_no_intervening_nn() {
    let t = thermo_aligned("ACGAAATG", "CACGT");
    assert!(close(t.dh, -25.1, 0.05), "dh = {}", t.dh);
    assert!(close(t.ds, -80.995, 0.05), "ds = {}", t.ds);
}

/// Compare a duplex containing a single bulge against the same top strand
/// paired with its true complement. The bulge penalty is the +4.0 kcal/mol
/// Table 4 length-1 value *plus* the loss of one paired NN stack — since the
/// bulged version effectively shortens the duplex by one base pair. Expected
/// total destabilization is in the 4-7 kcal/mol range depending on the
/// specific stacks involved.
#[test]
fn bulge_destabilizes_relative_to_perfect() {
    let top = Sequence::parse("ACGTG").unwrap();

    // Perfect: bottom = rev_comp(top) = "CACGT"
    let perfect = Duplex::perfect(top.clone());
    let t_perfect = nnm_dna::duplex_thermo(&perfect).unwrap();

    // Bulged: replace bottom with a 4-mer so top[3]=T has no partner
    let bot = Sequence::parse("CCGT").unwrap();
    let aligned = AlignedDuplex::align(top, bot).unwrap();
    let t_bulged = aligned_duplex_thermo(&aligned).unwrap();

    let penalty = t_bulged.delta_g_37() - t_perfect.delta_g_37();
    assert!(
        penalty > 4.0 && penalty < 7.0,
        "ΔΔG penalty = {penalty} kcal/mol — expected in 4-7 range"
    );
}

/// Verify the ΔS-from-ΔG derivation is internally consistent: ΔG°₃₇
/// computed via ΔH − T·ΔS exactly equals the table-derived value at T_REF.
#[test]
fn bulge_dg_round_trip() {
    let t = thermo_aligned("ACGTG", "CCGT");
    let dg_via_struct = t.delta_g_37();
    let dg_via_formula = t.dh - T_REF * t.ds / 1000.0;
    assert!(close(dg_via_struct, dg_via_formula, 1e-9));
}
