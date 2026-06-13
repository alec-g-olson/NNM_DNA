//! Integration tests for duplexes containing internal single-base mismatches.
//!
//! All 8 mismatch types are supported. Parameter sources:
//! - G·T: Allawi & SantaLucia (1997), Biochemistry 36:10581
//! - G·A: Allawi & SantaLucia (1998), Biochemistry 37:2170
//! - C·T: Allawi & SantaLucia (1998), NAR 26:2694
//! - A·C: Allawi & SantaLucia (1998), Biochemistry 37:9435 (pH 7.0)
//! - A·A, C·C, G·G, T·T: Peyret et al. (1999), Biochemistry 38:3468

use n3dna::{duplex_thermo, Base, Duplex, DuplexError, Sequence};

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

/// Worked example — single internal G·T mismatch.
///
/// Duplex:
///   5'-CGTGC-3'
///   3'-GTACG-5'   (bottom written 5'→3' = "GCATG")
///
/// Pairing:
///   pos 0: C·G WC
///   pos 1: G·T mismatch
///   pos 2: T·A WC
///   pos 3: G·C WC
///   pos 4: C·G WC
///
/// NN stacks:
///   0–1: 5'-CG-3'/3'-GT-5'  →  G·T-mm table "CG/GT" =  (-4.1, -11.7)
///   1–2: 5'-GT-3'/3'-TA-5'  →  G·T-mm, strand-flip of "AT/TG" =  (-2.5, -8.3)
///   2–3: 5'-TG-3'/3'-AC-5'  →  WC, nn_params(T,G) = (-8.5, -22.7)
///   3–4: 5'-GC-3'/3'-CG-5'  →  WC, nn_params(G,C) = (-9.8, -24.4)
///
///   ΔH (stacks) = -4.1 + -2.5 + -8.5 + -9.8  = -24.9
///   ΔS (stacks) = -11.7 + -8.3 + -22.7 + -24.4 = -67.1
///
/// Init: both ends C·G → 2·(0.1, -2.8) = (0.2, -5.6)
///
/// Total:
///   ΔH = -24.7 kcal/mol
///   ΔS = -72.7 cal/(mol·K)
///   ΔG°37 = -24.7 - 310.15·(-72.7/1000) = -2.15 kcal/mol
///
/// Compare against the same top strand paired with its true complement:
/// perfect "CGTGC/GCACG" gives ΔG°37 ≈ -5.37 — a ~3.2 kcal/mol penalty
/// for the single mismatch, which is in the expected range.
#[test]
fn single_gt_mismatch_cgtgc() {
    let dx = Duplex::new(
        Sequence::parse("CGTGC").unwrap(),
        Sequence::parse("GCATG").unwrap(),
    )
    .unwrap();
    assert!(!dx.is_perfect());
    assert_eq!(dx.pair_at(1), (Base::G, Base::T));

    let t = duplex_thermo(&dx).unwrap();
    assert!(close(t.dh, -24.7, 0.05), "dh = {}", t.dh);
    assert!(close(t.ds, -72.7, 0.05), "ds = {}", t.ds);
    assert!(close(t.delta_g_37(), -2.15, 0.05), "ΔG37 = {}", t.delta_g_37());
    assert!(!t.self_comp);
}

#[test]
fn mismatch_destabilizes_relative_to_perfect() {
    // Same top strand, two different bottoms: WC vs single-mismatch.
    let top = Sequence::parse("CGTGC").unwrap();

    let perfect = Duplex::perfect(top.clone());
    let mismatched =
        Duplex::new(top, Sequence::parse("GCATG").unwrap()).unwrap();

    let t_perfect = duplex_thermo(&perfect).unwrap();
    let t_mismatched = duplex_thermo(&mismatched).unwrap();

    // The mismatch must be destabilizing: ΔG°37 less negative.
    assert!(t_mismatched.delta_g_37() > t_perfect.delta_g_37());
    let penalty = t_mismatched.delta_g_37() - t_perfect.delta_g_37();
    assert!(penalty > 2.0 && penalty < 5.0, "ΔΔG penalty = {penalty}");
}

/// Hand-calculated A·A mismatch — previously failed with UnsupportedMismatch,
/// now resolves via Peyret 1999 parameters.
///
/// Duplex:
///   5'-A C A T-3'
///   3'-T G A A-5'   (bottom 5'→3' = "AAGT")
///
/// Pairings: A·T, C·G, A·A (mm), T·A.
///
/// NN stacks:
///   0-1: 5'-AC-3'/3'-TG-5'   WC, nn_params(A,C) = canonical GT/CA = (-8.4, -22.4)
///   1-2: 5'-CA-3'/3'-GA-5'   mm at pos 2: A·A. "CA/GA" → (-0.9, -4.2)  (Peyret 1999)
///   2-3: 5'-AT-3'/3'-AA-5'   mm at pos 1: A·A. Strand-flip of "AT/AA" → "AA/TA" → (1.2, 1.7)
///
/// ΔH stacks = -8.4 + -0.9 + 1.2  = -8.1
/// ΔS stacks = -22.4 + -4.2 + 1.7 = -24.9
///
/// Init: A·T (both ends, AT init): 2·(+2.3, +4.1) = (+4.6, +8.2)
/// Total ΔH = -8.1 + 4.6 = -3.5
/// Total ΔS = -24.9 + 8.2 = -16.7
#[test]
fn aa_mismatch_acat_now_supported() {
    let dx = Duplex::new(
        Sequence::parse("ACAT").unwrap(),
        Sequence::parse("AAGT").unwrap(),
    )
    .unwrap();
    let t = duplex_thermo(&dx).unwrap();
    assert!(close(t.dh, -3.5, 0.05), "dh = {}", t.dh);
    assert!(close(t.ds, -16.7, 0.05), "ds = {}", t.ds);
}

/// Hand-calculated C·T mismatch.
///
/// Duplex:
///   5'-A C T G-3'
///   3'-T A G C-5'   (bottom 5'→3' = "CGAT")
///
/// Pairings: A·T, C·A (mm — wait that's A·C), T·G WC? No.
///
/// Let me redesign:
///   5'-A C T G-3'
///   3'-T G G C-5'   (bottom 5'→3' = "CGGT")
///
/// Pairings: A·T, C·G, T·G (mm — that's G·T not C·T), G·C. No.
///
/// To get a clean C·T mismatch (C on top, T on bottom):
///   5'-A C T-3'                (3-mer)
///   3'-T T A-5'   (bot 5'→3' = "ATT")
///
/// Pairings:
///   pos 0: A·T WC
///   pos 1: C·T mismatch (C on top, T on bottom)
///   pos 2: T·A WC
///
/// NN stacks:
///   0-1: 5'-AC-3'/3'-TT-5'  → "AC/TT" = (0.7, 0.2)  (C·T table, mm at pos 2)
///   1-2: 5'-CT-3'/3'-TA-5'  → strand-flip of "AT/TC" = (-1.2, -6.2)
///
/// ΔH stacks = 0.7 + -1.2 = -0.5
/// ΔS stacks = 0.2 + -6.2 = -6.0
/// Init: A·T both ends = 2·(+2.3, +4.1) = (+4.6, +8.2)
/// Total ΔH = -0.5 + 4.6 = 4.1
/// Total ΔS = -6.0 + 8.2 = 2.2
#[test]
fn ct_mismatch_now_supported() {
    let dx = Duplex::new(
        Sequence::parse("ACT").unwrap(),
        Sequence::parse("ATT").unwrap(),
    )
    .unwrap();
    let t = duplex_thermo(&dx).unwrap();
    assert_eq!(dx.pair_at(1), (Base::C, Base::T));
    assert!(close(t.dh, 4.1, 0.05), "dh = {}", t.dh);
    assert!(close(t.ds, 2.2, 0.05), "ds = {}", t.ds);
}

#[test]
fn terminal_mismatch_rejected_at_construction() {
    // Top A at pos 0 should pair with T at bottom-end, but bottom 5'→3' is
    // "AAAA" → bottom[3] = A, so pos 0 is A·A mismatch at the terminus.
    let err = Duplex::new(
        Sequence::parse("AAAA").unwrap(),
        Sequence::parse("AAAA").unwrap(),
    )
    .unwrap_err();
    assert_eq!(err, DuplexError::TerminalMismatch);
}

#[test]
fn length_mismatch_rejected() {
    let err = Duplex::new(
        Sequence::parse("ACGT").unwrap(),
        Sequence::parse("ACG").unwrap(),
    )
    .unwrap_err();
    assert_eq!(err, DuplexError::LengthMismatch { top: 4, bottom: 3 });
}
