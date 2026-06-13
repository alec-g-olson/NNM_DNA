//! End-to-end thermo tests for duplexes containing internal loops.
//!
//! Internal loops can't be inferred from raw sequence by our heuristic NW
//! aligner (with the current scoring, "1 bulge + N mismatches" always beats
//! "two opposite-strand bulges" cost-wise — see `tests/align.rs` notes), so
//! these tests construct `AlignedDuplex` via `from_columns` directly.
//!
//! v1: internal-loop length penalty uses the bulge column from SantaLucia &
//! Hicks 2004 Table 4 as a fallback. Will be replaced when we extract the
//! proper internal-loop column.

use n3dna::{aligned_duplex_thermo, AlignedDuplex, AlignedPosition, Base, Sequence};

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

/// Asymmetric 2×1 loop (total length 3) embedded in a 6 bp duplex.
///
/// 5'-A C A A T G-3'   (6 bases on top)
/// 3'-T G C   A C-5'   (5 bases on bottom; bot 5'→3' = "CACGT")
///
/// Manual alignment:
///   col 0: Pair(A, T)   WC
///   col 1: Pair(C, G)   WC
///   col 2: BulgeTop(A)
///   col 3: BulgeTop(A)
///   col 4: BulgeBottom(C)
///   col 5: Pair(T, A)   WC
///   col 6: Pair(G, C)   WC
///
/// Features:
///   WcStem [A, C], InternalLoop (L=3), WcStem [T, G]
///
/// Init: A (AT: +2.3, +4.1) + G (GC: +0.1, -2.8) = (+2.4, +1.3)
/// WcStem [A,C] NN: AC = (-8.4, -22.4)
/// Loop L=3: ΔH=0, ΔS = -3.1·1000/310.15 = -9.995
/// WcStem [T,G] NN: TG = canonical CA/GT = (-8.5, -22.7)
///
/// ΔH = +2.4 - 8.4 + 0 - 8.5 = -14.5 kcal/mol
/// ΔS = +1.3 - 22.4 - 9.995 - 22.7 = -53.795 cal/(mol·K)
#[test]
fn asymmetric_loop_2_top_1_bottom() {
    use AlignedPosition::*;
    use Base::*;
    let top = Sequence::parse("ACAATG").unwrap();
    let bot = Sequence::parse("CACGT").unwrap();
    let cols = vec![
        Pair { top: A, bottom: T },
        Pair { top: C, bottom: G },
        BulgeTop(A),
        BulgeTop(A),
        BulgeBottom(C),
        Pair { top: T, bottom: A },
        Pair { top: G, bottom: C },
    ];
    let aligned = AlignedDuplex::from_columns(top, bot, cols).unwrap();
    let t = aligned_duplex_thermo(&aligned).unwrap();
    assert!(close(t.dh, -14.5, 0.05), "dh = {}", t.dh);
    assert!(close(t.ds, -53.795, 0.05), "ds = {}", t.ds);
}

/// Symmetric 2×2 loop (tandem mismatch, total length 4) in a 4 bp duplex.
///
/// 5'-A C G T-3'
/// 3'-T A T A-5'   (bot 5'→3' = "ATAT")
///
/// Pairings:  pos 0: A·T WC,  pos 1: C·A mismatch,  pos 2: G·T mismatch,  pos 3: T·A WC
///
/// Manual columns:
///   col 0: Pair(A, T)  WC
///   col 1: Pair(C, A)  non-WC (C·A mismatch)
///   col 2: Pair(G, T)  non-WC (G·T mismatch)
///   col 3: Pair(T, A)  WC
///
/// Two consecutive non-WC Pairs → classified as InternalLoop (length 4).
///
/// Init: A (AT) + T (AT) = 2·(+2.3, +4.1) = (+4.6, +8.2)
/// WcStem [A]: no stacks
/// Loop L=4: ΔH=0, ΔS = -3.2·1000/310.15 = -10.318
/// WcStem [T]: no stacks
///
/// ΔH = +4.6 + 0 = +4.6 kcal/mol  (positive — unstable)
/// ΔS = +8.2 - 10.318 = -2.118 cal/(mol·K)
#[test]
fn symmetric_2x2_loop_tandem_mismatch() {
    use AlignedPosition::*;
    use Base::*;
    let top = Sequence::parse("ACGT").unwrap();
    let bot = Sequence::parse("ATAT").unwrap();
    let cols = vec![
        Pair { top: A, bottom: T },
        Pair { top: C, bottom: A },
        Pair { top: G, bottom: T },
        Pair { top: T, bottom: A },
    ];
    let aligned = AlignedDuplex::from_columns(top, bot, cols).unwrap();
    let t = aligned_duplex_thermo(&aligned).unwrap();
    assert!(close(t.dh, 4.6, 0.05), "dh = {}", t.dh);
    assert!(close(t.ds, -2.118, 0.05), "ds = {}", t.ds);
}

/// Sanity: an asymmetric loop is destabilizing relative to a perfect duplex
/// of the same paired-region length. Compares 2×1-loop duplex against a
/// 4 bp WC duplex of the same flanking stems.
#[test]
fn loop_destabilizes_relative_to_perfect_stems() {
    use AlignedPosition::*;
    use Base::*;

    // 4 bp perfect duplex with same flanking stems: ACTG/CAGT
    let perfect = n3dna::Duplex::perfect(Sequence::parse("ACTG").unwrap());
    let t_perfect = n3dna::duplex_thermo(&perfect).unwrap();

    // Same flanking stems + 2×1 loop in middle
    let top = Sequence::parse("ACAATG").unwrap();
    let bot = Sequence::parse("CACGT").unwrap();
    let cols = vec![
        Pair { top: A, bottom: T },
        Pair { top: C, bottom: G },
        BulgeTop(A),
        BulgeTop(A),
        BulgeBottom(C),
        Pair { top: T, bottom: A },
        Pair { top: G, bottom: C },
    ];
    let aligned = AlignedDuplex::from_columns(top, bot, cols).unwrap();
    let t_loop = aligned_duplex_thermo(&aligned).unwrap();

    assert!(
        t_loop.delta_g_37() > t_perfect.delta_g_37(),
        "loop ΔG ({}) should be less negative than perfect ΔG ({})",
        t_loop.delta_g_37(),
        t_perfect.delta_g_37(),
    );
}
