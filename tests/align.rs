
//! Integration tests for the Needleman-Wunsch aligner. Each test uses a
//! hand-designed top/bottom strand pair where the expected alignment is
//! unambiguous given the scoring constants in `align.rs`.
//!
//! Reminder on the antiparallel convention: `AlignedDuplex::align(top, bottom)`
//! runs NW on `top` against `bottom.reverse_complement()`, so a duplex with WC
//! pairs at every position has `top == bottom.reverse_complement()`.

use nnm_dna::{AlignError, AlignedDuplex, AlignedPosition, Base, Sequence};

#[test]
fn wc_only_returns_diagonal_no_gaps() {
    // 5'-CGTTGA-3' / 3'-GCAACT-5'  (bottom 5'→3' = "TCAACG")
    let top = Sequence::parse("CGTTGA").unwrap();
    let bot = Sequence::parse("TCAACG").unwrap();
    let dx = AlignedDuplex::align(top, bot).unwrap();
    assert_eq!(dx.columns().len(), 6);
    for (i, col) in dx.columns().iter().enumerate() {
        assert!(col.is_wc_pair(), "column {i} is not WC: {col:?}");
    }
}

#[test]
fn palindrome_aligns_as_diagonal() {
    // CGCGCG is self-complementary
    let top = Sequence::parse("CGCGCG").unwrap();
    let bot = Sequence::parse("CGCGCG").unwrap();
    let dx = AlignedDuplex::align(top, bot).unwrap();
    assert_eq!(dx.columns().len(), 6);
    assert!(dx.is_self_complementary());
    for col in dx.columns() {
        assert!(col.is_wc_pair());
    }
}

#[test]
fn single_top_bulge_recovered() {
    // Designed so that top has one extra base in the middle.
    //   top:    5'-A C G T G-3'   (5 bases; T at pos 3 is bulged)
    //   bottom: 3'-T G C   C-5'   (4 bases; bottom 5'→3' = "CCGT")
    //
    // Pairing should be:
    //   pos 0: A·T WC,  pos 1: C·G WC,  pos 2: G·C WC,
    //   pos 3: BulgeTop(T),
    //   pos 4: G·C WC.
    let top = Sequence::parse("ACGTG").unwrap();
    let bot = Sequence::parse("CCGT").unwrap();
    let dx = AlignedDuplex::align(top, bot).unwrap();
    let cols = dx.columns();
    assert_eq!(cols.len(), 5, "got {} columns", cols.len());

    // The bulged base should be the T from top.
    let bulge_count = cols
        .iter()
        .filter(|c| matches!(c, AlignedPosition::BulgeTop(Base::T)))
        .count();
    assert_eq!(bulge_count, 1, "expected one BulgeTop(T), got cols: {cols:?}");

    // Termini are WC.
    assert!(cols.first().unwrap().is_wc_pair());
    assert!(cols.last().unwrap().is_wc_pair());
}

#[test]
fn single_bottom_bulge_recovered() {
    // Mirror of the previous test — bulge on the bottom strand.
    //   top:    5'-A C   G T-3'   (4 bases)
    //   bottom: 3'-T G C C A-5'   (5 bases; C at the middle is bulged)
    //
    // bottom 5'→3' = "ACCGT". rc("ACCGT") = "ACGGT".
    let top = Sequence::parse("ACGT").unwrap();
    let bot = Sequence::parse("ACCGT").unwrap();
    let dx = AlignedDuplex::align(top, bot).unwrap();
    let cols = dx.columns();
    assert_eq!(cols.len(), 5);

    let bulge_count = cols
        .iter()
        .filter(|c| matches!(c, AlignedPosition::BulgeBottom(_)))
        .count();
    assert_eq!(bulge_count, 1, "expected one BulgeBottom, got: {cols:?}");

    assert!(cols.first().unwrap().is_wc_pair());
    assert!(cols.last().unwrap().is_wc_pair());
}

#[test]
fn multi_base_bulge_recovered() {
    // 3-base bulge on top.
    //   top:    5'-A C G A A A T G-3'   (8 bases; AAA at positions 3,4,5 bulged)
    //   bottom: 3'-T G C       A C-5'   (5 bases; 5'→3' = "CACGT")
    let top = Sequence::parse("ACGAAATG").unwrap();
    let bot = Sequence::parse("CACGT").unwrap();
    let dx = AlignedDuplex::align(top, bot).unwrap();
    let cols = dx.columns();
    assert_eq!(cols.len(), 8);

    let bulge_count = cols
        .iter()
        .filter(|c| matches!(c, AlignedPosition::BulgeTop(_)))
        .count();
    assert_eq!(bulge_count, 3, "expected three consecutive BulgeTop, got: {cols:?}");

    // The bulges should be contiguous, not split.
    let positions: Vec<usize> = cols
        .iter()
        .enumerate()
        .filter(|(_, c)| matches!(c, AlignedPosition::BulgeTop(_)))
        .map(|(i, _)| i)
        .collect();
    assert_eq!(positions.len(), 3);
    assert_eq!(positions[1], positions[0] + 1);
    assert_eq!(positions[2], positions[1] + 1);
}

#[test]
fn asymmetric_internal_loop_via_from_columns() {
    // With our scoring (gap-open = +3.0), the DP prefers "1 bulge + N
    // mismatches" over "2 bulges on opposite strands" because opening a
    // second gap run is expensive. So loops with gaps on BOTH strands aren't
    // produced by the aligner from raw sequences — they're constructed
    // manually (e.g. when a caller has external alignment information) or
    // emerge when the surrounding stems anchor strongly enough to make the
    // double-gap alignment locally optimal in longer molecules.
    //
    // Here we exercise the data path by constructing a 2×1 loop manually
    // and verifying the structure round-trips cleanly.
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
    let dx = AlignedDuplex::from_columns(top, bot, cols).unwrap();
    assert_eq!(dx.columns().len(), 7);
    let top_bulges = dx
        .columns()
        .iter()
        .filter(|c| matches!(c, BulgeTop(_)))
        .count();
    let bot_bulges = dx
        .columns()
        .iter()
        .filter(|c| matches!(c, BulgeBottom(_)))
        .count();
    assert_eq!(top_bulges, 2);
    assert_eq!(bot_bulges, 1);
}

#[test]
fn align_rejects_empty_inputs() {
    // Sequence::parse won't allow empty/short input, so we can't easily
    // hit EmptyStrand via the public API. Instead, exercise the constructor
    // with explicit empty columns.
    let top = Sequence::parse("AA").unwrap();
    let bot = Sequence::parse("TT").unwrap();
    assert_eq!(
        AlignedDuplex::from_columns(top, bot, vec![]),
        Err(AlignError::EmptyStrand)
    );
}
