//! Structural feature extraction from an `AlignedDuplex`.
//!
//! Walks the alignment column-by-column and groups consecutive columns into
//! features:
//!
//! - `WcStem` — a run of one or more consecutive Watson-Crick pairs. NN
//!   stacks are summed inside the run; stack boundaries between adjacent
//!   `WcStem` features are *not* crossed (something — a mismatch, bulge, or
//!   loop — sits between them).
//! - `Mismatch` — a single mismatch column flanked by WC pairs.
//! - `Bulge` — a run of gap-on-one-strand columns (no opposing base at all).
//! - `InternalLoop` — a run where both strands have at least one unpaired
//!   base, possibly with non-WC `Pair` columns interleaved.
//!
//! Termini are validated as WC pairs; non-WC terminal columns produce
//! `AlignError::TerminalMismatch` / `TerminalBulge`.
//!
//! Tandem mismatches (two consecutive non-WC `Pair` columns with no gaps)
//! are classified as `InternalLoop { top_bases: [a, b], bottom_bases: [c, d] }`
//! per the SantaLucia & Hicks 2004 treatment of length-2 internal loops.

use crate::align::{AlignError, AlignedDuplex, AlignedPosition};
use crate::base::Base;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Feature {
    /// Run of ≥ 1 consecutive WC pairs. `start_idx` is the column index of
    /// the first pair in the run.
    WcStem {
        tops: Vec<Base>,
        bottoms: Vec<Base>,
        start_idx: usize,
    },
    /// A single non-WC `Pair` column flanked by WC pairs on both sides.
    Mismatch {
        top: Base,
        bottom: Base,
        idx: usize,
    },
    /// A bulge: a run of gap-on-one-strand columns. `bases` reads 5'→3'
    /// along whichever strand carries the unpaired bases. `left_close` and
    /// `right_close` are the column indices of the flanking WC pairs.
    Bulge {
        on_top: bool,
        bases: Vec<Base>,
        left_close: usize,
        right_close: usize,
    },
    /// An internal loop: a run with unpaired bases on both strands.
    InternalLoop {
        top_bases: Vec<Base>,
        bottom_bases: Vec<Base>,
        left_close: usize,
        right_close: usize,
    },
}

pub fn extract_features(aligned: &AlignedDuplex) -> Result<Vec<Feature>, AlignError> {
    let cols = aligned.columns();
    if cols.len() < 2 {
        return Err(AlignError::EmptyStrand);
    }
    // Termini must already be WC by `AlignedDuplex` construction, but check
    // defensively since `from_columns` could be bypassed by future callers.
    if !cols.first().unwrap().is_wc_pair() {
        return match cols.first().unwrap() {
            AlignedPosition::Pair { .. } => Err(AlignError::TerminalMismatch),
            _ => Err(AlignError::TerminalBulge),
        };
    }
    if !cols.last().unwrap().is_wc_pair() {
        return match cols.last().unwrap() {
            AlignedPosition::Pair { .. } => Err(AlignError::TerminalMismatch),
            _ => Err(AlignError::TerminalBulge),
        };
    }

    let mut features = Vec::new();
    let mut i = 0;
    while i < cols.len() {
        if cols[i].is_wc_pair() {
            let stem_start = i;
            let mut tops = Vec::new();
            let mut bottoms = Vec::new();
            while i < cols.len() && cols[i].is_wc_pair() {
                if let AlignedPosition::Pair { top, bottom } = cols[i] {
                    tops.push(top);
                    bottoms.push(bottom);
                }
                i += 1;
            }
            features.push(Feature::WcStem {
                tops,
                bottoms,
                start_idx: stem_start,
            });
        } else {
            let run_start = i;
            let mut top_bases = Vec::new();
            let mut bottom_bases = Vec::new();
            let mut has_top_gap = false;
            let mut has_bot_gap = false;
            let mut has_non_wc_pair = false;
            while i < cols.len() && !cols[i].is_wc_pair() {
                match cols[i] {
                    AlignedPosition::Pair { top, bottom } => {
                        has_non_wc_pair = true;
                        top_bases.push(top);
                        bottom_bases.push(bottom);
                    }
                    AlignedPosition::BulgeTop(b) => {
                        has_top_gap = true;
                        top_bases.push(b);
                    }
                    AlignedPosition::BulgeBottom(b) => {
                        has_bot_gap = true;
                        bottom_bases.push(b);
                    }
                }
                i += 1;
            }
            let run_len = i - run_start;
            let left_close = run_start - 1; // termini are WC → run_start ≥ 1
            let right_close = i; // termini are WC → i < cols.len() at this point

            if run_len == 1 && has_non_wc_pair && !has_top_gap && !has_bot_gap {
                features.push(Feature::Mismatch {
                    top: top_bases[0],
                    bottom: bottom_bases[0],
                    idx: run_start,
                });
            } else if has_top_gap && !has_bot_gap && !has_non_wc_pair {
                features.push(Feature::Bulge {
                    on_top: true,
                    bases: top_bases,
                    left_close,
                    right_close,
                });
            } else if has_bot_gap && !has_top_gap && !has_non_wc_pair {
                features.push(Feature::Bulge {
                    on_top: false,
                    bases: bottom_bases,
                    left_close,
                    right_close,
                });
            } else {
                features.push(Feature::InternalLoop {
                    top_bases,
                    bottom_bases,
                    left_close,
                    right_close,
                });
            }
        }
    }

    Ok(features)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::Sequence;
    use Base::*;

    fn wc(top: Base) -> AlignedPosition {
        AlignedPosition::Pair {
            top,
            bottom: top.complement(),
        }
    }

    #[test]
    fn all_wc_emits_one_stem() {
        let cols = vec![wc(A), wc(C), wc(G), wc(T)];
        let dx = AlignedDuplex::from_columns(
            Sequence::parse("ACGT").unwrap(),
            Sequence::parse("ACGT").unwrap(),
            cols,
        )
        .unwrap();
        let features = extract_features(&dx).unwrap();
        assert_eq!(features.len(), 1);
        match &features[0] {
            Feature::WcStem { tops, bottoms, start_idx } => {
                assert_eq!(tops, &vec![A, C, G, T]);
                assert_eq!(bottoms, &vec![T, G, C, A]);
                assert_eq!(*start_idx, 0);
            }
            _ => panic!("expected WcStem, got {:?}", features[0]),
        }
    }

    #[test]
    fn bulge_splits_stems() {
        // Stem(2) | Bulge(top, 1) | Stem(2)
        let cols = vec![
            wc(A),
            wc(C),
            AlignedPosition::BulgeTop(G),
            wc(T),
            wc(A),
        ];
        let dx = AlignedDuplex::from_columns(
            Sequence::parse("ACGTA").unwrap(),
            Sequence::parse("TACGT").unwrap(),
            cols,
        )
        .unwrap();
        let features = extract_features(&dx).unwrap();
        assert_eq!(features.len(), 3);
        assert!(matches!(features[0], Feature::WcStem { .. }));
        match &features[1] {
            Feature::Bulge { on_top, bases, left_close, right_close } => {
                assert!(*on_top);
                assert_eq!(bases, &vec![G]);
                assert_eq!(*left_close, 1);
                assert_eq!(*right_close, 3);
            }
            _ => panic!("expected Bulge, got {:?}", features[1]),
        }
        assert!(matches!(features[2], Feature::WcStem { .. }));
    }

    #[test]
    fn isolated_mismatch_is_mismatch_feature() {
        // Stem(2) | Pair(G·T mismatch) | Stem(2)
        let cols = vec![
            wc(A),
            wc(C),
            AlignedPosition::Pair { top: G, bottom: T },
            wc(C),
            wc(A),
        ];
        let dx = AlignedDuplex::from_columns(
            Sequence::parse("ACGCA").unwrap(),
            Sequence::parse("TGTGT").unwrap(),
            cols,
        )
        .unwrap();
        let features = extract_features(&dx).unwrap();
        assert_eq!(features.len(), 3);
        match &features[1] {
            Feature::Mismatch { top, bottom, idx } => {
                assert_eq!(*top, G);
                assert_eq!(*bottom, T);
                assert_eq!(*idx, 2);
            }
            _ => panic!("expected Mismatch, got {:?}", features[1]),
        }
    }

    #[test]
    fn tandem_mismatch_becomes_internal_loop() {
        // Two consecutive non-WC pairs → InternalLoop length-2.
        let cols = vec![
            wc(A),
            AlignedPosition::Pair { top: G, bottom: T },
            AlignedPosition::Pair { top: A, bottom: A },
            wc(T),
        ];
        let dx = AlignedDuplex::from_columns(
            Sequence::parse("AGAT").unwrap(),
            Sequence::parse("ATTT").unwrap(),
            cols,
        )
        .unwrap();
        let features = extract_features(&dx).unwrap();
        assert_eq!(features.len(), 3);
        match &features[1] {
            Feature::InternalLoop { top_bases, bottom_bases, .. } => {
                assert_eq!(top_bases, &vec![G, A]);
                assert_eq!(bottom_bases, &vec![T, A]);
            }
            _ => panic!("expected InternalLoop, got {:?}", features[1]),
        }
    }

    #[test]
    fn mixed_bulge_and_pair_becomes_internal_loop() {
        // A run with both BulgeTop and BulgeBottom → InternalLoop.
        let cols = vec![
            wc(A),
            wc(C),
            AlignedPosition::BulgeTop(G),
            AlignedPosition::BulgeBottom(A),
            wc(T),
            wc(A),
        ];
        let dx = AlignedDuplex::from_columns(
            Sequence::parse("ACGTA").unwrap(),
            Sequence::parse("TATGT").unwrap(),
            cols,
        )
        .unwrap();
        let features = extract_features(&dx).unwrap();
        assert_eq!(features.len(), 3);
        match &features[1] {
            Feature::InternalLoop { top_bases, bottom_bases, .. } => {
                assert_eq!(top_bases, &vec![G]);
                assert_eq!(bottom_bases, &vec![A]);
            }
            _ => panic!("expected InternalLoop, got {:?}", features[1]),
        }
    }
}
