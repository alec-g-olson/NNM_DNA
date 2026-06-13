//! Duplex: two antiparallel-paired DNA strands.
//!
//! In this crate's notation, both strands are stored 5'→3'. The bottom strand
//! is interpreted as paired antiparallel with the top, so `top[i]` pairs with
//! `bottom[len - 1 - i]`.
//!
//! Current scope: equal-length strands (no overhangs), Watson-Crick at both
//! termini (no terminal mismatches). Internal positions may contain mismatches.

use crate::base::Base;
use crate::sequence::Sequence;
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Duplex {
    top: Sequence,
    bottom: Sequence,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DuplexError {
    LengthMismatch { top: usize, bottom: usize },
    TooShort(usize),
    TerminalMismatch,
}

impl fmt::Display for DuplexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DuplexError::LengthMismatch { top, bottom } => {
                write!(f, "strand length mismatch: top {top} bp, bottom {bottom} bp")
            }
            DuplexError::TooShort(n) => write!(
                f,
                "duplex too short for nearest-neighbor model (need >= 2 bp, got {n})"
            ),
            DuplexError::TerminalMismatch => {
                f.write_str("terminal mismatches are not supported")
            }
        }
    }
}

impl std::error::Error for DuplexError {}

impl Duplex {
    /// Perfect Watson-Crick duplex from a single strand. The bottom strand is
    /// derived by reverse-complementing `top`.
    pub fn perfect(top: Sequence) -> Self {
        let bottom = top.reverse_complement();
        Duplex { top, bottom }
    }

    /// Duplex from two strands, both written 5'→3'. They will be paired
    /// antiparallel: `top[i]` pairs with `bottom[len - 1 - i]`.
    pub fn new(top: Sequence, bottom: Sequence) -> Result<Self, DuplexError> {
        if top.len() != bottom.len() {
            return Err(DuplexError::LengthMismatch {
                top: top.len(),
                bottom: bottom.len(),
            });
        }
        if top.len() < 2 {
            return Err(DuplexError::TooShort(top.len()));
        }
        let dx = Duplex { top, bottom };
        let (t0, b0) = dx.pair_at(0);
        let (tn, bn) = dx.pair_at(dx.len() - 1);
        if b0 != t0.complement() || bn != tn.complement() {
            return Err(DuplexError::TerminalMismatch);
        }
        Ok(dx)
    }

    pub fn len(&self) -> usize {
        self.top.len()
    }

    pub fn is_empty(&self) -> bool {
        self.top.is_empty()
    }

    pub fn top(&self) -> &Sequence {
        &self.top
    }

    pub fn bottom(&self) -> &Sequence {
        &self.bottom
    }

    /// `(top_base, bottom_base)` at top-strand position `i`. `bottom_base` is
    /// the base actually paired (antiparallel) with `top[i]`, located at index
    /// `len - 1 - i` of the 5'→3' bottom strand.
    pub fn pair_at(&self, i: usize) -> (Base, Base) {
        let n = self.len();
        let t = self.top.as_slice()[i];
        let b = self.bottom.as_slice()[n - 1 - i];
        (t, b)
    }

    /// True if all positions are Watson-Crick.
    pub fn is_perfect(&self) -> bool {
        (0..self.len()).all(|i| {
            let (t, b) = self.pair_at(i);
            b == t.complement()
        })
    }

    /// True iff the duplex is invariant under 180° rotation (top↔bottom swap).
    /// Equivalently, both strands written 5'→3' are the same string. This is
    /// what triggers the symmetry correction in the Tm/ΔG formulas.
    ///
    /// Note: for a `Duplex::perfect(top)` this reduces to "top is its own
    /// reverse-complement," i.e. `top` is a palindromic sequence.
    pub fn is_self_complementary(&self) -> bool {
        self.top == self.bottom
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_matches_reverse_complement() {
        let top = Sequence::parse("ACGTACGT").unwrap();
        let dx = Duplex::perfect(top.clone());
        assert_eq!(dx.bottom(), &top.reverse_complement());
        assert!(dx.is_perfect());
    }

    #[test]
    fn new_validates_lengths() {
        let top = Sequence::parse("ACGT").unwrap();
        let bot = Sequence::parse("ACG").unwrap();
        assert!(matches!(
            Duplex::new(top, bot),
            Err(DuplexError::LengthMismatch { top: 4, bottom: 3 })
        ));
    }

    #[test]
    fn new_rejects_terminal_mismatch() {
        // top "ACGT", bottom 5'→3' "ACGA" → top[3]=T pairs with bottom[0]=A: T·A WC.
        // top[0]=A pairs with bottom[3]=A: A·A mismatch at 5' terminus → reject.
        let top = Sequence::parse("ACGT").unwrap();
        let bot = Sequence::parse("ACGA").unwrap();
        assert!(matches!(
            Duplex::new(top, bot),
            Err(DuplexError::TerminalMismatch)
        ));
    }

    #[test]
    fn pair_at_antiparallel() {
        // Standard WC duplex 5'-ACGT-3' / 3'-TGCA-5'. Bottom 5'→3' is "ACGT".
        let dx = Duplex::perfect(Sequence::parse("ACGT").unwrap());
        assert_eq!(dx.pair_at(0), (Base::A, Base::T));
        assert_eq!(dx.pair_at(1), (Base::C, Base::G));
        assert_eq!(dx.pair_at(2), (Base::G, Base::C));
        assert_eq!(dx.pair_at(3), (Base::T, Base::A));
    }

    #[test]
    fn detects_internal_mismatch() {
        // 5'-CGTGC-3' / 3'-GTACG-5'  (bottom 5'→3' = "GCATG")
        // pair_at(1) = (G, T) → G·T mismatch.
        let top = Sequence::parse("CGTGC").unwrap();
        let bot = Sequence::parse("GCATG").unwrap();
        let dx = Duplex::new(top, bot).unwrap();
        assert!(!dx.is_perfect());
        assert_eq!(dx.pair_at(1), (Base::G, Base::T));
        // Termini are still WC: C·G at both ends.
        assert_eq!(dx.pair_at(0), (Base::C, Base::G));
        assert_eq!(dx.pair_at(4), (Base::C, Base::G));
    }
}
