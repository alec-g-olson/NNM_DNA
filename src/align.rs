//! Needleman-Wunsch alignment of two DNA strands to recover the column-level
//! structure of a duplex containing bulges and/or internal loops.
//!
//! ## The antiparallel trick
//!
//! Two strands of a DNA duplex are antiparallel. To pose the alignment as a
//! standard NW problem on two same-direction strings, we align the top strand
//! (5'→3') against the **reverse complement** of the bottom strand (also
//! 5'→3'). In this projection, *identical bases* at corresponding positions
//! correspond to physical *Watson-Crick pairs* in the duplex — because RC
//! places `complement(bottom[n−1−j])` at position `j`, which is exactly the
//! base top must equal for a WC pair.
//!
//! The aligner therefore uses plain base equality as the match signal. The
//! `Pair` columns it emits store the *physically paired* bottom base
//! (`complement(b[j])`), so downstream code never needs to redo the RC
//! reasoning.
//!
//! ## Affine gap penalty (Gotoh)
//!
//! Three matrices M, X, Y track the best path ending in (mis)match, gap-on-`b`
//! (top-strand bulge), and gap-on-`a` (bottom-strand bulge) respectively.
//! Gap-open and gap-extend are separated so multi-base bulges stay coherent.
//!
//! ## Scoring convention
//!
//! All scores are in **kcal/mol-equivalent ΔG° at 37 °C, 1 M NaCl** and the
//! algorithm **minimizes** (low score = stable structure). Sign is opposite
//! of the traditional sequence-alignment convention. Constants are tunable
//! but exposed as `pub const` so callers can see the values in use.

use crate::base::Base;
use crate::sequence::Sequence;
use std::fmt;

/// Score for a Watson-Crick paired column. Negative because WC is favorable.
/// Roughly half the magnitude of an average WC NN stack ΔG°₃₇.
pub const SCORE_WC_MATCH: f64 = -1.5;

/// Score for a non-WC paired column (mismatched bases at the same column).
/// Positive — mismatches are tolerable but penalized.
pub const SCORE_MISMATCH: f64 = 1.0;

/// Cost of opening a new gap (first bulged base). Roughly half the
/// Table-4 length-1 bulge penalty so two flanking columns share it.
pub const SCORE_GAP_OPEN: f64 = 3.0;

/// Cost of extending an existing gap by one base. Much cheaper than open;
/// reflects that multi-base bulges pay a length penalty, not a per-base one.
pub const SCORE_GAP_EXTEND: f64 = 0.5;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignedPosition {
    /// Both strands have a base at this column. `bottom` is the *physically
    /// paired* bottom base (not its RC). May be WC (when `bottom == top.complement()`)
    /// or a mismatch.
    Pair { top: Base, bottom: Base },
    /// Top strand has a base; bottom is gapped at this column.
    BulgeTop(Base),
    /// Bottom strand has a base; top is gapped at this column. The stored
    /// base is the physical bottom-strand base (not its RC).
    BulgeBottom(Base),
}

impl AlignedPosition {
    pub fn is_pair(&self) -> bool {
        matches!(self, AlignedPosition::Pair { .. })
    }

    pub fn is_wc_pair(&self) -> bool {
        match self {
            AlignedPosition::Pair { top, bottom } => *bottom == top.complement(),
            _ => false,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlignedDuplex {
    top: Sequence,
    bottom: Sequence,
    columns: Vec<AlignedPosition>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlignError {
    EmptyStrand,
    /// The DP didn't find any finite-score path. Cannot happen for two
    /// non-empty `Sequence`s with the recommended scoring; reserved for
    /// future score schemes that could produce it.
    NoAlignment,
    /// The alignment's first or last column is a `BulgeTop` / `BulgeBottom`
    /// (no terminal base pair). Terminal overhangs are out of scope for v1.
    TerminalBulge,
    /// The alignment's first or last column is a `Pair` whose bases don't
    /// form a Watson-Crick pair. Terminal mismatches are out of scope for v1.
    TerminalMismatch,
}

impl fmt::Display for AlignError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AlignError::EmptyStrand => f.write_str("empty strand passed to aligner"),
            AlignError::NoAlignment => f.write_str("no valid alignment found"),
            AlignError::TerminalBulge => f.write_str("terminal bulge: duplex must start and end with a base pair"),
            AlignError::TerminalMismatch => f.write_str("terminal mismatch: duplex must start and end with a Watson-Crick pair"),
        }
    }
}

impl std::error::Error for AlignError {}

impl AlignedDuplex {
    /// Run Needleman-Wunsch on `top` (5'→3') against `bottom.reverse_complement()`
    /// (also 5'→3'). The result is validated to have WC pairs at both termini.
    pub fn align(top: Sequence, bottom: Sequence) -> Result<Self, AlignError> {
        if top.is_empty() || bottom.is_empty() {
            return Err(AlignError::EmptyStrand);
        }
        let rc = bottom.reverse_complement();
        let columns = needleman_wunsch(top.as_slice(), rc.as_slice())?;
        Self::from_columns(top, bottom, columns)
    }

    /// Construct an `AlignedDuplex` from a pre-computed column list. Useful
    /// for tests and for callers who computed the alignment another way.
    /// Validates termini.
    pub fn from_columns(
        top: Sequence,
        bottom: Sequence,
        columns: Vec<AlignedPosition>,
    ) -> Result<Self, AlignError> {
        if columns.is_empty() {
            return Err(AlignError::EmptyStrand);
        }
        for terminal in [*columns.first().unwrap(), *columns.last().unwrap()] {
            match terminal {
                AlignedPosition::Pair { top: t, bottom: b } => {
                    if b != t.complement() {
                        return Err(AlignError::TerminalMismatch);
                    }
                }
                _ => return Err(AlignError::TerminalBulge),
            }
        }
        Ok(AlignedDuplex { top, bottom, columns })
    }

    pub fn columns(&self) -> &[AlignedPosition] {
        &self.columns
    }

    pub fn top(&self) -> &Sequence {
        &self.top
    }

    pub fn bottom(&self) -> &Sequence {
        &self.bottom
    }

    /// Permissive self-complementarity: true iff the two strands written 5'→3'
    /// are identical. A bulged or mismatched duplex with `top == bottom` still
    /// has 180° rotational symmetry, so the symmetry correction applies.
    pub fn is_self_complementary(&self) -> bool {
        self.top == self.bottom
    }
}

#[derive(Clone, Copy)]
enum State {
    M,
    X,
    Y,
}

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

fn needleman_wunsch(a: &[Base], b: &[Base]) -> Result<Vec<AlignedPosition>, AlignError> {
    let m = a.len();
    let n = b.len();
    let inf = f64::INFINITY;

    // M[i][j] = best score aligning a[..i] and b[..j] ending in a (mis)match column
    // X[i][j] = best score ending in a gap on b  (a consumes, b doesn't) → BulgeTop
    // Y[i][j] = best score ending in a gap on a  (b consumes, a doesn't) → BulgeBottom
    let mut mm = vec![vec![inf; n + 1]; m + 1];
    let mut xx = vec![vec![inf; n + 1]; m + 1];
    let mut yy = vec![vec![inf; n + 1]; m + 1];

    mm[0][0] = 0.0;
    for i in 1..=m {
        xx[i][0] = SCORE_GAP_OPEN + (i - 1) as f64 * SCORE_GAP_EXTEND;
    }
    for j in 1..=n {
        yy[0][j] = SCORE_GAP_OPEN + (j - 1) as f64 * SCORE_GAP_EXTEND;
    }

    for i in 1..=m {
        for j in 1..=n {
            let s = if a[i - 1] == b[j - 1] {
                SCORE_WC_MATCH
            } else {
                SCORE_MISMATCH
            };
            let prev_diag = mm[i - 1][j - 1].min(xx[i - 1][j - 1]).min(yy[i - 1][j - 1]);
            mm[i][j] = s + prev_diag;

            xx[i][j] = (mm[i - 1][j] + SCORE_GAP_OPEN)
                .min(xx[i - 1][j] + SCORE_GAP_EXTEND)
                .min(yy[i - 1][j] + SCORE_GAP_OPEN);

            yy[i][j] = (mm[i][j - 1] + SCORE_GAP_OPEN)
                .min(yy[i][j - 1] + SCORE_GAP_EXTEND)
                .min(xx[i][j - 1] + SCORE_GAP_OPEN);
        }
    }

    if mm[m][n].is_infinite() {
        return Err(AlignError::NoAlignment);
    }

    // Traceback from M[m][n]. Termini must pair, so we enter from M.
    let mut columns = Vec::with_capacity(m + n);
    let mut i = m;
    let mut j = n;
    let mut state = State::M;

    while i > 0 || j > 0 {
        match state {
            State::M => {
                let s = if a[i - 1] == b[j - 1] {
                    SCORE_WC_MATCH
                } else {
                    SCORE_MISMATCH
                };
                let bot_base = b[j - 1].complement();
                columns.push(AlignedPosition::Pair {
                    top: a[i - 1],
                    bottom: bot_base,
                });
                let target = mm[i][j] - s;
                state = if approx_eq(mm[i - 1][j - 1], target) {
                    State::M
                } else if approx_eq(xx[i - 1][j - 1], target) {
                    State::X
                } else {
                    State::Y
                };
                i -= 1;
                j -= 1;
            }
            State::X => {
                columns.push(AlignedPosition::BulgeTop(a[i - 1]));
                let xv = xx[i][j];
                state = if approx_eq(mm[i - 1][j] + SCORE_GAP_OPEN, xv) {
                    State::M
                } else if approx_eq(xx[i - 1][j] + SCORE_GAP_EXTEND, xv) {
                    State::X
                } else {
                    State::Y
                };
                i -= 1;
            }
            State::Y => {
                let bot_base = b[j - 1].complement();
                columns.push(AlignedPosition::BulgeBottom(bot_base));
                let yv = yy[i][j];
                state = if approx_eq(mm[i][j - 1] + SCORE_GAP_OPEN, yv) {
                    State::M
                } else if approx_eq(yy[i][j - 1] + SCORE_GAP_EXTEND, yv) {
                    State::Y
                } else {
                    State::X
                };
                j -= 1;
            }
        }
    }

    columns.reverse();
    Ok(columns)
}

#[cfg(test)]
mod tests {
    use super::*;
    use Base::*;

    #[test]
    fn align_perfect_wc_returns_all_pairs() {
        let top = Sequence::parse("ACGT").unwrap();
        let bot = Sequence::parse("ACGT").unwrap(); // RC of "ACGT" = "ACGT"
        let dx = AlignedDuplex::align(top, bot).unwrap();
        assert_eq!(dx.columns().len(), 4);
        for col in dx.columns() {
            assert!(col.is_wc_pair());
        }
    }

    #[test]
    fn align_six_mer_no_bulges() {
        // CGTTGA / TCAACG → equal length WC duplex
        let top = Sequence::parse("CGTTGA").unwrap();
        let bot = Sequence::parse("TCAACG").unwrap();
        let dx = AlignedDuplex::align(top, bot).unwrap();
        assert_eq!(dx.columns().len(), 6);
        assert!(dx.columns().iter().all(|c| c.is_wc_pair()));
    }

    #[test]
    fn from_columns_rejects_terminal_bulge() {
        let top = Sequence::parse("AT").unwrap();
        let bot = Sequence::parse("AT").unwrap();
        let cols = vec![
            AlignedPosition::BulgeTop(A),
            AlignedPosition::Pair { top: T, bottom: A },
        ];
        assert_eq!(
            AlignedDuplex::from_columns(top, bot, cols),
            Err(AlignError::TerminalBulge)
        );
    }

    #[test]
    fn from_columns_rejects_terminal_mismatch() {
        let top = Sequence::parse("AT").unwrap();
        let bot = Sequence::parse("AT").unwrap();
        let cols = vec![
            AlignedPosition::Pair { top: A, bottom: A }, // A·A mismatch at terminus
            AlignedPosition::Pair { top: T, bottom: A },
        ];
        assert_eq!(
            AlignedDuplex::from_columns(top, bot, cols),
            Err(AlignError::TerminalMismatch)
        );
    }

    #[test]
    fn from_columns_accepts_internal_bulge() {
        // 5'-A C T-3'  with bulged C on top
        // 3'-T   A-5'   (bottom 5'→3' = "AT")
        let top = Sequence::parse("ACT").unwrap();
        let bot = Sequence::parse("AT").unwrap();
        let cols = vec![
            AlignedPosition::Pair { top: A, bottom: T },
            AlignedPosition::BulgeTop(C),
            AlignedPosition::Pair { top: T, bottom: A },
        ];
        let dx = AlignedDuplex::from_columns(top, bot, cols).unwrap();
        assert_eq!(dx.columns().len(), 3);
    }
}
