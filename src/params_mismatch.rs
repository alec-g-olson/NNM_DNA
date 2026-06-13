//! Nearest-neighbor parameters for internal single-base mismatches in DNA.
//!
//! Supports all 8 mismatch types:
//!
//! | Mismatch | Source paper |
//! |----------|--------------|
//! | G·T      | Allawi & SantaLucia (1997), *Biochemistry* 36:10581–10594 |
//! | G·A      | Allawi & SantaLucia (1998), *Biochemistry* 37:2170–2179   |
//! | C·T      | Allawi & SantaLucia (1998), *Nucl. Acids Res.* 26:2694–2701 |
//! | A·C      | Allawi & SantaLucia (1998), *Biochemistry* 37:9435–9444 (pH 7.0) |
//! | A·A      | Peyret, Seneviratne, Allawi & SantaLucia (1999), *Biochemistry* 38:3468–3477 |
//! | C·C      | (same as above) |
//! | G·G      | (same as above) |
//! | T·T      | (same as above) |
//!
//! Values reproduced from Biopython's `DNA_IMM1` table (BSD-licensed), which
//! compiles the above papers. All entries are at 1 M NaCl.
//!
//! ## Convention
//!
//! For a nearest-neighbor stack written
//!
//! ```text
//!   5'-t1 t2-3'
//!   3'-b1 b2-5'
//! ```
//!
//! `t1` is paired antiparallel with `b1`, and `t2` with `b2`. Exactly one of
//! the two pairs must be a non-Watson-Crick pair for these parameters to apply.
//!
//! ## Strand-flip symmetry
//!
//! A duplex rotated 180° is the same physical molecule and has the same NN
//! energy. The published tables list one representative per equivalence class;
//! this module canonicalizes the input under the mapping `(t1,t2,b1,b2)
//! ↔ (b2,b1,t2,t1)`.

use crate::base::Base;

/// Returns NN ΔH° (kcal/mol) and ΔS° (cal/mol·K) for a stack with exactly
/// one Watson-Crick pair and one mismatch.
///
/// Returns `None` if:
/// - both pairs are Watson-Crick (use [`crate::params::nn_params`] instead),
/// - both pairs are mismatches (tandem; out of scope for this table), or
/// - the input contains a base not covered by the unified parameter set
///   (currently this cannot happen — `Base` has only A, C, G, T).
pub fn mismatch_params(t1: Base, t2: Base, b1: Base, b2: Base) -> Option<(f64, f64)> {
    let p1_wc = b1 == t1.complement();
    let p2_wc = b2 == t2.complement();
    if p1_wc == p2_wc {
        return None;
    }
    canonical_lookup(t1, t2, b1, b2).or_else(|| canonical_lookup(b2, b1, t2, t1))
}

/// Direct table lookup; no canonicalization. Caller is responsible for
/// trying the strand-flipped form if the direct lookup misses.
fn canonical_lookup(t1: Base, t2: Base, b1: Base, b2: Base) -> Option<(f64, f64)> {
    use Base::*;
    match (t1, t2, b1, b2) {
        // ─────────────── G·T mismatches ────────────────────────────────────
        // Allawi & SantaLucia (1997), Biochemistry 36:10581
        (A, G, T, T) => Some((1.0, 0.9)),    // AG/TT
        (A, T, T, G) => Some((-2.5, -8.3)),  // AT/TG
        (C, G, G, T) => Some((-4.1, -11.7)), // CG/GT
        (C, T, G, G) => Some((-2.8, -8.0)),  // CT/GG
        (G, G, C, T) => Some((3.3, 10.4)),   // GG/CT
        (G, T, C, G) => Some((-4.4, -12.3)), // GT/CG
        (T, G, A, T) => Some((-0.1, -1.7)),  // TG/AT
        (T, T, A, G) => Some((-1.3, -5.3)),  // TT/AG

        // ─────────────── G·A mismatches ────────────────────────────────────
        // Allawi & SantaLucia (1998), Biochemistry 37:2170
        (A, A, T, G) => Some((-0.6, -2.3)),  // AA/TG
        (A, G, T, A) => Some((-0.7, -2.3)),  // AG/TA
        (C, A, G, G) => Some((-0.7, -2.3)),  // CA/GG
        (C, G, G, A) => Some((-4.0, -13.2)), // CG/GA
        (G, A, C, G) => Some((-0.6, -1.0)),  // GA/CG
        (G, G, C, A) => Some((0.5, 3.2)),    // GG/CA
        (T, A, A, G) => Some((0.7, 0.7)),    // TA/AG
        (T, G, A, A) => Some((3.0, 7.4)),    // TG/AA

        // ─────────────── C·T mismatches ────────────────────────────────────
        // Allawi & SantaLucia (1998), Nucl. Acids Res. 26:2694
        (A, C, T, T) => Some((0.7, 0.2)),    // AC/TT
        (A, T, T, C) => Some((-1.2, -6.2)),  // AT/TC
        (C, C, G, T) => Some((-0.8, -4.5)),  // CC/GT
        (C, T, G, C) => Some((-1.5, -6.1)),  // CT/GC
        (G, C, C, T) => Some((2.3, 5.4)),    // GC/CT
        (G, T, C, C) => Some((5.2, 13.5)),   // GT/CC
        (T, C, A, T) => Some((1.2, 0.7)),    // TC/AT
        (T, T, A, C) => Some((1.0, 0.7)),    // TT/AC

        // ─────────────── A·C mismatches ────────────────────────────────────
        // Allawi & SantaLucia (1998), Biochemistry 37:9435 (pH 7.0)
        (A, A, T, C) => Some((2.3, 4.6)),    // AA/TC
        (A, C, T, A) => Some((5.3, 14.6)),   // AC/TA
        (C, A, G, C) => Some((1.9, 3.7)),    // CA/GC
        (C, C, G, A) => Some((0.6, -0.6)),   // CC/GA
        (G, A, C, C) => Some((5.2, 14.2)),   // GA/CC
        (G, C, C, A) => Some((-0.7, -3.8)),  // GC/CA
        (T, A, A, C) => Some((3.4, 8.0)),    // TA/AC
        (T, C, A, A) => Some((7.6, 20.2)),   // TC/AA

        // ─────────────── A·A mismatches ────────────────────────────────────
        // Peyret, Seneviratne, Allawi & SantaLucia (1999), Biochemistry 38:3468
        (A, A, T, A) => Some((1.2, 1.7)),    // AA/TA
        (C, A, G, A) => Some((-0.9, -4.2)),  // CA/GA
        (G, A, C, A) => Some((-2.9, -9.8)),  // GA/CA
        (T, A, A, A) => Some((4.7, 12.9)),   // TA/AA

        // ─────────────── C·C mismatches ────────────────────────────────────
        // Peyret et al. (1999)
        (A, C, T, C) => Some((0.0, -4.4)),   // AC/TC
        (C, C, G, C) => Some((-1.5, -7.2)),  // CC/GC
        (G, C, C, C) => Some((3.6, 8.9)),    // GC/CC
        (T, C, A, C) => Some((6.1, 16.4)),   // TC/AC

        // ─────────────── G·G mismatches ────────────────────────────────────
        // Peyret et al. (1999)
        (A, G, T, G) => Some((-3.1, -9.5)),  // AG/TG
        (C, G, G, G) => Some((-4.9, -15.3)), // CG/GG
        (G, G, C, G) => Some((-6.0, -15.8)), // GG/CG
        (T, G, A, G) => Some((1.6, 3.6)),    // TG/AG

        // ─────────────── T·T mismatches ────────────────────────────────────
        // Peyret et al. (1999)
        (A, T, T, T) => Some((-2.7, -10.8)), // AT/TT
        (C, T, G, T) => Some((-5.0, -15.8)), // CT/GT
        (G, T, C, T) => Some((-2.2, -8.4)),  // GT/CT
        (T, T, A, T) => Some((0.2, -1.5)),   // TT/AT

        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use Base::*;

    #[test]
    fn perfect_wc_returns_none() {
        // Both WC at 5'-AG-3'/3'-TC-5'
        assert_eq!(mismatch_params(A, G, T, C), None);
    }

    #[test]
    fn tandem_mismatch_returns_none() {
        // Two G·T mismatches at 5'-GG-3'/3'-TT-5'
        assert_eq!(mismatch_params(G, G, T, T), None);
    }

    #[test]
    fn gt_at_right_position() {
        assert_eq!(mismatch_params(A, G, T, T), Some((1.0, 0.9)));
    }

    #[test]
    fn gt_at_left_position_via_strand_flip() {
        // 5'-GT-3'/3'-TA-5' is the strand-flip of 5'-AT-3'/3'-TG-5'
        assert_eq!(mismatch_params(G, T, T, A), Some((-2.5, -8.3)));
    }

    #[test]
    fn aa_mismatch_now_supported() {
        // 5'-AA-3'/3'-TA-5' — A·A mismatch at right (used to be unsupported)
        assert_eq!(mismatch_params(A, A, T, A), Some((1.2, 1.7)));
    }

    #[test]
    fn all_mismatch_types_resolve() {
        // One representative entry per mismatch type
        let cases = [
            ((A, G, T, T), (1.0, 0.9), "G·T"),
            ((A, A, T, G), (-0.6, -2.3), "G·A"),
            ((A, C, T, T), (0.7, 0.2), "C·T"),
            ((A, A, T, C), (2.3, 4.6), "A·C"),
            ((A, A, T, A), (1.2, 1.7), "A·A"),
            ((A, C, T, C), (0.0, -4.4), "C·C"),
            ((A, G, T, G), (-3.1, -9.5), "G·G"),
            ((A, T, T, T), (-2.7, -10.8), "T·T"),
        ];
        for ((t1, t2, b1, b2), expected, label) in cases {
            assert_eq!(
                mismatch_params(t1, t2, b1, b2),
                Some(expected),
                "{label} direct lookup failed"
            );
            // Strand-flip equivalent must give the same answer
            assert_eq!(
                mismatch_params(b2, b1, t2, t1),
                Some(expected),
                "{label} strand-flip lookup failed"
            );
        }
    }

    #[test]
    fn every_published_entry_resolves_via_strand_flip() {
        // For each of the 48 published canonical entries, the strand-flipped
        // input (b2, b1, t2, t1) must return the same parameters.
        let entries: &[((Base, Base, Base, Base), (f64, f64))] = &[
            // G·T
            ((A, G, T, T), (1.0, 0.9)),
            ((A, T, T, G), (-2.5, -8.3)),
            ((C, G, G, T), (-4.1, -11.7)),
            ((C, T, G, G), (-2.8, -8.0)),
            ((G, G, C, T), (3.3, 10.4)),
            ((G, T, C, G), (-4.4, -12.3)),
            ((T, G, A, T), (-0.1, -1.7)),
            ((T, T, A, G), (-1.3, -5.3)),
            // G·A
            ((A, A, T, G), (-0.6, -2.3)),
            ((A, G, T, A), (-0.7, -2.3)),
            ((C, A, G, G), (-0.7, -2.3)),
            ((C, G, G, A), (-4.0, -13.2)),
            ((G, A, C, G), (-0.6, -1.0)),
            ((G, G, C, A), (0.5, 3.2)),
            ((T, A, A, G), (0.7, 0.7)),
            ((T, G, A, A), (3.0, 7.4)),
            // C·T
            ((A, C, T, T), (0.7, 0.2)),
            ((A, T, T, C), (-1.2, -6.2)),
            ((C, C, G, T), (-0.8, -4.5)),
            ((C, T, G, C), (-1.5, -6.1)),
            ((G, C, C, T), (2.3, 5.4)),
            ((G, T, C, C), (5.2, 13.5)),
            ((T, C, A, T), (1.2, 0.7)),
            ((T, T, A, C), (1.0, 0.7)),
            // A·C
            ((A, A, T, C), (2.3, 4.6)),
            ((A, C, T, A), (5.3, 14.6)),
            ((C, A, G, C), (1.9, 3.7)),
            ((C, C, G, A), (0.6, -0.6)),
            ((G, A, C, C), (5.2, 14.2)),
            ((G, C, C, A), (-0.7, -3.8)),
            ((T, A, A, C), (3.4, 8.0)),
            ((T, C, A, A), (7.6, 20.2)),
            // A·A
            ((A, A, T, A), (1.2, 1.7)),
            ((C, A, G, A), (-0.9, -4.2)),
            ((G, A, C, A), (-2.9, -9.8)),
            ((T, A, A, A), (4.7, 12.9)),
            // C·C
            ((A, C, T, C), (0.0, -4.4)),
            ((C, C, G, C), (-1.5, -7.2)),
            ((G, C, C, C), (3.6, 8.9)),
            ((T, C, A, C), (6.1, 16.4)),
            // G·G
            ((A, G, T, G), (-3.1, -9.5)),
            ((C, G, G, G), (-4.9, -15.3)),
            ((G, G, C, G), (-6.0, -15.8)),
            ((T, G, A, G), (1.6, 3.6)),
            // T·T
            ((A, T, T, T), (-2.7, -10.8)),
            ((C, T, G, T), (-5.0, -15.8)),
            ((G, T, C, T), (-2.2, -8.4)),
            ((T, T, A, T), (0.2, -1.5)),
        ];
        assert_eq!(entries.len(), 48, "expected 48 entries total");
        for ((t1, t2, b1, b2), expected) in entries.iter().copied() {
            assert_eq!(
                mismatch_params(t1, t2, b1, b2),
                Some(expected),
                "direct lookup failed for ({t1:?}, {t2:?}, {b1:?}, {b2:?})"
            );
            assert_eq!(
                mismatch_params(b2, b1, t2, t1),
                Some(expected),
                "strand-flip lookup failed for ({t1:?}, {t2:?}, {b1:?}, {b2:?})"
            );
        }
    }
}
