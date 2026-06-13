//! SantaLucia 1998 unified nearest-neighbor parameters.
//!
//! Reference: SantaLucia, J. (1998). "A unified view of polymer, dumbbell, and
//! oligonucleotide DNA nearest-neighbor thermodynamics." PNAS 95(4):1460-1465.
//!
//! All ΔH° values are in kcal/mol. All ΔS° values are in cal/(mol·K).
//! Parameters apply at 1 M NaCl. Salt corrections are not included here.

use crate::base::Base;

/// Nearest-neighbor enthalpy/entropy contribution for the stack 5'-b1 b2-3'
/// paired with its complement 3'-b1' b2'-5'.
///
/// Returns `(ΔH°, ΔS°)` in `(kcal/mol, cal/(mol·K))`.
///
/// There are 10 unique parameters (not 16) because each NN is equivalent under
/// strand-flip to another. For instance, 5'-AA-3'/3'-TT-5' is the same duplex
/// as 5'-TT-3'/3'-AA-5' read upside-down, so AA and TT (as top-strand
/// dinucleotides) share a parameter.
pub const fn nn_params(b1: Base, b2: Base) -> (f64, f64) {
    use Base::*;
    match (b1, b2) {
        (A, A) | (T, T) => (-7.9, -22.2),
        (A, T)          => (-7.2, -20.4),
        (T, A)          => (-7.2, -21.3),
        (C, A) | (T, G) => (-8.5, -22.7),
        (G, T) | (A, C) => (-8.4, -22.4),
        (C, T) | (A, G) => (-7.8, -21.0),
        (G, A) | (T, C) => (-8.2, -22.2),
        (C, G)          => (-10.6, -27.2),
        (G, C)          => (-9.8, -24.4),
        (G, G) | (C, C) => (-8.0, -19.9),
    }
}

/// Initiation contribution applied once per duplex terminus whose terminal
/// base pair is G·C.
pub const INIT_GC: (f64, f64) = (0.1, -2.8);

/// Initiation contribution applied once per duplex terminus whose terminal
/// base pair is A·T.
pub const INIT_AT: (f64, f64) = (2.3, 4.1);

/// Symmetry correction applied once if the duplex is self-complementary.
pub const SYMMETRY_CORRECTION: (f64, f64) = (0.0, -1.4);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strand_flip_equivalence() {
        // The 6 NN pairs that share parameters with their reverse-complement.
        for (a, b) in [
            (Base::A, Base::A),
            (Base::C, Base::A),
            (Base::G, Base::T),
            (Base::C, Base::T),
            (Base::G, Base::A),
            (Base::G, Base::G),
        ] {
            let p1 = nn_params(a, b);
            // reverse-complement of (a, b) read as a dinucleotide is
            // (complement(b), complement(a)).
            let p2 = nn_params(b.complement(), a.complement());
            assert_eq!(p1, p2, "{a:?}{b:?} vs its reverse-complement should share parameters");
        }
    }

    #[test]
    fn palindromic_nn_pairs() {
        // AT, TA, CG, GC are palindromic dinucleotides — their reverse-complement
        // is themselves, so they sit alone in the 10-entry set.
        for (a, b) in [(Base::A, Base::T), (Base::T, Base::A), (Base::C, Base::G), (Base::G, Base::C)] {
            assert_eq!(nn_params(a, b), nn_params(b.complement(), a.complement()));
        }
    }
}
