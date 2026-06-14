//! N3DNA — Rust implementation of the SantaLucia 1998 unified nearest-neighbor
//! thermodynamic model for DNA duplex stability.
//!
//! Reference: SantaLucia, J. (1998). PNAS 95(4):1460-1465.
//!
//! # Scope
//!
//! - Perfect Watson-Crick duplexes (blunt-ended, equal-length strands).
//! - Internal single-base **G·T mismatches** (Allawi & SantaLucia 1997).
//!
//! Parameters apply at 1 M NaCl. Other mismatch types, dangling ends, and
//! salt correction are planned but not yet implemented.
//!
//! # Example
//!
//! ```
//! use nnm_dna::{Duplex, Sequence, duplex_thermo};
//!
//! // Perfect duplex from one strand.
//! let dx = Duplex::perfect(Sequence::parse("CGTTGA").unwrap());
//! let t = duplex_thermo(&dx).unwrap();
//! assert!((t.dh - (-41.2)).abs() < 0.05);
//!
//! // Duplex with one internal G·T mismatch.
//! //   5'-CGTGC-3'
//! //   3'-GTACG-5'    (bottom written 5'→3' = "GCATG")
//! let dx = Duplex::new(
//!     Sequence::parse("CGTGC").unwrap(),
//!     Sequence::parse("GCATG").unwrap(),
//! ).unwrap();
//! let t = duplex_thermo(&dx).unwrap();
//! assert!(!dx.is_perfect());
//! ```

pub mod align;
pub mod base;
pub mod duplex;
pub mod feature;
pub mod params;
pub mod params_loop;
pub mod params_mismatch;
pub mod sequence;
pub mod thermo;

pub use align::{AlignError, AlignedDuplex, AlignedPosition};
pub use base::Base;
pub use duplex::{Duplex, DuplexError};
pub use feature::{extract_features, Feature};
pub use sequence::{ParseError, Sequence};
pub use thermo::{
    aligned_duplex_thermo, duplex_thermo, Equilibrium, ReactionConditions, SaltCorrection, Thermo,
    ThermoError, R,
};
