//! Single DNA bases and their basic operations.

use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum Base {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Base {
    pub const fn complement(self) -> Self {
        match self {
            Base::A => Base::T,
            Base::T => Base::A,
            Base::C => Base::G,
            Base::G => Base::C,
        }
    }

    pub fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'A' => Some(Base::A),
            'C' => Some(Base::C),
            'G' => Some(Base::G),
            'T' => Some(Base::T),
            _ => None,
        }
    }

    pub const fn to_char(self) -> char {
        match self {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
        }
    }

    /// True if this base, when paired with its complement at a duplex terminus,
    /// contributes an A·T initiation term (vs G·C).
    pub const fn is_at(self) -> bool {
        matches!(self, Base::A | Base::T)
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Base::A => "A",
            Base::C => "C",
            Base::G => "G",
            Base::T => "T",
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_pairs() {
        assert_eq!(Base::A.complement(), Base::T);
        assert_eq!(Base::T.complement(), Base::A);
        assert_eq!(Base::C.complement(), Base::G);
        assert_eq!(Base::G.complement(), Base::C);
    }

    #[test]
    fn complement_is_involutive() {
        for b in [Base::A, Base::C, Base::G, Base::T] {
            assert_eq!(b.complement().complement(), b);
        }
    }

    #[test]
    fn parse_case_insensitive() {
        assert_eq!(Base::from_char('a'), Some(Base::A));
        assert_eq!(Base::from_char('T'), Some(Base::T));
        assert_eq!(Base::from_char('U'), None);
        assert_eq!(Base::from_char('N'), None);
    }
}
