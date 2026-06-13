//! Validated DNA sequences for the nearest-neighbor model.

use crate::base::Base;
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Sequence(Vec<Base>);

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseError {
    Empty,
    TooShort(usize),
    InvalidBase { ch: char, position: usize },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::Empty => f.write_str("empty sequence"),
            ParseError::TooShort(n) => {
                write!(f, "sequence too short for nearest-neighbor model (need >= 2 bases, got {n})")
            }
            ParseError::InvalidBase { ch, position } => {
                write!(f, "invalid base {ch:?} at position {position}")
            }
        }
    }
}

impl std::error::Error for ParseError {}

impl Sequence {
    pub fn parse(s: &str) -> Result<Self, ParseError> {
        if s.is_empty() {
            return Err(ParseError::Empty);
        }
        let mut bases = Vec::with_capacity(s.len());
        for (i, ch) in s.chars().enumerate() {
            let b = Base::from_char(ch).ok_or(ParseError::InvalidBase { ch, position: i })?;
            bases.push(b);
        }
        if bases.len() < 2 {
            return Err(ParseError::TooShort(bases.len()));
        }
        Ok(Sequence(bases))
    }

    pub fn from_bases(bases: Vec<Base>) -> Result<Self, ParseError> {
        if bases.len() < 2 {
            return Err(ParseError::TooShort(bases.len()));
        }
        Ok(Sequence(bases))
    }

    pub fn as_slice(&self) -> &[Base] {
        &self.0
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn reverse_complement(&self) -> Sequence {
        Sequence(self.0.iter().rev().map(|b| b.complement()).collect())
    }

    pub fn is_self_complementary(&self) -> bool {
        if self.0.len() % 2 != 0 {
            return false;
        }
        let n = self.0.len();
        (0..n / 2).all(|i| self.0[i] == self.0[n - 1 - i].complement())
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for b in &self.0 {
            write!(f, "{b}")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_valid() {
        let s = Sequence::parse("ACGT").unwrap();
        assert_eq!(s.len(), 4);
        assert_eq!(s.to_string(), "ACGT");
    }

    #[test]
    fn parse_lowercase() {
        let s = Sequence::parse("acgt").unwrap();
        assert_eq!(s.to_string(), "ACGT");
    }

    #[test]
    fn parse_rejects_invalid() {
        match Sequence::parse("ACNT") {
            Err(ParseError::InvalidBase { ch: 'N', position: 2 }) => {}
            other => panic!("unexpected: {other:?}"),
        }
    }

    #[test]
    fn parse_rejects_too_short() {
        assert!(matches!(Sequence::parse("A"), Err(ParseError::TooShort(1))));
        assert!(matches!(Sequence::parse(""), Err(ParseError::Empty)));
    }

    #[test]
    fn reverse_complement_basic() {
        let s = Sequence::parse("ACGT").unwrap();
        assert_eq!(s.reverse_complement().to_string(), "ACGT"); // palindrome
        let s = Sequence::parse("CGTTGA").unwrap();
        assert_eq!(s.reverse_complement().to_string(), "TCAACG");
    }

    #[test]
    fn self_complementary() {
        assert!(Sequence::parse("CGCGCG").unwrap().is_self_complementary());
        assert!(Sequence::parse("ACGT").unwrap().is_self_complementary());
        assert!(!Sequence::parse("CGTTGA").unwrap().is_self_complementary());
        // odd length cannot be self-complementary
        assert!(!Sequence::parse("ACG").unwrap().is_self_complementary());
    }
}
