//! Alleles.
//!
//! Ported from `htsjdk.variant.variantcontext.Allele` and `SimpleAllele` at htsjdk 4.2.0.
//!
//! The construction path decides what the file will say, and it is not a pass-through. Two
//! rules in `SimpleAllele`'s constructor change the bases before anything is written:
//!
//!   * plain bases are **uppercased in place**, so `acgt` is stored and printed as `ACGT`;
//!   * **symbolic** bases are not, so `<del>` keeps its case.
//!
//! and "symbolic" is decided by `wouldBeSymbolicAllele`, which is broader than it sounds. It
//! is true for anything starting with `<`, anything *ending* with `>`, anything containing
//! `[` or `]`, and anything starting or ending with `.`. So `AT>` is a symbolic allele and is
//! not uppercased, while `AT` is not and is. A port that uppercased unconditionally, or not at
//! all, would produce a valid VCF that differs from htsjdk's on exactly those inputs.

use std::fmt;

/// `.`
pub const NO_CALL_STRING: &str = ".";
/// `*`
pub const SPAN_DEL_STRING: &str = "*";

const SINGLE_BREAKEND_INDICATOR: u8 = b'.';
const BREAKEND_EXTENDING_RIGHT: u8 = b'[';
const BREAKEND_EXTENDING_LEFT: u8 = b']';
const SYMBOLIC_ALLELE_START: u8 = b'<';
const SYMBOLIC_ALLELE_END: u8 = b'>';

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AlleleError {
    NullAllele,
    NoCallAsReference,
    SymbolicAsReference,
    SpanDelAsReference,
    UnacceptableBases(String),
}

impl fmt::Display for AlleleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NullAllele => write!(f, "Null alleles are not supported"),
            Self::NoCallAsReference => {
                write!(f, "Cannot tag a NoCall allele as the reference allele")
            }
            Self::SymbolicAsReference => {
                write!(f, "Cannot tag a symbolic allele as the reference allele")
            }
            Self::SpanDelAsReference => write!(
                f,
                "Cannot tag a spanning deletions allele as the reference allele"
            ),
            Self::UnacceptableBases(b) => write!(f, "Unexpected base in allele bases '{b}'"),
        }
    }
}

/// An allele: bases plus a reference flag.
///
/// Equality follows htsjdk's: the reference flag is part of it. Two alleles with the same bases
/// and different reference status are different alleles, which matters because `VCFEncoder`
/// looks genotype alleles up in a map keyed by allele to find their GT index.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Allele {
    bases: Vec<u8>,
    is_ref: bool,
    is_no_call: bool,
    is_symbolic: bool,
}

/// `bases.length == 1 && bases[0] == '-'`, or empty.
fn would_be_null(bases: &[u8]) -> bool {
    (bases.len() == 1 && bases[0] == b'-') || bases.is_empty()
}

fn would_be_no_call(bases: &[u8]) -> bool {
    bases.len() == 1 && bases[0] == b'.'
}

fn would_be_star(bases: &[u8]) -> bool {
    bases.len() == 1 && bases[0] == b'*'
}

fn would_be_breakpoint(bases: &[u8]) -> bool {
    bases.len() > 1
        && bases
            .iter()
            .any(|&b| b == BREAKEND_EXTENDING_LEFT || b == BREAKEND_EXTENDING_RIGHT)
}

fn would_be_single_breakend(bases: &[u8]) -> bool {
    bases.len() > 1
        && (bases[0] == SINGLE_BREAKEND_INDICATOR
            || bases[bases.len() - 1] == SINGLE_BREAKEND_INDICATOR)
}

/// The test that decides whether the bases are uppercased.
///
/// Note the `bases.len() <= 1` guard: a single `<` is *not* symbolic, and so is rejected as an
/// unacceptable base rather than accepted as a symbolic allele.
pub fn would_be_symbolic(bases: &[u8]) -> bool {
    if bases.len() <= 1 {
        false
    } else {
        bases[0] == SYMBOLIC_ALLELE_START
            || bases[bases.len() - 1] == SYMBOLIC_ALLELE_END
            || would_be_breakpoint(bases)
            || would_be_single_breakend(bases)
    }
}

fn acceptable_bases(bases: &[u8], is_reference: bool) -> bool {
    if would_be_null(bases) {
        return false;
    }
    if would_be_no_call(bases) || would_be_symbolic(bases) {
        return true;
    }
    if would_be_star(bases) {
        return !is_reference;
    }
    bases.iter().all(|b| {
        matches!(
            b,
            b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' | b'N' | b'n'
        )
    })
}

impl Allele {
    /// `Allele.create(bases, isRef)`.
    ///
    /// The single-base fast path in htsjdk returns interned constants for `A C G T N` in both
    /// cases and for `.` and `*`; the result is the same object this builds, so the fast path
    /// is an allocation optimisation and not a behavioural branch. The one place it *is*
    /// behavioural is the error: htsjdk reports `Illegal base [X] seen in the allele` for a bad
    /// single base and `Unexpected base in allele bases 'XY'` for a longer one. Only the
    /// message differs, so both map here to `UnacceptableBases`.
    pub fn create(bases: &[u8], is_ref: bool) -> Result<Self, AlleleError> {
        if would_be_null(bases) {
            return Err(AlleleError::NullAllele);
        }
        if would_be_no_call(bases) {
            if is_ref {
                return Err(AlleleError::NoCallAsReference);
            }
            // htsjdk represents a no-call as *no bases*, not as a '.' base.
            return Ok(Self {
                bases: Vec::new(),
                is_ref: false,
                is_no_call: true,
                is_symbolic: false,
            });
        }
        if would_be_star(bases) && is_ref {
            return Err(AlleleError::SpanDelAsReference);
        }

        let is_symbolic = would_be_symbolic(bases);
        if is_symbolic && is_ref {
            return Err(AlleleError::SymbolicAsReference);
        }
        let bases: Vec<u8> = if is_symbolic {
            bases.to_vec()
        } else {
            bases.to_ascii_uppercase()
        };

        if !acceptable_bases(&bases, is_ref) {
            return Err(AlleleError::UnacceptableBases(
                String::from_utf8_lossy(&bases).into_owned(),
            ));
        }
        Ok(Self {
            bases,
            is_ref,
            is_no_call: false,
            is_symbolic,
        })
    }

    /// `Allele.create(s, isRef)` for a `&str`.
    pub fn from_str(s: &str, is_ref: bool) -> Result<Self, AlleleError> {
        Self::create(s.as_bytes(), is_ref)
    }

    /// `Allele.NO_CALL`.
    pub fn no_call() -> Self {
        Self {
            bases: Vec::new(),
            is_ref: false,
            is_no_call: true,
            is_symbolic: false,
        }
    }

    pub fn is_reference(&self) -> bool {
        self.is_ref
    }

    pub fn is_no_call(&self) -> bool {
        self.is_no_call
    }

    pub fn is_symbolic(&self) -> bool {
        self.is_symbolic
    }

    /// `getDisplayString`: the bases as written into REF and ALT.
    ///
    /// For a no-call this is the empty string, not `.`, because a no-call carries no bases.
    /// `getBaseString` is the one that substitutes `.`, and the encoder does not use it.
    pub fn display_string(&self) -> String {
        String::from_utf8_lossy(&self.bases).into_owned()
    }

    /// `getBaseString`: the segregating bases, with `.` substituted for a no-call.
    pub fn base_string(&self) -> String {
        if self.is_no_call {
            NO_CALL_STRING.to_string()
        } else if self.is_symbolic {
            // getBases() returns no bases for a symbolic allele, so the base string is empty
            // even though the display string is `<TAG>`.
            String::new()
        } else {
            self.display_string()
        }
    }

    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_bases_are_uppercased() {
        assert_eq!(
            Allele::from_str("acgt", true).unwrap().display_string(),
            "ACGT"
        );
        assert_eq!(Allele::from_str("a", false).unwrap().display_string(), "A");
    }

    /// The rule that catches a port out: symbolic bases keep their case, and "symbolic" is
    /// decided by the first and last byte, so it reaches further than `<...>`.
    #[test]
    fn symbolic_bases_keep_their_case() {
        assert_eq!(
            Allele::from_str("<del>", false).unwrap().display_string(),
            "<del>"
        );
        assert_eq!(
            Allele::from_str("a]chr2:456]", false)
                .unwrap()
                .display_string(),
            "a]chr2:456]",
            "a breakend contains ']' so it is symbolic, so it is not uppercased"
        );
        assert_eq!(
            Allele::from_str("at>", false).unwrap().display_string(),
            "at>",
            "ending in '>' is enough to be symbolic"
        );
        assert_eq!(
            Allele::from_str("at", false).unwrap().display_string(),
            "AT",
            "the same bases without the '>' are uppercased"
        );
    }

    #[test]
    fn a_no_call_has_no_bases() {
        let n = Allele::from_str(".", false).unwrap();
        assert!(n.is_no_call());
        assert_eq!(n.display_string(), "", "no bases, so nothing to display");
        assert_eq!(n.base_string(), ".");
        assert_eq!(n, Allele::no_call());
    }

    #[test]
    fn a_single_angle_bracket_is_not_symbolic_and_is_rejected() {
        assert!(!would_be_symbolic(b"<"));
        assert!(Allele::from_str("<", false).is_err());
    }

    #[test]
    fn the_reference_flag_is_part_of_identity() {
        let r = Allele::from_str("A", true).unwrap();
        let a = Allele::from_str("A", false).unwrap();
        assert_ne!(r, a, "same bases, different alleles");
    }

    #[test]
    fn what_cannot_be_the_reference() {
        assert_eq!(
            Allele::from_str(".", true),
            Err(AlleleError::NoCallAsReference)
        );
        assert_eq!(
            Allele::from_str("*", true),
            Err(AlleleError::SpanDelAsReference)
        );
        assert_eq!(
            Allele::from_str("<DEL>", true),
            Err(AlleleError::SymbolicAsReference)
        );
        assert_eq!(Allele::from_str("-", true), Err(AlleleError::NullAllele));
        assert_eq!(Allele::from_str("", false), Err(AlleleError::NullAllele));
    }

    #[test]
    fn a_span_del_is_allowed_as_an_alternate() {
        assert_eq!(Allele::from_str("*", false).unwrap().display_string(), "*");
    }

    #[test]
    fn iupac_ambiguity_codes_are_not_acceptable_bases() {
        assert!(Allele::from_str("M", false).is_err());
        assert!(Allele::from_str("ATM", false).is_err());
        assert!(Allele::from_str("ATN", false).is_ok(), "N is acceptable");
    }
}
