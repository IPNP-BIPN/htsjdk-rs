//! `GenotypesContext` and the laziness `LazyGenotypesContext` gives it.
//!
//! Ported from `htsjdk.variant.variantcontext.GenotypesContext` and
//! `htsjdk.variant.variantcontext.LazyGenotypesContext` at htsjdk 4.2.0, restricted to the one
//! behaviour that reaches the bytes: **a record nobody looked at is written from the file's own
//! text**.
//!
//! # Why this type exists at all
//!
//! `VCFCodec` does not parse the genotype columns of a line it decodes. It hands the substring to
//! a `LazyGenotypesContext` and returns, and `VCFEncoder` then has two paths:
//!
//! ```java
//! if (gc.isLazyWithData() && ((LazyGenotypesContext) gc).getUnparsedGenotypeData() instanceof String) {
//!     vcfOutput.append(VCFConstants.FIELD_SEPARATOR);
//!     vcfOutput.append(((LazyGenotypesContext) gc).getUnparsedGenotypeData().toString());
//! } else {
//!     ... calcVCFGenotypeKeys, the sort, the trailing-field trim ...
//! }
//! ```
//!
//! The verbatim branch copies the FORMAT column too, so an untouched record keeps **the file's own
//! key order**, while a touched one gets the recomputed order: `GT` first and the rest sorted. The
//! two differ on any file whose FORMAT is not already in that order, which is most of them:
//!
//! ```text
//! input                      GT:GQ:DP:XX  0/1:60:10:7
//! copied through             GT:GQ:DP:XX  0/1:60:10:7   <- the file's order, verbatim
//! genotypes touched          GT:DP:GQ:XX  0/1:10:60:7   <- recomputed and sorted
//! ```
//!
//! # Decoding is what drops the text, and it is not free to look
//!
//! `decode()` ends with `unparsedGenotypeData = null`, and every accessor that needs the list runs
//! it. So *reading* a genotype is enough to change what the next write produces: it is not the
//! mutation that recomputes the keys, it is the look. `size()` and `isEmpty()` are the exceptions,
//! answered from `nUnparsedGenotypes` without parsing, and they stay exceptions here.
//!
//! # What this port does differently, and why
//!
//! htsjdk parses on first access and can therefore throw from a getter. This port parses when the
//! file is read, keeps the text beside the parsed genotypes, and defers only the *drop*: the
//! refusal a malformed genotype produces stays where [`crate::reader::read_vcf`] has always put
//! it, at the read, while the byte-level behaviour above is reproduced exactly. The observable
//! difference between the two is which call reports a broken file, and a file that parses at all
//! is written identically either way.

use std::sync::atomic::{AtomicBool, Ordering};

use crate::variant::Genotype;

/// The genotypes of one record, with the file's own text when they came from a file.
///
/// Derefs to the genotype list, so a caller reads it the way it always did -- and that read is
/// what marks the context decoded, which is `LazyGenotypesContext.decode()` clearing its text.
#[derive(Debug)]
pub struct GenotypesContext {
    genotypes: Vec<Genotype>,
    /// `unparsedGenotypeData`, which is the FORMAT column and every sample column after it, tab
    /// separated, exactly as the line carried them.
    unparsed: Option<String>,
    /// `loaded`. Once set, the text is no longer the record's answer.
    decoded: AtomicBool,
}

impl GenotypesContext {
    /// A context a caller built, which has no file text and is therefore never lazy.
    pub fn new(genotypes: Vec<Genotype>) -> Self {
        Self {
            genotypes,
            unparsed: None,
            decoded: AtomicBool::new(true),
        }
    }

    /// A context read from a file: the parsed genotypes and the text they were parsed from.
    pub fn lazy(genotypes: Vec<Genotype>, unparsed: String) -> Self {
        Self {
            genotypes,
            unparsed: Some(unparsed),
            decoded: AtomicBool::new(false),
        }
    }

    /// `isLazyWithData()`: the text is still the record's answer, so the encoder writes it.
    pub fn unparsed(&self) -> Option<&str> {
        if self.decoded.load(Ordering::Relaxed) {
            None
        } else {
            self.unparsed.as_deref()
        }
    }

    /// `size()`, which htsjdk answers from the unparsed count rather than by parsing.
    ///
    /// It is an inherent method so that it wins over the `Deref`: asking how many genotypes there
    /// are must not be what decides the FORMAT order of the next write.
    pub fn len(&self) -> usize {
        self.genotypes.len()
    }

    /// `isEmpty()`, the other accessor that does not decode.
    pub fn is_empty(&self) -> bool {
        self.genotypes.is_empty()
    }

    /// `decode()`: the text stops being the answer, and every later write recomputes the keys.
    pub fn decode(&self) {
        self.decoded.store(true, Ordering::Relaxed);
    }

    /// The genotypes without marking the context decoded.
    ///
    /// For a caller that has to look at the list *as the encoder does*, which is what the verbatim
    /// branch itself does when it does not take it. Nothing outside this crate needs it.
    pub(crate) fn peek(&self) -> &[Genotype] {
        &self.genotypes
    }

    /// The list, taken out. Used where a caller rebuilds the record from its parts.
    pub fn into_vec(self) -> Vec<Genotype> {
        self.genotypes
    }
}

impl std::ops::Deref for GenotypesContext {
    type Target = Vec<Genotype>;

    /// Reading the list is `decode()`, and `decode()` drops the text. That is htsjdk's rule and
    /// not a convenience: `ensureSampleNameMap` and `ensureSampleOrdering` both call it, so any
    /// accessor that reaches a genotype has already changed what the next write produces.
    fn deref(&self) -> &Vec<Genotype> {
        self.decode();
        &self.genotypes
    }
}

impl std::ops::DerefMut for GenotypesContext {
    fn deref_mut(&mut self) -> &mut Vec<Genotype> {
        self.decode();
        &mut self.genotypes
    }
}

impl From<Vec<Genotype>> for GenotypesContext {
    fn from(genotypes: Vec<Genotype>) -> Self {
        Self::new(genotypes)
    }
}

impl Clone for GenotypesContext {
    fn clone(&self) -> Self {
        Self {
            genotypes: self.genotypes.clone(),
            unparsed: self.unparsed.clone(),
            decoded: AtomicBool::new(self.decoded.load(Ordering::Relaxed)),
        }
    }
}

/// Equality is over the genotypes alone.
///
/// Two records with the same genotypes are the same record whether or not one of them still
/// remembers the text it was read from, and a comparison must not be what decodes a context.
impl PartialEq for GenotypesContext {
    fn eq(&self, other: &Self) -> bool {
        self.genotypes == other.genotypes
    }
}

impl Default for GenotypesContext {
    fn default() -> Self {
        Self::new(Vec::new())
    }
}

impl<'a> IntoIterator for &'a GenotypesContext {
    type Item = &'a Genotype;
    type IntoIter = std::slice::Iter<'a, Genotype>;

    fn into_iter(self) -> Self::IntoIter {
        self.decode();
        self.genotypes.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genotypes() -> Vec<Genotype> {
        vec![Genotype::new("s0", Vec::new())]
    }

    #[test]
    fn a_built_context_is_never_lazy() {
        let context = GenotypesContext::new(genotypes());
        assert_eq!(context.unparsed(), None);
    }

    #[test]
    fn a_read_context_answers_with_its_text_until_something_looks() {
        let context = GenotypesContext::lazy(genotypes(), "GT\t0/1".to_string());
        assert_eq!(context.unparsed(), Some("GT\t0/1"));
        // `size()` does not decode, which is the optimisation htsjdk keeps for exactly this.
        assert_eq!(context.len(), 1);
        assert!(!context.is_empty());
        assert_eq!(context.unparsed(), Some("GT\t0/1"));
        // Reading the list does.
        let _ = context.first();
        assert_eq!(context.unparsed(), None);
    }

    #[test]
    fn a_clone_carries_the_state_it_was_cloned_in() {
        let context = GenotypesContext::lazy(genotypes(), "GT\t0/1".to_string());
        let before = context.clone();
        assert_eq!(before.unparsed(), Some("GT\t0/1"));
        let _ = context.iter();
        assert_eq!(context.clone().unparsed(), None);
    }
}
