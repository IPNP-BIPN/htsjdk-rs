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
//! # When the genotypes are parsed: at the read, or on first access
//!
//! htsjdk parses on first access and so can throw from a getter. Both are available here
//! (decision 0048):
//!
//!  * a reader left at its default, [`crate::reader::GenotypeParsing::Eager`], parses when the file
//!    is read, keeps the text beside the parsed genotypes and defers only the *drop*. A malformed
//!    genotype is refused at the read, and no accessor can fail;
//!  * a reader asked for [`crate::reader::GenotypeParsing::Deferred`] keeps the text and parses it
//!    the first time something looks, as `LazyGenotypesContext` does. A record nothing looks at is
//!    never parsed, so a malformed one that is only copied is copied, as Picard copies it. The
//!    refusal comes from the first look instead: [`GenotypesContext::try_get`] and its siblings
//!    return it, worded with the reader's line counter and last line at that moment, which is what
//!    htsjdk's message carries. `Deref` panics with it, so a caller that opts in reads through
//!    `try_*`.
//!
//! A file whose samples are not sorted is parsed at the read either way, as `parseVCFLine` does.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, OnceLock};

use crate::variant::Genotype;

/// A genotype block refused when it was finally parsed: the exception htsjdk's getter throws.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenotypeError {
    /// The Java class, `htsjdk.tribble.TribbleException` or one of its nested classes.
    pub class: &'static str,
    pub message: String,
}

impl std::fmt::Display for GenotypeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.class, self.message)
    }
}

impl std::error::Error for GenotypeError {}

impl From<crate::record_parse::RecordError> for GenotypeError {
    fn from(error: crate::record_parse::RecordError) -> Self {
        Self {
            class: error.class(),
            message: error.message().to_string(),
        }
    }
}

impl From<crate::bcf2::Bcf2Error> for GenotypeError {
    fn from(error: crate::bcf2::Bcf2Error) -> Self {
        Self {
            class: error.class,
            message: error.message,
        }
    }
}

/// The genotypes of one record, with the file's own text when they came from a file.
///
/// Derefs to the genotype list, so a caller reads it the way it always did -- and that read is
/// what marks the context decoded, which is `LazyGenotypesContext.decode()` clearing its text. A
/// deferred context that fails to parse makes that `Deref` panic; [`Self::try_get`] returns the
/// refusal instead.
#[derive(Debug)]
pub struct GenotypesContext {
    genotypes: Genotypes,
    /// `unparsedGenotypeData`, which is the FORMAT column and every sample column after it, tab
    /// separated, exactly as the line carried them.
    unparsed: Option<String>,
    /// The sample columns that text is in the order of, which is the header of the file it came
    /// from. Shared between every record of one file, so it costs one allocation per file.
    ///
    /// The text is only an answer for a writer whose header has the SAME samples in the SAME
    /// order. This is the port's rule, not htsjdk's: `VCFEncoder` copies the text whatever the
    /// writer's sample order (decision 0048 measured it). It holds for what it was written for,
    /// `MergeVcfs` with inputs whose sample columns are ordered differently, because htsjdk decodes
    /// those at the read.
    samples: Option<Arc<[String]>>,
    /// `loaded`. Once set, the text is no longer the record's answer.
    decoded: AtomicBool,
    /// `BCF2Codec.LazyData`: the genotype block of a BCF record, which `BCF2Writer` copies into a
    /// BCF whose header agrees with the one it was read under. It is never VCF text, so the VCF
    /// encoder re-encodes such a record.
    bcf: Option<crate::bcf2::LazyData>,
}

/// The list, or what will produce it.
#[derive(Debug, Clone)]
enum Genotypes {
    Parsed(Vec<Genotype>),
    /// `LazyGenotypesContext` before its first `decode()`.
    Deferred(Box<Deferred>),
}

#[derive(Debug, Clone)]
struct Deferred {
    parser: Parser,
    /// `nUnparsedGenotypes`: the header's sample count, which `size()` answers before a parse.
    count: usize,
    parsed: OnceLock<Result<Vec<Genotype>, GenotypeError>>,
}

/// The parser a deferred block is handed to: `LazyVCFGenotypesParser` or
/// `BCF2LazyGenotypesDecoder`.
#[derive(Debug, Clone)]
enum Parser {
    Vcf(crate::genotype_parse::DeferredVcfGenotypes),
    Bcf(crate::bcf2::codec::DeferredBcfGenotypes),
}

impl GenotypesContext {
    /// A context a caller built, which has no file text and is therefore never lazy.
    pub fn new(genotypes: Vec<Genotype>) -> Self {
        Self {
            genotypes: Genotypes::Parsed(genotypes),
            unparsed: None,
            samples: None,
            decoded: AtomicBool::new(true),
            bcf: None,
        }
    }

    /// A context read from a file: the parsed genotypes, the text they were parsed from, and the
    /// sample order that text is in.
    pub fn lazy(genotypes: Vec<Genotype>, unparsed: String, samples: Arc<[String]>) -> Self {
        Self {
            genotypes: Genotypes::Parsed(genotypes),
            unparsed: Some(unparsed),
            samples: Some(samples),
            decoded: AtomicBool::new(false),
            bcf: None,
        }
    }

    /// A context read from a file and not parsed yet: the text, the sample order it is in, and
    /// what will parse it on first access. The header's sample count is what [`Self::len`]
    /// answers until then.
    pub fn deferred_vcf(
        unparsed: String,
        samples: Arc<[String]>,
        parser: crate::genotype_parse::DeferredVcfGenotypes,
    ) -> Self {
        Self {
            genotypes: Genotypes::Deferred(Box::new(Deferred {
                parser: Parser::Vcf(parser),
                count: samples.len(),
                parsed: OnceLock::new(),
            })),
            unparsed: Some(unparsed),
            samples: Some(samples),
            decoded: AtomicBool::new(false),
            bcf: None,
        }
    }

    /// A context read from a BCF file: the decoded genotypes and the block they were decoded from.
    pub fn lazy_bcf(genotypes: Vec<Genotype>, data: crate::bcf2::LazyData) -> Self {
        Self {
            genotypes: Genotypes::Parsed(genotypes),
            unparsed: None,
            samples: None,
            decoded: AtomicBool::new(false),
            bcf: Some(data),
        }
    }

    /// A BCF genotype block not decoded yet, `count` samples long.
    pub(crate) fn deferred_bcf(
        data: crate::bcf2::LazyData,
        parser: crate::bcf2::codec::DeferredBcfGenotypes,
        count: usize,
    ) -> Self {
        Self {
            genotypes: Genotypes::Deferred(Box::new(Deferred {
                parser: Parser::Bcf(parser),
                count,
                parsed: OnceLock::new(),
            })),
            unparsed: None,
            samples: None,
            decoded: AtomicBool::new(false),
            bcf: Some(data),
        }
    }

    /// The BCF block while the context is still lazy, which is `isLazyWithData()` with
    /// `BCF2Codec.LazyData` behind it.
    pub fn bcf_unparsed(&self) -> Option<&crate::bcf2::LazyData> {
        if self.decoded.load(Ordering::Relaxed) {
            None
        } else {
            self.bcf.as_ref()
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

    /// The text, but only for a writer whose samples are the ones it was written under.
    ///
    /// A record read from one file and written under another file's header has columns in the
    /// wrong order, and copying them would put one sample's genotype in another's column. That is
    /// not a hypothetical: `MergeVcfs` accepts inputs whose sample columns are ordered differently
    /// and writes them under the SORTED order, and the reference reorders every column while a
    /// verbatim copy does not.
    pub fn unparsed_for(&self, samples: &[String]) -> Option<&str> {
        let text = self.unparsed()?;
        match &self.samples {
            Some(own) if own.as_ref() == samples => Some(text),
            _ => None,
        }
    }

    /// `size()`, which htsjdk answers from the unparsed count rather than by parsing.
    ///
    /// It is an inherent method so that it wins over the `Deref`: asking how many genotypes there
    /// are must not be what decides the FORMAT order of the next write, nor what parses a deferred
    /// block.
    pub fn len(&self) -> usize {
        match &self.genotypes {
            Genotypes::Parsed(genotypes) => genotypes.len(),
            Genotypes::Deferred(deferred) => match deferred.parsed.get() {
                Some(Ok(genotypes)) => genotypes.len(),
                _ => deferred.count,
            },
        }
    }

    /// `isEmpty()`, the other accessor that does not decode.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// `decode()`: the text stops being the answer, and every later write recomputes the keys.
    pub fn decode(&self) {
        self.decoded.store(true, Ordering::Relaxed);
    }

    /// The parse's answer, run once, without marking the context decoded.
    fn parsed(&self) -> Result<&Vec<Genotype>, &GenotypeError> {
        match &self.genotypes {
            Genotypes::Parsed(genotypes) => Ok(genotypes),
            Genotypes::Deferred(deferred) => deferred
                .parsed
                .get_or_init(|| match &deferred.parser {
                    Parser::Vcf(parser) => {
                        parser.parse(self.unparsed.as_deref().unwrap_or_default())
                    }
                    Parser::Bcf(parser) => parser.parse(
                        self.bcf
                            .as_ref()
                            .expect("a deferred BCF context keeps its block"),
                    ),
                })
                .as_ref(),
        }
    }

    /// `getGenotypes()`: the list, parsing a deferred block the first time, and `decode()`. A
    /// refused parse leaves the context as it was, still lazy, as htsjdk's `decode()` does when it
    /// throws before `loaded = true`.
    pub fn try_get(&self) -> Result<&[Genotype], GenotypeError> {
        self.try_get_vec().map(Vec::as_slice)
    }

    fn try_get_vec(&self) -> Result<&Vec<Genotype>, GenotypeError> {
        let genotypes = self.parsed().map_err(Clone::clone)?;
        self.decode();
        Ok(genotypes)
    }

    /// [`Self::try_get`], for a caller that changes the list.
    pub fn try_get_mut(&mut self) -> Result<&mut Vec<Genotype>, GenotypeError> {
        if matches!(self.genotypes, Genotypes::Deferred(_)) {
            self.parsed().map_err(Clone::clone)?;
            let parsed = std::mem::replace(&mut self.genotypes, Genotypes::Parsed(Vec::new()));
            if let Genotypes::Deferred(deferred) = parsed {
                let genotypes = deferred.parsed.into_inner().expect("parsed just above");
                self.genotypes = Genotypes::Parsed(genotypes.expect("parsed without refusal"));
            }
        }
        self.decode();
        match &mut self.genotypes {
            Genotypes::Parsed(genotypes) => Ok(genotypes),
            Genotypes::Deferred(_) => unreachable!("parsed just above"),
        }
    }

    /// The list, taken out. Used where a caller rebuilds the record from its parts.
    pub fn try_into_vec(mut self) -> Result<Vec<Genotype>, GenotypeError> {
        self.try_get_mut()?;
        match self.genotypes {
            Genotypes::Parsed(genotypes) => Ok(genotypes),
            Genotypes::Deferred(_) => unreachable!("parsed by try_get_mut"),
        }
    }

    /// The genotypes without marking the context decoded, or the refusal of a deferred block.
    pub fn try_peek(&self) -> Result<&[Genotype], GenotypeError> {
        self.parsed().map(Vec::as_slice).map_err(Clone::clone)
    }

    /// The genotypes without marking the context decoded.
    ///
    /// For a caller that has to look at the list *as the encoder does*, which is what the verbatim
    /// branch itself does when it does not take it. Nothing outside this crate needs it. Panics on
    /// a deferred block that does not parse; the encoder asks [`Self::try_peek`] first.
    pub fn peek(&self) -> &[Genotype] {
        match self.parsed() {
            Ok(genotypes) => genotypes,
            Err(error) => refused(error),
        }
    }

    /// The list, taken out. Panics on a deferred block that does not parse.
    pub fn into_vec(self) -> Vec<Genotype> {
        self.try_into_vec().unwrap_or_else(|error| refused(&error))
    }
}

/// The panic a deferred block's refusal becomes behind an accessor that cannot return it.
#[cold]
fn refused(error: &GenotypeError) -> ! {
    panic!("genotypes refused on first access ({error}); read them with try_get to handle it")
}

impl std::ops::Deref for GenotypesContext {
    type Target = Vec<Genotype>;

    /// Reading the list is `decode()`, and `decode()` drops the text. That is htsjdk's rule and
    /// not a convenience: `ensureSampleNameMap` and `ensureSampleOrdering` both call it, so any
    /// accessor that reaches a genotype has already changed what the next write produces.
    fn deref(&self) -> &Vec<Genotype> {
        self.try_get_vec().unwrap_or_else(|error| refused(&error))
    }
}

impl std::ops::DerefMut for GenotypesContext {
    fn deref_mut(&mut self) -> &mut Vec<Genotype> {
        match self.try_get_mut() {
            Ok(genotypes) => genotypes,
            Err(error) => refused(&error),
        }
    }
}

impl From<Vec<Genotype>> for GenotypesContext {
    fn from(genotypes: Vec<Genotype>) -> Self {
        Self::new(genotypes)
    }
}

impl FromIterator<Genotype> for GenotypesContext {
    /// A context collected from genotypes a caller built, which is never lazy: there is no file
    /// text behind it. `GenotypesContext.create(ArrayList<Genotype>)` is the same constructor.
    fn from_iter<I: IntoIterator<Item = Genotype>>(iter: I) -> Self {
        Self::new(iter.into_iter().collect())
    }
}

impl IntoIterator for GenotypesContext {
    type Item = Genotype;
    type IntoIter = std::vec::IntoIter<Genotype>;

    fn into_iter(self) -> Self::IntoIter {
        self.into_vec().into_iter()
    }
}

impl Clone for GenotypesContext {
    fn clone(&self) -> Self {
        Self {
            genotypes: self.genotypes.clone(),
            unparsed: self.unparsed.clone(),
            samples: self.samples.clone(),
            decoded: AtomicBool::new(self.decoded.load(Ordering::Relaxed)),
            bcf: self.bcf.clone(),
        }
    }
}

/// Equality is over the genotypes alone.
///
/// Two records with the same genotypes are the same record whether or not one of them still
/// remembers the text it was read from, and a comparison must not be what decodes a context. It
/// may parse a deferred one; two refusals are equal when their exceptions are.
impl PartialEq for GenotypesContext {
    fn eq(&self, other: &Self) -> bool {
        self.parsed() == other.parsed()
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
        std::ops::Deref::deref(self).iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reader::{read_vcf_with, GenotypeParsing};

    /// One good record and one whose first sample has a value for no key, read deferred.
    fn deferred_records() -> Vec<crate::variant::VariantContext> {
        let text = "##fileformat=VCFv4.2\n\
                    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
                    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">\n\
                    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
                    chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT:GQ\t0/1:60\t1/1:70\n\
                    chr1\t200\t.\tC\tG\t50\tPASS\t.\tGT:GQ\t0/1:60:99\t1/1:70\n";
        read_vcf_with(text, GenotypeParsing::Deferred)
            .expect("the sites read")
            .records
    }

    #[test]
    fn a_deferred_context_answers_its_size_without_parsing() {
        let records = deferred_records();
        let bad = &records[1].genotypes;
        assert_eq!((bad.len(), bad.is_empty()), (2, false));
        assert!(matches!(&bad.genotypes, Genotypes::Deferred(d) if d.parsed.get().is_none()));
    }

    #[test]
    fn a_refused_parse_leaves_the_record_lazy() {
        let records = deferred_records();
        let bad = &records[1].genotypes;
        let error = bad.try_get().expect_err("three values, two keys");
        assert_eq!(error.class, "htsjdk.tribble.TribbleException");
        assert!(
            error.message.contains("too many keys for the sample s0"),
            "{error}"
        );
        // `decode()` threw before `loaded = true`: the text is still what a write copies.
        assert_eq!(bad.unparsed(), Some("GT:GQ\t0/1:60:99\t1/1:70"));
        assert_eq!(bad.len(), 2);
        // And the refusal is the same every time it is asked.
        assert_eq!(bad.try_peek().unwrap_err(), error);
        let mut owned = bad.clone();
        assert_eq!(owned.try_get_mut().unwrap_err(), error);
        assert_eq!(bad.clone().try_into_vec().unwrap_err(), error);
    }

    #[test]
    fn a_good_deferred_context_parses_once_and_then_is_decoded() {
        let records = deferred_records();
        let good = &records[0].genotypes;
        assert_eq!(good.unparsed(), Some("GT:GQ\t0/1:60\t1/1:70"));
        let parsed = good.try_get().expect("well formed").to_vec();
        assert_eq!(parsed.len(), 2);
        assert_eq!(good.unparsed(), None, "a look decodes");
        assert_eq!(good.clone().try_into_vec().expect("well formed"), parsed);
        let mut owned = good.clone();
        owned.try_get_mut().expect("well formed").pop();
        assert_eq!(owned.len(), 1);
        assert_eq!(*good, GenotypesContext::new(parsed));
    }

    #[test]
    #[should_panic(expected = "genotypes refused on first access")]
    fn deref_panics_on_a_refused_block() {
        let records = deferred_records();
        let _ = records[1].genotypes.first();
    }

    fn samples() -> Arc<[String]> {
        Arc::from(vec!["s0".to_string()])
    }

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
        let context = GenotypesContext::lazy(genotypes(), "GT\t0/1".to_string(), samples());
        assert_eq!(context.unparsed(), Some("GT\t0/1"));
        // `size()` does not decode, which is the optimisation htsjdk keeps for exactly this.
        assert_eq!(context.len(), 1);
        assert!(!context.is_empty());
        assert_eq!(context.unparsed(), Some("GT\t0/1"));
        // Reading the list does.
        let _ = context.first();
        assert_eq!(context.unparsed(), None);
    }

    /// The text is not an answer for a writer whose columns are in a different order.
    #[test]
    fn the_text_is_only_the_answer_under_its_own_sample_order() {
        let context = GenotypesContext::lazy(genotypes(), "GT\t0/1".to_string(), samples());
        assert_eq!(context.unparsed_for(&["s0".to_string()]), Some("GT\t0/1"));
        assert_eq!(context.unparsed_for(&["s1".to_string()]), None);
        assert_eq!(
            context.unparsed_for(&["s0".to_string(), "s1".to_string()]),
            None
        );
        // And asking did not decode it: the question is about the header, not the genotypes.
        assert_eq!(context.unparsed(), Some("GT\t0/1"));
    }

    #[test]
    fn a_collected_context_is_never_lazy() {
        let context: GenotypesContext = genotypes().into_iter().collect();
        assert_eq!(context.unparsed(), None);
        assert_eq!(context.len(), 1);
    }

    #[test]
    fn a_clone_carries_the_state_it_was_cloned_in() {
        let context = GenotypesContext::lazy(genotypes(), "GT\t0/1".to_string(), samples());
        let before = context.clone();
        assert_eq!(before.unparsed(), Some("GT\t0/1"));
        let _ = context.iter();
        assert_eq!(context.clone().unparsed(), None);
    }
}
