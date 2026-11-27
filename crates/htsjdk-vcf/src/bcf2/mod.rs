//! BCF2, read and written as htsjdk does it.
//!
//! Ported from `htsjdk.variant.bcf2` (`BCF2Codec`, `BCF2Decoder`, `BCF2GenotypeFieldDecoders`,
//! `BCF2LazyGenotypesDecoder`, `BCF2Type`, `BCF2Utils`, `BCFVersion`) and from
//! `htsjdk.variant.variantcontext.writer` (`BCF2Writer`, `BCF2Encoder`, `BCF2FieldEncoder`,
//! `BCF2FieldWriter`, `BCF2FieldWriterManager`, `IntGenotypeFieldAccessors`) at htsjdk 4.2.0. The
//! 5.0.0 sources Picard 3.5.0 ships differ from these only in formatting.
//!
//! # What htsjdk's BCF is
//!
//! **BCF 2.1, uncompressed.** The writer emits `BCF\2\1`, an INT32 header length, the VCF header
//! text with a terminating `\0`, and then per record an INT32 sites length, an INT32 genotypes
//! length and the two blocks, with no BGZF around any of it. The reader accepts major version 2
//! and minor version 1 and nothing else, so a `bcftools` file (2.2, BGZF) is refused.
//!
//! # The quirks that reach the bytes
//!
//! - **INFO is in hash order.** `buildInfo` iterates the attribute `HashMap`, and the order of a
//!   `HashMap` depends on every map the attributes were copied through on the way: see
//!   [`java_hash`] and [`decode`].
//! - **A record from a VCF is `fullyDecode`d first**, which turns each value into the type its
//!   header line declares, looking up the FORMAT line before the INFO line for an INFO key,
//!   refuses a value count the header disagrees with, and silently drops an INFO key named `FT`.
//! - **GQ is capped at 99** by `IntGenotypeFieldAccessors.GQAccessor`, so the BCF of a VCF with a
//!   GQ of 120 holds 99.
//! - **Genotype blocks are copied** from a BCF record nobody decoded into a writer whose header
//!   lists the same ID lines in the same order, so the FORMAT block of such a record is never
//!   re-encoded.
//! - **An integer INFO list with a missing element is a `NullPointerException`**:
//!   `determineIntegerType(List<Integer>)` unboxes every element, and `.` decodes to `null`.
//! - **Strings are measured in UTF-16 units and written in UTF-8 bytes**, so a non-ASCII value is
//!   cut short.

pub mod codec;
pub mod decode;
pub mod encoder;
pub mod java_hash;
pub mod types;
pub mod writer;

/// An exception as htsjdk throws it: the class Picard reports and its message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bcf2Error {
    pub class: &'static str,
    pub message: String,
}

/// A deferred genotype block's refusal, thrown on through a writer that reached it.
impl From<crate::genotypes_context::GenotypeError> for Bcf2Error {
    fn from(error: crate::genotypes_context::GenotypeError) -> Self {
        Self::new(error.class, error.message)
    }
}

impl Bcf2Error {
    pub fn new(class: &'static str, message: impl Into<String>) -> Self {
        Self {
            class,
            message: message.into(),
        }
    }

    pub fn tribble(message: impl Into<String>) -> Self {
        Self::new("htsjdk.tribble.TribbleException", message)
    }

    /// `TribbleException.InvalidHeader`, which prefixes its message.
    pub fn invalid_header(message: impl AsRef<str>) -> Self {
        Self::new(
            "htsjdk.tribble.TribbleException$InvalidHeader",
            format!(
                "Your input file has a malformed header: {}",
                message.as_ref()
            ),
        )
    }

    pub fn illegal_state(message: impl Into<String>) -> Self {
        Self::new("java.lang.IllegalStateException", message)
    }

    pub fn illegal_argument(message: impl Into<String>) -> Self {
        Self::new("java.lang.IllegalArgumentException", message)
    }

    pub fn null_pointer(message: impl Into<String>) -> Self {
        Self::new("java.lang.NullPointerException", message)
    }

    /// A failed cast, in the JDK 17 wording.
    pub fn class_cast(from: &str, to: &str) -> Self {
        Self::new(
            "java.lang.ClassCastException",
            format!(
                "class {from} cannot be cast to class {to} ({from} and {to} are in module \
                 java.base of loader 'bootstrap')"
            ),
        )
    }

    /// `ArrayList.get` out of range.
    pub fn index_out_of_bounds(index: i64, length: usize) -> Self {
        Self::new(
            "java.lang.IndexOutOfBoundsException",
            format!("Index {index} out of bounds for length {length}"),
        )
    }
}

impl std::fmt::Display for Bcf2Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.class, self.message)
    }
}

impl std::error::Error for Bcf2Error {}

/// `BCF2Utils.makeDictionary`: `PASS`, then the ID of every FILTER, INFO and FORMAT line in the
/// header's input order, each once.
pub fn make_dictionary<'a>(lines: impl IntoIterator<Item = &'a crate::HeaderLine>) -> Vec<String> {
    let mut dict = vec!["PASS".to_string()];
    for line in lines {
        let id = match line {
            crate::HeaderLine::Compound { key, id, .. } if key == "INFO" || key == "FORMAT" => id,
            crate::HeaderLine::Filter { id, .. } => id,
            // The parser makes a FILTER line a structured one keyed `FILTER`.
            crate::HeaderLine::Structured { key, fields } if key == "FILTER" => {
                match fields.iter().find(|(k, _)| k == "ID") {
                    Some((_, id)) => id,
                    None => continue,
                }
            }
            _ => continue,
        };
        if !dict.contains(id) {
            dict.push(id.clone());
        }
    }
    dict
}

/// The `(key, ID)` of every `VCFIDHeaderLine`, in order: what
/// `BCF2Utils.headerLinesAreOrderedConsistently` compares, since a line's class follows from its
/// key.
pub fn id_lines<'a>(
    lines: impl IntoIterator<Item = &'a crate::HeaderLine>,
) -> Vec<(String, String)> {
    lines
        .into_iter()
        .filter_map(|line| match line {
            crate::HeaderLine::Compound { key, id, .. } => Some((key.clone(), id.clone())),
            crate::HeaderLine::Filter { id, .. } => Some(("FILTER".to_string(), id.clone())),
            crate::HeaderLine::Contig { fields, .. }
            | crate::HeaderLine::Structured { fields, .. } => fields
                .iter()
                .find(|(k, _)| k == "ID")
                .map(|(_, v)| (line.key().to_string(), v.clone())),
            crate::HeaderLine::Unstructured { .. } => None,
        })
        .collect()
}

/// What `BCF2Codec.LazyData` carries beside the bytes: the header the genotype block was written
/// under, reduced to what `headerLinesAreOrderedConsistently` reads from it.
#[derive(Debug, PartialEq, Eq)]
pub struct LazyHeader {
    /// `getSampleNamesInOrder`: the samples, sorted.
    pub samples_in_order: Vec<String>,
    pub id_lines: Vec<(String, String)>,
}

/// `BCF2Codec.LazyData`: a genotype block as it was read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LazyData {
    pub header: std::sync::Arc<LazyHeader>,
    pub n_genotype_fields: i32,
    pub bytes: std::sync::Arc<[u8]>,
}

/// `String.compareTo`, which is UTF-16 order.
pub fn java_compare(a: &str, b: &str) -> std::cmp::Ordering {
    a.encode_utf16().cmp(b.encode_utf16())
}

/// `getSampleNamesInOrder`.
pub fn sample_names_in_order(samples: &[String]) -> Vec<String> {
    let mut names = samples.to_vec();
    names.sort_by(|a, b| java_compare(a, b));
    names
}
