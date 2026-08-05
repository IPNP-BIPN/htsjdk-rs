//! Writing a VCF and its `.idx` in one pass, which is what every GATK tool that emits a VCF does.
//!
//! Ported from `htsjdk.variant.variantcontext.writer.IndexingVariantContextWriter`,
//! `VariantContextWriterBuilder` and `htsjdk.tribble.index.TribbleIndexCreator` at htsjdk 4.2.0.
//!
//! [`crate::vcf_file::write_vcf`] writes the file and
//! [`htsjdk_tribble::index_write`] builds an index from a feature list. What neither can see on
//! its own is the join, and the join is where the file's bytes and the index's numbers have to
//! agree. Four things decide that agreement and none is in either format.
//!
//! # Indexing is on by default
//!
//! `VariantContextWriterBuilder.DEFAULT_OPTIONS` is `EnumSet.of(Options.INDEX_ON_THE_FLY)`, so a
//! caller who asks for nothing gets an index — and a caller who asks for nothing *and* supplies no
//! dictionary gets an exception instead of a file. The default is the behaviour, not a convenience.
//!
//! # The recorded position is the one before the record
//!
//! `IndexingVariantContextWriter.add` hands the indexer `locationSource.getPosition()` and only
//! then does `VCFWriter.add` write the line, so a record's block starts where the record starts.
//! The final position, handed to `finalizeIndex`, is the whole file's length. Both are absolute in
//! the output stream, so **the header is counted**: the first record's position is the header's
//! length, and an index built from a feature list that forgot the header is uniformly off by it.
//!
//! # The sequence dictionary becomes properties, not a flag
//!
//! `setIndexSequenceDictionary` writes one `DICT:<contig>` property per sequence, in dictionary
//! order, and `flags` stays **zero**: `SEQUENCE_DICTIONARY_FLAG` is only read for version < 3.
//! Those properties are added to the creator, and `DynamicIndexCreator.finalizeIndex` copies its
//! own map into the chosen creator **before** appending the four statistics, so the `DICT:` entries
//! come first and the order is observable in the bytes.
//!
//! # The layout is not the caller's choice, and it is usually the tree
//!
//! The writer always uses a `DynamicIndexCreator` with `FOR_SEEK_TIME`, so which layout lands
//! beside a VCF depends on the variants in it. Measured on ordinary files: one record at position
//! 100 gets a **linear** index and two thousand records get an **interval tree**, as does a
//! header-only file, whose feature density is a division of zero by zero.

use htsjdk_tribble::index::{TribbleIndex, INTERVAL_TREE, LINEAR, MAGIC_NUMBER, VERSION};
use htsjdk_tribble::index_write::{BalanceApproach, BuiltIndex, DynamicIndexCreator, Feature};

use crate::encoder::EncodeError;
use crate::header::VcfHeader;
use crate::variant::VariantContext;
use crate::vcf_file::write_vcf;

/// One entry of the `SAMSequenceDictionary` the builder is handed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceEntry {
    pub name: String,
    pub length: i64,
}

/// A VCF and the index written beside it.
#[derive(Debug, Clone, PartialEq)]
pub struct IndexedVcf {
    pub text: String,
    pub index: TribbleIndex,
    /// The byte offset of each record in `text`, which is what the index recorded.
    pub record_positions: Vec<i64>,
}

/// What the writer refuses, with the class each refusal arrives as.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WriteError {
    /// `build()` with `INDEX_ON_THE_FLY` and no dictionary.
    IllegalArgument(String),
    /// Indexing to a stream. Documented as an `IllegalArgumentException` and measured as a
    /// **NullPointerException**: nothing checks, and `writeBasedOnFeaturePath` dereferences the
    /// null path at `close()` — after the VCF itself has been written.
    NullPointer(String),
    Encode(EncodeError),
}

impl WriteError {
    pub fn class(&self) -> &'static str {
        match self {
            WriteError::IllegalArgument(_) => "java.lang.IllegalArgumentException",
            WriteError::NullPointer(_) => "java.lang.NullPointerException",
            WriteError::Encode(_) => "htsjdk.tribble.TribbleException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            WriteError::IllegalArgument(message) | WriteError::NullPointer(message) => {
                message.clone()
            }
            WriteError::Encode(_) => "the record could not be encoded".to_string(),
        }
    }
}

/// `VariantContextWriterBuilder.build()` with the default options, then `writeHeader`, `add` per
/// record, and `close`.
///
/// `indexed_path` is the URI the index header records, which upstream comes from
/// `location.toAbsolutePath().toUri()`. It is a parameter rather than derived because this port
/// does not open files: the caller owns where the bytes go.
///
/// `dictionary` being `None` is the refusal, not a mode: indexing is on by default, so there is no
/// way to reach the writer without one except by turning indexing off, which is
/// [`crate::vcf_file::write_vcf`].
pub fn write_vcf_indexed(
    header: &VcfHeader,
    records: &[VariantContext],
    dictionary: Option<&[SequenceEntry]>,
    indexed_path: &str,
) -> Result<IndexedVcf, WriteError> {
    let Some(dictionary) = dictionary else {
        return Err(WriteError::IllegalArgument(
            "A reference dictionary is required for creating Tribble indices on the fly"
                .to_string(),
        ));
    };
    if indexed_path.is_empty() {
        // `writerName` returns the stream's `toString()` when there is no path, and the null path
        // survives all the way to `close()`, where `writeBasedOnFeaturePath` dereferences it. The
        // message is the JVM's helpful-NPE text, which is why it is measured rather than written.
        return Err(WriteError::NullPointer(
            "Cannot invoke \"java.nio.file.Path.toAbsolutePath()\" because \"featurePath\" is null"
                .to_string(),
        ));
    }

    let text = write_vcf(header, records).map_err(WriteError::Encode)?;
    let header_text = header.write();

    // The positions the indexer was handed: the offset of each record's first byte, absolute in
    // the output stream and therefore counting the header.
    let mut creator = DynamicIndexCreator::new(BalanceApproach::ForSeekTime);
    let mut record_positions = Vec::with_capacity(records.len());
    let mut at = header_text.len() as i64;
    for (record, line) in records.iter().zip(text[header_text.len()..].lines()) {
        record_positions.push(at);
        creator.add_feature(
            &Feature {
                contig: record.contig.clone(),
                start: record.start as i32,
                end: record.stop as i32,
            },
            at,
        );
        at += line.len() as i64 + 1;
    }

    // `setIndexSequenceDictionary` runs at `close()`, on the creator, so its properties are in the
    // creator's own map before `finalizeIndex` appends the statistics.
    let mut properties: Vec<(String, String)> = dictionary
        .iter()
        .map(|entry| (format!("DICT:{}", entry.name), entry.length.to_string()))
        .collect();
    properties.extend(creator.properties());

    let built = creator
        .finalize(text.len() as i64)
        // The linear creator's `finalFilePosition == 0` guard is unreachable here: the header is
        // always written, so the final position is at least its length.
        .expect("a VCF always has a header, so the final position is never zero");

    let (index_type, contigs, interval_contigs) = match built {
        BuiltIndex::Linear(contigs) => (LINEAR, contigs, Vec::new()),
        BuiltIndex::IntervalTree(intervals) => (INTERVAL_TREE, Vec::new(), intervals),
    };

    Ok(IndexedVcf {
        text,
        index: TribbleIndex {
            index_type,
            version: VERSION,
            indexed_path: indexed_path.to_string(),
            indexed_file_size: 0,
            indexed_file_timestamp: 0,
            indexed_file_md5: String::new(),
            // Zero, not `SEQUENCE_DICTIONARY_FLAG`: from version 3 the dictionary is properties.
            flags: 0,
            properties,
            contigs,
            interval_contigs,
        },
        record_positions,
    })
}

/// `AbstractIndex.MAGIC_NUMBER`, re-exported so a caller checking a written index does not have to
/// reach across two crates for it.
pub const INDEX_MAGIC_NUMBER: i32 = MAGIC_NUMBER;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele::Allele;
    use crate::header::{Cardinality, HeaderLine, LineType};
    use crate::variant::Value;

    fn header() -> VcfHeader {
        let mut h = VcfHeader::new();
        h.lines.push(HeaderLine::info(
            "DP",
            Cardinality::Fixed(1),
            LineType::Integer,
            "Depth",
        ));
        h.lines.push(HeaderLine::contig("chr1", 100000, 0));
        h.lines.push(HeaderLine::contig("chr2", 200000, 1));
        h
    }

    fn record(contig: &str, start: i64) -> VariantContext {
        let mut vc = VariantContext::new(
            contig,
            start,
            vec![
                Allele::from_str("A", true).unwrap(),
                Allele::from_str("T", false).unwrap(),
            ],
        );
        vc.attributes = vec![("DP".to_string(), Value::Str("10".to_string()))];
        vc
    }

    fn dictionary() -> Vec<SequenceEntry> {
        vec![
            SequenceEntry {
                name: "chr1".into(),
                length: 100000,
            },
            SequenceEntry {
                name: "chr2".into(),
                length: 200000,
            },
        ]
    }

    #[test]
    fn the_first_record_sits_at_the_header_s_length() {
        let h = header();
        let written = write_vcf_indexed(
            &h,
            &[record("chr1", 100)],
            Some(&dictionary()),
            "file:///x.vcf",
        )
        .expect("the fixture writes");
        assert_eq!(written.record_positions[0], h.write().len() as i64);
    }

    #[test]
    fn the_dictionary_properties_come_before_the_statistics() {
        let written = write_vcf_indexed(
            &header(),
            &[record("chr1", 100)],
            Some(&dictionary()),
            "file:///x.vcf",
        )
        .expect("the fixture writes");
        let keys: Vec<&str> = written
            .index
            .properties
            .iter()
            .map(|(key, _)| key.as_str())
            .collect();
        assert_eq!(
            keys,
            vec![
                "DICT:chr1",
                "DICT:chr2",
                "FEATURE_LENGTH_MEAN",
                "FEATURE_LENGTH_STD_DEV",
                "MEAN_FEATURE_VARIANCE",
                "FEATURE_COUNT"
            ]
        );
        assert_eq!(written.index.flags, 0);
    }

    /// Indexing is on by default, so no dictionary is a refusal rather than an unindexed file.
    #[test]
    fn no_dictionary_is_refused_before_anything_is_written() {
        let error = write_vcf_indexed(&header(), &[], None, "file:///x.vcf").expect_err("no dict");
        assert_eq!(error.class(), "java.lang.IllegalArgumentException");
        assert!(error
            .message()
            .starts_with("A reference dictionary is required"));
    }

    /// One record close to the origin scores below the tree; the same writer over many records
    /// does not.
    #[test]
    fn the_layout_depends_on_the_variants_and_not_on_the_caller() {
        let one = write_vcf_indexed(
            &header(),
            &[record("chr1", 100)],
            Some(&dictionary()),
            "file:///x.vcf",
        )
        .expect("writes");
        assert_eq!(one.index.index_type, LINEAR);

        let many: Vec<VariantContext> = (0..2000).map(|i| record("chr1", 100 + i * 5)).collect();
        let lots = write_vcf_indexed(&header(), &many, Some(&dictionary()), "file:///x.vcf")
            .expect("writes");
        assert_eq!(lots.index.index_type, INTERVAL_TREE);
    }

    /// A header-only file still gets an index, and its density is zero divided by zero.
    #[test]
    fn a_header_only_file_is_indexed_too() {
        let written = write_vcf_indexed(&header(), &[], Some(&dictionary()), "file:///x.vcf")
            .expect("writes");
        assert_eq!(written.index.index_type, INTERVAL_TREE);
        assert!(written.index.contigs.is_empty());
        assert!(written.index.interval_contigs.is_empty());
    }
}
