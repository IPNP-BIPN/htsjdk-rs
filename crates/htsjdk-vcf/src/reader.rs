//! Reading a whole VCF file: the header, then every line after it, through one stateful codec.
//!
//! Ported from `htsjdk.variant.vcf.VCFCodec`, `AbstractVCFCodec.parseHeaderFromLines` and
//! `htsjdk.variant.vcf.VCFHeader` at htsjdk 4.2.0, driven the way `AbstractFeatureReader` drives
//! them: `readActualHeader` consumes lines up to and including `#CHROM`, and every line after it
//! goes to `decode`.
//!
//! [`crate::header_parse::read_header_frame`] reads the frame, [`crate::record_parse::decode_line`]
//! reads one line, and [`crate::genotype_parse::parse_genotypes`] reads one genotype block. Each
//! has its own conformance suite. What this adds is the state carried **between** those calls, and
//! that state is where the surprises are, because none of it is visible in any single line.
//!
//! # The header you get back is not the header the file contains
//!
//! Two rewrites happen before a caller ever sees it, and neither is in the format:
//!
//!  * **the file's `##fileformat` line is deleted and a different one is put back.** The
//!    `VCFHeader` constructor calls `removeVCFVersionLines`, and `getMetaDataInInputOrder`
//!    prepends a synthesized line: `VCFv4.3` when the header's own version is 4.3 or later, and
//!    `VCFv4.2` otherwise. So a v4.0 file reads back claiming to be v4.2, and a file declaring the
//!    version twice reads back with neither of its own lines;
//!  * **standard INFO and FORMAT lines are replaced when they disagree with htsjdk's**, because
//!    `doOnTheFlyModifications` defaults to true. See [`crate::standard_header_lines`].
//!
//! And the repair rebuilds the header through the constructor that takes no version, re-attaching
//! it **only for 4.3 and later**. So below 4.3 the header has forgotten which version it is while
//! the codec still knows: [`VcfFile::codec_version`] and [`VcfFile::header_version`] are different
//! answers to the same question and a port that keeps one field cannot give both.
//!
//! # The line counter is shared and incremented in two places
//!
//! `lineNo` counts header lines, then one per record **parsed**. `decodeLine`'s column-count check
//! runs before `parseVCFLine`'s increment and `generateException` runs after it, so the same
//! malformed line reports two different numbers depending on which check refuses it: a short line
//! at file line 13 says "Line 12" and a bad INFO field on that same line says "line number 13".
//! A `#` line increments nothing at all, because `decodeLine` returns before reaching the counter.
//!
//! # A `#` line in the body is a silently dropped record
//!
//! `decodeLine` returns `null` for it rather than a record or a refusal. A reader that collects
//! decode results without checking loses records with nothing to show for it, which is why
//! [`VcfFile::skipped`] is a field rather than a filter.

use crate::genotype_parse::{parse_genotypes, GenotypeContext};
use crate::header::{HeaderLine, VcfHeader};
use crate::header_lines::{parse_meta_line, HeaderLineError};
use crate::header_parse::{read_header_frame, InvalidHeader, VcfVersion};
use crate::record_parse::{decode_line, RecordError};
use crate::variant::VariantContext;

/// What a whole-file read produced.
#[derive(Debug, Clone, PartialEq)]
pub struct VcfFile {
    pub header: VcfHeader,
    /// `codec.getVersion()`: the version the file declared, which the codec keeps.
    pub codec_version: VcfVersion,
    /// `header.getVCFHeaderVersion()`: `None` below 4.3, because the repair rebuilt the header
    /// through a constructor that takes no version and only re-attached it from 4.3 up.
    pub header_version: Option<VcfVersion>,
    pub records: Vec<VariantContext>,
    /// The zero-based index each dropped `#` line would have had among the records.
    pub skipped: Vec<usize>,
}

impl VcfFile {
    /// `VCFHeader.getMetaDataInInputOrder`, which prepends a `fileformat` line of its own.
    ///
    /// The prepended line is **not** the file's. It says 4.3 when the header's own version is 4.3
    /// or later and 4.2 in every other case, including for a file that declared 4.0.
    pub fn meta_data_in_input_order(&self) -> Vec<HeaderLine> {
        let version = match self.header_version {
            Some(v) if v.is_at_least(VcfVersion::Vcf4_3) => VcfVersion::Vcf4_3,
            _ => VcfVersion::Vcf4_2,
        };
        let mut lines = vec![HeaderLine::Unstructured {
            key: version.format_string().to_string(),
            value: version.version_string().to_string(),
        }];
        lines.extend(self.header.lines.iter().cloned());
        lines
    }
}

/// A read that failed, with everything it had managed to read first.
///
/// The partial records are kept because a file that fails at its last line is not the same file as
/// one that fails at its first, and a reader that returns nothing in both cases cannot say which
/// happened.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadFailure {
    pub error: ReadError,
    pub records: Vec<VariantContext>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadError {
    /// Anything `readActualHeader` or `parseHeaderFromLines` refused.
    Header(InvalidHeader),
    /// A typed header line that would not parse.
    HeaderLine(HeaderLineError),
    /// Anything a data line refused.
    Record(RecordError),
}

impl ReadError {
    /// The Java class name, as a dump reports it.
    pub fn class(&self) -> &'static str {
        match self {
            ReadError::Header(_) => "htsjdk.tribble.TribbleException$InvalidHeader",
            ReadError::HeaderLine(error) => error.class(),
            ReadError::Record(error) => error.class(),
        }
    }

    pub fn message(&self) -> String {
        match self {
            ReadError::Header(error) => error.message(),
            ReadError::HeaderLine(error) => error.message(),
            ReadError::Record(error) => error.message().to_string(),
        }
    }
}

/// The header, as the codec establishes it, and the counter it leaves behind.
struct EstablishedHeader {
    header: VcfHeader,
    codec_version: VcfVersion,
    header_version: Option<VcfVersion>,
    /// `lineNo` after `readActualHeader`: one per header line, `#CHROM` included.
    line_number: usize,
}

/// `VCFCodec.readActualHeader` followed by `parseHeaderFromLines`, `repairStandardHeaderLines` and
/// the `VCFHeader` constructor, in that order.
fn establish_header(text: &str) -> Result<EstablishedHeader, ReadError> {
    let frame = read_header_frame(text).map_err(ReadError::Header)?;

    let mut lines: Vec<HeaderLine> = Vec::new();
    let mut contig_index = 0;
    for line in &frame.meta_lines {
        match parse_meta_line(line, frame.version, contig_index) {
            Ok(parsed) => {
                if matches!(parsed, HeaderLine::Contig { .. }) {
                    contig_index += 1;
                }
                lines.push(parsed);
            }
            // The one silent drop: `parseHeaderFromLines`' fallback only builds a line when the
            // text contains an '=', so a `##` line without one is discarded. Every other failure
            // is a throw upstream and a refusal here.
            Err(HeaderLineError::IllegalArgument(message))
                if message.starts_with("no '=' in header line") => {}
            Err(other) => return Err(ReadError::HeaderLine(other)),
        }
    }

    // A `LinkedHashSet`, so a line repeated exactly collapses onto its first occurrence and one
    // repeated with a different value does not.
    let mut deduped: Vec<HeaderLine> = Vec::with_capacity(lines.len());
    for line in lines {
        if !deduped.iter().any(|kept| header_lines_equal(kept, &line)) {
            deduped.push(line);
        }
    }

    // `repairStandardHeaderLines`, then `removeVCFVersionLines` in the constructor it calls. The
    // repair runs first upstream, on a header that still holds the fileformat line; the order is
    // immaterial because a fileformat line is never a compound one.
    let repaired: Vec<HeaderLine> = deduped
        .iter()
        .map(crate::standard_header_lines::repair)
        .filter(|line| !is_format_string(line.key()))
        .collect();

    Ok(EstablishedHeader {
        header: VcfHeader {
            lines: repaired,
            samples: frame.samples,
        },
        codec_version: frame.version,
        // `repairStandardHeaderLines` propagates the version only from 4.3 up.
        header_version: if frame.version.is_at_least(VcfVersion::Vcf4_3) {
            Some(frame.version)
        } else {
            None
        },
        line_number: frame.meta_lines.len() + 1,
    })
}

/// `VCFHeaderVersion.isFormatString`.
fn is_format_string(key: &str) -> bool {
    key == "fileformat" || key == "format"
}

/// `VCFHeaderLine.equals` and its overrides, which is what the `LinkedHashSet` deduplicates by.
///
/// Compound lines compare on ID, count, type and description and **ignore any extra tags**, so two
/// `##INFO` lines differing only in a `Source` tag are one line to htsjdk. A dedup by rendered
/// string would keep both, which is why this is spelled out rather than left to the renderer.
fn header_lines_equal(a: &HeaderLine, b: &HeaderLine) -> bool {
    match (a, b) {
        (
            HeaderLine::Compound {
                key: ka,
                id: ia,
                number: na,
                line_type: ta,
                description: da,
                ..
            },
            HeaderLine::Compound {
                key: kb,
                id: ib,
                number: nb,
                line_type: tb,
                description: db,
                ..
            },
        ) => ka == kb && ia == ib && na == nb && ta == tb && da == db,
        // Everything else compares by key and rendered value, which for these variants is the
        // whole line.
        _ => a.key() == b.key() && a.render() == b.render(),
    }
}

/// Read a whole VCF file.
///
/// The genotypes are decoded here rather than deferred. `LazyGenotypesContext` defers them
/// upstream and the deferral is observable in *when* a failure appears, not in what the record
/// holds, so a reader that decodes eagerly answers the same questions in a different order.
pub fn read_vcf(text: &str) -> Result<VcfFile, ReadFailure> {
    let established = match establish_header(text) {
        Ok(established) => established,
        Err(error) => {
            return Err(ReadFailure {
                error,
                records: Vec::new(),
            })
        }
    };
    let EstablishedHeader {
        header,
        codec_version,
        header_version,
        mut line_number,
    } = established;

    let mut records = Vec::new();
    let mut skipped = Vec::new();

    // The body is every line after the `#CHROM` one. `read_header_frame` stops there without
    // saying where, so the count of header lines is what locates it, and that count is the same
    // number the codec left in `lineNo`.
    for line in text.lines().skip(line_number) {
        let decoded = match decode_line(line, &header, line_number, codec_version) {
            Ok(Some(decoded)) => decoded,
            // A `#` line: no record, no refusal, and no increment of the counter either.
            Ok(None) => {
                skipped.push(records.len());
                continue;
            }
            Err(error) => {
                return Err(ReadFailure {
                    error: ReadError::Record(error),
                    records,
                })
            }
        };
        // `parseVCFLine` incremented on entry, which is why every field-level message on this line
        // carries a number one higher than the column check would have reported for it.
        line_number += 1;

        let mut variant = decoded.variant;
        if let Some(block) = decoded.genotype_block {
            let context = GenotypeContext {
                site_parts: &split_site_parts(line),
                header: &header,
                version: codec_version,
                contig: &variant.contig,
                pos: variant.start,
                line_number,
            };
            match parse_genotypes(&block, &variant.alleles, &context) {
                Ok(genotypes) => variant.genotypes = genotypes,
                Err(error) => {
                    return Err(ReadFailure {
                        error: ReadError::Record(error),
                        records,
                    })
                }
            }
        }
        records.push(variant);
    }

    Ok(VcfFile {
        header,
        codec_version,
        header_version,
        records,
        skipped,
    })
}

/// The site columns as the splitter left them, which the genotype layer needs only because one of
/// its messages quotes them.
fn split_site_parts(line: &str) -> Vec<String> {
    crate::record_parse::split_condensed(
        line,
        '\t',
        crate::record_parse::NUM_STANDARD_FIELDS + 1,
        true,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "##fileformat=VCFv4.2\n\
                          ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
                          #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    #[test]
    fn a_file_is_its_header_then_one_record_per_line() {
        let file = read_vcf(&format!(
            "{HEADER}chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\nchr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n"
        ))
        .expect("the fixture reads");
        assert_eq!(file.records.len(), 2);
        assert_eq!(file.records[1].start, 200);
    }

    /// The dropped record has no other trace, which is the reason `skipped` exists.
    #[test]
    fn a_hash_line_in_the_body_is_dropped_without_a_refusal() {
        let file = read_vcf(&format!(
            "{HEADER}chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\n# a comment\nchr1\t200\t.\tC\tG\t50\tPASS\tDP=20\n"
        ))
        .expect("the fixture reads");
        assert_eq!(file.records.len(), 2);
        assert_eq!(file.skipped, vec![1]);
    }

    /// Three header lines, so the first data line is file line 4 and the column check reports 3.
    #[test]
    fn the_column_check_reports_the_line_before_the_one_that_failed() {
        let failure =
            read_vcf(&format!("{HEADER}chr1\t100\t.\tA\tT\t50\tPASS\n")).expect_err("too few");
        assert!(
            failure.error.message().starts_with("Line 3: "),
            "{}",
            failure.error.message()
        );
        assert!(failure.records.is_empty());
    }

    /// The same line, refused by a check that runs after the increment, reports 4.
    #[test]
    fn a_field_level_refusal_on_that_line_reports_the_line_itself() {
        let failure = read_vcf(&format!("{HEADER}chr1\t100\t.\tA\tT\t50\tPASS\tDP=1 0\n"))
            .expect_err("space");
        assert!(
            failure
                .error
                .message()
                .contains("at approximately line number 4:"),
            "{}",
            failure.error.message()
        );
    }

    #[test]
    fn the_records_read_before_a_refusal_are_kept() {
        let failure = read_vcf(&format!(
            "{HEADER}chr1\t100\t.\tA\tT\t50\tPASS\tDP=10\nchr1\t200\t.\tC\tG\t50\tPASS\n"
        ))
        .expect_err("the second line is short");
        assert_eq!(failure.records.len(), 1);
    }

    /// Below 4.3 the header has forgotten its version and the codec has not.
    #[test]
    fn the_two_versions_disagree_below_four_three() {
        let file = read_vcf(HEADER).expect("a header-only file reads");
        assert_eq!(file.codec_version, VcfVersion::Vcf4_2);
        assert_eq!(file.header_version, None);
    }

    /// A v4.0 file reads back declaring v4.2, because the line it declared was deleted.
    #[test]
    fn the_fileformat_line_handed_back_is_not_the_file_s() {
        let text = "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let file = read_vcf(text).expect("a v4.0 file reads");
        assert_eq!(file.codec_version, VcfVersion::Vcf4_0);
        assert_eq!(
            file.meta_data_in_input_order()[0].render(),
            "fileformat=VCFv4.2"
        );
        assert!(
            !file.header.lines.iter().any(|l| l.key() == "fileformat"),
            "the file's own line is removed from the stored metadata"
        );
    }

    /// The same data line, two files, two meanings.
    #[test]
    fn a_percent_escape_is_decoded_only_from_four_three() {
        let body = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=%341\n";
        let head = |version| {
            format!(
                "##fileformat={version}\n\
                 ##INFO=<ID=DP,Number=1,Type=String,Description=\"Depth\">\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            )
        };
        let at42 = read_vcf(&format!("{}{body}", head("VCFv4.2"))).expect("4.2 reads");
        let at43 = read_vcf(&format!("{}{body}", head("VCFv4.3"))).expect("4.3 reads");
        assert_eq!(at42.records[0].attributes[0].1.format().unwrap(), "%341");
        assert_eq!(at43.records[0].attributes[0].1.format().unwrap(), "41");
    }
}
