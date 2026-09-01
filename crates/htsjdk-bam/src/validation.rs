//! `htsjdk.samtools.SamFileValidator`: what a SAM file is wrong about, in htsjdk's own words.
//!
//! The validator is htsjdk's; `ValidateSamFile` is the Picard tool that drives it and decides
//! whether to print every error or a histogram of their types. `picard-rs` ported the two together
//! because there was nowhere here to put the first half. This is that half: the checks, and the
//! `SAMValidationError.toString` rendering, with the mode selection and the summary histogram left
//! to the tool.
//!
//! An error carries a severity, a type, and optionally a record number and a read name, and prints
//! as `SEVERITY::TYPE:Record N, Read name R, message`. Two details of that line are the reference's
//! and not a formatter's choice: the record number is omitted when it is not positive, and two of
//! the read-group messages end with a trailing space before the read name is appended by the caller
//! ("A platform (PL) attribute was not found for read group "), which is where htsjdk's own string
//! concatenation stops.
//!
//! What is covered is what a SAM text file with an optional reference can be asked: the header
//! checks, `SAMRecord.isValid()`'s reference-free subset, the mate-pair checks
//! (`PairEndInfo.validateMates`), the sort-order check (`SAMSortOrderChecker`), and the `NM` tag,
//! by presence without a reference and by value with one. What is deferred is listed in
//! `picard-rs`'s tool module, which this was extracted from, and none of it is reachable from a
//! text file without a dictionary.

use std::collections::HashMap;

use crate::fasta::{read_fasta, FastaError};
use crate::header::SamHeader;
use crate::md_nm::calculate_md_and_nm;
use crate::record::BamRecord;
use crate::tag::Tag;

/// Why a validation run could not start: the reference FASTA did not parse.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValidationSetupError {
    Fasta(String),
}

impl From<FastaError> for ValidationSetupError {
    fn from(e: FastaError) -> Self {
        ValidationSetupError::Fasta(format!("{e:?}"))
    }
}

/// One `SAMValidationError`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidationError {
    /// `WARNING` or `ERROR`.
    pub severity: &'static str,
    /// The `SAMValidationError.Type` name, which is what a summary histogram is keyed on.
    pub error_type: &'static str,
    /// `recordNumber`, printed only when positive.
    pub record_number: Option<i64>,
    pub read_name: Option<String>,
    pub message: String,
}

impl std::fmt::Display for ValidationError {
    /// `SAMValidationError.toString`, for the case where there is no source file to name.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}::{}:", self.severity, self.error_type)?;
        if let Some(n) = self.record_number {
            if n > 0 {
                write!(f, "Record {n}, ")?;
            }
        }
        if let Some(name) = &self.read_name {
            write!(f, "Read name {name}, ")?;
        }
        write!(f, "{}", self.message)
    }
}

/// Collects the errors a run finds, in the order it finds them.
struct Report {
    errors: Vec<ValidationError>,
}

impl Report {
    fn new() -> Self {
        Report { errors: Vec::new() }
    }

    fn add(
        &mut self,
        severity: &str,
        ty: &str,
        record_number: Option<i64>,
        read_name: Option<&str>,
        message: &str,
    ) {
        // The severity and type are `&'static str` in the errors this module builds; the signature
        // takes `&str` so the call sites read as they do in htsjdk.
        let severity: &'static str = if severity == "WARNING" {
            WARNING
        } else {
            ERROR
        };
        self.errors.push(ValidationError {
            severity,
            error_type: TYPES
                .iter()
                .copied()
                .find(|known| *known == ty)
                .unwrap_or_else(|| panic!("unknown validation type {ty}")),
            record_number,
            read_name: read_name.map(str::to_string),
            message: message.to_string(),
        });
    }
}

/// Every `SAMValidationError.Type` this module can produce, so a type name is a fact rather than a
/// string literal that a typo could invent.
const TYPES: [&str; 31] = [
    "DUPLICATE_READ_GROUP_ID",
    "INVALID_CIGAR",
    "INVALID_FLAG_FIRST_OF_PAIR",
    "INVALID_FLAG_MATE_NEG_STRAND",
    "INVALID_FLAG_MATE_UNMAPPED",
    "INVALID_FLAG_NOT_PRIM_ALIGNMENT",
    "INVALID_FLAG_PROPER_PAIR",
    "INVALID_FLAG_READ_UNMAPPED",
    "INVALID_FLAG_SECOND_OF_PAIR",
    "INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT",
    "INVALID_MAPPING_QUALITY",
    "INVALID_PLATFORM_VALUE",
    "INVALID_TAG_NM",
    "INVALID_VERSION_NUMBER",
    "MATES_ARE_SAME_END",
    "MATE_NOT_FOUND",
    "MISMATCH_FLAG_MATE_NEG_STRAND",
    "MISMATCH_FLAG_MATE_UNMAPPED",
    "MISMATCH_MATE_ALIGNMENT_START",
    "MISMATCH_MATE_CIGAR_STRING",
    "MISMATCH_MATE_REF_INDEX",
    "MISSING_PLATFORM_VALUE",
    "MISSING_READ_GROUP",
    "MISSING_SEQUENCE_DICTIONARY",
    "MISSING_TAG_NM",
    "MISSING_VERSION_NUMBER",
    "QUALITY_NOT_STORED",
    "READ_GROUP_NOT_FOUND",
    "RECORD_MISSING_READ_GROUP",
    "RECORD_OUT_OF_ORDER",
    "REF_SEQ_TOO_LONG_FOR_BAI",
];

const READ_PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const READ_UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const READ_REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const FIRST_OF_PAIR: u16 = 0x40;
const SECOND_OF_PAIR: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;

/// `SAMFileHeader.ACCEPTABLE_VERSIONS`.
const ACCEPTABLE_VERSIONS: [&str; 5] = ["1.0", "1.3", "1.4", "1.5", "1.6"];

/// `SAMReadGroupRecord.PlatformValue` (compared case-insensitively, as htsjdk uppercases first).
const PLATFORM_VALUES: [&str; 14] = [
    "BGI",
    "CAPILLARY",
    "DNBSEQ",
    "ELEMENT",
    "HELICOS",
    "ILLUMINA",
    "IONTORRENT",
    "LS454",
    "ONT",
    "OTHER",
    "PACBIO",
    "SINGULAR",
    "SOLID",
    "ULTIMA",
];

/// `GenomicIndexUtil.BIN_GENOMIC_SPAN` (512 MiB): a reference longer than this cannot be BAI-indexed.
const BIN_GENOMIC_SPAN: i64 = 512 * 1024 * 1024;

const WARNING: &str = "WARNING";
const ERROR: &str = "ERROR";

/// `SamFileValidator.validateHeader`, restricted to the reference-free checks in scope.
fn validate_header(header: &SamHeader, rep: &mut Report) {
    // Version: getVersion() is the @HD VN attribute. (A missing VN is a parser-level error here and
    // is deferred; if present it must be one of the acceptable versions.)
    match header.attributes.get("VN") {
        None => rep.add(
            ERROR,
            "MISSING_VERSION_NUMBER",
            None,
            None,
            "Header has no version number",
        ),
        Some(v) if !ACCEPTABLE_VERSIONS.contains(&v) => rep.add(
            ERROR,
            "INVALID_VERSION_NUMBER",
            None,
            None,
            &format!(
                "Header version: {v} does not match any of the acceptable versions: {}",
                ACCEPTABLE_VERSIONS.join(", ")
            ),
        ),
        Some(_) => {}
    }

    // Sequence dictionary: an empty one only arms a warning that fires on the first mapped read.
    if !header.sequences.is_empty() {
        let long: Vec<&str> = header
            .sequences
            .iter()
            .filter(|s| s.length as i64 > BIN_GENOMIC_SPAN)
            .map(|s| s.name.as_str())
            .collect();
        if !long.is_empty() {
            rep.add(
                WARNING,
                "REF_SEQ_TOO_LONG_FOR_BAI",
                None,
                None,
                &format!(
                    "Reference sequences are too long for BAI indexing: {}",
                    long.join(", ")
                ),
            );
        }
    }

    if header.read_groups.is_empty() {
        rep.add(
            ERROR,
            "MISSING_READ_GROUP",
            None,
            None,
            "Read groups is empty",
        );
    }

    // Read groups: duplicate id, then missing / invalid platform.
    let mut seen: Vec<&str> = Vec::new();
    for rg in &header.read_groups {
        let id = rg.id.as_str();
        if seen.contains(&id) {
            rep.add(
                ERROR,
                "DUPLICATE_READ_GROUP_ID",
                None,
                None,
                &format!("Duplicate read group id: {id}"),
            );
        } else {
            seen.push(id);
        }

        match rg.attributes.get("PL") {
            None | Some("") => rep.add(
                ERROR,
                "MISSING_PLATFORM_VALUE",
                None,
                Some(id),
                "A platform (PL) attribute was not found for read group ",
            ),
            Some(pl) if !PLATFORM_VALUES.contains(&pl.to_ascii_uppercase().as_str()) => rep.add(
                ERROR,
                "INVALID_PLATFORM_VALUE",
                None,
                Some(id),
                &format!(
                    "The platform (PL) attribute ({pl}) + was not one of the valid values for read group "
                ),
            ),
            Some(_) => {}
        }
    }
}

/// `SamFileValidator.validateReadGroup`: the record's read group is unknown if it has no `RG` tag or
/// the tag's id is not in the header.
fn read_group_present(header: &SamHeader, rec: &BamRecord) -> bool {
    match rec.tags.get(Tag::new(b"RG")) {
        Some(crate::tag::TagValue::Str(id)) => header.read_groups.iter().any(|rg| rg.id == *id),
        _ => false,
    }
}

/// `SamFileValidator.PairEndInfo`: the per-read view kept while waiting to meet a read's mate,
/// carrying both the read's own fields and what the read asserts about its mate.
struct PairEndInfo {
    read_alignment_start: i32,
    read_reference_index: i32,
    read_neg_strand: bool,
    read_unmapped: bool,
    read_cigar: String,
    mate_alignment_start: i32,
    mate_reference_index: i32,
    mate_neg_strand: bool,
    mate_unmapped: bool,
    mate_cigar: Option<String>,
    first_of_pair: bool,
    record_number: i64,
}

impl PairEndInfo {
    fn new(rec: &BamRecord, record_number: i64) -> Self {
        let mate_cigar = match rec.tags.get(Tag::new(b"MC")) {
            Some(crate::tag::TagValue::Str(s)) => Some(s.clone()),
            _ => None,
        };
        PairEndInfo {
            read_alignment_start: rec.alignment_start,
            read_reference_index: rec.reference_index,
            read_neg_strand: rec.flags & READ_REVERSE != 0,
            read_unmapped: rec.flags & READ_UNMAPPED != 0,
            read_cigar: rec.cigar.to_text(),
            mate_alignment_start: rec.mate_alignment_start,
            mate_reference_index: rec.mate_reference_index,
            mate_neg_strand: rec.flags & MATE_REVERSE != 0,
            mate_unmapped: rec.flags & MATE_UNMAPPED != 0,
            mate_cigar,
            first_of_pair: rec.flags & FIRST_OF_PAIR != 0,
            record_number,
        }
    }
}

/// `PairEndInfo.validateMateFields(end1, end2)`: the mate fields `end1` asserts must agree with
/// `end2`'s own fields. All errors carry `end1`'s record number.
fn validate_mate_fields(end1: &PairEndInfo, end2: &PairEndInfo, read_name: &str, rep: &mut Report) {
    let rn = Some(end1.record_number);
    if end1.mate_alignment_start != end2.read_alignment_start {
        rep.add(
            ERROR,
            "MISMATCH_MATE_ALIGNMENT_START",
            rn,
            Some(read_name),
            "Mate alignment does not match alignment start of mate",
        );
    }
    if end1.mate_neg_strand != end2.read_neg_strand {
        rep.add(
            ERROR,
            "MISMATCH_FLAG_MATE_NEG_STRAND",
            rn,
            Some(read_name),
            "Mate negative strand flag does not match read negative strand flag of mate",
        );
    }
    if end1.mate_reference_index != end2.read_reference_index {
        rep.add(
            ERROR,
            "MISMATCH_MATE_REF_INDEX",
            rn,
            Some(read_name),
            "Mate reference index (MRNM) does not match reference index of mate",
        );
    }
    if end1.mate_unmapped != end2.read_unmapped {
        rep.add(
            ERROR,
            "MISMATCH_FLAG_MATE_UNMAPPED",
            rn,
            Some(read_name),
            "Mate unmapped flag does not match read unmapped flag of mate",
        );
    }
    if let Some(mc) = &end1.mate_cigar {
        if mc != &end2.read_cigar {
            rep.add(
                ERROR,
                "MISMATCH_MATE_CIGAR_STRING",
                rn,
                Some(read_name),
                "Mate CIGAR string does not match CIGAR string of mate",
            );
        }
    }
}

/// `PairEndInfo.validateMates`: both directions, then the both-marked-same-end check (reported once,
/// against the first-seen read's record number).
fn validate_mates(first: &PairEndInfo, second: &PairEndInfo, read_name: &str, rep: &mut Report) {
    validate_mate_fields(first, second, read_name, rep);
    validate_mate_fields(second, first, read_name, rep);
    if first.first_of_pair == second.first_of_pair {
        let which = if first.first_of_pair {
            "first"
        } else {
            "second"
        };
        rep.add(
            ERROR,
            "MATES_ARE_SAME_END",
            Some(first.record_number),
            Some(read_name),
            &format!("Both mates are marked as {which} of pair"),
        );
    }
}

/// `SAMRecord.isValid`, restricted to the reference-free, dictionary-independent flag / mapping /
/// CIGAR / read-group checks, emitted in htsjdk's own order. Every error carries the record number
/// (`SamFileValidator` calls `setRecordNumber` on each). The mate-reference checks for unpaired
/// reads, the paired branch's reference/position checks, `INVALID_INSERT_SIZE`, the mapped read's
/// empty-dictionary / missing-reference-name checks, and `isValidReferenceIndexAndPosition` are
/// deferred (each needs the mate reference, the insert-size bound, or the sequence dictionary).
fn is_valid_record(header: &SamHeader, rec: &BamRecord, record_number: i64, rep: &mut Report) {
    let rn = Some(record_number);
    let name = Some(rec.read_name.as_str());
    let paired = rec.flags & READ_PAIRED != 0;
    let unmapped = rec.flags & READ_UNMAPPED != 0;

    if !paired {
        if rec.flags & PROPER_PAIR != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_PROPER_PAIR",
                rn,
                name,
                "Proper pair flag should not be set for unpaired read.",
            );
        }
        if rec.flags & MATE_UNMAPPED != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_MATE_UNMAPPED",
                rn,
                name,
                "Mate unmapped flag should not be set for unpaired read.",
            );
        }
        if rec.flags & MATE_REVERSE != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_MATE_NEG_STRAND",
                rn,
                name,
                "Mate negative strand flag should not be set for unpaired read.",
            );
        }
        if rec.flags & FIRST_OF_PAIR != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_FIRST_OF_PAIR",
                rn,
                name,
                "First of pair flag should not be set for unpaired read.",
            );
        }
        if rec.flags & SECOND_OF_PAIR != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_SECOND_OF_PAIR",
                rn,
                name,
                "Second of pair flag should not be set for unpaired read.",
            );
        }
    }

    if unmapped {
        if rec.flags & SECONDARY != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_NOT_PRIM_ALIGNMENT",
                rn,
                name,
                "Secondary alignment flag should not be set for unmapped read.",
            );
        }
        if rec.flags & SUPPLEMENTARY != 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT",
                rn,
                name,
                "Supplementary alignment flag should not be set for unmapped read.",
            );
        }
        if rec.mapping_quality != 0 {
            rep.add(
                ERROR,
                "INVALID_MAPPING_QUALITY",
                rn,
                name,
                "MAPQ should be 0 for unmapped read.",
            );
        }
    } else {
        if rec.cigar.elements.is_empty() {
            // (MAPQ >= 256 is unreachable: the field is a single byte.)
            rep.add(
                ERROR,
                "INVALID_CIGAR",
                rn,
                name,
                "CIGAR should have > zero elements for mapped read.",
            );
        }
        // `!hasReferenceName()`: a record whose RNAME is `*` while its unmapped flag is clear. The
        // type name is the reference's and it is the confusing one -- INVALID_FLAG_READ_UNMAPPED on
        // a record that is not flagged unmapped -- because the error is that the flag disagrees
        // with the reference name rather than that the flag is set.
        if rec.reference_index < 0 {
            rep.add(
                ERROR,
                "INVALID_FLAG_READ_UNMAPPED",
                rn,
                name,
                "Mapped read should have valid reference name",
            );
        }
    }

    // The RG ID, when present, must resolve in the header.
    if let Some(crate::tag::TagValue::Str(id)) = rec.tags.get(Tag::new(b"RG")) {
        if !header.read_groups.iter().any(|rg| rg.id == *id) {
            rep.add(
                ERROR,
                "READ_GROUP_NOT_FOUND",
                rn,
                name,
                &format!("RG ID on SAMRecord not found in header: {id}"),
            );
        }
    }
}

/// The header sort orders that carry a comparator (`SortOrder.getComparatorInstance`). `unsorted`,
/// `unknown`, a missing `SO`, and `duplicate` (deferred) get no order check.
enum SortOrder {
    Coordinate,
    Queryname,
    Unchecked,
}

impl SortOrder {
    fn of(header: &SamHeader) -> Self {
        match header.attributes.get("SO") {
            Some("coordinate") => SortOrder::Coordinate,
            Some("queryname") => SortOrder::Queryname,
            _ => SortOrder::Unchecked,
        }
    }

    fn name(&self) -> &'static str {
        match self {
            SortOrder::Coordinate => "coordinate",
            SortOrder::Queryname => "queryname",
            SortOrder::Unchecked => "unsorted",
        }
    }

    /// `SAMRecordCoordinateComparator` / `SAMRecordQueryNameComparator` `fileOrderCompare`. `prev` is
    /// in order iff this is `<= 0`.
    fn file_order_compare(&self, prev: &BamRecord, rec: &BamRecord) -> i32 {
        match self {
            SortOrder::Coordinate => {
                let (r1, r2) = (prev.reference_index, rec.reference_index);
                if r1 == -1 {
                    return if r2 == -1 { 0 } else { 1 };
                }
                if r2 == -1 {
                    return -1;
                }
                if r1 != r2 {
                    return r1 - r2;
                }
                prev.alignment_start - rec.alignment_start
            }
            // compareReadNames is String.compareTo, i.e. UTF-16 code-unit order, which equals Rust's
            // byte order for the ASCII read names in practice.
            SortOrder::Queryname => match prev.read_name.cmp(&rec.read_name) {
                std::cmp::Ordering::Less => -1,
                std::cmp::Ordering::Equal => 0,
                std::cmp::Ordering::Greater => 1,
            },
            SortOrder::Unchecked => 0,
        }
    }
}

/// `SamFileValidator.validateSamFile`, over records that are already parsed.
///
/// Without a reference each mapped read's `NM` tag is checked for presence alone
/// (`MISSING_TAG_NM`); with one, its value is recomputed and compared (`INVALID_TAG_NM`), which is
/// exactly what htsjdk skips when no reference is given.
pub fn validate(
    header: &SamHeader,
    records: &[BamRecord],
    fasta: Option<&[u8]>,
) -> Result<Vec<ValidationError>, ValidationSetupError> {
    let mut rep = Report::new();

    // The reference bases by contig name, resolved once, for the `NM` value check.
    let contigs = match fasta {
        Some(bytes) => read_fasta(bytes)?,
        None => Vec::new(),
    };
    let ref_by_name: HashMap<&str, &[u8]> = contigs
        .iter()
        .map(|c| (c.name.as_str(), c.bases.as_slice()))
        .collect();
    let have_reference = fasta.is_some();

    validate_header(header, &mut rep);

    // Armed by an empty dictionary, disarmed by the first mapped read it reports on.
    let mut dict_empty_pending = header.sequences.is_empty();

    // Sort-order checking: the comparator from the header's SO, and the previous record seen.
    let sort_order = SortOrder::of(header);
    let mut prev_record: Option<&BamRecord> = None;

    // Reads awaiting their mate. Keyed, as htsjdk's coordinate-sorted map is, by (reference bucket,
    // read name): a read is stored under the reference index it claims for its mate, and matched
    // when a later read on that reference arrives. A linear vector keeps a deterministic order for
    // the leftover `MATE_NOT_FOUND` pass. (Cross-reference pairing and multi-leftover ordering are
    // deferred; the covered corpus keeps every pair on one reference with at most one leftover.)
    let mut pending: Vec<(i32, String, PairEndInfo)> = Vec::new();

    for (i, rec) in records.iter().enumerate() {
        let record_number = (i + 1) as i64;
        let unmapped = rec.flags & READ_UNMAPPED != 0;

        // isValid(): the per-record flag / mapping / CIGAR / read-group checks, emitted first.
        is_valid_record(header, rec, record_number, &mut rep);

        // validateMateFields: only for paired, primary reads. (The MC-as-valid-cigar check is
        // deferred; the corpus MC tags are valid cigars.)
        if rec.flags & READ_PAIRED != 0 && rec.flags & (SECONDARY | SUPPLEMENTARY) == 0 {
            let found = pending.iter().position(|(bucket, name, _)| {
                *bucket == rec.reference_index && *name == rec.read_name
            });
            if let Some(pos) = found {
                let (_, _, first) = pending.remove(pos);
                let second = PairEndInfo::new(rec, record_number);
                validate_mates(&first, &second, &rec.read_name, &mut rep);
            } else {
                pending.push((
                    rec.mate_reference_index,
                    rec.read_name.clone(),
                    PairEndInfo::new(rec, record_number),
                ));
            }
        }

        // validateSortOrder: compare against the previous record under the header's comparator.
        if let Some(prev) = prev_record {
            if sort_order.file_order_compare(prev, rec) > 0 {
                rep.add(
                    ERROR,
                    "RECORD_OUT_OF_ORDER",
                    Some(record_number),
                    Some(&rec.read_name),
                    &format!(
                        "The record is out of [{}] order, prior read name [{}], prior coodinates [{}:{}]",
                        sort_order.name(),
                        prev.read_name,
                        prev.reference_index,
                        prev.alignment_start,
                    ),
                );
            }
        }
        prev_record = Some(rec);

        // validateReadGroup
        if !read_group_present(header, rec) {
            rep.add(
                WARNING,
                "RECORD_MISSING_READ_GROUP",
                None,
                Some(&rec.read_name),
                "A record is missing a read group",
            );
        }

        // validateNmTag: the tag must be present (MISSING_TAG_NM) and, when a reference is given,
        // must equal the value recomputed from the reference (INVALID_TAG_NM).
        if !unmapped {
            match rec.tags.get(Tag::new(b"NM")) {
                None => rep.add(
                    WARNING,
                    "MISSING_TAG_NM",
                    Some(record_number),
                    Some(&rec.read_name),
                    "NM tag (nucleotide differences) is missing",
                ),
                Some(crate::tag::TagValue::Int(tag_nm)) if have_reference => {
                    let name = &header.sequences[rec.reference_index as usize].name;
                    if let Some(ref_bases) = ref_by_name.get(name.as_str()) {
                        let (_, actual) = calculate_md_and_nm(
                            rec.alignment_start,
                            &rec.cigar,
                            &rec.read_bases,
                            ref_bases,
                        );
                        if *tag_nm != actual as i64 {
                            rep.add(
                                ERROR,
                                "INVALID_TAG_NM",
                                Some(record_number),
                                Some(&rec.read_name),
                                &format!(
                                    "NM tag (nucleotide differences) in file [{tag_nm}] does not match reality [{actual}]"
                                ),
                            );
                        }
                    }
                }
                Some(_) => {}
            }
        }

        // Empty dictionary reported once, on the first mapped read.
        if dict_empty_pending && !unmapped {
            rep.add(
                ERROR,
                "MISSING_SEQUENCE_DICTIONARY",
                None,
                None,
                "Sequence dictionary is empty",
            );
            dict_empty_pending = false;
        }

        // QUAL == '*' (no stored qualities).
        if rec.base_qualities.is_empty() {
            rep.add(
                WARNING,
                "QUALITY_NOT_STORED",
                Some(record_number),
                Some(&rec.read_name),
                "QUAL field is set to * (unspecified quality scores), this is allowed by the SAM \
                 specification but many tools expect reads to include qualities ",
            );
        }
    }

    // validateUnmatchedPairs: reads marked paired whose mate never arrived.
    for (_, name, _) in &pending {
        rep.add(
            ERROR,
            "MATE_NOT_FOUND",
            None,
            Some(name),
            "Mate not found for paired read",
        );
    }

    // The caller decides what to do with them: `ValidateSamFile` prints one line each in VERBOSE
    // mode and a histogram of their types in SUMMARY mode, and prints `No errors found` for an
    // empty list. Both of those are the tool's, not the validator's.
    Ok(rep.errors)
}
