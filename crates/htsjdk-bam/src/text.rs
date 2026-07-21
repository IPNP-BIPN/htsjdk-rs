//! The SAM text record writer.
//!
//! Ported from `htsjdk.samtools.SAMTextWriter.writeAlignment` and
//! `htsjdk.samtools.TextTagCodec.encode`.
//!
//! ## The counterpoint to decision 0008
//!
//! Decision 0008 records the care htsjdk takes choosing an integer tag's binary width: a ladder
//! that writes 200 as an unsigned byte and 300 as a *signed short*, where the obvious rule gets
//! it wrong. `TextTagCodec.encode` throws all of that away:
//!
//! ```java
//! case 'c': case 'C': case 's': case 'S': case 'I':
//!     tagType = 'i';
//! ```
//!
//! Every integer type collapses to `i` in the text form. So a BAM and the SAM it converts to
//! carry the same value and different type letters, and converting SAM back to BAM re-derives
//! the width from the value through the same ladder. The width is therefore **not information**
//! in the file; it is a function of the value, and the text form simply declines to store it.
//!
//! That is worth knowing before anyone tries to "round-trip" a BAM through SAM and compare
//! bytes: it works, but only because the ladder is deterministic.

use crate::record::BamRecord;
use crate::tag::{Tag, TagValue};

/// `SAMRecord.NO_ALIGNMENT_REFERENCE_NAME`.
pub const NO_ALIGNMENT_REFERENCE_NAME: &str = "*";
/// `SAMRecord.NULL_SEQUENCE_STRING`.
pub const NULL_SEQUENCE_STRING: &str = "*";
/// `SAMRecord.NULL_QUALS_STRING`.
pub const NULL_QUALS_STRING: &str = "*";

const FIELD_SEPARATOR: char = '\t';

/// `SAMUtils.phredToFastq`: Phred plus 33.
pub fn phred_to_fastq(quals: &[u8]) -> String {
    quals.iter().map(|&q| (q + 33) as char).collect()
}

/// `TextTagCodec.encode`, the text form of one tag.
///
/// Returns `None` for a value htsjdk refuses to write, matching the exception it throws.
pub fn encode_tag(tag: Tag, value: &TagValue) -> Option<String> {
    let name = tag.to_string();
    Some(match value {
        TagValue::Char(c) => format!("{name}:A:{}", *c as char),
        TagValue::Int(v) => {
            // The spec's range is [-2^31, 2^32), which is wider than either a signed or an
            // unsigned 32-bit integer alone.
            if *v < i32::MIN as i64 || *v > u32::MAX as i64 {
                return None;
            }
            format!("{name}:i:{v}")
        }
        // Java renders a float with `Float.toString`, which is not `{}` in Rust for every
        // value; the shared helper keeps the two together.
        TagValue::Float(f) => format!("{name}:f:{}", java_float_to_string(*f)),
        TagValue::Str(s) => format!("{name}:Z:{s}"),
        TagValue::Hex(bytes) => {
            let hex: String = bytes.iter().map(|b| format!("{b:02X}")).collect();
            format!("{name}:H:{hex}")
        }
        TagValue::ByteArray { values, unsigned } => {
            let letter = if *unsigned { 'C' } else { 'c' };
            let body: String = values
                .iter()
                .map(|v| {
                    // `widenToUnsigned` masks with 0xff before rendering, so an unsigned array
                    // prints 255 where the signed one prints -1 for the same byte.
                    if *unsigned {
                        format!(",{}", *v as u8)
                    } else {
                        format!(",{v}")
                    }
                })
                .collect();
            format!("{name}:B:{letter}{body}")
        }
        TagValue::ShortArray { values, unsigned } => {
            let letter = if *unsigned { 'S' } else { 's' };
            let body: String = values
                .iter()
                .map(|v| {
                    if *unsigned {
                        format!(",{}", *v as u16)
                    } else {
                        format!(",{v}")
                    }
                })
                .collect();
            format!("{name}:B:{letter}{body}")
        }
        TagValue::IntArray { values, unsigned } => {
            let letter = if *unsigned { 'I' } else { 'i' };
            let body: String = values
                .iter()
                .map(|v| {
                    if *unsigned {
                        format!(",{}", *v as u32)
                    } else {
                        format!(",{v}")
                    }
                })
                .collect();
            format!("{name}:B:{letter}{body}")
        }
        TagValue::FloatArray(values) => {
            let body: String = values
                .iter()
                .map(|v| format!(",{}", java_float_to_string(*v)))
                .collect();
            format!("{name}:B:f{body}")
        }
    })
}

/// `Float.toString`, which differs from Rust's `{}` in two visible ways.
///
/// Java always shows a decimal point (`1.0`, not `1`) and switches to scientific notation
/// outside `[1e-3, 1e7)`, with the exponent written as `E7` rather than `e7`.
fn java_float_to_string(f: f32) -> String {
    if f.is_nan() {
        return "NaN".to_string();
    }
    if f.is_infinite() {
        return if f < 0.0 { "-Infinity" } else { "Infinity" }.to_string();
    }
    let a = f.abs();
    if a != 0.0 && !(1e-3..1e7).contains(&a) {
        // Java's form is `d.dddEn`, with at least one fractional digit.
        let s = format!("{:e}", f);
        let (mantissa, exp) = s.split_once('e').unwrap();
        let mantissa = if mantissa.contains('.') {
            mantissa.to_string()
        } else {
            format!("{mantissa}.0")
        };
        return format!("{mantissa}E{exp}");
    }
    let s = format!("{f}");
    if s.contains('.') || s.contains('e') {
        s
    } else {
        format!("{s}.0")
    }
}

/// `SAMTextWriter.writeAlignment`, without the trailing newline.
///
/// `reference_name` and `mate_reference_name` are resolved by the caller from the header, since
/// a record carries indices rather than names.
pub fn write_alignment(
    rec: &BamRecord,
    reference_name: &str,
    mate_reference_name: &str,
) -> Option<String> {
    let mut out = String::new();
    out.push_str(&rec.read_name);
    out.push(FIELD_SEPARATOR);
    out.push_str(&rec.flags.to_string());
    out.push(FIELD_SEPARATOR);
    out.push_str(reference_name);
    out.push(FIELD_SEPARATOR);
    out.push_str(&rec.alignment_start.to_string());
    out.push(FIELD_SEPARATOR);
    out.push_str(&rec.mapping_quality.to_string());
    out.push(FIELD_SEPARATOR);
    out.push_str(&rec.cigar.to_text());
    out.push(FIELD_SEPARATOR);
    // htsjdk compares the two names by *reference identity*, with a comment noting the strings
    // are interned. Comparing by value gives the same answer for any record whose names came
    // from one header, which is every record read from a file.
    if reference_name == mate_reference_name && reference_name != NO_ALIGNMENT_REFERENCE_NAME {
        out.push('=');
    } else {
        out.push_str(mate_reference_name);
    }
    out.push(FIELD_SEPARATOR);
    out.push_str(&rec.mate_alignment_start.to_string());
    out.push(FIELD_SEPARATOR);
    out.push_str(&rec.inferred_insert_size.to_string());
    out.push(FIELD_SEPARATOR);
    if rec.read_bases.is_empty() {
        out.push_str(NULL_SEQUENCE_STRING);
    } else {
        out.extend(rec.read_bases.iter().map(|&b| b as char));
    }
    out.push(FIELD_SEPARATOR);
    if rec.base_qualities.is_empty() {
        out.push_str(NULL_QUALS_STRING);
    } else {
        out.push_str(&phred_to_fastq(&rec.base_qualities));
    }
    for (tag, value) in rec.tags.iter() {
        out.push(FIELD_SEPARATOR);
        out.push_str(&encode_tag(*tag, value)?);
    }
    Some(out)
}

/// Decodes bases the way a SAM reader must, which is **not** what a BAM decoder produces.
///
/// SAM text carries the bases verbatim, so lower case survives. A BAM cannot represent it: its
/// nibble encoding decodes to upper case only (decision 0008). Any code that folds case is
/// therefore reachable through this path and unreachable through the binary one.
pub fn sam_text_preserves_case(bases: &str) -> Vec<u8> {
    bases.bytes().collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cigar::{Cigar, CigarElement, Op};

    fn rec() -> BamRecord {
        BamRecord {
            read_name: "read1".into(),
            flags: 99,
            reference_index: 0,
            alignment_start: 100,
            mapping_quality: 60,
            cigar: Cigar::new(vec![CigarElement {
                length: 4,
                op: Op::M,
            }]),
            mate_reference_index: 0,
            mate_alignment_start: 300,
            inferred_insert_size: 250,
            read_bases: b"ACGT".to_vec(),
            base_qualities: vec![30, 31, 32, 33],
            tags: Default::default(),
        }
    }

    fn t(name: &str) -> Tag {
        Tag::new(name.as_bytes().try_into().unwrap())
    }

    #[test]
    fn the_eleven_mandatory_fields_are_tab_separated() {
        let line = write_alignment(&rec(), "chr1", "chr1").unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), 11);
        assert_eq!(f[0], "read1");
        assert_eq!(f[1], "99");
        assert_eq!(f[2], "chr1");
        assert_eq!(f[3], "100");
        assert_eq!(f[4], "60");
        assert_eq!(f[5], "4M");
        assert_eq!(f[6], "=", "a mate on the same reference is written as =");
        assert_eq!(f[9], "ACGT");
        assert_eq!(f[10], "?@AB", "Phred 30..33 at offset 33");
    }

    #[test]
    fn a_mate_on_another_reference_is_named_in_full() {
        assert_eq!(
            write_alignment(&rec(), "chr1", "chr2")
                .unwrap()
                .split('\t')
                .nth(6),
            Some("chr2")
        );
    }

    /// The `=` shorthand is suppressed when there is no reference at all, so an unplaced pair
    /// writes `*` twice rather than `*` and `=`.
    #[test]
    fn two_absent_references_are_not_abbreviated() {
        let line = write_alignment(&rec(), "*", "*").unwrap();
        assert_eq!(line.split('\t').nth(6), Some("*"));
    }

    #[test]
    fn absent_bases_and_qualities_are_stars() {
        let mut r = rec();
        r.read_bases = Vec::new();
        r.base_qualities = Vec::new();
        let line = write_alignment(&r, "chr1", "chr1").unwrap();
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!((f[9], f[10]), ("*", "*"));
    }

    /// The counterpoint to decision 0008: every integer width collapses to `i` in text.
    #[test]
    fn every_integer_tag_type_collapses_to_i() {
        for v in [-1i64, 0, 100, 200, 300, 40_000, 70_000, 3_000_000_000] {
            let mut r = rec();
            r.tags.insert(t("XI"), TagValue::Int(v));
            let line = write_alignment(&r, "chr1", "chr1").unwrap();
            assert!(
                line.ends_with(&format!("XI:i:{v}")),
                "value {v} must be written as :i:, got {line}"
            );
        }
    }

    /// The binary form would have chosen four different type letters for those same values.
    #[test]
    fn the_binary_form_would_have_distinguished_them() {
        use crate::tag::integer_type;
        let letters: Vec<char> = [200i64, 300, 40_000, 70_000]
            .iter()
            .map(|v| integer_type(*v).unwrap() as char)
            .collect();
        assert_eq!(letters, vec!['C', 's', 'S', 'i']);
    }

    #[test]
    fn an_integer_outside_the_spec_range_is_refused() {
        let mut r = rec();
        r.tags.insert(t("XI"), TagValue::Int(u32::MAX as i64 + 1));
        assert_eq!(write_alignment(&r, "chr1", "chr1"), None);
    }

    #[test]
    fn arrays_keep_their_element_type_letter() {
        let mut r = rec();
        r.tags.insert(
            t("XB"),
            TagValue::IntArray {
                values: vec![1, 2, 3],
                unsigned: false,
            },
        );
        assert!(write_alignment(&r, "chr1", "chr1")
            .unwrap()
            .ends_with("XB:B:i,1,2,3"));
    }

    /// An unsigned array is widened before rendering, so the same stored byte prints 255 in one
    /// form and -1 in the other.
    #[test]
    fn unsigned_arrays_are_widened_before_rendering() {
        let mk = |unsigned| {
            let mut r = rec();
            r.tags.insert(
                t("XB"),
                TagValue::ByteArray {
                    values: vec![-1, 0, 1],
                    unsigned,
                },
            );
            write_alignment(&r, "chr1", "chr1").unwrap()
        };
        assert!(mk(false).ends_with("XB:B:c,-1,0,1"));
        assert!(mk(true).ends_with("XB:B:C,255,0,1"));
    }

    /// `Float.toString` always shows a decimal point and uses `E` for the exponent.
    #[test]
    fn floats_are_rendered_the_java_way() {
        assert_eq!(java_float_to_string(1.0), "1.0");
        assert_eq!(java_float_to_string(-0.5), "-0.5");
        assert_eq!(java_float_to_string(f32::NAN), "NaN");
        assert_eq!(java_float_to_string(f32::INFINITY), "Infinity");
        assert_eq!(java_float_to_string(f32::NEG_INFINITY), "-Infinity");
        assert!(java_float_to_string(1e10).contains('E'));
        assert!(!java_float_to_string(1e10).contains('e'));
    }

    #[test]
    fn tags_are_written_in_the_same_order_as_the_binary_form() {
        let mut r = rec();
        for name in ["ZA", "AZ", "NM"] {
            r.tags.insert(t(name), TagValue::Int(1));
        }
        let line = write_alignment(&r, "chr1", "chr1").unwrap();
        let tags: Vec<&str> = line.split('\t').skip(11).collect();
        assert_eq!(tags, vec!["ZA:i:1", "NM:i:1", "AZ:i:1"]);
    }

    /// SAM text carries bases verbatim, so lower case survives here and cannot survive a BAM.
    #[test]
    fn sam_text_is_where_lower_case_bases_are_reachable() {
        assert_eq!(sam_text_preserves_case("acgt"), b"acgt");
        assert_eq!(
            crate::bases::compressed_bases_to_bytes(
                4,
                &crate::bases::bytes_to_compressed_bases(b"acgt").unwrap(),
                0
            ),
            b"ACGT",
            "the binary form uppercases, so the two paths genuinely differ"
        );
    }
}
