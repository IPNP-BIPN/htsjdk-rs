//! The SAM text record parser.
//!
//! Ported from `htsjdk.samtools.SAMLineParser.parseLine` and
//! `htsjdk.samtools.TextTagCodec.decode`, the inverse of [`crate::text`].
//!
//! This is where the width information that [`crate::text`] throws away is **re-derived**. A
//! text tag says `i` for every integer; the parser reads the value and the binary encoder then
//! runs `BinaryTagCodec`'s promotion ladder over it again (decision 0008). So a
//! BAM → SAM → BAM round trip reproduces the original bytes, and it does so because the ladder
//! is a deterministic function of the value, not because the width survived the trip.
//!
//! It is also, per decision 0002, the only path on which a lower-case base is reachable: BAM
//! stores bases as nibbles that decode to upper case, and SAM text carries them verbatim.

use crate::cigar::{Cigar, CigarElement, Op};
use crate::record::BamRecord;
use crate::tag::{Tag, TagValue, Tags};

/// `SAMLineParser.NUM_REQUIRED_FIELDS`.
pub const NUM_REQUIRED_FIELDS: usize = 11;

/// `htsjdk.samtools.ValidationStringency`.
///
/// This exists because htsjdk's **writer does not enforce what its reader checks**. A record
/// with no reference name but a non-zero MAPQ is written without complaint and then rejected by
/// `SAMLineParser` at the default stringency, which is `STRICT`. Found by feeding htsjdk's own
/// SAM output back to this parser.
///
/// So a port that is unconditionally strict cannot read every file htsjdk produces, and a port
/// that is unconditionally lenient cannot reproduce htsjdk's default behaviour. Both are needed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ValidationStringency {
    /// `ValidationStringency.DEFAULT_STRINGENCY`.
    #[default]
    Strict,
    /// Logs and continues. The port returns the record and drops the message, since there is
    /// no logger here to write it to.
    Lenient,
    /// Neither throws nor logs.
    Silent,
}

impl ValidationStringency {
    fn rejects(self) -> bool {
        matches!(self, ValidationStringency::Strict)
    }
}

/// Why a SAM line was rejected. Messages mirror htsjdk's wording.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseError {
    NotEnoughFields(usize),
    /// `"Empty field at position N (zero-based)"`.
    EmptyField(usize),
    BadInteger {
        field: &'static str,
        value: String,
    },
    /// A mandatory-field combination htsjdk rejects.
    Inconsistent(String),
    InvalidReadBase(char),
    BadCigar(String),
    BadTag(String),
}

/// `SAMLineParser.isValidReadBase`.
///
/// The accepted alphabet is the IUPAC codes in both cases, plus `.` and `=`. Note what is
/// **not** here: `n` and `N` are, but nothing outside the table is, so a stray character is a
/// parse error rather than being folded into `N` the way [`crate::bases`] folds `.`.
pub fn is_valid_read_base(base: u8) -> bool {
    matches!(
        base,
        b'a' | b'c'
            | b'm'
            | b'g'
            | b'r'
            | b's'
            | b'v'
            | b't'
            | b'w'
            | b'y'
            | b'h'
            | b'k'
            | b'd'
            | b'b'
            | b'n'
            | b'A'
            | b'C'
            | b'M'
            | b'G'
            | b'R'
            | b'S'
            | b'V'
            | b'T'
            | b'W'
            | b'Y'
            | b'H'
            | b'K'
            | b'D'
            | b'B'
            | b'N'
            | b'.'
            | b'='
    )
}

/// `SAMUtils.fastqToPhred`: subtract 33 from each character.
pub fn fastq_to_phred(s: &str) -> Vec<u8> {
    s.bytes().map(|b| b.wrapping_sub(33)).collect()
}

/// `TextCigarCodec.decode`.
pub fn parse_cigar(text: &str) -> Result<Cigar, ParseError> {
    if text == "*" {
        return Ok(Cigar::default());
    }
    let mut elements = Vec::new();
    let mut len: u32 = 0;
    let mut saw_digit = false;
    for b in text.bytes() {
        if b.is_ascii_digit() {
            len = len
                .checked_mul(10)
                .and_then(|v| v.checked_add((b - b'0') as u32))
                .ok_or_else(|| ParseError::BadCigar(text.to_string()))?;
            saw_digit = true;
        } else {
            let op = match b {
                b'M' => Op::M,
                b'I' => Op::I,
                b'D' => Op::D,
                b'N' => Op::N,
                b'S' => Op::S,
                b'H' => Op::H,
                b'P' => Op::P,
                b'=' => Op::Eq,
                b'X' => Op::X,
                _ => return Err(ParseError::BadCigar(text.to_string())),
            };
            if !saw_digit {
                return Err(ParseError::BadCigar(text.to_string()));
            }
            elements.push(CigarElement { length: len, op });
            len = 0;
            saw_digit = false;
        }
    }
    if saw_digit {
        return Err(ParseError::BadCigar(text.to_string()));
    }
    Ok(Cigar::new(elements))
}

/// `TextTagCodec.decode` for one `TAG:TYPE:VALUE` field.
///
/// Returns `None` for a `B` array whose element type has no value list, matching htsjdk's
/// handling of a bare `B:i` with nothing after it.
pub fn parse_tag(field: &str) -> Result<(Tag, TagValue), ParseError> {
    let bad = || ParseError::BadTag(field.to_string());
    let mut parts = field.splitn(3, ':');
    let name = parts.next().ok_or_else(bad)?;
    let ty = parts.next().ok_or_else(bad)?;
    // htsjdk allows a two-field tag, treating the value as empty.
    let value = parts.next().unwrap_or("");
    if name.len() != 2 {
        return Err(bad());
    }
    let tag = Tag::new(name.as_bytes().try_into().map_err(|_| bad())?);

    let v = match ty {
        "Z" => TagValue::Str(value.to_string()),
        "A" => {
            if value.len() != 1 {
                return Err(ParseError::BadTag(
                    "Tag of type A should have a single-character value".into(),
                ));
            }
            TagValue::Char(value.as_bytes()[0])
        }
        "i" => {
            let n: i64 = value.parse().map_err(|_| {
                ParseError::BadTag("Tag of type i should have signed decimal value".into())
            })?;
            // htsjdk accepts [-2^31, 2^32), which is wider than either 32-bit type alone, and
            // returns an Integer or a Long accordingly. Both land in the same Rust variant,
            // and the binary encoder re-derives the width from the value.
            if n < i32::MIN as i64 || n > u32::MAX as i64 {
                return Err(ParseError::BadTag(format!(
                    "Integer is out of range for both a 32-bit signed and unsigned integer: {value}"
                )));
            }
            TagValue::Int(n)
        }
        "f" => TagValue::Float(value.parse().map_err(|_| {
            ParseError::BadTag(
                "Tag of type f should have single-precision floating point value".into(),
            )
        })?),
        "H" => {
            if !value.len().is_multiple_of(2) {
                return Err(ParseError::BadTag(
                    "Tag of type H should have valid hex string with even number of digits".into(),
                ));
            }
            let bytes = (0..value.len() / 2)
                .map(|i| u8::from_str_radix(&value[i * 2..i * 2 + 2], 16))
                .collect::<Result<Vec<_>, _>>()
                .map_err(|_| ParseError::BadTag("bad hex".into()))?;
            TagValue::Hex(bytes)
        }
        "B" => return parse_array_tag(tag, value).map(|v| (tag, v)),
        _ => return Err(ParseError::BadTag(format!("Unrecognized tag type: {ty}"))),
    };
    Ok((tag, v))
}

/// `TextTagCodec.covertStringArrayToObject`.
///
/// The method name is misspelled in htsjdk. Noted rather than corrected, because a reader
/// grepping for it in the reference will find the typo and not the fix.
fn parse_array_tag(_tag: Tag, value: &str) -> Result<TagValue, ParseError> {
    let bad = |m: &str| ParseError::BadTag(m.to_string());
    let (element_type, rest) = match value.split_once(',') {
        Some((t, r)) => (t, r),
        None => (value, ""),
    };
    if element_type.len() != 1 {
        return Err(bad("Unrecognized element type for array tag value"));
    }
    let letter = element_type.as_bytes()[0];
    let unsigned = letter.is_ascii_uppercase();
    let items: Vec<&str> = if rest.is_empty() {
        Vec::new()
    } else {
        rest.split(',').collect()
    };

    Ok(match letter.to_ascii_lowercase() {
        b'f' => TagValue::FloatArray(
            items
                .iter()
                .map(|s| s.parse::<f32>())
                .collect::<Result<_, _>>()
                .map_err(|_| {
                    bad("Array tag of type f should have single-precision floating point value")
                })?,
        ),
        b'c' => TagValue::ByteArray {
            // An unsigned array was widened before rendering, so parsing it back reads values
            // up to 255 and narrows them into the same signed byte the binary form holds.
            values: items
                .iter()
                .map(|s| s.parse::<i64>().map(|v| v as i8))
                .collect::<Result<_, _>>()
                .map_err(|_| bad("bad byte array element"))?,
            unsigned,
        },
        b's' => TagValue::ShortArray {
            values: items
                .iter()
                .map(|s| s.parse::<i64>().map(|v| v as i16))
                .collect::<Result<_, _>>()
                .map_err(|_| bad("bad short array element"))?,
            unsigned,
        },
        b'i' => TagValue::IntArray {
            values: items
                .iter()
                .map(|s| s.parse::<i64>().map(|v| v as i32))
                .collect::<Result<_, _>>()
                .map_err(|_| bad("bad int array element"))?,
            unsigned,
        },
        _ => return Err(bad("Unrecognized element type for array tag value")),
    })
}

/// `SAMLineParser.parseLine`.
///
/// `resolve` maps a reference name to its dictionary index, returning `None` for a name the
/// header does not carry. Reference names are indices in a `BamRecord`, so the caller supplies
/// the header's view rather than the parser guessing.
pub fn parse_line<F>(line: &str, resolve: F) -> Result<BamRecord, ParseError>
where
    F: FnMut(&str) -> Option<i32>,
{
    parse_line_with(line, resolve, ValidationStringency::default())
}

/// `SAMLineParser.parseLine` at an explicit stringency.
///
/// Only the *consistency* checks are governed by it, matching `reportErrorParsingLine`. A
/// malformed integer or an unparseable tag is a fatal error at every stringency, because
/// `reportFatalErrorParsingLine` throws unconditionally.
pub fn parse_line_with<F>(
    line: &str,
    mut resolve: F,
    stringency: ValidationStringency,
) -> Result<BamRecord, ParseError>
where
    F: FnMut(&str) -> Option<i32>,
{
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < NUM_REQUIRED_FIELDS {
        return Err(ParseError::NotEnoughFields(fields.len()));
    }
    for (i, f) in fields.iter().enumerate() {
        if f.is_empty() {
            return Err(ParseError::EmptyField(i));
        }
    }

    let int_of = |s: &str, field: &'static str| -> Result<i32, ParseError> {
        s.parse().map_err(|_| ParseError::BadInteger {
            field,
            value: s.to_string(),
        })
    };

    let flags: u16 = fields[1].parse().map_err(|_| ParseError::BadInteger {
        field: "FLAG",
        value: fields[1].to_string(),
    })?;
    let read_unmapped = flags & crate::record::READ_UNMAPPED_FLAG != 0;

    // Consistency failures are reported through here so the stringency governs them in one
    // place, exactly as `reportErrorParsingLine` does.
    let inconsistent = |msg: &str| -> Result<(), ParseError> {
        if stringency.rejects() {
            Err(ParseError::Inconsistent(msg.to_string()))
        } else {
            Ok(())
        }
    };

    let rname = fields[2];
    let reference_index = if rname == "*" {
        if !read_unmapped {
            inconsistent("RNAME is not specified but flags indicate mapped")?;
        }
        -1
    } else {
        resolve(rname).unwrap_or(-1)
    };

    let pos = int_of(fields[3], "POS")?;
    let mapq = int_of(fields[4], "MAPQ")?;
    let cigar_text = fields[5];

    // The consistency rules are asymmetric: with a reference name, POS must be non-zero and a
    // mapped read must have a real CIGAR; without one, all three of POS, MAPQ and CIGAR must be
    // absent. Only the second group is checked in both directions.
    if rname != "*" {
        if pos == 0 {
            inconsistent("POS must be non-zero if RNAME is specified")?;
        }
        if !read_unmapped && cigar_text == "*" {
            inconsistent("CIGAR must not be '*' if RNAME is specified")?;
        }
    } else {
        if pos != 0 {
            inconsistent("POS must be zero if RNAME is not specified")?;
        }
        if mapq != 0 {
            inconsistent("MAPQ must be zero if RNAME is not specified")?;
        }
        if cigar_text != "*" {
            inconsistent("CIGAR must be '*' if RNAME is not specified")?;
        }
    }

    let mate_rname = fields[6];
    let mate_reference_index = if mate_rname == "*" {
        -1
    } else if mate_rname == "=" {
        // `=` means "the same reference as this record", which is why an unplaced record
        // cannot use it. At a lenient stringency the record survives with no mate reference,
        // which is what htsjdk's `setMateReferenceName` leaves it as.
        if reference_index < 0 {
            inconsistent("MRNM is '=', but RNAME is not set")?;
        }
        reference_index
    } else {
        resolve(mate_rname).unwrap_or(-1)
    };

    let mpos = int_of(fields[7], "MPOS")?;
    let isize_ = int_of(fields[8], "ISIZE")?;

    let seq = fields[9];
    let read_bases = if seq == "*" {
        Vec::new()
    } else {
        if let Some(bad) = seq.bytes().find(|b| !is_valid_read_base(*b)) {
            return Err(ParseError::InvalidReadBase(bad as char));
        }
        seq.as_bytes().to_vec()
    };

    let qual = fields[10];
    let base_qualities = if qual == "*" {
        Vec::new()
    } else {
        fastq_to_phred(qual)
    };

    let mut tags = Tags::new();
    for field in &fields[NUM_REQUIRED_FIELDS..] {
        let (tag, value) = parse_tag(field)?;
        tags.insert(tag, value);
    }

    Ok(BamRecord {
        read_name: fields[0].to_string(),
        flags,
        reference_index,
        alignment_start: pos,
        mapping_quality: mapq as u8,
        cigar: parse_cigar(cigar_text)?,
        mate_reference_index,
        mate_alignment_start: mpos,
        inferred_insert_size: isize_,
        read_bases,
        base_qualities,
        tags,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::text::write_alignment;

    fn resolver(name: &str) -> Option<i32> {
        match name {
            "chr1" => Some(0),
            "chr2" => Some(1),
            _ => None,
        }
    }

    const PLAIN: &str = "read1\t99\tchr1\t100\t60\t4M\t=\t300\t250\tACGT\t?@AB";

    #[test]
    fn a_plain_line_parses_to_its_fields() {
        let r = parse_line(PLAIN, resolver).unwrap();
        assert_eq!(r.read_name, "read1");
        assert_eq!(r.flags, 99);
        assert_eq!(r.reference_index, 0);
        assert_eq!(r.alignment_start, 100);
        assert_eq!(r.mapping_quality, 60);
        assert_eq!(r.cigar.to_text(), "4M");
        assert_eq!(r.mate_reference_index, 0, "= resolves to RNAME");
        assert_eq!(r.read_bases, b"ACGT");
        assert_eq!(r.base_qualities, vec![30, 31, 32, 33]);
    }

    #[test]
    fn a_line_round_trips_through_the_writer() {
        let r = parse_line(PLAIN, resolver).unwrap();
        assert_eq!(write_alignment(&r, "chr1", "chr1").unwrap(), PLAIN);
    }

    #[test]
    fn too_few_fields_is_refused() {
        assert_eq!(
            parse_line("a\t1\t*\t0\t0\t*", resolver),
            Err(ParseError::NotEnoughFields(6))
        );
    }

    /// An empty field is an error, not an empty value. The position is zero-based in the
    /// message, which is worth reproducing because it is what a user will read.
    #[test]
    fn an_empty_field_is_refused_with_its_zero_based_position() {
        let line = "read1\t99\tchr1\t100\t60\t4M\t=\t300\t250\t\t?@AB";
        assert_eq!(parse_line(line, resolver), Err(ParseError::EmptyField(9)));
    }

    #[test]
    fn a_mapped_read_without_a_reference_name_is_refused() {
        let line = "read1\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\t?@AB";
        assert!(matches!(
            parse_line(line, resolver),
            Err(ParseError::Inconsistent(_))
        ));
    }

    /// Without a reference name, all three of POS, MAPQ and CIGAR must be absent.
    #[test]
    fn an_unplaced_read_must_have_no_position_no_mapq_and_no_cigar() {
        for (line, what) in [
            ("read1\t4\t*\t5\t0\t*\t*\t0\t0\tACGT\t?@AB", "POS"),
            ("read1\t4\t*\t0\t7\t*\t*\t0\t0\tACGT\t?@AB", "MAPQ"),
            ("read1\t4\t*\t0\t0\t4M\t*\t0\t0\tACGT\t?@AB", "CIGAR"),
        ] {
            assert!(
                matches!(parse_line(line, resolver), Err(ParseError::Inconsistent(_))),
                "{what} must be refused"
            );
        }
        // And the consistent form is accepted.
        assert!(parse_line("read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t?@AB", resolver).is_ok());
    }

    #[test]
    fn an_equals_mate_without_a_reference_is_refused() {
        let line = "read1\t1\t*\t0\t0\t*\t=\t0\t0\tACGT\t?@AB";
        assert!(matches!(
            parse_line(line, resolver),
            Err(ParseError::Inconsistent(_))
        ));
    }

    #[test]
    fn absent_sequence_and_qualities_are_empty_not_stars() {
        let line = "read1\t99\tchr1\t100\t60\t4M\t=\t300\t250\t*\t*";
        let r = parse_line(line, resolver).unwrap();
        assert!(r.read_bases.is_empty());
        assert!(r.base_qualities.is_empty());
    }

    /// This is the path on which lower case is reachable. A BAM cannot carry it; SAM text can.
    #[test]
    fn lower_case_bases_survive_the_text_parser() {
        let line = "read1\t99\tchr1\t100\t60\t4M\t=\t300\t250\tacgt\t?@AB";
        assert_eq!(parse_line(line, resolver).unwrap().read_bases, b"acgt");
    }

    /// A character outside the IUPAC table is a parse error, not folded into `N`. That differs
    /// from `bases::base_to_nibble`, which folds `.` but refuses everything else too.
    #[test]
    fn an_invalid_base_is_refused_rather_than_folded() {
        let line = "read1\t99\tchr1\t100\t60\t4M\t=\t300\t250\tACZT\t?@AB";
        assert_eq!(
            parse_line(line, resolver),
            Err(ParseError::InvalidReadBase('Z'))
        );
        // `.` and `=` are in the table.
        for seq in ["AC.T", "AC=T"] {
            let l = format!("read1\t99\tchr1\t100\t60\t4M\t=\t300\t250\t{seq}\t?@AB");
            assert!(parse_line(&l, resolver).is_ok(), "{seq} must be accepted");
        }
    }

    /// The width the text form discarded is re-derived by the binary encoder. This is the
    /// property that makes a BAM → SAM → BAM round trip byte-identical.
    #[test]
    fn integer_widths_are_re_derived_from_the_value() {
        use crate::tag::integer_type;
        for (v, expected) in [(200i64, b'C'), (300, b's'), (40_000, b'S'), (70_000, b'i')] {
            let line = format!("read1\t99\tchr1\t100\t60\t4M\t=\t300\t250\tACGT\t?@AB\tXI:i:{v}");
            let r = parse_line(&line, resolver).unwrap();
            match r.tags.get(Tag::new(b"XI")) {
                Some(TagValue::Int(n)) => {
                    assert_eq!(*n, v);
                    assert_eq!(
                        integer_type(*n).unwrap(),
                        expected,
                        "value {v} must re-derive to '{}'",
                        expected as char
                    );
                }
                other => panic!("expected an int tag, got {other:?}"),
            }
        }
    }

    #[test]
    fn an_integer_outside_the_spec_range_is_refused() {
        let line = format!(
            "read1\t99\tchr1\t100\t60\t4M\t=\t300\t250\tACGT\t?@AB\tXI:i:{}",
            u32::MAX as i64 + 1
        );
        assert!(matches!(
            parse_line(&line, resolver),
            Err(ParseError::BadTag(_))
        ));
    }

    #[test]
    fn arrays_parse_back_with_their_element_type() {
        let line =
            "read1\t99\tchr1\t100\t60\t4M\t=\t300\t250\tACGT\t?@AB\tXB:B:i,1,2,3\tYB:B:C,255,0";
        let r = parse_line(line, resolver).unwrap();
        assert_eq!(
            r.tags.get(Tag::new(b"XB")),
            Some(&TagValue::IntArray {
                values: vec![1, 2, 3],
                unsigned: false
            })
        );
        // 255 narrows to the signed byte the binary form holds, and the unsigned flag records
        // how to render it again.
        assert_eq!(
            r.tags.get(Tag::new(b"YB")),
            Some(&TagValue::ByteArray {
                values: vec![-1, 0],
                unsigned: true
            })
        );
    }

    #[test]
    fn an_empty_array_is_accepted() {
        let line = "read1\t99\tchr1\t100\t60\t4M\t=\t300\t250\tACGT\t?@AB\tXB:B:i";
        let r = parse_line(line, resolver).unwrap();
        assert_eq!(
            r.tags.get(Tag::new(b"XB")),
            Some(&TagValue::IntArray {
                values: vec![],
                unsigned: false
            })
        );
    }

    #[test]
    fn a_cigar_without_a_length_is_refused() {
        for c in ["M", "4", "4Z", "4M2"] {
            assert!(
                parse_cigar(c).is_err(),
                "{c} must be refused, got {:?}",
                parse_cigar(c)
            );
        }
        assert_eq!(parse_cigar("*").unwrap(), Cigar::default());
    }

    /// The asymmetry that made the stringency model necessary: htsjdk's writer emits a record
    /// its own default-stringency reader rejects. Both behaviours are reproduced.
    #[test]
    fn a_record_htsjdk_writes_but_strictly_rejects_is_readable_leniently() {
        // RNAME '*' with MAPQ 60: htsjdk's SAMFileWriter writes this without complaint.
        let line = "unplaced\t4\t*\t0\t60\t*\t*\t0\t0\tACGT\t?@AB";

        assert!(
            matches!(parse_line(line, resolver), Err(ParseError::Inconsistent(_))),
            "STRICT is the default and must reject it"
        );

        for lenient in [ValidationStringency::Lenient, ValidationStringency::Silent] {
            let r = parse_line_with(line, resolver, lenient).unwrap();
            assert_eq!(r.mapping_quality, 60, "the value survives, unchanged");
            assert_eq!(r.reference_index, -1);
        }
    }

    /// A malformed integer is fatal at every stringency, because htsjdk's
    /// `reportFatalErrorParsingLine` throws unconditionally rather than consulting it.
    #[test]
    fn a_malformed_field_is_fatal_even_when_lenient() {
        let line = "read1\t99\tchr1\tNOTANUMBER\t60\t4M\t=\t300\t250\tACGT\t?@AB";
        for s in [
            ValidationStringency::Strict,
            ValidationStringency::Lenient,
            ValidationStringency::Silent,
        ] {
            assert!(
                matches!(
                    parse_line_with(line, resolver, s),
                    Err(ParseError::BadInteger { .. })
                ),
                "{s:?} must still refuse a malformed integer"
            );
        }
    }

    #[test]
    fn phred_conversion_is_the_inverse_of_the_writer() {
        use crate::text::phred_to_fastq;
        let quals: Vec<u8> = (0..60).collect();
        assert_eq!(fastq_to_phred(&phred_to_fastq(&quals)), quals);
    }
}
