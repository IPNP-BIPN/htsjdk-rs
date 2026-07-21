//! Conformance for the SAM text record writer against htsjdk's `SAMTextWriter`.
//!
//! Goldens from `tools/sam-text-conformance/SamTextDump.java` in the pinned oracle.
//!
//! The property this suite exists for is the counterpoint to decision 0008: `TextTagCodec`
//! collapses every integer tag width to `i`, discarding exactly what `BinaryTagCodec`'s
//! promotion ladder computes. The corpus contains the whole ladder so that is checked rather
//! than asserted.
//!
//! **One case is declared as diverging**, and it is the same cause decision 0011 recorded for
//! doubles, now confirmed for floats: Java 17's `Float.toString` is not the shortest
//! round-trip decimal. `Float.MIN_VALUE` renders as `1.4E-45` there and `1.0E-45` here, and
//! both parse back to the same subnormal. Closing it needs the `FloatingDecimal` port that
//! decision 0011 already lists as outstanding, so it is pinned here rather than patched.

// `3.14159f32` below is the literal the Java harness wrote, and the golden holds the bits of
// exactly that literal. Substituting `std::f32::consts::PI`, as clippy suggests, is a different
// number and fails against the golden.
#![allow(clippy::approx_constant)]

use std::io::Read;

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::text::write_alignment;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sam_text.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn goldens() -> Vec<(String, String)> {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let (name, line) = l.split_once('\t').expect("name TAB line");
            (name.to_string(), line.to_string())
        })
        .collect()
}

fn t(name: &str) -> Tag {
    Tag::new(name.as_bytes().try_into().unwrap())
}

fn cigar(text: &str) -> Cigar {
    if text == "*" {
        return Cigar::default();
    }
    let mut elements = Vec::new();
    let mut len = 0u32;
    for b in text.bytes() {
        if b.is_ascii_digit() {
            len = len * 10 + (b - b'0') as u32;
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
                _ => panic!("bad operator"),
            };
            elements.push(CigarElement { length: len, op });
            len = 0;
        }
    }
    Cigar::new(elements)
}

fn java_hash(s: &str) -> i32 {
    s.encode_utf16()
        .fold(0i32, |h, u| h.wrapping_mul(31).wrapping_add(u as i32))
}

/// Reference names by index, as the harness's header defines them.
fn name_of(index: i32) -> &'static str {
    match index {
        0 => "chr1",
        1 => "chr2",
        _ => "*",
    }
}

struct Case {
    name: String,
    record: BamRecord,
}

fn base() -> BamRecord {
    BamRecord {
        read_name: "read1".into(),
        flags: 99,
        reference_index: 0,
        alignment_start: 100,
        mapping_quality: 60,
        cigar: cigar("4M"),
        mate_reference_index: 0,
        mate_alignment_start: 300,
        inferred_insert_size: 250,
        read_bases: b"ACGT".to_vec(),
        base_qualities: vec![30, 31, 32, 33],
        tags: Default::default(),
    }
}

#[allow(clippy::vec_init_then_push)]
fn cases() -> Vec<Case> {
    let mut out: Vec<Case> = Vec::new();
    let mut push = |name: String, record: BamRecord| out.push(Case { name, record });

    push("plain".into(), base());

    let mut other = base();
    other.mate_reference_index = 1;
    push("mate_other_reference".into(), other);

    let mut unplaced = base();
    unplaced.reference_index = -1;
    unplaced.alignment_start = 0;
    unplaced.mate_reference_index = -1;
    unplaced.mate_alignment_start = 0;
    unplaced.cigar = Cigar::default();
    unplaced.flags |= 0x4;
    push("unplaced".into(), unplaced);

    let mut no_seq = base();
    no_seq.read_bases = Vec::new();
    no_seq.base_qualities = Vec::new();
    push("no_sequence".into(), no_seq);

    let mut no_qual = base();
    no_qual.base_qualities = Vec::new();
    push("no_quals".into(), no_qual);

    for v in [
        i32::MIN as i64,
        -32769,
        -32768,
        -129,
        -128,
        -1,
        0,
        1,
        127,
        128,
        200,
        255,
        256,
        300,
        32767,
        32768,
        65535,
        65536,
        2147483647,
        2147483648,
        4294967295,
    ] {
        let mut r = base();
        r.tags.insert(t("XI"), TagValue::Int(v));
        push(format!("int_{v}"), r);
    }

    let mut char_tag = base();
    char_tag.tags.insert(t("XA"), TagValue::Char(b'Z'));
    push("tag_char".into(), char_tag);

    for f in [
        0.0f32,
        -0.0,
        1.0,
        -1.0,
        0.5,
        3.14159,
        1e10,
        1e-10,
        f32::from_bits(1),
        f32::MAX,
        f32::NAN,
        f32::INFINITY,
        f32::NEG_INFINITY,
    ] {
        let mut r = base();
        r.tags.insert(t("XF"), TagValue::Float(f));
        push(format!("float_{}", f.to_bits() as i32), r);
    }

    for s in ["", "hello", "with space", "100"] {
        let mut r = base();
        r.tags.insert(t("XS"), TagValue::Str(s.to_string()));
        push(format!("str_{}_{}", s.len(), java_hash(s)), r);
    }

    let arrays: &[(&str, TagValue)] = &[
        (
            "arr_byte_signed",
            TagValue::ByteArray {
                values: vec![-1, 0, 1, 127],
                unsigned: false,
            },
        ),
        (
            "arr_byte_unsigned",
            TagValue::ByteArray {
                values: vec![-1, 0, 1, 127],
                unsigned: true,
            },
        ),
        (
            "arr_short_signed",
            TagValue::ShortArray {
                values: vec![-1, 0, 300, 32767],
                unsigned: false,
            },
        ),
        (
            "arr_short_unsigned",
            TagValue::ShortArray {
                values: vec![-1, 0, 300, 32767],
                unsigned: true,
            },
        ),
        (
            "arr_int_signed",
            TagValue::IntArray {
                values: vec![-1, 0, 100000, i32::MAX],
                unsigned: false,
            },
        ),
        (
            "arr_int_unsigned",
            TagValue::IntArray {
                values: vec![-1, 0, 100000, i32::MAX],
                unsigned: true,
            },
        ),
        ("arr_float", TagValue::FloatArray(vec![1.0, -2.5, 1e10])),
    ];
    for (name, value) in arrays {
        let mut r = base();
        r.tags.insert(t("XB"), value.clone());
        push((*name).into(), r);
    }

    let mut ordered = base();
    for name in ["ZA", "AZ", "NM", "MD", "AS"] {
        ordered.tags.insert(t(name), TagValue::Int(1));
    }
    push("tag_order".into(), ordered);

    for c in [
        "4M",
        "2M2I",
        "1S2M1S",
        "2H4M2H",
        "4=",
        "4X",
        "1M1I1D1N1S1H1P1=1X",
    ] {
        let mut r = base();
        r.cigar = cigar(c);
        let len = r.cigar.read_length() as usize;
        r.read_bases = vec![b'A'; len];
        r.base_qualities = vec![30; len];
        push(format!("cigar_{c}"), r);
    }

    for flags in [0u16, 4, 16, 99, 147, 1024, 2048, 4095] {
        let mut r = base();
        r.flags = flags;
        push(format!("flags_{flags}"), r);
    }

    out
}

/// Cases where the port knowingly differs, with the cause.
///
/// Pinned rather than excluded: a suite trimmed to what passes reports 100% and means nothing.
const KNOWN_DIVERGENCES: &[(&str, &str, &str)] = &[("float_1", "XF:f:1.0E-45", "XF:f:1.4E-45")];

/// The gate: every record renders to exactly htsjdk's line.
#[test]
fn every_record_renders_to_htsjdks_sam_line() {
    let cases = cases();
    let goldens = goldens();
    assert_eq!(cases.len(), goldens.len(), "case lists have drifted");

    let mut failures = Vec::new();
    let mut declared = 0usize;
    for (case, (golden_name, expected)) in cases.iter().zip(&goldens) {
        assert_eq!(&case.name, golden_name, "case lists out of step");
        let ours = write_alignment(
            &case.record,
            name_of(case.record.reference_index),
            name_of(case.record.mate_reference_index),
        );
        match ours {
            None => failures.push(format!("{}: refused to render", case.name)),
            Some(line) if &line != expected => {
                match KNOWN_DIVERGENCES.iter().find(|(n, _, _)| *n == case.name) {
                    Some((_, ours_expected, theirs_expected)) => {
                        // A declared divergence must still be exactly the recorded one; drifting
                        // to a different wrong answer would otherwise pass unnoticed.
                        assert!(
                            line.ends_with(ours_expected) && expected.ends_with(theirs_expected),
                            "{} drifted: {line}",
                            case.name
                        );
                        declared += 1;
                    }
                    None => failures.push(format!(
                        "{}:\n  ours:   {line}\n  htsjdk: {expected}",
                        case.name
                    )),
                }
            }
            Some(_) => {}
        }
    }
    assert_eq!(
        declared,
        KNOWN_DIVERGENCES.len(),
        "a declared divergence now matches; remove it and update decision 0011"
    );
    assert!(
        failures.is_empty(),
        "{} of {} lines diverge:\n{}",
        failures.len(),
        cases.len(),
        failures.join("\n")
    );
}

/// The property the suite exists for, asserted against the goldens rather than the port.
#[test]
fn htsjdk_really_does_collapse_every_integer_width_to_i() {
    let mut seen = 0;
    for (name, line) in goldens() {
        if !name.starts_with("int_") {
            continue;
        }
        seen += 1;
        let tag = line.rsplit('\t').next().unwrap();
        assert!(
            tag.starts_with("XI:i:"),
            "{name} rendered as {tag}, so the text form does distinguish widths after all"
        );
    }
    assert!(
        seen >= 20,
        "the corpus must span the whole ladder, saw {seen}"
    );
}
