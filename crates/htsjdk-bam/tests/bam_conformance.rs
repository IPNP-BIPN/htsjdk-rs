//! Conformance against htsjdk's own `BAMRecordCodec`.
//!
//! The goldens in `tests/data/bam_codec.txt` were produced by `tools/bam-conformance/`
//! running inside the pinned oracle container. Nothing in this file decides what the right
//! answer is; htsjdk did, and this asserts we reproduce it byte for byte.
//!
//! The case list is regenerated here in the same order as the Java harness emits it, and the
//! case *names* are compared as well as the bytes. That guards the one failure mode this
//! design has: if the two lists drifted out of step, index pairing would compare the wrong
//! records and could pass by coincidence.

// `3.14159f` below is the literal the Java harness wrote, and the golden holds the f32 bits of
// exactly that literal. Substituting `std::f32::consts::PI`, as clippy suggests, would be a
// different number and would fail against the golden.
#![allow(clippy::approx_constant)]

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::record::{BamRecord, READ_UNMAPPED_FLAG};
use htsjdk_bam::tag::{Tag, TagValue};

/// The corpus is gzipped: the long-CIGAR cases alone are 2.8 MB of hex, and it compresses to
/// 7 KB. Same convention as `jmath`'s conformance corpus.
fn golden_text() -> String {
    use std::io::Read;
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bam_codec.txt.gz");
    let f = std::fs::File::open(&p).expect("golden corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("golden corpus is gzip");
    s
}

/// One case: the name the Java harness gave it and the record we claim reproduces it.
struct Case {
    name: String,
    record: BamRecord,
}

fn base() -> BamRecord {
    BamRecord {
        read_name: "read1".into(),
        flags: 0,
        reference_index: 0,
        alignment_start: 100,
        mapping_quality: 60,
        cigar: cigar("4M"),
        mate_reference_index: -1,
        mate_alignment_start: 0,
        inferred_insert_size: 0,
        read_bases: b"ACGT".to_vec(),
        base_qualities: vec![30, 30, 30, 30],
        tags: Default::default(),
    }
}

/// `TextCigarCodec.decode`, enough of it for the fixtures.
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
                _ => panic!("bad cigar operator {}", b as char),
            };
            elements.push(CigarElement { length: len, op });
            len = 0;
        }
    }
    Cigar::new(elements)
}

fn tag(name: &str) -> Tag {
    Tag::new(name.as_bytes().try_into().unwrap())
}

/// `java.lang.String.hashCode`, needed only because the harness puts it in a case name.
fn java_hash(s: &str) -> i32 {
    s.encode_utf16()
        .fold(0i32, |h, u| h.wrapping_mul(31).wrapping_add(u as i32))
}

fn filled(n: usize, base_char: u8, qual: u8) -> (Vec<u8>, Vec<u8>) {
    (vec![base_char; n], vec![qual; n])
}

/// Rebuilds the harness's case list, in order.
fn cases() -> Vec<Case> {
    let mut out: Vec<Case> = Vec::new();
    let mut push = |name: String, record: BamRecord| out.push(Case { name, record });

    push("plain".into(), base());

    for start in [
        1, 2, 16384, 16385, 16386, 131072, 131073, 1048576, 8388608, 67108864, 100000000, 250000000,
    ] {
        let mut r = base();
        r.alignment_start = start;
        push(format!("start_{start}"), r);
    }

    // The reference span comes from a deletion, so the read stays four bases long whatever the
    // span. Same bin arithmetic as a 64 Mb match, four orders of magnitude less golden.
    for span in [1u32, 100, 16384, 131072, 1048576, 8388608, 67108864] {
        let mut r = base();
        r.alignment_start = 16300;
        r.cigar = cigar(&format!("4M{span}D"));
        push(format!("span_{span}"), r);
    }

    let mut unmapped = base();
    unmapped.flags |= READ_UNMAPPED_FLAG;
    unmapped.reference_index = -1;
    unmapped.alignment_start = 0;
    unmapped.cigar = Cigar::default();
    push("unmapped".into(), unmapped);

    let mut placed_unmapped = base();
    placed_unmapped.flags |= READ_UNMAPPED_FLAG;
    push("placed_unmapped".into(), placed_unmapped);

    for seq in [
        "A",
        "AC",
        "ACG",
        "ACGT",
        "ACGTN",
        "=ACMGRSVTWYHKDBN",
        "acgt",
        "ACGTacgtNn.=",
    ] {
        let mut r = base();
        r.read_bases = seq.as_bytes().to_vec();
        r.cigar = cigar(&format!("{}M", seq.len()));
        r.base_qualities = (0..seq.len()).map(|i| (i % 60) as u8).collect();
        push(format!("seq_{seq}"), r);
    }

    let mut no_quals = base();
    no_quals.base_qualities = Vec::new();
    push("no_quals".into(), no_quals);

    for (i, name) in ["a", "read/1", &"x".repeat(100), &"x".repeat(254)]
        .iter()
        .enumerate()
    {
        let mut r = base();
        r.read_name = (*name).to_string();
        push(format!("name_{i}_len{}", name.len()), r);
    }

    for c in [
        "4M",
        "2M2I",
        "2M2D2M",
        "1S2M1S",
        "2H4M2H",
        "4=",
        "4X",
        "1M1I1D1N1S1H1P1=1X",
    ] {
        let mut r = base();
        r.cigar = cigar(c);
        let (bases, quals) = filled(r.cigar.read_length() as usize, b'A', 30);
        r.read_bases = bases;
        r.base_qualities = quals;
        push(format!("cigar_{c}"), r);
    }

    // The integer promotion ladder: the point of the whole exercise.
    for v in [
        i32::MIN as i64,
        -32769,
        -32768,
        -32767,
        -129,
        -128,
        -127,
        -1,
        0,
        1,
        126,
        127,
        128,
        129,
        200,
        254,
        255,
        256,
        257,
        300,
        32766,
        32767,
        32768,
        32769,
        65534,
        65535,
        65536,
        65537,
        2147483646,
        2147483647,
        2147483648,
        4294967294,
        4294967295,
    ] {
        let mut r = base();
        r.tags.insert(tag("XI"), TagValue::Int(v));
        push(format!("int_{v}"), r);
    }

    // The declared Java box type must not influence the bytes; only the value may.
    for v in [100i32, 30000] {
        for (kind, value) in [
            ("byte", (v as i8) as i64),
            ("short", (v as i16) as i64),
            ("int", v as i64),
            ("long", v as i64),
        ] {
            let mut r = base();
            r.tags.insert(tag("XI"), TagValue::Int(value));
            push(format!("box_{kind}_{v}"), r);
        }
    }

    let mut char_tag = base();
    char_tag.tags.insert(tag("XA"), TagValue::Char(b'Z'));
    push("tag_char".into(), char_tag);

    let mut float_tag = base();
    float_tag.tags.insert(tag("XF"), TagValue::Float(1.5));
    push("tag_float".into(), float_tag);

    for f in [
        0.0f32,
        -0.0,
        1.0,
        -1.0,
        3.14159,
        f32::from_bits(1),
        f32::MAX,
        f32::NAN,
        f32::INFINITY,
    ] {
        let mut r = base();
        r.tags.insert(tag("XF"), TagValue::Float(f));
        push(format!("float_{}", f.to_bits() as i32), r);
    }

    for s in ["", "a", "hello world", "100", "é"] {
        let mut r = base();
        r.tags.insert(tag("XS"), TagValue::Str(s.to_string()));
        push(
            format!("tag_str_{}_{}", s.encode_utf16().count(), java_hash(s)),
            r,
        );
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
        ("arr_float", TagValue::FloatArray(vec![1.0, -2.5])),
        (
            "arr_empty",
            TagValue::IntArray {
                values: vec![],
                unsigned: false,
            },
        ),
    ];
    for (name, value) in arrays {
        let mut r = base();
        r.tags.insert(tag("XB"), value.clone());
        push((*name).into(), r);
    }

    let mut ordered = base();
    for t in ["ZA", "AZ", "NM", "MD", "AS", "XS", "SA", "aa", "Aa"] {
        ordered.tags.insert(tag(t), TagValue::Int(1));
    }
    push("tag_order".into(), ordered);

    let mut realistic = base();
    realistic.tags.insert(tag("NM"), TagValue::Int(2));
    realistic
        .tags
        .insert(tag("MD"), TagValue::Str("2A1".into()));
    realistic.tags.insert(tag("AS"), TagValue::Int(100));
    realistic.tags.insert(tag("XS"), TagValue::Int(-30));
    realistic
        .tags
        .insert(tag("RG"), TagValue::Str("rg1".into()));
    realistic
        .tags
        .insert(tag("PG"), TagValue::Str("MarkDuplicates".into()));
    push("realistic".into(), realistic);

    let long_cigar = |n: usize| {
        Cigar::new(
            (0..n)
                .map(|i| CigarElement {
                    length: 1,
                    op: if i % 2 == 0 { Op::M } else { Op::I },
                })
                .collect(),
        )
    };
    for n in [65535usize, 65536, 65537] {
        let c = long_cigar(n);
        let mut r = base();
        let (bases, quals) = filled(c.read_length() as usize, b'A', 30);
        r.cigar = c;
        r.read_bases = bases;
        r.base_qualities = quals;
        push(format!("longcigar_{n}"), r);
    }

    let c = long_cigar(65536);
    let mut with_tags = base();
    let (bases, quals) = filled(c.read_length() as usize, b'A', 30);
    with_tags.cigar = c;
    with_tags.read_bases = bases;
    with_tags.base_qualities = quals;
    with_tags.tags.insert(tag("AG"), TagValue::Int(1));
    with_tags.tags.insert(tag("CH"), TagValue::Int(1));
    with_tags.tags.insert(tag("NM"), TagValue::Int(1));
    push("longcigar_with_tags".into(), with_tags);

    for flags in [0u16, 1, 4, 16, 99, 147, 1024, 2048, 4095] {
        let mut r = base();
        r.flags = flags;
        push(format!("flags_{flags}"), r);
    }
    for isize_ in [i32::MIN, -1000, -1, 0, 1, 1000, i32::MAX] {
        let mut r = base();
        r.inferred_insert_size = isize_;
        push(format!("isize_{isize_}"), r);
    }
    let mut mated = base();
    mated.mate_reference_index = 1;
    mated.mate_alignment_start = 5000;
    mated.inferred_insert_size = 300;
    push("mated".into(), mated);

    for mq in [0u8, 1, 60, 254, 255] {
        let mut r = base();
        r.mapping_quality = mq;
        push(format!("mapq_{mq}"), r);
    }

    out
}

fn goldens() -> Vec<(String, String)> {
    golden_text()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let (name, hex) = l
                .split_once('\t')
                .expect("golden line must be name TAB hex");
            (name.to_string(), hex.to_string())
        })
        .collect()
}

fn to_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// The gate: every record must encode to exactly the bytes htsjdk produced.
#[test]
fn every_record_encodes_to_htsjdks_bytes() {
    let cases = cases();
    let goldens = goldens();
    assert_eq!(
        cases.len(),
        goldens.len(),
        "case list and golden list are different lengths; the two have drifted"
    );

    let mut failures = Vec::new();
    for (case, (golden_name, golden_hex)) in cases.iter().zip(&goldens) {
        assert_eq!(
            &case.name, golden_name,
            "case lists out of step: pairing would compare the wrong records"
        );
        let encoded = match case.record.encode() {
            Ok(b) => to_hex(&b),
            Err(e) => {
                failures.push(format!("{}: encode failed: {e:?}", case.name));
                continue;
            }
        };
        if &encoded != golden_hex {
            let at = encoded
                .as_bytes()
                .iter()
                .zip(golden_hex.as_bytes())
                .position(|(a, b)| a != b)
                .unwrap_or(encoded.len().min(golden_hex.len()));
            failures.push(format!(
                "{}: first differs at byte {} (hex char {at})\n  ours:   {}\n  htsjdk: {}",
                case.name,
                at / 2,
                &encoded[at.saturating_sub(8)..(at + 24).min(encoded.len())],
                &golden_hex[at.saturating_sub(8)..(at + 24).min(golden_hex.len())],
            ));
        }
    }

    assert!(
        failures.is_empty(),
        "{} of {} cases diverge from htsjdk:\n{}",
        failures.len(),
        cases.len(),
        failures.join("\n")
    );
}

/// Decoding htsjdk's bytes and re-encoding must return the same bytes.
///
/// Separate from the encode test on purpose: a decoder that drops a field the encoder then
/// re-derives would pass the encode test and lose data on any real file.
#[test]
fn htsjdks_bytes_survive_a_decode_encode_round_trip() {
    let mut failures = Vec::new();
    for (name, hex) in goldens() {
        let bytes: Vec<u8> = (0..hex.len() / 2)
            .map(|i| u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16).unwrap())
            .collect();
        let (record, used) = match BamRecord::decode(&bytes) {
            Ok(Some(v)) => v,
            other => {
                failures.push(format!("{name}: decode gave {other:?}"));
                continue;
            }
        };
        if used != bytes.len() {
            failures.push(format!("{name}: consumed {used} of {} bytes", bytes.len()));
            continue;
        }
        match record.encode() {
            Ok(again) if again == bytes => {}
            Ok(again) => failures.push(format!(
                "{name}: re-encode differs ({} vs {} bytes)",
                again.len(),
                bytes.len()
            )),
            Err(e) => failures.push(format!("{name}: re-encode failed: {e:?}")),
        }
    }
    assert!(
        failures.is_empty(),
        "{} round-trip failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// The goldens must actually exercise the promotion ladder. A conformance suite that happens
/// to contain only `i` tags would pass while proving nothing about the interesting part.
#[test]
fn the_goldens_cover_every_integer_tag_type() {
    let mut seen = std::collections::BTreeSet::new();
    for (name, hex) in goldens() {
        if !name.starts_with("int_") && !name.starts_with("box_") {
            continue;
        }
        // The XI tag is the last thing in these records: "5849" then the type byte.
        let at = hex.rfind("5849").expect("XI tag present");
        let ty = u8::from_str_radix(&hex[at + 4..at + 6], 16).unwrap();
        seen.insert(ty as char);
    }
    let expected: std::collections::BTreeSet<char> = "cCsSiI".chars().collect();
    assert_eq!(
        seen, expected,
        "the goldens must contain every integer tag type htsjdk can emit"
    );
}
