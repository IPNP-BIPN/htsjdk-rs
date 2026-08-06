//! Conformance for the read side of the tag codec, against `htsjdk.samtools.BinaryTagCodec`.
//!
//! Goldens from `tools/bam-conformance/BinaryTagCodecReadDump.java` in the pinned oracle.
//!
//! The rows that justify the suite:
//!
//! ```text
//! value  I-small  49414905000000        IA  Integer  false  5              49416305      changed
//! value  H        48684834383635364300  Hh  byte[]   false  [72,101,108]   486842630300000048656c  changed
//! value  B-C      4243424303000000ff007f  BC  byte[]  true  [-1,0,127]     4243424303000000ff007f  same
//! textrt XX:H:48656C  byte[]  XX:B:c,72,101,108  changed
//! ```
//!
//! Reading is not the inverse of writing. An `I` holding a small value comes back as an ordinary
//! integer and is rewritten four bytes shorter; an `H` comes back as the same `byte[]` a `B` array
//! does and is rewritten a byte longer, in the text codec as well as the binary one. The unsigned
//! flag of an array changes the type letter and nothing else, so a `C` array holding `0xFF` still
//! holds `-1`.

use std::io::Read;

use htsjdk_bam::tag::{
    java_char, Tag, TagValue, Tags, FIXED_BINARY_ARRAY_TAG_SIZE, FIXED_TAG_SIZE,
};
use htsjdk_bam::text::{encode_tag, java_float_to_string};
use htsjdk_bam::text_parse::parse_tag;

fn corpus() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bam_tag_read.txt.gz");
    let file = std::fs::File::open(&path).expect("corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("corpus is gzip");
    text
}

fn rows<'a>(corpus: &'a str, kind: &str) -> Vec<Vec<&'a str>> {
    let prefix = format!("{kind}\t");
    corpus
        .lines()
        .filter_map(|line| line.strip_prefix(&prefix))
        .map(|rest| rest.split('\t').collect())
        .collect()
}

fn unhex(text: &str) -> Vec<u8> {
    if text == "-" {
        return Vec::new();
    }
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).expect("hex"))
        .collect()
}

fn hex(bytes: &[u8]) -> String {
    if bytes.is_empty() {
        return "-".to_string();
    }
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// The dump's `escape`: anything outside printable ASCII as `\uXXXX`, over UTF-16 units, and an
/// empty string as `-`.
fn escape(text: &str) -> String {
    let escaped: String = text
        .encode_utf16()
        .map(|unit| {
            if (0x20..=0x7E).contains(&unit) {
                char::from_u32(u32::from(unit)).expect("ascii").to_string()
            } else {
                format!("\\u{unit:04X}")
            }
        })
        .collect();
    if escaped.is_empty() {
        "-".to_string()
    } else {
        escaped
    }
}

fn unescape(text: &str) -> String {
    if text == "-" {
        return String::new();
    }
    let mut units = Vec::new();
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '\\' && chars.peek() == Some(&'u') {
            chars.next();
            let digits: String = (0..4).filter_map(|_| chars.next()).collect();
            units.push(u16::from_str_radix(&digits, 16).expect("escape"));
        } else {
            let mut buffer = [0u16; 2];
            units.extend_from_slice(c.encode_utf16(&mut buffer));
        }
    }
    String::from_utf16(&units).expect("utf-16")
}

/// The Java class the value comes back as.
///
/// Every narrow integer widens to `Integer`; only an `I` above `Integer.MAX_VALUE` is a `Long`,
/// and the value alone says which, because that is all htsjdk keeps.
fn class_of(value: &TagValue) -> &'static str {
    match value {
        TagValue::Char(_) => "Character",
        TagValue::Int(v) if *v > i32::MAX as i64 => "Long",
        TagValue::Int(_) => "Integer",
        TagValue::Float(_) => "Float",
        TagValue::Str(_) => "String",
        TagValue::ByteArray { .. } => "byte[]",
        TagValue::ShortArray { .. } => "short[]",
        TagValue::IntArray { .. } => "int[]",
        TagValue::FloatArray(_) => "float[]",
    }
}

/// The dump's `describe`.
fn describe(value: &TagValue) -> String {
    fn list<T: std::fmt::Display>(values: &[T]) -> String {
        let body: Vec<String> = values.iter().map(|v| v.to_string()).collect();
        format!("[{}]", body.join(","))
    }
    match value {
        // `(char) aSignedByte`, which is not the byte for anything above 0x7F.
        TagValue::Char(c) => format!("U+{:04X}", java_char(*c)),
        TagValue::Int(v) => v.to_string(),
        TagValue::Float(f) => java_float_to_string(*f),
        TagValue::Str(s) => escape(s),
        TagValue::ByteArray { values, .. } => list(values),
        TagValue::ShortArray { values, .. } => list(values),
        TagValue::IntArray { values, .. } => list(values),
        TagValue::FloatArray(values) => {
            let body: Vec<String> = values.iter().map(|f| java_float_to_string(*f)).collect();
            format!("[{}]", body.join(","))
        }
    }
}

fn unsigned_of(value: &TagValue) -> bool {
    match value {
        TagValue::ByteArray { unsigned, .. }
        | TagValue::ShortArray { unsigned, .. }
        | TagValue::IntArray { unsigned, .. } => *unsigned,
        // htsjdk carries the flag on the node's class, and there is no unsigned float node.
        _ => false,
    }
}

/// The constants, measured rather than assumed.
#[test]
fn the_constants_are_the_reference_constants() {
    let corpus = corpus();
    let sizes = rows(&corpus, "sizes");
    assert_eq!(sizes.len(), 1);
    assert_eq!(sizes[0][0], FIXED_TAG_SIZE.to_string());
    assert_eq!(sizes[0][1], FIXED_BINARY_ARRAY_TAG_SIZE.to_string());
}

/// Every block the reference read comes back the same way here, and writing it back gives the
/// same bytes the reference's own writer gave.
#[test]
fn every_value_reads_and_rewrites_as_the_reference_did() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "value") {
        let (label, input, name, class, unsigned, value, output, stable) = (
            row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7],
        );
        let bytes = unhex(input);
        let tags = Tags::read(&bytes).unwrap_or_else(|e| panic!("{label}: {e}"));
        assert_eq!(tags.len(), 1, "{label}: one tag");
        let (tag, read) = tags.iter().next().expect("one tag");

        assert_eq!(tag.to_string(), name, "{label}: tag name");
        assert_eq!(class_of(read), class, "{label}: class");
        assert_eq!(unsigned_of(read).to_string(), unsigned, "{label}: unsigned");
        assert_eq!(describe(read), value, "{label}: value");

        let mut written = Vec::new();
        tags.write(&mut written)
            .unwrap_or_else(|e| panic!("{label}: {e:?}"));
        assert_eq!(hex(&written), output, "{label}: rewritten");
        assert_eq!(
            if written == bytes { "same" } else { "changed" },
            stable,
            "{label}: stability"
        );
        compared += 1;
    }
    assert_eq!(compared, 33, "values compared");
}

/// Exactly two of the thirty-three break the round trip, and both are forms htsjdk never writes.
#[test]
fn only_the_two_types_htsjdk_never_writes_change_on_the_way_back() {
    let corpus = corpus();
    let changed: Vec<&str> = rows(&corpus, "value")
        .into_iter()
        .filter(|row| row[7] == "changed")
        .map(|row| row[0])
        .collect();
    assert_eq!(
        changed,
        ["I-small", "I-int-max", "H", "H-empty", "H-lowercase"]
    );

    // The I that shrinks: four bytes of value become one, and the type letter with it.
    let small = rows(&corpus, "value")
        .into_iter()
        .find(|row| row[0] == "I-small")
        .expect("I-small");
    assert_eq!(unhex(small[1]).len(), unhex(small[6]).len() + 3);
    assert_eq!(unhex(small[1])[2], b'I');
    assert_eq!(unhex(small[6])[2], b'c');

    // The H that grows: a null-terminated hex string becomes a length-prefixed array.
    let hex_tag = rows(&corpus, "value")
        .into_iter()
        .find(|row| row[0] == "H")
        .expect("H");
    assert_eq!(unhex(hex_tag[1])[2], b'H');
    assert_eq!(unhex(hex_tag[6])[2], b'B');
    assert_eq!(unhex(hex_tag[6]).len(), unhex(hex_tag[1]).len() + 1);
}

/// The unsigned flag is the case of the type letter and nothing else: the elements stay signed.
#[test]
fn an_unsigned_array_differs_from_its_signed_twin_only_in_the_letter() {
    let corpus = corpus();
    let find = |label: &str| {
        rows(&corpus, "value")
            .into_iter()
            .find(|row| row[0] == label)
            .map(|row| (row[1].to_string(), row[4].to_string(), row[5].to_string()))
            .expect("row")
    };
    let (signed_in, signed_flag, signed_value) = find("B-c");
    let (unsigned_in, unsigned_flag, unsigned_value) = find("B-C");
    assert_eq!(signed_flag, "false");
    assert_eq!(unsigned_flag, "true");
    assert_eq!(signed_value, unsigned_value, "the same elements either way");
    // The bytes differ in two places: the tag's own second character, which is how the dump keeps
    // the two cases apart, and the element type letter, which is the whole of the flag.
    let (a, b) = (unhex(&signed_in), unhex(&unsigned_in));
    let differing: Vec<usize> = (0..a.len()).filter(|i| a[*i] != b[*i]).collect();
    assert_eq!(differing, [1, 3]);
    assert_eq!((a[1], b[1]), (b'c', b'C'), "the tag names");
    assert_eq!((a[3], b[3]), (b'c', b'C'), "the element type letter");
    assert_eq!(&a[4..], &b[4..], "and the elements themselves are the same");
}

/// A whole block comes back sorted by the packed short, so by the second character first.
#[test]
fn a_block_comes_back_in_the_order_the_reference_gives_it() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "order") {
        let (label, input, expect) = (row[0], row[1], row[2]);
        let tags = Tags::read(&unhex(input)).expect("reads");
        let mine: Vec<String> = tags
            .iter()
            .map(|(tag, value)| format!("{tag}={}", describe(value)))
            .collect();
        assert_eq!(mine.join(","), expect, "{label}");
        compared += 1;
    }
    assert_eq!(compared, 2, "orders compared");
    // Both files hold the same four tags in different orders and come back identically.
    let orders: Vec<String> = rows(&corpus, "order")
        .into_iter()
        .map(|row| {
            row[2]
                .split(',')
                .map(|entry| entry.split('=').next().expect("name").to_string())
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect();
    assert_eq!(orders[0], orders[1]);
    assert_eq!(orders[0], "ZA,MD,NM,AZ");
}

/// A repeated tag replaces rather than duplicates, and the last one in the block wins.
#[test]
fn a_repeated_tag_leaves_one_entry_and_the_last_value() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "dup") {
        let (label, input, entries, survivor) = (row[0], row[1], row[2], row[3]);
        let tags = Tags::read(&unhex(input)).expect("reads");
        assert_eq!(tags.len().to_string(), entries, "{label}: entries");
        let value = tags.get(Tag::new(b"NM")).expect("NM");
        assert_eq!(describe(value), survivor, "{label}: survivor");
        compared += 1;
    }
    assert_eq!(compared, 3, "duplicates compared");
}

/// Every way a block can be malformed, with the exception the reference raises.
#[test]
fn the_read_errors_are_the_reference_errors() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "readerr") {
        let (label, input, class, message) = (row[0], row[1], row[2], row[3]);
        let read = Tags::read(&unhex(input));
        if class == "none" {
            // An empty block is an empty list, not an error.
            let tags = read.unwrap_or_else(|e| panic!("{label}: {e}"));
            assert_eq!(
                message,
                format!("{} entries, first -", tags.len()),
                "{label}"
            );
        } else {
            let error = read.expect_err(label);
            assert_eq!(error.java_exception(), class, "{label}: exception");
            // `String.valueOf(null)` is the string "null", which is what an exception with no
            // message of its own leaves in the golden.
            let mine = error.message();
            assert_eq!(
                if mine.is_empty() { "null" } else { &mine },
                message,
                "{label}: message"
            );
        }
        compared += 1;
    }
    assert_eq!(compared, 10, "malformed blocks compared");
}

/// There is no in-memory `H`, so the text codec cannot write one back either.
#[test]
fn an_h_tag_re_encodes_as_the_array_it_decoded_into() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "textrt") {
        let (text, class, back, stable) = (row[0], row[1], row[2], row[3]);
        let (tag, value) = parse_tag(text).unwrap_or_else(|e| panic!("{text}: {e:?}"));
        assert_eq!(class_of(&value), class, "{text}: class");
        let mine = encode_tag(tag, &value).expect("encodes");
        assert_eq!(mine, back, "{text}: re-encoded");
        assert_eq!(
            if mine == text { "same" } else { "changed" },
            stable,
            "{text}: stability"
        );
        compared += 1;
    }
    assert_eq!(compared, 5, "text tags compared");
    // The two that change are the two H forms, and both become B arrays.
    let changed: Vec<&str> = rows(&corpus, "textrt")
        .into_iter()
        .filter(|row| row[3] == "changed")
        .map(|row| row[0])
        .collect();
    assert_eq!(changed, ["XX:H:48656C", "XX:H:"]);
}

/// The packed short, and the tag names htsjdk can pack but cannot name again.
#[test]
fn a_tag_packs_the_second_character_into_the_high_byte() {
    let corpus = corpus();
    let mut compared = 0;
    for row in rows(&corpus, "strtag") {
        let (name, packed, back) = (unescape(row[0]), row[1], row[2]);
        if packed == "-" {
            // A name that is not two characters, which `Tag::new` cannot be given at all.
            assert_ne!(name.encode_utf16().count(), 2, "{name}: length");
            compared += 1;
            continue;
        }
        let units: Vec<u16> = name.encode_utf16().collect();
        let tag = Tag::new(&[units[0] as u8, units[1] as u8]);
        assert_eq!(tag.0.to_string(), packed, "{name}: packed");

        if back == "ArrayIndexOutOfBoundsException" {
            // `makeStringTag` indexes `new String[Short.MAX_VALUE]` with the packed short, so a
            // tag whose second character is above 0x7F packs to a negative index and cannot be
            // named again. The port has no such cache and names it without complaint.
            assert!(tag.0 < 0, "{name}: the packed short is negative");
        } else {
            assert_eq!(escape(&tag.to_string()), back, "{name}: named back");
        }
        compared += 1;
    }
    assert_eq!(compared, 8, "tag names compared");
}
