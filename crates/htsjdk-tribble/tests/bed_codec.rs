//! Conformance for `BEDCodec`, against the oracle.
//!
//! Golden from `tools/tribble-conformance/BedCodecDump.java`.
//!
//! The golden corrected the port on five points, and every one of them is a default or an
//! exception rather than arithmetic:
//!
//! ```text
//! bed  chr1\t100\t200                       ONE  chr1:101-200||NaN|.|null|
//! bed  \strack\sname=whatever               ONE  E:java.lang.NumberFormatException:For input string: "track"
//! bed  ...\t300,0,0                         ONE  E:java.lang.IllegalArgumentException:Color parameter outside of expected range: Red
//! bed  ...\t3\t10,20\t0,50                  ONE  E:java.lang.NumberFormatException:Cannot parse null string
//! bed  ...\t0\t\t                           ONE  E:java.lang.ArrayIndexOutOfBoundsException:Index 0 out of bounds for length 0
//! ```
//!
//! An unnamed feature carries `""` and an unscored one carries `NaN`, so neither is absent; a
//! header prefix with a leading space is not a header, and the leading space is itself a
//! separator, so the line throws on `"track"` as a coordinate; a colour component of 300 fails the
//! whole line rather than clamping; a declared exon count larger than the lists parses a `null`;
//! and a declared count of zero indexes a zero-length array.

use std::io::Read;

use htsjdk_tribble::bed::{
    can_decode, decode, split_bed_line, BedError, BedFeature, StartOffset, Strand,
};

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bed_codec.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

/// The dump's `escape`: a tab is `\t`, a space is `\s`, everything else outside printable ASCII is
/// `\uXXXX`.
fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('s') => out.push(' '),
            Some('\\') => out.push('\\'),
            Some('u') => {
                let hex: String = (0..4).filter_map(|_| chars.next()).collect();
                let code = u32::from_str_radix(&hex, 16).expect("four hex digits");
                out.push(char::from_u32(code).expect("a character"));
            }
            other => panic!("unknown escape {other:?}"),
        }
    }
    out
}

fn escape(text: &str) -> String {
    text.chars()
        .map(|c| match c {
            '\t' => "\\t".to_string(),
            ' ' => "\\s".to_string(),
            '\\' => "\\\\".to_string(),
            c if (c as u32) < 0x20 || (c as u32) > 0x7e => format!("\\u{:04x}", c as u32),
            c => c.to_string(),
        })
        .collect()
}

/// `Float.toString`, in the two shapes this golden needs: an integral value prints with a trailing
/// `.0`, and the specials print by name.
fn java_float(value: f32) -> String {
    if value.is_nan() {
        return "NaN".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_string();
    }
    if value == value.trunc() && value.abs() < 1e7 {
        return format!("{value:.1}");
    }
    format!("{value}")
}

/// The dump's rendering of one feature.
fn show(feature: &Option<BedFeature>) -> String {
    let Some(feature) = feature else {
        return "null".to_string();
    };
    let color = match feature.color {
        None => "null".to_string(),
        Some((r, g, b)) => format!("{r},{g},{b}"),
    };
    let exons: Vec<String> = feature
        .exons
        .iter()
        .map(|exon| {
            format!(
                "{}-{}#{}#{}-{}",
                exon.start, exon.end, exon.number, exon.cd_start, exon.cd_end
            )
        })
        .collect();
    format!(
        "{}:{}-{}|{}|{}|{}|{}|{}",
        feature.contig,
        feature.start,
        feature.end,
        escape(&feature.name),
        java_float(feature.score),
        feature.strand.name(),
        color,
        exons.join(";")
    )
}

fn rendered(result: Result<Option<BedFeature>, BedError>) -> String {
    match result {
        Ok(feature) => show(&feature),
        Err(error) => format!("E:{}:{}", error.class(), error.message()),
    }
}

#[test]
fn every_line_decodes_as_the_reference_decodes_it() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("bed\t") else {
            continue;
        };
        let mut fields = rest.splitn(3, '\t');
        let escaped = fields.next().expect("a line");
        let offset = match fields.next().expect("an offset") {
            "ZERO" => StartOffset::Zero,
            "ONE" => StartOffset::One,
            other => panic!("unknown offset {other}"),
        };
        let expected = fields.next().expect("a result");
        let ours = rendered(decode(&unescape(escaped), offset));
        assert_eq!(ours, expected, "decoding {escaped:?} at {offset:?}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no bed rows");
    println!("{count} decodes identical");
}

#[test]
fn every_split_matches_the_reference() {
    let text = golden();
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("split\t") else {
            continue;
        };
        let mut fields = rest.splitn(3, '\t');
        let escaped = fields.next().expect("a line");
        let count: usize = fields.next().expect("a count").parse().expect("a number");
        let expected = fields.next().unwrap_or("");
        let line = unescape(escaped);
        let ours = split_bed_line(&line);
        assert_eq!(ours.len(), count, "field count for {escaped:?}");
        assert_eq!(
            ours.iter().map(|f| escape(f)).collect::<Vec<_>>().join("|"),
            expected,
            "fields for {escaped:?}"
        );
    }
}

#[test]
fn every_path_decodes_or_does_not_as_the_reference_says() {
    let text = golden();
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("candecode\t") else {
            continue;
        };
        let (path, expected) = rest.split_once('\t').expect("a path and an answer");
        assert_eq!(
            can_decode(path).to_string(),
            expected,
            "canDecode({path:?})"
        );
    }
}

/// The rows the golden corrected, kept as assertions so a later change cannot quietly undo them.
#[test]
fn the_five_the_golden_corrected() {
    let text = golden();
    let row = |line: &str, offset: &str| -> String {
        let needle = format!("bed\t{line}\t{offset}\t");
        text.lines()
            .find(|l| l.starts_with(&needle))
            .unwrap_or_else(|| panic!("no row for {line:?}"))[needle.len()..]
            .to_string()
    };

    // Name and score are defaults, not absences: "" and NaN.
    assert_eq!(row(r"chr1\t100\t200", "ONE"), "chr1:101-200||NaN|.|null|");
    // And the strand of a feature with no strand column prints as `.`, like an explicit one.
    assert_eq!(
        row(r"chr1\t100\t200\tname\t500\t.", "ONE"),
        row(r"chr1\t100\t200\tname\t500", "ONE").replace("|NaN|", "|500.0|")
    );

    // A leading space makes a header line data, and the space is itself a separator, so the
    // line throws on "track" as a coordinate rather than decoding to null.
    assert!(row(r"\strack\sname=whatever", "ONE").contains(r#"For input string: "track""#));
    assert_eq!(row("track\\sname=whatever", "ONE"), "null");

    // A colour component out of range fails the whole line; an unknown colour name is black.
    assert!(
        row(r"chr1\t100\t200\tname\t500\t+\t100\t200\t300,0,0", "ONE")
            .contains("Color parameter outside of expected range: Red")
    );
    assert!(row(r"chr1\t100\t200\tname\t500\t+\t100\t200\tnonsense", "ONE").contains("|0,0,0|"));

    // A declared exon count larger than the lists parses a null; a count of zero indexes a
    // zero-length array.
    assert!(row(
        r"chr1\t100\t200\tname\t500\t+\t100\t200\t255,0,0\t3\t10,20\t0,50",
        "ONE"
    )
    .contains("Cannot parse null string"));
    assert!(row(
        r"chr1\t100\t200\tname\t500\t+\t100\t200\t255,0,0\t0\t\t",
        "ONE"
    )
    .contains("Index 0 out of bounds for length 0"));

    // The exon numbering runs backwards on the negative strand, and the coding bounds are clamped
    // to each exon rather than to the feature.
    assert!(row(
        r"chr1\t100\t200\tname\t500\t-\t100\t200\t255,0,0\t2\t10,20\t0,50",
        "ONE"
    )
    .ends_with("101-110#2#101-110;151-170#1#151-170"));

    // A bad score returns the feature early, so the strand and colour columns beside it are lost.
    assert_eq!(
        row(
            r"chr1\t100\t200\tname\t.\t+\t100\t200\t255,0,0\t2\t10,10\t0,50",
            "ONE"
        ),
        "chr1:101-200|name|NaN|.|null|"
    );

    // Java's float parser takes a type suffix and the spelled-out infinity, and refuses `inf`,
    // which is where a Rust parser would have disagreed in both directions.
    assert!(row(r"chr1\t100\t200\tname\t1.5f\t+", "ONE").contains("|1.5|+|"));
    assert!(row(r"chr1\t100\t200\tname\tInfinity\t+", "ONE").contains("|Infinity|+|"));
    assert!(row(r"chr1\t100\t200\tname\tinf\t+", "ONE").contains("|NaN|.|"));

    // And the strand default is NONE, which shares its rendering with an explicit dot.
    assert_eq!(Strand::None.name(), ".");
}
