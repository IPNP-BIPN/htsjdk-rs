//! Conformance for whole SAM text files against htsjdk, and for the BAM → SAM conversion.
//!
//! Each case carries the SAM htsjdk wrote **and** the BAM it wrote from the same records. The
//! test that matters reads the BAM and writes the SAM, which is the conversion a user actually
//! performs, and compares against htsjdk's SAM byte for byte.
//!
//! That path exercises the reader, the header parser, the text writer and the tag re-derivation
//! together. None of the unit tests covers the seam between them.

use std::io::Read;

use htsjdk_bam::reader::BamReader;
use htsjdk_bam::sam_file::{read_sam_with, write_sam};
use htsjdk_bam::text_parse::ValidationStringency;

fn corpus() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sam_file.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s
}

fn rows(kind: &str) -> Vec<(String, String)> {
    corpus()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.splitn(3, '\t');
            let k = it.next()?;
            let name = it.next()?.to_string();
            let payload = it.next()?.to_string();
            (k == kind).then_some((name, payload))
        })
        .collect()
}

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

fn from_hex(hex: &str) -> Vec<u8> {
    (0..hex.len() / 2)
        .map(|i| u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16).unwrap())
        .collect()
}

/// The gate: convert htsjdk's BAM to SAM and get htsjdk's SAM.
#[test]
fn converting_htsjdks_bam_produces_htsjdks_sam() {
    let bams = rows("bam");
    let sams = rows("sam");
    assert_eq!(bams.len(), sams.len(), "corpus is inconsistent");
    assert!(bams.len() >= 5, "expected at least 5 cases");

    let mut failures = Vec::new();
    for ((name, bam_hex), (sam_name, sam_escaped)) in bams.iter().zip(&sams) {
        assert_eq!(name, sam_name, "case lists out of step");
        let expected = unescape(sam_escaped);

        let plain = htsjdk_bgzf::decompress_all(&from_hex(bam_hex)).unwrap();
        let reader = BamReader::new(&plain).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let header = reader.header.text.clone();
        let records: Vec<_> = reader
            .map(|r| r.unwrap_or_else(|e| panic!("{name}: {e:?}")))
            .collect();

        match write_sam(&header, &records) {
            None => failures.push(format!("{name}: refused to render")),
            Some(ours) if ours != expected => {
                let at = ours
                    .lines()
                    .zip(expected.lines())
                    .position(|(a, b)| a != b)
                    .unwrap_or(0);
                failures.push(format!(
                    "{name}: line {}\n  ours:   {:?}\n  htsjdk: {:?}",
                    at + 1,
                    ours.lines().nth(at).unwrap_or("<eof>"),
                    expected.lines().nth(at).unwrap_or("<eof>")
                ));
            }
            Some(_) => {}
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} conversions diverge:\n{}",
        failures.len(),
        bams.len(),
        failures.join("\n\n")
    );
}

/// And the reverse: htsjdk's SAM parses back to the records the BAM holds.
///
/// Read at `SILENT`, because one of the files htsjdk wrote cannot be read at the default
/// `STRICT`. See the comment inside.
#[test]
fn htsjdks_sam_parses_back_to_the_same_records() {
    let bams = rows("bam");
    let sams = rows("sam");

    for ((name, bam_hex), (_, sam_escaped)) in bams.iter().zip(&sams) {
        let plain = htsjdk_bgzf::decompress_all(&from_hex(bam_hex)).unwrap();
        let reader = BamReader::new(&plain).unwrap();
        let from_bam: Vec<_> = reader.map(|r| r.unwrap()).collect();

        // SILENT, and the reason is a finding rather than a convenience. The `unmapped` case
        // contains a record htsjdk's own writer produced and its own default-stringency reader
        // rejects: RNAME `*` with MAPQ 60. Reading htsjdk's output at STRICT is therefore
        // impossible for that file, which is exactly why the port models stringency at all.
        let (_, from_sam) = read_sam_with(&unescape(sam_escaped), ValidationStringency::Silent)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"));

        assert_eq!(
            from_sam.len(),
            from_bam.len(),
            "{name}: record count differs"
        );
        for (i, (a, b)) in from_sam.iter().zip(&from_bam).enumerate() {
            // The binary form re-derives every integer tag's width, so comparing the encoded
            // bytes is the strictest check available and the one that matters.
            assert_eq!(
                a.encode().unwrap(),
                b.encode().unwrap(),
                "{name}: record {i} differs after SAM round trip"
            );
        }
    }
}

/// The corpus must reach the sections and shapes the port claims to handle.
#[test]
fn the_corpus_covers_the_shapes_that_matter() {
    let names: Vec<String> = rows("sam").into_iter().map(|(n, _)| n).collect();
    for expected in [
        "header_only",
        "two_references",
        "full_header_and_tags",
        "unmapped",
        "many_records",
    ] {
        assert!(names.iter().any(|n| n == expected), "missing {expected}");
    }
    let (_, full) = rows("sam")
        .into_iter()
        .find(|(n, _)| n == "full_header_and_tags")
        .unwrap();
    for section in ["@HD", "@SQ", "@RG", "@PG", "@CO"] {
        assert!(
            full.contains(section),
            "the full header must contain {section}"
        );
    }
}
