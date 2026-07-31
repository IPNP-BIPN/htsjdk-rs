//! Conformance for `IntervalListCodec`, against the oracle.
//!
//! Golden from `tools/tribble-conformance/IntervalListCodecDump.java`.
//!
//! The golden corrected nothing this time, which is worth saying because the last codec it was
//! asked about corrected five things. What it did settle is that four different malformed lines
//! fail in four different ways, and only one of them is a bad line to the reference:
//!
//! ```text
//! interval  chr1\t1\t10\t+\t     two-contigs  E:...TribbleException:Invalid interval record contains 4 fields: ...
//! interval  chr1\t1\t10\t.\tname two-contigs  E:java.lang.IllegalArgumentException:Invalid strand field: .
//! interval  chrX\t1\t10\t+\tname two-contigs  null
//! interval  chr1\t1\t201\t+\tname two-contigs E:java.lang.IllegalArgumentException:interval with end: 201 ...
//! ```
//!
//! A trailing tab costs a field rather than emptying the name; a dot is a strand everywhere else
//! in Tribble and an error here; an unknown contig is dropped and the file still loads; and an
//! interval one base past a contig kills the read. The empty dictionary answers `null` to every
//! well-formed line, so a file whose header lists no contig loads as no intervals rather than as
//! an error.

use std::io::Read;

use htsjdk_bam::header::{SamHeader, SequenceRecord};
use htsjdk_tribble::interval_list::{
    can_decode, decode, IntervalListError, IntervalRecord, Strand,
};

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/interval_list_codec.txt.gz");
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

/// The dictionaries the dump decoded against, by the label it wrote.
fn dictionary(label: &str) -> Option<SamHeader> {
    match label {
        "two-contigs" => {
            let mut header = SamHeader::default();
            for (name, length) in [("chr1", 200), ("chr2", 200), ("chr3", 0)] {
                header.sequences.push(SequenceRecord::new(name, length));
            }
            Some(header)
        }
        "empty" => Some(SamHeader::default()),
        "null" => None,
        other => panic!("unknown dictionary {other}"),
    }
}

fn show(record: &Option<IntervalRecord>) -> String {
    let Some(record) = record else {
        return "null".to_string();
    };
    format!(
        "{}:{}-{}|{}|{}",
        escape(&record.contig),
        record.start,
        record.end,
        match record.strand {
            Strand::Positive => "+",
            Strand::Negative => "-",
        },
        escape(&record.name)
    )
}

fn rendered(result: Result<Option<IntervalRecord>, IntervalListError>) -> String {
    match result {
        Ok(record) => show(&record),
        Err(error) => format!("E:{}:{}", error.class(), error.message()),
    }
}

#[test]
fn every_line_decodes_as_the_reference_decodes_it() {
    let text = golden();
    let mut count = 0;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("interval\t") else {
            continue;
        };
        let mut fields = rest.splitn(3, '\t');
        let escaped = fields.next().expect("a line");
        let label = fields.next().expect("a dictionary");
        let expected = fields.next().expect("a result");
        let header = dictionary(label);
        let ours = rendered(decode(&unescape(escaped), header.as_ref()));
        assert_eq!(ours, expected, "decoding {escaped:?} against {label}");
        count += 1;
    }
    assert!(count > 0, "the golden carries no interval rows");
    println!("{count} decodes identical");
}

#[test]
fn every_path_decodes_or_does_not_as_the_reference_says() {
    let text = golden();
    let mut count = 0;
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
        count += 1;
    }
    assert!(count > 0, "the golden carries no candecode rows");
}
