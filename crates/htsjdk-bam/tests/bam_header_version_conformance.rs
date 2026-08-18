//! Conformance for which `@HD` version a written BAM carries, against htsjdk.
//!
//! Golden from `tools/bam-conformance/BamHeaderVersionDump.java` in the pinned oracle container,
//! with the JDK deflater pinned per the oracle contract.
//!
//! # What this suite is for
//!
//! Two header paths that no other fixture can tell apart, because every other fixture builds its
//! header with `new SAMFileHeader()` and is therefore already at the current version.
//!
//!  * **the ordinary writer replaces the version**, `SAMFileWriterImpl.writeHeader` encoding with
//!    `keepExistingVersionNumber = false`, so a BAM written from a `VN:1.5` header says `VN:1.6`;
//!  * **the reheader path keeps it**, the static `BAMFileWriter.writeHeader` passing true;
//!  * **the writer stamps the sort order**, `setHeader` writing `unsorted` onto a header that
//!    carries none, which is why a written file has an `SO` its input did not;
//!  * **and `VN` cannot be anywhere but first**, even parsed from text that puts `SO` before it.
//!    That row is why the attribute order needs no fixing.

use std::io::Read;

use htsjdk_bam::header::{SamHeader, SequenceRecord};
use htsjdk_bam::writer::BamWriter;

fn golden_text() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/bam_header_version.txt.gz");
    let file = std::fs::File::open(&path).expect("golden corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden corpus is gzip");
    text
}

fn golden(kind: &str, name: &str) -> String {
    let prefix = format!("{kind}\t{name}\t");
    golden_text()
        .lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{name}"))
        .to_string()
}

fn minimal() -> SamHeader {
    let mut header = SamHeader::new();
    header
        .sequences
        .push(SequenceRecord::new("chr1", 250_000_000));
    header
}

/// The four headers the dump builds, in its order.
fn cases() -> Vec<(&'static str, SamHeader)> {
    let old = {
        let mut header = minimal();
        header.attributes.set("VN", "1.5");
        header
    };
    // The parsed case: `@HD SO:coordinate VN:1.5`, which htsjdk comes back from with VN first,
    // because the constructor set it before the text was read.
    let parsed = {
        let mut header = minimal();
        header.set_sort_order("coordinate");
        header.attributes.set("VN", "1.5");
        header
    };
    let sorted = {
        let mut header = minimal();
        header.attributes.set("VN", "1.4");
        header.set_sort_order("coordinate");
        header
    };
    vec![
        ("current", minimal()),
        ("old", old),
        ("parsed", parsed),
        ("sorted", sorted),
    ]
}

#[test]
fn the_reheader_path_keeps_the_version() {
    for (name, header) in cases() {
        let expected = golden("header", &format!("{name}_kept"));
        assert_eq!(
            header.encode().replace('\n', "\\n"),
            expected,
            "{name}: encoded with keepExistingVersionNumber = true"
        );
    }
}

#[test]
fn a_written_bam_carries_the_current_version() {
    for (name, header) in cases() {
        let expected = golden("file", &format!("{name}_written"));
        let writer = BamWriter::new(Vec::new(), &header).expect("a writer");
        let bytes = writer.finish().expect("a complete file");
        let hex: String = bytes.iter().map(|b| format!("{b:02x}")).collect();
        assert_eq!(hex, expected, "{name}: the whole file");
    }
}
