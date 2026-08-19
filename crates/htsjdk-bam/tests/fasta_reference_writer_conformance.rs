//! Conformance for `FastaReferenceWriter` against htsjdk 4.2.0.
//!
//! Golden from `tools/bam-conformance/FastaReferenceWriterDump.java` in the pinned oracle
//! container.
//!
//! # What this suite is for
//!
//!  * **the index offset counts bytes, not bases**: `chr2` in the two-sequence case starts at 58,
//!    which is all of chr1 including its header and its newlines plus chr2's own header;
//!  * **bytes-per-line is bases-per-line plus one**, whatever the sequence's length;
//!  * **one newline per sequence, at the end**, so a length that is a multiple of the width gets no
//!    blank line;
//!  * **chunked appends are the same bytes**, because the line breaks come from a running count;
//!  * **the md5 is of the upper-cased bases**, while the FASTA keeps its case;
//!  * **and the nine refusals each have their own class and message**, with a tab in a description
//!    allowed where every other control character is not.

use htsjdk_bam::fasta_writer::{FastaOutputs, FastaReferenceWriter, FastaWriterError};
use std::io::Read;

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/fasta_reference_writer.txt.gz");
    let file = std::fs::File::open(&path).expect("golden corpus");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden corpus is gzip");
    text
}

fn row(text: &str, kind: &str, label: &str) -> String {
    let prefix = format!("{kind}\t{label}\t");
    text.lines()
        .find_map(|line| line.strip_prefix(prefix.as_str()))
        .unwrap_or_else(|| panic!("the golden carries {kind}/{label}"))
        .to_string()
}

/// The dump's escaping, so the port's bytes can be compared against the golden's text.
fn escape(text: &str) -> String {
    text.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\n', "\\n")
        .replace('\u{1}', "\\u0001")
}

fn check(text: &str, label: &str, outputs: &FastaOutputs) {
    assert_eq!(
        escape(&String::from_utf8_lossy(&outputs.fasta)),
        row(text, "fasta", label),
        "{label}: the FASTA"
    );
    assert_eq!(
        escape(&outputs.index),
        row(text, "fai", label),
        "{label}: the index"
    );
    assert_eq!(
        escape(&outputs.dictionary),
        row(text, "dict", label),
        "{label}: the dictionary"
    );
}

fn repeat(unit: &str, times: usize) -> Vec<u8> {
    unit.repeat(times).into_bytes()
}

#[test]
fn every_written_reference_matches_the_golden() {
    let text = golden();

    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    writer.append_bases(&repeat("ACGT", 25)).expect("bases");
    check(&text, "default-width", &writer.close().expect("closed"));

    let mut writer = FastaReferenceWriter::new(10, false).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    writer.append_bases(&repeat("ACGT", 10)).expect("bases");
    check(&text, "exact-multiple", &writer.close().expect("closed"));

    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    writer.append_bases(b"ACGT").expect("bases");
    check(&text, "short", &writer.close().expect("closed"));

    let mut writer = FastaReferenceWriter::new(12, false).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    writer.append_bases(&repeat("ACGT", 7)).expect("bases");
    writer
        .start_sequence_with("chr2", "the second one", 5)
        .expect("a sequence");
    writer.append_bases(&repeat("TTGC", 3)).expect("bases");
    check(&text, "two-sequences", &writer.close().expect("closed"));

    // The same sequence in three pieces that do not line up with the width.
    let mut writer = FastaReferenceWriter::new(12, false).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    writer.append_bases(b"ACG").expect("bases");
    writer.append_bases(b"TACGTACG").expect("bases");
    writer.append_bases(b"TACGTACGTACGTACGT").expect("bases");
    check(&text, "chunked", &writer.close().expect("closed"));

    let mut writer = FastaReferenceWriter::new(10, false).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    writer.append_bases(b"acgtRYKMSWacgtNNNN").expect("bases");
    check(&text, "mixed-case", &writer.close().expect("closed"));

    let mut writer = FastaReferenceWriter::new(10, true).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    writer.append_bases(b"acgtacgtac").expect("bases");
    check(&text, "md5", &writer.close().expect("closed"));

    let mut writer = FastaReferenceWriter::new(10, true).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    writer.append_bases(b"ACGTACGTAC").expect("bases");
    check(&text, "md5-uppercase", &writer.close().expect("closed"));
}

/// The chunked case and the two-sequence case share their first sequence's bytes exactly, which is
/// the claim that the line breaks do not follow the calls.
#[test]
fn chunked_appends_write_the_same_bytes() {
    let text = golden();
    let chunked = row(&text, "fasta", "chunked");
    let two = row(&text, "fasta", "two-sequences");
    assert!(
        two.starts_with(chunked.trim_end_matches("\\n")),
        "the first sequence of the two-sequence case is the chunked one"
    );
}

fn error_of(result: Result<(), FastaWriterError>) -> FastaWriterError {
    result.expect_err("the reference refuses this")
}

#[test]
fn every_refusal_matches_the_golden() {
    let text = golden();

    let expect = |label: &str, error: FastaWriterError| {
        let expected = row(&text, "error", label);
        let mine = format!("{}:{}", error.java_class(), escape(&error.message()));
        assert_eq!(mine, expected, "{label}");
    };

    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    expect("empty-name", error_of(writer.start_sequence("")));

    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    expect("blank-in-name", error_of(writer.start_sequence("chr 1")));

    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    expect(
        "control-in-name",
        error_of(writer.start_sequence("chr\u{1}")),
    );

    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    expect(
        "control-in-description",
        error_of(writer.start_sequence_with("chr1", "a\u{1}b", 60)),
    );

    // A sequence opened and left empty, which the next startSequence refuses.
    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    expect("no-bases", error_of(writer.start_sequence("chr2")));

    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    writer.append_bases(b"ACGT").expect("bases");
    expect("duplicate-name", error_of(writer.start_sequence("chr1")));

    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    expect("bad-base", error_of(writer.append_bases(b"ACGZ")));

    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    expect(
        "bases-before-sequence",
        error_of(writer.append_bases(b"ACGT")),
    );

    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    expect(
        "zero-width",
        error_of(writer.start_sequence_with("chr1", "", 0)),
    );

    // And closing with a sequence open and empty, which is the same refusal from the other side.
    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    writer.start_sequence("chr1").expect("a sequence");
    let error = writer.close().expect_err("no base was added");
    let expected = row(&text, "error", "close-with-no-bases");
    assert_eq!(
        format!("{}:{}", error.java_class(), escape(&error.message())),
        expected
    );
}

/// A tab is the one control character a description may hold, and the golden carries the file it
/// produces rather than an error.
#[test]
fn a_tab_in_a_description_is_written_rather_than_refused() {
    let text = golden();
    let mut writer = FastaReferenceWriter::new(60, false).expect("a writer");
    writer
        .start_sequence_with("chr1", "a\tb", 60)
        .expect("a tab is allowed");
    writer.append_bases(b"ACGT").expect("bases");
    check(
        &text,
        "tab-in-description",
        &writer.close().expect("closed"),
    );
}
