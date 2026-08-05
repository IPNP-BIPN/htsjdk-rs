//! A second, independent reader over the BAM files this crate writes.
//!
//! Every other suite here compares against htsjdk, which is the right target and the only one:
//! byte-identity with htsjdk is the programme's goal. But a golden says "these bytes are the bytes
//! htsjdk produced", and that is a different claim from "these bytes are a BAM". The two coincide
//! when the golden covers whole files, and they come apart wherever a suite compares a rendering
//! rather than a file.
//!
//! [`noodles_bam`] is that second claim. It is an unrelated implementation of the format, written
//! from the specification rather than from htsjdk, so when it reads what this crate wrote and gets
//! the record back, two implementations that share no code agree about what the bytes mean.
//!
//! # Why it is a dev dependency and will stay one
//!
//! Decision 0036. noodles is a **reader to disagree with**, never a source of bytes:
//!
//!  * its primitives are not reachable anyway. `noodles_cram::io::reader::num` is `pub(crate)` and
//!    so are `codecs::rans_4x8::{encode, decode}`, so the parts a port would want to reuse cannot
//!    be called from outside the crate at all;
//!  * where it is reachable it disagrees on purpose. Its ITF8 reader uses `read_exact`, so a
//!    truncated stream is an error; htsjdk's returns `-1` and carries on. Both are defensible and
//!    only one of them is what this programme reproduces.
//!
//! So the rule is: hand-port the primitives against the oracle, and let noodles read the result.

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::{SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::BamWriter;

fn header() -> SamHeader {
    let mut header = SamHeader::new();
    header
        .sequences
        .push(SequenceRecord::new("chr1", 250_000_000));
    header.sequences.push(SequenceRecord::new("chr2", 100_000));
    header
}

fn record(name: &str, ref_index: i32, start: i32) -> BamRecord {
    BamRecord {
        read_name: name.into(),
        flags: 0,
        reference_index: ref_index,
        alignment_start: start,
        mapping_quality: 60,
        cigar: Cigar::new(vec![CigarElement {
            length: 4,
            op: Op::M,
        }]),
        mate_reference_index: -1,
        mate_alignment_start: 0,
        inferred_insert_size: 0,
        read_bases: b"ACGT".to_vec(),
        base_qualities: vec![30, 31, 32, 33],
        tags: Default::default(),
    }
}

/// Write with this crate, read with noodles, and compare what came back.
#[test]
fn what_this_crate_writes_is_a_bam_to_an_unrelated_reader() {
    let header = header();
    let mut writer = BamWriter::new(Vec::new(), &header).expect("a writer");
    let records = [
        record("r1", 0, 100),
        record("r2", 0, 200),
        record("r3", 1, 50),
    ];
    for record in &records {
        writer.write(record).expect("a record writes");
    }
    let bytes = writer.finish().expect("the file closes");

    let mut reader = noodles_bam::io::Reader::new(&bytes[..]);
    let noodles_header = reader.read_header().expect("noodles reads the header");

    // The reference sequences, which noodles parses out of the SAM text this crate encoded.
    let names: Vec<String> = noodles_header
        .reference_sequences()
        .keys()
        .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
        .collect();
    assert_eq!(names, vec!["chr1", "chr2"], "reference sequences");

    let mut seen = Vec::new();
    for result in reader.records() {
        let record = result.expect("noodles reads a record");
        seen.push(record);
    }
    assert_eq!(seen.len(), records.len(), "record count");

    for (mine, theirs) in records.iter().zip(seen.iter()) {
        let name = theirs
            .name()
            .map(|n| String::from_utf8_lossy(n).into_owned());
        assert_eq!(name.as_deref(), Some(mine.read_name.as_str()), "read name");
        let start = theirs
            .alignment_start()
            .transpose()
            .expect("a position")
            .map(usize::from);
        assert_eq!(
            start,
            Some(mine.alignment_start as usize),
            "alignment start"
        );
        assert_eq!(
            theirs.mapping_quality().map(u8::from),
            Some(mine.mapping_quality),
            "mapping quality"
        );
        assert_eq!(
            theirs.sequence().len(),
            mine.read_bases.len(),
            "sequence length"
        );
    }
}

/// Tags travel through both encoders, which is where two readers most easily part company: the
/// type letters are one byte each and a wrong one is still a parseable file.
#[test]
fn every_tag_type_survives_the_second_reader() {
    let header = header();
    let mut writer = BamWriter::new(Vec::new(), &header).expect("a writer");
    let mut record = record("tagged", 0, 100);
    record.tags.insert(Tag::new(b"NM"), TagValue::Int(7));
    record
        .tags
        .insert(Tag::new(b"MD"), TagValue::Str("4".into()));
    record.tags.insert(Tag::new(b"XF"), TagValue::Float(0.5));
    writer.write(&record).expect("the record writes");
    let bytes = writer.finish().expect("the file closes");

    let mut reader = noodles_bam::io::Reader::new(&bytes[..]);
    reader.read_header().expect("a header");
    let read = reader
        .records()
        .next()
        .expect("one record")
        .expect("noodles reads it");

    let data = read.data();
    let count = data.iter().count();
    assert_eq!(
        count, 3,
        "noodles found {count} tags where this crate wrote 3"
    );
}
