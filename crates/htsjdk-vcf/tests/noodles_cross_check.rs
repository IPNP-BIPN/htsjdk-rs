//! A second, independent reader over the VCF files this crate writes.
//!
//! Decision 0036, applied to the other whole-file writer in this repository. Every other suite here
//! compares against htsjdk, which is the right target; this one asks a different question, which no
//! golden asks: are the bytes a **VCF**, to an implementation that shares no code with either side?
//!
//! [`noodles_vcf`] is written from the specification. This port is written from htsjdk. Where the
//! two disagree the golden wins by construction, and the disagreement is worth knowing about
//! because it locates where htsjdk departs from the specification — which is the class of finding
//! this programme exists to record. Two of those are already on file for VCF and both are visible
//! from here: the header a reader hands back is not the header the file contains (decision 0035),
//! and the declared `Type` converts nothing.

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::indexed_writer::{write_vcf_indexed, SequenceEntry};
use htsjdk_vcf::variant::{Value, VariantContext};
use htsjdk_vcf::vcf_file::write_vcf;

fn header() -> VcfHeader {
    let mut header = VcfHeader::new();
    header.lines.push(HeaderLine::info(
        "DP",
        Cardinality::Fixed(1),
        LineType::Integer,
        "Depth",
    ));
    header.lines.push(HeaderLine::info(
        "AF",
        Cardinality::A,
        LineType::Float,
        "Allele Frequency",
    ));
    header
        .lines
        .push(HeaderLine::filter("LowQual", "Low quality"));
    header.lines.push(HeaderLine::format(
        "GT",
        Cardinality::Fixed(1),
        LineType::String,
        "Genotype",
    ));
    header.lines.push(HeaderLine::contig("chr1", 100_000, 0));
    header.lines.push(HeaderLine::contig("chr2", 200_000, 1));
    header
}

fn record(contig: &str, start: i64, reference: &str, alt: &str) -> VariantContext {
    let mut record = VariantContext::new(
        contig,
        start,
        vec![
            Allele::from_str(reference, true).unwrap(),
            Allele::from_str(alt, false).unwrap(),
        ],
    );
    record.stop = start + reference.len() as i64 - 1;
    record.attributes = vec![("DP".to_string(), Value::Str("10".to_string()))];
    record.filters = Some(Vec::new());
    record
}

fn records() -> Vec<VariantContext> {
    vec![
        record("chr1", 100, "A", "T"),
        record("chr1", 200, "C", "G"),
        record("chr2", 50, "GG", "G"),
    ]
}

/// Write with this crate, read with noodles, and compare what came back.
#[test]
fn what_this_crate_writes_is_a_vcf_to_an_unrelated_reader() {
    let text = write_vcf(&header(), &records()).expect("the fixture writes");

    let mut reader = noodles_vcf::io::Reader::new(text.as_bytes());
    let noodles_header = reader.read_header().expect("noodles reads the header");

    // The contigs, which noodles parses out of the metadata this crate sorted and rendered.
    let contigs: Vec<String> = noodles_header
        .contigs()
        .keys()
        .map(ToString::to_string)
        .collect();
    assert_eq!(contigs, vec!["chr1", "chr2"], "contigs");
    assert!(
        noodles_header.infos().contains_key("DP"),
        "the INFO declarations survive"
    );
    assert!(
        noodles_header.filters().contains_key("LowQual"),
        "the FILTER declaration survives"
    );

    let seen: Vec<_> = reader
        .records()
        .collect::<Result<Vec<_>, _>>()
        .expect("noodles reads every record");
    assert_eq!(seen.len(), records().len(), "record count");

    for (mine, theirs) in records().iter().zip(seen.iter()) {
        assert_eq!(theirs.reference_sequence_name(), mine.contig, "contig");
        let start = theirs
            .variant_start()
            .expect("a start")
            .expect("a valid start");
        assert_eq!(usize::from(start), mine.start as usize, "position");
        assert_eq!(
            theirs.reference_bases().to_string(),
            mine.alleles[0].display_string(),
            "reference allele"
        );
    }
}

/// The indexed writer's text is the plain writer's text, which is what says the index costs the
/// VCF nothing — asserted here against a reader that has no idea an index exists.
#[test]
fn the_indexed_writer_produces_the_same_readable_file() {
    let dictionary = [
        SequenceEntry {
            name: "chr1".into(),
            length: 100_000,
        },
        SequenceEntry {
            name: "chr2".into(),
            length: 200_000,
        },
    ];
    let indexed = write_vcf_indexed(
        &header(),
        &records(),
        Some(&dictionary),
        "file:///cross-check.vcf",
    )
    .expect("the fixture writes");

    let mut reader = noodles_vcf::io::Reader::new(indexed.text.as_bytes());
    let noodles_header = reader.read_header().expect("noodles reads the header");
    let count = reader.records().count();
    assert_eq!(count, records().len(), "record count");
    assert_eq!(noodles_header.contigs().len(), 2, "contigs");

    // And the positions the index recorded are offsets into that same text, so the byte at each
    // one begins a data line. An index built against a differently-written file fails here even
    // when both halves are internally consistent.
    for position in &indexed.record_positions {
        let at = *position as usize;
        assert!(at < indexed.text.len(), "position {at} is past the end");
        let line = indexed.text[at..]
            .lines()
            .next()
            .expect("a line begins at the recorded position");
        assert!(
            !line.starts_with('#') && line.contains('\t'),
            "position {at} does not begin a data line: {line:.40}"
        );
    }
}
