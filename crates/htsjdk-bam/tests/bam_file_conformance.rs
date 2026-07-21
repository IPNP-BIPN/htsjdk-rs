//! Conformance for whole BAM files and SAM text headers, against htsjdk.
//!
//! This is the first test in the port that compares a **complete file**. It composes three
//! things that were verified separately: the BGZF writer, the record codec, and the header
//! encoder. Composition is where a port that passes every unit test still produces the wrong
//! file, because framing errors live in the seams rather than in the parts.
//!
//! Goldens from `tools/bam-conformance/BamFileDump.java` in the pinned oracle container, with
//! the JDK deflater pinned per the oracle contract.

use std::io::Read;

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::{ProgramRecord, ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_bam::writer::BamWriter;

fn golden_text() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bam_file.txt.gz");
    let f = std::fs::File::open(&p).expect("golden corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("golden corpus is gzip");
    s
}

/// `kind -> name -> payload`.
fn goldens(kind: &str) -> Vec<(String, String)> {
    golden_text()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.splitn(3, '\t');
            let k = it.next()?;
            let name = it.next()?;
            let payload = it.next()?;
            (k == kind).then(|| (name.to_string(), payload.to_string()))
        })
        .collect()
}

fn tag(name: &str) -> Tag {
    Tag::new(name.as_bytes().try_into().unwrap())
}

fn minimal() -> SamHeader {
    let mut h = SamHeader::new();
    h.sequences.push(SequenceRecord::new("chr1", 250_000_000));
    h
}

/// The header cases, in the order the Java harness emits them.
fn header_cases() -> Vec<(String, SamHeader)> {
    let mut out: Vec<(String, SamHeader)> = Vec::new();

    out.push(("empty_dict".into(), SamHeader::new()));
    out.push(("minimal".into(), minimal()));

    let mut sorted = minimal();
    sorted.set_sort_order("coordinate");
    out.push(("sorted".into(), sorted));

    let mut queryname = minimal();
    queryname.set_sort_order("queryname");
    queryname.set_group_order("query");
    out.push(("queryname_grouped".into(), queryname));

    let many_sq = {
        let mut h = SamHeader::new();
        let mut s1 = SequenceRecord::new("chr1", 250_000_000);
        s1.attributes.set("AS", "GRCh38");
        s1.attributes.set("M5", "d41d8cd98f00b204e9800998ecf8427e");
        s1.attributes.set("SP", "Homo sapiens");
        s1.attributes.set("UR", "file:/ref/chr1.fa");
        h.sequences.push(s1);
        h.sequences.push(SequenceRecord::new("chr2", 200_000_000));
        let mut s3 = SequenceRecord::new("chrM", 16571);
        s3.attributes.set("M5", "c68f52674c9fb33aef52dcf399755519");
        h.sequences.push(s3);
        h
    };
    out.push(("many_sq".into(), many_sq.clone()));

    let rg1 = {
        let mut rg = ReadGroup::new("rg1");
        rg.attributes.set("SM", "sample1");
        rg.attributes.set("LB", "lib1");
        rg.attributes.set("PL", "ILLUMINA");
        rg.attributes.set("PU", "unit1");
        rg
    };
    let mut with_rg = minimal();
    with_rg.read_groups.push(rg1.clone());
    let mut rg2 = ReadGroup::new("rg2");
    rg2.attributes.set("SM", "sample2");
    with_rg.read_groups.push(rg2);
    out.push(("read_groups".into(), with_rg));

    // The LinkedHashMap property: SM is overwritten after LB was set, and must stay first.
    let mut reorder = minimal();
    let mut rg = ReadGroup::new("rgx");
    rg.attributes.set("SM", "first");
    rg.attributes.set("LB", "lib");
    rg.attributes.set("SM", "second");
    rg.attributes.set("PL", "ILLUMINA");
    reorder.read_groups.push(rg);
    out.push(("attribute_overwrite_keeps_position".into(), reorder));

    let pg1 = {
        let mut pg = ProgramRecord::new("prog1");
        pg.attributes.set("PN", "MarkDuplicates");
        pg.attributes.set("VN", "3.4.0");
        pg.attributes.set("CL", "MarkDuplicates I=in.bam O=out.bam");
        pg
    };
    let mut with_pg = minimal();
    with_pg.programs.push(pg1.clone());
    let mut pg2 = ProgramRecord::new("prog2");
    pg2.attributes.set("PP", "prog1");
    with_pg.programs.push(pg2);
    with_pg.add_comment("a comment");
    with_pg.add_comment("another comment");
    out.push(("programs_and_comments".into(), with_pg));

    let mut full = many_sq;
    full.set_sort_order("coordinate");
    full.read_groups.push(rg1);
    full.programs.push(pg1);
    full.add_comment("full header");
    out.push(("full".into(), full));

    out
}

#[test]
fn every_header_encodes_to_htsjdks_text() {
    let cases = header_cases();
    let goldens = goldens("header");
    assert_eq!(cases.len(), goldens.len(), "case lists have drifted");

    let mut failures = Vec::new();
    for ((name, header), (golden_name, golden)) in cases.iter().zip(&goldens) {
        assert_eq!(name, golden_name, "case lists out of step");
        // The harness escapes newlines so each case stays on one line.
        let expected = golden.replace("\\n", "\n");
        let ours = header.encode();
        if ours != expected {
            failures.push(format!(
                "{name}:\n  ours:   {:?}\n  htsjdk: {:?}",
                ours, expected
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} headers diverge:\n{}",
        failures.len(),
        cases.len(),
        failures.join("\n")
    );
}

/// The one that would be lost by any sorted or append-on-overwrite attribute map.
#[test]
fn the_overwrite_case_is_actually_in_the_goldens() {
    let g = goldens("header");
    let (_, text) = g
        .iter()
        .find(|(n, _)| n == "attribute_overwrite_keeps_position")
        .expect("the overwrite case must be in the corpus");
    assert!(
        // Only newlines are escaped by the harness; tabs are literal.
        text.contains("ID:rgx\tSM:second\tLB:lib\tPL:ILLUMINA"),
        "htsjdk keeps an overwritten key in its original position; golden says: {text}"
    );
}

fn file_header() -> SamHeader {
    let mut h = minimal();
    h.set_sort_order("coordinate");
    h
}

fn record(name: &str, start: i32, cigar_len: u32, bases: Vec<u8>, quals: Vec<u8>) -> BamRecord {
    BamRecord {
        read_name: name.into(),
        flags: 0,
        reference_index: 0,
        alignment_start: start,
        mapping_quality: 60,
        cigar: Cigar::new(vec![CigarElement {
            length: cigar_len,
            op: Op::M,
        }]),
        mate_reference_index: -1,
        mate_alignment_start: 0,
        inferred_insert_size: 0,
        read_bases: bases,
        base_qualities: quals,
        tags: Default::default(),
    }
}

/// The file cases, in harness order.
fn file_cases() -> Vec<(String, SamHeader, Vec<BamRecord>)> {
    let mut out = Vec::new();

    let mut sorted = minimal();
    sorted.set_sort_order("coordinate");
    out.push(("file_empty".to_string(), sorted, Vec::new()));

    let mut r1 = record("read1", 100, 4, b"ACGT".to_vec(), vec![30, 30, 30, 30]);
    r1.tags.insert(tag("NM"), TagValue::Int(0));
    out.push(("file_one_record".to_string(), file_header(), vec![r1]));

    let many: Vec<BamRecord> = (0..500)
        .map(|i| {
            let bases: Vec<u8> = (0..10).map(|j| b"ACGT"[((i + j) % 4) as usize]).collect();
            let mut r = record(
                &format!("read{i}"),
                100 + i * 37,
                10,
                bases,
                vec![(20 + i % 20) as u8; 10],
            );
            r.tags.insert(tag("NM"), TagValue::Int((i % 5) as i64));
            r.tags.insert(tag("RG"), TagValue::Str("rg1".into()));
            r
        })
        .collect();
    out.push(("file_500_records".to_string(), file_header(), many));

    let lots: Vec<BamRecord> = (0..20000i32)
        .map(|i| {
            let bases: Vec<u8> = (0..50)
                .map(|j| b"ACGTN"[((i * 7 + j) % 5) as usize])
                .collect();
            let quals: Vec<u8> = (0..50).map(|j| ((i + j) % 60) as u8).collect();
            record(&format!("r{i}"), 1 + i * 11, 50, bases, quals)
        })
        .collect();
    out.push(("file_20000_records".to_string(), file_header(), lots));

    out
}

/// The gate for the whole write path: header, dictionary, records and BGZF framing together.
#[test]
fn every_file_is_byte_identical_to_htsjdks() {
    let cases = file_cases();
    let goldens = goldens("file");
    assert_eq!(cases.len(), goldens.len(), "case lists have drifted");

    let mut failures = Vec::new();
    for ((name, header, records), (golden_name, hex)) in cases.iter().zip(&goldens) {
        assert_eq!(name, golden_name, "case lists out of step");
        let expected: Vec<u8> = (0..hex.len() / 2)
            .map(|i| u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16).unwrap())
            .collect();

        let mut w = BamWriter::new(Vec::new(), header).unwrap();
        for r in records {
            w.write(r).unwrap();
        }
        let ours = w.finish().unwrap();

        if ours != expected {
            // Locate the divergence the way the differential harness would: in the BGZF
            // structure, not just as a byte offset.
            let at = ours
                .iter()
                .zip(&expected)
                .position(|(a, b)| a != b)
                .unwrap_or(ours.len().min(expected.len()));
            failures.push(format!(
                "{name}: {} vs {} bytes, first differs at {at}",
                ours.len(),
                expected.len()
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} files diverge from htsjdk:\n{}",
        failures.len(),
        cases.len(),
        failures.join("\n")
    );
}

/// The 20,000-record case must actually span several BGZF blocks, or it is not testing the
/// framing it was added to test.
#[test]
fn the_large_case_spans_multiple_bgzf_blocks() {
    let (_, hex) = goldens("file")
        .into_iter()
        .find(|(n, _)| n == "file_20000_records")
        .expect("the large case must exist");
    let bytes: Vec<u8> = (0..hex.len() / 2)
        .map(|i| u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16).unwrap())
        .collect();

    let mut blocks = 0;
    let mut off = 0usize;
    while off + 18 <= bytes.len() {
        let bsize = u16::from_le_bytes([bytes[off + 16], bytes[off + 17]]) as usize + 1;
        if bsize < 18 || off + bsize > bytes.len() {
            break;
        }
        blocks += 1;
        off += bsize;
    }
    assert!(
        blocks >= 3,
        "the large case must cross block boundaries to test the framing; got {blocks} blocks"
    );
}
