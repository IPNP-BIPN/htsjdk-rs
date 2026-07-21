//! Conformance for the BAI index, against htsjdk's `BAMIndexer`.
//!
//! Goldens from `tools/bam-conformance/BaiDump.java` in the pinned oracle container.
//!
//! Each case carries both the BAM and its BAI. The BAM is checked first, deliberately: an index
//! that matches for a file that does not is a coincidence, not a result, because the index is
//! made of virtual file pointers into that exact byte stream.

use std::io::Read;

use htsjdk_bam::cigar::{Cigar, CigarElement, Op};
use htsjdk_bam::header::{SamHeader, SequenceRecord};
use htsjdk_bam::record::{BamRecord, READ_UNMAPPED_FLAG};
use htsjdk_bam::writer::BamWriter;

fn golden_text() -> String {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/bai.txt.gz");
    let f = std::fs::File::open(&p).expect("golden corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("golden corpus is gzip");
    s
}

fn goldens(kind: &str) -> Vec<(String, Vec<u8>)> {
    golden_text()
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| {
            let mut it = l.splitn(3, '\t');
            let k = it.next()?;
            let name = it.next()?.to_string();
            let hex = it.next()?;
            if k != kind {
                return None;
            }
            let bytes = (0..hex.len() / 2)
                .map(|i| u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16).unwrap())
                .collect();
            Some((name, bytes))
        })
        .collect()
}

fn header(lengths: &[i32]) -> SamHeader {
    let mut h = SamHeader::new();
    for (i, &len) in lengths.iter().enumerate() {
        h.sequences
            .push(SequenceRecord::new(&format!("chr{}", i + 1), len));
    }
    h.set_sort_order("coordinate");
    h
}

/// `TextCigarCodec.decode`, for the shapes the harness uses.
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
                b'S' => Op::S,
                _ => panic!("unexpected operator {}", b as char),
            };
            elements.push(CigarElement { length: len, op });
            len = 0;
        }
    }
    Cigar::new(elements)
}

fn read(name: &str, reference: i32, start: i32, cigar_text: &str) -> BamRecord {
    let c = cigar(cigar_text);
    let len = c.read_length() as usize;
    BamRecord {
        read_name: name.into(),
        flags: 0,
        reference_index: reference,
        alignment_start: start,
        mapping_quality: 60,
        cigar: c,
        mate_reference_index: -1,
        mate_alignment_start: 0,
        inferred_insert_size: 0,
        read_bases: vec![b'A'; len],
        base_qualities: vec![30; len],
        tags: Default::default(),
    }
}

/// The case list, in harness order.
///
/// Built by successive `push` so the order matches the Java harness line by line and stays
/// checkable against it; clippy would prefer one `vec![]` literal, which would make the
/// correspondence harder to verify and is not worth it here.
#[allow(clippy::vec_init_then_push)]
fn cases() -> Vec<(String, SamHeader, Vec<BamRecord>)> {
    let mut out: Vec<(String, SamHeader, Vec<BamRecord>)> = Vec::new();

    out.push(("empty".into(), header(&[250_000_000]), vec![]));

    out.push((
        "one_read".into(),
        header(&[250_000_000]),
        vec![read("r1", 0, 100, "4M")],
    ));

    out.push((
        "same_window".into(),
        header(&[250_000_000]),
        vec![read("r1", 0, 100, "50M"), read("r2", 0, 200, "50M")],
    ));

    out.push((
        "sparse_windows".into(),
        header(&[250_000_000]),
        vec![
            read("r1", 0, 1, "50M"),
            read("r2", 0, 5 * 16384 + 1, "50M"),
            read("r3", 0, 40 * 16384 + 1, "50M"),
        ],
    ));

    out.push((
        "all_levels".into(),
        header(&[250_000_000]),
        vec![
            read("lvl5", 0, 1, "100M"),
            read("lvl4", 0, 100, "200000M"),
            read("lvl3", 0, 300000, "2000000M"),
            read("lvl2", 0, 3000000, "20000000M"),
            read("lvl1", 0, 30000000, "100000000M"),
        ],
    ));

    out.push((
        "multi_reference_with_gap".into(),
        header(&[250_000_000, 200_000_000, 100_000_000]),
        vec![read("a", 0, 100, "50M"), read("b", 2, 500, "50M")],
    ));

    let unmapped = {
        let mut placed = read("placed", 0, 100, "50M");
        placed.flags |= READ_UNMAPPED_FLAG;
        let mut unplaced = read("unplaced", 0, 100, "50M");
        unplaced.flags |= READ_UNMAPPED_FLAG;
        unplaced.reference_index = -1;
        unplaced.alignment_start = 0;
        unplaced.cigar = Cigar::default();
        vec![
            read("mapped", 0, 50, "50M"),
            placed,
            unplaced.clone(),
            unplaced,
        ]
    };
    out.push(("unmapped".into(), header(&[250_000_000]), unmapped));

    let many: Vec<BamRecord> = (0..20000)
        .map(|i| read(&format!("r{i}"), 0, 1 + i * 11, "50M"))
        .collect();
    out.push(("many_blocks".into(), header(&[250_000_000]), many));

    out.push((
        "window_boundary".into(),
        header(&[250_000_000]),
        vec![
            read("before", 0, 16384, "1M"),
            read("on", 0, 16385, "1M"),
            read("after", 0, 16386, "1M"),
        ],
    ));

    out
}

fn build(header: &SamHeader, records: &[BamRecord]) -> (Vec<u8>, Vec<u8>) {
    let mut w = BamWriter::new(Vec::new(), header).unwrap().with_index();
    for r in records {
        w.write(r).unwrap();
    }
    w.finish_with_index().unwrap()
}

/// The BAM must match before the index means anything.
#[test]
fn the_indexed_bam_files_are_themselves_byte_identical() {
    let cases = cases();
    let goldens = goldens("bam");
    assert_eq!(cases.len(), goldens.len(), "case lists have drifted");

    let mut failures = Vec::new();
    for ((name, header, records), (golden_name, expected)) in cases.iter().zip(&goldens) {
        assert_eq!(name, golden_name, "case lists out of step");
        let (bam, _) = build(header, records);
        if &bam != expected {
            failures.push(format!("{name}: {} vs {} bytes", bam.len(), expected.len()));
        }
    }
    assert!(
        failures.is_empty(),
        "BAM divergences:\n{}",
        failures.join("\n")
    );
}

#[test]
fn every_bai_is_byte_identical_to_htsjdks() {
    let cases = cases();
    let goldens = goldens("bai");
    assert_eq!(cases.len(), goldens.len(), "case lists have drifted");

    let mut failures = Vec::new();
    for ((name, header, records), (golden_name, expected)) in cases.iter().zip(&goldens) {
        assert_eq!(name, golden_name, "case lists out of step");
        let (_, bai) = build(header, records);
        if &bai != expected {
            let at = bai
                .iter()
                .zip(expected)
                .position(|(a, b)| a != b)
                .unwrap_or(bai.len().min(expected.len()));
            failures.push(format!(
                "{name}: {} vs {} bytes, first differs at {at}",
                bai.len(),
                expected.len()
            ));
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} indices diverge from htsjdk:\n{}",
        failures.len(),
        cases.len(),
        failures.join("\n")
    );
}

/// The corpus must actually contain the structures the port is claiming to reproduce, or a
/// passing run says nothing about them.
#[test]
fn the_goldens_exercise_the_structures_that_matter() {
    let by_name: std::collections::HashMap<String, Vec<u8>> = goldens("bai").into_iter().collect();

    // The empty case is the minimum: magic, n_ref, one null reference, no-coordinate count.
    let empty = &by_name["empty"];
    assert_eq!(&empty[0..4], b"BAI\x01");
    assert_eq!(empty.len(), 4 + 4 + 8 + 8);

    // The pseudo-bin must appear, and be numbered 37450.
    let one = &by_name["one_read"];
    let n_bin = i32::from_le_bytes(one[8..12].try_into().unwrap());
    assert_eq!(n_bin, 2, "one real bin plus the pseudo-bin");
    let pseudo_at = 12 + 8 + 16;
    assert_eq!(
        i32::from_le_bytes(one[pseudo_at..pseudo_at + 4].try_into().unwrap()),
        37450
    );

    // The sparse case must have back-filled windows: several consecutive equal entries.
    let sparse = &by_name["sparse_windows"];
    let n_bin = i32::from_le_bytes(sparse[8..12].try_into().unwrap()) as usize;
    // Skip the bins to reach n_intv. Each real bin is 4 + 4 + 16*n_chunk.
    let mut p = 12;
    for _ in 0..n_bin {
        p += 4;
        let n_chunk = i32::from_le_bytes(sparse[p..p + 4].try_into().unwrap()) as usize;
        p += 4 + 16 * n_chunk;
    }
    let n_intv = i32::from_le_bytes(sparse[p..p + 4].try_into().unwrap()) as usize;
    p += 4;
    let entries: Vec<i64> = (0..n_intv)
        .map(|i| i64::from_le_bytes(sparse[p + i * 8..p + i * 8 + 8].try_into().unwrap()))
        .collect();
    assert!(n_intv > 40, "the sparse case must span many windows");
    assert!(
        entries.windows(2).any(|w| w[0] == w[1] && w[0] != 0),
        "empty windows must have inherited a previous non-zero offset"
    );
    assert!(
        !entries.contains(&-1),
        "no window may be left uninitialised in a written index"
    );

    // The many-blocks case must produce more than one chunk somewhere, or block-boundary
    // behaviour is untested.
    let many = &by_name["many_blocks"];
    let n_bin = i32::from_le_bytes(many[8..12].try_into().unwrap()) as usize;
    let mut p = 12;
    let mut max_chunks = 0;
    for _ in 0..n_bin {
        p += 4;
        let n_chunk = i32::from_le_bytes(many[p..p + 4].try_into().unwrap()) as usize;
        max_chunks = max_chunks.max(n_chunk);
        p += 4 + 16 * n_chunk;
    }
    assert!(
        max_chunks > 1,
        "chunks must stop coalescing across BGZF blocks; largest chunk list was {max_chunks}"
    );
}
