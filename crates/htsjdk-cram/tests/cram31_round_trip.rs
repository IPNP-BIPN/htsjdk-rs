//! CRAM 3.1 written by [`htsjdk_cram::writer`] and read back by [`htsjdk_cram::reader`].
//!
//! The byte-for-byte claim against htsjdk 5.0.0 is measured downstream, by picard-rs, whose oracle
//! is Picard 3.5.0: this crate's oracle is htsjdk 4.2.0, which writes CRAM 3.0. What is checked
//! here is what a round trip must keep, and what it must change the way htsjdk does.

use htsjdk_bam::cigar::Cigar;
use htsjdk_bam::header::{ReadGroup, SamHeader, SequenceRecord};
use htsjdk_bam::record::BamRecord;
use htsjdk_bam::tag::{Tag, TagValue};
use htsjdk_cram::reader::{CramReader, ReferenceSource};
use htsjdk_cram::writer::{write_cram, write_cram_indexed};
use std::collections::HashMap;
use std::sync::Arc;

struct Memory(HashMap<String, Arc<Vec<u8>>>);

impl ReferenceSource for Memory {
    fn sequence_bases(&mut self, name: &str) -> Option<Arc<Vec<u8>>> {
        self.0.get(name).cloned()
    }
}

/// A deterministic contig: a linear congruential walk over ACGT.
fn contig(length: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..length)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 33) as usize % 4]
        })
        .collect()
}

fn reference() -> (SamHeader, Memory) {
    let mut header = SamHeader::new();
    header.set_sort_order("coordinate");
    let mut map = HashMap::new();
    for (i, (name, length)) in [("chr1", 30_000usize), ("chr2", 8_000)].iter().enumerate() {
        header
            .sequences
            .push(SequenceRecord::new(name, *length as i32));
        map.insert(name.to_string(), Arc::new(contig(*length, i as u64 + 7)));
    }
    for id in ["rg1", "rg2"] {
        let mut group = ReadGroup::new(id);
        group.attributes.set("SM", "sample");
        header.read_groups.push(group);
    }
    (header, Memory(map))
}

fn record(
    name: &str,
    reference: i32,
    start: i32,
    cigar: &str,
    bases: Vec<u8>,
    flags: u16,
) -> BamRecord {
    let mut r = BamRecord {
        read_name: name.to_string(),
        flags,
        reference_index: reference,
        alignment_start: start,
        mapping_quality: 37,
        cigar: htsjdk_bam::text_parse::parse_cigar(cigar).unwrap(),
        base_qualities: (0..bases.len()).map(|i| (i % 41) as u8).collect(),
        read_bases: bases,
        ..BamRecord::default()
    };
    r.tags
        .insert(Tag::new(b"RG"), TagValue::Str("rg1".to_string()));
    r
}

/// A mixture of everything the writer turns into read features, with tags of every type.
fn records(header: &SamHeader, source: &mut Memory, count: usize) -> Vec<BamRecord> {
    let chr1 = source.sequence_bases("chr1").unwrap();
    let mut out = Vec::new();
    for i in 0..count {
        let start = 1 + (i * 2) as i32 % 29_000;
        let mut bases = chr1[start as usize - 1..start as usize - 1 + 60].to_vec();
        let cigar = match i % 5 {
            0 => "60M",
            1 => {
                bases[10] = if bases[10] == b'A' { b'C' } else { b'A' };
                bases[20] = b'N';
                "60M"
            }
            2 => "5S50M5S",
            3 => {
                bases.truncate(57);
                "20M3D20M3I14M"
            }
            _ => "30M100N30M",
        };
        let mut r = record(&format!("read{i:06}"), 0, start, cigar, bases, 0);
        if i % 7 == 0 {
            r.tags.insert(Tag::new(b"XC"), TagValue::Char(b'q'));
        }
        if i % 3 == 0 {
            r.tags
                .insert(Tag::new(b"XI"), TagValue::Int(i as i64 * 997 - 40_000));
            r.tags
                .insert(Tag::new(b"XF"), TagValue::Float(i as f32 / 3.0));
            r.tags
                .insert(Tag::new(b"XZ"), TagValue::Str(format!("value{i}")));
        }
        if i % 11 == 0 {
            r.tags.insert(
                Tag::new(b"XB"),
                TagValue::ByteArray {
                    values: vec![1, -2, 3],
                    unsigned: false,
                },
            );
        }
        if i % 13 == 0 {
            // A wrong NM is kept verbatim.
            r.tags.insert(Tag::new(b"NM"), TagValue::Int(99));
        }
        out.push(r);
    }
    // Pairs on chr2, then unplaced unmapped reads last.
    for i in 0..40 {
        let start = 100 + i * 50;
        let chr2 = source.sequence_bases("chr2").unwrap();
        let mate = start + 120;
        let a = chr2[start as usize - 1..start as usize + 49].to_vec();
        let b = chr2[mate as usize - 1..mate as usize + 49].to_vec();
        let mut first = record(
            &format!("pair{i}"),
            1,
            start,
            "50M",
            a,
            0x1 | 0x2 | 0x20 | 0x40,
        );
        let mut second = record(
            &format!("pair{i}"),
            1,
            mate,
            "50M",
            b,
            0x1 | 0x2 | 0x10 | 0x80,
        );
        first.mate_reference_index = 1;
        first.mate_alignment_start = mate;
        second.mate_reference_index = 1;
        second.mate_alignment_start = start;
        // htslib's template length, which makes the writer attach the pair.
        first.inferred_insert_size = 170;
        second.inferred_insert_size = -170;
        out.push(first);
        out.push(second);
    }
    out.sort_by_key(|r| (r.reference_index, r.alignment_start));
    for i in 0..10 {
        let mut r = record(&format!("unmapped{i}"), -1, 0, "*", contig(40, i), 0x4);
        r.cigar = Cigar::default();
        r.mapping_quality = 0;
        r.tags
            .insert(Tag::new(b"RG"), TagValue::Str("rg2".to_string()));
        out.push(r);
    }
    let _ = header;
    out
}

fn read(cram: &[u8], source: &mut Memory) -> (SamHeader, Vec<BamRecord>) {
    let mut reader = CramReader::new(cram).unwrap();
    let records = reader.read_all(source).unwrap();
    (reader.header.clone(), records)
}

#[test]
fn records_survive_the_round_trip_and_md_nm_are_regenerated() {
    let (header, mut source) = reference();
    let written = records(&header, &mut source, 2_500);
    let cram = write_cram(&header, &written, "file:///test.cram", &mut source).unwrap();
    assert_eq!(&cram[..6], b"CRAM\x03\x01");
    let (read_header, read_back) = read(&cram, &mut source);
    assert_eq!(read_header.sequences, header.sequences);
    assert_eq!(read_back.len(), written.len());

    for (before, after) in written.iter().zip(&read_back) {
        assert_eq!(before.read_name, after.read_name);
        assert_eq!(before.flags, after.flags, "{}", before.read_name);
        assert_eq!(before.reference_index, after.reference_index);
        assert_eq!(before.alignment_start, after.alignment_start);
        assert_eq!(before.cigar, after.cigar, "{}", before.read_name);
        assert_eq!(before.read_bases, after.read_bases, "{}", before.read_name);
        assert_eq!(before.base_qualities, after.base_qualities);
        assert_eq!(before.mate_alignment_start, after.mate_alignment_start);
        assert_eq!(before.inferred_insert_size, after.inferred_insert_size);
        for (tag, value) in before.tags.iter() {
            assert_eq!(
                after.tags.get(*tag),
                Some(value),
                "{} {tag}",
                before.read_name
            );
        }
        // A mapped read without NM or MD gets both back.
        if before.flags & 0x4 == 0 {
            assert!(
                after.tags.get(Tag::new(b"MD")).is_some(),
                "{}",
                before.read_name
            );
            assert!(
                after.tags.get(Tag::new(b"NM")).is_some(),
                "{}",
                before.read_name
            );
        }
    }
}

#[test]
fn a_regenerated_md_and_nm_are_dropped_on_the_way_back_out() {
    let (header, mut source) = reference();
    let written = records(&header, &mut source, 300);
    let first = write_cram(&header, &written, "file:///a.cram", &mut source).unwrap();
    let (_, decoded) = read(&first, &mut source);
    let second = write_cram(&header, &decoded, "file:///a.cram", &mut source).unwrap();
    assert_eq!(first, second, "writing what was read gives the same file");
}

#[test]
fn many_containers_and_an_index() {
    let (header, mut source) = reference();
    let written = records(&header, &mut source, 25_000);
    let (cram, bai) =
        write_cram_indexed(&header, &written, "file:///big.cram", &mut source).unwrap();
    assert_eq!(&bai[..4], b"BAI\x01");
    let (_, read_back) = read(&cram, &mut source);
    assert_eq!(read_back.len(), written.len());
    assert!(read_back
        .iter()
        .zip(&written)
        .all(|(a, b)| a.read_name == b.read_name));
}

#[test]
fn a_missing_reference_fails_as_htsjdk_does() {
    let (header, mut source) = reference();
    let written = records(&header, &mut source, 10);
    let cram = write_cram(&header, &written, "file:///x.cram", &mut source).unwrap();
    let mut reader = CramReader::new(&cram).unwrap();
    let error = reader
        .read_all(&mut htsjdk_cram::reader::NoReference)
        .unwrap_err();
    assert_eq!(
        error.java_line(),
        "java.lang.IllegalArgumentException: Failure getting reference bases for sequence chr1"
    );
    // A writer made with no reference has the lazy source, which throws its own message.
    let error = write_cram(
        &header,
        &written,
        "file:///x.cram",
        &mut htsjdk_cram::reader::LazyReference,
    )
    .unwrap_err();
    assert_eq!(
        error.message(),
        "A reference must be supplied that includes the reference sequence for chr1."
    );
    // One given a reference that lacks the sequence refuses with the record.
    let error = write_cram(
        &header,
        &written,
        "file:///x.cram",
        &mut htsjdk_cram::reader::NoReference,
    )
    .unwrap_err();
    assert!(error.message().starts_with(
        "A reference must be supplied (reference sequence SAMSequenceRecord(name=chr1,length=30000,dict_index=0,"
    ));
}

/// A reader that hands out a few bytes per call, so the stream reader's buffer has to grow.
struct Trickle<'a>(&'a [u8]);

impl std::io::Read for Trickle<'_> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        let n = out.len().min(self.0.len()).min(7);
        out[..n].copy_from_slice(&self.0[..n]);
        self.0 = &self.0[n..];
        Ok(n)
    }
}

#[test]
fn a_streamed_cram_reads_as_the_bytes_do() {
    let (header, mut source) = reference();
    let written = records(&header, &mut source, 25_000);
    let cram = write_cram(&header, &written, "file:///s.cram", &mut source).unwrap();
    let (whole_header, whole) = read(&cram, &mut source);
    for input in [
        Box::new(&cram[..]) as Box<dyn std::io::Read>,
        Box::new(Trickle(&cram)),
    ] {
        let mut reader = CramReader::from_read(input).unwrap();
        assert_eq!(reader.header, whole_header);
        assert_eq!(reader.read_all(&mut source).unwrap(), whole);
    }
    // Cut short anywhere, the stream refuses where the bytes do, with the same records before it.
    for cut in [10, 30, cram.len() / 3, cram.len() / 2, cram.len() - 40] {
        let short = &cram[..cut];
        let whole = CramReader::new(short).map(|mut r| {
            let mut out = Vec::new();
            let end = loop {
                match r.next_container(&mut source) {
                    Ok(Some(batch)) => out.extend(batch),
                    other => break other.map(|_| ()),
                }
            };
            (out, end)
        });
        let streamed = CramReader::from_read(Trickle(short)).map(|mut r| {
            let mut out = Vec::new();
            let end = loop {
                match r.next_container(&mut source) {
                    Ok(Some(batch)) => out.extend(batch),
                    other => break other.map(|_| ()),
                }
            };
            (out, end)
        });
        assert_eq!(streamed.is_ok(), whole.is_ok(), "cut at {cut}");
        match (streamed, whole) {
            (Ok(a), Ok(b)) => assert!(a == b, "cut at {cut}"),
            (Err(a), Err(b)) => assert_eq!(a, b, "cut at {cut}"),
            _ => unreachable!(),
        }
    }
}

/// What a reader gives container by container: each batch, then how it stopped (the end, a
/// refusal or a panic), then what the next call gives after a refusal.
fn container_by_container<'a>(
    reader: Result<CramReader<'a>, htsjdk_cram::reader::CramError>,
    source: &mut Memory,
) -> Vec<String> {
    let mut reader = match reader {
        Ok(reader) => reader,
        Err(error) => return vec![format!("open: {error:?}")],
    };
    let mut seen = Vec::new();
    loop {
        let next = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            reader.next_container(source)
        }));
        match next {
            Ok(Ok(Some(batch))) => seen.push(format!("{batch:?}")),
            Ok(Ok(None)) => {
                seen.push("end".to_string());
                return seen;
            }
            Ok(Err(error)) => {
                seen.push(format!("refused: {error:?}"));
                let again = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    reader
                        .next_container(source)
                        .map(|batch| batch.map(|b| b.len()))
                }));
                seen.push(format!("then: {:?}", again.map_err(|_| "panic")));
                return seen;
            }
            Err(_) => {
                seen.push("panic".to_string());
                return seen;
            }
        }
    }
}

/// Containers decoded on many threads answer as on one: the same records container by
/// container, then the same end, refusal or panic, and the same answer to the call after it, in
/// memory and streamed, whole, cut short and with one bit changed.
#[test]
fn a_cram_read_on_many_threads_reads_as_on_one() {
    let (header, mut source) = reference();
    let written = records(&header, &mut source, 25_000);
    let cram = write_cram(&header, &written, "file:///p.cram", &mut source).unwrap();
    let pools: Vec<rayon::ThreadPool> = [1, 8]
        .into_iter()
        .map(|n| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
        })
        .collect();
    let hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let mut seed = 3u64;
    let mut next = move || {
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        (seed >> 33) as usize
    };
    let mut variants = vec![cram.clone()];
    for _ in 0..12 {
        variants.push(cram[..next() % cram.len()].to_vec());
    }
    for _ in 0..40 {
        let mut damaged = cram.clone();
        let at = next() % damaged.len();
        damaged[at] ^= 1 << (next() % 8);
        variants.push(damaged);
    }
    for (i, bytes) in variants.iter().enumerate() {
        let answers: Vec<[Vec<String>; 2]> = pools
            .iter()
            .map(|pool| {
                pool.install(|| {
                    [
                        container_by_container(CramReader::new(bytes), &mut source),
                        container_by_container(CramReader::from_read(Trickle(bytes)), &mut source),
                    ]
                })
            })
            .collect();
        assert!(answers[0] == answers[1], "variant {i}");
    }
    std::panic::set_hook(hook);
    // The whole file has three data containers, so eight threads read them together.
    assert!(written.len() > 2 * 10_000);
}
