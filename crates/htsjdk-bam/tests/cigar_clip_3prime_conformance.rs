//! `CigarUtil.softClip3PrimeEndOfRead` against the reference, as a record transform.
//!
//! Each row is a cigar, a strand, a start and a clip point, answered with the new cigar, the new
//! start, whether the record ended up unmapped, and whether `NM`/`MD`/`UQ` were dropped. The last
//! two are separate rules and the dump keeps them apart: the tags go whenever the reference length
//! changed, which is far more often than the record is unmapped.
//!
//! The golden is committed once CI has produced it; until then `CIGAR_CLIP_3PRIME_DUMP` names the
//! dump (decision 0008).

use std::io::Read;
use std::path::Path;

use htsjdk_bam::cigar::{soft_clip_3prime_end_of_read, Cigar, CigarElement, Op};

/// `TextCigarCodec.decode`, for the shapes this corpus uses.
fn parse_cigar(text: &str) -> Cigar {
    let mut elements = Vec::new();
    let mut length = 0u32;
    for c in text.chars() {
        if c.is_ascii_digit() {
            length = length * 10 + c.to_digit(10).expect("a digit");
        } else {
            let op = match c {
                'M' => Op::M,
                'I' => Op::I,
                'D' => Op::D,
                'N' => Op::N,
                'S' => Op::S,
                'H' => Op::H,
                'P' => Op::P,
                '=' => Op::Eq,
                'X' => Op::X,
                other => panic!("unknown cigar operator {other}"),
            };
            elements.push(CigarElement { length, op });
            length = 0;
        }
    }
    Cigar::new(elements)
}

fn show(cigar: &Cigar) -> String {
    cigar
        .elements
        .iter()
        .map(|e| format!("{}{}", e.length, e.op.to_char() as char))
        .collect()
}

#[test]
fn every_clip_matches_the_reference() {
    // The golden was produced by the pinned container on real x86-64 and is re-derived on every
    // run; `CIGAR_CLIP_3PRIME_DUMP` still overrides it, which is how a harness change is checked before CI
    // sees it.
    let dump = match std::env::var("CIGAR_CLIP_3PRIME_DUMP") {
        Ok(path) => {
            std::fs::read_to_string(path).expect("the dump named by CIGAR_CLIP_3PRIME_DUMP")
        }
        Err(_) => {
            let golden =
                Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/cigar_clip_3prime.txt.gz");
            let file = std::fs::File::open(&golden).expect("the committed golden");
            let mut text = String::new();
            flate2::read::GzDecoder::new(file)
                .read_to_string(&mut text)
                .expect("the golden decompresses");
            text
        }
    };

    let mut rows = 0;
    for line in dump.lines() {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        let ["clip3", cigar_text, strand, start, clip_from, expected_cigar, expected_start, expected_unmapped, expected_dropped] =
            fields.as_slice()
        else {
            panic!("unrecognized dump line: {line}");
        };
        let clipped = soft_clip_3prime_end_of_read(
            &parse_cigar(cigar_text),
            *strand == "-",
            start.parse().expect("a start"),
            clip_from.parse().expect("a clip point"),
        )
        .expect("the reference did not throw on this row");

        // An unmapped result writes `*` as its cigar and zero as its start, which is
        // SAMRecord.NO_ALIGNMENT_CIGAR and NO_ALIGNMENT_START.
        let (ours_cigar, ours_start) = match &clipped.cigar {
            None => ("*".to_string(), 0),
            Some(cigar) => (show(cigar), clipped.alignment_start),
        };
        let context = format!("{cigar_text} {strand} from {clip_from}");
        assert_eq!(ours_cigar, *expected_cigar, "cigar: {context}");
        assert_eq!(ours_start.to_string(), *expected_start, "start: {context}");
        assert_eq!(
            clipped.unmapped.to_string(),
            *expected_unmapped,
            "unmapped: {context}"
        );
        assert_eq!(
            clipped.invalidate_nm_md_uq.to_string(),
            *expected_dropped,
            "tags dropped: {context}"
        );
        rows += 1;
    }
    assert_eq!(rows, 120, "ten cigars, two strands, six clip points");
    println!("rows={rows}");
}
