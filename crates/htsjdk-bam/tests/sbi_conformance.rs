//! `SBIIndexWriter`'s bytes against the reference's, as a digest, a length and the two counters.
//!
//! The record stream is rebuilt here rather than read from a fixture: record `i` is written at the
//! virtual offset `(i << 16) | (i % 7)`, which is what `SbiDump` feeds it.
//!
//! The golden is committed and re-derived by the `sbi` suite on every run; the dump can
//! still be overridden with an environment variable while a harness change is being checked.

use std::io::Read;
use std::path::Path;

use htsjdk_bam::sbi::SbiIndexWriter;
use md5::{Digest, Md5};

#[test]
fn every_index_matches_the_reference() {
    // The golden was produced by the pinned container on real x86-64 and is re-derived on every
    // run; `SBI_DUMP` still overrides it, which is how a harness change is checked before CI
    // sees it.
    let dump = match std::env::var("SBI_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by SBI_DUMP"),
        Err(_) => {
            let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/sbi.txt.gz");
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
        let ["sbi", granularity, records, final_offset, file_length, expected_len, expected_md5, expected_records, expected_offsets] =
            fields.as_slice()
        else {
            panic!("unrecognized dump line: {line}");
        };
        let granularity: u64 = granularity.parse().expect("a granularity");
        let records: u64 = records.parse().expect("a record count");

        let mut writer = SbiIndexWriter::new(granularity);
        for i in 0..records {
            writer
                .process_record((i << 16) | (i % 7))
                .expect("the offsets increase");
        }
        let bytes = writer
            .finish(
                final_offset.parse().expect("a final offset"),
                file_length.parse().expect("a file length"),
                None,
                None,
            )
            .expect("nothing here is refused");

        let mut hasher = Md5::new();
        hasher.update(&bytes);
        let digest: String = hasher
            .finalize()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();

        let context = format!("granularity {granularity}, {records} records");
        assert_eq!(bytes.len().to_string(), *expected_len, "length: {context}");
        assert_eq!(digest, *expected_md5, "digest: {context}");
        assert_eq!(
            u64::from_le_bytes(bytes[44..52].try_into().unwrap()).to_string(),
            *expected_records,
            "record count: {context}"
        );
        assert_eq!(
            u64::from_le_bytes(bytes[60..68].try_into().unwrap()).to_string(),
            *expected_offsets,
            "offset count: {context}"
        );
        rows += 1;
    }
    assert_eq!(rows, 20, "four granularities against five record counts");
    println!("rows={rows}");
}
