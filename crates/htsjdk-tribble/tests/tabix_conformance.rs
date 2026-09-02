//! `TabixIndexCreator`'s `.tbi` against the reference's, over feature streams the dump carries.
//!
//! The creator is fed features and their file positions directly on both sides, so what is
//! compared is the index and not a codec in front of it. Each case is dumped twice: `body` is the
//! little-endian index, and `file` is the same index inside the BGZF stream that lands beside the
//! feature file, so the composition is pinned as well as the arithmetic.
//!
//! One refusal quotes the FEATURE's own `toString`, which is `VariantContext`'s in GATK and a
//! record's in the harness. It cannot be composed here, so [`FeatureRef::description`] carries it
//! and this test supplies the harness's text; the wrapper around it is what is compared.
//!
//! While the suite is `golden-pending` the dump is named by `TABIX_INDEX_DUMP`.

use htsjdk_tribble::tabix::{FeatureRef, TabixFormat, TabixIndexCreator};

/// One case's features, as `(contig, start, end, file position)`.
type Rows<'a> = &'a [(&'a str, i32, i32, i64)];

/// The `toString` of the harness's own record, which one refusal quotes.
fn description(row: &(&str, i32, i32, i64)) -> String {
    format!(
        "Row[contig={}, start={}, end={}, position={}]",
        row.0, row.1, row.2, row.3
    )
}

/// The cases the harness runs, in its own order. The features are the inputs; the golden is the
/// answer.
#[allow(clippy::type_complexity)]
const CASES: &[(&str, TabixFormat, i64, Rows<'static>)] = &[
    (
        "one-feature",
        TabixFormat::VCF,
        4096,
        &[("chr1", 100, 100, 512)],
    ),
    (
        "two-in-one-window",
        TabixFormat::VCF,
        8192,
        &[("chr1", 100, 100, 512), ("chr1", 200, 200, 1024)],
    ),
    (
        "gap-in-the-linear-index",
        TabixFormat::VCF,
        1 << 20,
        &[("chr1", 1, 1, 512), ("chr1", 200_000, 200_000, 65536)],
    ),
    (
        "spans-many-windows",
        TabixFormat::VCF,
        1 << 20,
        &[("chr1", 1, 100_000, 512), ("chr1", 150_000, 150_000, 70000)],
    ),
    (
        "unset-end",
        TabixFormat::VCF,
        8192,
        &[("chr1", 16385, 0, 512), ("chr1", 16385, 0, 1024)],
    ),
    (
        "two-contigs",
        TabixFormat::VCF,
        8192,
        &[("chr1", 100, 100, 512), ("chr2", 100, 100, 2048)],
    ),
    (
        "every-bin-level",
        TabixFormat::VCF,
        1 << 28,
        &[
            ("chr1", 1, 1000, 512),
            ("chr1", 20_000, 100_000, 1024),
            ("chr1", 200_000, 900_000, 2048),
            ("chr1", 1_000_000, 9_000_000, 4096),
            ("chr1", 10_000_000, 200_000_000, 8192),
        ],
    ),
    (
        "format-gff",
        TabixFormat::GFF,
        4096,
        &[("chr1", 100, 200, 512)],
    ),
    (
        "format-bed",
        TabixFormat::BED,
        4096,
        &[("chr1", 100, 200, 512)],
    ),
    (
        "format-sam",
        TabixFormat::SAM,
        4096,
        &[("chr1", 100, 200, 512)],
    ),
    (
        "format-psltbl",
        TabixFormat::PSLTBL,
        4096,
        &[("chr1", 100, 200, 512)],
    ),
    (
        "long-and-wide-names",
        TabixFormat::VCF,
        8192,
        &[
            ("a_very_long_contig_name_0123456789", 100, 100, 512),
            ("chré", 100, 100, 2048),
        ],
    ),
    ("no-features", TabixFormat::VCF, 0, &[]),
    (
        "virtual-offsets",
        TabixFormat::VCF,
        30000 << 16,
        &[
            ("chr1", 100, 100, (12345 << 16) | 678),
            ("chr1", 20_000, 20_000, (23456 << 16) | 90),
        ],
    ),
    (
        "sequence-out-of-order",
        TabixFormat::VCF,
        8192,
        &[
            ("chr1", 100, 100, 512),
            ("chr2", 100, 100, 1024),
            ("chr1", 100, 100, 2048),
        ],
    ),
    (
        "features-out-of-order",
        TabixFormat::VCF,
        8192,
        &[("chr1", 200, 200, 512), ("chr1", 100, 100, 1024)],
    ),
    (
        "position-did-not-advance",
        TabixFormat::VCF,
        8192,
        &[("chr1", 100, 100, 512), ("chr1", 200, 200, 512)],
    ),
    (
        "final-position-did-not-advance",
        TabixFormat::VCF,
        512,
        &[("chr1", 100, 100, 512)],
    ),
    (
        "equal-starts",
        TabixFormat::VCF,
        8192,
        &[("chr1", 100, 500, 512), ("chr1", 100, 200, 1024)],
    ),
];

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn field(dump: &str, kind: &str, name: &str) -> Option<String> {
    let prefix = format!("{kind}\t{name}\t");
    dump.lines()
        .find(|line| line.starts_with(&prefix))
        .map(|line| line[prefix.len()..].to_string())
}

/// What this port answered for one case: the index as body and as file, or the refusal.
fn build(
    format: TabixFormat,
    final_position: i64,
    rows: Rows<'_>,
) -> Result<(Vec<u8>, Vec<u8>), String> {
    let mut creator = TabixIndexCreator::new(format);
    for row in rows {
        let text = description(row);
        let feature = FeatureRef {
            contig: row.0,
            start: row.1,
            end: row.2,
            description: &text,
            // No sequence dictionary, which is what the harness's creator was given.
            sequence_length: 0,
        };
        if let Err(error) = creator.add_feature(feature, row.3) {
            return Err(format!("{}: {}", error.java_class(), error.message()));
        }
    }
    match creator.finish(final_position) {
        Ok(index) => Ok((index.write_body(), index.write())),
        Err(error) => Err(format!("{}: {}", error.java_class(), error.message())),
    }
}

#[test]
fn every_stream_indexes_as_the_reference_indexes_it() {
    let dump = match std::env::var("TABIX_INDEX_DUMP") {
        Ok(path) => std::fs::read_to_string(path).expect("the dump named by TABIX_INDEX_DUMP"),
        Err(_) => {
            println!(
                "skipped: the tabix-index golden is still pending. Run the suite and point \
                 TABIX_INDEX_DUMP at tools/conformance/pending/tabix-index.TabixIndexDump.txt"
            );
            return;
        }
    };

    let mut compared = 0;
    for (case, format, final_position, rows) in CASES {
        match build(*format, *final_position, rows) {
            Ok((body, file)) => {
                let expected = field(&dump, "body", case)
                    .unwrap_or_else(|| panic!("{case}: the golden has no body, so it refused"));
                assert_eq!(hex(&body), expected, "{case}");
                // The same index once block compressed, which is the file that lands on disk. It
                // depends on the deflate pin like every other BGZF golden.
                let expected = field(&dump, "file", case)
                    .unwrap_or_else(|| panic!("{case}: the golden has no file"));
                assert_eq!(hex(&file), expected, "{case}: block compressed");
            }
            Err(refusal) => {
                let expected = field(&dump, "error", case)
                    .unwrap_or_else(|| panic!("{case}: the golden has no error, so it succeeded"));
                assert_eq!(refusal, expected, "{case}");
            }
        }
        compared += 1;
    }
    assert_eq!(compared, CASES.len());

    // Every row of the dump is answered: a case added to the harness and not here would otherwise
    // pass unnoticed.
    let bodies = dump.lines().filter(|l| l.starts_with("body\t")).count();
    let errors = dump.lines().filter(|l| l.starts_with("error\t")).count();
    assert_eq!(
        bodies + errors,
        CASES.len(),
        "the dump carries a case this test does not"
    );
    println!("cases={compared}");
}
