//! How fast this crate writes a VCF, against `noodles-vcf` on the same records.
//!
//! Not a criterion benchmark and not a claim about either library's ceiling: one process, one
//! machine, wall clock, no warm-up beyond a discarded first round. It answers one question — is
//! reproducing htsjdk's bytes costing an order of magnitude — and it is ignored by default because
//! a timing in CI is a flaky test.
//!
//! Run it with:
//!
//! ```text
//! cargo test --release -p htsjdk-vcf --test write_speed -- --ignored --nocapture
//! ```

use std::time::Instant;

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{Value, VariantContext};
use htsjdk_vcf::vcf_file::write_vcf;

const RECORDS: usize = 50_000;
const ROUNDS: usize = 5;

fn header() -> VcfHeader {
    let mut h = VcfHeader::new();
    h.lines.push(HeaderLine::info(
        "DP",
        Cardinality::Fixed(1),
        LineType::Integer,
        "Depth",
    ));
    h.lines.push(HeaderLine::info(
        "AF",
        Cardinality::A,
        LineType::Float,
        "Allele Frequency",
    ));
    h.lines.push(HeaderLine::filter("LowQual", "Low quality"));
    h.lines.push(HeaderLine::contig("chr1", 250_000_000, 0));
    h
}

fn ours_records() -> Vec<VariantContext> {
    (0..RECORDS)
        .map(|i| {
            let mut r = VariantContext::new(
                "chr1",
                (i as i64) * 10 + 1,
                vec![
                    Allele::from_str("A", true).unwrap(),
                    Allele::from_str("T", false).unwrap(),
                ],
            );
            r.stop = (i as i64) * 10 + 1;
            r.log10_p_error = -5.0;
            r.filters = Some(Vec::new());
            r.attributes = vec![
                ("AF".to_string(), Value::Str("0.5".to_string())),
                ("DP".to_string(), Value::Str(format!("{}", 10 + i % 90))),
            ];
            r
        })
        .collect()
}

#[test]
#[ignore = "a timing, not an assertion"]
fn how_fast_is_this_crate_against_noodles() {
    use noodles_vcf::header::record::value::map::info::{Number, Type};
    use noodles_vcf::header::record::value::map::{Contig, Filter, Info};
    use noodles_vcf::header::record::value::Map;
    use noodles_vcf::header::FileFormat;
    use noodles_vcf::variant::io::Write as _;
    use noodles_vcf::variant::record_buf::info::field::Value as IV;
    use noodles_vcf::variant::record_buf::{AlternateBases, Filters};
    use noodles_vcf::variant::RecordBuf;

    let h = header();
    let records = ours_records();

    let nh = noodles_vcf::Header::builder()
        .set_file_format(FileFormat::new(4, 2))
        .add_filter("LowQual", Map::<Filter>::new("Low quality"))
        .add_info(
            "AF",
            Map::<Info>::new(Number::AlternateBases, Type::Float, "Allele Frequency"),
        )
        .add_info(
            "DP",
            Map::<Info>::new(Number::Count(1), Type::Integer, "Depth"),
        )
        .add_contig("chr1", {
            let mut c = Map::<Contig>::new();
            *c.length_mut() = Some(250_000_000);
            c
        })
        .build();
    let nrecords: Vec<RecordBuf> = (0..RECORDS)
        .map(|i| {
            RecordBuf::builder()
                .set_reference_sequence_name("chr1")
                .set_variant_start(noodles_core::Position::try_from(i * 10 + 1).unwrap())
                .set_reference_bases("A")
                .set_alternate_bases(AlternateBases::from(vec![String::from("T")]))
                .set_quality_score(50.0)
                .set_filters(Filters::pass())
                .set_info(
                    [
                        ("AF".to_string(), Some(IV::Float(0.5))),
                        ("DP".to_string(), Some(IV::Integer((10 + i % 90) as i32))),
                    ]
                    .into_iter()
                    .collect(),
                )
                .build()
        })
        .collect();

    let mut ours_best = f64::MAX;
    let mut theirs_best = f64::MAX;
    let mut ours_bytes = 0;
    let mut theirs_bytes = 0;

    for round in 0..ROUNDS {
        let t = Instant::now();
        let out = write_vcf(&h, &records).expect("ours writes");
        let elapsed = t.elapsed().as_secs_f64();
        ours_bytes = out.len();
        if round > 0 {
            ours_best = ours_best.min(elapsed);
        }

        let t = Instant::now();
        let mut w = noodles_vcf::io::Writer::new(Vec::with_capacity(ours_bytes));
        w.write_header(&nh).expect("header");
        for r in &nrecords {
            w.write_variant_record(&nh, r).expect("record");
        }
        let out = w.into_inner();
        let elapsed = t.elapsed().as_secs_f64();
        theirs_bytes = out.len();
        if round > 0 {
            theirs_best = theirs_best.min(elapsed);
        }
    }

    let mib = |bytes: usize, secs: f64| (bytes as f64 / (1024.0 * 1024.0)) / secs;
    println!("records            {RECORDS}");
    println!(
        "ours     {ours_best:>8.4} s  {:>8.1} MiB/s  {ours_bytes} bytes",
        mib(ours_bytes, ours_best)
    );
    println!(
        "noodles  {theirs_best:>8.4} s  {:>8.1} MiB/s  {theirs_bytes} bytes",
        mib(theirs_bytes, theirs_best)
    );
    println!(
        "ratio    {:.2}x  (>1 means this crate is slower)",
        ours_best / theirs_best
    );
}
