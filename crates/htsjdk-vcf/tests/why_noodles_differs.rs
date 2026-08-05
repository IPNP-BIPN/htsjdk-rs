//! Why noodles does not write what GATK writes, and how much of that is forcible.
//!
//! Not a conformance suite: a measurement, asserted, so the boundary between decision 0039's two
//! surfaces is a diff rather than an argument. It exists because "just use the crate for writing
//! too" is a reasonable thing to ask and the answer is specific rather than philosophical.
//!
//! On the simplest VCF there is, the two writers differ three ways. **Two are forcible through
//! noodles' API and one is not.**
//!
//! | | forcible |
//! |---|---|
//! | noodles writes the newest version it knows, `VCFv4.5`; htsjdk always writes `VCFv4.2` | **yes**, `set_file_format` |
//! | noodles keeps INFO fields in insertion order; htsjdk sorts them | **yes**, insert them sorted |
//! | noodles writes header lines grouped by kind; htsjdk sorts across kinds | **no** |
//!
//! The third is hard-coded in `noodles_vcf::io::writer::header::write_header`, which iterates
//! `infos()`, then `filters()`, then `formats()`, then `contigs()`. Insertion order decides the
//! order *within* a kind and nothing decides the order *between* kinds. htsjdk sorts every line by
//! its rendered string, so `FILTER=<...>` precedes `INFO=<...>` precedes `contig=<...>`, and no
//! sequence of `add_*` calls reaches that.
//!
//! And forcing it would not be enough even with a fork. Decision 0016 measured htsjdk's comparator
//! and found it **is not a total order**: it has a cycle, so on some headers the output depends on
//! the insertion order of a Java `TreeSet`. A correct implementation cannot reproduce that. Only a
//! port of the same broken comparator can, which is what `VcfHeader::write` is.
//!
//! None of this is a defect in noodles. The specification does not fix any of the three, so two
//! correct implementations differ — and this is the *text* format, the mildest case. A BAM adds a
//! deflate stream on top, where decision 0001 already measured two zlibs disagreeing.

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{Value, VariantContext};
use htsjdk_vcf::vcf_file::write_vcf;

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
    h.lines.push(HeaderLine::contig("chr1", 100_000, 0));
    h
}

fn record() -> VariantContext {
    let mut r = VariantContext::new(
        "chr1",
        100,
        vec![
            Allele::from_str("A", true).unwrap(),
            Allele::from_str("T", false).unwrap(),
        ],
    );
    r.stop = 100;
    r.id = "rs42".into();
    r.log10_p_error = -5.0;
    r.filters = Some(Vec::new());
    r.attributes = vec![
        ("DP".to_string(), Value::Str("10".to_string())),
        ("AF".to_string(), Value::Str("0.5".to_string())),
    ];
    r
}

#[test]
fn print_both_writers() {
    let ours = write_vcf(&header(), &[record()]).expect("ours writes");

    use noodles_vcf::header::record::value::map::info::{Number, Type};
    use noodles_vcf::header::record::value::map::{Contig, Filter, Info};
    use noodles_vcf::header::record::value::Map;
    use noodles_vcf::variant::record_buf::info::field::Value as IV;
    use noodles_vcf::variant::record_buf::{AlternateBases, Filters};
    use noodles_vcf::variant::RecordBuf;

    // FORCING ATTEMPT: insert in htsjdk's sorted order, and pin the version to 4.2.
    use noodles_vcf::header::FileFormat;
    let mut nh = noodles_vcf::Header::builder()
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
            *c.length_mut() = Some(100_000);
            c
        })
        .build();
    let _ = &mut nh;

    let nr = RecordBuf::builder()
        .set_reference_sequence_name("chr1")
        .set_variant_start(noodles_core::Position::try_from(100).unwrap())
        .set_ids([String::from("rs42")].into_iter().collect())
        .set_reference_bases("A")
        .set_alternate_bases(AlternateBases::from(vec![String::from("T")]))
        .set_quality_score(50.0)
        .set_filters(Filters::pass())
        .set_info(
            [
                ("AF".to_string(), Some(IV::Float(0.5))),
                ("DP".to_string(), Some(IV::Integer(10))),
            ]
            .into_iter()
            .collect(),
        )
        .build();

    let mut w = noodles_vcf::io::Writer::new(Vec::new());
    w.write_header(&nh).expect("header");
    use noodles_vcf::variant::io::Write as _;
    w.write_variant_record(&nh, &nr).expect("record");
    let theirs = String::from_utf8(w.into_inner()).expect("utf8");

    println!("=== OURS (htsjdk bytes, oracle-backed) ===\n{ours}");
    println!("=== NOODLES, forced as far as the API allows ===\n{theirs}");

    let differing: Vec<usize> = ours
        .lines()
        .zip(theirs.lines())
        .enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(i, (a, b))| {
            println!("line {i}\n  ours:    {a}\n  noodles: {b}");
            i
        })
        .collect();

    // The two forcible differences are gone: the version line and the INFO field order both match.
    assert_eq!(
        ours.lines().next(),
        theirs.lines().next(),
        "the version line is forcible with set_file_format"
    );
    assert_eq!(
        ours.lines().last(),
        theirs.lines().last(),
        "the INFO field order is forcible by inserting them sorted"
    );

    // The one that is not: header lines 1 to 3, where htsjdk sorts across kinds and noodles groups
    // by kind. If this ever becomes empty, noodles changed and this file should say so.
    assert_eq!(
        differing,
        vec![1, 2, 3],
        "exactly the header lines htsjdk sorts and noodles groups"
    );
    assert_ne!(ours, theirs, "and so the files are not identical");
}
