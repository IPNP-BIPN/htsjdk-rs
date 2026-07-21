//! Conformance for VCF data lines against htsjdk's `VCFEncoder`.
//!
//! Goldens from `tools/vcf-conformance/VcfRecordDump.java` in the pinned oracle.
//!
//! The case list below mirrors the harness case for case, in the same order, and the case
//! *names* are asserted as well as the values. Without that, a port could drift a case out of
//! the list and the suite would still pass on the ones that remain.

use std::io::Read;

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::encoder::VcfEncoder;
use htsjdk_vcf::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use htsjdk_vcf::variant::{format_vcf_double, Genotype, Value, VariantContext};

fn corpus() -> Vec<(String, String)> {
    let p = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/vcf_record.txt.gz");
    let f = std::fs::File::open(&p).expect("corpus");
    let mut s = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut s)
        .expect("corpus is gzip");
    s.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let (name, value) = l.split_once('\t').expect("name TAB value");
            (name.to_string(), unescape(value))
        })
        .collect()
}

fn unescape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('n') => out.push('\n'),
            Some('t') => out.push('\t'),
            Some('\\') => out.push('\\'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// Java's `Double.toHexString`, read back. Used only to recover the exact double a golden case
/// was generated from, so that the sweep cases need not be duplicated by hand in two languages
/// and drift apart.
fn parse_java_hex_double(s: &str) -> f64 {
    let (sign, rest) = match s.strip_prefix('-') {
        Some(r) => (-1.0, r),
        None => (1.0, s),
    };
    let rest = rest
        .strip_prefix("0x")
        .expect("Double.toHexString starts 0x");
    let (mantissa, exp) = rest.split_once('p').expect("Double.toHexString has p");
    let (int_part, frac) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let mut value = int_part.parse::<f64>().expect("hex integer part is 0 or 1");
    for (i, d) in frac.chars().enumerate() {
        let digit = d.to_digit(16).expect("hex digit") as f64;
        value += digit * 16f64.powi(-(i as i32 + 1));
    }
    sign * value * 2f64.powi(exp.parse().expect("binary exponent"))
}

fn header() -> VcfHeader {
    let mut h = VcfHeader::new();
    for k in ["AA", "DP", "ZZ", "MM", "AF", "STR", "LIST", "NEG"] {
        h.lines.push(HeaderLine::info(
            k,
            Cardinality::Fixed(1),
            LineType::String,
            k,
        ));
    }
    h.lines.push(HeaderLine::info(
        "FLAG",
        Cardinality::Fixed(0),
        LineType::Flag,
        "a flag",
    ));
    for (id, n, t, d) in [
        ("GT", Cardinality::Fixed(1), LineType::String, "Genotype"),
        (
            "GQ",
            Cardinality::Fixed(1),
            LineType::Integer,
            "Genotype Quality",
        ),
        ("DP", Cardinality::Fixed(1), LineType::Integer, "Depth"),
        ("AD", Cardinality::R, LineType::Integer, "Allele Depths"),
        ("PL", Cardinality::G, LineType::Integer, "Likelihoods"),
        (
            "FT",
            Cardinality::Fixed(1),
            LineType::String,
            "Genotype Filter",
        ),
        ("XX", Cardinality::Fixed(1), LineType::String, "extended"),
    ] {
        h.lines.push(HeaderLine::format(id, n, t, d));
    }
    for (id, d) in [
        ("q10", "Quality below 10"),
        ("aFilter", "a"),
        ("zFilter", "z"),
    ] {
        h.lines.push(HeaderLine::filter(id, d));
    }
    h.samples = vec!["SAMPLE1".to_string(), "SAMPLE2".to_string()];
    h
}

fn allele(s: &str, is_ref: bool) -> Allele {
    Allele::from_str(s, is_ref).unwrap()
}

/// The harness's `base(...)`: first string is the reference, the rest are alternates.
fn base(alleles: &[&str]) -> VariantContext {
    let mut list = vec![allele(alleles[0], true)];
    for a in &alleles[1..] {
        list.push(allele(a, false));
    }
    VariantContext::new("chr1", 100, list)
}

fn gt(sample: &str, alleles: &[Allele]) -> Genotype {
    Genotype::new(sample, alleles.to_vec())
}

/// Builds every case in the harness's order. `None` marks a case produced by a direct call to
/// `formatVCFDouble` rather than by encoding a record.
fn cases() -> Vec<(String, String)> {
    let h = header();
    let enc = VcfEncoder::new(&h);
    let mut out: Vec<(String, String)> = Vec::new();
    let mut push = |name: &str, vc: &VariantContext| {
        out.push((name.to_string(), enc.encode(vc).expect(name)));
    };

    push("minimal", &base(&["A", "T"]));
    let mut with_id = base(&["A", "T"]);
    with_id.id = "rs123".to_string();
    push("with_id", &with_id);
    push("no_alt", &base(&["A"]));
    push("multi_alt", &base(&["A", "T", "C", "GG"]));
    push("span_del", &base(&["A", "T", "*"]));
    push("symbolic", &base(&["A", "<DEL>"]));
    push("breakend", &base(&["A", "A]chr2:456]"]));
    push("lowercase_bases", &base(&["acgt", "a"]));
    push("symbolic_lowercase", &base(&["A", "<del>"]));
    push("ends_with_gt", &base(&["A", "at>"]));
    push("breakend_lowercase", &base(&["A", "a]chr2:456]"]));
    push("single_breakend", &base(&["A", ".a"]));

    for q in [
        0.0, 1.0, 10.0, 29.5, 30.0, 30.004, 30.005, 30.015, 0.125, 2.675, 1234.5678, 1e-3,
        99999.999,
    ] {
        let mut vc = base(&["A", "T"]);
        vc.log10_p_error = q / -10.0;
        // Java names the case after `Double.toString(q)`, which prints an integral double with
        // a trailing `.0` and small values in exponent form.
        push(&format!("qual_{}", java_double_to_string(q)), &vc);
    }

    push("filter_none", &base(&["A", "T"]));
    let mut passed = base(&["A", "T"]);
    passed.filters = Some(Vec::new());
    push("filter_pass", &passed);
    let mut one = base(&["A", "T"]);
    one.filters = Some(vec!["q10".to_string()]);
    push("filter_one", &one);
    let mut many = base(&["A", "T"]);
    many.filters = Some(
        ["zFilter", "q10", "aFilter"]
            .iter()
            .map(|s| s.to_string())
            .collect(),
    );
    push("filter_sorted", &many);

    let mut sorted = base(&["A", "T"]);
    sorted.attributes = vec![
        ("ZZ".to_string(), Value::Str("last".to_string())),
        ("AA".to_string(), Value::Str("first".to_string())),
        ("MM".to_string(), Value::Str("middle".to_string())),
    ];
    push("info_sorted", &sorted);

    let mut with_attr = |name: &str, key: &str, v: Value| {
        let mut vc = base(&["A", "T"]);
        vc.attributes = vec![(key.to_string(), v)];
        out.push((name.to_string(), enc.encode(&vc).expect(name)));
    };
    with_attr("info_flag_true", "FLAG", Value::Bool(true));
    with_attr("info_flag_false", "FLAG", Value::Bool(false));
    with_attr("info_empty_string", "STR", Value::Str(String::new()));
    with_attr("info_null", "STR", Value::Missing);
    with_attr("info_int", "DP", Value::Int(42));
    with_attr(
        "info_list",
        "LIST",
        Value::List(vec![Value::Int(1), Value::Int(2), Value::Int(3)]),
    );
    with_attr("info_empty_list", "LIST", Value::List(vec![]));
    with_attr(
        "info_int_array",
        "LIST",
        Value::List(vec![Value::Int(4), Value::Int(5)]),
    );
    with_attr("info_double", "AF", Value::Double(0.5));
    with_attr("info_double_small", "AF", Value::Double(0.001));
    with_attr("info_double_negative", "NEG", Value::Double(-0.5));
    with_attr(
        "info_double_list",
        "LIST",
        Value::List(vec![
            Value::Double(0.5),
            Value::Double(0.001),
            Value::Double(1e-30),
        ]),
    );

    let a = allele("A", true);
    let t = allele("T", false);
    let nc = Allele::no_call();

    let mut with_gts = |name: &str, gs: Vec<Genotype>| {
        let mut vc = base(&["A", "T"]);
        vc.genotypes = gs;
        out.push((name.to_string(), enc.encode(&vc).expect(name)));
    };

    with_gts(
        "gt_simple",
        vec![
            gt("SAMPLE1", &[a.clone(), t.clone()]),
            gt("SAMPLE2", &[a.clone(), a.clone()]),
        ],
    );
    let mut phased = gt("SAMPLE1", &[a.clone(), t.clone()]);
    phased.phased = true;
    with_gts(
        "gt_phased",
        vec![phased, gt("SAMPLE2", &[a.clone(), a.clone()])],
    );
    with_gts(
        "gt_nocall",
        vec![
            gt("SAMPLE1", &[nc.clone(), nc.clone()]),
            gt("SAMPLE2", &[a.clone(), t.clone()]),
        ],
    );
    with_gts(
        "gt_haploid",
        vec![
            gt("SAMPLE1", std::slice::from_ref(&t)),
            gt("SAMPLE2", &[a.clone(), t.clone()]),
        ],
    );
    with_gts(
        "gt_triploid",
        vec![
            gt("SAMPLE1", &[a.clone(), t.clone(), t.clone()]),
            gt("SAMPLE2", &[a.clone(), t.clone()]),
        ],
    );
    with_gts(
        "gt_one_sample_only",
        vec![gt("SAMPLE1", &[a.clone(), t.clone(), t.clone()])],
    );

    let mut g1 = gt("SAMPLE1", &[a.clone(), t.clone()]);
    g1.pl = Some(vec![0, 10, 100]);
    g1.ad = Some(vec![5, 6]);
    g1.dp = Some(11);
    g1.gq = Some(30);
    let mut g2 = gt("SAMPLE2", &[a.clone(), a.clone()]);
    g2.pl = Some(vec![0, 20, 200]);
    g2.ad = Some(vec![9, 0]);
    g2.dp = Some(9);
    g2.gq = Some(40);
    with_gts("gt_all_int_fields", vec![g1, g2]);

    let mut ragged = gt("SAMPLE1", &[a.clone(), t.clone()]);
    ragged.gq = Some(30);
    ragged.dp = Some(11);
    with_gts(
        "gt_ragged",
        vec![ragged, gt("SAMPLE2", &[a.clone(), a.clone()])],
    );

    let mut only_gq = gt("SAMPLE1", &[a.clone(), t.clone()]);
    only_gq.gq = Some(30);
    let mut only_dp = gt("SAMPLE2", &[a.clone(), a.clone()]);
    only_dp.dp = Some(9);
    with_gts("gt_trailing_stripped", vec![only_gq, only_dp]);

    let mut filtered = gt("SAMPLE1", &[a.clone(), t.clone()]);
    filtered.filters = Some("q10".to_string());
    with_gts(
        "gt_filtered",
        vec![filtered, gt("SAMPLE2", &[a.clone(), a.clone()])],
    );

    let mut extended = gt("SAMPLE1", &[a.clone(), t.clone()]);
    extended.extended = vec![("XX".to_string(), Value::Str("hello".to_string()))];
    with_gts(
        "gt_extended",
        vec![extended, gt("SAMPLE2", &[a.clone(), a.clone()])],
    );

    let mut extended_d = gt("SAMPLE1", &[a.clone(), t.clone()]);
    extended_d.extended = vec![("XX".to_string(), Value::Double(0.001))];
    with_gts(
        "gt_extended_double",
        vec![extended_d, gt("SAMPLE2", &[a.clone(), a.clone()])],
    );

    out
}

/// `Double.toString` for the handful of literals the harness names cases with. Not a general
/// implementation: it covers integral values and the one exponent-form value in the list, which
/// is all the case names need, and it says so rather than pretending to be more.
fn java_double_to_string(d: f64) -> String {
    if d == 0.001 {
        return "0.001".to_string();
    }
    if d == d.trunc() && d.abs() < 1e7 {
        return format!("{d:.1}");
    }
    format!("{d}")
}

#[test]
fn every_encoded_record_matches_the_oracle() {
    let golden = corpus();
    let record_cases = cases();

    // The record cases come first in the harness, then the formatVCFDouble sweep.
    assert!(
        golden.len() > record_cases.len(),
        "the golden must also carry the formatVCFDouble sweep"
    );
    let mut mismatches = Vec::new();
    for (i, (name, ours)) in record_cases.iter().enumerate() {
        let (gname, gvalue) = &golden[i];
        assert_eq!(
            gname, name,
            "case {i} is {gname} in the golden and {name} here: the lists have drifted"
        );
        if gvalue != ours {
            mismatches.push(format!("{name}\n  htsjdk: {gvalue}\n  ours  : {ours}"));
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} of {} record cases diverge:\n{}",
        mismatches.len(),
        record_cases.len(),
        mismatches.join("\n")
    );
}

#[test]
fn every_format_vcf_double_case_matches_the_oracle() {
    let golden = corpus();
    let sweep: Vec<_> = golden
        .iter()
        .filter_map(|(n, v)| n.strip_prefix("formatVCFDouble_").map(|h| (h, v)))
        .collect();
    assert!(!sweep.is_empty(), "the golden carries the sweep");

    let mut mismatches = Vec::new();
    for (hex, expected) in &sweep {
        let d = parse_java_hex_double(hex);
        let ours = format_vcf_double(d);
        if &&ours != expected {
            mismatches.push(format!(
                "{hex} ({d:e})\n  htsjdk: {expected}\n  ours  : {ours}"
            ));
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} of {} sweep cases diverge:\n{}",
        mismatches.len(),
        sweep.len(),
        mismatches.join("\n")
    );
}

/// The `%f` / `%e` agreement rate against the oracle, over the sweep in
/// `tools/vcf-conformance/JavaFormatSweep.java`.
///
/// Ignored by default because the sweep file is generated in the pinned container rather than
/// committed: 127,803 lines of formatted doubles is corpus, not source. CI runs it with
/// `JFORMAT_SWEEP` pointing at the freshly generated file.
///
/// The thresholds are the rates measured in decision 0017, minus nothing. They are assertions
/// and not targets: if the rate *falls*, something regressed; if it rises, htsjdk or the JDK
/// changed and the decision needs revisiting either way.
#[test]
#[ignore = "needs JFORMAT_SWEEP from the pinned oracle"]
fn the_jvm_format_model_holds_at_the_measured_rate() {
    let path = std::env::var("JFORMAT_SWEEP").expect("JFORMAT_SWEEP");
    let text = std::fs::read_to_string(&path).expect("sweep file");

    let (mut n, mut bad_f, mut bad_e, mut bad_small, mut n_small) = (0u64, 0u64, 0u64, 0u64, 0u64);
    let mut smallest_divergence = f64::INFINITY;
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 4 {
            continue;
        }
        let d = f64::from_bits(u64::from_str_radix(p[0], 16).expect("hex bits"));
        n += 1;
        let small = d.abs() < 1e15;
        if small {
            n_small += 1;
        }
        let mut diverged = false;
        if htsjdk_vcf::jformat::format_fixed(d, 2) != p[1] {
            bad_f += 1;
            diverged = true;
        }
        if htsjdk_vcf::jformat::format_fixed(d, 3) != p[2] {
            diverged = true;
        }
        if htsjdk_vcf::jformat::format_scientific(d, 3) != p[3] {
            bad_e += 1;
            diverged = true;
        }
        if diverged {
            smallest_divergence = smallest_divergence.min(d.abs());
            if small {
                bad_small += 1;
            }
        }
    }
    assert!(
        n > 100_000,
        "the sweep should cover the whole corpus, got {n}"
    );

    let rate = 100.0 * (n - bad_f) as f64 / n as f64;
    let rate_small = 100.0 * (n_small - bad_small) as f64 / n_small as f64;
    println!("n={n}  %.2f agreement {rate:.4}%  below 1e15 {rate_small:.6}%");
    println!("smallest diverging magnitude {smallest_divergence:e}");

    assert_eq!(
        bad_e, 0,
        "%.3e was exact in decision 0017 and must stay exact"
    );
    assert!(rate >= 99.85, "%.2f agreement fell to {rate:.4}%");
    assert!(
        rate_small >= 99.99,
        "agreement below 1e15 fell to {rate_small:.6}%"
    );
    assert!(
        smallest_divergence > 1e14,
        "a divergence appeared at {smallest_divergence:e}, below the 6.9e14 bound in decision 0017"
    );
}
