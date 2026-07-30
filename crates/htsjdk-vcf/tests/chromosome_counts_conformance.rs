//! Conformance for `AC`, `AF` and `AN` against the oracle.
//!
//! Golden from `tools/vcf-conformance/ChromosomeCountsDump.java`, which calls
//! `VariantContextUtils.calculateChromosomeCounts` on the same twenty fixtures this file rebuilds.
//!
//! Four rows carry the whole point of the suite:
//!
//! ```text
//! type   two-alts                  java.util.ArrayList
//! type   one-alt-het               java.lang.Integer
//! af     pedigree-one-founder      4602678819172646912   (0.5, while AC says 5)
//! af     pedigree-absent-founder   -2251799813685248     (NaN, from 0/0)
//! ```
//!
//! The first two are the type change: one alternate allele gives a scalar and two or more give a
//! list, and the two render identically. The third is `AF` not being `AC / AN`: the whole cohort
//! carries five alternate chromosomes out of six, but the single declared founder is a het, so `AF`
//! is 1/2 while `AC` alongside it is 5. The fourth is a division by zero that is not an error.
//!
//! # `AF` is compared as raw bits, except where it is a NaN
//!
//! `af` travels as `doubleToRawLongBits` because a decimal rendering of a division hides a
//! divergence in the last place. The one exception is `pedigree-absent-founder`, where the golden
//! records `0xFFF8000000000000` (a negative quiet NaN) because that is what `0.0 / 0.0` produces on
//! x86-64, while the same expression on aarch64 produces `0x7FF8000000000000`. The sign bit of a
//! generated NaN is architecture-defined and not a property of the port, so two NaNs are compared
//! as equal here. Every other value is compared bit for bit.

use std::io::Read;

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::chromosome_counts::{
    calculate_chromosome_counts, called_chr_count, ChromosomeCounts, Count, Frequency,
};
use htsjdk_vcf::variant::{Genotype, VariantContext};

fn golden() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/chromosome_counts.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn reference() -> Allele {
    Allele::from_str("A", true).expect("reference allele")
}

fn alt1() -> Allele {
    Allele::from_str("C", false).expect("alternate allele")
}

fn alt2() -> Allele {
    Allele::from_str("G", false).expect("alternate allele")
}

fn gt(sample: &str, alleles: &[Allele]) -> Genotype {
    Genotype::new(sample, alleles.to_vec())
}

fn filtered(mut genotype: Genotype, filter: &str) -> Genotype {
    genotype.filters = Some(filter.to_string());
    genotype
}

fn build(alleles: &[Allele], genotypes: Vec<Genotype>) -> VariantContext {
    let mut vc = VariantContext::new("chr1", 100, alleles.to_vec());
    vc.stop = 100;
    vc.genotypes = genotypes;
    vc
}

fn names(samples: &[&str]) -> Vec<String> {
    samples.iter().map(|s| (*s).to_string()).collect()
}

struct Case {
    label: &'static str,
    vc: VariantContext,
    remove_stale: bool,
    founders: Vec<String>,
}

/// The dump's fixtures, in its order. A pedigree is reused across three founder sets, exactly as
/// the dump reuses one `VariantContext`.
fn cases() -> Vec<Case> {
    let (r, a1, a2) = (reference(), alt1(), alt2());
    let no_call = Allele::no_call();

    let case = |label, vc, remove_stale, founders: &[&str]| Case {
        label,
        vc,
        remove_stale,
        founders: names(founders),
    };

    let pedigree = || {
        build(
            &[r.clone(), a1.clone()],
            vec![
                gt("founder", &[r.clone(), a1.clone()]),
                gt("child1", &[a1.clone(), a1.clone()]),
                gt("child2", &[a1.clone(), a1.clone()]),
            ],
        )
    };

    vec![
        case(
            "no-genotypes",
            build(&[r.clone(), a1.clone()], vec![]),
            true,
            &[],
        ),
        case(
            "one-alt-het",
            build(
                &[r.clone(), a1.clone()],
                vec![gt("s1", &[r.clone(), a1.clone()])],
            ),
            true,
            &[],
        ),
        case(
            "one-alt-hom-var",
            build(
                &[r.clone(), a1.clone()],
                vec![gt("s1", &[a1.clone(), a1.clone()])],
            ),
            true,
            &[],
        ),
        case(
            "one-alt-hom-ref",
            build(
                &[r.clone(), a1.clone()],
                vec![gt("s1", &[r.clone(), r.clone()])],
            ),
            true,
            &[],
        ),
        case(
            "two-alts",
            build(
                &[r.clone(), a1.clone(), a2.clone()],
                vec![
                    gt("s1", &[r.clone(), a1.clone()]),
                    gt("s2", &[a1.clone(), a2.clone()]),
                ],
            ),
            true,
            &[],
        ),
        case(
            "ref-only",
            build(
                std::slice::from_ref(&r),
                vec![gt("s1", &[r.clone(), r.clone()])],
            ),
            true,
            &[],
        ),
        case(
            "all-no-call-remove",
            build(
                &[r.clone(), a1.clone()],
                vec![gt("s1", &[no_call.clone(), no_call.clone()])],
            ),
            true,
            &[],
        ),
        case(
            "all-no-call-keep",
            build(
                &[r.clone(), a1.clone()],
                vec![gt("s1", &[no_call.clone(), no_call.clone()])],
            ),
            false,
            &[],
        ),
        case(
            "half-no-call",
            build(
                &[r.clone(), a1.clone()],
                vec![gt("s1", &[a1.clone(), no_call.clone()])],
            ),
            true,
            &[],
        ),
        case(
            "filtered-genotype",
            build(
                &[r.clone(), a1.clone()],
                vec![
                    filtered(gt("s1", &[a1.clone(), a1.clone()]), "LowGQ"),
                    gt("s2", &[r.clone(), a1.clone()]),
                ],
            ),
            true,
            &[],
        ),
        case(
            "empty-filter",
            build(
                &[r.clone(), a1.clone()],
                vec![
                    filtered(gt("s1", &[a1.clone(), a1.clone()]), ""),
                    gt("s2", &[r.clone(), a1.clone()]),
                ],
            ),
            true,
            &[],
        ),
        case(
            "all-filtered",
            build(
                &[r.clone(), a1.clone()],
                vec![filtered(gt("s1", &[a1.clone(), a1.clone()]), "LowGQ")],
            ),
            true,
            &[],
        ),
        case("pedigree-no-founders", pedigree(), true, &[]),
        case("pedigree-one-founder", pedigree(), true, &["founder"]),
        case(
            "pedigree-two-founders",
            pedigree(),
            true,
            &["founder", "child1"],
        ),
        case("pedigree-absent-founder", pedigree(), true, &["nobody"]),
        case(
            "haploid",
            build(
                &[r.clone(), a1.clone()],
                vec![gt("s1", std::slice::from_ref(&a1))],
            ),
            true,
            &[],
        ),
        case(
            "mixed-ploidy",
            build(
                &[r.clone(), a1.clone()],
                vec![
                    gt("s1", std::slice::from_ref(&a1)),
                    gt("s2", &[r.clone(), a1.clone()]),
                    gt("s3", &[r.clone(), a1.clone(), a1.clone()]),
                ],
            ),
            true,
            &[],
        ),
        case(
            "three-hets",
            build(
                &[r.clone(), a1.clone()],
                (0..3)
                    .map(|i| gt(&format!("m{i}"), &[r.clone(), a1.clone()]))
                    .collect(),
            ),
            true,
            &[],
        ),
    ]
}

/// The golden's rows of one kind, as `label -> value`.
fn rows(text: &str, kind: &str) -> Vec<(String, String)> {
    text.lines()
        .filter_map(|line| line.strip_prefix(&format!("{kind}\t")))
        .map(|rest| {
            let (label, value) = rest.split_once('\t').expect("a label and a value");
            (label.to_string(), value.to_string())
        })
        .collect()
}

fn value(text: &str, kind: &str, label: &str) -> String {
    rows(text, kind)
        .into_iter()
        .find(|(key, _)| key == label)
        .unwrap_or_else(|| panic!("no {kind} row for {label:?}"))
        .1
}

fn show_count(counts: &ChromosomeCounts) -> String {
    match &counts.allele_count {
        None => "absent".to_string(),
        Some(Count::One(value)) => value.to_string(),
        Some(Count::Many(values)) => values
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(","),
    }
}

/// The Java class of the `AC` value, which is the type change the enum exists to model.
fn show_type(counts: &ChromosomeCounts) -> &'static str {
    match &counts.allele_count {
        None => "absent",
        Some(Count::One(_)) => "java.lang.Integer",
        Some(Count::Many(_)) => "java.util.ArrayList",
    }
}

fn frequencies(counts: &ChromosomeCounts) -> Option<Vec<f64>> {
    match &counts.allele_frequency {
        None => None,
        Some(Frequency::One(value)) => Some(vec![*value]),
        Some(Frequency::Many(values)) => Some(values.clone()),
    }
}

/// Bit equality, with two NaNs equal whatever their sign. See the module note.
fn same_bits(ours: f64, golden: i64) -> bool {
    let expected = f64::from_bits(golden as u64);
    if ours.is_nan() && expected.is_nan() {
        return true;
    }
    ours.to_bits() as i64 == golden
}

#[test]
fn every_allele_number_matches_the_reference() {
    let text = golden();
    for case in cases() {
        let counts = calculate_chromosome_counts(&case.vc, case.remove_stale, &case.founders);
        let ours = counts
            .allele_number
            .map_or_else(|| "absent".to_string(), |an| an.to_string());
        assert_eq!(
            ours,
            value(&text, "an", case.label),
            "AN for {}",
            case.label
        );
    }
}

#[test]
fn every_allele_count_matches_the_reference() {
    let text = golden();
    for case in cases() {
        let counts = calculate_chromosome_counts(&case.vc, case.remove_stale, &case.founders);
        assert_eq!(
            show_count(&counts),
            value(&text, "ac", case.label),
            "AC for {}",
            case.label
        );
    }
}

/// The type of the `AC` value, which the golden reports as a Java class name. A port that always
/// wrote a vector passes every `ac` row above and fails here.
#[test]
fn the_allele_count_changes_type_with_the_alternate_count() {
    let text = golden();
    for case in cases() {
        let counts = calculate_chromosome_counts(&case.vc, case.remove_stale, &case.founders);
        assert_eq!(
            show_type(&counts),
            value(&text, "type", case.label),
            "the type of AC for {}",
            case.label
        );
    }
    // And the two shapes really are distinct in the golden, so the assertion above has teeth.
    assert_eq!(value(&text, "type", "one-alt-het"), "java.lang.Integer");
    assert_eq!(value(&text, "type", "two-alts"), "java.util.ArrayList");
}

#[test]
fn every_allele_frequency_matches_the_reference_bit_for_bit() {
    let text = golden();
    for case in cases() {
        let counts = calculate_chromosome_counts(&case.vc, case.remove_stale, &case.founders);
        let expected = value(&text, "af", case.label);
        match frequencies(&counts) {
            None => assert_eq!(expected, "absent", "AF for {}", case.label),
            Some(values) => {
                let golden: Vec<i64> = expected
                    .split(',')
                    .map(|bits| bits.parse().expect("raw bits"))
                    .collect();
                assert_eq!(values.len(), golden.len(), "AF arity for {}", case.label);
                for (ours, bits) in values.iter().zip(&golden) {
                    assert!(
                        same_bits(*ours, *bits),
                        "AF for {}: {ours} has bits {} where the reference had {bits}",
                        case.label,
                        ours.to_bits() as i64
                    );
                }
            }
        }
    }
}

/// `getCalledChrCount()` and `getCalledChrCount(founders)` on their own, which is where the two
/// denominators come from.
#[test]
fn every_called_chromosome_count_matches_the_reference() {
    let text = golden();
    for case in cases() {
        let expected = value(&text, "called", case.label);
        let (all, founders) = expected.split_once('\t').expect("two counts");
        assert_eq!(
            called_chr_count(&case.vc, &[]).to_string(),
            all,
            "the whole cohort's called count for {}",
            case.label
        );
        assert_eq!(
            called_chr_count(&case.vc, &case.founders).to_string(),
            founders,
            "the founders' called count for {}",
            case.label
        );
    }
}

/// The rows a port gets wrong by writing what the field names suggest.
#[test]
fn the_rows_that_the_field_names_get_wrong() {
    let text = golden();

    // AF is not AC over AN. The cohort carries 5 alternate chromosomes out of 6, and the single
    // founder is a het, so AF is 1/2 while the AC written beside it is still 5.
    assert_eq!(value(&text, "ac", "pedigree-one-founder"), "5");
    assert_eq!(value(&text, "an", "pedigree-one-founder"), "6");
    assert_eq!(
        f64::from_bits(
            value(&text, "af", "pedigree-one-founder")
                .parse::<i64>()
                .unwrap() as u64
        ),
        0.5
    );

    // A founder set naming nobody present divides zero by zero, which is a NaN and not an error.
    let vc = cases()
        .into_iter()
        .find(|case| case.label == "pedigree-absent-founder")
        .expect("the fixture");
    let counts = calculate_chromosome_counts(&vc.vc, true, &vc.founders);
    assert!(matches!(counts.allele_frequency, Some(Frequency::One(af)) if af.is_nan()));
    assert_eq!(counts.allele_count, Some(Count::One(5)));

    // AN is not the ploidy times the sample count: a filtered genotype is invisible to it, and so
    // is a no-call allele.
    assert_eq!(value(&text, "an", "all-filtered"), "absent");
    assert_eq!(value(&text, "an", "half-no-call"), "1");
    // An empty FT is not a filter, so this genotype counts and the one above does not.
    assert_eq!(value(&text, "an", "empty-filter"), "4");

    // And the same site gives three keys or none depending on removeStaleValues.
    assert_eq!(value(&text, "an", "all-no-call-remove"), "absent");
    assert_eq!(value(&text, "an", "all-no-call-keep"), "0");
    assert_eq!(value(&text, "af", "all-no-call-keep"), "0");
}
