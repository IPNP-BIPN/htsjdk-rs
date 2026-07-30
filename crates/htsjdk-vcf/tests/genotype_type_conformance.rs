//! Conformance for the genotype type ladder, against the oracle.
//!
//! Golden from `tools/vcf-conformance/GenotypeTypeDump.java`.
//!
//! The rows that a port gets wrong by believing the javadoc, or by writing the obvious Rust:
//!
//! ```text
//! type   het-non-ref     HET     (two alternates, no reference allele in the genotype)
//! type   het-by-ref-flag HET     (the reference A and a non-reference A, printing as A/A)
//! mono   no-genotypes    false   (alternate alleles and nobody called: NOT monomorphic)
//! order  supplementary   a,😀,￿  (UTF-16 code units, not code points)
//! ```

use std::io::Read;

use htsjdk_vcf::allele::Allele;
use htsjdk_vcf::genotype_type::{
    compare_sample_names, determine_type, genotypes_ordered_by_name, is_monomorphic_in_samples,
    is_polymorphic_in_samples,
};
use htsjdk_vcf::variant::{Genotype, VariantContext};

fn golden() -> String {
    let path =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/genotype_type.txt.gz");
    let file = std::fs::File::open(&path).expect("golden");
    let mut text = String::new();
    flate2::read::GzDecoder::new(file)
        .read_to_string(&mut text)
        .expect("golden is gzip");
    text
}

fn allele(bases: &str, is_ref: bool) -> Allele {
    Allele::from_str(bases, is_ref).expect("an allele")
}

fn reference() -> Allele {
    allele("A", true)
}

/// The same bases as the reference, not flagged reference: a different allele to `equals`.
fn reference_bases_as_alt() -> Allele {
    allele("A", false)
}

fn alt1() -> Allele {
    allele("C", false)
}

fn alt2() -> Allele {
    allele("G", false)
}

fn gt(sample: &str, alleles: Vec<Allele>) -> Genotype {
    Genotype::new(sample, alleles)
}

fn filtered(mut genotype: Genotype, filter: &str) -> Genotype {
    genotype.filters = Some(filter.to_string());
    genotype
}

fn build(alleles: Vec<Allele>, genotypes: Vec<Genotype>) -> VariantContext {
    let mut vc = VariantContext::new("chr1", 100, alleles);
    vc.stop = 100;
    vc.genotypes = genotypes;
    vc
}

/// The dump's `type` fixtures, by label.
fn type_cases() -> Vec<(&'static str, Vec<Allele>)> {
    let (r, ru, a1, a2) = (reference(), reference_bases_as_alt(), alt1(), alt2());
    let n = Allele::no_call();
    vec![
        ("no-alleles", vec![]),
        ("hom-ref", vec![r.clone(), r.clone()]),
        ("hom-var", vec![a1.clone(), a1.clone()]),
        ("het", vec![r.clone(), a1.clone()]),
        ("het-reversed", vec![a1.clone(), r.clone()]),
        ("het-non-ref", vec![a1.clone(), a2.clone()]),
        ("het-by-ref-flag", vec![r.clone(), ru.clone()]),
        ("hom-var-unflagged-ref", vec![ru.clone(), ru.clone()]),
        ("no-call", vec![n.clone(), n.clone()]),
        ("mixed-ref", vec![r.clone(), n.clone()]),
        ("mixed-alt", vec![a1.clone(), n.clone()]),
        ("haploid-ref", vec![r.clone()]),
        ("haploid-alt", vec![a1.clone()]),
        ("haploid-no-call", vec![n.clone()]),
        ("triploid-hom-ref", vec![r.clone(), r.clone(), r.clone()]),
        ("triploid-het", vec![r.clone(), r.clone(), a1.clone()]),
        ("triploid-hom-var", vec![a1.clone(), a1.clone(), a1.clone()]),
        ("triploid-two-alts", vec![a1.clone(), a1.clone(), a2]),
        ("triploid-mixed", vec![r, a1.clone(), n]),
        (
            "tetraploid-hom-var",
            vec![a1.clone(), a1.clone(), a1.clone(), a1],
        ),
    ]
}

/// The dump's `mono` fixtures, by label.
fn mono_cases() -> Vec<(&'static str, VariantContext)> {
    let (r, a1) = (reference(), alt1());
    let n = Allele::no_call();
    vec![
        ("no-genotypes", build(vec![r.clone(), a1.clone()], vec![])),
        (
            "ref-only-site",
            build(vec![r.clone()], vec![gt("s1", vec![r.clone(), r.clone()])]),
        ),
        ("ref-only-site-no-genotypes", build(vec![r.clone()], vec![])),
        (
            "all-hom-ref",
            build(
                vec![r.clone(), a1.clone()],
                vec![
                    gt("s1", vec![r.clone(), r.clone()]),
                    gt("s2", vec![r.clone(), r.clone()]),
                ],
            ),
        ),
        (
            "one-het",
            build(
                vec![r.clone(), a1.clone()],
                vec![
                    gt("s1", vec![r.clone(), r.clone()]),
                    gt("s2", vec![r.clone(), a1.clone()]),
                ],
            ),
        ),
        (
            "all-no-call",
            build(
                vec![r.clone(), a1.clone()],
                vec![gt("s1", vec![n.clone(), n.clone()])],
            ),
        ),
        (
            "half-no-call-ref",
            build(
                vec![r.clone(), a1.clone()],
                vec![gt("s1", vec![r.clone(), n.clone()])],
            ),
        ),
        (
            "half-no-call-alt",
            build(
                vec![r.clone(), a1.clone()],
                vec![gt("s1", vec![a1.clone(), n])],
            ),
        ),
        (
            "filtered-het",
            build(
                vec![r.clone(), a1.clone()],
                vec![
                    filtered(gt("s1", vec![r.clone(), a1.clone()]), "LowGQ"),
                    gt("s2", vec![r.clone(), r.clone()]),
                ],
            ),
        ),
        (
            "hom-var-only",
            build(vec![r, a1.clone()], vec![gt("s1", vec![a1.clone(), a1])]),
        ),
    ]
}

/// The dump's `order` fixtures, by label, in the order the dump declares them (which is not their
/// sorted order: that is the point).
fn order_cases() -> Vec<(&'static str, Vec<&'static str>)> {
    vec![
        ("ascii", vec!["b", "A", "a", "B"]),
        ("digits", vec!["10", "2", "1", "20"]),
        ("underscore-vs-letters", vec!["S_1", "SA", "Sa", "S1"]),
        (
            "shared-prefix",
            vec!["sample", "sample1", "sampl", "sample10", "sample2"],
        ),
        ("non-ascii", vec!["é", "z", "Z", "É"]),
        ("supplementary", vec!["😀", "\u{ffff}", "a"]),
    ]
}

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

/// The dump's `escape`: everything outside printable ASCII travels as `\uXXXX`, one escape per
/// UTF-16 code unit, so a supplementary character arrives as its surrogate pair.
fn unescape(text: &str) -> String {
    let mut units: Vec<u16> = Vec::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            let mut buffer = [0u16; 2];
            units.extend_from_slice(c.encode_utf16(&mut buffer));
            continue;
        }
        match chars.next() {
            Some('\\') => units.push(u16::from(b'\\')),
            Some('u') => {
                let hex: String = (0..4).filter_map(|_| chars.next()).collect();
                units.push(u16::from_str_radix(&hex, 16).expect("four hex digits"));
            }
            other => panic!("unknown escape {other:?}"),
        }
    }
    String::from_utf16(&units).expect("a valid UTF-16 sequence")
}

#[test]
fn every_genotype_gets_the_reference_type_and_predicates() {
    let text = golden();
    for (label, alleles) in type_cases() {
        let genotype = gt("s1", alleles);
        let expected = value(&text, "type", label);
        // `isHetNonRef` indexes allele 1, so the dump only asks it of a ploidy of two or more.
        let het_non_ref = genotype.alleles.len() >= 2 && genotype.is_het_non_ref();
        let ours = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            determine_type(&genotype).name(),
            genotype.is_called(),
            genotype.is_hom(),
            genotype.is_hom_ref(),
            genotype.is_hom_var(),
            genotype.is_het(),
            het_non_ref,
            genotype.is_mixed(),
            genotype.is_no_call(),
            genotype.is_available(),
        );
        assert_eq!(ours, expected, "the type ladder for {label}");
    }
}

#[test]
fn every_site_gets_the_reference_monomorphism() {
    let text = golden();
    for (label, vc) in mono_cases() {
        let expected = value(&text, "mono", label);
        let ours = format!(
            "{}\t{}\t{}\t{}",
            is_monomorphic_in_samples(&vc),
            is_polymorphic_in_samples(&vc),
            htsjdk_vcf::chromosome_counts::called_chr_count_for(&vc, vc.reference(), &[]),
            htsjdk_vcf::chromosome_counts::called_chr_count(&vc, &[]),
        );
        assert_eq!(ours, expected, "monomorphism for {label}");
    }
}

#[test]
fn every_sample_order_matches_collections_sort() {
    let text = golden();
    for (label, samples) in order_cases() {
        let vc = build(
            vec![reference(), alt1()],
            samples
                .iter()
                .map(|sample| gt(sample, vec![reference(), alt1()]))
                .collect(),
        );
        let ours: Vec<String> = genotypes_ordered_by_name(&vc)
            .into_iter()
            .map(|genotype| genotype.sample_name.clone())
            .collect();
        let expected: Vec<String> = value(&text, "order", label)
            .split(',')
            .map(unescape)
            .collect();
        assert_eq!(ours, expected, "the sample order for {label}");
    }
}

/// The rows a port gets wrong by believing the javadoc, or by writing the obvious Rust.
#[test]
fn the_rows_that_the_javadoc_gets_wrong() {
    let text = golden();

    // HET is "two called alleles that are not equal", not "one ref and one alt".
    assert_eq!(
        value(&text, "type", "het-non-ref").split('\t').next(),
        Some("HET")
    );
    // And equality is bases AND the reference flag, so this genotype prints A/A and is HET.
    assert_eq!(
        value(&text, "type", "het-by-ref-flag").split('\t').next(),
        Some("HET")
    );
    let by_flag = gt("s1", vec![reference(), reference_bases_as_alt()]);
    assert_eq!(
        by_flag.alleles[0].base_string(),
        by_flag.alleles[1].base_string()
    );
    assert!(by_flag.is_het());

    // A site with alternate alleles and no genotypes is NOT monomorphic, even though both called
    // counts are zero and every "all samples are hom-ref" reading says it should be.
    assert_eq!(value(&text, "mono", "no-genotypes"), "false\ttrue\t0\t0");
    // Whereas the same zero counts with a genotype present are monomorphic.
    assert_eq!(value(&text, "mono", "all-no-call"), "true\tfalse\t0\t0");

    // The order is UTF-16 code units, so the surrogate pair of U+1F600 sorts BELOW U+FFFF. Code
    // point order, which is what Rust's own `str` comparison gives, puts it above.
    assert_eq!(
        value(&text, "order", "supplementary"),
        "a,\\ud83d\\ude00,\\uffff"
    );
    assert_eq!(
        compare_sample_names("😀", "\u{ffff}"),
        std::cmp::Ordering::Less
    );
    assert!("😀" > "\u{ffff}", "Rust's own ordering disagrees");

    // Uppercase before lowercase, and "10" before "2", neither of which is a locale's collation.
    assert_eq!(value(&text, "order", "ascii"), "A,B,a,b");
    assert_eq!(value(&text, "order", "digits"), "1,10,2,20");
}
