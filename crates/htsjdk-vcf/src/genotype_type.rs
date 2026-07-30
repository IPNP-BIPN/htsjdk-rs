//! The genotype type ladder, ported from `Genotype.determineType`,
//! `VariantContext.isMonomorphicInSamples` and `GenotypesContext.iterateInSampleNameOrder`
//! (htsjdk 4.2.0).
//!
//! This is what every consumer downstream means by "het" and by "called", and three of the answers
//! contradict the javadoc sitting above them.
//!
//! # `HET` is "two called alleles that are not equal"
//!
//! The javadoc says "heterozygous, with at least one ref and at least one alt in any order". The
//! code says something else:
//!
//! ```java
//! return sawMultipleAlleles ? GenotypeType.HET : firstCallAllele.isReference() ? HOM_REF : HOM_VAR;
//! ```
//!
//! `C/G`, which holds no reference allele at all, is `HET`. That is exactly why `isHetNonRef`
//! exists as a separate question, and why a port that defined `HET` as "one ref and one alt" would
//! disagree on every multiallelic site.
//!
//! # Equality is bases **and** the reference flag
//!
//! `Allele.equals` compares both, so a genotype holding the reference `A` and a non-reference `A`
//! has two unequal alleles and comes back `HET`, while printing as `A/A`. The bases alone are not
//! the allele.
//!
//! # Ploidy is not two
//!
//! A haploid call is `HOM_REF` or `HOM_VAR` and never `HET`. A triploid `A/A/C` is `HET`. `MIXED`
//! is any genotype with at least one call and at least one no-call, whatever its ploidy. And the
//! empty allele list, which is the only way to reach `UNAVAILABLE`, is what `isAvailable` answers.
//!
//! # `isMonomorphicInSamples` is not "every genotype is hom-ref"
//!
//! ```java
//! monomorphic = !isVariant() || (hasGenotypes() && getCalledChrCount(getReference()) == getCalledChrCount());
//! ```
//!
//! A site carrying alternate alleles and **no genotypes** is therefore *not* monomorphic: the first
//! disjunct needs there to be no alternate allele and the second needs genotypes. Since this is the
//! guard on GATK's `SampleList` annotation, getting it backwards changes what a VCF says at every
//! sites-only record.
//!
//! # The name order is `String.compareTo`
//!
//! `getGenotypesOrderedByName` sorts the names with `Collections.sort`, which is UTF-16 code-unit
//! order: uppercase before lowercase, digits before letters, and `"10"` before `"2"`. It is not any
//! locale's collation, and it is not code-point order either. The golden settles the second half:
//!
//! ```text
//! order  supplementary  a,😀,￿
//! ```
//!
//! `U+1F600` is a supplementary character, so a port comparing code points (which is what Rust's
//! `str` ordering does) puts it **above** `U+FFFF`. Java compares UTF-16 code units, and the first
//! unit of the pair is the surrogate `0xD83D`, which is below `0xFFFF`. The two orders disagree,
//! and the reference's is the one that reaches the file.

use std::cmp::Ordering;

use crate::allele::Allele;
use crate::chromosome_counts::{called_chr_count, called_chr_count_for};
use crate::variant::{Genotype, VariantContext};

/// `htsjdk.variant.variantcontext.GenotypeType`, in declaration order. The order is observable:
/// `VariantContext.calculateGenotypeCounts` indexes an array by `ordinal()`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenotypeType {
    NoCall,
    HomRef,
    Het,
    HomVar,
    Unavailable,
    Mixed,
}

impl GenotypeType {
    /// The enum constant's name, as Java prints it.
    pub fn name(self) -> &'static str {
        match self {
            GenotypeType::NoCall => "NO_CALL",
            GenotypeType::HomRef => "HOM_REF",
            GenotypeType::Het => "HET",
            GenotypeType::HomVar => "HOM_VAR",
            GenotypeType::Unavailable => "UNAVAILABLE",
            GenotypeType::Mixed => "MIXED",
        }
    }

    /// `ordinal()`, which is the index `calculateGenotypeCounts` writes to.
    pub fn ordinal(self) -> usize {
        match self {
            GenotypeType::NoCall => 0,
            GenotypeType::HomRef => 1,
            GenotypeType::Het => 2,
            GenotypeType::HomVar => 3,
            GenotypeType::Unavailable => 4,
            GenotypeType::Mixed => 5,
        }
    }
}

/// `Genotype.determineType`, transcribed including the walk that only remembers the **first**
/// called allele: every later one is compared against that one, not against its predecessor.
pub fn determine_type(genotype: &Genotype) -> GenotypeType {
    if genotype.alleles.is_empty() {
        return GenotypeType::Unavailable;
    }

    let mut saw_no_call = false;
    let mut saw_multiple_alleles = false;
    let mut first_call_allele: Option<&Allele> = None;

    for allele in &genotype.alleles {
        if allele.is_no_call() {
            saw_no_call = true;
        } else if first_call_allele.is_none() {
            first_call_allele = Some(allele);
        } else if Some(allele) != first_call_allele {
            saw_multiple_alleles = true;
        }
    }

    if saw_no_call {
        return match first_call_allele {
            None => GenotypeType::NoCall,
            Some(_) => GenotypeType::Mixed,
        };
    }

    // The reference reaches an IllegalStateException here, which is unreachable: the list is not
    // empty and nothing in it was a no-call, so a called allele was seen.
    let first = first_call_allele.expect("a called allele, since the list is non-empty");

    if saw_multiple_alleles {
        GenotypeType::Het
    } else if first.is_reference() {
        GenotypeType::HomRef
    } else {
        GenotypeType::HomVar
    }
}

impl Genotype {
    pub fn genotype_type(&self) -> GenotypeType {
        determine_type(self)
    }

    /// `isCalled()`: neither `NO_CALL` nor `UNAVAILABLE`. A partially called genotype is `MIXED`,
    /// and therefore **called**.
    pub fn is_called(&self) -> bool {
        !matches!(
            self.genotype_type(),
            GenotypeType::NoCall | GenotypeType::Unavailable
        )
    }

    pub fn is_hom(&self) -> bool {
        self.is_hom_ref() || self.is_hom_var()
    }

    pub fn is_hom_ref(&self) -> bool {
        self.genotype_type() == GenotypeType::HomRef
    }

    pub fn is_hom_var(&self) -> bool {
        self.genotype_type() == GenotypeType::HomVar
    }

    pub fn is_het(&self) -> bool {
        self.genotype_type() == GenotypeType::Het
    }

    /// `isHetNonRef()`: het **and** the first two alleles are both non-reference. It indexes
    /// allele 1 unconditionally, so it can only be asked of a genotype with a ploidy of at least
    /// two; a haploid genotype is never `HET`, so the reference never reaches the index.
    pub fn is_het_non_ref(&self) -> bool {
        self.is_het()
            && self.alleles.len() >= 2
            && !self.alleles[0].is_reference()
            && !self.alleles[1].is_reference()
    }

    pub fn is_no_call(&self) -> bool {
        self.genotype_type() == GenotypeType::NoCall
    }

    pub fn is_mixed(&self) -> bool {
        self.genotype_type() == GenotypeType::Mixed
    }
}

/// `VariantContext.isMonomorphicInSamples`.
///
/// Note the shape: a site with alternate alleles and no genotypes is **not** monomorphic.
pub fn is_monomorphic_in_samples(vc: &VariantContext) -> bool {
    !vc.is_variant()
        || (!vc.genotypes.is_empty()
            && called_chr_count_for(vc, vc.reference(), &[]) == called_chr_count(vc, &[]))
}

/// `isPolymorphicInSamples()`.
pub fn is_polymorphic_in_samples(vc: &VariantContext) -> bool {
    !is_monomorphic_in_samples(vc)
}

/// `String.compareTo`: UTF-16 code units, then length. Not UTF-8 byte order, which disagrees on
/// every string containing a supplementary character, and not any locale's collation.
pub fn compare_sample_names(left: &str, right: &str) -> Ordering {
    for (a, b) in left.encode_utf16().zip(right.encode_utf16()) {
        if a != b {
            return a.cmp(&b);
        }
    }
    left.encode_utf16()
        .count()
        .cmp(&right.encode_utf16().count())
}

/// `GenotypesContext.getSampleNamesOrderedByName()`: the names, sorted with `Collections.sort`.
pub fn sample_names_ordered_by_name(vc: &VariantContext) -> Vec<String> {
    let mut names: Vec<String> = vc
        .genotypes
        .iter()
        .map(|genotype| genotype.sample_name.clone())
        .collect();
    names.sort_by(|a, b| compare_sample_names(a, b));
    names
}

/// `VariantContext.getGenotypesOrderedByName()`.
///
/// The reference iterates the sorted **names** and looks each one up, so a duplicated sample name
/// yields the same genotype twice and hides the other. That is transcribed rather than fixed.
pub fn genotypes_ordered_by_name(vc: &VariantContext) -> Vec<&Genotype> {
    sample_names_ordered_by_name(vc)
        .into_iter()
        .filter_map(|name| {
            vc.genotypes
                .iter()
                .find(|genotype| genotype.sample_name == name)
        })
        .collect()
}
