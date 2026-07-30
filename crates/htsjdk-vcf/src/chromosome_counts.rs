//! `AC`, `AF` and `AN`, ported from `VariantContext.getCalledChrCount` and
//! `VariantContextUtils.calculateChromosomeCounts` (htsjdk 4.2.0).
//!
//! GATK's `ChromosomeCounts` annotation is two lines that delegate straight here, so this is where
//! the three commonest INFO fields in any VCF are actually decided, and three of its behaviours are
//! not what the field names suggest.
//!
//! # `AC` and `AF` change **type** with the number of alternate alleles
//!
//! ```java
//! attributes.put(ALLELE_COUNT_KEY, alleleCounts.size() == 1 ? alleleCounts.get(0) : alleleCounts);
//! ```
//!
//! One alternate allele gives a scalar; two or more give a list. The encoded line differs by more
//! than a comma, because a scalar and a one-element list are written the same way but a consumer
//! that fetched `AC` gets an `Integer` in one case and an `ArrayList` in the other. That is why
//! [`Count`] is an enum rather than always a vector.
//!
//! # A filtered genotype contributes nothing, and a no-call contributes nothing
//!
//! `getCalledChrCount` skips a genotype whose `FT` is set, and within a kept genotype it counts
//! alleles that are not no-calls. So `AN` is not the ploidy times the sample count: a `./.` sample
//! and a sample filtered by `FT` are both invisible to it, for different reasons.
//!
//! # `AF` is not `AC / AN`
//!
//! The numerator is the founders' count for that allele and the denominator is the founders'
//! *total* called count, both computed over `founderIds`. With no founders declared, which is what
//! GATK's annotation passes, the two sets coincide and `AF` is `AC / AN` after all. With founders
//! declared they diverge, and the `AC` reported is still the **whole cohort's** while the `AF` is
//! the founders'. A port that computed `AF` from the `AC` it had just written would agree on every
//! unrelated cohort and disagree on every pedigree.
//!
//! There is also a division that cannot be reached: the `AN == 0` guard inside the loop is dead,
//! because `AN == 0` has already returned at the top of the function whenever stale values are
//! being removed, and when they are not the loop is still entered. It is transcribed anyway, since
//! a caller passing `removeStaleValues = false` reaches it.

use crate::allele::Allele;
use crate::variant::{Value, VariantContext};

/// `VCFConstants.ALLELE_NUMBER_KEY`.
pub const ALLELE_NUMBER_KEY: &str = "AN";
/// `VCFConstants.ALLELE_COUNT_KEY`.
pub const ALLELE_COUNT_KEY: &str = "AC";
/// `VCFConstants.ALLELE_FREQUENCY_KEY`.
pub const ALLELE_FREQUENCY_KEY: &str = "AF";

/// `ChromosomeCounts.keyNames`, in the order the annotation declares them, which is **not** the
/// order they are written in: the encoder sorts.
pub const KEY_NAMES: [&str; 3] = [ALLELE_NUMBER_KEY, ALLELE_COUNT_KEY, ALLELE_FREQUENCY_KEY];

/// `getCalledChrCount(sampleIds)`: chromosomes carrying any allele, no-calls excluded.
///
/// An empty `sample_ids` means **every** sample, not none. That reading is load-bearing: the
/// annotation passes an empty set and would otherwise report zero everywhere.
pub fn called_chr_count(vc: &VariantContext, sample_ids: &[String]) -> i32 {
    let mut n = 0;
    for genotype in selected(vc, sample_ids) {
        if genotype.is_filtered() {
            continue;
        }
        for allele in &genotype.alleles {
            if !allele.is_no_call() {
                n += 1;
            }
        }
    }
    n
}

/// `getCalledChrCount(a, sampleIds)`: chromosomes carrying one particular allele.
///
/// `Genotype.countAllele` compares with `Allele.equals`, which is bases **and** the reference flag,
/// so an allele that is the reference in one context and an alternate in another does not match
/// across them.
pub fn called_chr_count_for(vc: &VariantContext, allele: &Allele, sample_ids: &[String]) -> i32 {
    let mut n = 0;
    for genotype in selected(vc, sample_ids) {
        if genotype.is_filtered() {
            continue;
        }
        n += genotype
            .alleles
            .iter()
            .filter(|candidate| *candidate == allele)
            .count() as i32;
    }
    n
}

fn selected<'a>(
    vc: &'a VariantContext,
    sample_ids: &'a [String],
) -> impl Iterator<Item = &'a crate::variant::Genotype> {
    vc.genotypes
        .iter()
        .filter(move |genotype| sample_ids.is_empty() || sample_ids.contains(&genotype.sample_name))
}

/// An `AC` or `AF` value, which is a scalar for one alternate allele and a list for more.
#[derive(Debug, Clone, PartialEq)]
pub enum Count {
    One(i32),
    Many(Vec<i32>),
}

/// An `AF` value, same shape.
#[derive(Debug, Clone, PartialEq)]
pub enum Frequency {
    One(f64),
    Many(Vec<f64>),
}

/// What `calculateChromosomeCounts` decided.
///
/// The three fields are `Option` because the function **removes** keys as well as adding them: a
/// site with no alternate alleles has its `AC` and `AF` deleted from whatever the caller had, and a
/// site where nobody is called has all three deleted when stale values are being removed.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ChromosomeCounts {
    pub allele_number: Option<i32>,
    pub allele_count: Option<Count>,
    pub allele_frequency: Option<Frequency>,
    /// Whether the three keys should be **removed** from the caller's existing attributes. It is
    /// the same signal as "all three are `None`", named so a caller cannot mistake "not computed"
    /// for "computed as absent".
    pub remove_stale: bool,
}

/// `VariantContextUtils.calculateChromosomeCounts`.
///
/// `founder_ids` empty is the annotation's own call and means the whole cohort.
pub fn calculate_chromosome_counts(
    vc: &VariantContext,
    remove_stale_values: bool,
    founder_ids: &[String],
) -> ChromosomeCounts {
    let an = called_chr_count(vc, &[]);

    // Everyone is a no-call: the three keys are removed rather than written as zero.
    if an == 0 && remove_stale_values {
        return ChromosomeCounts {
            remove_stale: true,
            ..Default::default()
        };
    }

    if vc.genotypes.is_empty() {
        return ChromosomeCounts::default();
    }

    let mut counts = ChromosomeCounts {
        allele_number: Some(an),
        ..Default::default()
    };

    if vc.alternate_alleles().is_empty() {
        // No alternate allele: AC and AF are removed, and AN stays.
        return counts;
    }

    let total_founders = f64::from(called_chr_count(vc, founder_ids));
    let mut allele_counts = Vec::new();
    let mut allele_freqs = Vec::new();
    for allele in vc.alternate_alleles() {
        let founders_alt = called_chr_count_for(vc, allele, founder_ids);
        // The whole cohort's count, not the founders': the two differ as soon as founders are
        // declared, and only AF uses the founders'.
        allele_counts.push(called_chr_count_for(vc, allele, &[]));
        if an == 0 {
            allele_freqs.push(0.0);
        } else {
            allele_freqs.push(f64::from(founders_alt) / total_founders);
        }
    }

    counts.allele_count = Some(if allele_counts.len() == 1 {
        Count::One(allele_counts[0])
    } else {
        Count::Many(allele_counts)
    });
    counts.allele_frequency = Some(if allele_freqs.len() == 1 {
        Frequency::One(allele_freqs[0])
    } else {
        Frequency::Many(allele_freqs)
    });
    counts
}

impl ChromosomeCounts {
    /// The attributes as the annotation hands them on, in `ChromosomeCounts.keyNames` order.
    ///
    /// Empty for a site the function only removed keys from, which is how the annotation reports
    /// "nothing to say" as distinct from "zero".
    pub fn attributes(&self) -> Vec<(String, Value)> {
        let mut out = Vec::new();
        if let Some(an) = self.allele_number {
            out.push((ALLELE_NUMBER_KEY.to_string(), Value::Int(an as i64)));
        }
        if let Some(ac) = &self.allele_count {
            out.push((
                ALLELE_COUNT_KEY.to_string(),
                match ac {
                    Count::One(value) => Value::Int(*value as i64),
                    Count::Many(values) => Value::List(
                        values
                            .iter()
                            .map(|value| Value::Int(*value as i64))
                            .collect(),
                    ),
                },
            ));
        }
        if let Some(af) = &self.allele_frequency {
            out.push((
                ALLELE_FREQUENCY_KEY.to_string(),
                match af {
                    Frequency::One(value) => Value::Double(*value),
                    Frequency::Many(values) => {
                        Value::List(values.iter().map(|value| Value::Double(*value)).collect())
                    }
                },
            ));
        }
        out
    }
}
