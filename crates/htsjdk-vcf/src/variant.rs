//! Variant contexts and genotypes.
//!
//! Ported from `htsjdk.variant.variantcontext.VariantContext`, `Genotype`, `GenotypeBuilder`
//! and `CommonInfo` at htsjdk 4.2.0, restricted to what `VCFEncoder` reads. `VariantContext`
//! is 1873 lines of which the encoder touches a narrow band; the rest (type inference,
//! subsetting, allele trimming) is ported when a tool needs it, not speculatively.
//!
//! The one piece of real logic here is `calc_vcf_genotype_keys`, which decides the FORMAT
//! column. It collects keys into a `HashSet`, sorts them, and then puts `GT` in front
//! *afterwards* - so `GT` is first and everything else is in ASCII order, and `GT` does not
//! participate in the sort at all.

use std::collections::BTreeSet;

use crate::allele::Allele;
use crate::jformat::{format_fixed, format_scientific};

/// `CommonInfo.NO_LOG10_PERROR`, the sentinel for "no QUAL".
pub const NO_LOG10_PERROR: f64 = 1.0;

/// An INFO or FORMAT attribute value.
///
/// The variants mirror the `instanceof` ladder in `VCFEncoder.formatVCFField`, in its order,
/// because that order is observable: a value is tested against `Double` before `List`, so a
/// `List<Double>` formats each element through the double rules while a bare `Double` does too,
/// and anything not matched falls through to `toString()`.
///
/// `Float` is deliberately absent. htsjdk tests `instanceof Double` only, so a `java.lang.Float`
/// attribute takes the `toString()` path and prints Java's *float* representation rather than
/// `formatVCFDouble`. Nothing in the ported surface produces one; when something does, it needs
/// its own variant rather than being folded into `Double`, which would silently change the text.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Missing,
    Int(i64),
    Double(f64),
    Bool(bool),
    Str(String),
    List(Vec<Value>),
}

impl Value {
    /// `VCFEncoder.formatVCFField`. `None` means the field is dropped entirely, which only a
    /// `false` boolean produces.
    pub fn format(&self) -> Option<String> {
        match self {
            Value::Missing => Some(".".to_string()),
            Value::Double(d) => Some(format_vcf_double(*d)),
            // Empty string for true, dropped for false. The empty string is what makes a flag
            // print as a bare key with no '='.
            Value::Bool(b) => {
                if *b {
                    Some(String::new())
                } else {
                    None
                }
            }
            Value::List(items) => {
                if items.is_empty() {
                    // formatVCFField(null), i.e. ".", not an empty string.
                    return Some(".".to_string());
                }
                let parts: Vec<String> = items
                    .iter()
                    // A nested null formats as "."; a nested false boolean would append the
                    // literal "null" in Java, which no caller produces, so it is not modelled.
                    .map(|v| v.format().unwrap_or_else(|| ".".to_string()))
                    .collect();
                Some(parts.join(","))
            }
            Value::Int(i) => Some(i.to_string()),
            Value::Str(s) => Some(s.clone()),
        }
    }
}

/// `VCFEncoder.formatVCFDouble`.
///
/// The branch structure is `d < 1` and then `d < 0.01`, both **signed**. So the sign decides as
/// much as the magnitude does: every negative value, however large, falls through to the
/// exponent format, and `-1e10` prints as `-1.000e+10` while `1e10` prints as `10000000000.00`.
/// The `Math.abs(d) >= 1e-20` guard is the only place magnitude is taken without sign, and it
/// sends everything nearer zero than 1e-20, in either direction, to the literal string `0.00`.
pub fn format_vcf_double(d: f64) -> String {
    if d < 1.0 {
        if d < 0.01 {
            if d.abs() >= 1e-20 {
                format_scientific(d, 3)
            } else {
                "0.00".to_string()
            }
        } else {
            format_fixed(d, 3)
        }
    } else {
        format_fixed(d, 2)
    }
}

/// A genotype call for one sample.
#[derive(Debug, Clone, PartialEq)]
pub struct Genotype {
    pub sample_name: String,
    pub alleles: Vec<Allele>,
    pub phased: bool,
    pub gq: Option<i32>,
    pub dp: Option<i32>,
    pub ad: Option<Vec<i32>>,
    pub pl: Option<Vec<i32>>,
    /// `getFilters()`: already a joined string in htsjdk, not a set.
    pub filters: Option<String>,
    /// Insertion-ordered, because htsjdk's `GenotypeBuilder` uses a `LinkedHashMap`. The order
    /// does not reach the output - the FORMAT column sorts - but keeping it makes the two
    /// structures comparable when debugging a divergence.
    pub extended: Vec<(String, Value)>,
}

impl Genotype {
    pub fn new(sample_name: &str, alleles: Vec<Allele>) -> Self {
        Self {
            sample_name: sample_name.to_string(),
            alleles,
            phased: false,
            gq: None,
            dp: None,
            ad: None,
            pl: None,
            filters: None,
            extended: Vec::new(),
        }
    }

    /// `GenotypeBuilder.createMissing(sample, ploidy)`: `ploidy` no-calls.
    pub fn missing(sample_name: &str, ploidy: usize) -> Self {
        Self::new(sample_name, vec![Allele::no_call(); ploidy])
    }

    /// `isAvailable()`: `getType() != UNAVAILABLE`, and the type is unavailable exactly when
    /// there are no alleles. A genotype of no-calls *is* available.
    pub fn is_available(&self) -> bool {
        !self.alleles.is_empty()
    }

    pub fn is_filtered(&self) -> bool {
        self.filters.as_ref().is_some_and(|f| !f.is_empty())
    }

    pub fn ploidy(&self) -> usize {
        self.alleles.len()
    }

    pub fn get(&self, key: &str) -> Option<&Value> {
        self.extended.iter().find(|(k, _)| k == key).map(|(_, v)| v)
    }
}

/// One VCF data line's worth of variant.
#[derive(Debug, Clone, PartialEq)]
pub struct VariantContext {
    pub contig: String,
    pub start: i64,
    pub id: String,
    /// The reference allele first, then the alternates, as htsjdk stores them. The GT indices
    /// are positions in this list.
    pub alleles: Vec<Allele>,
    pub log10_p_error: f64,
    /// `None` means filters were never applied, which prints `.`; `Some(empty)` means they were
    /// applied and passed, which prints `PASS`. The distinction is invisible in the type system
    /// of most ports and is the difference between two different files.
    pub filters: Option<Vec<String>>,
    /// Insertion-ordered. The encoder sorts, so this order is not observable in the output; it
    /// is kept so that a divergence can be traced back to what the caller set.
    pub attributes: Vec<(String, Value)>,
    pub genotypes: Vec<Genotype>,
}

impl VariantContext {
    pub fn new(contig: &str, start: i64, alleles: Vec<Allele>) -> Self {
        Self {
            contig: contig.to_string(),
            start,
            id: ".".to_string(),
            alleles,
            log10_p_error: NO_LOG10_PERROR,
            filters: None,
            attributes: Vec::new(),
            genotypes: Vec::new(),
        }
    }

    pub fn reference(&self) -> &Allele {
        &self.alleles[0]
    }

    pub fn alternate_alleles(&self) -> &[Allele] {
        &self.alleles[1..]
    }

    /// `isVariant()`: `getType() != NO_VARIATION`, and the type is no-variation exactly when
    /// there is nothing but the reference allele.
    pub fn is_variant(&self) -> bool {
        self.alleles.len() > 1
    }

    pub fn has_log10_p_error(&self) -> bool {
        self.log10_p_error != NO_LOG10_PERROR
    }

    /// `getPhredScaledQual()`: `log10PError * -10`.
    pub fn phred_scaled_qual(&self) -> f64 {
        self.log10_p_error * -10.0
    }

    pub fn is_filtered(&self) -> bool {
        self.filters.as_ref().is_some_and(|f| !f.is_empty())
    }

    pub fn filters_were_applied(&self) -> bool {
        self.filters.is_some()
    }

    pub fn genotype(&self, sample: &str) -> Option<&Genotype> {
        self.genotypes.iter().find(|g| g.sample_name == sample)
    }

    /// `getMaxPloidy(defaultPloidy)`.
    pub fn max_ploidy(&self, default_ploidy: usize) -> usize {
        self.genotypes
            .iter()
            .map(|g| g.ploidy())
            .max()
            .unwrap_or(default_ploidy)
    }

    /// `calcVCFGenotypeKeys(header)`: the FORMAT column.
    ///
    /// The four integer fields and the genotype filter are added by *presence in any sample*,
    /// so one sample carrying DP puts DP in the FORMAT of every sample, and the samples without
    /// it print `.`. Then the whole set is sorted and `GT` is prepended, which is why the
    /// familiar `GT:AD:DP:GQ:PL` is alphabetical after the first field and not by any
    /// significance ordering.
    ///
    /// `header_has_genotyping_data` is the fallback's condition: a record where every sample is
    /// a no-call with no attributes still gets a `GT` column, because otherwise the sample
    /// columns would have no field at all.
    pub fn calc_vcf_genotype_keys(&self, header_has_genotyping_data: bool) -> Vec<String> {
        let mut keys: BTreeSet<String> = BTreeSet::new();
        let mut saw_good_gt = false;
        for g in &self.genotypes {
            for (k, _) in &g.extended {
                keys.insert(k.clone());
            }
            if g.is_available() {
                saw_good_gt = true;
            }
            if g.gq.is_some() {
                keys.insert("GQ".to_string());
            }
            if g.dp.is_some() {
                keys.insert("DP".to_string());
            }
            if g.ad.is_some() {
                keys.insert("AD".to_string());
            }
            if g.pl.is_some() {
                keys.insert("PL".to_string());
            }
            if g.is_filtered() {
                keys.insert("FT".to_string());
            }
        }

        let mut sorted: Vec<String> = keys.into_iter().collect();
        if saw_good_gt {
            sorted.insert(0, "GT".to_string());
        }
        if sorted.is_empty() && header_has_genotyping_data {
            return vec!["GT".to_string()];
        }
        sorted
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn allele(s: &str, is_ref: bool) -> Allele {
        Allele::from_str(s, is_ref).unwrap()
    }

    fn vc() -> VariantContext {
        VariantContext::new("chr1", 100, vec![allele("A", true), allele("T", false)])
    }

    /// GT is prepended after the sort, so it is first regardless of where it would sort.
    #[test]
    fn the_format_column_is_gt_then_ascii_order() {
        let mut v = vc();
        let mut g = Genotype::new("S1", vec![allele("A", true), allele("T", false)]);
        g.pl = Some(vec![0, 10, 100]);
        g.ad = Some(vec![5, 6]);
        g.dp = Some(11);
        g.gq = Some(30);
        v.genotypes.push(g);
        assert_eq!(
            v.calc_vcf_genotype_keys(true),
            ["GT", "AD", "DP", "GQ", "PL"]
        );
    }

    /// One sample carrying a field puts it in every sample's FORMAT.
    #[test]
    fn a_field_present_in_one_sample_is_in_the_format_for_all() {
        let mut v = vc();
        let mut g1 = Genotype::new("S1", vec![allele("A", true), allele("T", false)]);
        g1.dp = Some(11);
        v.genotypes.push(g1);
        v.genotypes.push(Genotype::new(
            "S2",
            vec![allele("A", true), allele("A", true)],
        ));
        assert_eq!(v.calc_vcf_genotype_keys(true), ["GT", "DP"]);
    }

    #[test]
    fn a_record_with_no_genotypes_still_gets_gt_if_the_header_has_samples() {
        assert_eq!(vc().calc_vcf_genotype_keys(true), ["GT"]);
        assert!(vc().calc_vcf_genotype_keys(false).is_empty());
    }

    /// The sign is tested, not the magnitude, so a large negative takes the exponent format.
    #[test]
    fn format_vcf_double_branches_on_the_signed_value() {
        assert_eq!(format_vcf_double(1e10), "10000000000.00");
        assert_eq!(format_vcf_double(-1e10), "-1.000e+10");
        assert_eq!(format_vcf_double(0.5), "0.500");
        assert_eq!(format_vcf_double(-0.5), "-5.000e-01");
        assert_eq!(format_vcf_double(0.001), "1.000e-03");
    }

    /// Everything nearer zero than 1e-20, in either direction, is the literal "0.00".
    #[test]
    fn very_small_values_collapse_to_a_literal() {
        assert_eq!(format_vcf_double(1e-21), "0.00");
        assert_eq!(format_vcf_double(-1e-21), "0.00");
        assert_eq!(format_vcf_double(0.0), "0.00");
        assert_eq!(format_vcf_double(-0.0), "0.00");
        assert_eq!(
            format_vcf_double(1e-20),
            "1.000e-20",
            "the boundary is inclusive"
        );
    }

    #[test]
    fn a_false_flag_is_dropped_and_a_true_flag_is_empty() {
        assert_eq!(Value::Bool(false).format(), None);
        assert_eq!(Value::Bool(true).format(), Some(String::new()));
    }

    #[test]
    fn an_empty_list_formats_as_missing_not_as_nothing() {
        assert_eq!(Value::List(vec![]).format().as_deref(), Some("."));
    }
}
