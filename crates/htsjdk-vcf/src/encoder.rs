//! Encoding a variant as a VCF data line.
//!
//! Ported from `htsjdk.variant.vcf.VCFEncoder` at htsjdk 4.2.0.
//!
//! Two orderings decide the bytes, and neither is stated in the VCF specification:
//!
//!   * **INFO** is written from `new TreeMap<>()`, so the keys come out in ASCII order of the
//!     key, discarding whatever order the caller set them in.
//!   * **FILTER** is written from `ParsingUtils.sortList`, so the filters are sorted too.
//!
//! Taken with the two already recorded, that makes four distinct ordering rules in one library:
//! SAM header attributes keep insertion order (decision 0009); VCF header lines sort by the
//! rendered text of the whole line (decision 0016); VCF INFO sorts by key alone; and the FORMAT
//! column sorts, except for `GT` which is prepended after the sort. A port that picked one
//! convention and applied it throughout would be wrong three times out of four, and every one
//! of those wrongs produces a file that reads back correctly.
//!
//! The third thing worth knowing is the **trailing-field trim**. Per sample, trailing fields
//! whose value is entirely `.` and `,` are removed, so two samples on the same record can end
//! up with different numbers of colons. `isMissingValue` counts those two characters and
//! compares to the length, which also makes the *empty string* a missing value.

use std::collections::BTreeMap;

use crate::allele::Allele;
use crate::header::{HeaderLine, VcfHeader};
use crate::variant::{Genotype, VariantContext};

/// Missing-from-header handling, mirroring `allowMissingFieldsInHeader`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MissingFields {
    /// `allowMissingFieldsInHeader = false`, htsjdk's default: refuse.
    #[default]
    Refuse,
    Allow,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EncodeError {
    /// `fieldIsMissingFromHeaderError`.
    MissingFromHeader {
        key: String,
        field: &'static str,
        contig: String,
        start: i64,
    },
    /// A genotype's allele is not one of the record's alleles, so it has no GT index.
    AlleleNotInVariant(String),
    /// `"GTs cannot be missing for some samples if they are available for others in the record"`.
    UnavailableGenotype(String),
}

pub struct VcfEncoder<'a> {
    header: &'a VcfHeader,
    missing_fields: MissingFields,
    /// `outputTrailingFormatFields`. htsjdk defaults it to false, i.e. trailing missing fields
    /// *are* trimmed.
    output_trailing_format_fields: bool,
}

impl<'a> VcfEncoder<'a> {
    pub fn new(header: &'a VcfHeader) -> Self {
        Self {
            header,
            missing_fields: MissingFields::default(),
            output_trailing_format_fields: false,
        }
    }

    pub fn with_missing_fields(mut self, m: MissingFields) -> Self {
        self.missing_fields = m;
        self
    }

    pub fn with_trailing_format_fields(mut self, yes: bool) -> Self {
        self.output_trailing_format_fields = yes;
        self
    }

    fn check_header(
        &self,
        present: bool,
        key: &str,
        field: &'static str,
        vc: &VariantContext,
    ) -> Result<(), EncodeError> {
        if present || self.missing_fields == MissingFields::Allow {
            Ok(())
        } else {
            Err(EncodeError::MissingFromHeader {
                key: key.to_string(),
                field,
                contig: vc.contig.clone(),
                start: vc.start,
            })
        }
    }

    /// `VCFEncoder.encode`.
    ///
    /// A convenience over [`Self::encode_into`] for a caller that wants one record on its own. The
    /// whole-file writer uses the other one, because a `String` per record is an allocation per
    /// record and then a copy of it.
    pub fn encode(&self, vc: &VariantContext) -> Result<String, EncodeError> {
        let mut out = String::new();
        self.encode_into(vc, &mut out)?;
        Ok(out)
    }

    /// `VCFEncoder.encode`, appending rather than returning.
    ///
    /// Byte for byte what [`Self::encode`] produced before this existed; the conformance suites are
    /// the gate on that and there are three over this function alone. What changed is only where
    /// the bytes land: straight into the caller's buffer, so a file of N records allocates once
    /// rather than N times.
    pub fn encode_into(&self, vc: &VariantContext, out: &mut String) -> Result<(), EncodeError> {
        use std::fmt::Write as _;

        out.push_str(&vc.contig);
        out.push('\t');
        // `write!` rather than `push_str(&to_string())`: the same digits, without the String.
        let _ = write!(out, "{}", vc.start);
        out.push('\t');
        out.push_str(&vc.id);
        out.push('\t');
        out.push_str(&vc.reference().display_string());
        out.push('\t');

        // ALT
        if vc.is_variant() {
            for (index, allele) in vc.alternate_alleles().iter().enumerate() {
                if index > 0 {
                    out.push(',');
                }
                out.push_str(&allele.display_string());
            }
        } else {
            out.push('.');
        }
        out.push('\t');

        // QUAL
        if vc.has_log10_p_error() {
            out.push_str(&format_qual_value(vc.phred_scaled_qual()));
        } else {
            out.push('.');
        }
        out.push('\t');

        self.append_filter_string(vc, out)?;
        out.push('\t');

        // INFO. The TreeMap is the ordering: keys come out in ASCII order. The keys are borrowed
        // from the record rather than cloned into the map, which is the same order over the same
        // strings.
        let mut info: BTreeMap<&str, String> = BTreeMap::new();
        for (key, value) in &vc.attributes {
            self.check_header(self.header.has_info_line(key), key, "INFO", vc)?;
            if let Some(text) = value.format() {
                info.insert(key.as_str(), text);
            }
        }
        self.write_info_string(&info, out);

        // FORMAT and the sample columns.
        let keys = vc.calc_vcf_genotype_keys(!self.header.samples.is_empty());
        if !keys.is_empty() {
            for key in &keys {
                self.check_header(self.header.has_format_line(key), key, "FORMAT", vc)?;
            }
            out.push('\t');
            for (index, key) in keys.iter().enumerate() {
                if index > 0 {
                    out.push(':');
                }
                out.push_str(key);
            }
            self.append_genotype_data(vc, &keys, out)?;
        }

        Ok(())
    }

    /// `getFilterString`.
    ///
    /// The three states are distinct in the output: filtered prints the sorted filters, applied
    /// and passed prints `PASS`, never applied prints `.`.
    fn append_filter_string(
        &self,
        vc: &VariantContext,
        out: &mut String,
    ) -> Result<(), EncodeError> {
        if vc.is_filtered() {
            // Sorted by reference: the strings stay where they are and only the pointers move.
            let mut filters: Vec<&str> = vc
                .filters
                .as_deref()
                .unwrap_or_default()
                .iter()
                .map(String::as_str)
                .collect();
            for f in &filters {
                self.check_header(self.header.has_filter_line(f), f, "FILTER", vc)?;
            }
            filters.sort_unstable();
            for (index, f) in filters.iter().enumerate() {
                if index > 0 {
                    out.push(';');
                }
                out.push_str(f);
            }
        } else if vc.filters_were_applied() {
            out.push_str("PASS");
        } else {
            out.push('.');
        }
        Ok(())
    }

    /// `writeInfoString`.
    ///
    /// A key whose formatted value is empty prints bare, with no `=`. That happens for a `true`
    /// boolean, whose formatted value is the empty string. The `getCount() != 0` guard in
    /// htsjdk is a second route to the same thing, but `Number=0` is rejected by
    /// `VCFCompoundHeaderLine.validate` for every type but `Flag`, so the guard is reachable
    /// only for lines that are flags already. Verified against the oracle, which throws
    /// `Invalid count number, with fixed count the number should be 1 or higher` on the attempt.
    fn write_info_string(&self, info: &BTreeMap<&str, String>, out: &mut String) {
        if info.is_empty() {
            out.push('.');
            return;
        }
        let mut first = true;
        for (key, value) in info {
            if first {
                first = false;
            } else {
                out.push(';');
            }
            out.push_str(key);
            if !value.is_empty() {
                out.push('=');
                out.push_str(value);
            }
        }
    }

    /// `appendGenotypeData`.
    ///
    /// The samples come from the **header**, in header order, not from the record. A record
    /// carrying genotypes for a subset of samples gets a synthesised no-call for the rest, at
    /// the record's own maximum ploidy - so a record with one triploid sample gives the absent
    /// sample `././.` rather than `./.`.
    fn append_genotype_data(
        &self,
        vc: &VariantContext,
        keys: &[String],
        out: &mut String,
    ) -> Result<(), EncodeError> {
        let ploidy = vc.max_ploidy(2);
        let allele_index = build_allele_indices(vc);
        let has_gt = keys.iter().any(|k| k == "GT");

        for sample in &self.header.samples {
            out.push('\t');
            let owned;
            let g = match vc.genotype(sample) {
                Some(g) => g,
                None => {
                    owned = Genotype::missing(sample, ploidy);
                    &owned
                }
            };

            let mut attrs: Vec<String> = Vec::with_capacity(keys.len());
            for field in keys {
                if field == "GT" {
                    if !g.is_available() {
                        return Err(EncodeError::UnavailableGenotype(sample.clone()));
                    }
                    // GT is written straight to the output and never enters `attrs`, which is
                    // why the separator logic below has to know whether GT was in the keys.
                    write_gt_field(&allele_index, g, out)?;
                    continue;
                }
                let value = if field == "FT" {
                    if g.is_filtered() {
                        g.filters.clone().unwrap_or_default()
                    } else {
                        "PASS".to_string()
                    }
                } else if let Some(ints) = int_field(g, field) {
                    match ints {
                        None => ".".to_string(),
                        Some(v) => v
                            .iter()
                            .map(|i| i.to_string())
                            .collect::<Vec<_>>()
                            .join(","),
                    }
                } else {
                    match g.get(field) {
                        Some(v) => match v.format() {
                            Some(s) => s,
                            None => continue,
                        },
                        // An absent extended attribute becomes the *string* ".", which then
                        // goes through formatVCFField unchanged.
                        None => ".".to_string(),
                    }
                };
                attrs.push(value);
            }

            if !self.output_trailing_format_fields {
                while attrs.last().is_some_and(|a| is_missing_value(a)) {
                    attrs.pop();
                }
            }

            for (i, a) in attrs.iter().enumerate() {
                if i > 0 || has_gt {
                    out.push(':');
                }
                out.push_str(a);
            }
        }
        Ok(())
    }
}

/// The four fields `IntGenotypeFieldAccessors` handles, which bypass `formatVCFField` entirely.
///
/// The outer `Option` says whether this is one of those fields at all; the inner one says
/// whether the sample has a value, since an accessor returning null prints `.`.
fn int_field(g: &Genotype, field: &str) -> Option<Option<Vec<i32>>> {
    match field {
        "GQ" => Some(g.gq.map(|v| vec![v])),
        "DP" => Some(g.dp.map(|v| vec![v])),
        "AD" => Some(g.ad.clone()),
        "PL" => Some(g.pl.clone()),
        _ => None,
    }
}

/// `buildAlleleStrings`: allele to its index in the record, with the no-call mapped to `.`.
///
/// The map is keyed by allele *including its reference flag*, so a genotype holding an allele
/// with the wrong reference flag does not resolve and the encode fails. That is htsjdk's
/// behaviour and not an accident of the port: `SimpleAllele.equals` compares `isRef`.
fn build_allele_indices(vc: &VariantContext) -> Vec<(Allele, String)> {
    let mut map = vec![(Allele::no_call(), ".".to_string())];
    for (i, a) in vc.alleles.iter().enumerate() {
        map.push((a.clone(), i.to_string()));
    }
    map
}

/// `writeGtField`.
fn write_gt_field(
    index: &[(Allele, String)],
    g: &Genotype,
    out: &mut String,
) -> Result<(), EncodeError> {
    let lookup = |a: &Allele| -> Result<&str, EncodeError> {
        index
            .iter()
            // A later entry wins, matching `HashMap.put` overwriting on a duplicate allele.
            .rev()
            .find(|(k, _)| k == a)
            .map(|(_, v)| v.as_str())
            .ok_or_else(|| EncodeError::AlleleNotInVariant(a.display_string()))
    };
    out.push_str(lookup(&g.alleles[0])?);
    for a in &g.alleles[1..] {
        out.push(if g.phased { '|' } else { '/' });
        out.push_str(lookup(a)?);
    }
    Ok(())
}

/// `formatQualValue`: `%.2f` with a trailing `.00` removed.
///
/// The trim is unconditional on the *text*, not on the value, so it fires for any quality whose
/// first two decimals round to zero. QUAL 0.001 and QUAL 0 both print `0`, and QUAL 100000.0004
/// prints `100000`. The format is lossy before the trim and the trim makes it look exact.
pub fn format_qual_value(qual: f64) -> String {
    let s = crate::jformat::format_fixed(qual, 2);
    match s.strip_suffix(".00") {
        Some(trimmed) => trimmed.to_string(),
        None => s,
    }
}

/// `isMissingValue`: the string consists only of `.` and `,`.
///
/// Counting rather than matching makes the **empty string** a missing value too, since zero
/// occurrences equals a length of zero. That is what lets a formatted `true` flag be trimmed
/// off the end of a sample's fields.
pub fn is_missing_value(s: &str) -> bool {
    s.bytes().filter(|&b| b == b'.' || b == b',').count() == s.len()
}

impl VcfHeader {
    /// INFO and FORMAT share the `Compound` variant, distinguished by their `key`, exactly as
    /// `VCFCompoundHeaderLine` is shared between `VCFInfoHeaderLine` and `VCFFormatHeaderLine`.
    fn has_compound(&self, kind: &str, wanted: &str) -> bool {
        self.lines.iter().any(
            |l| matches!(l, HeaderLine::Compound { key, id, .. } if key == kind && id == wanted),
        )
    }

    pub fn has_info_line(&self, id: &str) -> bool {
        self.has_compound("INFO", id)
    }

    pub fn has_format_line(&self, id: &str) -> bool {
        self.has_compound("FORMAT", id)
    }

    /// Whether the header declares this filter.
    ///
    /// TWO SHAPES CARRY A FILTER LINE. A header built in code holds
    /// [`HeaderLine::Filter`]; a header PARSED from a file holds a [`HeaderLine::Structured`] whose
    /// key is `FILTER`, because [`crate::header_parse`] types only INFO, FORMAT and contig lines.
    /// Checking one shape and not the other made every parsed record carrying a filter other than
    /// `PASS` unencodable: the round trip refused a file it had just read, and htsjdk does not.
    pub fn has_filter_line(&self, wanted: &str) -> bool {
        self.lines.iter().any(|l| match l {
            HeaderLine::Filter { id, .. } => id == wanted,
            HeaderLine::Structured { key, fields } => {
                key == "FILTER"
                    && fields
                        .iter()
                        .any(|(name, value)| name == "ID" && value == wanted)
            }
            _ => false,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_qual_trim_hides_a_lossy_format() {
        assert_eq!(format_qual_value(30.0), "30");
        assert_eq!(format_qual_value(0.0), "0");
        assert_eq!(
            format_qual_value(0.001),
            "0",
            "not distinguishable from zero"
        );
        assert_eq!(format_qual_value(29.5), "29.50");
        assert_eq!(format_qual_value(1234.5678), "1234.57");
    }

    /// The empty string counts as missing, because the test counts characters.
    #[test]
    fn missing_values_are_dots_and_commas_and_nothing() {
        assert!(is_missing_value("."));
        assert!(is_missing_value(".,."));
        assert!(is_missing_value(","));
        assert!(is_missing_value(""));
        assert!(!is_missing_value("0"));
        assert!(!is_missing_value(".0"));
    }
}

#[cfg(test)]
mod filter_line_tests {
    use crate::reader::read_vcf;
    use crate::vcf_file::write_vcf;

    const FILE: &str = "##fileformat=VCFv4.2\n\
                        ##FILTER=<ID=LowQual,Description=\"Low quality\">\n\
                        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
                        ##contig=<ID=chr1,length=240>\n\
                        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n\
                        chr1\t20\t.\tT\tG\t.\tLowQual\t.\tGT\t0/1\n";

    #[test]
    fn a_parsed_filter_line_is_a_declaration() {
        let file = read_vcf(FILE).expect("the file reads");
        assert!(
            file.header.has_filter_line("LowQual"),
            "a FILTER line the parser typed as Structured still declares its filter"
        );
        assert!(!file.header.has_filter_line("Absent"));
    }

    #[test]
    fn a_record_carrying_a_declared_filter_writes_back_out() {
        let file = read_vcf(FILE).expect("the file reads");
        let written = write_vcf(&file.header, &file.records).expect("the record encodes");
        assert_eq!(written, FILE, "the round trip is the file it read");
    }
}
