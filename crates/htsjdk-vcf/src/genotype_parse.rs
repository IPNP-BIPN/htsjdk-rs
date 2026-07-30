//! Reading the genotype columns of a VCF data line.
//!
//! Ported from `htsjdk.variant.vcf.AbstractVCFCodec.createGenotypeMap`, `parseGenotypeAlleles` and
//! `oneAllele` (htsjdk 4.2.0), on top of [`crate::record_parse`], which stops at the site columns.
//!
//! # The GT separators are three characters, not two
//!
//! `VCFConstants.PHASING_TOKENS` is `"/|\\"`, and the split goes through a `StringTokenizer`, which
//! **skips empty tokens**. Two consequences a `split('/')` port does not have: a backslash
//! separates alleles like a slash does, and `0//1`, `/0/1` and `0/1/` all yield the same two
//! alleles as `0/1`. Nothing rejects them.
//!
//! # A malformed AD or PL is silently dropped
//!
//! `decodeInts` catches `NumberFormatException` and returns **null**, which the builder stores as
//! "no AD" rather than raising. So `AD=1,x` produces a genotype with no AD at all, while `DP=x`,
//! two lines further down the same method, throws an uncaught `NumberFormatException`. Same kind of
//! malformed integer, two different outcomes, decided by which key it was under.
//!
//! # GT has to be first, and before 4.1 it has to exist
//!
//! A `GT` anywhere but position 0 is a refusal. A record with no `GT` at all is fine from VCF 4.1
//! and a refusal before it, so the same file is valid or not depending on a header line.
//!
//! # `GL` becomes `PL`, unless a `PL` was seen first
//!
//! A `GL` field is converted through `GenotypeLikelihoods.fromGLField().getAsPLs()`, which lives in
//! [`crate::genotype_likelihoods`] because it is floating-point work with its own rounding rather
//! than parsing. The gate on it is record-wide and not per sample: `plIsSet` is set by the first
//! `PL` in *any* sample and suppresses the conversion in every sample after it, so a record whose
//! first sample has a `PL` and whose second has only a `GL` leaves the second sample without one.

use crate::allele::Allele;
use crate::header::VcfHeader;
use crate::header_parse::VcfVersion;
use crate::record_parse::{split_condensed, RecordError, MISSING_VALUE, NUM_STANDARD_FIELDS};
use crate::variant::{Genotype, Value};

/// `VCFConstants.PHASING_TOKENS`.
const PHASING_TOKENS: [char; 3] = ['/', '|', '\\'];

/// Everything the genotype layer needs from the record around it.
///
/// `site_parts` is there for one reason: the "too many keys" message quotes the **site** columns
/// rather than the genotype ones, so reproducing it needs the record's own tokens.
#[derive(Clone, Copy)]
pub struct GenotypeContext<'a> {
    pub site_parts: &'a [String],
    pub header: &'a VcfHeader,
    pub version: VcfVersion,
    pub contig: &'a str,
    pub pos: i64,
    pub line_number: usize,
}

/// `AbstractVCFCodec.createGenotypeMap`.
///
/// `block` is `parts[8]`: the FORMAT column and every sample column, joined by tabs, exactly as the
/// site decoder left it.
pub fn parse_genotypes(
    block: &str,
    alleles: &[Allele],
    context: &GenotypeContext<'_>,
) -> Result<Vec<Genotype>, RecordError> {
    let GenotypeContext {
        site_parts,
        header,
        version,
        contig,
        pos,
        line_number,
    } = *context;
    // The array is sized from the header, and this split does **not** condense, so an extra column
    // is a count mismatch rather than a joined last field.
    let expected = crate::record_parse::column_count(header) - NUM_STANDARD_FIELDS;
    let parts = split_condensed(block, '\t', expected, false);
    if parts.len() != expected {
        return Err(malformed(
            line_number,
            &format!(
                "there are {} genotypes while the header requires that {} genotypes be present \
                 for all records at {contig}:{pos}",
                parts.len() as i64 - 1,
                expected as i64 - 1
            ),
        ));
    }

    let keys: Vec<&str> = parts[0].split(':').collect();
    let mut genotypes = Vec::with_capacity(parts.len() - 1);
    // `PlIsSet` is per record, not per sample: a PL in one sample suppresses the GL conversion in
    // every sample after it.
    let mut pl_is_set = false;

    for (offset, column) in parts[1..].iter().enumerate() {
        let sample = header
            .samples
            .get(offset)
            .cloned()
            .unwrap_or_else(|| String::from("<missing>"));
        let values: Vec<&str> = column.split(':').collect();

        if keys.len() < values.len() {
            // The message is wrong upstream and reproduced as it is: `values` is
            // `parts[genotypeOffset]` over the **site** columns, not over the genotype ones, so
            // the first sample's failure quotes the record's POS. `keys` is `parts[8]`, the whole
            // genotype block rather than the FORMAT column.
            let quoted_values = site_parts.get(offset + 1).cloned().unwrap_or_default();
            return Err(malformed(
                line_number,
                &format!(
                    "There are too many keys for the sample {sample}, keys = {}, values = \
                     {quoted_values}",
                    site_parts
                        .get(NUM_STANDARD_FIELDS)
                        .cloned()
                        .unwrap_or_default()
                ),
            ));
        }

        let mut genotype = Genotype::new(&sample, Vec::new());
        let mut gt_position: i64 = -1;

        for (index, key) in keys.iter().enumerate() {
            let missing = index >= values.len();
            if *key == "GT" {
                gt_position = index as i64;
                continue;
            }
            if missing {
                // A key with no value at all is skipped, which is not the same as a value of ".":
                // both leave the field unset here, but only the second one had to be looked at.
                continue;
            }
            let value = values[index];
            if *key == "FT" {
                // The genotype filter goes through the record's filter parser, so `PASS` and `.`
                // mean what they mean on a record.
                if let Some(filters) = parse_genotype_filters(value, line_number)? {
                    genotype.filters = Some(filters.join(";"));
                }
                continue;
            }
            if value == MISSING_VALUE {
                continue;
            }
            match *key {
                "GQ" => {
                    // The VCF 3 encoding of a missing GQ, tested on the **string**: "-1.0"
                    // takes the other path, rounds to -1, and is then indistinguishable from
                    // absent because -1 is `Genotype`'s own sentinel for "no GQ".
                    if value == "-1" {
                        genotype.gq = None;
                    } else {
                        let parsed: f64 = value.parse().map_err(|_| {
                            RecordError::NumberFormat(format!("For input string: \"{value}\""))
                        })?;
                        // `Math.round`, which is not `floor(x + 0.5)` however much its javadoc
                        // says so: see [`crate::genotype_likelihoods::java_round`]. It still
                        // rounds half **up** rather than half away from zero, so -1.5 gives -1.
                        genotype.gq = Some(crate::genotype_likelihoods::java_round(parsed) as i32);
                    }
                }
                "AD" => genotype.ad = decode_ints(value),
                "PL" => {
                    genotype.pl = decode_ints(value);
                    pl_is_set = true;
                }
                "GL" => {
                    // `plIsSet` is checked, not `genotype.pl`: a PL seen in an *earlier sample*
                    // suppresses this sample's conversion, which is the record-wide flag doing
                    // something a per-sample one would not.
                    if !pl_is_set {
                        genotype.pl = crate::genotype_likelihoods::gl_field_to_pls(value).map_err(
                            |error| match error {
                                crate::genotype_likelihoods::LikelihoodError::PartialMissing => {
                                    RecordError::Tribble(error.message())
                                }
                                other => RecordError::NumberFormat(other.message()),
                            },
                        )?;
                    }
                }
                "DP" => {
                    // Not `decodeInts`: this one is a bare `Integer.parseInt` and its failure is
                    // not caught anywhere.
                    genotype.dp = Some(value.parse::<i32>().map_err(|_| {
                        RecordError::NumberFormat(format!("For input string: \"{value}\""))
                    })?);
                }
                other => genotype
                    .extended
                    .push((other.to_string(), Value::Str(value.to_string()))),
            }
        }

        if !version.is_at_least(VcfVersion::Vcf4_1) && gt_position == -1 {
            return Err(malformed(
                line_number,
                "Unable to find the GT field for the record; the GT field is required before VCF4.1",
            ));
        }
        if gt_position > 0 {
            return Err(malformed(
                line_number,
                &format!(
                    "Saw GT field at position {gt_position}, but it must be at the first position \
                     for genotypes when present"
                ),
            ));
        }

        if gt_position == 0 {
            let gt = values.first().copied().unwrap_or("");
            genotype.alleles = parse_genotype_alleles(gt, alleles, contig, pos)?;
            genotype.phased = gt.contains('|');
        }

        genotypes.push(genotype);
    }

    Ok(genotypes)
}

/// `parseGenotypeAlleles`, whose separators are `/`, `|` **and** `\`, and whose tokenizer drops
/// empty tokens, so repeated or leading separators change nothing.
fn parse_genotype_alleles(
    gt: &str,
    alleles: &[Allele],
    contig: &str,
    pos: i64,
) -> Result<Vec<Allele>, RecordError> {
    let mut result = Vec::new();
    for token in gt.split(PHASING_TOKENS).filter(|t| !t.is_empty()) {
        result.push(one_allele(token, alleles, contig, pos)?);
    }
    Ok(result)
}

/// `oneAllele`.
fn one_allele(
    index: &str,
    alleles: &[Allele],
    contig: &str,
    pos: i64,
) -> Result<Allele, RecordError> {
    if index == MISSING_VALUE {
        return Ok(Allele::no_call());
    }
    let position: usize = match index.parse::<i64>() {
        Ok(value) if value >= 0 => value as usize,
        _ => {
            return Err(RecordError::InternalCodec(format!(
                "The following invalid GT allele index was encountered in the file: {index}"
            )))
        }
    };
    alleles.get(position).cloned().ok_or_else(|| {
        // The builder wraps the failure with the position, which is how a caller learns *where* the
        // undefined index was without the codec carrying a line number into the allele layer.
        let _ = (contig, pos);
        RecordError::InternalCodec(format!(
            "The allele with index {index} is not defined in the REF/ALT columns in the record"
        ))
    })
}

/// `decodeInts`, which returns **null** when any element fails to parse: the whole field is
/// dropped, silently, rather than the record being refused.
fn decode_ints(text: &str) -> Option<Vec<i32>> {
    text.split(',')
        .map(|item| item.parse::<i32>().ok())
        .collect()
}

/// The `FT` field, which reuses the record's filter rules.
fn parse_genotype_filters(
    text: &str,
    line_number: usize,
) -> Result<Option<Vec<String>>, RecordError> {
    crate::record_parse::parse_filters_public(text, line_number)
}

fn malformed(line_number: usize, message: &str) -> RecordError {
    RecordError::Tribble(format!(
        "The provided VCF file is malformed at approximately line number {line_number}: {message}"
    ))
}
