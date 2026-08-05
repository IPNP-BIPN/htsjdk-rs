//! Reading one VCF data line into a `VariantContext`, sites only.
//!
//! Ported from `htsjdk.variant.vcf.AbstractVCFCodec` (`decodeLine`, `parseVCFLine`, `parseInfo`,
//! `parseAlleles`, `checkAllele`, `parseQual`) and `htsjdk.tribble.util.ParsingUtils.split`
//! (htsjdk 4.2.0).
//!
//! Genotypes are not here. `LazyGenotypesContext` defers them until something asks, and the
//! deferral is itself observable, so it gets its own slice; this decodes the eight site columns and
//! keeps the genotype block as the raw text the splitter produced.
//!
//! # The splitter is not `String.split`
//!
//! `ParsingUtils.split(line, parts, '\t', true)` fills a **fixed-size** array, and both of its
//! peculiarities change what a record means:
//!
//!  * **a leading delimiter is skipped, not honoured.** `if (end == 0) { start = 1; ... }`, so a
//!    line beginning with a tab loses its first character instead of producing an empty first
//!    column. A record whose CHROM is empty therefore parses as a record whose CHROM is the second
//!    column's text;
//!  * **trailing columns are condensed into the last slot.** The array is
//!    `min(header.getColumnCount(), 9)` long, so every genotype column after the first arrives
//!    joined back together with tabs in `parts[8]`. That is why the count check compares against 9
//!    and not against the number of samples: the record is *always* nine tokens once a header has
//!    samples, whatever the sample count.
//!
//! # `END` decides the stop, and it is trusted
//!
//! When the INFO field carries `END`, the stop comes from it and the reference allele's length is
//! ignored, so a record can end before it starts. Only a value that fails to parse is refused.
//!
//! # A flag written as `KEY=0` disappears
//!
//! `parseInfo` looks the key up in the header, and if it is declared `Flag` and the value is the
//! string `0`, the attribute is **skipped entirely** rather than stored as false. The same file
//! read against a header that does not declare `KEY` keeps `KEY=0` as a string, so the header
//! changes the record.

use crate::allele::Allele;
use crate::header::{HeaderLine, LineType, VcfHeader};
use crate::header_parse::VcfVersion;
use crate::text_transformer::TextTransformer;
use crate::variant::{Value, VariantContext, NO_LOG10_PERROR};

/// `VCFConstants.MISSING_VALUE_v4`.
pub const MISSING_VALUE: &str = ".";
/// The eight columns every record has.
pub const NUM_STANDARD_FIELDS: usize = 8;

/// The failures the record decoder produces, kept apart by their Java class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RecordError {
    /// `TribbleException`, which `generateException` raises with a line number in front.
    Tribble(String),
    /// `TribbleException.InternalCodecException`, which the genotype layer raises.
    InternalCodec(String),
    /// `java.lang.NumberFormatException`, which `parseQual` lets through uncaught: a malformed
    /// QUAL is not a malformed VCF, it is a raw JDK failure with no line number attached.
    NumberFormat(String),
}

impl RecordError {
    pub fn class(&self) -> &'static str {
        match self {
            RecordError::Tribble(_) => "htsjdk.tribble.TribbleException",
            RecordError::InternalCodec(_) => {
                "htsjdk.tribble.TribbleException$InternalCodecException"
            }
            RecordError::NumberFormat(_) => "java.lang.NumberFormatException",
        }
    }

    pub fn message(&self) -> &str {
        match self {
            RecordError::Tribble(message)
            | RecordError::InternalCodec(message)
            | RecordError::NumberFormat(message) => message,
        }
    }
}

impl std::fmt::Display for RecordError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.class(), self.message())
    }
}

impl std::error::Error for RecordError {}

/// What one decoded data line yields.
#[derive(Debug, Clone, PartialEq)]
pub struct DecodedRecord {
    pub variant: VariantContext,
    /// `parts[8]` as the splitter produced it: every genotype column joined back with tabs. Left
    /// undecoded on purpose, since `LazyGenotypesContext` is its own slice.
    pub genotype_block: Option<String>,
}

/// `ParsingUtils.split(aString, tokens, delim, condenseTrailingTokens)`.
///
/// Returns the tokens, at most `max_tokens` of them. Reproduced rather than replaced by
/// `str::split` because two of its behaviours differ from every ordinary splitter: a delimiter at
/// position zero is skipped, and trailing tokens are condensed into the last slot.
pub fn split_condensed(
    text: &str,
    delimiter: char,
    max_tokens: usize,
    condense: bool,
) -> Vec<String> {
    let bytes: Vec<char> = text.chars().collect();
    let index_of = |from: usize| -> Option<usize> {
        (from..bytes.len()).find(|&position| bytes[position] == delimiter)
    };

    let mut tokens: Vec<String> = Vec::new();
    let mut start = 0usize;
    let mut end = index_of(0);

    // The leading-delimiter case: the first character is stepped over rather than producing an
    // empty first token.
    if end == Some(0) {
        if bytes.len() > 1 {
            start = 1;
            end = index_of(start);
        } else {
            return tokens;
        }
    }

    let take = |from: usize, to: usize| bytes[from..to].iter().collect::<String>();

    let Some(mut end) = end else {
        tokens.push(take(start, bytes.len()));
        return tokens;
    };

    // `end > 0` upstream, and `end` can only be zero at the very start, which is handled above.
    while end > 0 && tokens.len() < max_tokens {
        tokens.push(take(start, end));
        start = end + 1;
        match index_of(start) {
            Some(next) => end = next,
            None => break,
        }
    }

    if condense && tokens.len() == max_tokens {
        let tail = take(start, bytes.len());
        let last = tokens.len() - 1;
        tokens[last] = format!("{}{delimiter}{tail}", tokens[last]);
    } else if tokens.len() < max_tokens {
        tokens.push(take(start, bytes.len()));
    }

    tokens
}

/// `VCFHeader.getColumnCount`.
pub fn column_count(header: &VcfHeader) -> usize {
    NUM_STANDARD_FIELDS
        + if header.samples.is_empty() {
            0
        } else {
            header.samples.len() + 1
        }
}

/// `AbstractVCFCodec.decodeLine`.
///
/// `Ok(None)` is the header-line case: a line starting with `#` is skipped rather than refused,
/// because the codec is handed lines by a reader that has already passed the header once.
pub fn decode_line(
    line: &str,
    header: &VcfHeader,
    line_number: usize,
    version: VcfVersion,
) -> Result<Option<DecodedRecord>, RecordError> {
    if line.starts_with('#') {
        return Ok(None);
    }

    let has_genotypes = !header.samples.is_empty();
    let max_tokens = column_count(header).min(NUM_STANDARD_FIELDS + 1);
    let parts = split_condensed(line, '\t', max_tokens, true);

    let expected = if has_genotypes {
        NUM_STANDARD_FIELDS + 1
    } else {
        NUM_STANDARD_FIELDS
    };
    if parts.len() != expected {
        // The number in the message is **not** the number that was just checked. Upstream the
        // check is `header.hasGenotypingData() ? 9 : 8` and the message is
        // `header == null ? 8 : 9`, so a sites-only file, whose records are checked against 8,
        // is told it was expecting 9. A port that formats `expected` into the message agrees with
        // the check and disagrees with htsjdk.
        return Err(RecordError::Tribble(format!(
            "Line {line_number}: there aren't enough columns for line {line} (we expected {} \
             tokens, and saw {} )",
            NUM_STANDARD_FIELDS + 1,
            parts.len()
        )));
    }

    // `parseVCFLine` increments the line counter on entry, so every field error below reports one
    // more than the column check above does.
    parse_vcf_line(&parts, header, line_number + 1, version).map(Some)
}

/// `AbstractVCFCodec.parseVCFLine`.
fn parse_vcf_line(
    parts: &[String],
    header: &VcfHeader,
    line_number: usize,
    version: VcfVersion,
) -> Result<DecodedRecord, RecordError> {
    let contig = parts[0].clone();

    let pos: i64 = parts[1].parse().map_err(|_| {
        generate_exception(
            line_number,
            &format!(
                "{} is not a valid start position in the VCF format",
                parts[1]
            ),
        )
    })?;

    // Empty and `.` are different: one is a refusal, the other is "no ID".
    let id = if parts[2].is_empty() {
        return Err(generate_exception(
            line_number,
            "The VCF specification requires a valid ID field",
        ));
    } else {
        parts[2].clone()
    };

    // The reference is upper-cased before anything looks at it, so a lower-case REF is not an
    // error, it is silently rewritten.
    let reference = parts[3].to_uppercase();
    let alts = &parts[4];

    let log10_p_error = parse_qual(&parts[5])?;
    let filters = parse_filters(&parts[6], line_number)?;
    let attributes = parse_info(&parts[7], header, line_number, version)?;

    // `END` wins over the reference allele's length, and nothing checks that the result is sane.
    let stop = match attributes.iter().find(|(key, _)| key == "END") {
        Some((_, value)) => value_text(value).parse::<i64>().map_err(|_| {
            generate_exception(line_number, "the END value in the INFO field is not valid")
        })?,
        None => pos + reference.len() as i64 - 1,
    };

    let alleles = parse_alleles(&reference, alts, line_number)?;

    let mut variant = VariantContext::new(&contig, pos, alleles);
    variant.id = id;
    variant.stop = stop;
    variant.log10_p_error = log10_p_error;
    variant.filters = filters;
    variant.attributes = attributes;

    Ok(DecodedRecord {
        variant,
        genotype_block: parts.get(NUM_STANDARD_FIELDS).cloned(),
    })
}

/// `generateException`, whose wording is not the column check's.
///
/// Two different messages carry the line number in this one decoder: `generateException` writes
/// "The provided VCF file is malformed at approximately line number N" and the column check writes
/// "Line N". They are also **different numbers**: `parseVCFLine` increments the counter on entry,
/// so a field error reports the record's own line while the column check, which runs before that
/// increment, reports the line before it.
fn generate_exception(line_number: usize, message: &str) -> RecordError {
    RecordError::Tribble(format!(
        "The provided VCF file is malformed at approximately line number {line_number}: {message}"
    ))
}

/// `AbstractVCFCodec.parseQual`.
///
/// Three answers, not two: `.` is "no quality", a VCF 3 style `-1` is *also* "no quality" (within
/// an epsilon, so `-1.0` and `-0.9999999` both qualify), and anything else is divided by -10.
fn parse_qual(text: &str) -> Result<f64, RecordError> {
    if text == MISSING_VALUE {
        return Ok(NO_LOG10_PERROR);
    }
    let value: f64 = text.parse().map_err(|_| {
        // `VCFUtils.parseVcfDouble` throws a NumberFormatException, which nothing catches: it
        // reaches the caller as a raw JDK failure with no line number and no "malformed VCF".
        RecordError::NumberFormat(format!("For input string: \"{text}\""))
    })?;
    // `VCFConstants.MISSING_QUALITY_v3_DOUBLE` is -1 and the epsilon is 1e-6.
    if value < 0.0 && (value - -1.0).abs() < 1e-6 {
        return Ok(NO_LOG10_PERROR);
    }
    Ok(value / -10.0)
}

/// `VCFCodec.parseFilters`.
///
/// `None` is "no filters were applied", which is a different file from `Some(vec![])`, "filters
/// were applied and passed".
pub fn parse_filters_public(
    text: &str,
    line_number: usize,
) -> Result<Option<Vec<String>>, RecordError> {
    parse_filters(text, line_number)
}

fn parse_filters(text: &str, line_number: usize) -> Result<Option<Vec<String>>, RecordError> {
    if text == MISSING_VALUE {
        return Ok(None);
    }
    if text == "PASS" {
        return Ok(Some(Vec::new()));
    }
    if text == "0" {
        return Err(generate_exception(
            line_number,
            "0 is an invalid filter name in vcf4",
        ));
    }
    if text.is_empty() {
        return Err(generate_exception(
            line_number,
            &format!("The VCF specification requires a valid filter status: filter was {text}"),
        ));
    }
    Ok(Some(text.split(';').map(str::to_string).collect()))
}

/// `AbstractVCFCodec.parseInfo`.
fn parse_info(
    text: &str,
    header: &VcfHeader,
    line_number: usize,
    version: VcfVersion,
) -> Result<Vec<(String, Value)>, RecordError> {
    let mut attributes: Vec<(String, Value)> = Vec::new();
    // The transformer runs on the value **before** the Flag test below, so under 4.3 a declared
    // Flag written as `DB=%30` is dropped exactly as `DB=0` is. The key is never transformed.
    let transformer = TextTransformer::for_version(version);

    if text.is_empty() {
        return Err(generate_exception(
            line_number,
            "The VCF specification requires a valid (non-zero length) info field",
        ));
    }
    if text == MISSING_VALUE {
        return Ok(attributes);
    }
    if text.contains('\t') || text.contains(' ') {
        return Err(generate_exception(
            line_number,
            &format!(
                "The VCF specification does not allow for whitespace in the INFO field. Offending \
                 field value was \"{text}\""
            ),
        ));
    }

    for field in text.split(';') {
        let (key, value) = match field.find('=') {
            Some(equals) => {
                let key = field[..equals].to_string();
                let raw = &field[equals + 1..];
                let parts: Vec<&str> = raw.split(',').collect();
                if parts.len() == 1 {
                    let decoded = transformer.decode(parts[0]);
                    // A declared Flag written as `KEY=0` is dropped rather than stored.
                    if info_type(header, &key) == Some(LineType::Flag) && decoded == "0" {
                        continue;
                    }
                    (key, Value::Str(decoded))
                } else {
                    (
                        key,
                        Value::List(
                            parts
                                .iter()
                                .map(|p| Value::Str(transformer.decode(p)))
                                .collect(),
                        ),
                    )
                }
            }
            None => {
                let key = field.to_string();
                // A bare key whose header type is not Flag becomes the *string* ".", not a flag
                // and not an absence.
                match info_type(header, &key) {
                    // `VCFConstants.MISSING_VALUE_v4` itself, not a substring that reads the same.
                    // The two are indistinguishable in every rendering and are *not*
                    // indistinguishable to `getAttributeAsInt`, which tests the reference; see
                    // [`crate::attributes`].
                    Some(line_type) if line_type != LineType::Flag => (key, Value::Missing),
                    _ => (key, Value::Bool(true)),
                }
            }
        };
        // `key=` with nothing after it is assigned the missing-value **constant**, which is a
        // different object from a `.` written in the file even though both render as ".".
        let value = match &value {
            Value::Str(text) if text.is_empty() => Value::Missing,
            _ => value,
        };
        put(&mut attributes, key, value);
    }

    Ok(attributes)
}

/// The attributes go into a `HashMap`, so a repeated key keeps the **last** value. The order is
/// not observable: the encoder sorts before writing.
fn put(attributes: &mut Vec<(String, Value)>, key: String, value: Value) {
    if let Some(slot) = attributes.iter_mut().find(|(existing, _)| *existing == key) {
        slot.1 = value;
    } else {
        attributes.push((key, value));
    }
}

/// The declared type of an INFO key, which is what makes a header change a record.
fn info_type(header: &VcfHeader, key: &str) -> Option<LineType> {
    header.lines.iter().find_map(|line| match line {
        HeaderLine::Compound {
            key: line_key,
            id,
            line_type,
            ..
        } if line_key == "INFO" && id == key => Some(*line_type),
        _ => None,
    })
}

fn value_text(value: &Value) -> String {
    match value {
        Value::Str(text) => text.clone(),
        other => other.format().unwrap_or_default(),
    }
}

/// `AbstractVCFCodec.parseAlleles`.
///
/// A no-call alternate (`.`) is checked and then **not added**, so `A` with an ALT of `.` yields a
/// one-allele record rather than a two-allele one.
fn parse_alleles(
    reference: &str,
    alts: &str,
    line_number: usize,
) -> Result<Vec<Allele>, RecordError> {
    check_allele(reference, true, line_number)?;
    let mut alleles = vec![Allele::from_str(reference, true)
        .map_err(|e| generate_exception(line_number, &format!("{e:?}")))?];

    // `alts.indexOf(',') == -1` takes a different path, but both end in the same per-allele check.
    for alt in alts.split(',') {
        check_allele(alt, false, line_number)?;
        let allele = Allele::from_str(alt, false)
            .map_err(|e| generate_exception(line_number, &format!("{e:?}")))?;
        if !allele.is_no_call() {
            alleles.push(allele);
        }
    }

    Ok(alleles)
}

/// `AbstractVCFCodec.checkAllele`, whose message depends on *why* the allele is unacceptable.
fn check_allele(allele: &str, is_ref: bool, line_number: usize) -> Result<(), RecordError> {
    if allele.is_empty() {
        return Err(generate_exception(
            line_number,
            "empty alleles are not permitted in VCF records",
        ));
    }

    if crate::allele::would_be_symbolic(allele.as_bytes()) {
        if is_ref {
            return Err(generate_exception(
                line_number,
                &format!("Symbolic alleles not allowed as reference allele: {allele}"),
            ));
        }
        return Ok(());
    }

    // VCF 3 wrote insertions and deletions with a leading D or I, and the codec refuses them by
    // name rather than failing on the bases.
    let first = allele.as_bytes()[0];
    if first == b'D' || first == b'I' {
        return Err(generate_exception(
            line_number,
            "Insertions/Deletions are not supported when reading 3.x VCF's. Please convert your \
             file to VCF4 using VCFTools, available at http://vcftools.sourceforge.net/index.html",
        ));
    }

    if !crate::allele::acceptable_bases(allele.as_bytes(), is_ref) {
        return Err(generate_exception(
            line_number,
            &bad_allele_bases_text(allele),
        ));
    }

    if is_ref && allele == MISSING_VALUE {
        return Err(generate_exception(
            line_number,
            "The reference allele cannot be missing",
        ));
    }

    Ok(())
}

/// `generateExceptionTextForBadAlleleBases`.
fn bad_allele_bases_text(allele: &str) -> String {
    if allele.is_empty() {
        return "empty alleles are not permitted in VCF records".to_string();
    }
    if allele.contains('[') || allele.contains(']') || allele.contains(':') || allele.contains('.')
    {
        return "VCF support for complex rearrangements with breakends has not yet been implemented"
            .to_string();
    }
    format!("unparsable vcf record with allele {allele}")
}
