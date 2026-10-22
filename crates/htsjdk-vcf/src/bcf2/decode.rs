//! `VariantContext.fullyDecode`, which `BCF2Writer.add` runs on every record a VCF reader made.
//!
//! A VCF codec leaves INFO and FORMAT values as the strings it read. BCF2 needs them typed, so
//! the writer decodes each through its header line first:
//!
//! - the line is looked up **as a FORMAT line first and an INFO line second**, whatever the
//!   field is (`VariantContextUtils.getMetaDataForField`), so an INFO key that is also a FORMAT
//!   ID is typed by the FORMAT line;
//! - a key with neither line is `Fully decoding VariantContext requires header line for all
//!   fields, but none was found for KEY`;
//! - a value is split on `,` with `String.split`, which drops trailing empty strings, and each
//!   piece goes through `Integer.valueOf`, `VCFUtils.parseVcfDouble`, or stays a string; `.` is
//!   `null`;
//! - at the site level, and only there, a value count other than the header's is
//!   `TribbleException.InvalidHeader`, except for unbounded and Flag lines;
//! - an INFO key named `FT` is dropped, because the loop that decodes INFO is the one that
//!   decodes genotype attributes and skips the genotype filter.
//!
//! The decoded map is a `new HashMap<>(10)` and is then copied by `VariantContextBuilder.attributes`
//! into another, so the result's attribute order is the hash order those two maps give, starting
//! from the order the codec's maps gave: see [`super::java_hash`].

use super::java_hash::{order_into, order_of, spread, string_hash, Fill};
use super::Bcf2Error;
use crate::attributes::parse_vcf_double;
use crate::genotypes_context::GenotypesContext;
use crate::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use crate::variant::{Genotype, Value, VariantContext};

/// A FORMAT or INFO line, as `VCFCompoundHeaderLine`.
#[derive(Debug, Clone, Copy)]
pub struct Compound<'a> {
    pub line: &'a HeaderLine,
    pub id: &'a str,
    pub number: Cardinality,
    pub line_type: LineType,
}

/// The position of the first INFO and FORMAT line of each ID, built once per header: a
/// GRCh38 header carries thousands of contig lines, and every field of every record looks here.
#[derive(Debug, Clone, Default)]
pub struct FieldIndex {
    /// Line key (`INFO`, `FORMAT`), then ID, so a lookup borrows both and allocates nothing.
    first: std::collections::HashMap<String, std::collections::HashMap<String, usize>>,
}

impl FieldIndex {
    pub fn new(header: &VcfHeader) -> Self {
        let mut first: std::collections::HashMap<String, std::collections::HashMap<String, usize>> =
            std::collections::HashMap::new();
        for (i, line) in header.lines.iter().enumerate() {
            if let HeaderLine::Compound { key, id, .. } = line {
                first
                    .entry(key.clone())
                    .or_default()
                    .entry(id.clone())
                    .or_insert(i);
            }
        }
        Self { first }
    }

    /// The first `key` line with this ID in `header`, which is the one `VCFHeader`'s map keeps.
    /// `header` must be the one the index was built from.
    pub fn compound<'a>(
        &self,
        header: &'a VcfHeader,
        key: &str,
        wanted: &str,
    ) -> Option<Compound<'a>> {
        let &i = self.first.get(key)?.get(wanted)?;
        match &header.lines[i] {
            line @ HeaderLine::Compound {
                id,
                number,
                line_type,
                ..
            } => Some(Compound {
                line,
                id,
                number: *number,
                line_type: *line_type,
            }),
            _ => None,
        }
    }
}

/// `VariantContextUtils.getMetaDataForField`: FORMAT, then INFO.
pub fn meta_data_for_field<'a>(
    header: &'a VcfHeader,
    index: &FieldIndex,
    field: &str,
) -> Result<Compound<'a>, Bcf2Error> {
    index
        .compound(header, "FORMAT", field)
        .or_else(|| index.compound(header, "INFO", field))
        .ok_or_else(|| {
            Bcf2Error::tribble(format!(
                "Fully decoding VariantContext requires header line for all fields, but none was \
                 found for {field}"
            ))
        })
}

/// `GenotypeLikelihoods.numLikelihoods`: `C(alleles + ploidy - 1, ploidy)`.
pub fn num_likelihoods(n_alleles: usize, ploidy: usize) -> i32 {
    if n_alleles == 0 {
        return 0;
    }
    let mut acc: u128 = 1;
    for k in 1..=ploidy as u128 {
        acc = acc * (n_alleles as u128 - 1 + k) / k;
    }
    acc as i32
}

/// `VCFCompoundHeaderLine.getCount(vc)`.
pub fn count_for(number: Cardinality, vc: &VariantContext) -> i32 {
    match number {
        Cardinality::Fixed(n) => n,
        Cardinality::Unbounded => -1,
        Cardinality::A => vc.alleles.len() as i32 - 1,
        Cardinality::R => vc.alleles.len() as i32,
        Cardinality::G => num_likelihoods(vc.alleles.len(), max_ploidy(vc)),
    }
}

/// `getMaxPloidy(2)`, without marking a lazy context decoded.
pub fn max_ploidy(vc: &VariantContext) -> usize {
    vc.genotypes
        .peek()
        .iter()
        .map(Genotype::ploidy)
        .max()
        .unwrap_or(2)
}

/// `String.split(",")`: trailing empty strings are removed, a leading one is kept.
pub fn java_split_comma(s: &str) -> Vec<&str> {
    let mut parts: Vec<&str> = s.split(',').collect();
    while parts.last() == Some(&"") {
        parts.pop();
    }
    parts
}

/// `Integer.valueOf`: an optional sign and ASCII digits, within `int`.
pub fn java_parse_int(s: &str) -> Option<i32> {
    let digits = s.strip_prefix(['+', '-']).unwrap_or(s);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    s.parse::<i32>().ok()
}

fn type_name(t: LineType) -> &'static str {
    match t {
        LineType::Integer => "Integer",
        LineType::Float => "Float",
        LineType::String => "String",
        LineType::Character => "Character",
        LineType::Flag => "Flag",
    }
}

/// `decodeOne`.
fn decode_one(
    field: &str,
    string: &str,
    format: &Compound,
    vc: &VariantContext,
) -> Result<Value, Bcf2Error> {
    if string == "." {
        return Ok(Value::Missing);
    }
    let failure = || {
        Bcf2Error::tribble(format!(
            "Could not decode field {field} with value {string} of declared type {}",
            type_name(format.line_type)
        ))
    };
    match format.line_type {
        LineType::Character | LineType::String => Ok(Value::Str(string.to_string())),
        LineType::Flag => {
            if string.eq_ignore_ascii_case("true") || string == "1" {
                Ok(Value::Bool(true))
            } else {
                Err(Bcf2Error::tribble(format!(
                    "VariantContext FLAG fields {field} cannot contain false values as seen at \
                     {}:{}",
                    vc.contig, vc.start
                )))
            }
        }
        LineType::Integer => java_parse_int(string)
            .map(|i| Value::Int(i as i64))
            .ok_or_else(failure),
        LineType::Float => parse_vcf_double(string)
            .map(Value::Double)
            .map_err(|_| failure()),
    }
}

/// `decodeValue`.
fn decode_value(
    field: &str,
    value: &Value,
    format: &Compound,
    vc: &VariantContext,
) -> Result<Value, Bcf2Error> {
    let text = match value {
        Value::Str(s) => Some(s.as_str()),
        // The codec stores a bare non-Flag key, and `key=`, as the constant ".".
        Value::Missing => Some("."),
        _ => None,
    };
    if let Some(string) = text {
        if field == "PL" {
            // `GenotypeLikelihoods.fromPLField`, an object no BCF2 encoder accepts.
            return Err(Bcf2Error::class_cast(
                "htsjdk.variant.variantcontext.GenotypeLikelihoods",
                "java.lang.Integer",
            ));
        }
        if string.contains(',') {
            return java_split_comma(string)
                .into_iter()
                .map(|piece| decode_one(field, piece, format, vc))
                .collect::<Result<Vec<_>, _>>()
                .map(Value::List);
        }
        return decode_one(field, string, format, vc);
    }
    match value {
        Value::List(items) if matches!(items.first(), Some(Value::Str(_))) => items
            .iter()
            .map(|item| match item {
                Value::Str(s) => decode_one(field, s, format, vc),
                Value::Missing => decode_one(field, ".", format, vc),
                other => Ok(other.clone()),
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Value::List),
        Value::List(items) if items.is_empty() => Err(Bcf2Error::index_out_of_bounds(0, 0)),
        other => Ok(other.clone()),
    }
}

/// `fullyDecodeAttributes`, over attributes already in the source map's iteration order.
fn decode_attributes(
    attributes: &[(String, Value)],
    header: &VcfHeader,
    index: &FieldIndex,
    lenient: bool,
    vc: &VariantContext,
) -> Result<Vec<(String, Value)>, Bcf2Error> {
    let mut decoded = Vec::with_capacity(attributes.len());
    for (field, value) in attributes {
        if field == "FT" {
            continue;
        }
        let format = meta_data_for_field(header, index, field)?;
        let value = decode_value(field, value, &format, vc)?;
        if value != Value::Missing
            && !lenient
            && format.number != Cardinality::Unbounded
            && format.line_type != LineType::Flag
        {
            let observed = match &value {
                Value::List(items) => items.len() as i32,
                _ => 1,
            };
            let expected = count_for(format.number, vc);
            if observed != expected {
                return Err(Bcf2Error::invalid_header(format!(
                    "Discordant field size detected for field {field} at {}:{}.  Field had \
                     {observed} values but the header says this should have {expected} values \
                     based on header record {}",
                    vc.contig,
                    vc.start,
                    format.line.render()
                )));
            }
        }
        decoded.push((field.clone(), value));
    }
    Ok(decoded)
}

/// Reorder `attributes` through the given maps.
pub fn reorder(attributes: &[(String, Value)], fills: &[Fill]) -> Vec<(String, Value)> {
    let keys: Vec<&str> = attributes.iter().map(|(k, _)| k.as_str()).collect();
    order_of(&keys, fills)
        .into_iter()
        .map(|i| attributes[i].clone())
        .collect()
}

/// [`reorder`] of attributes the caller is done with, moved rather than copied, and for an INFO
/// column of the usual size without a single allocation: the order is worked out in stack
/// buffers and applied in place.
pub(crate) fn reorder_owned(
    mut attributes: Vec<(String, Value)>,
    fills: &[Fill],
) -> Vec<(String, Value)> {
    const SMALL: usize = 16;
    let n = attributes.len();
    if n > SMALL {
        let keys: Vec<&str> = attributes.iter().map(|(k, _)| k.as_str()).collect();
        let order = order_of(&keys, fills);
        let mut slots: Vec<Option<(String, Value)>> = attributes.into_iter().map(Some).collect();
        return order
            .into_iter()
            .map(|i| slots[i].take().expect("each key once"))
            .collect();
    }
    let mut keys = [""; SMALL];
    let mut hashes = [0u32; SMALL];
    for (i, (key, _)) in attributes.iter().enumerate() {
        keys[i] = key.as_str();
        hashes[i] = spread(string_hash(key));
    }
    let mut order = [0usize; SMALL];
    let kept = order_into(&keys[..n], &hashes[..n], fills, &mut order);
    // `new[i] = old[order[i]]` in place: follow each source past the slots already filled.
    for i in 0..kept {
        let mut from = order[i];
        while from < i {
            from = order[from];
        }
        attributes.swap(i, from);
    }
    attributes.truncate(kept);
    attributes
}

/// The maps a VCF line's INFO passes through before a tool sees it: `parseInfo`'s
/// `new HashMap<>()`, then `VariantContextBuilder.attributes`' copy.
pub const VCF_CODEC_INFO: [Fill; 2] = [Fill::Put, Fill::PutAll];

/// `vc.fullyDecode(header, false)` on a record the VCF codec made, whose attributes are in file
/// order. The result's attributes are in the order `BCF2Writer.buildInfo` will meet them.
pub fn fully_decode(vc: &VariantContext, header: &VcfHeader) -> Result<VariantContext, Bcf2Error> {
    fully_decode_indexed(vc, header, &FieldIndex::new(header))
}

/// [`fully_decode`] with the header's [`FieldIndex`] already built.
pub fn fully_decode_indexed(
    vc: &VariantContext,
    header: &VcfHeader,
    index: &FieldIndex,
) -> Result<VariantContext, Bcf2Error> {
    let source = reorder(&vc.attributes, &VCF_CODEC_INFO);
    let decoded = decode_attributes(&source, header, index, false, vc)?;
    let attributes = reorder(&decoded, &[Fill::PutSized(10), Fill::PutAll]);

    // `fullyDecodeGenotypes` reads the list, which is the decode of a lazy context.
    let mut genotypes = Vec::with_capacity(vc.genotypes.len());
    for g in vc.genotypes.iter() {
        let mut g = g.clone();
        g.extended = decode_attributes(&g.extended, header, index, true, vc)?;
        genotypes.push(g);
    }

    Ok(VariantContext {
        contig: vc.contig.clone(),
        start: vc.start,
        stop: vc.stop,
        id: vc.id.clone(),
        alleles: vc.alleles.clone(),
        log10_p_error: vc.log10_p_error,
        filters: vc.filters.clone(),
        attributes,
        genotypes: GenotypesContext::new(genotypes),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The moving reorder against the copying one, through both of its paths and with repeats.
    #[test]
    fn reorder_owned_orders_as_reorder_does() {
        let mut seed = 5u64;
        let mut next = move |bound: u64| {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (seed >> 33) % bound
        };
        for _ in 0..300 {
            let n = next(30) as usize;
            let attributes: Vec<(String, Value)> = (0..n)
                .map(|i| (format!("K{}", next(40)), Value::Int(i as i64)))
                .collect();
            for fills in [
                &[Fill::PutSized(n), Fill::PutAll][..],
                &[Fill::Put, Fill::PutSized(10), Fill::PutAll][..],
            ] {
                assert_eq!(
                    reorder_owned(attributes.clone(), fills),
                    reorder(&attributes, fills),
                    "{attributes:?}"
                );
            }
        }
    }

    #[test]
    fn likelihood_counts() {
        assert_eq!(num_likelihoods(2, 2), 3);
        assert_eq!(num_likelihoods(3, 2), 6);
        assert_eq!(num_likelihoods(4, 2), 10);
        assert_eq!(num_likelihoods(2, 1), 2);
        assert_eq!(num_likelihoods(1, 5), 1);
        assert_eq!(num_likelihoods(3, 3), 10);
    }

    #[test]
    fn split_drops_trailing_empties_only() {
        assert_eq!(java_split_comma("1,"), ["1"]);
        assert_eq!(java_split_comma(",1"), ["", "1"]);
        assert_eq!(java_split_comma("1,,2"), ["1", "", "2"]);
        assert!(java_split_comma(",").is_empty());
    }

    #[test]
    fn integers_are_javas() {
        assert_eq!(java_parse_int("+5"), Some(5));
        assert_eq!(java_parse_int("-0"), Some(0));
        assert_eq!(java_parse_int("1.0"), None);
        assert_eq!(java_parse_int("2147483648"), None);
        assert_eq!(java_parse_int("+"), None);
    }
}
