//! `BCF2Codec`, `BCF2Decoder`, `BCF2GenotypeFieldDecoders` and `BCF2LazyGenotypesDecoder`.
//!
//! A record's genotypes are decoded as it is read, and kept beside the block they came from:
//! a BCF writer whose header agrees copies the block rather than re-encoding, exactly as it
//! does upstream for a block nobody decoded. Upstream also decodes at once when the header's
//! samples are not sorted, and then there is no block to copy; that is reproduced too.

use std::sync::Arc;

use super::decode::{meta_data_for_field, FieldIndex};
use super::java_hash::Fill;
use super::types::{decode_size, decode_type_id, Bcf2Type, Cursor, OVERFLOW_ELEMENT_MARKER};
use super::{id_lines, make_dictionary, sample_names_in_order, Bcf2Error, LazyData, LazyHeader};
use crate::allele::Allele;
use crate::genotypes_context::GenotypesContext;
use crate::header::{HeaderLine, LineType, VcfHeader};
use crate::variant::{Genotype, Value, VariantContext, NO_LOG10_PERROR};

/// `MAX_HEADER_SIZE`.
const MAX_HEADER_SIZE: i32 = 0x0800_0000;

/// A decoded value: `null`, an `Integer`, a `Double`, a `String`, or a `List` of them.
#[derive(Debug, Clone, PartialEq)]
enum Typed {
    Null,
    Int(i32),
    Double(f64),
    Str(String),
    List(Vec<Typed>),
}

impl Typed {
    fn class(&self) -> &'static str {
        match self {
            Typed::Null => "null",
            Typed::Int(_) => "java.lang.Integer",
            Typed::Double(_) => "java.lang.Double",
            Typed::Str(_) => "java.lang.String",
            // `Arrays.asList` for an exploded string, an `ArrayList` otherwise.
            Typed::List(items) if matches!(items.first(), Some(Typed::Str(_))) => {
                "java.util.Arrays$ArrayList"
            }
            Typed::List(_) => "java.util.ArrayList",
        }
    }

    fn into_value(self) -> Value {
        match self {
            Typed::Null => Value::Missing,
            Typed::Int(i) => Value::Int(i as i64),
            Typed::Double(d) => Value::Double(d),
            Typed::Str(s) => Value::Str(s),
            Typed::List(items) => Value::List(items.into_iter().map(Typed::into_value).collect()),
        }
    }
}

/// `BCF2Decoder` over one block.
struct Decoder<'a> {
    cursor: Cursor<'a>,
}

fn type_for(descriptor: u8) -> Result<Bcf2Type, Bcf2Error> {
    let id = decode_type_id(descriptor);
    Bcf2Type::from_id(id)
        .ok_or_else(|| Bcf2Error::null_pointer(format!("Cannot decode BCF2 type ID {id}")))
}

impl<'a> Decoder<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self {
            cursor: Cursor::new(bytes),
        }
    }

    fn read_type_descriptor(&mut self) -> u8 {
        self.cursor.byte()
    }

    fn decode_int(&mut self, t: Bcf2Type) -> Result<i32, Bcf2Error> {
        t.read(&mut self.cursor)
    }

    fn decode_int_or(&mut self, descriptor: u8, missing: i32) -> Result<i32, Bcf2Error> {
        let t = type_for(descriptor)?;
        let i = self.decode_int(t)?;
        Ok(if i == t.missing_bytes() { missing } else { i })
    }

    fn decode_number_of_elements(&mut self, descriptor: u8) -> Result<i32, Bcf2Error> {
        if decode_size(descriptor) == OVERFLOW_ELEMENT_MARKER {
            let inner = self.read_type_descriptor();
            self.decode_int_or(inner, -1)
        } else {
            Ok(decode_size(descriptor))
        }
    }

    fn decode_typed_value(&mut self) -> Result<Typed, Bcf2Error> {
        let descriptor = self.read_type_descriptor();
        let size = self.decode_number_of_elements(descriptor)?;
        self.decode_typed_value_sized(descriptor, size)
    }

    fn decode_typed_value_sized(&mut self, descriptor: u8, size: i32) -> Result<Typed, Bcf2Error> {
        if size == 0 {
            return Ok(Typed::Null);
        }
        let t = type_for(descriptor)?;
        if t == Bcf2Type::Char {
            return Ok(self.decode_literal_string(size));
        }
        if size == 1 {
            return self.decode_single_value(t);
        }
        let mut values = Vec::new();
        for _ in 0..size.max(0) {
            match self.decode_single_value(t)? {
                Typed::Null => continue,
                v => values.push(v),
            }
        }
        Ok(if values.is_empty() {
            Typed::Null
        } else {
            Typed::List(values)
        })
    }

    fn decode_single_value(&mut self, t: Bcf2Type) -> Result<Typed, Bcf2Error> {
        let value = self.decode_int(t)?;
        if value == t.missing_bytes() {
            return Ok(Typed::Null);
        }
        Ok(match t {
            Bcf2Type::Int8 | Bcf2Type::Int16 | Bcf2Type::Int32 => Typed::Int(value),
            Bcf2Type::Float => Typed::Double(f32::from_bits(value as u32) as f64),
            Bcf2Type::Char => Typed::Int(value & 0xFF),
            Bcf2Type::Missing => {
                return Err(Bcf2Error::tribble(format!(
                    "BCF2 codec doesn't know how to decode type {t}"
                )))
            }
        })
    }

    /// `decodeLiteralString`: up to the first `\0`, exploded when it starts with a comma.
    fn decode_literal_string(&mut self, size: i32) -> Typed {
        let bytes = self.cursor.bytes(size.max(0) as usize);
        let good = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
        if good == 0 {
            return Typed::Null;
        }
        let s = String::from_utf8_lossy(&bytes[..good]).into_owned();
        if let Some(rest) = s.strip_prefix(',') {
            let mut parts: Vec<&str> = rest.split(',').collect();
            while parts.last() == Some(&"") {
                parts.pop();
            }
            Typed::List(
                parts
                    .into_iter()
                    .map(|p| Typed::Str(p.to_string()))
                    .collect(),
            )
        } else {
            Typed::Str(s)
        }
    }

    /// `decodeIntArray(size, type)`: `null` when the first value is missing, cut at the first
    /// missing value after it.
    fn decode_int_array(&mut self, size: i32, t: Bcf2Type) -> Result<Option<Vec<i32>>, Bcf2Error> {
        if size == 0 {
            return Ok(None);
        }
        let first = self.decode_int(t)?;
        if first == t.missing_bytes() {
            for _ in 1..size {
                self.decode_int(t)?;
            }
            return Ok(None);
        }
        let mut ints = vec![first];
        for i in 1..size {
            let v = self.decode_int(t)?;
            if v == t.missing_bytes() {
                for _ in i + 1..size {
                    self.decode_int(t)?;
                }
                return Ok(Some(ints));
            }
            ints.push(v);
        }
        Ok(Some(ints))
    }
}

/// What a whole-file read produced.
#[derive(Debug, Clone)]
pub struct Bcf2File {
    /// The header as `VCFCodec.readActualHeader` parsed it from the embedded text.
    pub header: VcfHeader,
    /// The version the embedded header declared, which its codec keeps.
    pub codec_version: crate::header_parse::VcfVersion,
    /// `header.getVCFHeaderVersion()`, as [`crate::reader::VcfFile`] has it.
    pub header_version: Option<crate::header_parse::VcfVersion>,
    pub records: Vec<VariantContext>,
}

/// A read that stopped: the exception, and the records decoded before it.
#[derive(Debug, Clone)]
pub struct Bcf2ReadFailure {
    pub error: Bcf2Error,
    pub records: Vec<VariantContext>,
    /// Whether the header was what failed, which `TribbleIndexedFeatureReader` wraps differently.
    pub in_header: bool,
}

impl Bcf2ReadFailure {
    /// The exception as `VCFFileReader` surfaces it for `source`, the file's URI.
    ///
    /// `TribbleIndexedFeatureReader.readHeader` wraps anything the header throws in a
    /// `MalformedFeatureFile` naming the file. `readNextRecord` names the file on a
    /// `TribbleException` and lets every other exception through as it was.
    pub fn reader_exception(&self, source: &str) -> Bcf2Error {
        if self.in_header {
            return Bcf2Error::new(
                "htsjdk.tribble.TribbleException$MalformedFeatureFile",
                format!(
                    "Unable to parse header with error: {}, for input source: {source}",
                    self.error.message
                ),
            );
        }
        if self
            .error
            .class
            .starts_with("htsjdk.tribble.TribbleException")
        {
            return Bcf2Error::new(
                self.error.class,
                format!("{}, for input source: {source}", self.error.message),
            );
        }
        self.error.clone()
    }
}

/// `BCF2Codec` after `readHeader`.
struct Codec {
    header: VcfHeader,
    fields: FieldIndex,
    codec_version: crate::header_parse::VcfVersion,
    header_version: Option<crate::header_parse::VcfVersion>,
    contig_names: Vec<String>,
    dictionary: Vec<String>,
    lazy_header: Arc<LazyHeader>,
    samples_sorted: bool,
    record_no: i32,
    pos: i32,
}

impl Codec {
    /// `error(message)`: a `TribbleException` naming the record and position.
    fn error(&self, message: impl AsRef<str>) -> Bcf2Error {
        Bcf2Error::tribble(format!(
            "{}, at record {} with position {}:",
            message.as_ref(),
            self.record_no,
            self.pos
        ))
    }

    fn dictionary_string(&self, offset: i32) -> Result<&str, Bcf2Error> {
        self.dictionary
            .get(usize::try_from(offset).unwrap_or(usize::MAX))
            .map(String::as_str)
            .ok_or_else(|| Bcf2Error::index_out_of_bounds(offset as i64, self.dictionary.len()))
    }

    /// `getDictionaryString()`: the offset read as an `Integer`.
    fn read_dictionary_string(&self, d: &mut Decoder) -> Result<String, Bcf2Error> {
        match d.decode_typed_value()? {
            Typed::Int(offset) => Ok(self.dictionary_string(offset)?.to_string()),
            Typed::Null => Err(Bcf2Error::null_pointer(
                "Cannot invoke \"java.lang.Integer.intValue()\" because the return value of \
                 \"htsjdk.variant.bcf2.BCF2Decoder.decodeTypedValue()\" is null",
            )),
            other => Err(Bcf2Error::class_cast(other.class(), "java.lang.Integer")),
        }
    }

    /// The sites block: the record without its genotypes, and the genotype field count.
    fn decode_sites(&mut self, sites: &[u8]) -> Result<(VariantContext, i32), Bcf2Error> {
        let mut d = Decoder::new(sites);

        // `decodeSiteLoc`.
        let contig_offset = d.decode_int(Bcf2Type::Int32)?;
        let contig = self
            .contig_names
            .get(usize::try_from(contig_offset).unwrap_or(usize::MAX))
            .cloned()
            .ok_or_else(|| {
                Bcf2Error::index_out_of_bounds(contig_offset as i64, self.contig_names.len())
            })?;
        self.pos = d.decode_int(Bcf2Type::Int32)?.wrapping_add(1);
        let ref_length = d.decode_int(Bcf2Type::Int32)?;
        let start = self.pos as i64;
        let stop = self.pos as i64 + ref_length as i64 - 1;

        // `decodeSitesExtendedInfo`.
        let log10_p_error = match d.decode_single_value(Bcf2Type::Float)? {
            Typed::Double(q) => q / -10.0,
            _ => NO_LOG10_PERROR,
        };
        let n_allele_info = d.decode_int(Bcf2Type::Int32)?;
        let n_format_samples = d.decode_int(Bcf2Type::Int32)?;
        let n_alleles = n_allele_info >> 16;
        let n_info = n_allele_info & 0x0000_FFFF;
        let n_format_fields = n_format_samples >> 24;
        let n_samples = n_format_samples & 0x000F_FFFF;
        if self.header.samples.len() as i32 != n_samples {
            return Err(self.error(format!(
                "Reading BCF2 files with different numbers of samples per record is not \
                 currently supported.  Saw {} samples in header but have a record with \
                 {n_samples} samples",
                self.header.samples.len()
            )));
        }

        let id = match d.decode_typed_value()? {
            Typed::Null => ".".to_string(),
            Typed::Str(s) => s,
            other => return Err(Bcf2Error::class_cast(other.class(), "java.lang.String")),
        };

        let mut alleles = Vec::with_capacity(n_alleles.max(0) as usize);
        for i in 0..n_alleles {
            let bases = match d.decode_typed_value()? {
                Typed::Str(s) => s,
                Typed::Null => {
                    return Err(Bcf2Error::null_pointer(
                        "Cannot invoke \"String.getBytes()\" because \"bases\" is null",
                    ))
                }
                other => return Err(Bcf2Error::class_cast(other.class(), "java.lang.String")),
            };
            let allele = Allele::from_str(&bases, i == 0)
                .map_err(|e| Bcf2Error::illegal_argument(e.to_string()))?;
            alleles.push(allele);
        }

        // `decodeFilter`.
        let filters = match d.decode_typed_value()? {
            Typed::Null => None,
            Typed::Int(offset) => {
                let f = self.dictionary_string(offset)?;
                if f == "PASS" {
                    Some(Vec::new())
                } else {
                    Some(vec![f.to_string()])
                }
            }
            Typed::List(offsets) => {
                let mut filters: Vec<String> = Vec::new();
                for offset in offsets {
                    let Typed::Int(offset) = offset else {
                        return Err(Bcf2Error::class_cast(offset.class(), "java.lang.Integer"));
                    };
                    let f = self.dictionary_string(offset)?.to_string();
                    if !filters.contains(&f) {
                        filters.push(f);
                    }
                }
                Some(filters)
            }
            other => return Err(Bcf2Error::class_cast(other.class(), "java.lang.Integer")),
        };

        // `decodeInfo`: a `HashMap(nInfo)` and then the builder's copy.
        let mut attributes: Vec<(String, Value)> = Vec::new();
        for _ in 0..n_info {
            let key = self.read_dictionary_string(&mut d)?;
            let mut value = d.decode_typed_value()?.into_value();
            let meta = meta_data_for_field(&self.header, &self.fields, &key)?;
            if meta.line_type == LineType::Flag {
                value = Value::Bool(true);
            }
            match attributes.iter_mut().find(|(k, _)| *k == key) {
                Some(slot) => slot.1 = value,
                None => attributes.push((key, value)),
            }
        }
        let attributes = if attributes.is_empty() {
            attributes
        } else {
            let n = n_info.max(0) as usize;
            super::decode::reorder(&attributes, &[Fill::PutSized(n), Fill::PutAll])
        };

        // `SitesInfoForDecoding.isValid`.
        if n_format_fields < 0 || n_samples < 0 || alleles.is_empty() || !alleles[0].is_reference()
        {
            let shown: Vec<String> = alleles
                .iter()
                .map(|a| {
                    format!(
                        "{}{}",
                        a.display_string(),
                        if a.is_reference() { "*" } else { "" }
                    )
                })
                .collect();
            return Err(self.error(format!(
                "Sites info is malformed: nFormatFields = {n_format_fields}, nSamples = \
                 {n_samples}, alleles = [{}]",
                shown.join(", ")
            )));
        }

        let vc = VariantContext {
            contig,
            start,
            stop,
            id,
            alleles,
            log10_p_error,
            filters,
            attributes,
            genotypes: GenotypesContext::default(),
        };
        validate(&vc)?;
        Ok((vc, n_format_fields))
    }

    /// `createLazyGenotypesDecoder`, on the genotype block read after the sites were decoded.
    fn attach_genotypes(
        &self,
        mut vc: VariantContext,
        n_format_fields: i32,
        genotypes: &[u8],
    ) -> Result<VariantContext, Bcf2Error> {
        if !self.header.samples.is_empty() {
            let decoded = self.decode_genotypes(&vc, n_format_fields, genotypes)?;
            vc.genotypes = if self.samples_sorted {
                GenotypesContext::lazy_bcf(
                    decoded,
                    LazyData {
                        header: self.lazy_header.clone(),
                        n_genotype_fields: n_format_fields,
                        bytes: Arc::from(genotypes),
                    },
                )
            } else {
                GenotypesContext::new(decoded)
            };
        }
        Ok(vc)
    }

    /// `BCF2LazyGenotypesDecoder.parse`.
    fn decode_genotypes(
        &self,
        vc: &VariantContext,
        n_fields: i32,
        bytes: &[u8],
    ) -> Result<Vec<Genotype>, Bcf2Error> {
        let mut d = Decoder::new(bytes);
        let mut gbs: Vec<Genotype> = self
            .header
            .samples
            .iter()
            .map(|s| Genotype::new(s, Vec::new()))
            .collect();
        for _ in 0..n_fields {
            let field = self.read_dictionary_string(&mut d)?;
            let descriptor = d.read_type_descriptor();
            let n = d.decode_number_of_elements(descriptor)?;
            match field.as_str() {
                "GT" => decode_gt(&vc.alleles, &mut d, descriptor, n, &mut gbs)?,
                "FT" => {
                    for gb in gbs.iter_mut() {
                        gb.filters = match d.decode_typed_value_sized(descriptor, n)? {
                            Typed::Null => None,
                            Typed::Str(s) => Some(s),
                            other => {
                                return Err(Bcf2Error::tribble(format!(
                                    "BUG: expected encoding of field {field} inconsistent with \
                                     the value observed in the decoded value ({})",
                                    other.class()
                                )))
                            }
                        };
                    }
                }
                "DP" | "GQ" => {
                    for gb in gbs.iter_mut() {
                        let v = d.decode_int_or(descriptor, -1)?;
                        let v = if v == -1 { None } else { Some(v) };
                        if field == "DP" {
                            gb.dp = v;
                        } else {
                            gb.gq = v;
                        }
                    }
                }
                "AD" | "PL" => {
                    let t = type_for(descriptor)?;
                    for gb in gbs.iter_mut() {
                        let v = d.decode_int_array(n, t)?;
                        if field == "AD" {
                            gb.ad = v;
                        } else {
                            gb.pl = v;
                        }
                    }
                }
                _ => {
                    for gb in gbs.iter_mut() {
                        let value = match d.decode_typed_value_sized(descriptor, n)? {
                            Typed::Null => continue,
                            Typed::List(mut items) if items.len() == 1 => items.remove(0),
                            v => v,
                        };
                        let value = value.into_value();
                        match gb.extended.iter_mut().find(|(k, _)| *k == field) {
                            Some(slot) => slot.1 = value,
                            None => gb.extended.push((field.clone(), value)),
                        }
                    }
                }
            }
        }
        Ok(gbs)
    }
}

/// `GTDecoder`.
fn decode_gt(
    site_alleles: &[Allele],
    d: &mut Decoder,
    descriptor: u8,
    n: i32,
    gbs: &mut [Genotype],
) -> Result<(), Bcf2Error> {
    let t = type_for(descriptor)?;
    let allele = |encoded: i32| -> Result<Allele, Bcf2Error> {
        let offset = encoded >> 1;
        if offset == 0 {
            Ok(Allele::no_call())
        } else {
            site_alleles
                .get(usize::try_from(offset - 1).unwrap_or(usize::MAX))
                .cloned()
                .ok_or_else(|| {
                    Bcf2Error::index_out_of_bounds((offset - 1) as i64, site_alleles.len())
                })
        }
    };
    if site_alleles.len() == 2 && n == 2 {
        // `fastBiallelicDiploidDecode`: the phase comes from the second value even when it is
        // missing, and a missing first value leaves the genotype with no alleles at all.
        for gb in gbs.iter_mut() {
            let a1 = d.decode_int(t)?;
            let a2 = d.decode_int(t)?;
            gb.alleles = if a1 == t.missing_bytes() {
                Vec::new()
            } else if a2 == t.missing_bytes() {
                vec![allele(a1)?]
            } else {
                vec![allele(a1)?, allele(a2)?]
            };
            gb.phased = (a2 & 0x01) == 1;
        }
    } else {
        for gb in gbs.iter_mut() {
            match d.decode_int_array(n, t)? {
                None => gb.alleles = Vec::new(),
                Some(encoded) => {
                    gb.alleles = encoded
                        .iter()
                        .map(|&e| allele(e))
                        .collect::<Result<_, _>>()?;
                    let phase_from = if encoded.len() > 1 {
                        encoded[1]
                    } else {
                        encoded[0]
                    };
                    gb.phased = (phase_from & 0x01) == 1;
                }
            }
        }
    }
    Ok(())
}

/// The checks `VariantContext`'s constructor makes on a decoded record.
fn validate(vc: &VariantContext) -> Result<(), Bcf2Error> {
    let shown = || {
        let alleles: Vec<String> = vc
            .alleles
            .iter()
            .map(|a| {
                format!(
                    "{}{}",
                    a.display_string(),
                    if a.is_reference() { "*" } else { "" }
                )
            })
            .collect();
        format!("[{}]", alleles.join(", "))
    };
    for (i, a) in vc.alleles.iter().enumerate() {
        if vc.alleles[..i].contains(a) {
            return Err(Bcf2Error::illegal_argument(format!(
                "Duplicate allele added to VariantContext: {}{}",
                a.display_string(),
                if a.is_reference() { "*" } else { "" }
            )));
        }
    }
    if vc.alleles.iter().filter(|a| a.is_reference()).count() > 1 {
        return Err(Bcf2Error::illegal_argument(format!(
            "BUG: Received two reference tagged alleles in VariantContext {}",
            shown()
        )));
    }
    if vc.log10_p_error > 0.0 && vc.log10_p_error != NO_LOG10_PERROR {
        return Err(Bcf2Error::illegal_argument(format!(
            "BUG: log10PError cannot be > 0 : {}",
            vc.log10_p_error
        )));
    }
    if vc.log10_p_error.is_infinite() {
        return Err(Bcf2Error::illegal_argument(
            "BUG: log10PError should not be Infinity",
        ));
    }
    if vc.log10_p_error.is_nan() {
        return Err(Bcf2Error::illegal_argument(
            "BUG: log10PError should not be NaN",
        ));
    }
    Ok(())
}

/// `readHeader`.
fn read_header(bytes: &[u8]) -> Result<(Codec, usize), Bcf2Error> {
    let bare =
        |message: &str| Bcf2Error::tribble(format!("{message}, at record 0 with position 0:"));
    if bytes.len() < 3 || &bytes[..3] != b"BCF" {
        return Err(bare(
            "Input stream does not contain a BCF encoded file; BCF magic header info not found",
        ));
    }
    let major = bytes.get(3).map_or(-1, |&b| b as i32);
    let minor = bytes.get(4).map_or(-1, |&b| b as i32);
    if major != 2 {
        return Err(bare(&format!(
            "BCF2Codec can only process BCF2 files, this file has major version {major}"
        )));
    }
    if minor != 1 {
        return Err(bare(&format!(
            "BCF2Codec can only process BCF2 files with minor version = 1 but this file has minor \
             version {minor}"
        )));
    }
    let mut cursor = Cursor::new(bytes.get(5..).unwrap_or(&[]));
    let size = Bcf2Type::Int32.read(&mut cursor)?;
    if size <= 0 || size > MAX_HEADER_SIZE {
        return Err(bare(&format!(
            "BCF2 header has invalid length: {size} must be >= 0 and < {MAX_HEADER_SIZE}"
        )));
    }
    let start = 9;
    let end = start + size as usize;
    if end > bytes.len() {
        return Err(bare(&format!(
            "Couldn't read all of the bytes specified in the header length = {size}"
        )));
    }
    let text_bytes = &bytes[start..end];
    let cut = text_bytes
        .iter()
        .position(|&b| b == 0)
        .unwrap_or(text_bytes.len());
    let text = String::from_utf8_lossy(&text_bytes[..cut]);
    let file = crate::reader::read_vcf(&text)
        .map_err(|failure| Bcf2Error::new(failure.error.class(), failure.error.message()))?;
    let (header, codec_version, header_version) =
        (file.header, file.codec_version, file.header_version);

    let mut contig_names = Vec::new();
    for line in &header.lines {
        if let HeaderLine::Contig { fields, .. } = line {
            let id = fields
                .iter()
                .find(|(k, _)| k == "ID")
                .map(|(_, v)| v.clone());
            match id {
                Some(id) if !id.is_empty() => contig_names.push(id),
                _ => {
                    return Err(bare(&format!(
                        "found a contig with an invalid ID {}",
                        line.render()
                    )))
                }
            }
        }
    }
    if contig_names.is_empty() {
        return Err(bare("Didn't find any contig lines in BCF2 file header"));
    }
    let dictionary = make_dictionary(&header.lines);
    let samples_sorted = header
        .samples
        .windows(2)
        .all(|w| super::java_compare(&w[0], &w[1]) != std::cmp::Ordering::Greater);
    let lazy_header = Arc::new(LazyHeader {
        samples_in_order: sample_names_in_order(&header.samples),
        id_lines: id_lines(&header.lines),
    });
    Ok((
        Codec {
            fields: FieldIndex::new(&header),
            header,
            codec_version,
            header_version,
            contig_names,
            dictionary,
            lazy_header,
            samples_sorted,
            record_no: 0,
            pos: 0,
        },
        end,
    ))
}

/// The header of a BCF file alone, as `VCFFileReader.getFileHeader` reads it.
pub fn read_bcf_header(bytes: &[u8]) -> Result<VcfHeader, Bcf2Error> {
    read_header(bytes).map(|(codec, _)| codec.header)
}

/// Read a whole BCF file, uncompressed, as `VCFFileReader` iterates it.
pub fn read_bcf(bytes: &[u8]) -> Result<Bcf2File, Bcf2ReadFailure> {
    let (mut codec, mut at) = read_header(bytes).map_err(|error| Bcf2ReadFailure {
        error,
        records: Vec::new(),
        in_header: true,
    })?;
    let mut records = Vec::new();
    while at < bytes.len() {
        codec.record_no += 1;
        let result = (|| {
            let mut cursor = Cursor::new(&bytes[at..]);
            let sites_size = Bcf2Type::Int32.read(&mut cursor)?;
            let genotypes_size = Bcf2Type::Int32.read(&mut cursor)?;
            let block = |offset: usize, size: i32| -> Result<&[u8], Bcf2Error> {
                if size < 0 {
                    return Err(Bcf2Error::tribble(format!("Invalid block size {size}")));
                }
                let begin = (offset).min(bytes.len());
                let end = offset + size as usize;
                if end > bytes.len() {
                    return Err(Bcf2Error::tribble(format!(
                        "Failed to read next complete record: expected {size} bytes but read only \
                         {} after 0 iterations",
                        bytes.len() - begin
                    )));
                }
                Ok(&bytes[begin..end])
            };
            let sites_at = (at + 8).min(bytes.len());
            let sites = block(sites_at, sites_size)?;
            let genotypes_at = sites_at + sites_size as usize;
            // The sites block is decoded before the genotype block is read.
            let (vc, n_format_fields) = codec.decode_sites(sites)?;
            let genotypes = block(genotypes_at, genotypes_size)?;
            let vc = codec.attach_genotypes(vc, n_format_fields, genotypes)?;
            Ok((vc, genotypes_at + genotypes_size as usize))
        })();
        match result {
            Ok((vc, next)) => {
                records.push(vc);
                at = next;
            }
            Err(error) => {
                return Err(Bcf2ReadFailure {
                    error,
                    records,
                    in_header: false,
                })
            }
        }
    }
    Ok(Bcf2File {
        header: codec.header,
        codec_version: codec.codec_version,
        header_version: codec.header_version,
        records,
    })
}
