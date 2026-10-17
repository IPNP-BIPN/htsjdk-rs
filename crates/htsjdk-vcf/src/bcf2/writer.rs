//! `BCF2Writer`, with the field encoders and writers it builds from the header.
//!
//! The writer does not own a stream. [`Bcf2Writer::header_bytes`] is what `writeHeader` writes
//! and [`Bcf2Writer::encode`] is what one `add` writes, so a caller that indexes on the fly knows
//! each record's offset before its bytes, which is where `IndexingVariantContextWriter.add` takes
//! it.

use std::collections::HashMap;

use super::decode::{count_for, fully_decode_indexed, max_ploidy, Compound, FieldIndex};
use super::encoder::Bcf2Encoder;
use super::java_hash::{order_through, Fill};
use super::types::{
    determine_integer_type, determine_integer_type_of_array, determine_integer_type_of_list,
    max_integer_type, Bcf2Type, MAX_ALLELES_IN_GENOTYPES,
};
use super::{id_lines, make_dictionary, sample_names_in_order, Bcf2Error, LazyData, LazyHeader};
use crate::allele::Allele;
use crate::header::{Cardinality, HeaderLine, LineType, VcfHeader};
use crate::variant::{Genotype, Value, VariantContext};

/// `BCF2Writer.MAJOR_VERSION` and `MINOR_VERSION`.
pub const MAJOR_VERSION: u8 = 2;
pub const MINOR_VERSION: u8 = 1;

/// `VCFConstants.MAX_GENOTYPE_QUAL`, the cap `GQAccessor` applies.
const MAX_GENOTYPE_QUAL: i32 = 99;

/// Where a record's attributes came from, which decides what `add` does before it encodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Origin {
    /// The VCF codec: not fully decoded, attributes in file order, filters in file order. `add`
    /// runs `fullyDecode` and reproduces the hash orders of the codec's maps.
    VcfCodec,
    /// The BCF2 codec: fully decoded, attributes and filters already in the order the decoded
    /// record iterates them.
    Bcf2Codec,
}

/// The value a field encoder is handed: an attribute, an accessor's `int[]`, or a filter string.
#[derive(Debug, Clone, Copy)]
enum Obj<'a> {
    Null,
    V(&'a Value),
    Ints(&'a [i32]),
    Text(&'a str),
}

impl<'a> Obj<'a> {
    fn of(value: Option<&'a Value>) -> Self {
        match value {
            None | Some(Value::Missing) => Obj::Null,
            Some(v) => Obj::V(v),
        }
    }
}

fn class_of(v: &Value) -> &'static str {
    match v {
        Value::Missing => "null",
        Value::Int(_) => "java.lang.Integer",
        Value::Double(_) => "java.lang.Double",
        Value::Bool(_) => "java.lang.Boolean",
        Value::Str(_) => "java.lang.String",
        Value::List(_) => "java.util.ArrayList",
    }
}

/// `BCF2Utils.toList(Integer.class, o)`, read element by element as `Integer`.
fn int_list(o: Obj) -> Result<Vec<Option<i32>>, Bcf2Error> {
    let one = |v: &Value| match v {
        Value::Missing => Ok(None),
        Value::Int(i) => Ok(Some(*i as i32)),
        other => Err(Bcf2Error::class_cast(class_of(other), "java.lang.Integer")),
    };
    match o {
        Obj::Null => Ok(Vec::new()),
        Obj::Ints(ints) => Ok(ints.iter().map(|&i| Some(i)).collect()),
        Obj::Text(_) => Err(Bcf2Error::class_cast(
            "java.lang.String",
            "java.lang.Integer",
        )),
        Obj::V(Value::List(items)) => items.iter().map(one).collect(),
        Obj::V(v) => Ok(vec![one(v)?]),
    }
}

/// `BCF2Utils.toList(Double.class, o)`.
fn double_list(o: Obj) -> Result<Vec<Option<f64>>, Bcf2Error> {
    let one = |v: &Value| match v {
        Value::Missing => Ok(None),
        Value::Double(d) => Ok(Some(*d)),
        other => Err(Bcf2Error::class_cast(class_of(other), "java.lang.Double")),
    };
    match o {
        Obj::Null => Ok(Vec::new()),
        Obj::Ints(_) => Err(Bcf2Error::class_cast("[I", "java.lang.Double")),
        Obj::Text(_) => Err(Bcf2Error::class_cast(
            "java.lang.String",
            "java.lang.Double",
        )),
        Obj::V(Value::List(items)) => items.iter().map(one).collect(),
        Obj::V(v) => Ok(vec![one(v)?]),
    }
}

/// `StringOrCharacter.javaStringToBCF2String`.
fn bcf2_string(o: Obj) -> Result<String, Bcf2Error> {
    match o {
        Obj::Null => Ok(String::new()),
        Obj::Text(s) => Ok(s.to_string()),
        Obj::Ints(_) => Err(Bcf2Error::class_cast("[I", "java.lang.String")),
        Obj::V(Value::Str(s)) => Ok(s.clone()),
        Obj::V(Value::List(items)) => {
            let strings = items
                .iter()
                .map(|v| match v {
                    Value::Missing => Ok(None),
                    Value::Str(s) => Ok(Some(s.as_str())),
                    other => Err(Bcf2Error::class_cast(class_of(other), "java.lang.String")),
                })
                .collect::<Result<Vec<_>, _>>()?;
            collapse_string_list(&strings)
        }
        Obj::V(other) => Err(Bcf2Error::class_cast(class_of(other), "java.lang.String")),
    }
}

/// `BCF2Utils.collapseStringList`: one string as itself, more behind a leading comma.
fn collapse_string_list(strings: &[Option<&str>]) -> Result<String, Bcf2Error> {
    match strings {
        [] => Ok(String::new()),
        [one] => one.map(str::to_string).ok_or_else(|| {
            Bcf2Error::null_pointer("Cannot invoke \"String.length()\" because \"s\" is null")
        }),
        many => {
            let mut b = String::new();
            for s in many.iter().flatten() {
                b.push(',');
                b.push_str(s);
            }
            Ok(b)
        }
    }
}

/// `String.length()`, in UTF-16 units.
fn java_length(s: &str) -> usize {
    s.encode_utf16().count()
}

/// The subclasses of `BCF2FieldEncoder`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EncoderKind {
    StringOrCharacter,
    Flag,
    Float { atomic: bool },
    IntArray,
    AtomicInt,
    GenericInts,
}

/// `BCF2FieldEncoder`.
#[derive(Debug, Clone)]
struct FieldEncoder {
    field: String,
    number: Cardinality,
    kind: EncoderKind,
    dictionary_offset: i32,
    dictionary_offset_type: Bcf2Type,
}

impl FieldEncoder {
    fn new(
        line: &Compound,
        dict: &HashMap<String, i32>,
        kind: EncoderKind,
    ) -> Result<Self, Bcf2Error> {
        let offset = *dict.get(line.id).ok_or_else(|| {
            Bcf2Error::illegal_state(format!(
                "Format error: could not find string {} in header as required by BCF",
                line.id
            ))
        })?;
        if kind == EncoderKind::Flag && line.number != Cardinality::Fixed(0) {
            return Err(Bcf2Error::illegal_state(format!(
                "Flag encoder only supports atomic flags for field {}",
                line.id
            )));
        }
        Ok(Self {
            field: line.id.to_string(),
            number: line.number,
            kind,
            dictionary_offset: offset,
            dictionary_offset_type: determine_integer_type(offset)?,
        })
    }

    fn static_type(&self) -> Option<Bcf2Type> {
        match self.kind {
            EncoderKind::StringOrCharacter => Some(Bcf2Type::Char),
            EncoderKind::Flag => Some(Bcf2Type::Int8),
            EncoderKind::Float { .. } => Some(Bcf2Type::Float),
            _ => None,
        }
    }

    fn write_field_key(&self, encoder: &mut Bcf2Encoder) -> Result<(), Bcf2Error> {
        encoder.encode_typed_int_as(self.dictionary_offset, self.dictionary_offset_type)
    }

    fn has_constant_num_elements(&self) -> bool {
        self.kind != EncoderKind::StringOrCharacter && matches!(self.number, Cardinality::Fixed(_))
    }

    fn has_value_determined_num_elements(&self) -> bool {
        self.kind == EncoderKind::StringOrCharacter || self.number == Cardinality::Unbounded
    }

    fn has_context_determined_num_elements(&self) -> bool {
        !self.has_constant_num_elements() && !self.has_value_determined_num_elements()
    }

    /// `numElements()`: the header's count, except that a Flag writes one value for its zero.
    fn num_elements(&self) -> i32 {
        match (self.kind, self.number) {
            (EncoderKind::Flag, _) => 1,
            (_, Cardinality::Fixed(n)) => n,
            _ => -1,
        }
    }

    fn num_elements_from_value(&self, value: Obj) -> Result<i32, Bcf2Error> {
        Ok(match (self.kind, value) {
            (_, Obj::Null) => 0,
            (EncoderKind::StringOrCharacter, v) => java_length(&bcf2_string(v)?) as i32,
            (_, Obj::Ints(ints)) => ints.len() as i32,
            (EncoderKind::IntArray, Obj::V(v)) => {
                return Err(Bcf2Error::class_cast(class_of(v), "[I"))
            }
            (_, Obj::V(Value::List(items))) => items.len() as i32,
            _ => 1,
        })
    }

    /// `numElements(vc, value)`.
    fn num_elements_for(&self, vc: &VariantContext, value: Obj) -> Result<i32, Bcf2Error> {
        if self.has_constant_num_elements() {
            Ok(self.num_elements())
        } else if self.has_context_determined_num_elements() {
            Ok(count_for(self.number, vc))
        } else {
            self.num_elements_from_value(value)
        }
    }

    /// `getType(value)`: the static type, or the dynamic one an integer encoder computes.
    fn get_type(&self, value: Obj) -> Result<Option<Bcf2Type>, Bcf2Error> {
        if let Some(t) = self.static_type() {
            return Ok(Some(t));
        }
        if let Obj::Null = value {
            return Ok(Some(Bcf2Type::Int8));
        }
        let t = match self.kind {
            EncoderKind::IntArray => match value {
                Obj::Ints(ints) => determine_integer_type_of_array(ints)?,
                Obj::V(v) => return Err(Bcf2Error::class_cast(class_of(v), "[I")),
                _ => return Err(Bcf2Error::class_cast("java.lang.String", "[I")),
            },
            EncoderKind::AtomicInt => match value {
                Obj::V(Value::Int(i)) => determine_integer_type(*i as i32)?,
                Obj::V(v) => return Err(Bcf2Error::class_cast(class_of(v), "java.lang.Integer")),
                _ => return Err(Bcf2Error::class_cast("[I", "java.lang.Integer")),
            },
            _ => {
                // `for (final int value : values)` unboxes each element, and `.` decoded to null.
                let ints = int_list(value)?;
                if ints.iter().any(Option::is_none) {
                    return Err(Bcf2Error::null_pointer(
                        "Cannot invoke \"java.lang.Integer.intValue()\" because the return value \
                         of \"java.util.Iterator.next()\" is null",
                    ));
                }
                determine_integer_type_of_list(&ints.into_iter().flatten().collect::<Vec<_>>())?
            }
        };
        Ok(Some(t))
    }

    /// `encodeValue(encoder, value, type, minValues)`.
    fn encode_value(
        &self,
        encoder: &mut Bcf2Encoder,
        value: Obj,
        t: Option<Bcf2Type>,
        min_values: i32,
    ) -> Result<(), Bcf2Error> {
        let t = t.ok_or_else(|| Bcf2Error::null_pointer("type is null"))?;
        let mut count = 0;
        match self.kind {
            EncoderKind::StringOrCharacter => {
                let s = bcf2_string(value)?;
                let size = java_length(&s).max(min_values.max(0) as usize);
                encoder.encode_raw_string(s.as_bytes(), size);
                return Ok(());
            }
            EncoderKind::Flag => return encoder.encode_raw_int(1, Bcf2Type::Int8),
            EncoderKind::Float { atomic: true } => match value {
                Obj::Null => {}
                Obj::V(Value::Double(d)) => {
                    encoder.encode_raw_float(*d)?;
                    count += 1;
                }
                Obj::V(v) => return Err(Bcf2Error::class_cast(class_of(v), "java.lang.Double")),
                _ => {
                    return Err(Bcf2Error::class_cast(
                        "java.lang.String",
                        "java.lang.Double",
                    ))
                }
            },
            EncoderKind::Float { atomic: false } => {
                for d in double_list(value)?.into_iter().flatten() {
                    encoder.encode_raw_float(d)?;
                    count += 1;
                }
            }
            EncoderKind::IntArray => match value {
                Obj::Null => {}
                Obj::Ints(ints) => {
                    for &i in ints {
                        encoder.encode_raw_int(i, t)?;
                        count += 1;
                    }
                }
                Obj::V(v) => return Err(Bcf2Error::class_cast(class_of(v), "[I")),
                Obj::Text(_) => return Err(Bcf2Error::class_cast("java.lang.String", "[I")),
            },
            EncoderKind::AtomicInt => match value {
                Obj::Null => {}
                Obj::V(Value::Int(i)) => {
                    encoder.encode_raw_int(*i as i32, t)?;
                    count += 1;
                }
                Obj::V(v) => return Err(Bcf2Error::class_cast(class_of(v), "java.lang.Integer")),
                _ => return Err(Bcf2Error::class_cast("[I", "java.lang.Integer")),
            },
            EncoderKind::GenericInts => {
                for i in int_list(value)?.into_iter().flatten() {
                    encoder.encode_raw_int(i, t)?;
                    count += 1;
                }
            }
        }
        while count < min_values {
            encoder.encode_raw_missing_value(t)?;
            count += 1;
        }
        Ok(())
    }
}

/// `IntGenotypeFieldAccessors`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Accessor {
    Gq,
    Dp,
    Ad,
    Pl,
}

impl Accessor {
    fn for_field(field: &str) -> Option<Accessor> {
        match field {
            "GQ" => Some(Accessor::Gq),
            "DP" => Some(Accessor::Dp),
            "AD" => Some(Accessor::Ad),
            "PL" => Some(Accessor::Pl),
            _ => None,
        }
    }

    /// `getValues(g)`. The atomic two answer `null` for -1, which is also what an absent value
    /// reads as, after `GQAccessor`'s `Math.min(gq, 99)`.
    fn values(self, g: &Genotype) -> Option<Vec<i32>> {
        let atomic = |v: Option<i32>| match v {
            None | Some(-1) => None,
            Some(v) => Some(vec![v]),
        };
        match self {
            Accessor::Gq => atomic(g.gq.map(|q| q.min(MAX_GENOTYPE_QUAL))),
            Accessor::Dp => atomic(g.dp),
            Accessor::Ad => g.ad.clone(),
            Accessor::Pl => g.pl.clone(),
        }
    }
}

/// The subclasses of `BCF2FieldWriter.GenotypesWriter`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GenotypesKind {
    StaticallyTyped,
    IntegerTyped,
    Igf(Accessor),
    Ft,
    Gt,
}

/// `BCF2FieldWriter.GenotypesWriter`: the per-record state lives here as it does in Java.
#[derive(Debug, Clone)]
struct GenotypesWriter {
    encoder: FieldEncoder,
    kind: GenotypesKind,
    n_values_per_genotype: i32,
    encoding_type: Option<Bcf2Type>,
}

impl GenotypesWriter {
    fn new(encoder: FieldEncoder, kind: GenotypesKind) -> Self {
        let n_values_per_genotype = if encoder.has_constant_num_elements() {
            encoder.num_elements()
        } else {
            -1
        };
        let encoding_type = match kind {
            GenotypesKind::StaticallyTyped | GenotypesKind::Ft => encoder.static_type(),
            _ => None,
        };
        Self {
            encoder,
            kind,
            n_values_per_genotype,
            encoding_type,
        }
    }

    /// The value `addGenotype` and `numElements(vc, g)` read for one genotype.
    fn value<'a>(&self, g: &'a Genotype, scratch: &'a Option<Vec<i32>>) -> Obj<'a> {
        match self.kind {
            GenotypesKind::Igf(_) => match scratch {
                Some(ints) => Obj::Ints(ints),
                None => Obj::Null,
            },
            GenotypesKind::Ft => match &g.filters {
                Some(f) => Obj::Text(f),
                None => Obj::Null,
            },
            _ => Obj::of(g.get(&self.encoder.field)),
        }
    }

    fn accessor_values(&self, g: &Genotype) -> Option<Vec<i32>> {
        match self.kind {
            GenotypesKind::Igf(accessor) => accessor.values(g),
            _ => None,
        }
    }

    /// `computeMaxSizeOfGenotypeFieldFromValues`, from -1.
    fn max_size_from_values(&self, vc: &VariantContext) -> Result<i32, Bcf2Error> {
        let mut size = -1;
        for g in vc.genotypes.iter() {
            let scratch = self.accessor_values(g);
            let n = match self.kind {
                GenotypesKind::Igf(_) => scratch.as_ref().map_or(0, |v| v.len() as i32),
                _ => self.encoder.num_elements_for(vc, self.value(g, &scratch))?,
            };
            size = size.max(n);
        }
        Ok(size)
    }

    fn start(
        &mut self,
        encoder: &mut Bcf2Encoder,
        vc: &VariantContext,
        gt: &mut GtState,
    ) -> Result<(), Bcf2Error> {
        match self.kind {
            GenotypesKind::IntegerTyped => {
                let mut values = Vec::new();
                for g in vc.genotypes.iter() {
                    values.extend(
                        int_list(Obj::of(g.get(&self.encoder.field)))?
                            .into_iter()
                            .flatten(),
                    );
                }
                self.encoding_type = Some(determine_integer_type_of_list(&values)?);
            }
            GenotypesKind::Igf(accessor) => {
                let mut t = Bcf2Type::Int8;
                for g in vc.genotypes.iter() {
                    let values = accessor.values(g);
                    let value = values.as_deref().map_or(Obj::Null, Obj::Ints);
                    let values_type = self.encoder.get_type(value)?.unwrap_or(Bcf2Type::Int8);
                    t = max_integer_type(t, values_type)?;
                    if t == Bcf2Type::Int32 {
                        break;
                    }
                }
                self.encoding_type = Some(t);
            }
            GenotypesKind::Gt => {
                if vc.alleles.len() > MAX_ALLELES_IN_GENOTYPES {
                    return Err(Bcf2Error::illegal_state(format!(
                        "Current BCF2 encoder cannot handle sites with > {MAX_ALLELES_IN_GENOTYPES} \
                         alleles, but you have {} at {}:{}",
                        vc.alleles.len(),
                        vc.contig,
                        vc.start
                    )));
                }
                self.encoding_type = Some(Bcf2Type::Int8);
                gt.build(vc);
                self.n_values_per_genotype = max_ploidy(vc) as i32;
            }
            _ => {}
        }
        // `GenotypesWriter.start`.
        self.encoder.write_field_key(encoder)?;
        if !self.encoder.has_constant_num_elements() {
            self.n_values_per_genotype = if self.encoder.has_context_determined_num_elements() {
                count_for(self.encoder.number, vc)
            } else {
                self.max_size_from_values(vc)?
            };
        }
        let t = self
            .encoding_type
            .ok_or_else(|| Bcf2Error::null_pointer("encodingType is null"))?;
        encoder.encode_type(self.n_values_per_genotype, t)
    }

    fn add_genotype(
        &self,
        encoder: &mut Bcf2Encoder,
        g: &Genotype,
        gt: &GtState,
    ) -> Result<(), Bcf2Error> {
        if self.kind == GenotypesKind::Gt {
            return gt.add_genotype(encoder, g, self.n_values_per_genotype);
        }
        let scratch = self.accessor_values(g);
        self.encoder.encode_value(
            encoder,
            self.value(g, &scratch),
            self.encoding_type,
            self.n_values_per_genotype,
        )
    }
}

/// `GTWriter`'s allele map.
#[derive(Debug, Default)]
struct GtState {
    alleles: Vec<Allele>,
}

impl GtState {
    fn build(&mut self, vc: &VariantContext) {
        self.alleles = vc.alleles.clone();
    }

    /// `getAlleleOffset`: the reference, the first alternate, the no-call, then the map of the
    /// others, in which a repeated allele keeps its last index.
    fn offset(&self, a: &Allele) -> Result<i32, Bcf2Error> {
        if self.alleles.first() == Some(a) {
            return Ok(0);
        }
        if self.alleles.get(1) == Some(a) {
            return Ok(1);
        }
        if a.is_no_call() {
            return Ok(-1);
        }
        self.alleles
            .iter()
            .enumerate()
            .skip(2)
            .rev()
            .find(|(_, b)| *b == a)
            .map(|(i, _)| i as i32)
            .ok_or_else(|| {
                Bcf2Error::illegal_state(format!(
                    "BUG: Couldn't find allele offset for allele {}{}",
                    a.display_string(),
                    if a.is_reference() { "*" } else { "" }
                ))
            })
    }

    fn add_genotype(
        &self,
        encoder: &mut Bcf2Encoder,
        g: &Genotype,
        n_values: i32,
    ) -> Result<(), Bcf2Error> {
        for i in 0..n_values.max(0) as usize {
            match g.alleles.get(i) {
                Some(a) => {
                    let offset = self.offset(a)?;
                    let phase = if g.phased && i != 0 { 1 } else { 0 };
                    encoder.encode_raw_int(((offset + 1) << 1) | phase, Bcf2Type::Int8)?;
                }
                None => encoder.encode_raw_missing_value(Bcf2Type::Int8)?,
            }
        }
        Ok(())
    }
}

/// `BCF2Writer`, after `setHeader`.
#[derive(Debug)]
pub struct Bcf2Writer {
    /// `new VCFHeader(header.getMetaDataInSortedOrder(), samples)`.
    header: VcfHeader,
    do_not_write_genotypes: bool,
    contig_dictionary: HashMap<String, i32>,
    string_dictionary: HashMap<String, i32>,
    site_writers: HashMap<String, FieldEncoder>,
    genotypes_writers: HashMap<String, GenotypesWriter>,
    gt: GtState,
    encoder: Bcf2Encoder,
    /// The header of the last lazy block seen, and whether its blocks can be copied.
    last_lazy: Option<(std::sync::Arc<LazyHeader>, bool)>,
    own_lazy_header: LazyHeader,
    fields: FieldIndex,
}

impl Bcf2Writer {
    /// `setHeader`.
    pub fn new(header: &VcfHeader, do_not_write_genotypes: bool) -> Result<Self, Bcf2Error> {
        let sorted = VcfHeader {
            lines: header.sorted_lines().into_iter().cloned().collect(),
            samples: if do_not_write_genotypes {
                Vec::new()
            } else {
                header.samples.clone()
            },
        };

        let contigs: Vec<&HeaderLine> = sorted
            .lines
            .iter()
            .filter(|l| matches!(l, HeaderLine::Contig { .. }))
            .collect();
        if contigs.is_empty() {
            return Err(Bcf2Error::illegal_state(
                "Cannot write BCF2 file with missing contig lines",
            ));
        }
        let mut contig_dictionary = HashMap::new();
        for (offset, line) in contigs.iter().enumerate() {
            if let HeaderLine::Contig { fields, .. } = line {
                let id = fields
                    .iter()
                    .find(|(k, _)| k == "ID")
                    .map(|(_, v)| v.clone())
                    .unwrap_or_default();
                contig_dictionary.insert(id, offset as i32);
            }
        }

        let dict = make_dictionary(&sorted.lines);
        let string_dictionary: HashMap<String, i32> = dict
            .iter()
            .enumerate()
            .map(|(i, s)| (s.clone(), i as i32))
            .collect();

        // `BCF2FieldWriterManager.setup`, over the IDs `VCFHeader` keeps, the first line of each.
        let fields = FieldIndex::new(&sorted);
        let mut site_writers = HashMap::new();
        let mut genotypes_writers = HashMap::new();
        for line in &sorted.lines {
            let HeaderLine::Compound { key, id, .. } = line else {
                continue;
            };
            let Some(c) = fields.compound(&sorted, key, id) else {
                continue;
            };
            if !std::ptr::eq(c.line, line) {
                continue;
            }
            if key == "INFO" {
                let kind = basic_kind(&c)?;
                site_writers.insert(id.clone(), FieldEncoder::new(&c, &string_dictionary, kind)?);
            } else if key == "FORMAT" {
                let accessor = Accessor::for_field(id);
                let kind = if accessor.is_some() {
                    EncoderKind::IntArray
                } else if id == "GT" {
                    EncoderKind::GenericInts
                } else {
                    basic_kind(&c)?
                };
                let encoder = FieldEncoder::new(&c, &string_dictionary, kind)?;
                let writer_kind = if id == "GT" {
                    GenotypesKind::Gt
                } else if id == "FT" {
                    GenotypesKind::Ft
                } else if let Some(accessor) = accessor {
                    GenotypesKind::Igf(accessor)
                } else if c.line_type == LineType::Integer {
                    GenotypesKind::IntegerTyped
                } else {
                    GenotypesKind::StaticallyTyped
                };
                genotypes_writers.insert(id.clone(), GenotypesWriter::new(encoder, writer_kind));
            }
        }

        let own_lazy_header = LazyHeader {
            samples_in_order: sample_names_in_order(&sorted.samples),
            id_lines: id_lines(&sorted.lines),
        };
        Ok(Self {
            header: sorted,
            do_not_write_genotypes,
            contig_dictionary,
            string_dictionary,
            site_writers,
            genotypes_writers,
            gt: GtState::default(),
            encoder: Bcf2Encoder::new(),
            last_lazy: None,
            own_lazy_header,
            fields,
        })
    }

    /// The header the writer writes and encodes against.
    pub fn header(&self) -> &VcfHeader {
        &self.header
    }

    /// `writeHeader`: the magic, the INT32 length of the text and its `\0`, and the text.
    pub fn header_bytes(&self) -> Vec<u8> {
        let mut text = self.header.write().into_bytes();
        text.push(0);
        let mut out = vec![b'B', b'C', b'F', MAJOR_VERSION, MINOR_VERSION];
        out.extend_from_slice(&(text.len() as i32).to_le_bytes());
        out.extend_from_slice(&text);
        out
    }

    /// `add(vc)`: the bytes of one record, both block lengths and both blocks.
    pub fn encode(&mut self, vc: &VariantContext, origin: Origin) -> Result<Vec<u8>, Bcf2Error> {
        let stripped;
        let vc = if self.do_not_write_genotypes {
            stripped = VariantContext {
                genotypes: Default::default(),
                ..vc.clone()
            };
            &stripped
        } else {
            vc
        };
        let decoded;
        let (vc, filters) = match origin {
            Origin::VcfCodec => {
                decoded = fully_decode_indexed(vc, &self.header, &self.fields)?;
                (&decoded, vcf_codec_filters(vc))
            }
            Origin::Bcf2Codec => (vc, bcf2_codec_filters(vc)),
        };
        let sites = self.build_sites_data(vc, filters.as_deref())?;
        let genotypes = self.build_samples_data(vc)?;
        let mut out = Vec::with_capacity(8 + sites.len() + genotypes.len());
        out.extend_from_slice(&(sites.len() as i32).to_le_bytes());
        out.extend_from_slice(&(genotypes.len() as i32).to_le_bytes());
        out.extend_from_slice(&sites);
        out.extend_from_slice(&genotypes);
        Ok(out)
    }

    /// `getLazyData`, with `canSafelyWriteRawGenotypesBytes` cached per header; a block that
    /// cannot be copied is decoded, which is the side effect the Java warns about.
    fn lazy_data(&mut self, vc: &VariantContext) -> Option<LazyData> {
        let data = vc.genotypes.bcf_unparsed()?.clone();
        let consistent = match &self.last_lazy {
            Some((header, consistent)) if std::sync::Arc::ptr_eq(header, &data.header) => {
                *consistent
            }
            _ => {
                let consistent = ordered_consistently(&self.own_lazy_header, &data.header);
                self.last_lazy = Some((data.header.clone(), consistent));
                consistent
            }
        };
        if consistent {
            Some(data)
        } else {
            vc.genotypes.decode();
            None
        }
    }

    fn build_sites_data(
        &mut self,
        vc: &VariantContext,
        filters: Option<&[String]>,
    ) -> Result<Vec<u8>, Bcf2Error> {
        let contig_index = *self.contig_dictionary.get(&vc.contig).ok_or_else(|| {
            Bcf2Error::null_pointer(
                "Cannot invoke \"java.lang.Integer.intValue()\" because the return value of \
                 \"java.util.Map.get(Object)\" is null",
            )
        })?;
        let e = &mut self.encoder;
        e.encode_raw_int(contig_index, Bcf2Type::Int32)?;
        e.encode_raw_int((vc.start - 1) as i32, Bcf2Type::Int32)?;
        e.encode_raw_int((vc.stop - vc.start + 1) as i32, Bcf2Type::Int32)?;
        if vc.has_log10_p_error() {
            e.encode_raw_float(vc.phred_scaled_qual() as f32 as f64)?;
        } else {
            e.encode_raw_missing_value(Bcf2Type::Float)?;
        }

        let n_alleles = vc.alleles.len() as i32;
        let n_info = vc.attributes.len() as i32;
        let n_format_fields = match self.lazy_data(vc) {
            Some(data) => data.n_genotype_fields,
            None => vc
                .calc_vcf_genotype_keys(!self.header.samples.is_empty())
                .len() as i32,
        };
        let n_samples = self.header.samples.len() as i32;
        let e = &mut self.encoder;
        e.encode_raw_int((n_alleles << 16) | (n_info & 0x0000_FFFF), Bcf2Type::Int32)?;
        e.encode_raw_int(
            (n_format_fields << 24) | (n_samples & 0x000F_FFFF),
            Bcf2Type::Int32,
        )?;

        e.encode_typed_string(vc.id.as_bytes())?;
        for allele in &vc.alleles {
            e.encode_typed_string(allele.display_string().as_bytes())?;
        }

        // `buildFilter`.
        match filters {
            Some(f) if !f.is_empty() => self.encode_strings_by_ref(f)?,
            Some(_) => self.encode_strings_by_ref(&["PASS".to_string()])?,
            None => self.encoder.encode_typed_missing(Bcf2Type::Int8)?,
        }

        // `buildInfo`.
        for (field, value) in &vc.attributes {
            let Some(writer) = self.site_writers.get(field) else {
                return Err(unexpected_field(vc, field, "INFO"));
            };
            writer.write_field_key(&mut self.encoder)?;
            let raw = Obj::of(Some(value));
            let t = writer.get_type(raw)?;
            match raw {
                Obj::Null => self.encoder.encode_type(0, t.unwrap_or(Bcf2Type::Int8))?,
                _ => {
                    let count = writer.num_elements_for(vc, raw)?;
                    let t = t.ok_or_else(|| Bcf2Error::null_pointer("type is null"))?;
                    self.encoder.encode_type(count, t)?;
                    writer.encode_value(&mut self.encoder, raw, Some(t), count)?;
                }
            }
        }
        Ok(self.encoder.record_bytes())
    }

    fn encode_strings_by_ref(&mut self, strings: &[String]) -> Result<(), Bcf2Error> {
        let mut offsets = Vec::with_capacity(strings.len());
        for s in strings {
            let offset = *self.string_dictionary.get(s).ok_or_else(|| {
                Bcf2Error::illegal_state(format!(
                    "Format error: could not find string {s} in header as required by BCF"
                ))
            })?;
            offsets.push(offset);
        }
        let t = determine_integer_type_of_list(&offsets)?;
        self.encoder.encode_typed_ints(&offsets, t)
    }

    fn build_samples_data(&mut self, vc: &VariantContext) -> Result<Vec<u8>, Bcf2Error> {
        if let Some(data) = self.lazy_data(vc) {
            return Ok(data.bytes.to_vec());
        }
        let fields = vc.calc_vcf_genotype_keys(!self.header.samples.is_empty());
        for field in &fields {
            let Some(mut writer) = self.genotypes_writers.remove(field) else {
                return Err(unexpected_field(vc, field, "FORMAT"));
            };
            let result = (|| {
                writer.start(&mut self.encoder, vc, &mut self.gt)?;
                for name in &self.header.samples {
                    let missing;
                    let g = match vc.genotype(name) {
                        Some(g) => g,
                        None => {
                            if writer.n_values_per_genotype < 0 {
                                return Err(Bcf2Error::illegal_argument(format!(
                                    "List length = {}",
                                    writer.n_values_per_genotype
                                )));
                            }
                            missing =
                                Genotype::missing(name, writer.n_values_per_genotype as usize);
                            &missing
                        }
                    };
                    writer.add_genotype(&mut self.encoder, g, &self.gt)?;
                }
                Ok(())
            })();
            self.genotypes_writers.insert(field.clone(), writer);
            result?;
        }
        Ok(self.encoder.record_bytes())
    }
}

/// The encoder kind `createFieldEncoder` picks from a line's type alone.
fn basic_kind(c: &Compound) -> Result<EncoderKind, Bcf2Error> {
    Ok(match c.line_type {
        LineType::Character | LineType::String => EncoderKind::StringOrCharacter,
        LineType::Flag => EncoderKind::Flag,
        LineType::Float => EncoderKind::Float {
            atomic: c.number == Cardinality::Fixed(1),
        },
        LineType::Integer => {
            if c.number == Cardinality::Fixed(1) {
                EncoderKind::AtomicInt
            } else {
                EncoderKind::GenericInts
            }
        }
    })
}

fn unexpected_field(vc: &VariantContext, field: &str, field_type: &str) -> Bcf2Error {
    Bcf2Error::illegal_state(format!(
        "Found field {field} in the {field_type} fields of VariantContext at {}:{} from Unknown \
         but this hasn't been defined in the VCFHeader",
        vc.contig, vc.start
    ))
}

/// `BCF2Utils.headerLinesAreOrderedConsistently(outputHeader, genotypesBlockHeader)`: the same
/// samples, and every ID line of the block's header at the same place in the output's.
pub fn ordered_consistently(output: &LazyHeader, block: &LazyHeader) -> bool {
    output.samples_in_order == block.samples_in_order
        && block.id_lines.len() <= output.id_lines.len()
        && block
            .id_lines
            .iter()
            .zip(&output.id_lines)
            .all(|(a, b)| a == b)
}

/// The filters of a VCF codec record in the order its `HashSet` iterates them:
/// `new HashSet<>(filters)` sizes its map to `max((int) (n / .75f) + 1, 16)`.
fn vcf_codec_filters(vc: &VariantContext) -> Option<Vec<String>> {
    let filters = vc.filters.as_ref()?;
    let mut distinct: Vec<String> = Vec::new();
    for f in filters {
        if !distinct.contains(f) {
            distinct.push(f.clone());
        }
    }
    let capacity = ((filters.len() as f32 / 0.75) as usize + 1).max(16);
    let order = order_through(&distinct, &[Fill::PutSized(capacity)]);
    Some(order.into_iter().map(|i| distinct[i].clone()).collect())
}

/// The filters of a BCF codec record: a `LinkedHashSet`, in decode order.
fn bcf2_codec_filters(vc: &VariantContext) -> Option<Vec<String>> {
    let filters = vc.filters.as_ref()?;
    let mut distinct: Vec<String> = Vec::new();
    for f in filters {
        if !distinct.contains(f) {
            distinct.push(f.clone());
        }
    }
    Some(distinct)
}
