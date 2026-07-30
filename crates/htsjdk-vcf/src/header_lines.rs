//! Turning one `##` line into a typed header line.
//!
//! Ported from `htsjdk.variant.vcf.VCFCompoundHeaderLine`, `VCFInfoHeaderLine`,
//! `VCFFormatHeaderLine`, `VCFFilterHeaderLine`, `VCFSimpleHeaderLine` and `VCFContigHeaderLine`
//! (htsjdk 4.2.0), on top of the scanner in [`crate::header_parse`].
//!
//! The scanner says which pairs a line carries. This says what they mean, and it is where the
//! refusals live. There are **four** different failures, and which one a file gets is a property of
//! the field that is wrong:
//!
//! | wrong field | failure |
//! |---|---|
//! | `Number=x` | `java.lang.NumberFormatException`, uncaught: nothing wraps `Integer.parseInt` |
//! | `Number=-1` | `TribbleException$InvalidHeader`, "Count < 0 for fixed size VCF header field" |
//! | `Number=0` on a non-`Flag` | `java.lang.IllegalArgumentException` from `validate()` |
//! | `Type=integer` | plain `TribbleException`, "not a valid type ... types are case-sensitive" |
//!
//! Two of those are unchecked Java exceptions that no `catch` in the codec touches, so a malformed
//! `Number` does not produce "malformed header": it produces a `NumberFormatException` carrying the
//! offending string. A port that funnelled every failure into one error type would be answering a
//! different question from the reference.
//!
//! # `Flag` is not symmetric
//!
//! `INFO` allows `Type=Flag` and `FORMAT` does not, so the same line is valid under one key and an
//! `IllegalArgumentException` under the other. And a `Flag` with a non-zero count is **silently
//! rewritten to 0** rather than refused, which is a value change rather than an error.
//!
//! # `Source` and `Version` are version-gated
//!
//! They are read only for VCF 4.2 and later. Under 4.0 or 4.1 they are not recommended tags either,
//! so the tag-order validation rejects them outright: the same line is a valid 4.2 header line and
//! an invalid 4.1 one.

use crate::header::{Cardinality, HeaderLine, LineType};
use crate::header_parse::{parse_structured_value, InvalidHeader, VcfVersion};

/// The four ways a line fails, kept apart because upstream throws four different classes and a
/// caller can tell them apart.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HeaderLineError {
    /// `TribbleException.InvalidHeader`, which prefixes its message.
    InvalidHeader(String),
    /// A plain `TribbleException`, which does not.
    Tribble(String),
    /// `java.lang.IllegalArgumentException`.
    IllegalArgument(String),
    /// `java.lang.NumberFormatException`, from an unguarded `Integer.parseInt`.
    NumberFormat(String),
    /// `java.lang.NullPointerException`, from an unguarded `mapping.get("Number").equals(...)`.
    /// The message is the JVM's helpful-NPE text, which is why it is measured rather than written.
    NullPointer(String),
}

impl HeaderLineError {
    /// The Java class name, as a dump reports it.
    pub fn class(&self) -> &'static str {
        match self {
            HeaderLineError::InvalidHeader(_) => "htsjdk.tribble.TribbleException$InvalidHeader",
            HeaderLineError::Tribble(_) => "htsjdk.tribble.TribbleException",
            HeaderLineError::IllegalArgument(_) => "java.lang.IllegalArgumentException",
            HeaderLineError::NumberFormat(_) => "java.lang.NumberFormatException",
            HeaderLineError::NullPointer(_) => "java.lang.NullPointerException",
        }
    }

    /// What `getMessage()` returns.
    pub fn message(&self) -> String {
        match self {
            HeaderLineError::InvalidHeader(reason) => InvalidHeader(reason.clone()).message(),
            HeaderLineError::Tribble(message)
            | HeaderLineError::IllegalArgument(message)
            | HeaderLineError::NumberFormat(message)
            | HeaderLineError::NullPointer(message) => message.clone(),
        }
    }
}

impl From<InvalidHeader> for HeaderLineError {
    fn from(error: InvalidHeader) -> Self {
        HeaderLineError::InvalidHeader(error.0)
    }
}

impl std::fmt::Display for HeaderLineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.class(), self.message())
    }
}

impl std::error::Error for HeaderLineError {}

impl VcfVersion {
    /// `VCFHeaderVersion.isAtLeastAsRecentAs`, which compares **ordinals**, so the order the enum
    /// is declared in is the order of the versions.
    pub fn is_at_least(self, other: VcfVersion) -> bool {
        self.ordinal() >= other.ordinal()
    }

    fn ordinal(self) -> u8 {
        match self {
            VcfVersion::Vcf3_2 => 0,
            VcfVersion::Vcf3_3 => 1,
            VcfVersion::Vcf4_0 => 2,
            VcfVersion::Vcf4_1 => 3,
            VcfVersion::Vcf4_2 => 4,
            VcfVersion::Vcf4_3 => 5,
            VcfVersion::Vcf4_4 => 6,
        }
    }
}

/// `VCFHeaderLine.UNBOUND_DESCRIPTION`, substituted when a compound line carries no `Description`.
/// `ALLOW_UNBOUND_DESCRIPTIONS` is `true`, so a missing description is a default rather than a
/// refusal.
pub const UNBOUND_DESCRIPTION: &str = "Not provided in original VCF header";

/// The JVM's helpful-NullPointerException text for the unguarded `Number` dereference. Measured
/// from the pinned oracle rather than written from the source: the wording is the JVM's, not
/// htsjdk's, so it is an observable of the pinned image.
pub const NUMBER_IS_NULL: &str = "MEASURE-ME";

/// `AbstractVCFCodec.parseHeaderFromLines`, for one line: which prefix decides which type.
///
/// `contig_index` is the running counter the codec keeps, incremented for every `##contig` line and
/// for no other, so it counts contigs rather than lines.
pub fn parse_meta_line(
    line: &str,
    version: VcfVersion,
    contig_index: i32,
) -> Result<HeaderLine, HeaderLineError> {
    if let Some(value) = line.strip_prefix("##INFO=") {
        return compound("INFO", value, version, true);
    }
    if let Some(value) = line.strip_prefix("##FORMAT=") {
        return compound("FORMAT", value, version, false);
    }
    if let Some(value) = line.strip_prefix("##FILTER=") {
        return filter(value);
    }
    if let Some(value) = line.strip_prefix("##contig=") {
        return contig(value, contig_index);
    }
    // The fallback: a `##key=value` line, kept verbatim. A line with no `=` at all is dropped
    // upstream rather than kept, which is why this returns an error the caller reads as "skip".
    let rest = &line[2..];
    match rest.find('=') {
        Some(equals) => Ok(HeaderLine::Unstructured {
            key: rest[..equals].to_string(),
            value: rest[equals + 1..].to_string(),
        }),
        None => Err(HeaderLineError::IllegalArgument(format!(
            "no '=' in header line {line}"
        ))),
    }
}

/// `VCFCompoundHeaderLine(String line, VCFHeaderVersion version, SupportedHeaderLineType lineType)`.
fn compound(
    key: &str,
    value: &str,
    version: VcfVersion,
    allow_flag: bool,
) -> Result<HeaderLine, HeaderLineError> {
    let expected = ["ID", "Number", "Type", "Description"];
    // Recommended tags exist only from 4.2, so under 4.1 a `Source` tag is simply unexpected.
    let recommended: &[&str] = if version.is_at_least(VcfVersion::Vcf4_2) {
        &["Source", "Version"]
    } else {
        &[]
    };
    let mapping = parse_structured_value(value, Some(&expected), recommended)?;
    let get = |name: &str| {
        mapping
            .iter()
            .find(|(k, _)| k == name)
            .map(|(_, v)| v.clone())
    };

    let name = get("ID");
    // `mapping.get("Number")` is dereferenced without a null check upstream, so a line with no
    // `Number` throws a NullPointerException. The tag-order validation cannot prevent it: it only
    // checks the tags that are present, so `<ID=A>` reaches this line.
    let Some(number) = get("Number") else {
        return Err(HeaderLineError::NullPointer(NUMBER_IS_NULL.to_string()));
    };

    let mut count: i32 = -1;
    let cardinality = match number.as_str() {
        "A" => Cardinality::A,
        "R" => Cardinality::R,
        "G" => Cardinality::G,
        // `.` for VCF 4 and up, `-1` before it.
        "." if version.is_at_least(VcfVersion::Vcf4_0) => Cardinality::Unbounded,
        "-1" if !version.is_at_least(VcfVersion::Vcf4_0) => Cardinality::Unbounded,
        other => {
            count = other.parse::<i32>().map_err(|_| {
                // Nothing catches this upstream, and the message is Java's own.
                HeaderLineError::NumberFormat(format!("For input string: \"{other}\""))
            })?;
            Cardinality::Fixed(count)
        }
    };

    if count < 0 && matches!(cardinality, Cardinality::Fixed(_)) {
        return Err(HeaderLineError::InvalidHeader(format!(
            "Count < 0 for fixed size VCF header field {}",
            name.clone().unwrap_or_default()
        )));
    }

    let type_name = get("Type");
    let line_type = match type_name.as_deref() {
        Some("Integer") => LineType::Integer,
        Some("Float") => LineType::Float,
        Some("String") => LineType::String,
        Some("Character") => LineType::Character,
        Some("Flag") => LineType::Flag,
        other => {
            // `VCFHeaderLineType.valueOf` throws, and the catch turns it into a plain
            // TribbleException without the "malformed header" prefix.
            return Err(HeaderLineError::Tribble(format!(
                "{} is not a valid type in the VCF specification (note that types are \
                 case-sensitive)",
                other.unwrap_or("null")
            )));
        }
    };

    if line_type == LineType::Flag && !allow_flag {
        return Err(HeaderLineError::IllegalArgument(format!(
            "Flag is an unsupported type for this kind of field at line - {value}"
        )));
    }

    let description = get("Description").unwrap_or_else(|| UNBOUND_DESCRIPTION.to_string());

    // `validate()`.
    let mut cardinality = cardinality;
    if line_type != LineType::Flag && matches!(cardinality, Cardinality::Fixed(n) if n <= 0) {
        return Err(HeaderLineError::IllegalArgument(format!(
            "Invalid count number, with fixed count the number should be 1 or higher: key={key} \
             name={} type={} desc={description} lineType={key} count={count}",
            name.clone().unwrap_or_default(),
            type_render(line_type)
        )));
    }
    let Some(name) = name else {
        return Err(HeaderLineError::IllegalArgument(format!(
            "Invalid VCFCompoundHeaderLine: key={key} name=null type={} desc={description} \
             lineType={key}",
            type_render(line_type)
        )));
    };
    if name.contains('<') || name.contains('>') {
        return Err(HeaderLineError::IllegalArgument(
            "VCFHeaderLine: ID cannot contain angle brackets".to_string(),
        ));
    }
    if name.contains('=') {
        return Err(HeaderLineError::IllegalArgument(
            "VCFHeaderLine: ID cannot contain an equals sign".to_string(),
        ));
    }
    // A Flag with any other count is rewritten rather than refused.
    if line_type == LineType::Flag {
        cardinality = Cardinality::Fixed(0);
    }

    // Everything after Description, in order, which is where Source and Version land when the
    // version admits them.
    let extra: Vec<(String, String)> = mapping
        .iter()
        .filter(|(k, _)| !expected.contains(&k.as_str()))
        .filter(|(k, _)| {
            version.is_at_least(VcfVersion::Vcf4_2) || (k != "Source" && k != "Version")
        })
        .cloned()
        .collect();

    Ok(HeaderLine::Compound {
        key: key.to_string(),
        id: name,
        number: cardinality,
        line_type,
        description,
        extra,
    })
}

fn type_render(line_type: LineType) -> &'static str {
    match line_type {
        LineType::Integer => "Integer",
        LineType::Float => "Float",
        LineType::String => "String",
        LineType::Character => "Character",
        LineType::Flag => "Flag",
    }
}

/// `VCFFilterHeaderLine(String line, VCFHeaderVersion version)`, whose expected tags are `ID` and
/// `Description` and which admits no recommended tags at any version.
fn filter(value: &str) -> Result<HeaderLine, HeaderLineError> {
    let mapping = parse_structured_value(value, Some(&["ID", "Description"]), &[])?;
    let fields = simple_fields("FILTER", &mapping)?;
    let id = fields
        .iter()
        .find(|(k, _)| k == "ID")
        .map(|(_, v)| v.clone())
        .expect("simple_fields refuses a line with no ID");
    let description = fields
        .iter()
        .find(|(k, _)| k == "Description")
        .map(|(_, v)| v.clone())
        .unwrap_or_default();
    Ok(HeaderLine::Filter { id, description })
}

/// `VCFContigHeaderLine(String line, VCFHeaderVersion version, String key, int contigIndex)`, whose
/// expected tag order is `null`: a contig line may carry its fields in any order.
fn contig(value: &str, contig_index: i32) -> Result<HeaderLine, HeaderLineError> {
    let mapping = parse_structured_value(value, None, &[])?;
    let fields = simple_fields("contig", &mapping)?;
    if contig_index < 0 {
        return Err(HeaderLineError::Tribble(
            "The contig index is less than zero.".to_string(),
        ));
    }
    Ok(HeaderLine::Contig {
        index: contig_index,
        fields,
    })
}

/// `VCFSimpleHeaderLine.initialize`, shared by every structured line that is not compound.
fn simple_fields(
    key: &str,
    mapping: &[(String, String)],
) -> Result<Vec<(String, String)>, HeaderLineError> {
    let name = mapping.iter().find(|(k, _)| k == "ID").map(|(_, v)| v);
    // `name == null || genericFields.isEmpty()`, and the message names the key rather than the
    // line, so two different malformed lines produce the same text.
    if name.is_none() || mapping.is_empty() {
        return Err(HeaderLineError::IllegalArgument(format!(
            "Invalid VCFSimpleHeaderLine: key={key} name={}",
            name.cloned().unwrap_or_else(|| "null".to_string())
        )));
    }
    let name = name.expect("checked above");
    if name.contains('<') || name.contains('>') {
        return Err(HeaderLineError::IllegalArgument(
            "VCFHeaderLine: ID cannot contain angle brackets".to_string(),
        ));
    }
    if name.contains('=') {
        return Err(HeaderLineError::IllegalArgument(
            "VCFHeaderLine: ID cannot contain an equals sign".to_string(),
        ));
    }
    Ok(mapping.to_vec())
}
