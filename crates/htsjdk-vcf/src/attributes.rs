//! The typed attribute accessors, which are the layer between the strings a VCF is parsed into and
//! the numbers an annotation does arithmetic on.
//!
//! Ported from `htsjdk.variant.variantcontext.CommonInfo` and
//! `htsjdk.variant.vcf.VCFUtils.parseVcfDouble` at htsjdk 4.2.0.
//!
//! # The declared Type does not decide what a record holds
//!
//! Measured, every `Type` in the format stores a **String**, or a list of Strings when `Number` is
//! not 1, with the single exception of `Flag`, which stores a boolean. `Integer`, `Float`,
//! `Character` and `String` are indistinguishable in a decoded record. So the header's Type does
//! not convert anything: it tells a caller which accessor to reach for, and the conversion happens
//! there, once per call, against a default the caller supplies.
//!
//! # The missing-value test is a reference comparison, and only one accessor has it
//!
//! `getAttributeAsInt` reads `x == VCFConstants.MISSING_VALUE_v4` — with `==`, on a String. It is
//! true only when the stored object *is* that constant. The codec assigns the constant for a bare
//! key and for `KEY=`, so those return the default; a value written `KEY=.` arrives as a substring
//! of the line, a different reference, and reaches `Integer.parseInt(".")` instead.
//!
//! ```text
//! I1=.      asInt  NumberFormatException      asString "."
//! I1=       asInt  -1 (the default)           asString "."
//! I1        asInt  -1 (the default)           asString "."
//! ```
//!
//! Three spellings of missing, identical in every rendering, and two outcomes: a number and an
//! exception. [`Value::Missing`] is that constant here and `Value::Str(".")` is the substring; the
//! distinction exists in this port for exactly this reason.
//!
//! **`getAttributeAsDouble` does not have the test at all.** Only null is checked, so all three
//! spellings reach `parseVcfDouble` and throw. The two accessors disagree about what missing means
//! for the same field.
//!
//! # `parseVcfDouble` is not `Double.parseDouble`
//!
//! On failure it retries against `^(?<sign>[-+]?)((?<inf>(INF|INFINITY))|(?<nan>NAN))$`, case
//! insensitive, so `inf`, `-INF`, `Infinity` and `nan` all parse. And `Double.parseDouble` itself
//! is wider than the format: it accepts surrounding whitespace, a trailing `f` or `d`, and hex
//! float literals, so `1f`, `1d`, `0x1p3` and `" 1"` are all numbers to a VCF reader.
//!
//! # A scalar read as a list is a list of one
//!
//! `getAttributeAsList` wraps a non-list in a singleton and answers an absent key with an **empty**
//! list rather than a list holding the default. So a caller cannot tell "one value" from "one value
//! because there was only one" without asking twice.

use crate::variant::Value;

/// `VCFConstants.MISSING_VALUE_v4`.
pub const MISSING_VALUE: &str = ".";

/// How an accessor fails. Both are unchecked upstream and neither is caught anywhere.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AttributeError {
    /// `Integer.parseInt` or `Double.parseDouble` on text that is not a number.
    NumberFormat(String),
    /// A scalar accessor reaching a list: `(String) x` on an `ArrayList`.
    ClassCast(String),
}

impl AttributeError {
    pub fn class(&self) -> &'static str {
        match self {
            AttributeError::NumberFormat(_) => "java.lang.NumberFormatException",
            AttributeError::ClassCast(_) => "java.lang.ClassCastException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            AttributeError::NumberFormat(message) | AttributeError::ClassCast(message) => {
                message.clone()
            }
        }
    }
}

/// The JVM's own wording, which a dump reports and a caller sees. The **stored** class is in it,
/// so a list and a flag reaching the same accessor produce two different messages.
fn class_cast(stored: &str) -> AttributeError {
    AttributeError::ClassCast(format!(
        "class {stored} cannot be cast to class java.lang.String ({stored} and java.lang.String \
         are in module java.base of loader 'bootstrap')"
    ))
}

/// The Java class a stored value reports.
fn java_class(value: &Value) -> &'static str {
    match value {
        Value::List(_) => "java.util.ArrayList",
        Value::Bool(_) => "java.lang.Boolean",
        Value::Int(_) => "java.lang.Integer",
        Value::Double(_) => "java.lang.Double",
        _ => "java.lang.String",
    }
}

fn number_format(text: &str) -> AttributeError {
    AttributeError::NumberFormat(format!("For input string: \"{text}\""))
}

/// `String.valueOf(Object)` over a stored attribute, which is what every renderer here agrees on.
fn value_of(value: &Value) -> String {
    match value {
        Value::Str(text) => text.clone(),
        Value::Missing => MISSING_VALUE.to_string(),
        Value::Bool(flag) => flag.to_string(),
        Value::Int(n) => n.to_string(),
        Value::Double(d) => crate::variant::format_vcf_double(*d),
        // `AbstractList.toString`, which separates with ", " and not ",".
        Value::List(items) => format!(
            "[{}]",
            items.iter().map(value_of).collect::<Vec<_>>().join(", ")
        ),
    }
}

/// `getAttributeAsString(key, default)`.
pub fn as_string(value: Option<&Value>, default: &str) -> String {
    match value {
        None => default.to_string(),
        Some(value) => value_of(value),
    }
}

/// `getAttributeAsInt(key, default)`.
///
/// The missing test is the reference comparison described in the module doc, so only
/// [`Value::Missing`] takes the default and a `Value::Str(".")` reaches the parser.
pub fn as_int(value: Option<&Value>, default: i32) -> Result<i32, AttributeError> {
    match value {
        None | Some(Value::Missing) => Ok(default),
        Some(Value::Int(n)) => Ok(*n as i32),
        // Anything that is not a String reaches `(String) x` and fails the cast rather than the
        // parse, so a flag and a list produce two different exceptions with two different classes.
        Some(other @ (Value::List(_) | Value::Bool(_) | Value::Double(_))) => {
            Err(class_cast(java_class(other)))
        }
        Some(other) => {
            let text = value_of(other);
            text.parse::<i32>().map_err(|_| number_format(&text))
        }
    }
}

/// `getAttributeAsDouble(key, default)`.
///
/// **No missing test.** A `Value::Missing` is not null, so it reaches the parser and throws, which
/// is the disagreement with [`as_int`] the module doc names.
pub fn as_double(value: Option<&Value>, default: f64) -> Result<f64, AttributeError> {
    match value {
        None => Ok(default),
        Some(Value::Double(d)) => Ok(*d),
        Some(Value::Int(n)) => Ok(*n as f64),
        Some(other @ (Value::List(_) | Value::Bool(_))) => Err(class_cast(java_class(other))),
        Some(other) => parse_vcf_double(&value_of(other)),
    }
}

/// `getAttributeAsBoolean(key, default)`, which is `Boolean.valueOf`: `"true"` ignoring case, and
/// **false for everything else**, `"1"` and `"TRUE "` included, with nothing reporting that the
/// text was not a boolean.
pub fn as_boolean(value: Option<&Value>, default: bool) -> Result<bool, AttributeError> {
    match value {
        None => Ok(default),
        Some(Value::Bool(flag)) => Ok(*flag),
        Some(other @ (Value::List(_) | Value::Int(_) | Value::Double(_))) => {
            Err(class_cast(java_class(other)))
        }
        Some(other) => Ok(value_of(other).eq_ignore_ascii_case("true")),
    }
}

/// `getAttributeAsList(key)`: a list as itself, a scalar as a singleton, an absent key as empty.
pub fn as_list(value: Option<&Value>) -> Vec<Value> {
    match value {
        None => Vec::new(),
        Some(Value::List(items)) => items.clone(),
        Some(other) => vec![other.clone()],
    }
}

/// `getAttributeAsStringList(key, default)`.
pub fn as_string_list(value: Option<&Value>, default: &str) -> Vec<String> {
    // The transformer's `x == null` arm is unreachable from a decoded record: a list holds parsed
    // tokens and never a null. The default is carried anyway, because a caller supplies it.
    let _ = default;
    as_list(value).iter().map(value_of).collect()
}

/// `getAttributeAsIntList(key, default)`.
///
/// Each element gets the same reference-comparison missing test as [`as_int`], so a list written
/// `IA=1,.,3` throws on its second element: the tokens of a list are always substrings.
pub fn as_int_list(value: Option<&Value>, default: i32) -> Result<Vec<i32>, AttributeError> {
    as_list(value)
        .iter()
        .map(|item| match item {
            Value::Missing => Ok(default),
            // `x instanceof Number` is the only pre-cast test, so a boolean element reaches
            // `(String) x` and fails the cast rather than the parse.
            Value::Int(n) => Ok(*n as i32),
            other @ (Value::Bool(_) | Value::Double(_) | Value::List(_)) => {
                Err(class_cast(java_class(other)))
            }
            other => {
                let text = value_of(other);
                text.parse::<i32>().map_err(|_| number_format(&text))
            }
        })
        .collect()
}

/// `getAttributeAsDoubleList(key, default)`.
///
/// Unlike [`as_double`], this one **does** have the missing test, because the transformer written
/// for the list case carries it and the scalar accessor was written separately. Two accessors over
/// the same conversion, one with the test and one without.
pub fn as_double_list(value: Option<&Value>, default: f64) -> Result<Vec<f64>, AttributeError> {
    as_list(value)
        .iter()
        .map(|item| match item {
            Value::Missing => Ok(default),
            Value::Double(d) => Ok(*d),
            Value::Int(n) => Ok(*n as f64),
            other @ (Value::Bool(_) | Value::List(_)) => Err(class_cast(java_class(other))),
            other => parse_vcf_double(&value_of(other)),
        })
        .collect()
}

/// `VCFUtils.parseVcfDouble`.
///
/// `Double.parseDouble` first, then the infinity-or-NaN pattern. Both halves are wider than the
/// format: the first accepts whitespace, a trailing `f`/`d` and hex float literals, and the second
/// accepts `inf`, `infinity` and `nan` in any case with an optional sign. `-nan` is `NaN`, not a
/// negative one, because the sign group is read only on the infinity branch.
pub fn parse_vcf_double(text: &str) -> Result<f64, AttributeError> {
    if let Some(value) = java_parse_double(text) {
        return Ok(value);
    }
    let trimmed = text;
    let (sign, rest) = match trimmed.as_bytes().first() {
        Some(b'-') => (-1.0, &trimmed[1..]),
        Some(b'+') => (1.0, &trimmed[1..]),
        _ => (1.0, trimmed),
    };
    if rest.eq_ignore_ascii_case("inf") || rest.eq_ignore_ascii_case("infinity") {
        return Ok(sign * f64::INFINITY);
    }
    if rest.eq_ignore_ascii_case("nan") {
        return Ok(f64::NAN);
    }
    Err(if text.is_empty() {
        AttributeError::NumberFormat("empty String".to_string())
    } else {
        number_format(text)
    })
}

/// `Double.parseDouble`, whose grammar is not Rust's.
///
/// The differences that reach a VCF, measured: leading and trailing whitespace are **trimmed**, a
/// single trailing `f`, `F`, `d` or `D` is **allowed**, and a hexadecimal float literal such as
/// `0x1p3` parses to 8. Rust's `f64::from_str` accepts none of those and accepts `inf`/`NaN`, which
/// Java's does not — so those are stripped out here and left to the caller's pattern, which is
/// where upstream handles them.
fn java_parse_double(text: &str) -> Option<f64> {
    let trimmed = text.trim_matches(|c: char| c <= ' ');
    if trimmed.is_empty() {
        return None;
    }
    // Rust accepts these spellings and Java does not; letting them through here would make the
    // pattern below unreachable and change which branch a value took, which is observable in
    // nothing today and would be a silent divergence tomorrow.
    let bare = trimmed.trim_start_matches(['-', '+']);
    if bare.eq_ignore_ascii_case("inf")
        || bare.eq_ignore_ascii_case("infinity")
        || bare.eq_ignore_ascii_case("nan")
    {
        return None;
    }

    let body = match trimmed.as_bytes().last() {
        Some(b'f' | b'F' | b'd' | b'D') if !is_hex_literal(trimmed) => {
            &trimmed[..trimmed.len() - 1]
        }
        _ => trimmed,
    };
    if is_hex_literal(body) {
        return parse_hex_float(body);
    }
    body.parse::<f64>().ok()
}

fn is_hex_literal(text: &str) -> bool {
    let bare = text.trim_start_matches(['-', '+']);
    bare.len() > 2 && (bare.starts_with("0x") || bare.starts_with("0X"))
}

/// Java's hexadecimal floating-point literal: `0x<hex>[.<hex>]p<decimal exponent>`.
///
/// The binary exponent is mandatory, which is what makes `0x1` a refusal and `0x1p3` a number.
fn parse_hex_float(text: &str) -> Option<f64> {
    let (sign, rest) = match text.as_bytes().first() {
        Some(b'-') => (-1.0, &text[1..]),
        Some(b'+') => (1.0, &text[1..]),
        _ => (1.0, text),
    };
    let rest = rest.get(2..)?;
    let (mantissa, exponent) = rest.split_once(['p', 'P'])?;
    let (whole, fraction) = match mantissa.split_once('.') {
        Some((whole, fraction)) => (whole, fraction),
        None => (mantissa, ""),
    };
    if whole.is_empty() && fraction.is_empty() {
        return None;
    }
    let mut value = 0.0f64;
    for digit in whole.chars() {
        value = value * 16.0 + f64::from(digit.to_digit(16)?);
    }
    let mut scale = 1.0 / 16.0;
    for digit in fraction.chars() {
        value += f64::from(digit.to_digit(16)?) * scale;
        scale /= 16.0;
    }
    let exponent: i32 = exponent.parse().ok()?;
    Some(sign * value * 2f64.powi(exponent))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The three spellings of missing, and the two outcomes.
    #[test]
    fn only_the_constant_takes_the_default() {
        assert_eq!(as_int(Some(&Value::Missing), -1), Ok(-1));
        assert_eq!(
            as_int(Some(&Value::Str(".".into())), -1),
            Err(number_format("."))
        );
        assert_eq!(as_string(Some(&Value::Missing), "D"), ".");
        assert_eq!(as_string(Some(&Value::Str(".".into())), "D"), ".");
    }

    /// And the accessor that has no test at all.
    #[test]
    fn as_double_refuses_every_spelling_of_missing() {
        assert!(as_double(Some(&Value::Missing), -1.0).is_err());
        assert!(as_double(Some(&Value::Str(".".into())), -1.0).is_err());
        assert_eq!(as_double(None, -1.0), Ok(-1.0));
    }

    /// The list accessor has the test the scalar one lacks.
    #[test]
    fn the_list_accessors_disagree_with_the_scalar_ones() {
        assert_eq!(as_double_list(Some(&Value::Missing), -1.0), Ok(vec![-1.0]));
        assert!(as_double(Some(&Value::Missing), -1.0).is_err());
    }

    #[test]
    fn a_scalar_accessor_on_a_list_is_a_cast_failure_and_not_a_parse_failure() {
        let list = Value::List(vec![Value::Str("1".into()), Value::Str("2".into())]);
        assert_eq!(
            as_int(Some(&list), -1).unwrap_err().class(),
            "java.lang.ClassCastException"
        );
        // A flag casts too, and names its own class rather than the list one.
        assert!(as_int(Some(&Value::Bool(true)), -1)
            .unwrap_err()
            .message()
            .starts_with("class java.lang.Boolean"));
    }

    #[test]
    fn parse_vcf_double_is_wider_than_the_format_at_both_ends() {
        assert_eq!(parse_vcf_double("1f"), Ok(1.0));
        assert_eq!(parse_vcf_double("0x1p3"), Ok(8.0));
        assert_eq!(parse_vcf_double(" 1"), Ok(1.0));
        assert_eq!(parse_vcf_double("inf"), Ok(f64::INFINITY));
        assert_eq!(parse_vcf_double("-INF"), Ok(f64::NEG_INFINITY));
        assert!(parse_vcf_double("nan").unwrap().is_nan());
        // `-nan` is NaN and not a signed one: the sign is read on the infinity branch only.
        assert!(parse_vcf_double("-nan").unwrap().is_nan());
        assert_eq!(
            parse_vcf_double(""),
            Err(AttributeError::NumberFormat("empty String".into()))
        );
    }

    /// `Boolean.valueOf`, which never fails and is almost always false.
    #[test]
    fn as_boolean_is_true_only_for_the_word() {
        assert_eq!(
            as_boolean(Some(&Value::Str("true".into())), false),
            Ok(true)
        );
        assert_eq!(
            as_boolean(Some(&Value::Str("TRUE".into())), false),
            Ok(true)
        );
        assert_eq!(as_boolean(Some(&Value::Str("1".into())), false), Ok(false));
        assert_eq!(as_boolean(Some(&Value::Bool(true)), false), Ok(true));
    }

    /// A scalar is a list of one and an absent key is empty, not a list of the default.
    #[test]
    fn a_scalar_read_as_a_list_is_a_list_of_one() {
        assert_eq!(as_list(Some(&Value::Str("x".into()))).len(), 1);
        assert!(as_list(None).is_empty());
        assert!(as_int_list(None, -1).unwrap().is_empty());
    }
}
