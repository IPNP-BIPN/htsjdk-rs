//! `GenotypeLikelihoods` and the double parser it stands on, ported from
//! `htsjdk.variant.variantcontext.GenotypeLikelihoods` and `htsjdk.variant.vcf.VCFUtils`
//! (htsjdk 4.2.0).
//!
//! This is what the record decoder needs in order to accept a `GL` field: a genotype carrying `GL`
//! and no `PL` gets its `PL` computed here, and every downstream tool reads the `PL`. Until this
//! landed, `genotype_parse` refused a `GL` key rather than approximating it.
//!
//! # `Double.parseDouble` is not `str::parse::<f64>`
//!
//! The two agree on ordinary decimals and disagree in four places, all of which appear in real VCFs
//! and all of which the reference accepts:
//!
//! * a trailing type suffix: `1.5f`, `1.5F`, `1.5d`, `1.5D` are doubles in Java and errors in Rust;
//! * hexadecimal floating point: `0x1p3` is 8.0 in Java and an error in Rust;
//! * leading and trailing whitespace, which Java trims and Rust refuses;
//! * `Infinity` spelled out, which Java accepts and Rust does not (Rust wants `inf`).
//!
//! On top of that, `VCFUtils.parseVcfDouble` catches the failure and retries against
//! `^(?<sign>[-+]?)((?<inf>(INF|INFINITY))|(?<nan>NAN))$`, case-insensitive, so `inf`, `-INFINITY`
//! and `nan` parse too. The wrapper exists because those spellings are in the wild; the four
//! behaviours above come free with `Double.parseDouble` and are the easier ones to miss.
//!
//! # The conversion re-normalises, and the clamp happens before the rounding
//!
//! ```java
//! pls[i] = (int) Math.round(Math.min(-10 * (GLs[i] - adjust), MAX_PL));
//! ```
//!
//! `adjust` is the **maximum** of the likelihoods, so the best genotype always comes out at PL 0
//! whatever the input scale. `MAX_PL` is `Integer.MAX_VALUE`, the clamp is applied to the `double`
//! **before** rounding, and `Math.round` returns a `long` that is then truncated to `int`. Each of
//! those three steps changes the answer at the extremes, and a port that reordered them would agree
//! on every ordinary genotype.

/// `GenotypeLikelihoods.MAX_PL`.
pub const MAX_PL: f64 = i32::MAX as f64;

/// `VCFConstants.MISSING_VALUE_v4`.
pub const MISSING_VALUE: &str = ".";

/// What the GL field can fail with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LikelihoodError {
    /// `TribbleException("partial missing values for GL field")`.
    PartialMissing,
    /// `Double.parseDouble` refusing a value, after `parseVcfDouble` has tried its pattern.
    NumberFormat(String),
}

impl LikelihoodError {
    pub fn class(&self) -> &'static str {
        match self {
            LikelihoodError::PartialMissing => "htsjdk.tribble.TribbleException",
            LikelihoodError::NumberFormat(_) => "java.lang.NumberFormatException",
        }
    }

    pub fn message(&self) -> String {
        match self {
            LikelihoodError::PartialMissing => "partial missing values for GL field".to_string(),
            LikelihoodError::NumberFormat(value) => format!("For input string: \"{value}\""),
        }
    }
}

/// `Double.parseDouble`, including the three spellings Rust's parser refuses.
///
/// Written out rather than delegated because the differences are not edge cases in this context: a
/// GL field is written by whatever produced the VCF, and `1.5f` or a leading space is exactly the
/// kind of thing a hand-edited file carries.
pub fn parse_java_double(text: &str) -> Option<f64> {
    // Java trims the ASCII control characters and space, which is `String.trim()`, not Unicode
    // whitespace: a non-breaking space is *not* trimmed and the parse fails.
    let trimmed: &str = {
        let bytes = text.as_bytes();
        let start = bytes.iter().position(|byte| *byte > 0x20);
        let end = bytes.iter().rposition(|byte| *byte > 0x20);
        match (start, end) {
            (Some(start), Some(end)) => &text[start..=end],
            _ => "",
        }
    };
    if trimmed.is_empty() {
        return None;
    }

    // A single trailing type suffix, which Java allows on both decimal and hexadecimal literals.
    let body = match trimmed.chars().last() {
        Some('f') | Some('F') | Some('d') | Some('D') => {
            let stripped = &trimmed[..trimmed.len() - 1];
            // `Infinity` ends in a `y`, but `NaN` does not end in a suffix and `1d` does. A body
            // left empty, or one that is only a sign, is not a number.
            if stripped.is_empty() || stripped == "-" || stripped == "+" {
                return None;
            }
            stripped
        }
        _ => trimmed,
    };

    // Java spells the infinities out and is case-sensitive about it here; the case-insensitive
    // spellings are `parseVcfDouble`'s job, not this function's.
    match body {
        "Infinity" | "+Infinity" => return Some(f64::INFINITY),
        "-Infinity" => return Some(f64::NEG_INFINITY),
        "NaN" | "+NaN" | "-NaN" => return Some(f64::NAN),
        _ => {}
    }

    if let Some(value) = parse_hexadecimal(body) {
        return Some(value);
    }

    // Rust's parser accepts `inf`, `infinity` and `nan` case-insensitively, which Java does not, so
    // they are refused here and left to `parse_vcf_double` to accept where the reference does.
    let lowered = body.to_ascii_lowercase();
    let lowered = lowered
        .strip_prefix('-')
        .or_else(|| lowered.strip_prefix('+'))
        .unwrap_or(&lowered);
    if lowered == "inf" || lowered == "infinity" || lowered == "nan" {
        return None;
    }
    body.parse::<f64>().ok()
}

/// Java's hexadecimal floating-point literal: `0x1.8p3`, where the `p` exponent is mandatory.
fn parse_hexadecimal(body: &str) -> Option<f64> {
    let (sign, rest) = match body.strip_prefix('-') {
        Some(rest) => (-1.0, rest),
        None => (1.0, body.strip_prefix('+').unwrap_or(body)),
    };
    let rest = rest
        .strip_prefix("0x")
        .or_else(|| rest.strip_prefix("0X"))?;
    let (mantissa, exponent) = rest.split_once(['p', 'P'])?;
    let exponent: i32 = exponent.parse().ok()?;

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
    let mut scale = 1.0f64 / 16.0;
    for digit in fraction.chars() {
        value += f64::from(digit.to_digit(16)?) * scale;
        scale /= 16.0;
    }
    Some(sign * value * 2f64.powi(exponent))
}

/// `VCFUtils.parseVcfDouble`: `Double.parseDouble`, then the infinity-or-NaN pattern.
///
/// The pattern is tried **only** after the plain parse has failed, so a string both parsers accept
/// takes the plain parse's answer. It matters for nothing today and is the reference's order.
pub fn parse_vcf_double(text: &str) -> Option<f64> {
    if let Some(value) = parse_java_double(text) {
        return Some(value);
    }
    let (sign, body) = match text.strip_prefix('-') {
        Some(body) => (-1.0, body),
        None => (1.0, text.strip_prefix('+').unwrap_or(text)),
    };
    match body.to_ascii_uppercase().as_str() {
        "INF" | "INFINITY" => Some(sign * f64::INFINITY),
        // The sign is captured and then ignored for NaN, exactly as upstream ignores it.
        "NAN" => Some(f64::NAN),
        _ => None,
    }
}

/// `GenotypeLikelihoods.parseDeprecatedGLString`.
///
/// `None` is the reference's `null`: either the whole field is `.`, or every element is. A field
/// with **some** elements missing is the refusal, and the count is what decides which of the three
/// it is.
pub fn parse_gl_field(text: &str) -> Result<Option<Vec<f64>>, LikelihoodError> {
    if text == MISSING_VALUE {
        return Ok(None);
    }
    // `String.split(",")` on the whole field. A trailing empty element is dropped by Java's split,
    // which is why "1,2," has three values here and two there.
    let parts = java_split_on_comma(text);
    let mut values = vec![0.0f64; parts.len()];
    let mut missing = 0usize;
    for (index, part) in parts.iter().enumerate() {
        if *part == MISSING_VALUE {
            missing += 1;
        } else {
            values[index] = parse_vcf_double(part)
                .ok_or_else(|| LikelihoodError::NumberFormat((*part).to_string()))?;
        }
    }
    if missing == 0 {
        Ok(Some(values))
    } else if missing == values.len() {
        Ok(None)
    } else {
        Err(LikelihoodError::PartialMissing)
    }
}

/// `String.split(",")`: trailing empty strings are removed, leading and interior ones are kept.
fn java_split_on_comma(text: &str) -> Vec<&str> {
    let mut parts: Vec<&str> = text.split(',').collect();
    while parts.len() > 1 && parts.last() == Some(&"") {
        parts.pop();
    }
    if parts.len() == 1 && parts[0].is_empty() {
        parts.clear();
    }
    parts
}

/// `GenotypeLikelihoods.GLsToPLs`.
///
/// Three things in one line, and their order is the behaviour: the likelihoods are shifted so the
/// largest becomes zero, the shifted value is clamped to `Integer.MAX_VALUE` **as a double**, and
/// only then rounded. `Math.round` is `floor(x + 0.5)`, so it rounds half **up** rather than half
/// away from zero, and it answers 0 for a NaN.
pub fn gls_to_pls(likelihoods: &[f64]) -> Vec<i32> {
    let adjust = max_pl(likelihoods);
    likelihoods
        .iter()
        .map(|likelihood| {
            let scaled = (-10.0 * (likelihood - adjust)).min(MAX_PL);
            java_round(scaled) as i32
        })
        .collect()
}

/// `GenotypeLikelihoods.maxPL`: a plain maximum starting from negative infinity.
///
/// `Math.max` propagates NaN, so a single NaN in the field makes the adjustment NaN and every PL
/// comes out 0. That is the reference's answer and not an accident this port smooths over.
fn max_pl(likelihoods: &[f64]) -> f64 {
    let mut adjust = f64::NEG_INFINITY;
    for likelihood in likelihoods {
        adjust = java_math_max(adjust, *likelihood);
    }
    adjust
}

/// `Math.max(double, double)`, which returns NaN when either argument is NaN.
///
/// Rust's `f64::max` returns the **other** operand instead, so it cannot be used here: one NaN in a
/// GL field would silently disappear and the PLs would be computed from the remaining values.
fn java_math_max(left: f64, right: f64) -> f64 {
    if left.is_nan() || right.is_nan() {
        return f64::NAN;
    }
    if left > right {
        left
    } else if right > left {
        right
    } else if left == 0.0 && right == 0.0 {
        // Math.max(-0.0, 0.0) is 0.0.
        if left.is_sign_positive() {
            left
        } else {
            right
        }
    } else {
        left
    }
}

/// `Math.round(double)`: `floor(x + 0.5)` as a `long`, with the saturating edges Java specifies.
///
/// The half-up rule is the visible part: -1.5 rounds to -1, not to -2. The saturation matters at
/// the other end, where an infinity becomes `Long.MAX_VALUE` and truncating that to `int` gives
/// -1, which is why the clamp to `MAX_PL` has to happen first.
pub fn java_round(value: f64) -> i64 {
    if value.is_nan() {
        return 0;
    }
    let shifted = (value + 0.5).floor();
    if shifted >= i64::MAX as f64 {
        i64::MAX
    } else if shifted <= i64::MIN as f64 {
        i64::MIN
    } else {
        shifted as i64
    }
}

/// `GenotypeLikelihoods.fromGLField().getAsPLs()`, which is what the record decoder calls.
///
/// `None` where the reference returns `null`: `getAsPLs` propagates the null vector rather than
/// producing an empty array.
pub fn gl_field_to_pls(text: &str) -> Result<Option<Vec<i32>>, LikelihoodError> {
    Ok(parse_gl_field(text)?.map(|likelihoods| gls_to_pls(&likelihoods)))
}

/// `GenotypeLikelihoods.PLsToGLs`: the other direction, a plain division by -10.
pub fn pls_to_gls(pls: &[i32]) -> Vec<f64> {
    pls.iter().map(|pl| f64::from(*pl) / -10.0).collect()
}
