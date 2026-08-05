//! Reading a VCF header, as opposed to writing one.
//!
//! Ported from `htsjdk.variant.vcf.VCFHeaderLineTranslator`'s `VCF4Parser`, from
//! `htsjdk.variant.vcf.VCFCodec.readActualHeader` and from the `#CHROM` handling of
//! `htsjdk.variant.vcf.AbstractVCFCodec.parseHeaderFromLines` (htsjdk 4.2.0).
//!
//! This is the frame of the header: which version the file declares, which lines are header lines,
//! and which samples the column line names. Turning an `##INFO=<...>` line into a typed compound
//! header line is the next layer and is deliberately not here, so this slice can be measured on its
//! own. The dump feeds only unstructured meta lines for the same reason.
//!
//! # The state machine swallows characters
//!
//! `VCF4Parser.parseLine` is a hand-written scanner, and its `switch` falls through:
//!
//! ```java
//! case ('<') : if (index == 0) break;  // no break when index != 0: falls into '>'
//! case ('>') : if (index == valueLine.length()-1) ret.put(key,builder.toString().trim()); break;
//! case ('=') : key = builder.toString().trim(); builder = new StringBuilder(); break;
//! case (',') : ret.put(key,builder.toString().trim()); builder = new StringBuilder(); break;
//! default: builder.append(c);
//! ```
//!
//! Three consequences, none of them in the VCF specification:
//!
//!  * an unquoted `<` anywhere but position 0 falls into the `>` case and is **dropped**, and if it
//!    is the last character it also closes the entry;
//!  * an unquoted `>` anywhere is **dropped**, and only the one at the very end stores the pending
//!    entry;
//!  * therefore a line that does not end in `>` **loses its last field entirely**, silently.
//!
//! # Quotes are a toggle, not a delimiter
//!
//! The `c == '"'` test comes before the `inQuote` test, so a quote anywhere flips the state,
//! including one in the middle of an unquoted value. Inside a quote, `\"` yields `"`, `\\` yields
//! `\`, and a backslash before anything else is **kept along with the character**, so `\n` stays
//! two characters. An unclosed quote is a refusal rather than a truncation.
//!
//! # The version line is split on every `=`
//!
//! `readActualHeader` does `line.substring(2).split("=")` and only records a version when the
//! result has exactly two fields. `##fileformat=VCFv4.2=x` therefore parses as three fields, no
//! version is recorded, and the failure surfaces later as "we never saw a header line specifying
//! VCF version" rather than as a complaint about the line itself.

use std::fmt;

/// `VCFHeaderVersion`.
///
/// `VCF3_2`'s version string is `VCRv3.2` upstream, letters transposed. It is reproduced because a
/// file carrying that exact string is recognised and one carrying `VCFv3.2` is not.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VcfVersion {
    Vcf3_2,
    Vcf3_3,
    Vcf4_0,
    Vcf4_1,
    Vcf4_2,
    Vcf4_3,
    Vcf4_4,
}

impl VcfVersion {
    /// `VCFHeaderVersion.getVersionString`.
    pub fn version_string(self) -> &'static str {
        match self {
            VcfVersion::Vcf3_2 => "VCRv3.2",
            VcfVersion::Vcf3_3 => "VCFv3.3",
            VcfVersion::Vcf4_0 => "VCFv4.0",
            VcfVersion::Vcf4_1 => "VCFv4.1",
            VcfVersion::Vcf4_2 => "VCFv4.2",
            VcfVersion::Vcf4_3 => "VCFv4.3",
            VcfVersion::Vcf4_4 => "VCFv4.4",
        }
    }

    /// `VCFHeaderVersion.getFormatString`: `format` for 3.2, `fileformat` for everything else.
    pub fn format_string(self) -> &'static str {
        match self {
            VcfVersion::Vcf3_2 => "format",
            _ => "fileformat",
        }
    }

    /// `VCFHeaderVersion.toHeaderVersion`.
    pub fn from_version_string(version: &str) -> Option<VcfVersion> {
        [
            VcfVersion::Vcf3_2,
            VcfVersion::Vcf3_3,
            VcfVersion::Vcf4_0,
            VcfVersion::Vcf4_1,
            VcfVersion::Vcf4_2,
            VcfVersion::Vcf4_3,
            VcfVersion::Vcf4_4,
        ]
        .into_iter()
        .find(|candidate| candidate.version_string() == version)
    }

    /// `VCFHeaderVersion.isFormatString`.
    pub fn is_format_string(format: &str) -> bool {
        format == "format" || format == "fileformat"
    }

    /// The enum constant's own name, which is what `String.format("%s", version)` produces.
    ///
    /// Not the same text as [`Self::version_string`], and the difference is observable: the
    /// writer's refusal of a 4.3 header says `VCF4_3` where the file said `VCFv4.3`.
    pub fn constant_name(self) -> &'static str {
        match self {
            VcfVersion::Vcf3_2 => "VCF3_2",
            VcfVersion::Vcf3_3 => "VCF3_3",
            VcfVersion::Vcf4_0 => "VCF4_0",
            VcfVersion::Vcf4_1 => "VCF4_1",
            VcfVersion::Vcf4_2 => "VCF4_2",
            VcfVersion::Vcf4_3 => "VCF4_3",
            VcfVersion::Vcf4_4 => "VCF4_4",
        }
    }
}

/// `TribbleException.InvalidHeader`, carrying the message because the message is the behaviour: it
/// names which tag was out of order and what was expected.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InvalidHeader(pub String);

/// The prefix `TribbleException.InvalidHeader` puts in front of every message it is given.
pub const INVALID_HEADER_PREFIX: &str = "Your input file has a malformed header: ";

impl InvalidHeader {
    /// What `getMessage()` returns: the reason with the exception's own prefix in front.
    pub fn message(&self) -> String {
        format!("{INVALID_HEADER_PREFIX}{}", self.0)
    }
}

impl fmt::Display for InvalidHeader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message())
    }
}

impl std::error::Error for InvalidHeader {}

/// `VCF4Parser.parseLine`.
///
/// The returned pairs are in insertion order, and a repeated key **overwrites the value while
/// keeping the original position**, because upstream collects into a `LinkedHashMap`.
///
/// `expected_tag_order` of `None` is Java's `null`: no validation at all. `Some(&[])` is an empty
/// list, which validates and accepts everything, and is not the same thing.
pub fn parse_structured_value(
    value_line: &str,
    expected_tag_order: Option<&[&str]>,
    recommended_tags: &[&str],
) -> Result<Vec<(String, String)>, InvalidHeader> {
    let mut pairs: Vec<(String, String)> = Vec::new();
    let mut builder = String::new();
    let mut key = String::new();
    let mut in_quote = false;
    let mut escape = false;

    // `valueLine.length()` counts UTF-16 units upstream, and `toCharArray` walks them, so the
    // "last character" test is over the same units. Indices here are over `char`s, which agree for
    // everything outside the supplementary planes; a header line carrying an astral character
    // would need the UTF-16 count instead.
    let chars: Vec<char> = value_line.chars().collect();
    let last = chars.len().saturating_sub(1);

    for (index, &c) in chars.iter().enumerate() {
        if c == '"' {
            // Checked before `in_quote`, so a quote toggles the state wherever it appears.
            if escape {
                builder.push(c);
                escape = false;
            } else {
                in_quote = !in_quote;
            }
        } else if in_quote {
            if escape {
                // Only `\` and `"` are escapable; anything else keeps the backslash.
                if c == '\\' {
                    builder.push(c);
                } else {
                    builder.push('\\');
                    builder.push(c);
                }
                escape = false;
            } else if c != '\\' {
                builder.push(c);
            } else {
                escape = true;
            }
        } else {
            escape = false;
            match c {
                // The fall-through: a `<` away from position 0 is handled as a `>`.
                '<' if index == 0 => {}
                '<' | '>' => {
                    if index == last {
                        put(&mut pairs, &key, builder.trim());
                        builder = String::new();
                    }
                }
                '=' => {
                    key = builder.trim().to_string();
                    builder = String::new();
                }
                ',' => {
                    put(&mut pairs, &key, builder.trim());
                    builder = String::new();
                }
                other => builder.push(other),
            }
        }
    }

    if in_quote {
        return Err(InvalidHeader(format!(
            "Unclosed quote in header line value {value_line}"
        )));
    }

    if let Some(expected) = expected_tag_order {
        if pairs.is_empty() && !expected.is_empty() {
            return Err(InvalidHeader(format!(
                "Header with no tags is not supported when there are expected tags in line \
                 {value_line}"
            )));
        }
        for (index, (tag, _)) in pairs.iter().enumerate() {
            if index >= expected.len() {
                continue;
            }
            if expected[index] == tag {
                continue;
            }
            return Err(if let Some(at) = expected.iter().position(|e| e == tag) {
                InvalidHeader(format!(
                    "Tag {tag} in wrong order (was #{}, expected #{}) in line {value_line}",
                    index + 1,
                    at + 1
                ))
            } else if recommended_tags.contains(&tag.as_str()) {
                InvalidHeader(format!(
                    "Recommended tag {tag} must be listed after all expected tags in line \
                     {value_line}"
                ))
            } else {
                InvalidHeader(format!("Unexpected tag {tag} in line {value_line}"))
            });
        }
    }

    Ok(pairs)
}

/// `LinkedHashMap.put`: replace in place, keep the original position.
fn put(pairs: &mut Vec<(String, String)>, key: &str, value: &str) {
    if let Some(slot) = pairs.iter_mut().find(|(k, _)| k == key) {
        slot.1 = value.to_string();
    } else {
        pairs.push((key.to_string(), value.to_string()));
    }
}

/// What `readActualHeader` establishes before any line is turned into a typed header line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HeaderFrame {
    /// The version the file declared, which is not optional: a file that never declares one is a
    /// refusal rather than a default.
    pub version: VcfVersion,
    /// The `##` lines, in file order, with the `##` still on them, exactly as the codec collects
    /// them before handing them on.
    pub meta_lines: Vec<String>,
    /// The sample names from the `#CHROM` line, deduplicated with the first occurrence kept.
    pub samples: Vec<String>,
}

impl HeaderFrame {
    /// The keys of the metadata `VCFHeader` ends up holding, in input order.
    ///
    /// Three upstream steps land here and none of them is in the file format:
    ///
    ///  * a `##` line with no `=` at all is **dropped**, because `parseHeaderFromLines` only adds
    ///    one when `str.indexOf('=') != -1`;
    ///  * the lines go into a `LinkedHashSet`, and `VCFHeaderLine.equals` is over key **and**
    ///    value, so a line repeated exactly collapses and one repeated with a different value does
    ///    not;
    ///  * the `fileformat` line **stays**. `VCFHeader.removeVCFVersionLines` exists, but the
    ///    constructor the codec calls does not reach it, so the header carries the version both as
    ///    a field and as the line it came from. Measured rather than assumed: the golden's
    ///    `minimal` frame lists `fileformat` among its keys. Nothing appears twice on output,
    ///    because the writer skips every line whose key is a format string and writes its own.
    pub fn meta_keys(&self) -> Vec<String> {
        let mut seen: Vec<(String, String)> = Vec::new();
        for line in &self.meta_lines {
            let rest = &line[2..];
            let Some(equals) = rest.find('=') else {
                continue;
            };
            let (key, value) = (rest[..equals].to_string(), rest[equals + 1..].to_string());
            if !seen
                .iter()
                .any(|pair| *pair == (key.clone(), value.clone()))
            {
                seen.push((key, value));
            }
        }
        seen.into_iter().map(|(key, _)| key).collect()
    }
}

/// `VCFHeader.HEADER_FIELDS`, the eight mandatory columns in their required order.
pub const HEADER_FIELDS: [&str; 8] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

/// `VCFCodec.readActualHeader` plus the `#CHROM` handling of `parseHeaderFromLines`.
///
/// Every refusal upstream is a refusal here, with its message, because the messages distinguish
/// cases a caller can act on: "never saw a header line specifying VCF version" and "never saw the
/// required CHROM header line" are different failures of the same file.
pub fn read_header_frame(text: &str) -> Result<HeaderFrame, InvalidHeader> {
    let mut meta_lines: Vec<String> = Vec::new();
    let mut version: Option<VcfVersion> = None;

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("##") {
            // `split("=")` with no limit, and a version is recorded only when that yields exactly
            // two fields. Java's split also drops trailing empty strings, so `##fileformat=` is one
            // field and records nothing.
            let fields = java_split_on_equals(rest);
            if fields.len() == 2 && VcfVersion::is_format_string(&fields[0]) {
                match VcfVersion::from_version_string(&fields[1]) {
                    None => {
                        return Err(InvalidHeader(format!(
                            "{} is not a supported version",
                            fields[1]
                        )));
                    }
                    Some(found) => {
                        version = Some(found);
                        // `Defaults.OPTIMISTIC_VCF_4_4` is off in the pinned oracle, so 4.4 is not
                        // downgraded to 4.3 and falls through to the codec's version check below.
                        if !matches!(
                            found,
                            VcfVersion::Vcf4_0
                                | VcfVersion::Vcf4_1
                                | VcfVersion::Vcf4_2
                                | VcfVersion::Vcf4_3
                        ) {
                            return Err(InvalidHeader(format!(
                                "This codec is strictly for VCFv4 and does not support {}",
                                fields[1]
                            )));
                        }
                    }
                }
            }
            meta_lines.push(line.to_string());
        } else if line.starts_with('#') {
            let Some(version) = version else {
                return Err(InvalidHeader(
                    "We never saw a header line specifying VCF version".to_string(),
                ));
            };
            let samples = parse_column_line(line)?;
            return Ok(HeaderFrame {
                version,
                meta_lines,
                samples,
            });
        } else {
            return Err(InvalidHeader(
                "We never saw the required CHROM header line (starting with one #) for the input \
                 VCF file"
                    .to_string(),
            ));
        }
    }

    Err(InvalidHeader(
        "We never saw the required CHROM header line (starting with one #) for the input VCF file"
            .to_string(),
    ))
}

/// The `#CHROM` branch of `parseHeaderFromLines`.
fn parse_column_line(line: &str) -> Result<Vec<String>, InvalidHeader> {
    // `substring(1)` then split on tab. Java's split drops trailing empty fields, so a line ending
    // in tabs has fewer columns than it looks like it has.
    let columns = java_split(&line[1..], '\t');
    if columns.len() < HEADER_FIELDS.len() {
        return Err(InvalidHeader(format!(
            "there are not enough columns present in the header line: {line}"
        )));
    }

    for (index, field) in HEADER_FIELDS.iter().enumerate() {
        let seen = &columns[index];
        // Upstream this is `HEADER_FIELDS.valueOf(seen)`, so an unknown name and a known name in
        // the wrong place are different messages.
        if !HEADER_FIELDS.contains(&seen.as_str()) {
            return Err(InvalidHeader(format!(
                "unknown column name '{seen}'; it does not match a legal column header name."
            )));
        }
        if field != seen {
            return Err(InvalidHeader(format!(
                "we were expecting column name '{field}' but we saw '{seen}'"
            )));
        }
    }

    let mut index = HEADER_FIELDS.len();
    let mut saw_format = false;
    if index < columns.len() {
        if columns[index] != "FORMAT" {
            return Err(InvalidHeader(format!(
                "we were expecting column name 'FORMAT' but we saw '{}'",
                columns[index]
            )));
        }
        saw_format = true;
        index += 1;
    }

    // A `LinkedHashSet`, so a repeated sample name collapses onto its first occurrence rather than
    // producing two columns' worth of one name.
    let mut samples: Vec<String> = Vec::new();
    while index < columns.len() {
        if !samples.contains(&columns[index]) {
            samples.push(columns[index].clone());
        }
        index += 1;
    }

    if saw_format && samples.is_empty() {
        return Err(InvalidHeader(
            "The FORMAT field was provided but there is no genotype/sample data".to_string(),
        ));
    }

    Ok(samples)
}

/// `String.split(String)` with no limit: trailing empty fields are dropped, leading and interior
/// ones are kept, and a string with no separator at all yields itself.
fn java_split(text: &str, separator: char) -> Vec<String> {
    let mut parts: Vec<String> = text.split(separator).map(str::to_string).collect();
    while parts.len() > 1 && parts.last().is_some_and(String::is_empty) {
        parts.pop();
    }
    // Java returns `[""]` for an empty input rather than an empty array, which the loop above
    // already preserves.
    parts
}

fn java_split_on_equals(text: &str) -> Vec<String> {
    java_split(text, '=')
}
