//! `BEDCodec`, ported from `htsjdk.tribble.bed.BEDCodec` (htsjdk 4.2.0).
//!
//! This is what GATK's `-L regions.bed` runs through, so its coordinate convention is the one a
//! whole interval argument inherits. Five of its decisions are not what "a BED parser" suggests.
//!
//! # The start is shifted, and the shift is a constructor argument
//!
//! ```java
//! int start = Integer.parseInt(tokens[1]) + startOffsetValue;
//! ```
//!
//! `StartOffset.ONE` is the default, so a BED file's 0-based start becomes 1-based on the way in.
//! A caller that constructs the codec with `StartOffset.ZERO` gets the file's own numbers, and
//! nothing downstream can tell which was used: the feature carries the shifted value and not the
//! shift. GATK builds the default, so `-L` sees 1-based starts.
//!
//! # A two-token line is a **point**, not a zero-length interval
//!
//! ```java
//! int end = start;
//! if (tokenCount > 2) { end = Integer.parseInt(tokens[2]); }
//! ```
//!
//! With no end column the end is the *shifted* start, so `chr1 10` becomes `chr1:11-11` under the
//! default offset. That is one base, at a coordinate the file never mentions.
//!
//! # The separator is a tab **or** a run of spaces
//!
//! ```java
//! Pattern.compile("\\t|( +)")
//! ```
//!
//! One tab splits, and so does any number of consecutive spaces, but a single space between two
//! tabs is a run of one and splits too. The limit is `-1`, so trailing empty fields are kept and a
//! line ending in a tab has one more token than it looks like.
//!
//! # A blank or header line is `null`, not an error
//!
//! `decode` returns null for an empty line and for one starting with `#`, `track` or `browser`.
//! The prefixes are matched on the raw line, so ` track` with a leading space is data.
//!
//! # A bad score is not a bad line
//!
//! ```java
//! try { float score = Float.parseFloat(tokens[4]); feature.setScore(score); }
//! catch (NumberFormatException e) { return feature; }
//! ```
//!
//! The feature is returned **early**, so a line whose score is `.` keeps its name and loses its
//! strand, its colour and its exons even though those columns are present and well formed. Every
//! other malformed number in the line throws instead.

/// `BEDCodec.StartOffset`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StartOffset {
    /// `StartOffset.ZERO`: the file's own numbers.
    Zero,
    /// `StartOffset.ONE`, the default, which is what GATK constructs.
    One,
}

impl StartOffset {
    pub fn value(self) -> i32 {
        match self {
            StartOffset::Zero => 0,
            StartOffset::One => 1,
        }
    }
}

/// `Strand`, as `FullBEDFeature` stores it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Positive,
    Negative,
    None,
}

impl Strand {
    /// `Strand.toString()`, which is what a dump can compare.
    pub fn name(self) -> &'static str {
        match self {
            Strand::Positive => "+",
            Strand::Negative => "-",
            Strand::None => "NONE",
        }
    }
}

/// One exon, as `FullBEDFeature.addExon` records it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Exon {
    pub start: i32,
    pub end: i32,
    pub cd_start: i32,
    pub cd_end: i32,
    pub number: i32,
}

/// `FullBEDFeature`, with only the fields the codec sets.
#[derive(Debug, Clone, PartialEq)]
pub struct BedFeature {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub name: Option<String>,
    pub score: Option<f32>,
    pub strand: Option<Strand>,
    /// `ParsingUtils.parseColor`'s answer, as `(r, g, b)`.
    pub color: Option<(u8, u8, u8)>,
    pub exons: Vec<Exon>,
}

/// What the codec throws, which is Java's own parse failure with its own message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NumberFormat {
    pub input: String,
}

impl NumberFormat {
    pub fn class(&self) -> &'static str {
        "java.lang.NumberFormatException"
    }

    pub fn message(&self) -> String {
        format!("For input string: \"{}\"", self.input)
    }
}

/// `Integer.parseInt`, which takes a leading `+` and refuses everything else Rust refuses.
fn parse_int(text: &str) -> Result<i32, NumberFormat> {
    text.parse::<i32>().map_err(|_| NumberFormat {
        input: text.to_string(),
    })
}

/// `SPLIT_PATTERN.split(line, -1)`: a tab, or a run of one or more spaces, with trailing empty
/// fields kept.
///
/// A run of spaces is one separator however long it is, so `"a   b"` is two fields; a tab beside a
/// space is two separators, so `"a \tb"` has an empty field between them.
pub fn split_bed_line(line: &str) -> Vec<&str> {
    let bytes = line.as_bytes();
    let mut fields = Vec::new();
    let mut start = 0usize;
    let mut index = 0usize;
    while index < bytes.len() {
        if bytes[index] == b'\t' {
            fields.push(&line[start..index]);
            index += 1;
            start = index;
        } else if bytes[index] == b' ' {
            fields.push(&line[start..index]);
            while index < bytes.len() && bytes[index] == b' ' {
                index += 1;
            }
            start = index;
        } else {
            index += 1;
        }
    }
    fields.push(&line[start..]);
    fields
}

/// `isBEDHeaderLine`: matched on the raw line, so a leading space makes it data.
pub fn is_header_line(line: &str) -> bool {
    line.starts_with('#') || line.starts_with("track") || line.starts_with("browser")
}

/// `String.trim().isEmpty()`, where `trim` drops everything at or below `U+0020`.
fn is_blank(line: &str) -> bool {
    line.bytes().all(|b| b <= b' ')
}

/// `BEDCodec.decode(String)`.
///
/// `Ok(None)` is the codec's `null`: a blank line, a header line, or a line with fewer than two
/// fields. None of the three is an error, and a reader that treated them as one would stop at the
/// first comment.
pub fn decode(line: &str, start_offset: StartOffset) -> Result<Option<BedFeature>, NumberFormat> {
    if is_blank(line) {
        return Ok(None);
    }
    if is_header_line(line) {
        return Ok(None);
    }
    decode_tokens(&split_bed_line(line), start_offset)
}

/// `BEDCodec.decode(String[])`.
pub fn decode_tokens(
    tokens: &[&str],
    start_offset: StartOffset,
) -> Result<Option<BedFeature>, NumberFormat> {
    if tokens.len() < 2 {
        return Ok(None);
    }

    let start = parse_int(tokens[1])? + start_offset.value();
    // With no end column the end is the shifted start: one base, at a coordinate the file does not
    // contain.
    let end = if tokens.len() > 2 {
        parse_int(tokens[2])?
    } else {
        start
    };

    let mut feature = BedFeature {
        contig: tokens[0].to_string(),
        start,
        end,
        name: None,
        score: None,
        strand: None,
        color: None,
        exons: Vec::new(),
    };

    if tokens.len() > 3 {
        // `replaceAll("\"", "")`: every quote anywhere, not a surrounding pair.
        feature.name = Some(tokens[3].replace('"', ""));
    }

    if tokens.len() > 4 {
        match parse_float(tokens[4]) {
            Some(score) => feature.score = Some(score),
            // The early return: everything after column five is dropped, however well formed.
            None => return Ok(Some(feature)),
        }
    }

    if tokens.len() > 5 {
        let trimmed = tokens[5].trim_matches(|c: char| c <= ' ');
        let strand = trimmed.chars().next().unwrap_or(' ');
        feature.strand = Some(match strand {
            '-' => Strand::Negative,
            '+' => Strand::Positive,
            _ => Strand::None,
        });
    }

    if tokens.len() > 8 {
        feature.color = parse_color(tokens[8]);
    }

    if tokens.len() > 11 {
        create_exons(start, tokens, &mut feature, start_offset)?;
    }

    Ok(Some(feature))
}

/// `Float.parseFloat`, whose failure the codec catches.
///
/// Java's parser takes a trailing `f`/`d` and hexadecimal floats and refuses the `inf`/`nan`
/// spellings Rust takes, exactly as `Double.parseDouble` does. Only whether it *succeeds* reaches
/// the feature here, plus the value when it does.
fn parse_float(text: &str) -> Option<f32> {
    let trimmed = text.trim_matches(|c: char| c <= ' ');
    if trimmed.is_empty() {
        return None;
    }
    let lowered = trimmed.to_ascii_lowercase();
    // Java takes a type suffix; Rust does not, and Rust takes `inf`/`nan`, which Java does not.
    let body = lowered
        .strip_suffix('f')
        .or_else(|| lowered.strip_suffix('d'))
        .unwrap_or(&lowered);
    if body.is_empty() || body.contains("inf") || body.contains("nan") {
        return None;
    }
    body.parse::<f32>().ok()
}

/// `ParsingUtils.parseColor`, in the two shapes a BED file uses.
///
/// A comma-separated triple is RGB; a `#rrggbb` string is hex. Anything else falls through the
/// named-colour table and, failing that, comes back **black** rather than absent, which is why an
/// unparsable colour is not distinguishable from a black one.
fn parse_color(text: &str) -> Option<(u8, u8, u8)> {
    if text.contains(',') {
        let parts: Vec<&str> = text.split(',').collect();
        if parts.len() < 3 {
            return Some((0, 0, 0));
        }
        let mut rgb = [0i32; 3];
        for (slot, part) in rgb.iter_mut().zip(&parts) {
            match part.trim().parse::<i32>() {
                Ok(value) => *slot = value,
                Err(_) => return Some((0, 0, 0)),
            }
        }
        if rgb.iter().any(|v| !(0..=255).contains(v)) {
            // `new Color(r, g, b)` throws for a component out of range, and the caller catches it.
            return Some((0, 0, 0));
        }
        return Some((rgb[0] as u8, rgb[1] as u8, rgb[2] as u8));
    }
    if let Some(hex) = text.strip_prefix('#') {
        if hex.len() == 6 {
            if let Ok(value) = u32::from_str_radix(hex, 16) {
                return Some((
                    ((value >> 16) & 0xFF) as u8,
                    ((value >> 8) & 0xFF) as u8,
                    (value & 0xFF) as u8,
                ));
            }
        }
    }
    Some((0, 0, 0))
}

/// `createExons`.
///
/// The exon numbering runs backwards on the negative strand, and the loop is skipped entirely
/// unless the two comma-separated lists have the same length, so a line with mismatched counts
/// yields a feature with no exons rather than an error.
fn create_exons(
    start: i32,
    tokens: &[&str],
    feature: &mut BedFeature,
    start_offset: StartOffset,
) -> Result<(), NumberFormat> {
    let cd_start = parse_int(tokens[6])? + start_offset.value();
    let cd_end = parse_int(tokens[7])?;
    let exon_count = parse_int(tokens[9])?;

    let sizes = split_fixed(tokens[10], exon_count as usize);
    let starts = split_fixed(tokens[11], exon_count as usize);
    if starts.len() != sizes.len() {
        return Ok(());
    }

    let negative = feature.strand == Some(Strand::Negative);
    let mut number = if negative { exon_count } else { 1 };
    for (offset, size) in starts.iter().zip(&sizes) {
        let exon_start = start + parse_int(offset)?;
        let exon_end = exon_start + parse_int(size)? - 1;
        feature.exons.push(Exon {
            start: exon_start,
            end: exon_end,
            cd_start,
            cd_end,
            number,
        });
        if negative {
            number -= 1;
        } else {
            number += 1;
        }
    }
    Ok(())
}

/// `ParsingUtils.split(string, array, delim)`: fills a fixed-size array, so a list longer than the
/// declared exon count is truncated and a shorter one leaves nulls the caller then parses.
///
/// The port returns only the fields that exist, and the caller's length test then behaves as the
/// reference's does for the mismatched case.
fn split_fixed(text: &str, expected: usize) -> Vec<&str> {
    text.split(',').take(expected).collect()
}

/// `canDecode(path)`: a block-compressed extension is stripped first, then the name is lowercased
/// and tested for `.bed`.
///
/// So `regions.BED.gz` decodes and `regions.bed.gz.gz` does not, because only one extension comes
/// off.
pub fn can_decode(path: &str) -> bool {
    let to_decode = if has_block_compressed_extension(path) {
        match path.rfind('.') {
            Some(index) => &path[..index],
            None => path,
        }
    } else {
        path
    };
    to_decode.to_lowercase().ends_with(".bed")
}

/// `IOUtil.hasBlockCompressedExtension`.
fn has_block_compressed_extension(path: &str) -> bool {
    let lowered = path.to_lowercase();
    lowered.ends_with(".gz")
        || lowered.ends_with(".gzip")
        || lowered.ends_with(".bgz")
        || lowered.ends_with(".bgzf")
}
