//! The version-dependent rewrite every attribute value goes through on the way in.
//!
//! Ported from `htsjdk.variant.vcf.VCFPercentEncodedTextTransformer` and
//! `AbstractVCFCodec.getTextTransformerForVCFVersion` at htsjdk 4.2.0.
//!
//! From VCF 4.3 an attribute value may carry percent-encoded characters, and htsjdk decodes them
//! on read. Below 4.3 it does not, so **the same data line means different things in two files
//! that differ only in their `##fileformat` line**, and nothing on the line itself says which.
//! The transformer is chosen once, when the header is read, and then applied to every INFO value
//! and every genotype value in the file.
//!
//! # The decoder is `Integer.parseInt(s, 16)`, and that is not the same as "two hex digits"
//!
//! `parseInt` accepts a sign, so `%+1` parses as 1 and decodes to U+0001, a character no
//! percent-encoding can name. `%-1` parses as -1, `Character.toChars(-1)` throws, the throw is
//! caught, and the text is left alone. Neither is in any specification; both are measured.
//!
//! The catch is the general shape here: **a sequence that fails to decode is emitted verbatim**
//! rather than refused, and the `%` is *not* skipped, so the next character is examined again as a
//! possible start. That is why `%%41` decodes to `%A`: the first `%` fails on `%4`, is written
//! out, and the second one succeeds on `41`.
//!
//! The guard is `(i + 2) < length`, so the last two characters of a string can never begin an
//! escape: `%4` at the end stays `%4`, and so does a trailing `%`.

/// Which transformer a version selects. `VCFPassThruTextTransformer` for anything below 4.3.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TextTransformer {
    /// Below VCF 4.3: the text is the text.
    PassThru,
    /// VCF 4.3 and later: percent escapes are decoded.
    PercentEncoded,
}

impl TextTransformer {
    /// `getTextTransformerForVCFVersion`.
    pub fn for_version(version: crate::header_parse::VcfVersion) -> Self {
        if version.is_at_least(crate::header_parse::VcfVersion::Vcf4_3) {
            TextTransformer::PercentEncoded
        } else {
            TextTransformer::PassThru
        }
    }

    /// `decodeText(String)`.
    pub fn decode(self, raw: &str) -> String {
        match self {
            TextTransformer::PassThru => raw.to_string(),
            TextTransformer::PercentEncoded => decode_percent_encoded_chars(raw),
        }
    }

    /// `decodeText(List<String>)`, which is the single-string one mapped over the list.
    pub fn decode_all<'a>(self, raw: impl Iterator<Item = &'a str>) -> Vec<String> {
        raw.map(|part| self.decode(part)).collect()
    }
}

/// `VCFPercentEncodedTextTransformer.decodePercentEncodedChars`.
///
/// Indexing is over UTF-16 code units upstream and over bytes here. The two agree because the
/// only inputs that decode are ASCII escapes, and a non-ASCII character is copied through
/// untouched by both; working in bytes keeps the `(i + 2) < length` guard on the same units the
/// Java one uses for the inputs that reach it.
pub fn decode_percent_encoded_chars(raw: &str) -> String {
    // The whole function is skipped when there is no '%' at all, which is not just an
    // optimisation: it is the reason a string with no escapes is returned identical rather than
    // rebuilt.
    if !raw.contains('%') {
        return raw.to_string();
    }

    let bytes = raw.as_bytes();
    let mut out = String::with_capacity(raw.len());
    let mut i = 0usize;
    while i < bytes.len() {
        let c = bytes[i];
        if c == b'%' && i + 2 < bytes.len() {
            // `get` rather than a slice: the two bytes after a '%' may land inside a multi-byte
            // character, where a slice panics. Upstream indexes UTF-16 units and always gets a
            // substring, which `parseInt` then refuses, so `None` is the same answer.
            match parse_int_radix16(raw.get(i + 1..i + 3).unwrap_or("")) {
                // `Character.toChars` refuses a negative value and anything above the last code
                // point, and the throw is caught, so the escape is emitted verbatim. Two hex
                // digits with a sign cannot exceed 0xff, so only the negative arm is reachable.
                // A lone surrogate would be `None` here: `Character.toChars` produces one char and
                // htsjdk keeps it, where Rust has no such `char`. Unreachable from two hex digits,
                // which cannot exceed 0xff.
                Some(value) if (0..=0x10_FFFF).contains(&value) => {
                    if let Some(decoded) = char::from_u32(value as u32) {
                        out.push(decoded);
                        i += 3;
                        continue;
                    }
                }
                _ => {}
            }
            // The '%' is written and **not** skipped over, so the next byte is examined as a
            // possible escape start in its own right.
            out.push('%');
            i += 1;
        } else {
            // Copy the whole UTF-8 sequence, not the byte: `raw` is a `str`, so a multi-byte
            // character copied byte by byte would still be valid, but pushing the char keeps the
            // loop honest about what a position means.
            let ch = raw[i..].chars().next().expect("i is a char boundary");
            out.push(ch);
            i += ch.len_utf8();
        }
    }
    out
}

/// `Integer.parseInt(s, 16)`, which is not `i64::from_str_radix`: it accepts a leading `+`.
///
/// Rust's own `from_str_radix` accepts `+` too, so this is a thin wrapper, but the difference is
/// worth naming because the *sign* is the whole reason `%+1` and `%-1` behave as they do, and a
/// port written from the specification would reject both.
fn parse_int_radix16(text: &str) -> Option<i64> {
    i64::from_str_radix(text, 16).ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header_parse::VcfVersion;

    #[test]
    fn the_transformer_is_chosen_by_version_and_the_cut_is_at_four_three() {
        assert_eq!(
            TextTransformer::for_version(VcfVersion::Vcf4_2),
            TextTransformer::PassThru
        );
        assert_eq!(
            TextTransformer::for_version(VcfVersion::Vcf4_3),
            TextTransformer::PercentEncoded
        );
    }

    #[test]
    fn a_failed_escape_is_emitted_and_the_percent_is_re_examined() {
        assert_eq!(decode_percent_encoded_chars("%%41"), "%A");
        assert_eq!(decode_percent_encoded_chars("%4G"), "%4G");
    }

    /// The two characters at the end of a string can never start an escape.
    #[test]
    fn the_guard_is_strictly_less_than() {
        assert_eq!(decode_percent_encoded_chars("%41"), "A");
        assert_eq!(decode_percent_encoded_chars("x%4"), "x%4");
        assert_eq!(decode_percent_encoded_chars("x%"), "x%");
    }

    #[test]
    fn the_sign_is_accepted_by_parse_int_and_that_decides_two_cases() {
        assert_eq!(decode_percent_encoded_chars("%+1"), "\u{1}");
        assert_eq!(decode_percent_encoded_chars("%-1"), "%-1");
    }
}
