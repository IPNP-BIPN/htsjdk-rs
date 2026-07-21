//! Whole VCF files: header plus data lines.
//!
//! Ported from `htsjdk.variant.variantcontext.writer.VCFWriter` at htsjdk 4.2.0.
//!
//! The layer is thin, and its two decisions are both invisible in the format specification.
//!
//! **The version line is a constant, not a property of the header.** `VCFWriter` writes
//! `VERSION_LINE`, built from `VCFHeaderVersion.VCF4_2`, and *skips* any `fileformat` line the
//! header carries:
//!
//! ```java
//! for (final VCFHeaderLine line : header.getMetaDataInSortedOrder()) {
//!     if (VCFHeaderVersion.isFormatString(line.getKey())) continue;
//! ```
//!
//! So a header holding `##fileformat=VCFv4.3` writes `##fileformat=VCFv4.2` anyway. The writer
//! declares the version it can produce, not the version it was handed.
//!
//! **Every record is followed by a newline, and the header supplies its own.** There is no
//! separator between the two, so a writer that added one produces a file with a blank line that
//! most readers tolerate and no byte comparison accepts. Same shape as the SAM file layer.

use crate::encoder::{EncodeError, VcfEncoder};
use crate::header::VcfHeader;
use crate::variant::VariantContext;

/// `VCFWriter.VERSION_LINE`.
pub const VERSION_LINE: &str = "##fileformat=VCFv4.2";

/// Writes a complete VCF file.
pub fn write_vcf(header: &VcfHeader, records: &[VariantContext]) -> Result<String, EncodeError> {
    let encoder = VcfEncoder::new(header);
    let mut out = header.write();
    for record in records {
        out.push_str(&encoder.encode(record)?);
        out.push('\n');
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele::Allele;
    use crate::header::{Cardinality, HeaderLine, LineType};
    use crate::variant::{Genotype, Value};

    fn header() -> VcfHeader {
        let mut h = VcfHeader::new();
        h.lines.push(HeaderLine::info(
            "DP",
            Cardinality::Fixed(1),
            LineType::Integer,
            "Total Depth",
        ));
        h.lines.push(HeaderLine::format(
            "GT",
            Cardinality::Fixed(1),
            LineType::String,
            "Genotype",
        ));
        h.lines.push(HeaderLine::contig("chr1", 1000, 0));
        h.samples = vec!["s1".to_string()];
        h
    }

    fn record() -> VariantContext {
        let mut vc = VariantContext::new(
            "chr1",
            100,
            vec![
                Allele::from_str("A", true).unwrap(),
                Allele::from_str("T", false).unwrap(),
            ],
        );
        vc.attributes = vec![("DP".to_string(), Value::Int(42))];
        vc.genotypes = vec![Genotype::new(
            "s1",
            vec![
                Allele::from_str("A", true).unwrap(),
                Allele::from_str("T", false).unwrap(),
            ],
        )];
        vc
    }

    #[test]
    fn a_file_is_the_header_then_one_line_per_record() {
        let h = header();
        let text = write_vcf(&h, &[record(), record()]).unwrap();
        assert!(text.starts_with(VERSION_LINE));
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 7, "4 metadata + 1 column line + 2 records");
        assert_eq!(lines[5], "chr1\t100\t.\tA\tT\t.\t.\tDP=42\tGT\t0/1");
        assert_eq!(lines[5], lines[6]);
    }

    /// No blank line between the header and the first record: the header's own trailing newline
    /// is the only one.
    #[test]
    fn the_first_record_follows_the_header_directly() {
        let h = header();
        let text = write_vcf(&h, &[record()]).unwrap();
        let header_len = h.write().len();
        assert_eq!(&text[..header_len], h.write());
        assert!(
            text[header_len..].starts_with("chr1\t"),
            "a blank line crept in: {:?}",
            &text[header_len..header_len + 10]
        );
    }

    #[test]
    fn the_file_ends_with_a_newline() {
        let text = write_vcf(&header(), &[record()]).unwrap();
        assert!(text.ends_with('\n'));
        assert!(!text.ends_with("\n\n"), "no trailing blank line");
    }

    #[test]
    fn a_header_only_file_has_no_data_lines() {
        let h = header();
        let text = write_vcf(&h, &[]).unwrap();
        assert_eq!(text, h.write());
    }

    /// The version line is what the writer can produce, not what the header says. A header
    /// carrying a different `fileformat` is overwritten and its own line is dropped.
    #[test]
    fn the_version_line_is_the_writers_not_the_headers() {
        let mut h = header();
        h.lines.push(HeaderLine::Unstructured {
            key: "fileformat".to_string(),
            value: "VCFv4.3".to_string(),
        });
        let text = write_vcf(&h, &[]).unwrap();
        assert!(text.starts_with(VERSION_LINE));
        assert!(
            !text.contains("VCFv4.3"),
            "the header's own fileformat line must be skipped, not written twice"
        );
    }
}
