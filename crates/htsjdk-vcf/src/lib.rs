//! Port of htsjdk's VCF support.
//!
//! Ported from htsjdk 4.2.0 `htsjdk.variant.vcf` and
//! `htsjdk.variant.variantcontext.writer`.

pub mod allele;
pub mod chromosome_counts;
pub mod encoder;
pub mod genotype_likelihoods;
pub mod genotype_parse;
pub mod header;
pub mod header_lines;
pub mod header_parse;
pub mod jformat;
pub mod record_parse;
pub mod variant;
pub mod vcf_file;

pub use header::{Cardinality, HeaderLine, LineType, VcfHeader};
