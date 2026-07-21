//! Port of htsjdk's VCF support.
//!
//! Ported from htsjdk 4.2.0 `htsjdk.variant.vcf` and
//! `htsjdk.variant.variantcontext.writer`.

pub mod allele;
pub mod encoder;
pub mod header;
pub mod jformat;
pub mod variant;
pub mod vcf_file;

pub use header::{Cardinality, HeaderLine, LineType, VcfHeader};
