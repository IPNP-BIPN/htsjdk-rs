//! Port of htsjdk's VCF support.
//!
//! Ported from htsjdk 4.2.0 `htsjdk.variant.vcf` and
//! `htsjdk.variant.variantcontext.writer`.

pub mod header;

pub use header::{Cardinality, HeaderLine, LineType, VcfHeader};
