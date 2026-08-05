//! Port of htsjdk's CRAM support, from `htsjdk.samtools.cram` (htsjdk 4.2.0).
//!
//! Milestone H.3, a sub-project on its own: 169 Java files covering a container model, a dozen
//! encodings, codec negotiation and reference-based compression. It is built from the bottom, and
//! [`varint`] is the bottom: every structure above it is a run of ITF8s, so nothing above can be
//! checked until those are.
pub mod block;
pub mod container;
pub mod encoding_map;
pub mod preservation_map;
pub mod rans;
pub mod rans_order1;
pub mod slice_header;
pub mod substitution_matrix;
pub mod tag_encoding_map;
pub mod varint;
