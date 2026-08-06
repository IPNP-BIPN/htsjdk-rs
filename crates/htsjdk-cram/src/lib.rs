//! Port of htsjdk's CRAM support, from `htsjdk.samtools.cram` (htsjdk 4.2.0).
//!
//! Milestone H.3, a sub-project on its own: 169 Java files covering a container model, a dozen
//! encodings, codec negotiation and reference-based compression. It is built from the bottom, and
//! [`varint`] is the bottom: every structure above it is a run of ITF8s, so nothing above can be
//! checked until those are.
pub mod bit_stream;
pub mod block;
pub mod compression_header;
pub mod container;
pub mod core_codecs;
pub mod encoding_factory;
pub mod encoding_map;
pub mod external_codecs;
pub mod golomb;
pub mod huffman;
pub mod preservation_map;
pub mod rans;
pub mod rans_order1;
pub mod read_features;
pub mod record_flags;
pub mod slice_blocks;
pub mod slice_header;
pub mod substitution_matrix;
pub mod tag_encoding_map;
pub mod varint;
