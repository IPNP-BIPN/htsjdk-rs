//! Port of htsjdk's CRAM support, from `htsjdk.samtools.cram`, measured against htsjdk 5.0.0. The
//! 3.0 codecs were transcribed from 4.2.0, CRAM 3.1 and the writer from 5.0.0 (decision 0047).
//!
//! Milestone H.3, a sub-project on its own: 169 Java files covering a container model, a dozen
//! encodings, codec negotiation and reference-based compression. It is built from the bottom, and
//! [`varint`] is the bottom: every structure above it is a run of ITF8s, so nothing above can be
//! checked until those are.
pub mod bai;
pub mod bit_stream;
pub mod block;
pub mod codecs;
pub mod compression_header;
pub mod container;
pub mod core_codecs;
pub mod crai;
pub mod encoding_factory;
pub mod encoding_map;
pub mod external_codecs;
pub mod file;
pub mod fqzcomp;
pub mod golomb;
pub mod gzip;
pub mod huffman;
pub mod name_tokeniser;
pub mod negotiation;
pub mod preservation_map;
pub mod range;
pub mod rans;
pub mod rans_nx16;
pub mod rans_order1;
pub mod read_features;
pub mod reader;
pub mod record_flags;
pub mod record_read;
pub mod record_write;
pub mod slice_blocks;
pub mod slice_header;
pub mod substitution_matrix;
pub mod tag_encoding_map;
pub mod varint;
pub mod writer;
