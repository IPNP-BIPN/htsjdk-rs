//! Byte-level port of htsjdk's BAM record codec.
//!
//! Ported from htsjdk 4.2.0, symbol by symbol, from the pinned clone. Every module names the
//! Java class and method it comes from.
//!
//! The BAM *format* is public and short. Reimplementing it from the specification produces a
//! file that `samtools` reads, that carries the right reads, and that is not htsjdk's file.
//! The gap between those two is what this crate is for, and it is made of choices the
//! specification explicitly leaves open:
//!
//! - the indexing `bin`, which readers scanning linearly never check ([`bin`]);
//! - the width chosen for each integer tag, which is picked from the value by a ladder that is
//!   not the obvious one ([`tag`]);
//! - the order tags are written in, which sorts on the tag's *second* character ([`tag`]);
//! - the trailing nibble of an odd-length read, which is `=` and not `N` ([`bases`]);
//! - what an absent quality string becomes ([`record`]);
//! - where a CIGAR too long for the format goes, and what stands in for it ([`record`]).

pub mod alignment_block;
pub mod bases;
pub mod bin;
pub mod build_index;
pub mod cigar;
pub mod coordinate;
pub mod fasta;
pub mod fastq;
pub mod gather;
pub mod header;
pub mod index;
pub mod interval;
pub mod md_nm;
pub mod murmur3;
pub mod overlap;
pub mod pair;
pub mod query_name;
pub mod read_group_checksum;
pub mod reader;
pub mod record;
pub mod reheader;
pub mod sam_file;
pub mod sequence;
pub mod tag;
pub mod text;
pub mod text_parse;
pub mod writer;

pub use bin::compute_indexing_bin;
pub use build_index::{build_bam_index, BuildIndexError};
pub use cigar::{Cigar, CigarElement, Op};
pub use gather::gather_bam_files;
pub use header::SamHeader;
pub use reader::BamReader;
pub use record::{BamRecord, DecodeError, EncodeError};
pub use reheader::{block_copy, reheader_bam, ReheaderError};
pub use tag::{Tag, TagValue, Tags};
pub use writer::write_bam_header_block;
pub use writer::BamWriter;
