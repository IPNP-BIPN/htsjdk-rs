//! The CRAM slice header: the last frame before CRAM becomes reads.
//!
//! Ported from `htsjdk.samtools.cram.structure.Slice`'s stream constructor and
//! `createSliceHeaderBlockContent` at htsjdk 4.2.0.
//!
//! The compression header's three maps are pinned. A container's blocks then hold one or more
//! slices, and each slice begins with a `MAPPED_SLICE` block whose raw content is this header.
//!
//! # The block count does not count the header block
//!
//! `getNumberOfBlocks` returns `1 + numberOfExternalBlocks`: the core block plus the externals.
//! Measured on five files, it equals **exactly** the number of blocks that follow the header. A
//! reader that counts the header block among them reads one block too few and stops before the
//! last one.
//!
//! # Six tags ride in the header, and four of them digest nothing
//!
//! Measured: `B1` and `S1` carry a SHA-1, `B5` and `S5` a SHA-512, and on an unmapped slice all
//! four are the digest of the **empty string**, byte for byte identical in every file. That is 168
//! bytes of constant per slice. Only `BD` and `SD`, four bytes apiece, move with the reads, and
//! they do not move when only the tags change.
//!
//! # The tag section has no length
//!
//! It is read with `readFully` to the end of the block, so the slice header block's own length is
//! the only thing that delimits it. A header with no tags is indistinguishable from one whose tags
//! are zero bytes long. This module therefore carries the section as **opaque bytes**: reproducing
//! it is what byte-identity needs, and decoding it is `BinaryTagCodec`, a slice of its own.
//!
//! # The MD5 is sixteen zeroes when there is none
//!
//! `createSliceHeaderBlockContent` writes `new byte[16]` when the reference MD5 is null, so the
//! field is always present and its emptiness is a value rather than an absence.
//!
//! # An absent embedded reference is -1
//!
//! Written as an ITF8, which puts it in the **five-byte** form. The commonest value of this field
//! is also its longest encoding.
//!
//! # The alignment context carries magic numbers
//!
//! A reference id of -1 is unmapped-unplaced and -2 is multiple-reference, and both force the
//! start and the span to 0.

use crate::varint::{
    read_unsigned_itf8, read_unsigned_ltf8, write_unsigned_itf8, write_unsigned_ltf8, RuntimeEof,
};

/// `Slice.MD5_BYTE_SIZE`.
pub const MD5_BYTE_SIZE: usize = 16;
/// `Slice.EMBEDDED_REFERENCE_ABSENT_CONTENT_ID`.
pub const EMBEDDED_REFERENCE_ABSENT: i32 = -1;
/// `AlignmentContext.NO_ALIGNMENT_START` and `NO_ALIGNMENT_SPAN`, which are the same value.
pub const NO_ALIGNMENT: i32 = 0;
/// `ReferenceContext`: unmapped and unplaced.
pub const UNMAPPED_UNPLACED_ID: i32 = -1;
/// `ReferenceContext`: more than one reference in this slice.
pub const MULTIPLE_REFERENCE_ID: i32 = -2;

/// One slice header, as its own block's raw content.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SliceHeader {
    /// A sequence index, or [`UNMAPPED_UNPLACED_ID`], or [`MULTIPLE_REFERENCE_ID`].
    pub reference_context_id: i32,
    pub alignment_start: i32,
    pub alignment_span: i32,
    pub record_count: i32,
    pub global_record_counter: i64,
    /// The core block plus the externals, and **not** the header block.
    pub block_count: i32,
    /// The external content ids this slice uses, which are the data series ids plus one per tag.
    pub content_ids: Vec<i32>,
    /// [`EMBEDDED_REFERENCE_ABSENT`] when there is none.
    pub embedded_reference_content_id: i32,
    /// Always present, sixteen zeroes when there is no reference.
    pub reference_md5: [u8; MD5_BYTE_SIZE],
    /// The BAM-encoded tags, carried opaquely: the section has no length of its own and decoding
    /// it is a separate slice.
    pub tags: Vec<u8>,
}

impl SliceHeader {
    /// Is this slice unmapped and unplaced?
    pub fn is_unmapped_unplaced(&self) -> bool {
        self.reference_context_id == UNMAPPED_UNPLACED_ID
    }

    /// Does it span more than one reference?
    pub fn is_multiple_reference(&self) -> bool {
        self.reference_context_id == MULTIPLE_REFERENCE_ID
    }

    /// Read the header from the raw content of a `MAPPED_SLICE` block.
    ///
    /// The tag section is whatever is left, which is why this takes the block's content rather
    /// than a stream: nothing inside the header says where it ends.
    pub fn read(content: &[u8]) -> Result<Self, RuntimeEof> {
        let mut at = 0usize;
        let itf8 = |at: &mut usize| -> Result<i32, RuntimeEof> {
            let (value, consumed) = read_unsigned_itf8(&content[(*at).min(content.len())..])?;
            *at += consumed;
            Ok(value)
        };

        let reference_context_id = itf8(&mut at)?;
        let alignment_start = itf8(&mut at)?;
        let alignment_span = itf8(&mut at)?;
        let record_count = itf8(&mut at)?;
        let (global_record_counter, consumed) =
            read_unsigned_ltf8(&content[at.min(content.len())..])?;
        at += consumed;
        let block_count = itf8(&mut at)?;

        // `CramIntArray`: a count, then that many values.
        let content_id_count = itf8(&mut at)?;
        let mut content_ids = Vec::with_capacity(content_id_count.max(0) as usize);
        for _ in 0..content_id_count.max(0) {
            content_ids.push(itf8(&mut at)?);
        }
        let embedded_reference_content_id = itf8(&mut at)?;

        let reference_md5: [u8; MD5_BYTE_SIZE] = content
            .get(at..at + MD5_BYTE_SIZE)
            .ok_or(RuntimeEof)?
            .try_into()
            .expect("sixteen bytes");
        at += MD5_BYTE_SIZE;

        Ok(Self {
            reference_context_id,
            alignment_start,
            alignment_span,
            record_count,
            global_record_counter,
            block_count,
            content_ids,
            embedded_reference_content_id,
            reference_md5,
            tags: content[at.min(content.len())..].to_vec(),
        })
    }

    /// The header's own bytes, which are the raw content of its block.
    ///
    /// `major` decides whether the tag section is written at all: below version 3 it is neither
    /// written nor read, so the same slice is shorter in a 2.1 file by however many bytes its tags
    /// would have occupied.
    pub fn write(&self, major: u8) -> Vec<u8> {
        let mut out = Vec::new();
        push_itf8(self.reference_context_id, &mut out);
        push_itf8(self.alignment_start, &mut out);
        push_itf8(self.alignment_span, &mut out);
        push_itf8(self.record_count, &mut out);
        out.extend_from_slice(&write_unsigned_ltf8(self.global_record_counter).0);
        push_itf8(self.block_count, &mut out);

        push_itf8(self.content_ids.len() as i32, &mut out);
        for id in &self.content_ids {
            push_itf8(*id, &mut out);
        }
        push_itf8(self.embedded_reference_content_id, &mut out);
        out.extend_from_slice(&self.reference_md5);

        if major >= 3 {
            out.extend_from_slice(&self.tags);
        }
        out
    }
}

fn push_itf8(value: i32, out: &mut Vec<u8>) {
    out.extend_from_slice(&write_unsigned_itf8(value).0);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unmapped() -> SliceHeader {
        SliceHeader {
            reference_context_id: UNMAPPED_UNPLACED_ID,
            alignment_start: NO_ALIGNMENT,
            alignment_span: NO_ALIGNMENT,
            record_count: 4,
            global_record_counter: 0,
            block_count: 27,
            content_ids: (1..=26).collect(),
            embedded_reference_content_id: EMBEDDED_REFERENCE_ABSENT,
            reference_md5: [0u8; MD5_BYTE_SIZE],
            tags: vec![b'B', b'D', b'B', b'c', 4, 0, 0, 0, 1, 2, 3, 4],
        }
    }

    #[test]
    fn a_header_round_trips() {
        let header = unmapped();
        assert_eq!(SliceHeader::read(&header.write(3)), Ok(header));
    }

    /// Below version 3 the tags are not written, so the same header is shorter and reads back
    /// without them.
    #[test]
    fn the_tag_section_is_version_gated() {
        let header = unmapped();
        let three = header.write(3);
        let two = header.write(2);
        assert_eq!(three.len(), two.len() + header.tags.len());
        assert!(SliceHeader::read(&two).expect("parses").tags.is_empty());
    }

    /// The commonest value of the embedded reference field is also its longest encoding.
    #[test]
    fn an_absent_embedded_reference_takes_five_bytes() {
        let (bytes, _) = write_unsigned_itf8(EMBEDDED_REFERENCE_ABSENT);
        assert_eq!(bytes.len(), 5);
    }

    /// The MD5 is a value, not an absence: sixteen zeroes are still sixteen bytes.
    #[test]
    fn the_md5_is_always_sixteen_bytes() {
        let header = unmapped();
        assert_eq!(header.reference_md5.len(), MD5_BYTE_SIZE);
        let written = header.write(3);
        // It sits between the embedded reference id and the tags.
        let read = SliceHeader::read(&written).expect("parses");
        assert_eq!(read.reference_md5, [0u8; MD5_BYTE_SIZE]);
    }

    #[test]
    fn the_two_magic_reference_contexts_are_recognised() {
        let mut header = unmapped();
        assert!(header.is_unmapped_unplaced());
        assert!(!header.is_multiple_reference());
        header.reference_context_id = MULTIPLE_REFERENCE_ID;
        assert!(header.is_multiple_reference());
    }
}
