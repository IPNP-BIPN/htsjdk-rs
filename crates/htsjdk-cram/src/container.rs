//! The CRAM file definition and the container header: the first structures built out of
//! [`crate::varint`].
//!
//! Ported from `htsjdk.samtools.cram.structure.CramHeader`,
//! `htsjdk.samtools.cram.structure.ContainerHeader` and `htsjdk.samtools.cram.io.CramInt` at
//! htsjdk 4.2.0.
//!
//! A CRAM file is a 26-byte definition, then containers. The definition is fixed-width; the
//! container header is not, and it mixes little-endian `int32`s with the ITF8s and LTF8s underneath.
//!
//! # The file id is padded to exactly 20 bytes and truncated in silence
//!
//! `CramHeader` fills the array with zeros and copies `min(id.length(), 20)`. A shorter id is
//! zero-padded rather than terminated, and a longer one **loses its tail with no error**. Measured:
//! `far-too-long-to-fit-in-twenty` is written as `far-too-long-to-fit-`.
//!
//! # The checksum covers the header, not the container, and is little-endian
//!
//! The CRC32 is computed over the container header's own bytes, up to but not including the
//! checksum, and written **little-endian** — the opposite of the CRC in a BGZF block's gzip
//! trailer, which is the neighbouring format in this repository.
//!
//! And **it is absent below version 3**: `cramVersion.getMajor() >= 3` guards both the read and the
//! write, so the same container is four bytes shorter in a 2.1 file and a reader that always
//! consumes them is off by four from the first container onwards.
//!
//! # Every file has at least two containers, and the last one is a magic number
//!
//! Measured on four files: a SAM header container, then a data container per block of records, then
//! an **EOF container**. A file with no records has only the first and the last.
//!
//! The EOF container's `alignmentStart` is **4542278**, which is `0x454F46`, which is `EOF` in
//! ASCII. A magic number in a coordinate field, and a port that validated coordinates would refuse
//! the end of every well-formed CRAM.

use crate::varint::{read_unsigned_itf8, read_unsigned_ltf8, RuntimeEof};

/// `CramHeader.MAGIC`.
pub const MAGIC: &[u8; 4] = b"CRAM";
/// `CramHeader.CRAM_ID_LENGTH`.
pub const ID_LENGTH: usize = 20;
/// `CramHeader.CRAM_HEADER_LENGTH`: magic, two version bytes, then the id.
pub const FILE_DEFINITION_LENGTH: usize = 4 + 2 + ID_LENGTH;

/// The `alignmentStart` of the EOF container, which is `EOF` read as a big-endian integer.
pub const EOF_ALIGNMENT_START: i32 = 0x0045_4F46;

/// `CramHeader`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileDefinition {
    pub major: u8,
    pub minor: u8,
    /// Exactly [`ID_LENGTH`] bytes, zero-padded.
    pub id: [u8; ID_LENGTH],
}

/// What a malformed file definition is refused with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContainerError {
    /// The first four bytes are not `CRAM`.
    BadMagic { found: [u8; 4] },
    /// Fewer than [`FILE_DEFINITION_LENGTH`] bytes.
    Truncated,
    /// `CramVersions.isSupportedVersion` answers false. Measured: 2.1 and 3.0 are true, 3.1 and
    /// 4.0 are false, and the message is the one htsjdk throws.
    UnsupportedVersion { major: u8, minor: u8 },
}

impl ContainerError {
    pub fn message(&self) -> String {
        match self {
            ContainerError::BadMagic { .. } => "not a CRAM file".to_string(),
            ContainerError::Truncated => "the file definition is truncated".to_string(),
            ContainerError::UnsupportedVersion { major, minor } => {
                format!("CRAM version {major}.{minor} is not supported")
            }
        }
    }
}

/// `CramVersions.isSupportedVersion`, measured in decision 0038 rather than read from a list.
pub fn is_supported_version(major: u8, minor: u8) -> bool {
    matches!((major, minor), (2, 1) | (3, 0))
}

impl FileDefinition {
    /// `CramHeader`'s constructor: pad to 20 bytes, truncate without complaint.
    pub fn new(major: u8, minor: u8, id: &str) -> Self {
        let mut padded = [0u8; ID_LENGTH];
        let bytes = id.as_bytes();
        let take = bytes.len().min(ID_LENGTH);
        padded[..take].copy_from_slice(&bytes[..take]);
        Self {
            major,
            minor,
            id: padded,
        }
    }

    /// Read the 26 bytes at the head of a CRAM file.
    ///
    /// The version is **not** checked here: htsjdk reads the definition and refuses the version
    /// later, which is why a 3.1 file gets past this and fails with a version message rather than
    /// a magic one. [`is_supported_version`] is the separate test.
    pub fn read(bytes: &[u8]) -> Result<Self, ContainerError> {
        if bytes.len() < FILE_DEFINITION_LENGTH {
            return Err(ContainerError::Truncated);
        }
        let magic: [u8; 4] = bytes[..4].try_into().expect("four bytes");
        if &magic != MAGIC {
            return Err(ContainerError::BadMagic { found: magic });
        }
        let mut id = [0u8; ID_LENGTH];
        id.copy_from_slice(&bytes[6..FILE_DEFINITION_LENGTH]);
        Ok(Self {
            major: bytes[4],
            minor: bytes[5],
            id,
        })
    }

    pub fn write(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(FILE_DEFINITION_LENGTH);
        out.extend_from_slice(MAGIC);
        out.push(self.major);
        out.push(self.minor);
        out.extend_from_slice(&self.id);
        out
    }
}

/// `ContainerHeader`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContainerHeader {
    /// The size of the blocks that follow, **not** including this header.
    pub blocks_byte_size: i32,
    /// `ReferenceContext`: a sequence index, or -1 for unmapped-unplaced, or -2 for multiple.
    pub reference_context_id: i32,
    pub alignment_start: i32,
    pub alignment_span: i32,
    pub record_count: i32,
    pub global_record_counter: i64,
    pub base_count: i64,
    pub block_count: i32,
    /// Byte offsets of the slices within the container's blocks.
    pub landmarks: Vec<i32>,
    /// Zero below version 3, where the field is not written at all.
    pub checksum: i32,
    /// How many bytes this header occupied, which the caller needs to find the blocks.
    pub byte_length: usize,
}

impl ContainerHeader {
    /// Is this the container that ends the file?
    ///
    /// The test is the magic in the coordinate field. htsjdk compares whole containers against a
    /// constant; this compares the field that carries the constant's signature, which is the same
    /// answer on every file either produces.
    pub fn is_eof(&self) -> bool {
        self.alignment_start == EOF_ALIGNMENT_START && self.record_count == 0
    }

    /// `ContainerHeader(CRAMVersion, InputStream)`.
    ///
    /// `major` decides whether the checksum is there. Passing the wrong one does not fail: it reads
    /// four bytes of the next block as a checksum, or leaves four bytes of checksum to be read as a
    /// block. That is the shape of every version-dependent field and the reason this takes the
    /// version rather than guessing.
    pub fn read(bytes: &[u8], major: u8) -> Result<Self, RuntimeEof> {
        let mut at = 0usize;

        let blocks_byte_size = read_int32(bytes, &mut at);
        let reference_context_id = itf8(bytes, &mut at)?;
        let alignment_start = itf8(bytes, &mut at)?;
        let alignment_span = itf8(bytes, &mut at)?;
        let record_count = itf8(bytes, &mut at)?;
        let global_record_counter = ltf8(bytes, &mut at)?;
        let base_count = ltf8(bytes, &mut at)?;
        let block_count = itf8(bytes, &mut at)?;

        // `CramIntArray`: a count, then that many values, all ITF8.
        let landmark_count = itf8(bytes, &mut at)?;
        let mut landmarks = Vec::with_capacity(landmark_count.max(0) as usize);
        for _ in 0..landmark_count.max(0) {
            landmarks.push(itf8(bytes, &mut at)?);
        }

        let checksum = if major >= 3 {
            read_int32(bytes, &mut at)
        } else {
            0
        };

        Ok(Self {
            blocks_byte_size,
            reference_context_id,
            alignment_start,
            alignment_span,
            record_count,
            global_record_counter,
            base_count,
            block_count,
            landmarks,
            checksum,
            byte_length: at,
        })
    }
}

/// `CramInt.readInt32`: four bytes, **little-endian**, and each one read as a signed byte shifted
/// into place, so a short slice yields whatever the missing bytes read as rather than an error.
fn read_int32(bytes: &[u8], at: &mut usize) -> i32 {
    let byte = |offset: usize| -> i32 {
        match bytes.get(*at + offset) {
            Some(value) => i32::from(*value),
            // `InputStream.read()` past the end, which the arithmetic uses as it stands.
            None => -1,
        }
    };
    let value = byte(0) | (byte(1) << 8) | (byte(2) << 16) | (byte(3) << 24);
    *at += 4;
    value
}

fn itf8(bytes: &[u8], at: &mut usize) -> Result<i32, RuntimeEof> {
    let (value, consumed) = read_unsigned_itf8(&bytes[(*at).min(bytes.len())..])?;
    *at += consumed;
    Ok(value)
}

fn ltf8(bytes: &[u8], at: &mut usize) -> Result<i64, RuntimeEof> {
    let (value, consumed) = read_unsigned_ltf8(&bytes[(*at).min(bytes.len())..])?;
    *at += consumed;
    Ok(value)
}

/// Walk a whole CRAM's container headers, skipping the blocks between them.
pub fn container_headers(
    cram: &[u8],
) -> Result<(FileDefinition, Vec<ContainerHeader>), ContainerError> {
    let definition = FileDefinition::read(cram)?;
    let mut headers = Vec::new();
    let mut at = FILE_DEFINITION_LENGTH;
    while at < cram.len() {
        let Ok(header) = ContainerHeader::read(&cram[at..], definition.major) else {
            break;
        };
        at += header.byte_length + header.blocks_byte_size.max(0) as usize;
        let eof = header.is_eof();
        headers.push(header);
        if eof {
            break;
        }
    }
    Ok((definition, headers))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_id_is_padded_and_truncated_without_complaint() {
        assert_eq!(FileDefinition::new(3, 0, "").id, [0u8; ID_LENGTH]);
        let short = FileDefinition::new(3, 0, "short");
        assert_eq!(&short.id[..5], b"short");
        assert_eq!(&short.id[5..], &[0u8; 15]);
        let long = FileDefinition::new(3, 0, "far-too-long-to-fit-in-twenty");
        assert_eq!(&long.id[..], b"far-too-long-to-fit-");
    }

    #[test]
    fn a_file_definition_round_trips() {
        let definition = FileDefinition::new(3, 0, "abc");
        assert_eq!(FileDefinition::read(&definition.write()), Ok(definition));
    }

    #[test]
    fn the_magic_is_checked_and_the_version_is_not() {
        assert!(matches!(
            FileDefinition::read(b"NOPE\x03\x00                    "),
            Err(ContainerError::BadMagic { .. })
        ));
        // 3.1 reads fine here; it is refused by the version test, separately.
        let three_one = FileDefinition::read(b"CRAM\x03\x01                    ").expect("reads");
        assert_eq!((three_one.major, three_one.minor), (3, 1));
        assert!(!is_supported_version(3, 1));
        assert!(is_supported_version(3, 0));
        assert!(is_supported_version(2, 1));
    }

    /// The magic hidden in a coordinate.
    #[test]
    fn the_eof_container_is_recognised_by_a_number_that_spells_eof() {
        assert_eq!(EOF_ALIGNMENT_START, 0x0045_4F46);
        assert_eq!(&EOF_ALIGNMENT_START.to_be_bytes()[1..], b"EOF");
    }
}
