//! BGZF virtual file pointers.
//!
//! Ported from htsjdk 4.2.0
//! `src/main/java/htsjdk/samtools/util/BlockCompressedFilePointerUtil.java`.
//!
//! A virtual file pointer packs the byte offset of a block's first byte in the compressed
//! stream (48 bits) with an offset inside that block's uncompressed data (16 bits).

const SHIFT_AMOUNT: u32 = 16;
const OFFSET_MASK: u64 = 0xffff;
const ADDRESS_MASK: u64 = 0xFFFF_FFFF_FFFF;

pub const MAX_BLOCK_ADDRESS: u64 = ADDRESS_MASK;
pub const MAX_OFFSET: u16 = u16::MAX;

/// Error from constructing an out-of-range virtual file pointer.
#[derive(Debug, PartialEq, Eq)]
pub enum VfpError {
    BlockAddressTooLarge(u64),
    BlockOffsetTooLarge(u32),
}

impl std::fmt::Display for VfpError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::BlockAddressTooLarge(a) => write!(f, "blockAddress {a} too large."),
            Self::BlockOffsetTooLarge(o) => write!(f, "blockOffset {o} too large."),
        }
    }
}

impl std::error::Error for VfpError {}

/// Packs a block address and an in-block offset into a virtual file pointer.
///
/// Java takes signed values here and rejects negatives; the Rust signature makes that
/// unrepresentable, so only the upper bounds remain as runtime checks.
pub fn make_file_pointer(block_address: u64, block_offset: u32) -> Result<u64, VfpError> {
    if block_offset > MAX_OFFSET as u32 {
        return Err(VfpError::BlockOffsetTooLarge(block_offset));
    }
    if block_address > MAX_BLOCK_ADDRESS {
        return Err(VfpError::BlockAddressTooLarge(block_address));
    }
    Ok(block_address << SHIFT_AMOUNT | block_offset as u64)
}

pub fn block_address(vfp: u64) -> u64 {
    (vfp >> SHIFT_AMOUNT) & ADDRESS_MASK
}

pub fn block_offset(vfp: u64) -> u16 {
    (vfp & OFFSET_MASK) as u16
}

/// True when the two pointers are in the same block or in consecutive blocks.
pub fn are_in_same_or_adjacent_blocks(vfp1: u64, vfp2: u64) -> bool {
    let b1 = block_address(vfp1);
    let b2 = block_address(vfp2);
    b1 == b2 || b1 + 1 == b2
}

/// Moves a pointer by whole blocks, keeping the in-block offset.
pub fn shift(vfp: u64, offset: u64) -> Result<u64, VfpError> {
    make_file_pointer(block_address(vfp) + offset, block_offset(vfp) as u32)
}

/// Renders a pointer the way htsjdk's `asString` does.
pub fn as_string(vfp: u64) -> String {
    format!(
        "{}(0x{:x}): (block address: {}, offset: {})",
        vfp,
        vfp,
        block_address(vfp),
        block_offset(vfp)
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_address_and_offset() {
        for &(addr, off) in &[(0u64, 0u32), (1, 1), (12345, 65535), (MAX_BLOCK_ADDRESS, 0)] {
            let vfp = make_file_pointer(addr, off).unwrap();
            assert_eq!(block_address(vfp), addr);
            assert_eq!(block_offset(vfp) as u32, off);
        }
    }

    #[test]
    fn rejects_out_of_range() {
        assert_eq!(
            make_file_pointer(0, 65536),
            Err(VfpError::BlockOffsetTooLarge(65536))
        );
        assert_eq!(
            make_file_pointer(MAX_BLOCK_ADDRESS + 1, 0),
            Err(VfpError::BlockAddressTooLarge(MAX_BLOCK_ADDRESS + 1))
        );
    }

    #[test]
    fn adjacency() {
        let a = make_file_pointer(10, 5).unwrap();
        assert!(are_in_same_or_adjacent_blocks(
            a,
            make_file_pointer(10, 900).unwrap()
        ));
        assert!(are_in_same_or_adjacent_blocks(
            a,
            make_file_pointer(11, 0).unwrap()
        ));
        assert!(!are_in_same_or_adjacent_blocks(
            a,
            make_file_pointer(12, 0).unwrap()
        ));
    }

    /// Matches the format string in htsjdk's `asString`.
    #[test]
    fn as_string_format() {
        let vfp = make_file_pointer(3, 7).unwrap();
        assert_eq!(
            as_string(vfp),
            "196615(0x30007): (block address: 3, offset: 7)"
        );
    }
}
