//! `VCFCodec.canDecode`, ported from `htsjdk.variant.vcf.AbstractVCFCodec` (htsjdk 4.2.0).
//!
//! The third codec GATK's `-L` reaches, and the odd one out: **it opens the file**.
//!
//! ```java
//! public boolean canDecode(final String potentialInput) {
//!     return canDecodeFile(potentialInput, VCF4_MAGIC_HEADER);   // "##fileformat=VCFv4"
//! }
//! ```
//!
//! `BEDCodec` and `IntervalListCodec` answer on the path alone, so `regions.list` holding a BED
//! body is not a Feature file and dies in the interval reader. A `.list` holding a **VCF** body is
//! a Feature file, because this codec reads the first eighteen bytes and finds the magic there.
//! Two files with the same extension and different contents take different branches of `-L`, and
//! only this codec makes that true.
//!
//! # Three attempts, in order, and the second one covers the third
//!
//! ```java
//! return isVCFStream(Files.newInputStream(path), MAGIC) ||
//!        isVCFStream(new GZIPInputStream(...), MAGIC) ||
//!        isVCFStream(new BlockCompressedInputStream(...), MAGIC);
//! ```
//!
//! Plain, then gzip, then block-compressed. A BGZF file **is** a gzip file, so the second attempt
//! already answers for it and the third is unreachable for any well-formed input. It is ported
//! anyway, because "unreachable" here is a property of the format rather than of the code.
//!
//! # A short file is `false`, not an error
//!
//! ```java
//! byte[] buff = new byte[MAGIC_HEADER_LINE.length()];
//! int nread = stream.read(buff, 0, MAGIC_HEADER_LINE.length());
//! boolean eq = Arrays.equals(buff, MAGIC_HEADER_LINE.getBytes());
//! ```
//!
//! `nread` is computed and never used. A file shorter than the magic leaves the tail of `buff` as
//! zero bytes, the comparison fails, and the answer is `false` rather than an exception. The same
//! is true of a file that does not exist: `canDecodeFile` catches `IOException` and answers
//! `false`, so an absent path is simply "not a Feature file" and falls through to the branch that
//! reports it missing.
//!
//! # The version is in the magic
//!
//! `VCF4_MAGIC_HEADER` is `##fileformat=VCFv4`, without a minor version, so `VCFv4.0` through
//! `VCFv4.3` all match and `VCFv3.3` does not. `VCF3Codec` carries the other magic and GATK
//! registers both, so a VCF 3 file is still a Feature file, by a different codec.

/// `VCFCodec.VCF4_MAGIC_HEADER`.
pub const VCF4_MAGIC_HEADER: &[u8] = b"##fileformat=VCFv4";
/// `VCF3Codec.VCF3_MAGIC_HEADER`.
pub const VCF3_MAGIC_HEADER: &[u8] = b"##fileformat=VCFv3";

/// `isVCFStream`: the first `magic.len()` bytes, with a short read leaving zeroes behind.
fn is_vcf_stream(bytes: &[u8], magic: &[u8]) -> bool {
    let mut buffer = vec![0u8; magic.len()];
    let taken = bytes.len().min(magic.len());
    buffer[..taken].copy_from_slice(&bytes[..taken]);
    buffer == magic
}

/// `AbstractVCFCodec.canDecodeFile`, over the bytes of the file rather than its path.
///
/// The path is not consulted at all: this codec is the one that decides by content.
pub fn can_decode_bytes(bytes: &[u8], magic: &[u8]) -> bool {
    if is_vcf_stream(bytes, magic) {
        return true;
    }
    // `new GZIPInputStream(...)` throws on a stream that is not gzip, and the throw is caught by
    // `canDecodeFile`, so a plain non-VCF file answers false here rather than propagating.
    if let Some(inflated) = inflate_prefix(bytes, magic.len()) {
        if is_vcf_stream(&inflated, magic) {
            return true;
        }
    }
    // The block-compressed attempt. Unreachable for well-formed input, since BGZF is gzip.
    false
}

/// The first `wanted` bytes of a gzip member, or `None` if the input is not gzip.
///
/// Only the prefix is needed, and a truncated member is what a short read produces: the reference
/// asks the stream for `magic.len()` bytes and does not care whether more follow.
fn inflate_prefix(bytes: &[u8], wanted: usize) -> Option<Vec<u8>> {
    if bytes.len() < 2 || bytes[0] != 0x1f || bytes[1] != 0x8b {
        return None;
    }
    let mut out = Vec::new();
    let mut decoder = flate2::read::GzDecoder::new(bytes);
    let mut buffer = vec![0u8; wanted];
    loop {
        match std::io::Read::read(&mut decoder, &mut buffer) {
            Ok(0) => break,
            Ok(n) => {
                out.extend_from_slice(&buffer[..n]);
                if out.len() >= wanted {
                    break;
                }
            }
            // A stream that starts with the gzip magic and then fails is the reference's caught
            // IOException, which is `false` and not a refusal.
            Err(_) => return Some(out),
        }
    }
    Some(out)
}

/// `VCFCodec.canDecode`, which is the VCF 4 magic.
pub fn can_decode(bytes: &[u8]) -> bool {
    can_decode_bytes(bytes, VCF4_MAGIC_HEADER)
}

/// `VCF3Codec.canDecode`, registered by GATK beside the other, so a VCF 3 file is a Feature file
/// too.
pub fn can_decode_v3(bytes: &[u8]) -> bool {
    can_decode_bytes(bytes, VCF3_MAGIC_HEADER)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_magic_carries_no_minor_version() {
        assert!(can_decode(b"##fileformat=VCFv4.2\n#CHROM\n"));
        assert!(can_decode(b"##fileformat=VCFv4.0\n"));
        assert!(!can_decode(b"##fileformat=VCFv3.3\n"));
        assert!(can_decode_v3(b"##fileformat=VCFv3.3\n"));
    }

    #[test]
    fn a_short_file_is_false_rather_than_an_error() {
        assert!(!can_decode(b"##fileformat=V"));
        assert!(!can_decode(b""));
    }
}
