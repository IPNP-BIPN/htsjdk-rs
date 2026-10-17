//! A CRAM GZIP block as htsjdk 5.0.0 writes one: `GZIPExternalCompressor` over `GzipCodec`.
//!
//! The codec frames the deflate stream itself: ten header bytes with no flags, no time and an OS
//! of 0xFF, then the raw deflate, then the CRC-32 and the length. The deflate is whatever
//! `new DeflaterFactory().makeDeflater(level, true)` hands back, and in Picard 3.5.0 that is
//! libdeflate v1.25 (jlibdeflate 0.1.0) for every run: Picard's `USE_JDK_DEFLATER` and GKL only
//! replace the factory `BlockCompressedOutputStream` uses, and this codec builds its own.

/// `GzipCodec.compress(input, GZIP)` at `level`.
pub fn compress(input: &[u8], level: i32) -> Vec<u8> {
    let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
    // jlibdeflate writes nothing at all for no input, so an empty block is the header and a
    // trailer of zeroes: zlib refuses it, `GzipCodec.decompress` never inflates it.
    if input.is_empty() {
        out.extend_from_slice(&[0; 8]);
        return out;
    }
    // `compress(ByteBuffer)` sizes its output generously, so the one-shot compressor always fits.
    let capacity = input.len() + (input.len() / 16000 + 1) * 5 + 256;
    let mut deflated = vec![0u8; capacity];
    // SAFETY: the compressor is allocated, used once with valid buffers and freed here.
    let written = unsafe {
        let compressor = libdeflate_sys::libdeflate_alloc_compressor(level as std::ffi::c_int);
        assert!(
            !compressor.is_null(),
            "libdeflate could not allocate a compressor"
        );
        let written = libdeflate_sys::libdeflate_deflate_compress(
            compressor,
            input.as_ptr().cast(),
            input.len(),
            deflated.as_mut_ptr().cast(),
            capacity,
        );
        libdeflate_sys::libdeflate_free_compressor(compressor);
        written
    };
    assert!(written != 0, "a GZIP block's output was sized to fit");
    out.extend_from_slice(&deflated[..written]);
    out.extend_from_slice(&crate::compression_header::crc32(input).to_le_bytes());
    out.extend_from_slice(&(input.len() as u32).to_le_bytes());
    out
}

/// `GzipCodec.decompress`: the header parsed by its flags, then the raw deflate inflated until the
/// trailer's length is reached. A length of zero inflates nothing, which is what lets htsjdk read
/// the empty blocks it writes.
pub fn decompress(input: &[u8]) -> Result<Vec<u8>, String> {
    const FEXTRA: u8 = 4;
    const FNAME: u8 = 8;
    const FCOMMENT: u8 = 16;
    const FHCRC: u8 = 2;
    if input.len() < 18 {
        return Err("Input too small to be a valid GZIP block".into());
    }
    if input[0] != 0x1f || input[1] != 0x8b || input[2] != 8 {
        return Err("Invalid GZIP header".into());
    }
    let flags = input[3];
    let mut at = 10;
    if flags & FEXTRA != 0 {
        let xlen = u16::from_le_bytes([input[at], input[at + 1]]) as usize;
        at += 2 + xlen;
    }
    for flag in [FNAME, FCOMMENT] {
        if flags & flag != 0 {
            while *input.get(at).ok_or("truncated")? != 0 {
                at += 1;
            }
            at += 1;
        }
    }
    if flags & FHCRC != 0 {
        at += 2;
    }
    let end = input.len() - 8;
    if at > end {
        return Err("Invalid GZIP block: no room for deflated data and trailer".into());
    }
    let size = u32::from_le_bytes(input[end + 4..end + 8].try_into().unwrap()) as usize;
    let mut out = vec![0u8; size];
    if size > 0 {
        let mut inflater = flate2::Decompress::new(false);
        inflater
            .decompress(&input[at..end], &mut out, flate2::FlushDecompress::Finish)
            .map_err(|e| format!("Error inflating GZIP data: {e}"))?;
        if inflater.total_out() as usize != size {
            return Err(format!(
                "GZIP ISIZE mismatch: expected {size}, got {}",
                inflater.total_out()
            ));
        }
    }
    Ok(out)
}
