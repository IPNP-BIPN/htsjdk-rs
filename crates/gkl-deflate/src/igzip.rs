//! Levels 1 and 2, through the library GKL links rather than through a translation of it.
//!
//! Decision 0034 is why this is a binding and not a port. ISA-L ships two implementations of its
//! own compressor: readable C in `igzip_base.c`, and hand-written SIMD kernels. **They do not
//! agree.** On the same fixtures the C build produces 19749 bytes where the assembly build
//! produces 19044, and GKL ships the assembly. So translating the readable version would have
//! produced a confident wrong answer, and the only cheap exact route is to link the library.
//!
//! That is the trade this repository already makes for the JDK deflater: decision 0001 pins
//! `flate2` to a **vendored C zlib**, not to a Rust reimplementation, for exactly the same reason.
//!
//! ## The configuration is not a guess
//!
//! Decision 0031 found it by trying, because the disassembly is misleading here:
//!
//! ```text
//! level           1                       (2 gives 19141 where GKL gives 19044)
//! level_buf_size  ISAL_DEF_LVL1_DEFAULT   (the 0x141D0 in the disassembly gives 63373, not 63311)
//! end_of_stream   1
//! ```
//!
//! Java levels 1 and 2 both land here, because GKL does not pass the level through to ISA-L.
//!
//! ## What guards it
//!
//! ISA-L falls back to that same base C when it is built without an assembler, when the CPU
//! reports no SSE4.2 (decision 0033), and on any architecture that has no kernels at all. In every
//! one of those states it returns valid deflate that decompresses correctly, so a round-trip check
//! passes and a length check passes; the bytes are simply not GKL's.
//!
//! Nothing in the linked library reports which state it is in. So this module asks it, once, on a
//! 2 KB input GKL has already answered: [`Self::usable`] compresses [`igzip_canary::INPUT`] and
//! compares it with [`igzip_canary::EXPECTED`] byte for byte. If they disagree, every call refuses
//! rather than returning a wrong answer that decompresses.
//!
//! That check is why this is safe to link rather than merely convenient. It is also why running
//! the crate on this project's own development machines is honest: Apple Silicon has no x86
//! kernels, the canary fails there, and levels 1 and 2 say so instead of quietly differing.

use std::os::raw::c_int;
use std::sync::OnceLock;

use isal_sys::igzip_lib as isal;

use crate::igzip_canary;

/// Whether this build and this host reproduce GKL's igzip.
///
/// Computed once. The cost is one 2 KB compression, which is less than the cost of finding out
/// from a failing byte comparison three layers downstream.
pub fn usable() -> bool {
    static USABLE: OnceLock<bool> = OnceLock::new();
    *USABLE.get_or_init(|| deflate_unchecked(&igzip_canary::INPUT) == igzip_canary::EXPECTED)
}

/// Compress `data` exactly as `IntelDeflater` does at Java levels 1 and 2.
///
/// # Panics
///
/// If ISA-L refuses the input. `isal_deflate_stateless` can return `STATELESS_OVERFLOW` when the
/// output buffer is too small to hold even a stored block; the buffer here is sized from ISA-L's
/// own worst case, so a refusal means an assumption broke rather than that the caller was unlucky,
/// and it should be loud.
pub fn deflate(data: &[u8]) -> Vec<u8> {
    assert!(
        usable(),
        "this build of ISA-L does not reproduce GKL's igzip. It falls back to its readable C when \
         built without an assembler, on a CPU without SSE4.2, or on a non-x86 architecture, and \
         that C finds different matches (decision 0034). Levels 1 and 2 refuse rather than return \
         valid deflate that is not GKL's."
    );
    deflate_unchecked(data)
}

fn deflate_unchecked(data: &[u8]) -> Vec<u8> {
    // ISA-L's own stored-block worst case: a 5-byte header per 65535-byte block, plus the data.
    // Sized from the rule rather than from a guess, so a large incompressible input cannot quietly
    // land in the overflow path.
    const TYPE0_BLK_HDR_LEN: usize = 5;
    const TYPE0_MAX_BLK_LEN: usize = 65535;
    let worst_case =
        data.len() + TYPE0_BLK_HDR_LEN * (data.len().div_ceil(TYPE0_MAX_BLK_LEN) + 1) + 64;

    let mut out = vec![0u8; worst_case];
    let mut level_buf = vec![0u8; isal::ISAL_DEF_LVL1_DEFAULT as usize];
    // `next_in` is `*mut u8` in the binding although ISA-L only reads it, so the input is copied
    // rather than cast away from a shared reference.
    let mut input = data.to_vec();

    let written = unsafe {
        let mut stream: isal::isal_zstream = std::mem::zeroed();
        isal::isal_deflate_stateless_init(&mut stream);
        stream.level = 1;
        stream.level_buf = level_buf.as_mut_ptr();
        stream.level_buf_size = level_buf.len() as u32;
        stream.next_in = input.as_mut_ptr();
        stream.avail_in = input.len() as u32;
        stream.next_out = out.as_mut_ptr();
        stream.avail_out = out.len() as u32;
        stream.end_of_stream = 1;
        stream.flush = isal::NO_FLUSH as u16;

        let rc = isal::isal_deflate_stateless(&mut stream);
        assert_eq!(
            rc,
            isal::COMP_OK as c_int,
            "isal_deflate_stateless refused a {}-byte input with {rc}",
            data.len()
        );
        (out.len() - stream.avail_out as usize) as usize
    };
    out.truncate(written);
    out
}
