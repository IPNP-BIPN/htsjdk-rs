//! ISA-L itself, linked, for the `isal` feature's cross-check of the port in `igzip`.
//!
//! This used to be how levels 1 and 2 were answered. It reproduces GKL only where ISA-L runs its
//! x86 kernels: built with nasm, on a CPU with SSE4.2. Anywhere else the library falls back to its
//! readable C, which makes different bytes (decision 0034), so [`usable`] checks the canary first.
//! The configuration is GKL's: the Java level, `ISAL_DEF_LVL1_DEFAULT` or `ISAL_DEF_LVL2_DEFAULT`,
//! `end_of_stream` set (decisions 0031 and 0043).

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
    *USABLE.get_or_init(|| deflate_unchecked(&igzip_canary::INPUT, 1) == igzip_canary::EXPECTED)
}

/// Compress `data` exactly as `IntelDeflater` does at Java level `level`, which is 1 or 2.
///
/// # Panics
///
/// If ISA-L refuses the input. `isal_deflate_stateless` can return `STATELESS_OVERFLOW` when the
/// output buffer is too small to hold even a stored block; the buffer here is sized from ISA-L's
/// own worst case, so a refusal means an assumption broke rather than that the caller was unlucky,
/// and it should be loud.
pub fn deflate(data: &[u8], level: usize) -> Vec<u8> {
    assert!(
        usable(),
        "this build of ISA-L does not reproduce GKL's igzip. It falls back to its readable C when \
         built without an assembler, on a CPU without SSE4.2, or on a non-x86 architecture, and \
         that C finds different matches (decision 0034). Levels 1 and 2 refuse rather than return \
         valid deflate that is not GKL's."
    );
    deflate_unchecked(data, level)
}

fn deflate_unchecked(data: &[u8], level: usize) -> Vec<u8> {
    // ISA-L's own stored-block worst case: a 5-byte header per 65535-byte block, plus the data.
    // Sized from the rule rather than from a guess, so a large incompressible input cannot quietly
    // land in the overflow path.
    const TYPE0_BLK_HDR_LEN: usize = 5;
    const TYPE0_MAX_BLK_LEN: usize = 65535;
    let worst_case =
        data.len() + TYPE0_BLK_HDR_LEN * (data.len().div_ceil(TYPE0_MAX_BLK_LEN) + 1) + 64;

    let mut out = vec![0u8; worst_case];
    // GKL passes the Java level through: 1 is ISA-L level 1 and 2 is ISA-L level 2, each with its
    // own default token buffer (decision 0043, which corrects 0031).
    let buffer = if level == 2 {
        isal::ISAL_DEF_LVL2_DEFAULT
    } else {
        isal::ISAL_DEF_LVL1_DEFAULT
    };
    let mut level_buf = vec![0u8; buffer as usize];
    // `next_in` is `*mut u8` in the binding although ISA-L only reads it, so the input is copied
    // rather than cast away from a shared reference.
    let mut input = data.to_vec();

    let written = unsafe {
        let mut stream: isal::isal_zstream = std::mem::zeroed();
        isal::isal_deflate_stateless_init(&mut stream);
        stream.level = level as u32;
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
