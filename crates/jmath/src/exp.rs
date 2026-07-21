//! `java.lang.Math.exp`, reproducing HotSpot's x86 intrinsic.
//!
//! Ported from openjdk/jdk `src/hotspot/cpu/x86/macroAssembler_x86_exp.cpp` at `jdk-17-ga`
//! (`MacroAssembler::fast_exp`), which is Intel LIBM's table-driven exp.
//!
//! Unlike `log` and `log10`, `Math.exp` is **not** correctly rounded (measured 99.9833%
//! against MPFR), so rounding the true result is the wrong target and the algorithm itself
//! has to be reproduced. See
//! `docs/decisions/0006-correct-rounding-is-the-target-for-log-and-log10.md`.
//!
//! # Why this is transcribed operation by operation
//!
//! The original is SSE2 packed-double code, and the two lanes of each register carry
//! *different* polynomial terms that are combined at specific points. The result depends on
//! that exact association and on the exact order of the final accumulation, because every
//! intermediate rounds to `f64`. Restructuring the arithmetic into a more natural Horner form
//! would be algebraically identical and give different bits.
//!
//! ```text
//!   exp(x) = 2^n * T[j] * (1 + P(r))
//!   x = m*ln2/64 + r,   m = n*64 + j,   j in [0,64)
//! ```
//!
//! `clippy::excessive_precision` is allowed module-wide: the literals are transcribed verbatim
//! from HotSpot's `_cv` table and its algorithm comment, digits included. Trimming them to what
//! `f64` can hold would make the constants harder to diff against the reference, which is the
//! only check that matters here.
#![allow(clippy::excessive_precision)]

use crate::exp_table::{T_HI, T_LO};

// These are HotSpot's constants, transcribed from `_cv`, not the mathematically nearest
// doubles. `LN2_64_HI` in particular has its low mantissa bits deliberately cleared
// (`0x3f862e42fefa0000`) so that `m * LN2_64_HI` is exact for the integers `m` that occur,
// which is what makes the two-step argument reduction work. Replacing any of these with a
// `std::f64::consts` value, or with a more accurate literal, changes the result bits.
/// 64 / ln(2)
const INV_LN2_64: f64 = 92.33248261689366;
/// ln(2)/64, high part, low mantissa bits cleared on purpose.
const LN2_64_HI: f64 = 0.010830424696223417;
/// ln(2)/64, the remainder the high part cannot represent.
const LN2_64_LO: f64 = 2.572804622327669e-14;

/// The classic round-to-nearest-integer trick: adding 1.5 * 2^52 forces the fraction out.
const SHIFTER: f64 = 6_755_399_441_055_744.0;

// Minimax coefficients. Their pairing into SIMD lanes is preserved below, because it
// determines which products are formed and therefore which roundings occur.
const C_HALF: f64 = 0.4999999999999999; // cv[48], the r^2/2 term
const C4: f64 = 0.001388874731341548; // lane 0, multiplies r
const C5: f64 = 0.04166666666952449; // lane 1, multiplies r
const C2: f64 = 0.008333368243348653; // lane 0 constant
const C3: f64 = 0.16666666666523347; // lane 1 constant

/// Overflow and underflow thresholds quoted in the HotSpot algorithm description.
const OVERFLOW_THRESHOLD: f64 = 709.782712893383973096;
const UNDERFLOW_THRESHOLD: f64 = -745.133219101941108420;

pub fn exp(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == f64::INFINITY {
        return f64::INFINITY;
    }
    if x == f64::NEG_INFINITY {
        return 0.0;
    }
    if x > OVERFLOW_THRESHOLD {
        return f64::INFINITY;
    }
    if x < UNDERFLOW_THRESHOLD {
        return 0.0;
    }

    // The intrinsic takes its fast path only for |x| in roughly [2^-54, 2^10); outside it the
    // result is 1.0 to within a rounding.
    let ax = x.abs();
    if ax < 5.551_115_123_125_783e-17 {
        return 1.0 + x;
    }

    // Argument reduction. `m` is the nearest integer to x * 64/ln2, recovered from the low
    // bits of the shifted double exactly as the original does.
    let shifted = x * INV_LN2_64 + SHIFTER;
    let m = shifted.to_bits() as u32 as i32;
    let fm = shifted - SHIFTER;

    let j = (m & 63) as usize;
    let n = m >> 6;

    // r = (x - m*hi) - m*lo, in this order: the two subtractions do not commute in f64.
    let r = (x - fm * LN2_64_HI) - fm * LN2_64_LO;

    // Two independent polynomial lanes, exactly as the packed registers hold them.
    let r2 = r * r;
    let r3 = r * r2;
    let lane0_poly = C2 + C4 * r; // xmm5 low
    let lane1_poly = C3 + C5 * r; // xmm5 high
    let lane0 = (r3 * r2) * lane0_poly; // mulsd made this r^5 in the low lane only
    let lane1 = r3 * lane1_poly;
    let half_r2 = r2 * C_HALF;

    // Reassemble the table entry with its runtime exponent, which the original does with a
    // single `por` into the mantissa-only T_HI.
    let scaled_exponent = ((n as i64 + 1023) as u64) << 52;
    let t = f64::from_bits(T_HI[j] | scaled_exponent);

    // Accumulation order is load-bearing: each step rounds. Written to mirror the original's
    // `addsd` sequence rather than to read naturally.
    let mut acc = r + T_LO[j];
    acc += lane0;
    acc += lane1;
    acc += half_r2;

    // A separate multiply and add, not a fused multiply-add: the original emits `mulsd`
    // followed by `addsd`, so the product rounds before the sum.
    acc * t + t
}
