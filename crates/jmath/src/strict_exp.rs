//! `StrictMath.exp`.
//!
//! Ported from `fdlibm/e_exp.c` (FDLIBM 5.3, Sun Microsystems), which is what
//! `java.lang.StrictMath.exp` is specified to be.
//!
//! ## Why this file exists, and why it is not the file decision 0014 removed
//!
//! Decision 0014 withdrew `jmath::math::exp` because it was an operation-by-operation
//! transcription of HotSpot's x86 intrinsic, whose source is **GPL2 with no Classpath Exception**
//! and therefore cannot ship in an MIT crate. That is still true, and nothing here changes it.
//!
//! This is a different source with a different licence. FDLIBM is Sun's freely-distributable math
//! library, and its notice — preserved below, as it requires — grants exactly what a port needs:
//!
//! ```text
//! Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//!
//! Developed at SunPro, a Sun Microsystems, Inc. business.
//! Permission to use, copy, modify, and distribute this software
//! is freely granted, provided that this notice is preserved.
//! ```
//!
//! `java.lang.StrictMath` is **specified** to be FDLIBM: "To help ensure the portability of Java
//! programs, the definitions of some of the numeric functions in this package require that they
//! produce the same results as certain published algorithms ... the algorithms contained in the
//! `fdlibm` package". So this is not an approximation of `StrictMath.exp`; it is the algorithm
//! `StrictMath.exp` is defined as, and the corpus checks that claim on every point rather than
//! taking the specification's word for it.
//!
//! ## What it does and does not buy
//!
//! It buys `StrictMath.exp` exactly. It does **not** buy `Math.exp`: `Math` is free to use a
//! faster intrinsic, and on x86-64 HotSpot does. The two are different functions and decision 0005
//! is about exactly that. What this file makes possible is a *measurement*: whether FDLIBM is
//! closer to `Math.exp` than the system libm the corpus currently routes `exp` through, and by how
//! much. `tests/conformance.rs` reports both rates, so the cost of not having the intrinsic stops
//! being a single number of unknown provenance.
//!
//! Every call site that reaches `Math.exp` still cannot be reproduced bit-for-bit. What changes is
//! that a call site reaching `StrictMath.exp` now can, and that the gap for the others is
//! quantified against the best permissively-licensed implementation rather than against whatever
//! libm the host happens to ship.
//!
//! ## The algorithm, as FDLIBM states it
//!
//! ```text
//! 1. Argument reduction: reduce x to r so that |r| <= 0.5*ln2.
//!    Given x, find r and integer k such that
//!         x = k*ln2 + r,  |r| <= 0.5*ln2.
//!    (inexact, but the error is < 2**-60)
//!
//! 2. Approximation of exp(r) by a special rational function on [0, 0.34658]:
//!         R(r**2) = r*(exp(r)+1)/(exp(r)-1) = 2 + r*r/6 - r**4/360 + ...
//!    so that
//!         exp(r) = 1 + r + r*R/(2 - R)
//!
//! 3. Scale back: exp(x) = 2**k * exp(r)
//! ```
//!
//! The reduction is written as two subtractions against a split `ln2` (`ln2HI` and `ln2LO`)
//! because a single `f64` cannot hold `ln2` to the precision the error bound needs. That split is
//! the whole reason step 1 is exact enough to be worth doing.

/// `one`.
const ONE: f64 = 1.0;
/// `halF[2]`, indexed by the sign of x. The negative half is stored rather than negated at use.
const HALF: [f64; 2] = [0.5, -0.5];
/// `huge`, used only to raise the overflow flag; its value never reaches the result.
const HUGE: f64 = 1.0e300;
/// `twom1000`, `2**-1000`, which the subnormal path multiplies back in.
const TWOM1000: f64 = 9.332_636_185_032_189e-302;
/// `o_threshold`: above this, `exp` overflows.
const O_THRESHOLD: f64 = 7.097_827_128_933_84e2;
/// `u_threshold`: below this, `exp` underflows to zero.
const U_THRESHOLD: f64 = -7.451_332_191_019_411e2;
/// `ln2HI[2]`: the leading half of `ln2`, exact in an `f64`.
const LN2_HI: [f64; 2] = [6.931_471_803_691_238e-1, -6.931_471_803_691_238e-1];
/// `ln2LO[2]`: the trailing half, which is what makes the reduction accurate to 2**-60.
const LN2_LO: [f64; 2] = [1.908_214_929_270_587_7e-10, -1.908_214_929_270_587_7e-10];
/// `invln2`.
///
/// `clippy::approx_constant` flags this as an approximation of `std::f64::consts::LOG2_E` and the
/// advice is refused here for the same reason [`crate::log`] refuses it: the value is FDLIBM's
/// literal, and substituting a constant from elsewhere would make the port depend on two sources
/// agreeing to the last bit rather than on the one it cites. They do agree today. That is not the
/// same as being allowed to assume it.
#[allow(clippy::approx_constant)]
const INVLN2: f64 = 1.442_695_040_888_963_4;
/// `P1` to `P5`, the coefficients of the rational approximation of step 2.
const P1: f64 = 1.666_666_666_666_660_2e-1;
const P2: f64 = -2.777_777_777_701_559_3e-3;
const P3: f64 = 6.613_756_321_437_934e-5;
const P4: f64 = -1.653_390_220_546_525_2e-6;
const P5: f64 = 4.138_136_797_057_238_6e-8;

/// `StrictMath.exp(x)`, which is FDLIBM's `__ieee754_exp`.
///
/// The structure is the C function's, including the order of the floating-point operations, which
/// is what makes the result reproducible: rearranging `y = one - ((lo - (x*c)/(2.0 - c)) - hi)`
/// into anything algebraically equivalent changes the rounding and therefore the answer.
pub fn exp(x: f64) -> f64 {
    let bits = x.to_bits();
    // `hx`: the high word with the sign cleared; `xsb`: the sign bit, used as an index into the
    // paired constants above.
    let xsb = ((bits >> 63) & 1) as usize;
    let hx = ((bits >> 32) & 0x7fff_ffff) as u32;

    // Filter out non-finite arguments.
    if hx >= 0x4086_2E42 {
        // |x| >= 709.78...
        if hx >= 0x7ff0_0000 {
            let lx = (bits & 0xffff_ffff) as u32;
            if ((hx & 0xf_ffff) | lx) != 0 {
                // NaN: FDLIBM returns x+x, which quiets a signalling NaN.
                return x + x;
            }
            // +-Inf: exp(+inf) = inf, exp(-inf) = 0.
            return if xsb == 0 { x } else { 0.0 };
        }
        if x > O_THRESHOLD {
            // Overflow. `huge*huge` is how the C raises the flag; the value is +Inf either way.
            return HUGE * HUGE;
        }
        if x < U_THRESHOLD {
            // Underflow. Likewise `twom1000*twom1000`, which is +0.0.
            return TWOM1000 * TWOM1000;
        }
    }

    // Argument reduction.
    if hx > 0x3fd6_2e42 {
        // |x| > 0.5*ln2
        let (hi, lo, k);
        if hx < 0x3FF0_A2B2 {
            // |x| < 1.5*ln2: k is +-1 and the multiply is avoided.
            hi = x - LN2_HI[xsb];
            lo = LN2_LO[xsb];
            k = 1 - (xsb as i32) - (xsb as i32);
        } else {
            k = (INVLN2 * x + HALF[xsb]) as i32;
            let t = f64::from(k);
            // The two subtractions are what the split constant is for.
            hi = x - t * LN2_HI[0];
            lo = t * LN2_LO[0];
        }
        // Not folded into the polynomial: the C keeps this as its own rounding step.
        let x = hi - lo;
        finish(x, hi, lo, k)
    } else if hx < 0x3e30_0000 {
        // |x| < 2**-28. FDLIBM writes `if(huge+x>one) return one+x;`, whose only purpose is to
        // raise the inexact flag: the comparison is true for every x in range, so the result is
        // `1+x` unconditionally. Kept as one return rather than two identical ones, with the
        // reason here, because Rust has no way to raise the flag anyway.
        ONE + x
    } else {
        // |x| <= 0.5*ln2: no reduction, and `hi`/`lo` are unused by the k == 0 path.
        finish(x, 0.0, 0.0, 0)
    }
}

/// Steps 2 and 3: the rational approximation, then the scaling by `2**k`.
///
/// Split out of [`exp`] because the C reaches it from two places with different `hi`/`lo`, and
/// duplicating it would be two chances to get the operation order wrong.
fn finish(x: f64, hi: f64, lo: f64, k: i32) -> f64 {
    let t = x * x;
    let c = x - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
    // `y = one - ((lo - (x*c)/(2.0 - c)) - hi)` when k != 0, and the shorter form when k == 0.
    // Written exactly as the C writes them, parentheses included.
    if k == 0 {
        return ONE - ((x * c) / (c - 2.0) - x);
    }
    let y = ONE - ((lo - (x * c) / (2.0 - c)) - hi);

    if k >= -1021 {
        // `SET_HIGH_WORD(y, HIGH_WORD(y) + (k<<20))`: scaling by adding to the exponent field,
        // which is exact and, unlike a multiply by `2**k`, cannot round twice.
        let bits = y.to_bits();
        let high = ((bits >> 32) as u32).wrapping_add((k as u32) << 20);
        f64::from_bits((u64::from(high) << 32) | (bits & 0xffff_ffff))
    } else {
        // The result is subnormal: add 1000 to the exponent first, then scale back down by
        // `2**-1000`, so the intermediate stays normal and only the final multiply rounds.
        let bits = y.to_bits();
        let high = ((bits >> 32) as u32).wrapping_add(((k + 1000) as u32) << 20);
        let scaled = f64::from_bits((u64::from(high) << 32) | (bits & 0xffff_ffff));
        scaled * TWOM1000
    }
}
