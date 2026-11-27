//! `StrictMath.log1p`, which is also `Math.log1p`.
//!
//! Ported from `fdlibm/s_log1p.c` (FDLIBM 5.3, Sun Microsystems), which is what
//! `java.lang.StrictMath.log1p` is specified to be. Same source family and same licence as
//! [`crate::strict_math::exp`] (decision 0025):
//!
//! ```text
//! Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//!
//! Developed at SunSoft, a Sun Microsystems, Inc. business.
//! Permission to use, copy, modify, and distribute this
//! software is freely granted, provided that this notice
//! is preserved.
//! ```
//!
//! # Why this one is exact against `Math` too
//!
//! `exp` and `pow` have a HotSpot intrinsic behind `Math`, so FDLIBM is only `StrictMath` for
//! them. `log1p` has none: `Math.log1p` is declared as `return StrictMath.log1p(x);` and is not an
//! intrinsic candidate, so on every JVM and every architecture the two are one function. The
//! corpus checks that instead of citing it: its `Math` and `StrictMath` columns agree on all
//! 44,996 `log1p` points, and this port matches both (decision 0051).
//!
//! # The algorithm, as FDLIBM states it
//!
//! ```text
//! 1. Argument reduction: find k and f such that 1+x = 2^k * (1+f), sqrt(2)/2 < 1+f < sqrt(2).
//!    If k != 0, f may not be exact: with u = 1+x rounded, c = (1+x)-u, and
//!    log(1+x) - log(u) ~ c/u, so log(u) is computed and c/u added back.
//! 2. With s = f/(2+f), log(1+f) = 2s + s*R, R a degree-14 polynomial in s on [0, 0.1716];
//!    log1p(f) = f - (hfsq - s*(hfsq+R)), hfsq = f*f/2.
//! 3. log1p(x) = k*ln2_hi + (f - (hfsq - (s*(hfsq+R) + k*ln2_lo))).
//! ```
//!
//! As in [`crate::strict_exp`], the order of the floating-point operations is the C's, parentheses
//! included: rewriting any of the returns into an algebraically equal form changes the rounding.

// The literals are FDLIBM's to the digit, as decimals that round to its hex values.
#![allow(clippy::excessive_precision)]

/// `ln2_hi`, `3fe62e42 fee00000`: `n*ln2_hi` is exact for `|n| < 2000`.
const LN2_HI: f64 = 6.931_471_803_691_238_164_90e-1;
/// `ln2_lo`, `3dea39ef 35793c76`.
const LN2_LO: f64 = 1.908_214_929_270_587_700_02e-10;
/// `two54`, `2**54`, used only to build `-inf` for `log1p(-1)`.
const TWO54: f64 = 1.801_439_850_948_198_4e16;
/// `Lp1` to `Lp7`, the polynomial of step 2.
const LP1: f64 = 6.666_666_666_666_735_130e-1;
const LP2: f64 = 3.999_999_999_940_941_908e-1;
const LP3: f64 = 2.857_142_874_366_239_149e-1;
const LP4: f64 = 2.222_219_843_214_978_396e-1;
const LP5: f64 = 1.818_357_216_161_805_012e-1;
const LP6: f64 = 1.531_383_769_920_937_332e-1;
const LP7: f64 = 1.479_819_860_511_658_591e-1;

/// The high 32 bits of `x`, as the C's signed `int` (`__HI`).
fn high(x: f64) -> i32 {
    (x.to_bits() >> 32) as u32 as i32
}

/// `x` with its high word replaced and its low word kept (`__HI(x) = h`).
fn with_high(x: f64, h: i32) -> f64 {
    f64::from_bits((u64::from(h as u32) << 32) | (x.to_bits() & 0xffff_ffff))
}

/// `StrictMath.log1p(x)`, FDLIBM's `log1p`.
pub fn log1p(x: f64) -> f64 {
    let hx = high(x);
    let ax = hx & 0x7fff_ffff;

    // `k`, `f`, `hu` and `c` are assigned on every path that reads them; the C leaves `f` and `c`
    // uninitialised on the k == 0 path, where neither is read before it is set.
    let mut k: i32 = 1;
    let mut f = 0.0;
    let mut hu: i32 = 0;
    let mut c = 0.0;
    if hx < 0x3FDA_827A {
        // x < 0.41422
        if ax >= 0x3ff0_0000 {
            // x <= -1.0
            if x == -1.0 {
                // `-two54/zero`: -inf.
                return -TWO54 / 0.0;
            }
            // `(x-x)/(x-x)`: NaN, the hardware's default one, as the C produces it.
            #[allow(clippy::eq_op)]
            return (x - x) / (x - x);
        }
        if ax < 0x3e20_0000 {
            // |x| < 2**-29. `two54+x>zero` only raises inexact and is always true here.
            if ax < 0x3c90_0000 {
                // |x| < 2**-54
                return x;
            }
            return x - x * x * 0.5;
        }
        // `(int)0xbfd2bec3`: a negative high word, so this reads -0.2929 < x < 0.41422.
        if hx > 0 || hx <= 0xbfd2_bec3_u32 as i32 {
            k = 0;
            f = x;
            hu = 1;
        }
    }
    if hx >= 0x7ff0_0000 {
        return x + x;
    }
    if k != 0 {
        let mut u;
        if hx < 0x4340_0000 {
            u = 1.0 + x;
            hu = high(u);
            k = (hu >> 20) - 1023;
            // The correction term.
            c = if k > 0 { 1.0 - (u - x) } else { x - (u - 1.0) };
            c /= u;
        } else {
            u = x;
            hu = high(u);
            k = (hu >> 20) - 1023;
            c = 0.0;
        }
        hu &= 0x000f_ffff;
        if hu < 0x6a09e {
            // Normalise u.
            u = with_high(u, hu | 0x3ff0_0000);
        } else {
            // Normalise u/2.
            k += 1;
            u = with_high(u, hu | 0x3fe0_0000);
            hu = (0x0010_0000 - hu) >> 2;
        }
        f = u - 1.0;
    }
    let hfsq = 0.5 * f * f;
    let kf = f64::from(k);
    if hu == 0 {
        // |f| < 2**-20
        if f == 0.0 {
            if k == 0 {
                return 0.0;
            }
            c += kf * LN2_LO;
            return kf * LN2_HI + c;
        }
        let r = hfsq * (1.0 - 0.666_666_666_666_666_66 * f);
        if k == 0 {
            return f - r;
        }
        return kf * LN2_HI - ((r - (kf * LN2_LO + c)) - f);
    }
    let s = f / (2.0 + f);
    let z = s * s;
    let r = z * (LP1 + z * (LP2 + z * (LP3 + z * (LP4 + z * (LP5 + z * (LP6 + z * LP7))))));
    if k == 0 {
        return f - (hfsq - s * (hfsq + r));
    }
    kf * LN2_HI - ((hfsq - (s * (hfsq + r) + (kf * LN2_LO + c))) - f)
}
