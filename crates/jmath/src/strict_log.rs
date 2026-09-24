//! `StrictMath.log`.
//!
//! Ported from `fdlibm/e_log.c` (FDLIBM 5.3, Sun Microsystems), which is what
//! `java.lang.StrictMath.log` is specified to be. Same standing as [`crate::strict_exp`]: a
//! permissively licensed source, preserved notice, and a claim checked on every point of the
//! corpus rather than taken from the specification.
//!
//! ```text
//! Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//!
//! Developed at SunPro, a Sun Microsystems, Inc. business.
//! Permission to use, copy, modify, and distribute this software
//! is freely granted, provided that this notice is preserved.
//! ```
//!
//! ## Why it is needed when `Math.log` is already exact
//!
//! [`crate::math::log`] is correctly rounded, because `Math.log` measured correctly rounded
//! (decision 0006). `StrictMath.log` is not: it is FDLIBM's algorithm, within one ulp, and it
//! differs from the correctly rounded answer on 186 of the corpus's 44,996 points. A call site
//! that reaches `StrictMath.log` therefore needs this function and not the other. The one GATK
//! needs first is `java.util.Random.nextGaussian`, whose polar method takes `StrictMath.log` of a
//! uniform draw: `QualByDepth` jitters every QD above 35 with it, so a GenotypeGVCFs QD is only
//! reproducible through here.
//!
//! ## The algorithm, as FDLIBM states it
//!
//! ```text
//! 1. Argument reduction: find k and f such that
//!         x = 2^k * (1+f),
//!    where  sqrt(2)/2 < 1+f < sqrt(2) .
//!
//! 2. Approximation of log(1+f).
//!    Let s = f/(2+f) ; based on log(1+f) = log(1+s) - log(1-s)
//!              = 2s + 2/3 s**3 + 2/5 s**5 + .....,
//!              = 2s + s*R
//!    We use a special Remez algorithm on [0,0.1716] to generate
//!    a polynomial of degree 14 to approximate R. The maximum error
//!    of this polynomial approximation is bounded by 2**-58.45.
//!
//! 3. Finally, log(x) = k*ln2 + log(1+f).
//!                    = k*ln2_hi+(f-(hfsq-(s*(hfsq+R)+k*ln2_lo)))
//!    Here ln2 is split into two floating point numbers:
//!              ln2_hi + ln2_lo,
//!    where n*ln2_hi is always exact for |n| < 2000.
//! ```
//!
//! `clippy::excessive_precision` is allowed module-wide: the constants carry FDLIBM's printed
//! digits, which say more than an `f64` holds, and trimming them would make the port cite a
//! value its source never wrote. They round to the same bits either way.

#![allow(clippy::excessive_precision)]

/// `ln2_hi`: the leading half of `ln2`, whose low bits are zero so that `k * ln2_hi` is exact.
///
/// `clippy::approx_constant` flags it as an approximation of `LN_2`, which it deliberately is:
/// see [`crate::strict_exp`] for why FDLIBM's literals are kept rather than substituted.
#[allow(clippy::approx_constant)]
const LN2_HI: f64 = 6.931_471_803_691_238_164_90e-1;
/// `ln2_lo`: the trailing half.
const LN2_LO: f64 = 1.908_214_929_270_587_700_02e-10;
/// `two54`, which scales a subnormal argument up into the normal range.
const TWO54: f64 = 1.801_439_850_948_198_400_00e16;
/// `Lg1` to `Lg7`, the Remez coefficients of step 2.
const LG1: f64 = 6.666_666_666_666_735_130e-1;
const LG2: f64 = 3.999_999_999_940_941_908e-1;
const LG3: f64 = 2.857_142_874_366_239_149e-1;
const LG4: f64 = 2.222_219_843_214_978_396e-1;
const LG5: f64 = 1.818_357_216_161_805_012e-1;
const LG6: f64 = 1.531_383_769_920_937_332e-1;
const LG7: f64 = 1.479_819_860_511_658_591e-1;

fn high_word(x: f64) -> i32 {
    (x.to_bits() >> 32) as u32 as i32
}

fn low_word(x: f64) -> u32 {
    x.to_bits() as u32
}

fn with_high_word(x: f64, high: i32) -> f64 {
    f64::from_bits((u64::from(high as u32) << 32) | u64::from(low_word(x)))
}

/// `StrictMath.log(x)`, which is FDLIBM's `__ieee754_log`.
///
/// The operations are the C function's, in its order and with its parentheses: the result is
/// reproducible only because nothing is rearranged into an algebraically equal form.
pub fn log(x: f64) -> f64 {
    let mut x = x;
    let mut hx = high_word(x);
    let lx = low_word(x);
    let mut k: i32 = 0;

    if hx < 0x0010_0000 {
        // x < 2**-1022
        if ((hx & 0x7fff_ffff) as u32 | lx) == 0 {
            // log(+-0) = -inf. The C writes `-two54/zero`.
            return f64::NEG_INFINITY;
        }
        if hx < 0 {
            // log(-#) = NaN. The C writes `(x-x)/zero`, kept as written rather than as `f64::NAN`
            // because the NaN an FPU produces for 0/0 need not have the constant's bits.
            #[allow(clippy::eq_op)]
            let zero = x - x;
            return zero / 0.0;
        }
        // Subnormal: scale up.
        k -= 54;
        x *= TWO54;
        hx = high_word(x);
    }
    if hx >= 0x7ff0_0000 {
        return x + x;
    }
    k += (hx >> 20) - 1023;
    hx &= 0x000f_ffff;
    let i = (hx + 0x95f64) & 0x10_0000;
    // Normalise x or x/2.
    x = with_high_word(x, hx | (i ^ 0x3ff0_0000));
    k += i >> 20;
    let f = x - 1.0;
    if (0x000f_ffff & (2 + hx)) < 3 {
        // |f| < 2**-20
        if f == 0.0 {
            if k == 0 {
                return 0.0;
            }
            let dk = f64::from(k);
            return dk * LN2_HI + dk * LN2_LO;
        }
        let r = f * f * (0.5 - 0.333_333_333_333_333_33 * f);
        if k == 0 {
            return f - r;
        }
        let dk = f64::from(k);
        return dk * LN2_HI - ((r - dk * LN2_LO) - f);
    }
    let s = f / (2.0 + f);
    let dk = f64::from(k);
    let z = s * s;
    let mut i = hx - 0x6147a;
    let w = z * z;
    let j = 0x6b851 - hx;
    let t1 = w * (LG2 + w * (LG4 + w * LG6));
    let t2 = z * (LG1 + w * (LG3 + w * (LG5 + w * LG7)));
    i |= j;
    let r = t2 + t1;
    if i > 0 {
        let hfsq = 0.5 * f * f;
        if k == 0 {
            f - (hfsq - s * (hfsq + r))
        } else {
            dk * LN2_HI - ((hfsq - (s * (hfsq + r) + dk * LN2_LO)) - f)
        }
    } else if k == 0 {
        f - s * (f - r)
    } else {
        dk * LN2_HI - ((s * (f - r) - dk * LN2_LO) - f)
    }
}
