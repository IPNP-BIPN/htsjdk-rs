//! `StrictMath.expm1`, which is also `Math.expm1`.
//!
//! Ported from `fdlibm/s_expm1.c` (FDLIBM 5.3, Sun Microsystems), which is what
//! `java.lang.StrictMath.expm1` is specified to be. Same standing as [`crate::strict_log1p`]:
//!
//! ```text
//! Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
//!
//! Permission to use, copy, modify, and distribute this
//! software is freely granted, provided that this notice
//! is preserved.
//! ```
//!
//! `Math.expm1` is `return StrictMath.expm1(x);` with no intrinsic behind it, so this is exact
//! against `Math` as well, and the corpus's two columns agree on all 44,996 `expm1` points
//! (decision 0051).
//!
//! # The algorithm, as FDLIBM states it
//!
//! ```text
//! 1. Argument reduction: x = k*ln2 + r, |r| <= 0.5*ln2, with r kept as hi-lo and the rounding
//!    error of hi-lo carried in c.
//! 2. On [0, 0.5*ln2], a rational approximation in r*r/2 with coefficients Q1..Q5:
//!    expm1(r) = r - (r*e - r*r/2), e built from R1(r*r/2).
//! 3. Scale back: expm1(x) = 2^k * (expm1(r) + 1) - 1, written seven ways by the size of k so the
//!    final subtraction loses nothing.
//! ```
//!
//! The operation order is the C's, parentheses included, for the reason [`crate::strict_exp`]
//! gives: an algebraically equal rewrite rounds differently.

// The literals are FDLIBM's to the digit, as decimals that round to its hex values.
#![allow(clippy::excessive_precision)]

const ONE: f64 = 1.0;
/// `huge`, which only raises the overflow and inexact flags.
const HUGE: f64 = 1.0e300;
/// `tiny`, likewise.
const TINY: f64 = 1.0e-300;
/// `o_threshold`, `0x40862E42 0xFEFA39EF`.
const O_THRESHOLD: f64 = 7.097_827_128_933_839_730_96e2;
/// `ln2_hi`, `0x3fe62e42 0xfee00000`.
const LN2_HI: f64 = 6.931_471_803_691_238_164_90e-1;
/// `ln2_lo`, `0x3dea39ef 0x35793c76`.
const LN2_LO: f64 = 1.908_214_929_270_587_700_02e-10;
/// `invln2`, `0x3ff71547 0x652b82fe`. FDLIBM's literal rather than `LOG2_E`, as in
/// [`crate::strict_exp`].
#[allow(clippy::approx_constant)]
const INVLN2: f64 = 1.442_695_040_888_963_387_00;
/// `Q1` to `Q5`, the scaled coefficients of step 2.
const Q1: f64 = -3.333_333_333_333_313_164_28e-2;
const Q2: f64 = 1.587_301_587_254_814_601_65e-3;
const Q3: f64 = -7.936_507_578_674_879_424_73e-5;
const Q4: f64 = 4.008_217_827_329_362_395_52e-6;
const Q5: f64 = -2.010_992_181_836_243_713_26e-7;

/// `y` with `k << 20` added to its high word (`__HI(y) += (k<<20)`), which adds `k` to the
/// exponent without rounding.
fn add_to_exponent(y: f64, k: i32) -> f64 {
    let bits = y.to_bits();
    let high = ((bits >> 32) as u32).wrapping_add((k << 20) as u32);
    f64::from_bits((u64::from(high) << 32) | (bits & 0xffff_ffff))
}

/// `1.0` with its high word replaced (`t = one; __HI(t) = h`): the low word of `1.0` is zero.
fn from_high(h: u32) -> f64 {
    f64::from_bits(u64::from(h) << 32)
}

/// `StrictMath.expm1(x)`, FDLIBM's `expm1`.
pub fn expm1(x: f64) -> f64 {
    let bits = x.to_bits();
    let xsb = (bits >> 63) != 0;
    let hx = ((bits >> 32) as u32) & 0x7fff_ffff;
    let lx = (bits & 0xffff_ffff) as u32;

    // Huge and non-finite arguments.
    if hx >= 0x4043_687A {
        // |x| >= 56*ln2
        if hx >= 0x4086_2E42 {
            // |x| >= 709.78...
            if hx >= 0x7ff0_0000 {
                if ((hx & 0xfffff) | lx) != 0 {
                    // NaN
                    return x + x;
                }
                // expm1(+-inf) = {inf, -1}
                return if xsb { -1.0 } else { x };
            }
            if x > O_THRESHOLD {
                return HUGE * HUGE;
            }
        }
        if xsb {
            // x < -56*ln2: -1 with inexact. `x+tiny<0.0` only raises the flag and is always
            // true here.
            return TINY - ONE;
        }
    }

    // Argument reduction.
    let mut x = x;
    let k: i32;
    let mut c = 0.0;
    if hx > 0x3fd6_2e42 {
        // |x| > 0.5 ln2
        let (hi, lo);
        if hx < 0x3FF0_A2B2 {
            // and |x| < 1.5 ln2
            if xsb {
                hi = x + LN2_HI;
                lo = -LN2_LO;
                k = -1;
            } else {
                hi = x - LN2_HI;
                lo = LN2_LO;
                k = 1;
            }
        } else {
            // C's double-to-int conversion truncates toward zero, as `as` does.
            k = (INVLN2 * x + if xsb { -0.5 } else { 0.5 }) as i32;
            let t = f64::from(k);
            // t*ln2_hi is exact here.
            hi = x - t * LN2_HI;
            lo = t * LN2_LO;
        }
        x = hi - lo;
        c = (hi - x) - lo;
    } else if hx < 0x3c90_0000 {
        // |x| < 2**-54: return x, with the C's inexact-raising expression kept because it is
        // also the value: `x - (t-(huge+x))` with `t = huge+x`.
        let t = HUGE + x;
        return x - (t - (HUGE + x));
    } else {
        k = 0;
    }

    // x is now in the primary range.
    let hfx = 0.5 * x;
    let hxs = x * hfx;
    let r1 = ONE + hxs * (Q1 + hxs * (Q2 + hxs * (Q3 + hxs * (Q4 + hxs * Q5))));
    let t = 3.0 - r1 * hfx;
    let mut e = hxs * ((r1 - t) / (6.0 - x * t));
    if k == 0 {
        // c is 0.
        return x - (x * e - hxs);
    }
    e = x * (e - c) - c;
    e -= hxs;
    if k == -1 {
        return 0.5 * (x - e) - 0.5;
    }
    if k == 1 {
        if x < -0.25 {
            return -2.0 * (e - (x + 0.5));
        }
        return ONE + 2.0 * (x - e);
    }
    if k <= -2 || k > 56 {
        // Suffices to return exp(x)-1.
        let y = ONE - (e - x);
        return add_to_exponent(y, k) - ONE;
    }
    if k < 20 {
        // t = 1-2^-k
        let t = from_high(0x3ff0_0000 - (0x0020_0000_u32 >> k));
        let y = t - (e - x);
        add_to_exponent(y, k)
    } else {
        // t = 2^-k
        let t = from_high(((0x3ff - k) as u32) << 20);
        let mut y = x - (e + t);
        y += ONE;
        add_to_exponent(y, k)
    }
}
