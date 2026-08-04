//! `StrictMath.pow`.
//!
//! Ported from `fdlibm/e_pow.c` (FDLIBM 5.3, Sun Microsystems), which is what
//! `java.lang.StrictMath.pow` is specified to be. Same source and same licence as
//! [`crate::strict_math::exp`], and for the same reason: decision 0014 withdrew the HotSpot
//! transcription because `java.base` and HotSpot are GPL2, and FDLIBM is not.
//!
//! ```text
//! Copyright (C) 2004 by Sun Microsystems, Inc. All rights reserved.
//!
//! Permission to use, copy, modify, and distribute this
//! software is freely granted, provided that this notice
//! is preserved.
//! ```
//!
//! # Why this file exists
//!
//! `Math.pow` is the one intrinsic decision 0007 deferred, and the reason was specific: HotSpot's
//! `pow` uses `rcpps`, the packed approximate reciprocal, at six sites, and unlike `log` it does
//! not refine the approximation away. The hazard that raised — that `pow`'s bits might depend on
//! the silicon — was tested and did not materialise: 404,883 points regenerated on real AMD EPYC
//! matched the emulated corpus exactly.
//!
//! What was never produced is a **bound**. The record holds a rate — 99.9378% agreement with the
//! host libm, against `exp`'s 99.9711% — and a rate says how often, not how far. Decision 0025 did
//! for `exp` exactly what this file exists to make possible for `pow`: port the permissively
//! licensed implementation, then measure the distance to the intrinsic rather than assume it.
//!
//! # What it buys and what it does not
//!
//! It buys `StrictMath.pow` exactly, which is a stronger claim than closeness: `StrictMath` is
//! *specified* to be FDLIBM, so a divergence against it would mean this port is wrong rather than
//! that the algorithms differ. The corpus checks that on every point instead of taking the
//! specification's word for it.
//!
//! It does not buy `Math.pow`, and cannot. `Math` is free to use an intrinsic and on x86-64 HotSpot
//! does, so every call site reaching `Math.pow` stays unreproducible bit for bit. What changes is
//! that the gap is measured against the best permissively licensed implementation rather than
//! against whatever libm the host happens to ship.
//!
//! # The algorithm, as FDLIBM states it
//!
//! ```text
//!                   n
//! Method:  Let x = 2  * (1+f)
//!   1. Compute and return log2(x) in two pieces:
//!         log2(x) = w1 + w2,
//!      where w1 has 53-24 = 29 bit trailing zeros.
//!   2. Perform y*log2(x) = n+y' by simulating muti-precision arithmetic, where |y'| <= 0.5.
//!   3. Return x**y = 2**n * exp(y'*log2)
//! ```
//!
//! The two-piece split is the whole point, and it is why nothing here may be simplified. `__LO(t1)
//! = 0` clears the low word of a double so the high part carries exactly 29 trailing zero bits;
//! the tail is then recovered as a separate term. Rewriting `t2 = v-(t1-u)` into anything
//! algebraically equal destroys the split and the result changes.
//!
//! # The nineteen special cases are the specification
//!
//! FDLIBM enumerates them, and they are not all what an implementation would produce by accident:
//! `(anything)**0` is `1` including `NaN**0`; `(+-1)**(+-Inf)` is `NaN` rather than `1`;
//! `(-anything)**(non-integer)` is `NaN` while `(-anything)**(integer)` carries the sign of
//! `(-1)**integer`. Each is transcribed rather than derived, and [`pow`]'s tests walk them.

// The constants below are FDLIBM's decimal literals, transcribed. Several carry more digits than
// are needed to name their double, and clippy is right that the tail is redundant — but the
// redundant digits are the source's, not this port's, and shortening them would make the file a
// paraphrase of fdlibm rather than a transcription of it. The values are identical either way; the
// provenance is not.
#![allow(clippy::excessive_precision)]

/// The high 32 bits of a double, as the C's signed `int`.
fn high(x: f64) -> i32 {
    (x.to_bits() >> 32) as u32 as i32
}

/// The low 32 bits, as the C's `unsigned`.
fn low(x: f64) -> u32 {
    x.to_bits() as u32
}

/// `__HI(x) = value`.
fn with_high(x: f64, value: i32) -> f64 {
    f64::from_bits((u64::from(value as u32) << 32) | (x.to_bits() & 0xffff_ffff))
}

/// `__LO(x) = 0`, which is how FDLIBM takes the leading half of a double.
fn clear_low(x: f64) -> f64 {
    f64::from_bits(x.to_bits() & 0xffff_ffff_0000_0000)
}

/// `2**k`, exactly, for `k` in the normal exponent range.
fn two_to(k: i32) -> f64 {
    f64::from_bits(((0x3ff + k) as u64) << 52)
}

/// `scalbn(x, n)`: `x * 2**n`, with the scaling split so only the final multiply can round.
///
/// Reached from one place, the subnormal-output branch at the end of [`pow`]. A single
/// `x * 2**n` would be wrong there because `2**n` is not representable.
fn scalbn(mut x: f64, mut n: i32) -> f64 {
    if n > 1023 {
        x *= two_to(1023);
        n -= 1023;
        if n > 1023 {
            x *= two_to(1023);
            n -= 1023;
            if n > 1023 {
                n = 1023;
            }
        }
    } else if n < -1022 {
        // 2**-969 keeps the intermediate normal; two steps cover the whole subnormal range.
        x *= two_to(-1022) * two_to(53);
        n += 969;
        if n < -1022 {
            x *= two_to(-1022) * two_to(53);
            n += 969;
            if n < -1022 {
                n = -1022;
            }
        }
    }
    x * two_to(n)
}

const BP: [f64; 2] = [1.0, 1.5];
/// `dp_h[1]`, the leading half of `log2(1.5)`, with its low word already zero.
const DP_H: [f64; 2] = [0.0, 5.849_624_872_207_641_6e-1];
/// `dp_l[1]`, the tail the split leaves behind.
const DP_L: [f64; 2] = [0.0, 1.350_039_202_129_748_9e-8];
const TWO53: f64 = 9_007_199_254_740_992.0;
/// `huge` and `tiny` exist to raise the overflow and underflow flags; their products are the
/// returned infinity and zero.
const HUGE: f64 = 1.0e300;
const TINY: f64 = 1.0e-300;
/// `L1` to `L6`: the polynomial for `(3/2)*(log(x)-2s-2/3*s**3)`.
const L1: f64 = 5.999_999_999_999_946_5e-1;
const L2: f64 = 4.285_714_285_785_502e-1;
const L3: f64 = 3.333_333_298_183_774_3e-1;
const L4: f64 = 2.727_281_238_085_340_1e-1;
const L5: f64 = 2.306_607_457_755_617_6e-1;
const L6: f64 = 2.069_750_178_003_384_2e-1;
/// `P1` to `P5`: the polynomial of step 3, not the same coefficients as `exp`'s.
const P1: f64 = 1.666_666_666_666_660_2e-1;
const P2: f64 = -2.777_777_777_701_559_3e-3;
const P3: f64 = 6.613_756_321_437_934e-5;
const P4: f64 = -1.653_390_220_546_525_2e-6;
const P5: f64 = 4.138_136_797_057_238_5e-8;
#[allow(clippy::approx_constant)]
const LG2: f64 = 6.931_471_805_599_453e-1;
/// `lg2_h`, `ln2` truncated to its leading half, and `lg2_l`, the tail.
const LG2_H: f64 = 6.931_471_824_645_996e-1;
const LG2_L: f64 = -1.904_654_299_957_768e-9;
/// `ovt`: `-(1024-log2(ovfl+.5ulp))`, the overflow threshold's tail.
const OVT: f64 = 8.008_566_259_537_294e-17;
/// `cp` = `2/(3*ln2)`, with its leading half `cp_h` and tail `cp_l`.
const CP: f64 = 9.617_966_939_259_756e-1;
const CP_H: f64 = 9.617_967_009_544_373e-1;
const CP_L: f64 = -7.028_461_650_952_758e-9;
/// `ivln2` = `1/ln2`, split the same way.
#[allow(clippy::approx_constant)]
const IVLN2: f64 = 1.442_695_040_888_963_4;
#[allow(clippy::approx_constant)]
const IVLN2_H: f64 = 1.442_695_021_629_333_5;
const IVLN2_L: f64 = 1.925_962_991_126_617_5e-8;

/// `StrictMath.pow(x, y)`, which is FDLIBM's `__ieee754_pow`.
///
/// The structure is the C function's, including the order of every floating-point operation. The
/// two-piece splits (`t1`/`t2`, `p_h`/`p_l`, `z_h`/`z_l`) are what make the result reproducible,
/// and each one depends on a low word having been cleared at exactly the point the C clears it.
///
/// Two lint families are silenced here rather than worked around. `eq_op` fires on `(x-x)/(x-x)`
/// and `z-z`, which are how FDLIBM produces a NaN *and* raises the invalid flag; writing
/// `f64::NAN` instead would be a different program in C and is kept faithful here. `approx_constant`
/// fires on `lg2`, which is `ln 2` — the same refusal [`crate::strict_exp`] makes, for the same
/// reason: the value is fdlibm's literal, and substituting one from elsewhere would make the port
/// depend on two sources agreeing to the last bit rather than on the one it cites.
#[allow(
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::eq_op,
    clippy::approx_constant
)]
pub fn pow(x: f64, y: f64) -> f64 {
    let (hx, lx) = (high(x), low(x));
    let (hy, ly) = (high(y), low(y));
    let ix = hx & 0x7fff_ffff;
    let iy = hy & 0x7fff_ffff;

    // y == 0: x**0 is 1 for every x, NaN included. Case 1 of the nineteen.
    if (iy as u32 | ly) == 0 {
        return 1.0;
    }

    // Either argument NaN: return x+y, which quiets a signalling NaN and picks up whichever is
    // NaN. This runs *after* the y == 0 test, which is why NaN**0 is 1.
    if ix > 0x7ff0_0000
        || (ix == 0x7ff0_0000 && lx != 0)
        || iy > 0x7ff0_0000
        || (iy == 0x7ff0_0000 && ly != 0)
    {
        return x + y;
    }

    // Whether y is an integer, and if so its parity, which only matters when x is negative:
    //   0 = not an integer, 1 = odd integer, 2 = even integer.
    let mut yisint = 0i32;
    if hx < 0 {
        if iy >= 0x4340_0000 {
            // |y| >= 2**53: every such double is an even integer.
            yisint = 2;
        } else if iy >= 0x3ff0_0000 {
            let k = (iy >> 20) - 0x3ff;
            if k > 20 {
                let j = ly >> (52 - k);
                if (j << (52 - k)) == ly {
                    yisint = 2 - ((j & 1) as i32);
                }
            } else if ly == 0 {
                let j = iy >> (20 - k);
                if (j << (20 - k)) == iy {
                    yisint = 2 - (j & 1);
                }
            }
        }
    }

    // The special values of y, cases 3 to 9 plus the two shortcuts.
    if ly == 0 {
        if iy == 0x7ff0_0000 {
            // y is +-Inf.
            if ((ix as u32).wrapping_sub(0x3ff0_0000) | lx) == 0 {
                // (+-1)**(+-Inf) is NaN, not 1. Case 9.
                return y - y;
            } else if ix >= 0x3ff0_0000 {
                // (|x|>1)**(+-Inf) is Inf or 0.
                return if hy >= 0 { y } else { 0.0 };
            }
            // (|x|<1)**(-Inf, +Inf) is Inf, 0.
            return if hy < 0 { -y } else { 0.0 };
        }
        if iy == 0x3ff0_0000 {
            // y is +-1.
            return if hy < 0 { 1.0 / x } else { x };
        }
        if hy == 0x4000_0000 {
            // y is 2. Exact, and the reason `Histogram.getStandardDeviation` never needed this
            // function ported at all.
            return x * x;
        }
        if hy == 0x3fe0_0000 && hx >= 0 {
            // y is 0.5 and x >= +0. `sqrt` is correctly rounded, so this is exact.
            return x.sqrt();
        }
    }

    let mut ax = x.abs();

    // The special values of x: +-0, +-Inf, +-1.
    if lx == 0 && (ix == 0x7ff0_0000 || ix == 0 || ix == 0x3ff0_0000) {
        let mut z = ax;
        if hy < 0 {
            z = 1.0 / z;
        }
        if hx < 0 {
            if ((ix - 0x3ff0_0000) | yisint) == 0 {
                // (-1)**(non-integer) is NaN.
                z = (z - z) / (z - z);
            } else if yisint == 1 {
                // (x<0)**(odd) = -(|x|**odd).
                z = -z;
            }
        }
        return z;
    }

    // `n` is 0 for x >= 0 and 1 for x < 0, from the arithmetic shift of the sign bit.
    let mut n = (hx >> 31) + 1;

    // (x<0)**(non-integer) is NaN. Case 19.
    if (n | yisint) == 0 {
        return (x - x) / (x - x);
    }

    // `s` carries the sign of the result: -1 only for a negative base raised to an odd integer.
    let mut s = 1.0f64;
    if (n | (yisint - 1)) == 0 {
        s = -1.0;
    }

    let (t1, t2);
    if iy > 0x41e0_0000 {
        // |y| > 2**31: the result overflows or underflows unless x is extremely close to 1.
        if iy > 0x43f0_0000 {
            // |y| > 2**64 must over/underflow.
            if ix <= 0x3fef_ffff {
                return if hy < 0 { HUGE * HUGE } else { TINY * TINY };
            }
            if ix >= 0x3ff0_0000 {
                return if hy > 0 { HUGE * HUGE } else { TINY * TINY };
            }
        }
        if ix < 0x3fef_ffff {
            return if hy < 0 {
                s * HUGE * HUGE
            } else {
                s * TINY * TINY
            };
        }
        if ix > 0x3ff0_0000 {
            return if hy > 0 {
                s * HUGE * HUGE
            } else {
                s * TINY * TINY
            };
        }
        // |1-x| is tiny, so log(x) is computed by its series rather than by the general path.
        let t = ax - 1.0; // t has 20 trailing zeros
        let w = (t * t) * (0.5 - t * (0.333_333_333_333_333_33 - t * 0.25));
        let u = IVLN2_H * t; // ivln2_h has 21 significant bits
        let v = t * IVLN2_L - w * IVLN2;
        let first = u + v;
        t1 = clear_low(first);
        t2 = v - (t1 - u);
    } else {
        let mut n_local = 0i32;
        let mut ix = ix;
        // Take care of a subnormal x by scaling it up and paying for it in the exponent.
        if ix < 0x0010_0000 {
            ax *= TWO53;
            n_local -= 53;
            ix = high(ax);
        }
        n_local += (ix >> 20) - 0x3ff;
        let j = ix & 0x000f_ffff;
        // Normalise `ix` into [1, 2) and choose which of the two intervals x falls in.
        ix = j | 0x3ff0_0000;
        let k;
        if j <= 0x0003_988E {
            // |x| < sqrt(3/2)
            k = 0;
        } else if j < 0x000B_B67A {
            // |x| < sqrt(3)
            k = 1;
        } else {
            k = 0;
            n_local += 1;
            ix -= 0x0010_0000;
        }
        ax = with_high(ax, ix);

        // ss = s_h + s_l = (x-1)/(x+1) or (x-1.5)/(x+1.5).
        let u = ax - BP[k];
        let v = 1.0 / (ax + BP[k]);
        let ss = u * v;
        let s_h = clear_low(ss);
        // t_h is the high half of ax+bp[k], built directly in the exponent field.
        let t_h = with_high(
            0.0,
            ((ix >> 1) | 0x2000_0000) + 0x0008_0000 + ((k as i32) << 18),
        );
        let t_l = ax - (t_h - BP[k]);
        let s_l = v * ((u - s_h * t_h) - s_h * t_l);

        // log(ax), as a polynomial in ss**2.
        let s2 = ss * ss;
        let mut r = s2 * s2 * (L1 + s2 * (L2 + s2 * (L3 + s2 * (L4 + s2 * (L5 + s2 * L6)))));
        r += s_l * (s_h + ss);
        let s2 = s_h * s_h;
        let t_h = clear_low(3.0 + s2 + r);
        let t_l = r - ((t_h - 3.0) - s2);

        // u+v = ss*(1+...)
        let u = s_h * t_h;
        let v = s_l * t_h + t_l * ss;
        let p_h = clear_low(u + v);
        let p_l = v - (p_h - u);
        let z_h = CP_H * p_h;
        let z_l = CP_L * p_h + p_l * CP + DP_L[k];

        // log2(ax) = n + dp_h[k] + z_h + z_l, split into a head with 29 zero low bits and a tail.
        let t = f64::from(n_local);
        let first = ((z_h + z_l) + DP_H[k]) + t;
        t1 = clear_low(first);
        t2 = z_l - (((t1 - t) - DP_H[k]) - z_h);
    }

    // Split y into y1 + (y - y1) and form (y1+y2)*(t1+t2).
    let y1 = clear_low(y);
    let p_l = (y - y1) * t1 + y * t2;
    let mut p_h = y1 * t1;
    let z = p_l + p_h;
    let j = high(z);
    let i = low(z) as i32;

    if j >= 0x4090_0000 {
        // z >= 1024
        if (j.wrapping_sub(0x4090_0000) | i) != 0 {
            return s * HUGE * HUGE;
        }
        if p_l + OVT > z - p_h {
            return s * HUGE * HUGE;
        }
    } else if (j & 0x7fff_ffff) >= 0x4090_cc00 {
        // z <= -1075. The comparison the C writes as an unsigned subtraction, kept unsigned here
        // because `0xc090cc00` does not fit in a signed int and the C promotes.
        if ((j as u32).wrapping_sub(0xc090_cc00) | (i as u32)) != 0 {
            return s * TINY * TINY;
        }
        if p_l <= z - p_h {
            return s * TINY * TINY;
        }
    }

    // Compute 2**(p_h+p_l).
    let i = j & 0x7fff_ffff;
    let mut k = (i >> 20) - 0x3ff;
    n = 0;
    if i > 0x3fe0_0000 {
        // |z| > 0.5, so set n = [z+0.5] and reduce p_h by it.
        n = j + (0x0010_0000 >> (k + 1));
        k = ((n & 0x7fff_ffff) >> 20) - 0x3ff;
        let t = with_high(0.0, n & !(0x000f_ffff >> k));
        n = ((n & 0x000f_ffff) | 0x0010_0000) >> (20 - k);
        if j < 0 {
            n = -n;
        }
        p_h -= t;
    }
    let t = clear_low(p_l + p_h);
    let u = t * LG2_H;
    let v = (p_l - (t - p_h)) * LG2 + t * LG2_L;
    let mut z = u + v;
    let w = v - (z - u);
    let t = z * z;
    let t1 = z - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
    let r = (z * t1) / (t1 - 2.0) - (w + z * w);
    z = 1.0 - (r - z);
    let j = high(z) + (n << 20);
    if (j >> 20) <= 0 {
        // The result is subnormal, so the exponent cannot simply be added into the field.
        z = scalbn(z, n);
    } else {
        z = with_high(z, j);
    }
    s * z
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The nineteen special cases FDLIBM enumerates, in its order. Several are not what an
    /// implementation would produce by accident.
    #[test]
    fn the_special_cases_are_the_specification() {
        // 1. anything ** 0 is 1, NaN included.
        assert_eq!(pow(2.0, 0.0), 1.0);
        assert_eq!(pow(f64::NAN, 0.0), 1.0);
        assert_eq!(pow(f64::INFINITY, 0.0), 1.0);
        assert_eq!(pow(0.0, 0.0), 1.0);
        // 2. anything ** 1 is itself.
        assert_eq!(pow(-3.5, 1.0), -3.5);
        // 3, 4. NaN propagates otherwise.
        assert!(pow(2.0, f64::NAN).is_nan());
        assert!(pow(f64::NAN, 2.0).is_nan());
        // 5 to 8. |x| either side of 1, raised to +-Inf.
        assert_eq!(pow(2.0, f64::INFINITY), f64::INFINITY);
        assert_eq!(pow(2.0, f64::NEG_INFINITY), 0.0);
        assert_eq!(pow(0.5, f64::INFINITY), 0.0);
        assert_eq!(pow(0.5, f64::NEG_INFINITY), f64::INFINITY);
        // 9. +-1 ** +-Inf is NaN, not 1.
        assert!(pow(1.0, f64::INFINITY).is_nan());
        assert!(pow(-1.0, f64::NEG_INFINITY).is_nan());
        // 10 to 14. The zeros, where the parity of an integer exponent carries the sign.
        assert_eq!(pow(0.0, 3.0), 0.0);
        assert_eq!(pow(-0.0, 2.0), 0.0);
        assert_eq!(pow(-0.0, 3.0), -0.0);
        assert_eq!(pow(0.0, -3.0), f64::INFINITY);
        assert_eq!(pow(-0.0, -3.0), f64::NEG_INFINITY);
        // 15 to 17. The infinities.
        assert_eq!(pow(f64::INFINITY, 2.0), f64::INFINITY);
        assert_eq!(pow(f64::INFINITY, -2.0), 0.0);
        assert_eq!(pow(f64::NEG_INFINITY, 3.0), f64::NEG_INFINITY);
        // 18, 19. A negative base is fine for an integer exponent and NaN otherwise.
        assert_eq!(pow(-2.0, 3.0), -8.0);
        assert_eq!(pow(-2.0, 2.0), 4.0);
        assert!(pow(-2.0, 0.5).is_nan());
    }

    /// `pow(integer, integer)` returns the exact integer when it is representable, which FDLIBM
    /// states as an accuracy guarantee rather than as a side effect.
    #[test]
    fn integer_powers_are_exact() {
        for base in 2..=10i32 {
            let mut expected = 1.0f64;
            for exponent in 0..=10 {
                assert_eq!(
                    pow(f64::from(base), f64::from(exponent)),
                    expected,
                    "{base}**{exponent}"
                );
                expected *= f64::from(base);
            }
        }
        assert_eq!(pow(10.0, 22.0), 1e22);
    }

    /// The shortcuts the special-value block takes, which must agree with the general path's own
    /// answer rather than merely being fast.
    #[test]
    fn the_shortcuts_agree_with_the_arithmetic() {
        for x in [0.5, 1.5, 3.0, 7.25, 1e-8, 1e8] {
            // y == 2 returns x*x directly.
            assert_eq!(pow(x, 2.0), x * x);
            // y == 0.5 returns sqrt, which is correctly rounded.
            assert_eq!(pow(x, 0.5), x.sqrt());
            // y == -1 returns the reciprocal.
            assert_eq!(pow(x, -1.0), 1.0 / x);
        }
    }

    #[test]
    fn overflow_and_underflow_saturate() {
        assert_eq!(pow(10.0, 400.0), f64::INFINITY);
        assert_eq!(pow(10.0, -400.0), 0.0);
        assert_eq!(pow(2.0, 1024.0), f64::INFINITY);
        // A subnormal result, which is the only path that reaches `scalbn`.
        assert!(pow(2.0, -1060.0) > 0.0);
        assert!(pow(2.0, -1060.0) < f64::MIN_POSITIVE);
    }

    /// A subnormal base, which the general path scales by 2**53 before taking its logarithm.
    #[test]
    fn a_subnormal_base_is_scaled_rather_than_flushed() {
        let subnormal = f64::from_bits(1);
        assert!(pow(subnormal, 1.0) == subnormal);
        assert!(pow(subnormal, 0.5) > 0.0);
        assert_eq!(pow(subnormal, 2.0), 0.0);
    }
}
