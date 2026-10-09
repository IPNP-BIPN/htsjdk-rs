//! Correctly-rounded `log` and `log10`, matching `java.lang.Math`.
//!
//! `Math.log` and `Math.log10` were measured to be **correctly rounded** on every point of the
//! conformance corpus, so the target here is not "reproduce HotSpot's algorithm" but simply
//! "round the true result once". That is a much smaller job than porting the intrinsic, and it
//! is why these two land before `exp` and `pow`, which are *not* correctly rounded and do need
//! the algorithm ported.
//!
//! See `docs/decisions/0006-correct-rounding-is-the-target-for-log-and-log10.md`.
//!
//! `clippy::approx_constant` is allowed module-wide. It flags the `hi` halves of the
//! double-double constants as approximations of `std::f64::consts::LN_2` and `LOG10_E`, and
//! taking that advice would defeat the entire point: these are pairs precisely so they carry
//! ~53 bits *more* than a single `f64`, and `hi` alone is deliberately only the leading half.
//! Substituting the std constant silently discards the precision that makes correct rounding
//! possible, and the corpus would fail on exactly the hard-to-round points.
//!
//! `clippy::excessive_precision` is allowed for the same reason: the pairs were emitted at
//! 400-bit precision and the digits record what was generated, not what `f64` can hold.
#![allow(clippy::approx_constant, clippy::excessive_precision)]

use crate::dd::{self, DoubleDouble};

// ln(2) and log10(e) to ~106 bits, generated at 400-bit precision.
const LN2: DoubleDouble = DoubleDouble::new(6.93147180559945286e-01, 2.31904681384629956e-17);
const LOG10_E: DoubleDouble = DoubleDouble::new(4.34294481903251817e-01, 1.09831965021676507e-17);

/// `ln(x)` in double-double, for finite positive normal-or-subnormal `x`.
///
/// Argument reduction puts the mantissa in `[sqrt(1/2), sqrt(2))` so that `s = (m-1)/(m+1)`
/// satisfies `|s| <= 0.1716`. The atanh series in `s^2` then converges by a factor of ~0.029
/// per term, so 22 terms carry it past 106 bits.
fn ln_dd(x: f64) -> DoubleDouble {
    // Decompose x = m * 2^e with m in [1, 2).
    let mut e = ((x.to_bits() >> 52) & 0x7ff) as i64 - 1023;
    let mut m = f64::from_bits((x.to_bits() & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    if x.to_bits() >> 52 & 0x7ff == 0 {
        // Subnormal: scale into the normal range and correct the exponent afterwards.
        let scaled = x * f64::from_bits(0x4350_0000_0000_0000); // 2^54
        e = ((scaled.to_bits() >> 52) & 0x7ff) as i64 - 1023 - 54;
        m = f64::from_bits((scaled.to_bits() & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    }

    // Centre the mantissa on 1 to shrink |s|.
    if m > std::f64::consts::SQRT_2 {
        m *= 0.5;
        e += 1;
    }

    // s = (m - 1) / (m + 1), computed in double-double because the subtraction cancels.
    let num = dd::add_f64(DoubleDouble::from_f64(m), -1.0);
    let den = dd::add_f64(DoubleDouble::from_f64(m), 1.0);
    let s = dd::div(num, den);
    let s2 = dd::mul(s, s);

    // atanh(s) = s + s^3/3 + s^5/5 + ...  and  ln(m) = 2 * atanh(s)
    //
    // The coefficients are divided, not multiplied by a precomputed reciprocal. `1.0/3.0` is
    // not representable, so multiplying by it injects a relative error of about 2^-53 into the
    // largest term, which lands near 2^-60 overall. That is precisely the range where the
    // hard-to-round corpus cases live (they need 62 to 64 bits to resolve), so the reciprocal
    // form fails exactly the points this function exists to get right, and passes everywhere
    // else. Dividing by the exactly-representable odd integer keeps full precision.
    let mut sum = DoubleDouble::from_f64(1.0);
    let mut term = DoubleDouble::from_f64(1.0);
    for k in 1..=22 {
        term = dd::mul(term, s2);
        let denom = DoubleDouble::from_f64((2 * k + 1) as f64);
        sum = dd::add(sum, dd::div(term, denom));
    }
    let ln_m = dd::mul_f64(dd::mul(s, sum), 2.0);

    dd::add(ln_m, dd::mul_f64(LN2, e as f64))
}

/// Intervals the fast phase splits `[1, 2)` into, by the mantissa's leading bits.
const TABLE_BITS: u32 = 7;

/// Per interval `i`: `r_i`, an `f64` near `1 / c_i` where `c_i` is the interval's centre (exactly
/// `1.0` for the interval that starts at one, so a value just above one is not cancelled against a
/// table entry), and `-ln(r_i)` in double-double, computed once by [`ln_dd`] itself.
fn table() -> &'static [(f64, DoubleDouble); 1 << TABLE_BITS] {
    static TABLE: std::sync::OnceLock<[(f64, DoubleDouble); 1 << TABLE_BITS]> =
        std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        let mut table = [(1.0, DoubleDouble::from_f64(0.0)); 1 << TABLE_BITS];
        for (i, entry) in table.iter_mut().enumerate().skip(1) {
            let centre = 1.0 + (i as f64 + 0.5) / f64::from(1u32 << TABLE_BITS);
            let r = 1.0 / centre;
            let ln_r = ln_dd(r);
            *entry = (r, DoubleDouble::new(-ln_r.hi, -ln_r.lo));
        }
        table
    })
}

/// `ln(x)` in double-double with an error bound, for a positive normal `x`: the fast phase.
///
/// `x = 2^e * m`, `m` in `[1, 2)`; with `r` the table's entry for `m`'s interval, `z = m * r - 1`
/// is formed EXACTLY (the product by `two_prod`, the subtraction by Sterbenz, since `m * r` is
/// within 2^-7 of one), so `|z| < 2^-7` carries no error. Then
/// `ln(x) = e * ln(2) - ln(r) + ln(1 + z)`, the first two in double-double from constants accurate
/// to about 2^-100, and `ln(1 + z) = z - z^2/2 + z^3 * tail(z)`, the first two terms in
/// double-double and the tail (`1/3 - z/4 + ... - z^7/10`) in plain `f64`.
///
/// The bound returned is far wider than the error it covers. The truncated series leaves
/// `|z|^11/11 < 2^-80 |z|`; the tail's coefficients and Horner evaluation are relative errors of
/// a few 2^-53 on a value below `|z|^3 / 3 < 2^-22 |z|`, under 2^-70 |z| together, and using
/// `z.hi` for `z` in the tail moves it by `2^-53` relative; the double-double sums and constants
/// contribute about 2^-100 of the larger magnitudes. The bound claims 2^-60 |z| and 2^-90 of the
/// rest.
fn ln_fast(x: f64) -> (DoubleDouble, f64) {
    let bits = x.to_bits();
    let (e, minus_ln_r, z) = if (x - 1.0).abs() < 1.0 / f64::from(1u32 << TABLE_BITS) {
        // Near one, on either side, `z = x - 1` exactly (Sterbenz) and nothing else: below one
        // the decomposition would be `2^-1 * m` with `m` near two, and the result, of the order
        // of `z`, would be the difference of two terms near `ln(2)`.
        (
            0,
            DoubleDouble::from_f64(0.0),
            DoubleDouble::from_f64(x - 1.0),
        )
    } else {
        let e = ((bits >> 52) & 0x7ff) as i64 - 1023;
        let m = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
        let index = ((bits >> (52 - TABLE_BITS)) & ((1 << TABLE_BITS) - 1)) as usize;
        let (r, minus_ln_r) = table()[index];
        let (p, p_err) = dd::two_prod(m, r);
        let (z_hi, z_lo) = dd::two_sum(p - 1.0, p_err);
        (e, minus_ln_r, DoubleDouble::new(z_hi, z_lo))
    };
    let z_hi = z.hi;

    let mut tail = -1.0 / 10.0;
    for k in (3..=9).rev() {
        let coefficient = if k % 2 == 1 { 1.0 } else { -1.0 } / k as f64;
        tail = tail * z_hi + coefficient;
    }
    let z3_tail = z_hi * z_hi * z_hi * tail;
    let half_z2 = dd::mul_f64(dd::mul(z, z), 0.5);
    let ln_1pz = dd::add_f64(dd::sub(z, half_z2), z3_tail);

    let head = dd::add(dd::mul_f64(LN2, e as f64), minus_ln_r);
    let result = dd::add(head, ln_1pz);
    let bound = z_hi.abs() * f64::from_bits(0x3c30_0000_0000_0000) // 2^-60
        + (head.hi.abs() + z_hi.abs()) * f64::from_bits(0x3a50_0000_0000_0000); // 2^-90
    (result, bound)
}

/// The nearest `f64` to a value known to lie within `bound` of `value`, if that interval does not
/// straddle a rounding boundary, so the answer is the correctly rounded one whatever the value
/// is; `None` when it might, and the slow phase has to decide.
///
/// The interval tested is twice the bound. That is what makes the answer [`ln_dd`]'s and not only
/// the true value's rounding: the true value is within `bound` of `value`, `ln_dd`'s is within
/// about 2^-100 of the true value, so both lie inside `2 * bound`, and every point of it rounds to
/// the same `f64`.
fn round_if_certain(value: DoubleDouble, bound: f64) -> Option<f64> {
    let down = value.hi + (value.lo - 2.0 * bound);
    let up = value.hi + (value.lo + 2.0 * bound);
    (down == up).then_some(down)
}

/// `ln_dd(x).to_f64()`, by way of the fast phase whenever its answer is certain.
///
/// The double-double series costs 22 divisions per call, and the speed baseline found it under
/// 94% of `ModelSegments` and most of `VariantRecalibrator` (Milestone S, #112). The fast phase is
/// a table, an exact reduction and a short polynomial; [`round_if_certain`] returns its answer
/// only where that answer is provably the one `ln_dd` rounds to, and falls back otherwise, so the
/// function's output is unchanged by construction.
fn ln_rounded(x: f64) -> f64 {
    if x.to_bits() >> 52 & 0x7ff != 0 {
        let (value, bound) = ln_fast(x);
        if let Some(answer) = round_if_certain(value, bound) {
            return answer;
        }
    }
    ln_dd(x).to_f64()
}

/// `java.lang.Math.log`, correctly rounded.
pub fn log(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    if x == 1.0 {
        return 0.0;
    }
    ln_rounded(x)
}

/// `java.lang.Math.log10`, correctly rounded.
pub fn log10(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    if x == 1.0 {
        return 0.0;
    }

    // Exact powers of ten must return exact integers. The series route would land within one
    // ulp but not necessarily *on* the integer, and log10 of a power of ten showing up as
    // 2.9999999999999996 is the kind of thing that survives all the way into a report column.
    if x.fract() == 0.0 && x > 0.0 && x <= 1e22 {
        let mut p = 1.0f64;
        for k in 0..=22 {
            if p == x {
                return k as f64;
            }
            p *= 10.0;
        }
    }

    if x.to_bits() >> 52 & 0x7ff != 0 {
        let (value, bound) = ln_fast(x);
        // `LOG10_E` is below one and carries about 2^-106 of its own, so the bound scales down
        // and gains a term for the product.
        let scaled = dd::mul(value, LOG10_E);
        if let Some(answer) = round_if_certain(
            scaled,
            bound + scaled.hi.abs() * f64::from_bits(0x3a50_0000_0000_0000),
        ) {
            return answer;
        }
    }
    dd::mul(ln_dd(x), LOG10_E).to_f64()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A xorshift stream, so the sample is the same on every run and every host.
    fn stream(seed: u64) -> impl FnMut() -> u64 {
        let mut state = seed;
        move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        }
    }

    /// Positive normal doubles from every region the fast phase treats differently: any
    /// exponent, just either side of one, either side of every table boundary, and either side
    /// of every power of two.
    fn inputs(count: usize) -> Vec<f64> {
        let mut next = stream(0x9e37_79b9_7f4a_7c15);
        let mut out = Vec::with_capacity(count);
        while out.len() < count {
            let raw = next();
            let x = match raw % 4 {
                0 => f64::from_bits(raw >> 2 | 0x0010_0000_0000_0000) % f64::MAX,
                1 => {
                    1.0 + ((raw >> 12) as f64 / (1u64 << 52) as f64 - 0.5)
                        * 2f64.powi(-((raw % 40) as i32))
                }
                2 => {
                    let boundary = 1.0 + ((raw >> 8) % 128) as f64 / 128.0;
                    let offset = ((raw >> 20) % 2048) as f64 - 1024.0;
                    (boundary + offset * f64::EPSILON) * 2f64.powi((raw >> 40) as i32 % 64 - 32)
                }
                _ => {
                    let power = 2f64.powi((raw >> 8) as i32 % 2000 - 1000);
                    f64::from_bits(
                        power
                            .to_bits()
                            .wrapping_add((raw >> 32) % 4096)
                            .wrapping_sub(2048),
                    )
                }
            };
            if x.is_finite() && x > 0.0 && x.to_bits() >> 52 & 0x7ff != 0 {
                out.push(x);
            }
        }
        out
    }

    #[test]
    fn the_fast_phase_answers_exactly_what_the_series_does() {
        let mut fast = 0usize;
        let sample = inputs(2_000_000);
        for &x in &sample {
            let (value, bound) = ln_fast(x);
            if let Some(answer) = round_if_certain(value, bound) {
                fast += 1;
                assert_eq!(answer.to_bits(), ln_dd(x).to_f64().to_bits(), "log({x:e})");
            }
            let scaled = dd::mul(value, LOG10_E);
            let bound10 = bound + scaled.hi.abs() * f64::from_bits(0x3a50_0000_0000_0000);
            if let Some(answer) = round_if_certain(scaled, bound10) {
                let slow = dd::mul(ln_dd(x), LOG10_E).to_f64();
                assert_eq!(answer.to_bits(), slow.to_bits(), "log10({x:e})");
            }
        }
        // The fallback is for the points the fast phase cannot decide; if it were most of them,
        // the fast phase would be a cost rather than a saving.
        assert!(
            fast * 100 > sample.len() * 99,
            "fast phase decided only {fast} of {}",
            sample.len()
        );
    }

    #[test]
    fn the_error_bound_holds_against_the_series() {
        // Stronger than agreement after rounding: the fast phase's double-double lies within its
        // claimed bound of the series' value, which is about 2^-100 from the true one.
        for x in inputs(500_000) {
            let (value, bound) = ln_fast(x);
            let reference = ln_dd(x);
            let difference = dd::sub(value, reference);
            assert!(
                difference.hi.abs() <= bound,
                "log({x:e}): off by {:e}, bound {bound:e}",
                difference.hi
            );
        }
    }
}

#[cfg(test)]
mod long_tests {
    use super::*;

    /// The agreement test at a hundred times the size, run by hand (`--ignored`) and recorded
    /// in the commit that introduced the fast phase, rather than on every CI run. It also
    /// reports the share the fast phase decided and the time per call of each phase.
    #[test]
    #[ignore]
    fn two_hundred_million_points_agree() {
        let sample_size = 200_000_000;
        let mut next = tests_stream();
        let (mut fast, mut checked) = (0u64, 0u64);
        while checked < sample_size {
            let raw = next();
            let x = f64::from_bits(raw >> 1);
            if !(x.is_finite() && x > 0.0 && x.to_bits() >> 52 & 0x7ff != 0) {
                continue;
            }
            // Half the sample within 2^-7 of one, where the fast phase works hardest.
            let x = if raw & 1 == 0 {
                x
            } else {
                1.0 + (x.fract() - 0.5) / 64.0
            };
            checked += 1;
            let (value, bound) = ln_fast(x);
            if let Some(answer) = round_if_certain(value, bound) {
                fast += 1;
                assert_eq!(answer.to_bits(), ln_dd(x).to_f64().to_bits(), "log({x:e})");
            }
        }
        eprintln!("{checked} points, {fast} decided by the fast phase, none differing");
        let xs: Vec<f64> = (1..=1_000_000).map(|i| f64::from(i) * 0.731).collect();
        let start = std::time::Instant::now();
        let slow: f64 = xs.iter().map(|&x| ln_dd(x).to_f64()).sum();
        let slow_ns = start.elapsed().as_nanos() as f64 / xs.len() as f64;
        let start = std::time::Instant::now();
        let quick: f64 = xs.iter().map(|&x| ln_rounded(x)).sum();
        let quick_ns = start.elapsed().as_nanos() as f64 / xs.len() as f64;
        assert_eq!(slow.to_bits(), quick.to_bits());
        eprintln!("ln_dd {slow_ns:.1} ns per call, ln_rounded {quick_ns:.1} ns per call");
    }

    fn tests_stream() -> impl FnMut() -> u64 {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        }
    }
}
