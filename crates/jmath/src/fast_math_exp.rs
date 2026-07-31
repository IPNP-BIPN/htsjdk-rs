//! `FastMath.exp` and the tables under it, ported from `org.apache.commons.math3.util.FastMath`
//! and `FastMathCalc` (commons-math3 3.5), which are Apache 2.0. See
//! `docs/decisions/0023-commons-math3-is-portable-where-the-jdk-is-not.md`.
//!
//! This is the exponential GATK's `MathUtils` and commons-math3's own `Gamma`, `Erf` and
//! `NormalDistribution` reach. It is **not** [`crate::math::exp`], which is `java.lang.Math.exp`,
//! a HotSpot intrinsic whose port was withdrawn in decision 0014 for being a GPL2 transcription.
//! Two exponentials, two licences, two different answers on the same input: a call site that
//! names `FastMath` must reach this one.
//!
//! # The tables are computed, not transcribed
//!
//! ```java
//! private static final boolean RECOMPUTE_TABLES_AT_RUNTIME = false;
//! ...
//! EXP_INT_TABLE_A = FastMathLiteralArrays.loadExpIntA();
//! ```
//!
//! The reference ships 6,175 lines of literal doubles and can regenerate them with `FastMathCalc`
//! instead; the flag choosing between the two is a constant. This port takes the **computing**
//! branch, because 3,550 transcribed literals are 3,550 chances to typo one and no reader would
//! catch it. That makes the two branches' agreement a claim rather than an assumption, so the
//! golden carries every table entry and the test compares all of them.
//!
//! # `FastMathCalc` is double-double arithmetic with its own conventions
//!
//! `split` puts the high 22 bits in `[0]` and the rest in `[1]`, and `resplit` renormalises. They
//! look like [`crate::dd`] and are not interchangeable with it: `split`'s threshold is the literal
//! `8e298`, and the large-value branch multiplies by `9.31322574615478515625E-10` rather than
//! scaling by a power of two. Reproducing the arithmetic means reproducing those constants.
//!
//! # The argument reduction has three special cases below zero
//!
//! ```java
//! if (x < -746d)      return 0.0;                                  // underflows
//! if (intVal < -709)  return exp(x + 40.19140625, ...) / 285040095144011776.0;
//! if (intVal == -709) return exp(x + 1.494140625, ...) / 4.455505956692756620;
//! ```
//!
//! Two of the three are recursive calls with a shifted argument and a magic divisor, which exist
//! so a subnormal result keeps its precision. A port that folded them into the general path would
//! agree everywhere except on the subnormal range, which is exactly where an exponential is
//! interesting to a likelihood.

// The literals below are the reference's, character for character, including the ones clippy
// would rather see written as a constant or shortened. `2.718281828459045` is the high half of a
// split representation of e and not `f64::consts::E`'s role, and the trailing digits of
// `9.31322574615478515625E-10` and `4.455505956692756620` are what the reference wrote: rounding
// them here would be a silent edit to an algorithm whose output is compared bit for bit.
#![allow(clippy::approx_constant, clippy::excessive_precision)]

use std::sync::OnceLock;

const HEX_40000000: f64 = 1073741824.0;
const EXP_INT_TABLE_MAX_INDEX: usize = 750;
const EXP_INT_TABLE_LEN: usize = EXP_INT_TABLE_MAX_INDEX * 2;
const EXP_FRAC_TABLE_LEN: usize = 1025;

/// `FastMathCalc.FACT`, the factorials `slowexp` divides by.
const FACT: [f64; 20] = [
    1.0,
    1.0,
    2.0,
    6.0,
    24.0,
    120.0,
    720.0,
    5040.0,
    40320.0,
    362880.0,
    3628800.0,
    39916800.0,
    479001600.0,
    6227020800.0,
    87178291200.0,
    1307674368000.0,
    20922789888000.0,
    355687428096000.0,
    6402373705728000.0,
    121645100408832000.0,
];

/// `FastMathCalc.split`: the high 22 bits in `[0]`, the remainder in `[1]`.
fn split(d: f64, out: &mut [f64; 2]) {
    if d < 8e298 && d > -8e298 {
        let a = d * HEX_40000000;
        out[0] = (d + a) - a;
        out[1] = d - out[0];
    } else {
        let a = d * 9.313_225_746_154_785_2e-10;
        out[0] = (d + a - d) * HEX_40000000;
        out[1] = d - out[0];
    }
}

/// `FastMathCalc.resplit`: renormalise a split in place.
fn resplit(a: &mut [f64; 2]) {
    let c = a[0] + a[1];
    let d = -(c - a[0] - a[1]);

    if c < 8e298 && c > -8e298 {
        let z = c * HEX_40000000;
        a[0] = (c + z) - z;
        a[1] = c - a[0] + d;
    } else {
        let z = c * 9.313_225_746_154_785_2e-10;
        a[0] = (c + z - c) * HEX_40000000;
        a[1] = c - a[0] + d;
    }
}

fn split_mult(a: &[f64; 2], b: &[f64; 2], ans: &mut [f64; 2]) {
    ans[0] = a[0] * b[0];
    ans[1] = a[0] * b[1] + a[1] * b[0] + a[1] * b[1];
    resplit(ans);
}

fn split_add(a: &[f64; 2], b: &[f64; 2], ans: &mut [f64; 2]) {
    ans[0] = a[0] + b[0];
    ans[1] = a[1] + b[1];
    resplit(ans);
}

/// `FastMathCalc.splitReciprocal`, including its two refinement passes.
///
/// The comment on the loop in the reference reads *"this may be overkill, probably once is
/// enough"*, and it runs twice anyway. Running it once produces different bits.
fn split_reciprocal(input: &mut [f64; 2], result: &mut [f64; 2]) {
    let b = 1.0 / 4194304.0;
    let a = 1.0 - b;

    if input[0] == 0.0 {
        input[0] = input[1];
        input[1] = 0.0;
    }

    result[0] = a / input[0];
    result[1] = (b * input[0] - a * input[1]) / (input[0] * input[0] + input[0] * input[1]);

    if result[1].is_nan() {
        result[1] = 0.0;
    }

    resplit(result);

    for _ in 0..2 {
        let mut err = 1.0
            - result[0] * input[0]
            - result[0] * input[1]
            - result[1] * input[0]
            - result[1] * input[1];
        err *= result[0] + result[1];
        result[1] += err;
    }
}

/// `FastMathCalc.quadMult`: `(a[0] + a[1]) * (b[0] + b[1])` in extended precision.
fn quad_mult(a: &[f64; 2], b: &[f64; 2], result: &mut [f64; 2]) {
    let mut xs = [0.0f64; 2];
    let mut ys = [0.0f64; 2];
    let mut zs = [0.0f64; 2];

    split(a[0], &mut xs);
    split(b[0], &mut ys);
    split_mult(&xs, &ys, &mut zs);
    result[0] = zs[0];
    result[1] = zs[1];

    let accumulate = |zs: &[f64; 2], result: &mut [f64; 2]| {
        let mut tmp = result[0] + zs[0];
        result[1] -= tmp - result[0] - zs[0];
        result[0] = tmp;
        tmp = result[0] + zs[1];
        result[1] -= tmp - result[0] - zs[1];
        result[0] = tmp;
    };

    split(b[1], &mut ys);
    split_mult(&xs, &ys, &mut zs);
    accumulate(&zs, result);

    split(a[1], &mut xs);
    split(b[0], &mut ys);
    split_mult(&xs, &ys, &mut zs);
    accumulate(&zs, result);

    // The reference's comment here says `a[1] * b[0]` a second time; the code is `a[1] * b[1]`.
    split(a[1], &mut xs);
    split(b[1], &mut ys);
    split_mult(&xs, &ys, &mut zs);
    accumulate(&zs, result);
}

/// `FastMathCalc.slowexp`: `exp(x)` for `x` in `[0, 1]`, by the Taylor series in split form.
fn slow_exp(x: f64, result: &mut [f64; 2]) {
    let mut xs = [0.0f64; 2];
    let mut ys = [0.0f64; 2];
    let mut facts = [0.0f64; 2];
    let mut as_ = [0.0f64; 2];
    split(x, &mut xs);
    ys[0] = 0.0;
    ys[1] = 0.0;

    for index in (0..FACT.len()).rev() {
        split_mult(&xs, &ys, &mut as_);
        ys[0] = as_[0];
        ys[1] = as_[1];

        split(FACT[index], &mut as_);
        let mut divisor = as_;
        split_reciprocal(&mut divisor, &mut facts);

        split_add(&ys, &facts, &mut as_);
        ys[0] = as_[0];
        ys[1] = as_[1];
    }

    result[0] = ys[0];
    result[1] = ys[1];
}

/// `FastMathCalc.expint`: `exp(p)` for an integer `p`, by binary exponentiation of `e` in split
/// form.
fn exp_int(mut p: i32, result: &mut [f64; 2]) {
    let mut xs = [2.718_281_828_459_045, 1.445_646_891_729_250_2e-16];
    let mut ys = [0.0f64; 2];
    let mut as_ = [0.0f64; 2];
    split(1.0, &mut ys);

    while p > 0 {
        if p & 1 != 0 {
            quad_mult(&ys, &xs, &mut as_);
            ys[0] = as_[0];
            ys[1] = as_[1];
        }
        let squared = xs;
        quad_mult(&squared, &squared, &mut as_);
        xs[0] = as_[0];
        xs[1] = as_[1];
        p >>= 1;
    }

    result[0] = ys[0];
    result[1] = ys[1];
    resplit(result);
}

/// The four tables, computed once, exactly as `RECOMPUTE_TABLES_AT_RUNTIME` computes them.
pub struct Tables {
    pub exp_int_a: Vec<f64>,
    pub exp_int_b: Vec<f64>,
    pub exp_frac_a: Vec<f64>,
    pub exp_frac_b: Vec<f64>,
}

static TABLES: OnceLock<Tables> = OnceLock::new();

pub fn tables() -> &'static Tables {
    TABLES.get_or_init(|| {
        let mut exp_int_a = vec![0.0; EXP_INT_TABLE_LEN];
        let mut exp_int_b = vec![0.0; EXP_INT_TABLE_LEN];
        let mut tmp = [0.0f64; 2];
        let mut recip = [0.0f64; 2];

        for i in 0..EXP_INT_TABLE_MAX_INDEX {
            exp_int(i as i32, &mut tmp);
            exp_int_a[i + EXP_INT_TABLE_MAX_INDEX] = tmp[0];
            exp_int_b[i + EXP_INT_TABLE_MAX_INDEX] = tmp[1];

            if i != 0 {
                let mut input = tmp;
                split_reciprocal(&mut input, &mut recip);
                exp_int_a[EXP_INT_TABLE_MAX_INDEX - i] = recip[0];
                exp_int_b[EXP_INT_TABLE_MAX_INDEX - i] = recip[1];
            }
        }

        let mut exp_frac_a = vec![0.0; EXP_FRAC_TABLE_LEN];
        let mut exp_frac_b = vec![0.0; EXP_FRAC_TABLE_LEN];
        let factor = 1.0 / (EXP_FRAC_TABLE_LEN - 1) as f64;
        for i in 0..EXP_FRAC_TABLE_LEN {
            slow_exp(i as f64 * factor, &mut tmp);
            exp_frac_a[i] = tmp[0];
            exp_frac_b[i] = tmp[1];
        }

        Tables {
            exp_int_a,
            exp_int_b,
            exp_frac_a,
            exp_frac_b,
        }
    })
}

/// `FastMath.exp(double)`.
pub fn exp(x: f64) -> f64 {
    exp_with_extra(x, 0.0, None)
}

/// `FastMath.exp(double, double, double[])`, the form `pow` calls with an extra term and a
/// high-precision output.
///
/// Kept public because `pow` and `expm1` need exactly this entry point, and because the extra
/// term changes the order of the final additions rather than being an optimisation.
pub fn exp_with_extra(x: f64, extra: f64, mut hi_prec: Option<&mut [f64; 2]>) -> f64 {
    let tables = tables();
    let mut int_val = x as i32;

    if x < 0.0 {
        // Not tested against `intVal`: the reference notes that converting a large negative
        // double may be affected by a JIT bug, so the comparison is on `x` itself.
        if x < -746.0 {
            if let Some(hi) = hi_prec.as_deref_mut() {
                hi[0] = 0.0;
                hi[1] = 0.0;
            }
            return 0.0;
        }
        if int_val < -709 {
            // The shift keeps a subnormal result precise; the divisor is exp(40.19140625).
            let mut inner = [0.0f64; 2];
            let result = exp_with_extra(
                x + 40.19140625,
                extra,
                hi_prec.as_deref_mut().map(|_| &mut inner),
            ) / 285040095144011776.0;
            if let Some(hi) = hi_prec.as_deref_mut() {
                hi[0] = inner[0] / 285040095144011776.0;
                hi[1] = inner[1] / 285040095144011776.0;
            }
            return result;
        }
        if int_val == -709 {
            let mut inner = [0.0f64; 2];
            let result = exp_with_extra(
                x + 1.494140625,
                extra,
                hi_prec.as_deref_mut().map(|_| &mut inner),
            ) / 4.455505956692756620;
            if let Some(hi) = hi_prec.as_deref_mut() {
                hi[0] = inner[0] / 4.455505956692756620;
                hi[1] = inner[1] / 4.455505956692756620;
            }
            return result;
        }
        int_val -= 1;
    } else if int_val > 709 {
        if let Some(hi) = hi_prec.as_deref_mut() {
            hi[0] = f64::INFINITY;
            hi[1] = 0.0;
        }
        return f64::INFINITY;
    }

    let index = (EXP_INT_TABLE_MAX_INDEX as i32 + int_val) as usize;
    let int_part_a = tables.exp_int_a[index];
    let int_part_b = tables.exp_int_b[index];

    let int_frac = ((x - int_val as f64) * 1024.0) as i32 as usize;
    let frac_part_a = tables.exp_frac_a[int_frac];
    let frac_part_b = tables.exp_frac_b[int_frac];

    // The subtraction is done last on purpose: doing it earlier loses precision.
    let epsilon = x - (int_val as f64 + int_frac as f64 / 1024.0);

    // The Remez polynomial for exp(epsilon) - 1 on [0, 2^-10].
    let mut z = 0.04168701738764507;
    z = z * epsilon + 0.1666666505023083;
    z = z * epsilon + 0.5000000000042687;
    z = z * epsilon + 1.0;
    z = z * epsilon + -3.940510424527919e-20;

    let temp_a = int_part_a * frac_part_a;
    let temp_b = int_part_a * frac_part_b + int_part_b * frac_part_a + int_part_b * frac_part_b;

    // "Order of operations is important. For accuracy add by increasing size."
    let temp_c = temp_b + temp_a;
    let result = if extra != 0.0 {
        temp_c * extra * z + temp_c * extra + temp_c * z + temp_b + temp_a
    } else {
        temp_c * z + temp_b + temp_a
    };

    if let Some(hi) = hi_prec {
        hi[0] = temp_a;
        hi[1] = temp_c * extra * z + temp_c * extra + temp_c * z + temp_b;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_tables_have_the_shape_the_reference_declares() {
        let tables = tables();
        assert_eq!(tables.exp_int_a.len(), 1500);
        assert_eq!(tables.exp_frac_a.len(), 1025);
        // exp(0) is 1, exactly, in both halves of the split.
        assert_eq!(tables.exp_int_a[EXP_INT_TABLE_MAX_INDEX], 1.0);
        assert_eq!(tables.exp_frac_a[0], 1.0);
    }

    #[test]
    fn the_ordinary_range_is_close_to_the_true_exponential() {
        // Exactness is the golden's job; this only catches a port that is wrong by an order of
        // magnitude, which a table indexing slip would be.
        for x in [0.0, 1.0, -1.0, 10.0, -10.0, 0.5, 700.0] {
            let ours = exp(x);
            assert!(
                (ours - x.exp()).abs() <= 1e-9 * x.exp().abs().max(1e-300),
                "exp({x}) = {ours}"
            );
        }
    }

    #[test]
    fn the_three_special_cases_below_zero() {
        assert_eq!(exp(-747.0), 0.0);
        assert!(exp(-710.0) > 0.0 && exp(-710.0).is_subnormal() || exp(-710.0) > 0.0);
        assert_eq!(exp(710.0), f64::INFINITY);
    }
}
