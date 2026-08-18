//! Java-semantics floating point.
//!
//! Shared by the whole program, not specific to htsjdk: Picard's
//! `DuplicationMetrics.estimateLibrarySize` and GATK's `RecalDatum` both depend on it. It lives
//! here because `htsjdk-rs` is the root of the dependency chain.
//!
//! # There is no single "Java math"
//!
//! Three implementations coexist in the reference and disagree pairwise by up to 1 ULP:
//!
//! | | portable? | used by |
//! |---|---|---|
//! | [`math`] (`java.lang.Math`) | **no**, HotSpot intrinsic | `RecalDatum.log10`, `MathUtils.pow` |
//! | [`strict_math`] (`java.lang.StrictMath`) | yes, fdlibm | rarely in GATK |
//! | [`fast_math`] (`commons-math3 FastMath`) | yes, pure Java | `MathUtils` |
//!
//! A ported call site must name the same one the reference used. That is why there is no
//! top-level `jmath::exp`: a single blessed entry point is precisely the API shape that lets a
//! call site bind silently to the wrong implementation.
//!
//! See `docs/decisions/0005-java-math-has-three-implementations.md` for the measurement.

pub mod beta;
pub mod binomial;
pub mod brent;
pub mod combinatorics;
pub mod continued_fraction;
pub mod dd;
pub mod fast_math_exp;
pub mod fast_math_log;
pub mod fast_math_tables;
pub mod gamma;
mod log;
pub mod normal;
pub mod percentile;
pub mod saddle_point;
pub mod strict_exp;
pub mod strict_pow;

/// `java.lang.Math`. Platform-specific HotSpot intrinsics; the target for most GATK call sites.
pub mod math {
    /// IEEE-754 mandates a correctly-rounded square root, so every implementation agrees
    /// exactly and Rust's is already bit-identical. Verified over the whole corpus.
    #[inline]
    pub fn sqrt(x: f64) -> f64 {
        x.sqrt()
    }

    /// `Math.log`. Measured to be correctly rounded, so the port rounds the true result once
    /// rather than reproducing HotSpot's intrinsic.
    #[inline]
    pub fn log(x: f64) -> f64 {
        crate::log::log(x)
    }

    /// `Math.log10`. Correctly rounded, as above.
    #[inline]
    pub fn log10(x: f64) -> f64 {
        crate::log::log10(x)
    }

    /// `Math.round(double)`, which is **not** `floor(x + 0.5)` however much its own javadoc says
    /// so.
    ///
    /// The javadoc has said `floor(x + 0.5)` since Java 1.0 and the implementation stopped doing
    /// that in Java 7. The two disagree on exactly one class of input, and the conformance golden
    /// in `htsjdk-vcf` carries the witness:
    ///
    /// ```text
    /// round  <0.49999999999999994>  0
    /// ```
    ///
    /// `0.49999999999999994` is the double immediately below a half. Adding `0.5` to it rounds
    /// **up** to exactly `1.0`, so the arithmetic version rounds twice and answers 1 where the
    /// correct half-up answer is 0. JDK-8010430 replaced the arithmetic with bit manipulation that
    /// cannot round twice, and this is that code: take the significand, shift it into place by the
    /// unbiased exponent, add one and shift once more, so the half-up decision is made on the exact
    /// bits and never on a sum that has already been rounded.
    ///
    /// The half-up rule itself is unchanged and is what makes the two look alike on ordinary
    /// input: -1.5 rounds to -1, not to -2.
    ///
    /// It lives here rather than beside its callers because it is a `java.lang.Math` function and
    /// a second copy of it would be a second definition of what rounding means.
    pub fn round(value: f64) -> i64 {
        /// `DoubleConsts.SIGNIFICAND_WIDTH`.
        const SIGNIFICAND_WIDTH: i64 = 53;
        /// `DoubleConsts.EXP_BIAS`.
        const EXP_BIAS: i64 = 1023;
        const EXP_BIT_MASK: i64 = 0x7FF0_0000_0000_0000u64 as i64;
        const SIGNIF_BIT_MASK: i64 = 0x000F_FFFF_FFFF_FFFF;

        let long_bits = value.to_bits() as i64;
        let biased_exp = (long_bits & EXP_BIT_MASK) >> (SIGNIFICAND_WIDTH - 1);
        let shift = (SIGNIFICAND_WIDTH - 2 + EXP_BIAS) - biased_exp;
        if (shift & -64) == 0 {
            // shift is in [0, 64): the value is representable and the significand can be shifted.
            let mut r = (long_bits & SIGNIF_BIT_MASK) | (SIGNIF_BIT_MASK + 1);
            if long_bits < 0 {
                r = -r;
            }
            ((r >> shift) + 1) >> 1
        } else {
            // Too large, too small, NaN or infinite: `(long) someDouble`, which saturates at both
            // ends and answers zero for NaN. Rust's `as` has had exactly those semantics since
            // 1.45.
            value as i64
        }
    }

    // `Math.exp` is WITHDRAWN, not missing. It was implemented as an operation-by-operation
    // transcription of HotSpot's x86 intrinsic, and that source file is GPL2 *only*, with no
    // Classpath Exception, so the transcription could not be published under this crate's MIT
    // licence. Removed in decision 0014. `Math.exp` now has the same status as `Math.pow`:
    // unported, with the reason recorded rather than the gap left unexplained.
    //
    // `strict_math::exp` is NOT a substitute for it, and decision 0025 measured why: FDLIBM agrees
    // with `Math.exp` on 98.6443% of the corpus where the system libm agrees on 99.9711%. The
    // permissive implementation is the *worse* stand-in for the intrinsic, and it would have been
    // easy to assume the opposite. What both agree on is the size of the gap: no divergence
    // exceeds 1 ulp.
}

/// `java.lang.StrictMath`. fdlibm, portable by specification.
pub mod strict_math {
    /// See [`crate::math::sqrt`]: exact in every implementation.
    #[inline]
    pub fn sqrt(x: f64) -> f64 {
        x.sqrt()
    }

    /// `StrictMath.exp`, exact over every point of the conformance corpus.
    ///
    /// Portable where `Math.exp` is not, and for a licence reason rather than a difficulty one:
    /// `StrictMath` is *specified* to be FDLIBM, and FDLIBM's notice grants the right to
    /// translate it. See [`crate::strict_exp`] and decision 0025.
    #[inline]
    pub fn exp(x: f64) -> f64 {
        crate::strict_exp::exp(x)
    }

    /// `StrictMath.pow`, exact over every point of the conformance corpus.
    ///
    /// Same standing as [`exp`], and the same limitation: this is not `Math.pow`, which is a
    /// HotSpot intrinsic that decision 0007 deferred. What it makes possible is the measurement
    /// that decision 0007 never had — a distance to the intrinsic rather than a rate of
    /// agreement. See [`crate::strict_pow`].
    #[inline]
    pub fn pow(x: f64, y: f64) -> f64 {
        crate::strict_pow::pow(x, y)
    }
}

/// `org.apache.commons.math3.util.FastMath`. Pure Java, portable, and a third set of answers.
pub mod fast_math {
    /// See [`crate::math::sqrt`]: exact in every implementation.
    #[inline]
    pub fn sqrt(x: f64) -> f64 {
        x.sqrt()
    }

    /// The four exponential tables, computed as `FastMathCalc` computes them.
    ///
    /// Public so the conformance suite can compare every entry against the literals the reference
    /// ships, which is what makes "computed rather than transcribed" a measured claim.
    pub fn exp_tables() -> &'static crate::fast_math_exp::Tables {
        crate::fast_math_exp::tables()
    }

    /// The tables the reference's *other* branch would compute, which differ. See decision 0024.
    pub fn recomputed_exp_tables() -> &'static crate::fast_math_exp::Tables {
        crate::fast_math_exp::recomputed_tables()
    }

    /// `FastMath.exp`, table-driven and pure Java, which is **not** `java.lang.Math.exp`.
    #[inline]
    pub fn exp(x: f64) -> f64 {
        crate::fast_math_exp::exp(x)
    }

    /// `FastMath.log`, table-driven and **not** correctly rounded, unlike `java.lang.Math.log`.
    #[inline]
    pub fn log(x: f64) -> f64 {
        crate::fast_math_log::log(x)
    }

    /// `FastMath.log1p`, which is `log(1 + x)` computed two ways around |x| = 1e-6.
    #[inline]
    pub fn log1p(x: f64) -> f64 {
        crate::fast_math_log::log1p(x)
    }

    /// `FastMath.round(double)`, which is literally `(long) floor(x + 0.5)`.
    ///
    /// That is the definition [`crate::math::round`] **stopped** using in Java 7, and the two are
    /// still one apart on `0.49999999999999994`. Both live in this crate because both are reached
    /// from ported code: `MathUtils.median` ends in this one, `htsjdk` ends in the other.
    ///
    /// The cast is Java's `(long)` narrowing, which clamps rather than wrapping and answers 0 for
    /// `NaN`.
    #[inline]
    pub fn round(x: f64) -> i64 {
        let shifted = (x + 0.5).floor();
        if shifted.is_nan() {
            return 0;
        }
        if shifted >= i64::MAX as f64 {
            return i64::MAX;
        }
        if shifted <= i64::MIN as f64 {
            return i64::MIN;
        }
        shifted as i64
    }
}
