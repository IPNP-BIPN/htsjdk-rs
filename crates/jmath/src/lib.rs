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

pub mod dd;
mod log;

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

    // `Math.exp` is WITHDRAWN, not missing. It was implemented as an operation-by-operation
    // transcription of HotSpot's x86 intrinsic, and that source file is GPL2 *only*, with no
    // Classpath Exception, so the transcription could not be published under this crate's MIT
    // licence. Removed in decision 0014. `Math.exp` now has the same status as `Math.pow`:
    // unported, with the reason recorded rather than the gap left unexplained.
}

/// `java.lang.StrictMath`. fdlibm, portable by specification.
pub mod strict_math {
    /// See [`crate::math::sqrt`]: exact in every implementation.
    #[inline]
    pub fn sqrt(x: f64) -> f64 {
        x.sqrt()
    }
}

/// `org.apache.commons.math3.util.FastMath`. Pure Java, portable, and a third set of answers.
pub mod fast_math {
    /// See [`crate::math::sqrt`]: exact in every implementation.
    #[inline]
    pub fn sqrt(x: f64) -> f64 {
        x.sqrt()
    }
}
