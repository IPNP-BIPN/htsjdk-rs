//! `Gamma` and `Erf`, ported from `org.apache.commons.math3.special` (commons-math3 3.5), which
//! is Apache 2.0. See `docs/decisions/0023-commons-math3-is-portable-where-the-jdk-is-not.md`.
//!
//! `Erf` is the whole of `NormalDistribution.cumulativeProbability`, which is the whole of the
//! normal approximation in GATK's `MannWhitneyU`, which is what every rank-sum annotation reports.
//! So this file is four layers below `MQRankSum` and decides its last digit.
//!
//! # The error function is a gamma integral, and its tolerance is an argument
//!
//! ```java
//! Gamma.regularizedGammaP(0.5, x * x, 1.0e-15, 10000)
//! ```
//!
//! `erf` is not a polynomial approximation here: it is an incomplete gamma function evaluated by
//! series or by continued fraction, stopped by a **relative** tolerance of 1e-15 or by 10,000
//! iterations. The iteration count is part of the result: a port that iterated to convergence
//! would agree almost everywhere and differ wherever the reference stopped first.
//!
//! # Which of the two representations runs is decided by `x >= a + 1`
//!
//! `regularizedGammaP` delegates to `regularizedGammaQ` when `x >= a + 1`, and `Q` delegates back
//! to `P` when `x < a + 1`, so the pair is one function with a switch in the middle. Around
//! `x = a + 1` the two sides are computed by different algorithms, and their last bits differ.
//!
//! # `logGamma` has four branches and only two of them are reachable from here
//!
//! `erf` always asks for `a = 0.5`, which takes the `x <= 2.5` branch into `logGamma1p`, and that
//! is a rational approximation with 36 transcribed constants. The Lanczos branch, for `x > 8`, is
//! ported too because `MannWhitneyU` reaches it through the binomial coefficient.
//!
//! # `abs(x) > 40` is a shortcut, not a limit
//!
//! ```java
//! if (FastMath.abs(x) > 40) { return x > 0 ? 1 : -1; }
//! ```
//!
//! Beyond 40 the reference returns exactly ±1 without computing anything, so `erf(40.0)` and
//! `erf(40.000001)` are computed two different ways and the discontinuity is at a value the
//! author chose.

#![allow(clippy::excessive_precision)]

use crate::fast_math;

const INV_GAMMA1P_M1_A0: f64 = 0.611609510448141581788E-08;
const INV_GAMMA1P_M1_A1: f64 = 0.624730830116465516210E-08;
const INV_GAMMA1P_M1_B1: f64 = 0.203610414066806987300E+00;
const INV_GAMMA1P_M1_B2: f64 = 0.266205348428949217746E-01;
const INV_GAMMA1P_M1_B3: f64 = 0.493944979382446875238E-03;
const INV_GAMMA1P_M1_B4: f64 = -0.851419432440314906588E-05;
const INV_GAMMA1P_M1_B5: f64 = -0.643045481779353022248E-05;
const INV_GAMMA1P_M1_B6: f64 = 0.992641840672773722196E-06;
const INV_GAMMA1P_M1_B7: f64 = -0.607761895722825260739E-07;
const INV_GAMMA1P_M1_B8: f64 = 0.195755836614639731882E-09;
const INV_GAMMA1P_M1_P0: f64 = 0.6116095104481415817861E-08;
const INV_GAMMA1P_M1_P1: f64 = 0.6871674113067198736152E-08;
const INV_GAMMA1P_M1_P2: f64 = 0.6820161668496170657918E-09;
const INV_GAMMA1P_M1_P3: f64 = 0.4686843322948848031080E-10;
const INV_GAMMA1P_M1_P4: f64 = 0.1572833027710446286995E-11;
const INV_GAMMA1P_M1_P5: f64 = -0.1249441572276366213222E-12;
const INV_GAMMA1P_M1_P6: f64 = 0.4343529937408594255178E-14;
const INV_GAMMA1P_M1_Q1: f64 = 0.3056961078365221025009E+00;
const INV_GAMMA1P_M1_Q2: f64 = 0.5464213086042296536016E-01;
const INV_GAMMA1P_M1_Q3: f64 = 0.4956830093825887312020E-02;
const INV_GAMMA1P_M1_Q4: f64 = 0.2692369466186361192876E-03;
const INV_GAMMA1P_M1_C: f64 = -0.422784335098467139393487909917598E+00;
const INV_GAMMA1P_M1_C0: f64 = 0.577215664901532860606512090082402E+00;
const INV_GAMMA1P_M1_C1: f64 = -0.655878071520253881077019515145390E+00;
const INV_GAMMA1P_M1_C2: f64 = -0.420026350340952355290039348754298E-01;
const INV_GAMMA1P_M1_C3: f64 = 0.166538611382291489501700795102105E+00;
const INV_GAMMA1P_M1_C4: f64 = -0.421977345555443367482083012891874E-01;
const INV_GAMMA1P_M1_C5: f64 = -0.962197152787697356211492167234820E-02;
const INV_GAMMA1P_M1_C6: f64 = 0.721894324666309954239501034044657E-02;
const INV_GAMMA1P_M1_C7: f64 = -0.116516759185906511211397108401839E-02;
const INV_GAMMA1P_M1_C8: f64 = -0.215241674114950972815729963053648E-03;
const INV_GAMMA1P_M1_C9: f64 = 0.128050282388116186153198626328164E-03;
const INV_GAMMA1P_M1_C10: f64 = -0.201348547807882386556893914210218E-04;
const INV_GAMMA1P_M1_C11: f64 = -0.125049348214267065734535947383309E-05;
const INV_GAMMA1P_M1_C12: f64 = 0.113302723198169588237412962033074E-05;
const INV_GAMMA1P_M1_C13: f64 = -0.205633841697760710345015413002057E-06;
const LANCZOS: [f64; 15] = [
    0.99999999999999709182,
    57.156235665862923517,
    -59.597960355475491248,
    14.136097974741747174,
    -0.49191381609762019978,
    0.33994649984811888699e-4,
    0.46523628927048575665e-4,
    -0.98374475304879564677e-4,
    0.15808870322491248884e-3,
    -0.21026444172410488319e-3,
    0.21743961811521264320e-3,
    -0.16431810653676389022e-3,
    0.84418223983852743293e-4,
    -0.26190838401581408670e-4,
    0.36899182659531622704e-5,
];

/// `Gamma.LANCZOS_G`.
const LANCZOS_G: f64 = 607.0 / 128.0;

/// `Gamma.HALF_LOG_2_PI`, which the reference computes with `FastMath.log`, not `Math.log`.
fn half_log_2_pi() -> f64 {
    0.5 * fast_math::log(2.0 * std::f64::consts::PI)
}

/// What the reference throws rather than answering.
#[derive(Debug, Clone, PartialEq)]
pub enum GammaError {
    /// `NumberIsTooSmallException`: an argument below -0.5 in `invGamma1pm1`/`logGamma1p`.
    TooSmall { value: f64, bound: f64 },
    /// `NumberIsTooLargeException`: an argument above 1.5 in the same pair.
    TooLarge { value: f64, bound: f64 },
    /// `MaxCountExceededException`: the series or the continued fraction ran out of iterations.
    MaxCountExceeded { max: i32 },
    /// `ConvergenceException` from the continued fraction, which distinguishes an infinite value
    /// from a NaN one.
    ContinuedFractionDiverged { infinite: bool, x: f64 },
    /// `StackOverflowError`: `digamma` and `trigamma` in 3.5 have no NaN guard, so `NaN` and
    /// `-Infinity` recurse forever. Not an `Exception`, so a caller that catches one does not catch
    /// this.
    NonTerminating,
}

/// `Gamma.lanczos`, summed from the far end so the small terms are added first.
pub fn lanczos(x: f64) -> f64 {
    let mut sum = 0.0;
    for i in (1..LANCZOS.len()).rev() {
        sum += LANCZOS[i] / (x + i as f64);
    }
    sum + LANCZOS[0]
}

/// `Gamma.invGamma1pm1`: `1 / Gamma(1 + x) - 1`, for `x` in `[-0.5, 1.5]`.
pub fn inv_gamma1pm1(x: f64) -> Result<f64, GammaError> {
    if x < -0.5 {
        return Err(GammaError::TooSmall {
            value: x,
            bound: -0.5,
        });
    }
    if x > 1.5 {
        return Err(GammaError::TooLarge {
            value: x,
            bound: 1.5,
        });
    }

    let t = if x <= 0.5 { x } else { (x - 0.5) - 0.5 };
    let ret = if t < 0.0 {
        let a = INV_GAMMA1P_M1_A0 + t * INV_GAMMA1P_M1_A1;
        let mut b = INV_GAMMA1P_M1_B8;
        b = INV_GAMMA1P_M1_B7 + t * b;
        b = INV_GAMMA1P_M1_B6 + t * b;
        b = INV_GAMMA1P_M1_B5 + t * b;
        b = INV_GAMMA1P_M1_B4 + t * b;
        b = INV_GAMMA1P_M1_B3 + t * b;
        b = INV_GAMMA1P_M1_B2 + t * b;
        b = INV_GAMMA1P_M1_B1 + t * b;
        b = 1.0 + t * b;

        let mut c = INV_GAMMA1P_M1_C13 + t * (a / b);
        c = INV_GAMMA1P_M1_C12 + t * c;
        c = INV_GAMMA1P_M1_C11 + t * c;
        c = INV_GAMMA1P_M1_C10 + t * c;
        c = INV_GAMMA1P_M1_C9 + t * c;
        c = INV_GAMMA1P_M1_C8 + t * c;
        c = INV_GAMMA1P_M1_C7 + t * c;
        c = INV_GAMMA1P_M1_C6 + t * c;
        c = INV_GAMMA1P_M1_C5 + t * c;
        c = INV_GAMMA1P_M1_C4 + t * c;
        c = INV_GAMMA1P_M1_C3 + t * c;
        c = INV_GAMMA1P_M1_C2 + t * c;
        c = INV_GAMMA1P_M1_C1 + t * c;
        c = INV_GAMMA1P_M1_C + t * c;
        if x > 0.5 {
            t * c / x
        } else {
            // `(c + 0.5) + 0.5` rather than `c + 1.0`: the two are different doubles.
            x * ((c + 0.5) + 0.5)
        }
    } else {
        let mut p = INV_GAMMA1P_M1_P6;
        p = INV_GAMMA1P_M1_P5 + t * p;
        p = INV_GAMMA1P_M1_P4 + t * p;
        p = INV_GAMMA1P_M1_P3 + t * p;
        p = INV_GAMMA1P_M1_P2 + t * p;
        p = INV_GAMMA1P_M1_P1 + t * p;
        p = INV_GAMMA1P_M1_P0 + t * p;

        let mut q = INV_GAMMA1P_M1_Q4;
        q = INV_GAMMA1P_M1_Q3 + t * q;
        q = INV_GAMMA1P_M1_Q2 + t * q;
        q = INV_GAMMA1P_M1_Q1 + t * q;
        q = 1.0 + t * q;

        let mut c = INV_GAMMA1P_M1_C13 + (p / q) * t;
        c = INV_GAMMA1P_M1_C12 + t * c;
        c = INV_GAMMA1P_M1_C11 + t * c;
        c = INV_GAMMA1P_M1_C10 + t * c;
        c = INV_GAMMA1P_M1_C9 + t * c;
        c = INV_GAMMA1P_M1_C8 + t * c;
        c = INV_GAMMA1P_M1_C7 + t * c;
        c = INV_GAMMA1P_M1_C6 + t * c;
        c = INV_GAMMA1P_M1_C5 + t * c;
        c = INV_GAMMA1P_M1_C4 + t * c;
        c = INV_GAMMA1P_M1_C3 + t * c;
        c = INV_GAMMA1P_M1_C2 + t * c;
        c = INV_GAMMA1P_M1_C1 + t * c;
        c = INV_GAMMA1P_M1_C0 + t * c;

        if x > 0.5 {
            (t / x) * ((c - 0.5) - 0.5)
        } else {
            x * c
        }
    };
    Ok(ret)
}

/// `Gamma.logGamma1p`.
pub fn log_gamma1p(x: f64) -> Result<f64, GammaError> {
    if x < -0.5 {
        return Err(GammaError::TooSmall {
            value: x,
            bound: -0.5,
        });
    }
    if x > 1.5 {
        return Err(GammaError::TooLarge {
            value: x,
            bound: 1.5,
        });
    }
    Ok(-fast_math::log1p(inv_gamma1pm1(x)?))
}

/// `Gamma.gamma(x)` for `abs(x) <= 20`, which is the whole of what a beta function reaches.
///
/// The reference has a third arm for `abs(x) > 20`, built on `FastMath.pow` and `FastMath.exp`
/// around the Lanczos series. Nothing measured here reaches it — `Beta.logBeta` calls this only
/// with `a < 1` and `b < 10`, so the largest argument is `a + b < 11` — and it refuses rather than
/// guessing, on the same rule as the treeified bucket in the hash-order port.
pub fn gamma(x: f64) -> Result<f64, GammaError> {
    // A non-positive integer is a pole.
    if x == x.round() && x <= 0.0 {
        return Ok(f64::NAN);
    }
    let absolute = x.abs();
    if absolute > 20.0 {
        return Err(GammaError::TooLarge {
            value: x,
            bound: 20.0,
        });
    }
    if x >= 1.0 {
        // Gamma(x) = (x - 1) * ... * (x - n) * Gamma(x - n), reduced until the argument is in the
        // range `invGamma1pm1` accepts.
        let mut prod = 1.0;
        let mut t = x;
        while t > 2.5 {
            t -= 1.0;
            prod *= t;
        }
        return Ok(prod / (1.0 + inv_gamma1pm1(t - 1.0)?));
    }
    // Gamma(x) = Gamma(x + n + 1) / [x * (x + 1) * ... * (x + n)], reduced upwards instead.
    let mut prod = x;
    let mut t = x;
    while t < -0.5 {
        t += 1.0;
        prod *= t;
    }
    Ok(1.0 / (prod * (1.0 + inv_gamma1pm1(t)?)))
}

/// `Gamma.logGamma`, all four branches.
pub fn log_gamma(x: f64) -> f64 {
    if x.is_nan() || x <= 0.0 {
        return f64::NAN;
    }
    if x < 0.5 {
        return log_gamma1p(x).expect("x < 0.5 is in range") - fast_math::log(x);
    }
    if x <= 2.5 {
        // `(x - 0.5) - 0.5`, again not `x - 1.0`.
        return log_gamma1p((x - 0.5) - 0.5).expect("x <= 2.5 is in range");
    }
    if x <= 8.0 {
        let n = (x - 1.5).floor() as i32;
        let mut prod = 1.0;
        for i in 1..=n {
            prod *= x - i as f64;
        }
        return log_gamma1p(x - (n as f64 + 1.0)).expect("the reduction lands in range")
            + fast_math::log(prod);
    }
    let sum = lanczos(x);
    let tmp = x + LANCZOS_G + 0.5;
    ((x + 0.5) * fast_math::log(tmp)) - tmp + half_log_2_pi() + fast_math::log(sum / x)
}

/// `ContinuedFraction.evaluate`, the modified Lentz algorithm the gamma tail uses.
///
/// `Precision.equals(value, 0.0, small)` is an absolute comparison against 1e-50, not an exact
/// one, so a term that is merely tiny is replaced by 1e-50 as if it were zero.
fn continued_fraction(
    a: impl Fn(i32, f64) -> f64,
    b: impl Fn(i32, f64) -> f64,
    x: f64,
    epsilon: f64,
    max_iterations: i32,
) -> Result<f64, GammaError> {
    let small = 1e-50;
    let mut h_prev = a(0, x);
    if (h_prev - 0.0).abs() <= small {
        h_prev = small;
    }

    let mut n = 1;
    let mut d_prev = 0.0;
    let mut c_prev = h_prev;
    let mut h_n = h_prev;

    while n < max_iterations {
        let an = a(n, x);
        let bn = b(n, x);

        let mut d_n = an + bn * d_prev;
        if (d_n - 0.0).abs() <= small {
            d_n = small;
        }
        let mut c_n = an + bn / c_prev;
        if (c_n - 0.0).abs() <= small {
            c_n = small;
        }

        d_n = 1.0 / d_n;
        let delta_n = c_n * d_n;
        h_n = h_prev * delta_n;

        if h_n.is_infinite() {
            return Err(GammaError::ContinuedFractionDiverged { infinite: true, x });
        }
        if h_n.is_nan() {
            return Err(GammaError::ContinuedFractionDiverged { infinite: false, x });
        }

        if (delta_n - 1.0).abs() < epsilon {
            break;
        }

        d_prev = d_n;
        c_prev = c_n;
        h_prev = h_n;
        n += 1;
    }

    if n >= max_iterations {
        return Err(GammaError::MaxCountExceeded {
            max: max_iterations,
        });
    }
    Ok(h_n)
}

/// `Gamma.regularizedGammaP(a, x, epsilon, maxIterations)`.
pub fn regularized_gamma_p(
    a: f64,
    x: f64,
    epsilon: f64,
    max_iterations: i32,
) -> Result<f64, GammaError> {
    if a.is_nan() || x.is_nan() || a <= 0.0 || x < 0.0 {
        return Ok(f64::NAN);
    }
    if x == 0.0 {
        return Ok(0.0);
    }
    if x >= a + 1.0 {
        // The pair is one function with a switch in the middle: either side of `a + 1` is
        // computed by a different algorithm.
        return Ok(1.0 - regularized_gamma_q(a, x, epsilon, max_iterations)?);
    }
    let mut n = 0.0;
    let mut an = 1.0 / a;
    let mut sum = an;
    while (an / sum).abs() > epsilon && n < max_iterations as f64 && sum < f64::INFINITY {
        n += 1.0;
        an *= x / (a + n);
        sum += an;
    }
    if n >= max_iterations as f64 {
        return Err(GammaError::MaxCountExceeded {
            max: max_iterations,
        });
    }
    if sum.is_infinite() {
        return Ok(1.0);
    }
    Ok(fast_math::exp(-x + (a * fast_math::log(x)) - log_gamma(a)) * sum)
}

/// `Gamma.regularizedGammaQ(a, x, epsilon, maxIterations)`.
pub fn regularized_gamma_q(
    a: f64,
    x: f64,
    epsilon: f64,
    max_iterations: i32,
) -> Result<f64, GammaError> {
    if a.is_nan() || x.is_nan() || a <= 0.0 || x < 0.0 {
        return Ok(f64::NAN);
    }
    if x == 0.0 {
        return Ok(1.0);
    }
    if x < a + 1.0 {
        return Ok(1.0 - regularized_gamma_p(a, x, epsilon, max_iterations)?);
    }
    let value = continued_fraction(
        |n, x| ((2.0 * n as f64) + 1.0) - a + x,
        |n, _| n as f64 * (a - n as f64),
        x,
        epsilon,
        max_iterations,
    )?;
    let ret = 1.0 / value;
    Ok(fast_math::exp(-x + (a * fast_math::log(x)) - log_gamma(a)) * ret)
}

/// `Erf.erf(x)`, which is `regularizedGammaP(0.5, x*x)` with the reference's own tolerances.
pub fn erf(x: f64) -> f64 {
    if x.abs() > 40.0 {
        return if x > 0.0 { 1.0 } else { -1.0 };
    }
    let ret = regularized_gamma_p(0.5, x * x, 1.0e-15, 10000).unwrap_or(f64::NAN);
    if x < 0.0 {
        -ret
    } else {
        ret
    }
}

/// `Erf.erfc(x)`.
pub fn erfc(x: f64) -> f64 {
    if x.abs() > 40.0 {
        return if x > 0.0 { 0.0 } else { 2.0 };
    }
    let ret = regularized_gamma_q(0.5, x * x, 1.0e-15, 10000).unwrap_or(f64::NAN);
    if x < 0.0 {
        2.0 - ret
    } else {
        ret
    }
}

/// `Erf.erfInv(x)`, the inverse error function.
///
/// A rational approximation in three ranges, and unlike `erf` it touches no gamma function at
/// all: the two are inverses in mathematics and unrelated in code, so a round trip through them
/// does not return the input bit for bit.
pub fn erf_inv(x: f64) -> f64 {
    // Beyond the branches below the reference answers the infinities, and `|x| > 1` is NaN.
    let w = -fast_math::log((1.0 - x) * (1.0 + x));
    let mut p;

    if w < 6.25 {
        let w = w - 3.125;
        p = -3.6444120640178196996e-21;
        p = -1.685059138182016589e-19 + p * w;
        p = 1.2858480715256400167e-18 + p * w;
        p = 1.115787767802518096e-17 + p * w;
        p = -1.333171662854620906e-16 + p * w;
        p = 2.0972767875968561637e-17 + p * w;
        p = 6.6376381343583238325e-15 + p * w;
        p = -4.0545662729752068639e-14 + p * w;
        p = -8.1519341976054721522e-14 + p * w;
        p = 2.6335093153082322977e-12 + p * w;
        p = -1.2975133253453532498e-11 + p * w;
        p = -5.4154120542946279317e-11 + p * w;
        p = 1.051212273321532285e-09 + p * w;
        p = -4.1126339803469836976e-09 + p * w;
        p = -2.9070369957882005086e-08 + p * w;
        p = 4.2347877827932403518e-07 + p * w;
        p = -1.3654692000834678645e-06 + p * w;
        p = -1.3882523362786468719e-05 + p * w;
        p = 0.0001867342080340571352 + p * w;
        p = -0.00074070253416626697512 + p * w;
        p = -0.0060336708714301490533 + p * w;
        p = 0.24015818242558961693 + p * w;
        p = 1.6536545626831027356 + p * w;
    } else if w < 16.0 {
        let w = fast_math::sqrt(w) - 3.25;
        p = 2.2137376921775787049e-09;
        p = 9.0756561938885390979e-08 + p * w;
        p = -2.7517406297064545428e-07 + p * w;
        p = 1.8239629214389227755e-08 + p * w;
        p = 1.5027403968909827627e-06 + p * w;
        p = -4.013867526981545969e-06 + p * w;
        p = 2.9234449089955446044e-06 + p * w;
        p = 1.2475304481671778723e-05 + p * w;
        p = -4.7318229009055733981e-05 + p * w;
        p = 6.8284851459573175448e-05 + p * w;
        p = 2.4031110387097893999e-05 + p * w;
        p = -0.0003550375203628474796 + p * w;
        p = 0.00095328937973738049703 + p * w;
        p = -0.0016882755560235047313 + p * w;
        p = 0.0024914420961078508066 + p * w;
        p = -0.0037512085075692412107 + p * w;
        p = 0.005370914553590063617 + p * w;
        p = 1.0052589676941592334 + p * w;
        p = 3.0838856104922207635 + p * w;
    } else if !w.is_infinite() {
        let w = fast_math::sqrt(w) - 5.0;
        p = -2.7109920616438573243e-11;
        p = -2.5556418169965252055e-10 + p * w;
        p = 1.5076572693500548083e-09 + p * w;
        p = -3.7894654401267369937e-09 + p * w;
        p = 7.6157012080783393804e-09 + p * w;
        p = -1.4960026627149240478e-08 + p * w;
        p = 2.9147953450901080826e-08 + p * w;
        p = -6.7711997758452339498e-08 + p * w;
        p = 2.2900482228026654717e-07 + p * w;
        p = -9.9298272942317002539e-07 + p * w;
        p = 4.5260625972231537039e-06 + p * w;
        p = -1.9681778105531670567e-05 + p * w;
        p = 7.5995277030017761139e-05 + p * w;
        p = -0.00021503011930044477347 + p * w;
        p = -0.00013871931833623122026 + p * w;
        p = 1.0103004648645343977 + p * w;
        p = 4.8499064014085844221 + p * w;
    } else {
        // `w` is infinite, which means `x` was exactly ±1.
        p = f64::INFINITY;
    }

    p * x
}

/// `Gamma.GAMMA`, the Euler-Mascheroni constant as commons-math3 transcribes it.
const GAMMA: f64 = 0.577_215_664_901_532_9;
/// `Gamma.S_LIMIT`, below which `digamma` uses "method 5 from Bernardo AS103".
const S_LIMIT: f64 = 1e-5;
/// `Gamma.C_LIMIT`, above which the asymptotic expansion is used.
const C_LIMIT: f64 = 49.0;

/// `Gamma.digamma`: the logarithmic derivative of the gamma function.
///
/// Three branches, and the middle one is a **recursion**: below 49 it calls itself at `x + 1` and
/// subtracts `1 / x`, so a digamma of 0.1 unwinds forty-nine frames before any arithmetic that is
/// not a subtraction happens, and the result is a sum of forty-nine reciprocals in a fixed order.
/// That order is the result: rearranging the sum would change the last bits.
///
/// # There is no `NaN` guard in 3.5, so two inputs do not terminate
///
/// A later commons-math3 returns `x` unchanged for a `NaN` or an infinity. **3.5 does not have that
/// line.** `NaN` fails both branch tests, so the middle branch recurses on `NaN + 1`, which is
/// `NaN`, forever; `-Infinity` recurses on `-Infinity + 1`, which is `-Infinity`, forever. Both
/// raise a `StackOverflowError`, which is an `Error` and not an `Exception`, so a caller catching
/// `Exception` does not catch it. The golden records both, and this port refuses them rather than
/// inventing the guard the version in use does not have.
///
/// `+Infinity` does terminate: it clears the asymptotic branch on the first test.
///
/// A negative integer, where the true function has a pole, recurses up through zero and comes back
/// finite.
pub fn digamma(x: f64) -> Result<f64, GammaError> {
    if x.is_nan() || x == f64::NEG_INFINITY {
        return Err(GammaError::NonTerminating);
    }
    Ok(digamma_terminating(x))
}

fn digamma_terminating(x: f64) -> f64 {
    if x > 0.0 && x <= S_LIMIT {
        // "use method 5 from Bernardo AS103, accurate to O(x)".
        return -GAMMA - 1.0 / x;
    }
    if x >= C_LIMIT {
        // "use method 4 (accurate to O(1/x^8))":
        //            1       1        1         1
        // log(x) -  --- - ------ + ------- - -------
        //           2 x   12 x^2   120 x^4   252 x^6
        let inv = 1.0 / (x * x);
        return crate::fast_math::log(x)
            - 0.5 / x
            - inv * ((1.0 / 12.0) + inv * (1.0 / 120.0 - inv / 252.0));
    }
    digamma_terminating(x + 1.0) - 1.0 / x
}

/// `Gamma.trigamma`: the derivative of [`digamma`], with the same three-branch shape.
///
/// # The asymptotic branch contradicts its own comment, and the comment is the correct one
///
/// ```java
/// //  1    1      1       1       1
/// //  - + ---- + ---- - ----- + -----
/// //  x      2      3       5       7
/// //      2 x    6 x    30 x    42 x
/// return 1 / x + inv / 2 + inv / x * (1.0 / 6 - inv * (1.0 / 30 + inv / 42));
/// ```
///
/// Expanding the code gives `1/(6x^3) - 1/(30x^5) - 1/(42x^7)`, whose last sign is not the comment's
/// and not the Bernoulli series'. The true asymptotic expansion alternates, so the comment is right
/// and the code is wrong, and the error is real rather than a last-bit one: at `x = 1`, which
/// recurses up to the asymptotic branch, the reference answers `1.6449340668481562` where the true
/// value is `pi^2 / 6 = 1.6449340668482264`, eleven digits in.
///
/// The port follows the **code**, because the code is what produced the golden. The comment-faithful
/// version was written first and the golden refused it, which is the only way this would have been
/// found: both forms look correct, and only one of them is the reference.
pub fn trigamma(x: f64) -> Result<f64, GammaError> {
    if x.is_nan() || x == f64::NEG_INFINITY {
        return Err(GammaError::NonTerminating);
    }
    Ok(trigamma_terminating(x))
}

fn trigamma_terminating(x: f64) -> f64 {
    if x > 0.0 && x <= S_LIMIT {
        return 1.0 / (x * x);
    }
    if x >= C_LIMIT {
        let inv = 1.0 / (x * x);
        //  1    1      1       1       1
        //  - + ---- + ---- - ----- + -----
        //  x      2      3       5       7
        //      2 x    6 x    30 x    42 x
        // `+ inv / 42`, as the code has it and not as the comment above it has it.
        return 1.0 / x + inv / 2.0 + inv / x * (1.0 / 6.0 - inv * (1.0 / 30.0 + inv / 42.0));
    }
    trigamma_terminating(x + 1.0) + 1.0 / (x * x)
}
