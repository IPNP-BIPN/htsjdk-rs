//! Double-double arithmetic: an unevaluated sum `hi + lo` carrying roughly 106 bits.
//!
//! Used to compute a function to well beyond `f64` precision so the final rounding to `f64` is
//! the *only* rounding, which is what "correctly rounded" means. See
//! `docs/decisions/0006-correct-rounding-is-the-target-for-log-and-log10.md`.
//!
//! All products rely on `f64::mul_add` being a genuine fused multiply-add. On any target where
//! it is emulated in software the results stay correct but slow.

/// An unevaluated sum where `|lo| <= ulp(hi) / 2`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoubleDouble {
    pub hi: f64,
    pub lo: f64,
}

impl DoubleDouble {
    pub const ZERO: Self = Self { hi: 0.0, lo: 0.0 };

    #[inline]
    pub const fn new(hi: f64, lo: f64) -> Self {
        Self { hi, lo }
    }

    #[inline]
    pub const fn from_f64(x: f64) -> Self {
        Self { hi: x, lo: 0.0 }
    }

    /// Rounds to the nearest `f64`. This is the single rounding the whole exercise protects.
    #[inline]
    pub fn to_f64(self) -> f64 {
        self.hi + self.lo
    }
}

/// Exact sum of two `f64`. Knuth's algorithm; no assumption about relative magnitude.
#[inline]
pub fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    (s, err)
}

/// Exact sum when `|a| >= |b|` is already known. Cheaper than [`two_sum`].
#[inline]
pub fn quick_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let err = b - (s - a);
    (s, err)
}

/// Exact product, using FMA to recover the rounding error.
#[inline]
pub fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    let err = a.mul_add(b, -p);
    (p, err)
}

#[inline]
pub fn add(a: DoubleDouble, b: DoubleDouble) -> DoubleDouble {
    let (s, e) = two_sum(a.hi, b.hi);
    let (s2, e2) = two_sum(a.lo, b.lo);
    let (h, l) = quick_two_sum(s, e + s2);
    let (h2, l2) = quick_two_sum(h, l + e2);
    DoubleDouble::new(h2, l2)
}

#[inline]
pub fn add_f64(a: DoubleDouble, b: f64) -> DoubleDouble {
    let (s, e) = two_sum(a.hi, b);
    let (h, l) = quick_two_sum(s, e + a.lo);
    DoubleDouble::new(h, l)
}

#[inline]
pub fn sub(a: DoubleDouble, b: DoubleDouble) -> DoubleDouble {
    add(a, DoubleDouble::new(-b.hi, -b.lo))
}

#[inline]
pub fn mul(a: DoubleDouble, b: DoubleDouble) -> DoubleDouble {
    let (p, e) = two_prod(a.hi, b.hi);
    let e = e + (a.hi * b.lo + a.lo * b.hi);
    let (h, l) = quick_two_sum(p, e);
    DoubleDouble::new(h, l)
}

#[inline]
pub fn mul_f64(a: DoubleDouble, b: f64) -> DoubleDouble {
    let (p, e) = two_prod(a.hi, b);
    let (h, l) = quick_two_sum(p, e + a.lo * b);
    DoubleDouble::new(h, l)
}

#[inline]
pub fn div(a: DoubleDouble, b: DoubleDouble) -> DoubleDouble {
    // One Newton correction on the f64 quotient recovers the extra ~53 bits.
    let q1 = a.hi / b.hi;
    let r = sub(a, mul_f64(b, q1));
    let q2 = r.hi / b.hi;
    let r2 = sub(r, mul_f64(b, q2));
    let q3 = r2.hi / b.hi;
    let (h, l) = quick_two_sum(q1, q2);
    add_f64(DoubleDouble::new(h, l), q3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_sum_is_exact() {
        let (s, e) = two_sum(1.0, 1e-30);
        assert_eq!(s, 1.0);
        assert_eq!(e, 1e-30);
    }

    #[test]
    fn two_prod_recovers_the_rounding_error() {
        let a = 1.0 + f64::EPSILON;
        let (p, e) = two_prod(a, a);
        // The exact product needs more than 53 bits, so the error term must be non-zero.
        assert_ne!(e, 0.0);
        assert_eq!(p, a * a);
    }

    #[test]
    fn multiplication_beats_plain_f64() {
        // 1/3 squared: the double-double result must round-trip more accurately than f64.
        let third = DoubleDouble::new(1.0 / 3.0, 1.849_396_431_339_59e-17);
        let sq = mul(third, third);
        let naive = (1.0f64 / 3.0) * (1.0 / 3.0);
        let exact = 1.0f64 / 9.0;
        assert!((sq.to_f64() - exact).abs() <= (naive - exact).abs());
    }
}
