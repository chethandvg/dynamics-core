// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! Exact binary scaling, for the two places where plain `f64` arithmetic would
//! underflow on the way to a result that fits.
//!
//! - The threshold behind [`crate::pulses_to_fraction`]: does `n` pulses reach
//!   the fraction `p`, that is, is 1 − e^(−y) ≥ p with y = n·T/τ? The decision
//!   is made in double-double arithmetic: a value is the unevaluated sum
//!   `hi + lo` of two doubles with |lo| ≤ ulp(hi)/2, about 106 significant bits.
//!   The inputs are split into significand and exponent first, so nothing
//!   underflows however small p or T/τ is.
//! - [`mul_exp_neg`], the product a·e^(−x) behind [`crate::relax`] and
//!   [`crate::steady_state_trough`].

use std::cmp::Ordering;

/// A double-double number `hi + lo`, normalised so that `hi` is `hi + lo` rounded.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Dd {
    hi: f64,
    lo: f64,
}

/// ln 2 to about 107 bits.
const LN2: Dd = Dd {
    hi: std::f64::consts::LN_2,
    lo: 2.319_046_813_846_299_6e-17,
};

/// a + b exactly (Knuth's two-sum).
fn two_sum(a: f64, b: f64) -> Dd {
    let s = a + b;
    let v = s - a;
    Dd {
        hi: s,
        lo: (a - (s - v)) + (b - v),
    }
}

/// a + b exactly, when |a| ≥ |b| (Dekker's fast two-sum).
fn fast_two_sum(a: f64, b: f64) -> Dd {
    let s = a + b;
    Dd {
        hi: s,
        lo: b - (s - a),
    }
}

/// a · b exactly, with one fused multiply-add.
fn two_prod(a: f64, b: f64) -> Dd {
    let p = a * b;
    Dd {
        hi: p,
        lo: a.mul_add(b, -p),
    }
}

impl Dd {
    fn neg(self) -> Self {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }

    fn add(self, o: Self) -> Self {
        let s = two_sum(self.hi, o.hi);
        let t = two_sum(self.lo, o.lo);
        let u = fast_two_sum(s.hi, s.lo + t.hi);
        fast_two_sum(u.hi, u.lo + t.lo)
    }

    fn add_f64(self, b: f64) -> Self {
        let s = two_sum(self.hi, b);
        fast_two_sum(s.hi, s.lo + self.lo)
    }

    fn mul(self, o: Self) -> Self {
        let p = two_prod(self.hi, o.hi);
        fast_two_sum(p.hi, p.lo + (self.hi * o.lo + self.lo * o.hi))
    }

    fn mul_f64(self, b: f64) -> Self {
        let p = two_prod(self.hi, b);
        fast_two_sum(p.hi, p.lo + self.lo * b)
    }

    fn div_f64(self, b: f64) -> Self {
        let q = self.hi / b;
        let p = two_prod(q, b);
        let s = two_sum(self.hi, -p.hi);
        fast_two_sum(q, (s.hi + (s.lo - p.lo + self.lo)) / b)
    }

    /// · 2ᵏ, exact while the parts stay normal.
    fn scale(self, k: i64) -> Self {
        let f = pow2(k);
        Self {
            hi: self.hi * f,
            lo: self.lo * f,
        }
    }

    /// Compares `hi + lo` with `x` exactly: because `hi` is `hi + lo` rounded,
    /// the sum is on the same side of the double `x` as `hi`, unless `hi = x`.
    fn cmp_f64(self, x: f64) -> Ordering {
        match self.hi.partial_cmp(&x) {
            Some(Ordering::Equal) => self.lo.partial_cmp(&0.0).unwrap_or(Ordering::Equal),
            other => other.unwrap_or(Ordering::Equal),
        }
    }
}

/// e^w − 1 for |w| ≤ 0.75, to about 2⁻¹⁰⁰ relative: the Taylor series to the
/// 10th power at z = w/2¹⁰, then ten doublings e^(2z) − 1 = E·(E + 2).
fn expm1(w: Dd) -> Dd {
    const HALVINGS: i64 = 10;
    let z = w.scale(-HALVINGS);
    let mut t = Dd { hi: 1.0, lo: 0.0 };
    for k in (2..=10).rev() {
        t = t.mul(z).div_f64(f64::from(k)).add_f64(1.0);
    }
    let mut e = t.mul(z);
    for _ in 0..HALVINGS {
        e = e.mul(e.add_f64(2.0));
    }
    e
}

/// 2ᵏ, exact for every `k`: 0 below the `f64` range, infinity above it.
fn pow2(k: i64) -> f64 {
    match k {
        -1022..=1023 => f64::from_bits(((k + 1023) as u64) << 52),
        -1074..=-1023 => f64::from_bits(1 << (k + 1074)),
        ..-1074 => 0.0,
        _ => f64::INFINITY,
    }
}

/// Splits a finite x > 0 exactly into (m, e) with x = m·2ᵉ and 1 ≤ m < 2.
fn split(x: f64) -> (f64, i64) {
    const FRACTION: u64 = (1 << 52) - 1;
    let bits = x.to_bits();
    match (bits >> 52) as i64 {
        0 => {
            let (m, e) = split(x * pow2(64)); // subnormal: scale up exactly
            (m, e - 64)
        }
        biased => (
            f64::from_bits((bits & FRACTION) | (1023 << 52)),
            biased - 1023,
        ),
    }
}

/// a·e^(−x) for finite a ≥ 0 and x ≥ 0, without forming e^(−x) on its own, which
/// leaves the normal range past x ≈ 708 even when the product is far inside it.
/// With a = m·2ᵏ (1 ≤ m < 2) and x = j·ln 2 + r (|r| ≤ ln 2/2), it is
/// m·e^(−r)·2^(k−j): r is within about 2⁻⁵³ of its true value and the power of 2
/// is exact, so the result is as accurate as a·e^(−x) is for small x.
pub(crate) fn mul_exp_neg(a: f64, x: f64) -> f64 {
    if a == 0.0 || x > 1500.0 {
        return 0.0; // past x = 1500 the product is below 2⁻¹¹⁴⁰ for every finite a
    }
    let (m, k) = split(a);
    let j = (x / LN2.hi).round();
    let r = (-j).mul_add(LN2.hi, x) - j * LN2.lo;
    let v = m * (-r).exp(); // in [0.7, 2.9), so v·2ˢ cannot overflow: it is at most a
    let s = k - j as i64;
    if s >= -1022 {
        v * pow2(s)
    } else {
        v * pow2(s + 64) * pow2(-64) // exact, then one rounding into the subnormal range
    }
}

/// The question "do n pulses reach the fraction p?" for one p, T and τ.
pub(crate) struct Threshold {
    p: f64,
    p_parts: (f64, i64),
    period: (f64, i64),
    tau: (f64, i64),
}

impl Threshold {
    /// For checked arguments: 0 < p < 1, period and tau finite and positive.
    pub(crate) fn new(p: f64, period: f64, tau: f64) -> Self {
        Self {
            p,
            p_parts: split(p),
            period: split(period),
            tau: split(tau),
        }
    }

    /// −ln(1 − p)·τ/T in floating point: within a few units in the last place of
    /// the answer, 0 if it is below 2⁻⁶⁰ and infinity if it is above 2⁶⁰.
    pub(crate) fn estimate(&self) -> f64 {
        let (l, el) = split(-(-self.p).ln_1p());
        let k = el - self.period.1 + self.tau.1;
        match k {
            ..-60 => 0.0,
            61.. => f64::INFINITY,
            _ => l * self.tau.0 / self.period.0 * pow2(k),
        }
    }

    /// Whether n pulses reach the fraction p: 1 − e^(−y) ≥ p, y = n·T/τ.
    pub(crate) fn reached(&self, n: i64) -> bool {
        // y = big·2ˢ, with big = n·mT/mτ in [n/2, 2n) and no risk of underflow.
        let big = two_prod(n as f64, self.period.0).div_f64(self.tau.0);
        let s = self.period.1 - self.tau.1;
        let e = split(big.hi).1 + s; // 2ᵉ ≤ y < 2ᵉ⁺¹
        if e >= 6 {
            return true; // y ≥ 64, and −ln(1 − p) ≤ 53·ln 2 < 37 for every f64 p < 1
        }
        if e < -61 {
            // y < 2⁻⁶⁰: 1 − e^(−y) = y·(1 − y/2) to within y²/6 < 2⁻¹²¹ relative,
            // compared with p after scaling both by 2⁻ˢ.
            let k = self.p_parts.1 - s;
            return match k {
                57.. => false, // p·2⁻ˢ ≥ 2⁵⁷ > big
                ..-3 => true,  // p·2⁻ˢ < 1/8 < big·(1 − y/2)
                _ => {
                    let half_y = 0.5 * big.hi * pow2(s);
                    let lhs = big.add_f64(-big.hi * half_y);
                    lhs.cmp_f64(self.p_parts.0 * pow2(k)) != Ordering::Less
                }
            };
        }
        let y = big.scale(s); // 2⁻⁶² ≤ y < 64: every part is normal
        let p = self.p;
        if p < 0.5 {
            // An answer near p < 1/2 needs y < ln 2: compare 1 − e^(−y) with p.
            y.hi >= 0.75 || expm1(y.neg()).neg().cmp_f64(p) != Ordering::Less
        } else {
            // Compare e^(−y) with 1 − p, which is exact here, as 2⁻ᵏ·e^(−r), r = y − k·ln 2.
            if y.hi >= 40.0 {
                return true; // e^(−40) < 2⁻⁵³ ≤ 1 − p
            }
            if y.hi < 0.6 {
                return false; // e^(−0.6) > 1/2 ≥ 1 − p
            }
            let k = (y.hi / LN2.hi).round();
            let r = y.add(LN2.mul_f64(-k));
            let e_neg_r = expm1(r.neg()).add_f64(1.0);
            e_neg_r.cmp_f64((1.0 - p) * pow2(k as i64)) != Ordering::Greater
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn value(x: Dd) -> f64 {
        x.hi + x.lo
    }

    #[test]
    fn double_double_keeps_what_f64_rounds_away() {
        let third = Dd { hi: 1.0, lo: 0.0 }.div_f64(3.0);
        let back = third.mul_f64(3.0);
        assert_eq!(back.cmp_f64(1.0), Ordering::Equal);
        assert_eq!(back.hi, 1.0);
        let tiny = Dd { hi: 1.0, lo: 0.0 }.add_f64(1e-30);
        assert_eq!(tiny.cmp_f64(1.0), Ordering::Greater);
        assert_eq!(tiny.neg().add_f64(1.0).hi, -1e-30);
    }

    #[test]
    fn expm1_matches_known_digits() {
        // e^(-1/2) - 1 = -0.39346934028736657639620046500882, which is
        // -0.3934693402873666 + -6.593178415491414e-19 as a double-double.
        let e = expm1(Dd { hi: -0.5, lo: 0.0 });
        assert_eq!(e.hi, -0.393_469_340_287_366_6);
        assert!((e.lo - -6.593_178_415_491_414e-19).abs() < 1e-32); // 2^-106 of 0.39
        // e^(1e-20) - 1 = 1e-20 + 5e-41: the second term survives in lo.
        let small = expm1(Dd { hi: 1e-20, lo: 0.0 });
        assert_eq!(small.hi, 1e-20);
        assert!((small.lo - 5e-41).abs() < 1e-55);
        assert!((value(expm1(Dd { hi: 0.75, lo: 0.0 })) - 0.75f64.exp_m1()).abs() < 4e-16);
    }

    #[test]
    fn split_and_pow2_are_exact() {
        assert_eq!(split(1.0), (1.0, 0));
        assert_eq!(split(600.0), (1.171875, 9));
        assert_eq!(split(f64::from_bits(1)), (1.0, -1074));
        assert_eq!(
            (pow2(-1074), pow2(-1075), pow2(1024)),
            (f64::from_bits(1), 0.0, f64::INFINITY)
        );
    }

    #[test]
    fn mul_exp_neg_matches_the_plain_product_where_that_is_exact() {
        for x in [0.0f64, 1e-300, 0.3, 1.0, 40.0, 700.0] {
            for a in [1.0, 2.5, 1e-300, 1e300] {
                let (want, got) = (a * (-x).exp(), mul_exp_neg(a, x));
                assert!((got - want).abs() <= 4e-16 * want, "{a}·e^-{x}");
            }
        }
        // A subnormal result rounds to the nearest subnormal.
        assert_eq!(mul_exp_neg(f64::from_bits(3), 1e-3), f64::from_bits(3));
        assert_eq!(mul_exp_neg(1.0, 745.0), (-745.0f64).exp());
        assert_eq!(mul_exp_neg(0.0, 1.0), 0.0);
    }

    #[test]
    fn an_exact_tie_is_not_reached() {
        // p = T/τ = 2⁻³⁰⁰: one pulse gives 1 − e^(−p) < p, two give more than p.
        let p = pow2(-300);
        let goal = Threshold::new(p, p, 1.0);
        assert_eq!(goal.estimate(), 1.0);
        assert!(!goal.reached(1));
        assert!(goal.reached(2));
    }
}
