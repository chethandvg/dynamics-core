// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! # Dynamics core
//!
//! Leaky integrators and periodic pulses, with Lean 4 proofs of when repeated pulses add up and when they are forgotten.
//!
//! Part of [Frequencies of Life](https://chethandvg.github.io/frequencies-of-life/). Each public function corresponds to
//! claims in `claims.yaml`; the README lists what is proved and how.
//!
//! A linear leaky integrator is the simplest way for a cell to "count" pulses.
//! Its state `A ≥ 0` decays with time constant `τ` between pulses,
//! `A(t + Δt) = A(t)·e^(−Δt/τ)`, and each pulse adds `h` at once. With a pulse
//! every `T`, one period leaves the residual factor `r = e^(−T/τ)`, so the peaks
//! follow `A(n+1) = A(n)·r + h` and rise towards `h / (1 − r)`. The summation
//! ratio `1 / (1 − r)` says how much repeated pulses add up: 1 means not at all.
//!
//! ```
//! use dynamics_core::{pulses_to_fraction, summation_ratio};
//!
//! // Spikes 10 minutes apart, read by an integrator with τ = 75 s: r = e^-8.
//! let ratio = summation_ratio(600.0, 75.0)?;
//! assert!(1.0 < ratio && ratio < 1.0004);
//! // Pulses every τ/2 instead: the fifth peak is the first past 90 % of the steady state.
//! assert_eq!(pulses_to_fraction(0.9, 37.5, 75.0)?, 5);
//! # Ok::<(), dynamics_core::DynamicsError>(())
//! ```
//!
//! Times are in any one unit (seconds, say). Time constants, periods and pulse
//! sizes must be finite and greater than 0, elapsed times and states finite and
//! at least 0; anything else is a [`DynamicsError`]. So is a result too large for
//! an `f64`, while one too small for it rounds to 0. T/τ and Δt/τ are rounded to
//! doubles first, so a ratio below the normal range (under 2.2·10⁻³⁰⁸) carries
//! fewer digits, and a T/τ that rounds to 0 makes [`steady_state_peak`] and
//! [`steady_state_trough`] an error even where h·τ/T would fit.

use std::f64::consts::LN_2;
use std::fmt;

mod exact;

/// The largest count [`pulses_to_fraction`] returns: 2⁵³, up to which every
/// whole number is exact in an `f64` (and so in JavaScript).
pub const MAX_PULSES: i64 = 1 << 53;

const POSITIVE: &str = "a finite number greater than 0";
const NON_NEGATIVE: &str = "a finite number of at least 0";
const FRACTION: &str = "a number strictly between 0 and 1";
const WHOLE: &str = "a whole number of at least 1";
const FINITE: &str = "a finite number";
const COUNTABLE: &str = "a count of at most 2^53 pulses";

/// Below this many time constants, `n` pulses add up to `n·h` to within 2⁻⁶¹.
const NO_DECAY: f64 = 1.0 / (1u64 << 60) as f64;

/// An argument (or the result) outside the range the function accepts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DynamicsError {
    name: &'static str,
    value: f64,
    requirement: &'static str,
}

impl DynamicsError {
    /// The parameter (or `"result"`) that was out of range.
    pub fn name(&self) -> &'static str {
        self.name
    }

    /// The offending value.
    pub fn value(&self) -> f64 {
        self.value
    }
}

impl fmt::Display for DynamicsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} must be {}, got {}",
            self.name, self.requirement, self.value
        )
    }
}

impl std::error::Error for DynamicsError {}

/// `rhythm-units` checks the same thing under the same parameter names.
impl From<rhythm_units::RhythmError> for DynamicsError {
    fn from(e: rhythm_units::RhythmError) -> Self {
        Self {
            name: e.name(),
            value: e.value(),
            requirement: POSITIVE,
        }
    }
}

/// The result type of this crate.
pub type Result<T> = std::result::Result<T, DynamicsError>;

fn require(ok: bool, name: &'static str, value: f64, requirement: &'static str) -> Result<f64> {
    if ok {
        Ok(value)
    } else {
        Err(DynamicsError {
            name,
            value,
            requirement,
        })
    }
}

fn positive(name: &'static str, value: f64) -> Result<f64> {
    require(value.is_finite() && value > 0.0, name, value, POSITIVE)
}

fn non_negative(name: &'static str, value: f64) -> Result<f64> {
    require(value.is_finite() && value >= 0.0, name, value, NON_NEGATIVE)
}

fn finite_result(value: f64) -> Result<f64> {
    require(value.is_finite(), "result", value, FINITE)
}

/// T/τ, the number of time constants in one period, from checked arguments.
fn time_constants(period: f64, tau: f64) -> Result<f64> {
    Ok(positive("period", period)? / positive("tau", tau)?)
}

/// Δt/τ, the number of time constants elapsed, from checked arguments.
fn elapsed(dt: f64, tau: f64) -> Result<f64> {
    Ok(non_negative("dt", dt)? / positive("tau", tau)?)
}

/// 1 − e^(−y) for y ≥ 0, accurate for small y too.
fn one_minus_exp_neg(y: f64) -> f64 {
    -(-y).exp_m1()
}

/// The fraction of the state left after `dt` (≥ 0) with time constant `tau`: e^(−dt/τ).
pub fn decay_factor(dt: f64, tau: f64) -> Result<f64> {
    Ok((-elapsed(dt, tau)?).exp())
}

/// The state `a` (≥ 0) after relaxing for `dt` (≥ 0): a·e^(−dt/τ). Relaxing for
/// dt₁ and then dt₂ is relaxing for dt₁ + dt₂ (`dynamics:RELAX-COMPOSE`).
///
/// The product never forms e^(−dt/τ) on its own, so a large state keeps its
/// digits after more than 708 time constants, where that factor alone falls
/// below the normal `f64` range.
pub fn relax(a: f64, dt: f64, tau: f64) -> Result<f64> {
    let a = non_negative("a", a)?;
    Ok(exact::mul_exp_neg(a, elapsed(dt, tau)?))
}

/// The time in which the state halves: τ·ln 2.
pub fn half_life(tau: f64) -> Result<f64> {
    Ok(positive("tau", tau)? * LN_2)
}

/// The peak right after the `n`-th pulse (n ≥ 1), starting from 0, for pulses
/// of size `h` every `period`: h·(1 − rⁿ)/(1 − r) with r = e^(−T/τ)
/// (`dynamics:PEAK-AFTER-N`).
///
/// Both 1 − rⁿ and 1 − r come from `exp_m1`, so the result stays accurate when
/// the period is much shorter than τ.
pub fn peak_after_pulses(h: f64, period: f64, tau: f64, n: i64) -> Result<f64> {
    let h = positive("h", h)?;
    let x = time_constants(period, tau)?;
    let n = require(n >= 1, "n", n as f64, WHOLE)?;
    // Without decay, n pulses add up to n·h; using that below 2⁻⁶⁰ also keeps a
    // subnormal x from losing precision.
    let sum = if n * x < NO_DECAY {
        n
    } else {
        one_minus_exp_neg(n * x) / one_minus_exp_neg(x)
    };
    finite_result(h * sum)
}

/// A pulse count given as an `f64`, as JavaScript numbers are: a whole number from
/// 1 to [`MAX_PULSES`], up to which every count is exact, converted to `i64`.
/// Anything else is a [`DynamicsError`] for `n`, never a truncated count.
pub fn pulse_count(n: f64) -> Result<i64> {
    require(n.fract() == 0.0 && n >= 1.0, "n", n, WHOLE)?;
    Ok(require(n <= MAX_PULSES as f64, "n", n, COUNTABLE)? as i64)
}

/// How much repeated pulses add up: 1/(1 − r) with r = e^(−T/τ). It is the steady
/// peak in units of one pulse; 1 means each pulse is forgotten before the next.
/// It falls towards 1 as the period grows (`dynamics:SUMMATION-MONOTONE`).
pub fn summation_ratio(period: f64, tau: f64) -> Result<f64> {
    steady_state_peak(1.0, period, tau)
}

/// [`summation_ratio`] for pulses at `hz` pulses per unit of time, with the
/// period from [`rhythm_units::period_from_hz`].
pub fn summation_ratio_at_hz(hz: f64, tau: f64) -> Result<f64> {
    summation_ratio(rhythm_units::period_from_hz(hz)?, tau)
}

/// The peak that periodic pulses approach, h/(1 − r): the limit of
/// [`peak_after_pulses`] as n grows (`dynamics:STEADY-STATE`).
pub fn steady_state_peak(h: f64, period: f64, tau: f64) -> Result<f64> {
    finite_result(positive("h", h)? / one_minus_exp_neg(time_constants(period, tau)?))
}

/// The steady state just before each pulse, h·r/(1 − r): the steady peak relaxed
/// for one period (`dynamics:STEADY-TROUGH`).
///
/// Below T/τ = 709 it is h/(e^(T/τ) − 1), from `exp_m1` and one division, so a
/// subnormal h is not rounded onto the subnormal grid before it is divided.
/// Above, where 1 − r is 1 in `f64`, it is h·r formed as in [`relax`], so a large
/// h keeps its digits although r alone is below the normal `f64` range.
pub fn steady_state_trough(h: f64, period: f64, tau: f64) -> Result<f64> {
    let h = positive("h", h)?;
    let x = time_constants(period, tau)?;
    finite_result(if x < 709.0 {
        h / x.exp_m1()
    } else {
        exact::mul_exp_neg(h, x)
    })
}

/// The smallest number of pulses n ≥ 1 whose peak reaches the fraction `p`
/// (0 < p < 1) of the steady-state peak: the smallest n with 1 − rⁿ ≥ p.
///
/// A floating-point formula such as ⌈ln(1 − p)/ln r⌉ is off by one whenever the
/// answer sits on a boundary. Here a first guess is checked against the
/// definition in double-double arithmetic (about 106 bits), so the answer could
/// only be wrong if 1 − e^(−n·T/τ) and p agreed to about 2⁻¹⁰⁰ relative; no input
/// on which it fails is known (the tie p = T/τ = 2⁻³⁰⁰ agrees far more closely and
/// is decided correctly). Counts above [`MAX_PULSES`] are an error.
pub fn pulses_to_fraction(p: f64, period: f64, tau: f64) -> Result<i64> {
    let p = require(p > 0.0 && p < 1.0, "p", p, FRACTION)?;
    let goal = exact::Threshold::new(p, positive("period", period)?, positive("tau", tau)?);
    let too_many = |value| DynamicsError {
        name: "result",
        value,
        requirement: COUNTABLE,
    };
    let guess = goal.estimate().ceil().max(1.0);
    if guess > MAX_PULSES as f64 {
        return Err(too_many(guess));
    }
    let mut n = guess as i64;
    while n > 1 && goal.reached(n - 1) {
        n -= 1;
    }
    while !goal.reached(n) {
        n += 1;
        if n > MAX_PULSES {
            return Err(too_many(n as f64));
        }
    }
    Ok(n)
}

/// The shape of one calcium spike `s` (≥ 0) after it starts, with rise and decay
/// time constants: (1 − e^(−s/rise))·e^(−s/decay). It is 0 at the start, peaks
/// below 1, and decays back to 0.
pub fn spike_shape(s: f64, rise: f64, decay: f64) -> Result<f64> {
    let s = non_negative("s", s)?;
    let (rise, decay) = (positive("rise", rise)?, positive("decay", decay)?);
    Ok(one_minus_exp_neg(s / rise) * (-s / decay).exp())
}

/// When [`spike_shape`] peaks: rise·ln(1 + decay/rise) (`dynamics:SPIKE-PEAK`).
/// It always lies between 0 and `decay`.
pub fn spike_peak_time(rise: f64, decay: f64) -> Result<f64> {
    let (rise, decay) = (positive("rise", rise)?, positive("decay", decay)?);
    let z = decay / rise;
    Ok(if z.is_infinite() {
        rise * (decay.ln() - rise.ln()) // 1 + z overflows; ln(1 + z) = ln z here
    } else if z < f64::MIN_POSITIVE {
        decay // rise·ln(1 + z) = decay·(1 − z/2 + …) = decay in f64
    } else {
        rise * z.ln_1p()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn errors_name_the_argument_and_the_requirement() {
        let e = relax(-1.0, 1.0, 1.0).unwrap_err();
        assert_eq!((e.name(), e.value()), ("a", -1.0));
        assert_eq!(
            e.to_string(),
            "a must be a finite number of at least 0, got -1"
        );
        assert_eq!(
            pulses_to_fraction(1.0, 1.0, 1.0).unwrap_err().to_string(),
            "p must be a number strictly between 0 and 1, got 1"
        );
        assert_eq!(peak_after_pulses(1.0, 1.0, 1.0, 0).unwrap_err().name(), "n");
        assert_eq!(summation_ratio_at_hz(0.0, 75.0).unwrap_err().name(), "hz");
        assert_eq!(
            summation_ratio_at_hz(1e-310, 75.0).unwrap_err().name(),
            "result"
        );
    }

    #[test]
    fn pulse_counts_from_f64_are_whole_or_rejected() {
        assert_eq!(pulse_count(2.0), Ok(2));
        assert_eq!(pulse_count(MAX_PULSES as f64), Ok(MAX_PULSES));
        for bad in [2.5, 0.0, -3.0, 0.5, f64::NAN, f64::INFINITY] {
            let e = pulse_count(bad).unwrap_err();
            assert_eq!(e.name(), "n");
            assert!(
                e.to_string()
                    .starts_with("n must be a whole number of at least 1")
            );
        }
        assert_eq!(
            pulse_count(2f64.powi(60)).unwrap_err().to_string(),
            "n must be a count of at most 2^53 pulses, got 1152921504606847000"
        );
    }

    #[test]
    fn products_with_an_underflowing_factor_keep_their_digits() {
        // e^-750 and e^-720 alone are 0 and subnormal; the products are normal doubles.
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-15 * b;
        assert!(close(
            relax(1e300, 750.0, 1.0).unwrap(),
            1.901_684_963_475_006_5e-26
        ));
        assert!(close(
            relax(1e10, 720.0, 1.0).unwrap(),
            2.032_230_802_424_293e-303
        ));
        assert!(close(
            steady_state_trough(1e300, 750.0, 1.0).unwrap(),
            1.901_684_963_475_006_5e-26
        ));
        // A subnormal h is divided once, not rounded onto the subnormal grid first.
        assert!(close(
            steady_state_trough(2.3e-316, 1e-8, 1.0).unwrap(),
            2.299_999_975_620_625_6e-308
        ));
        assert_eq!(relax(f64::MAX, 1600.0, 1.0).unwrap(), 0.0);
        // The summation ratio alone overflows here; the steady peak does not.
        assert!(steady_state_peak(1e-20, 1e-320, 1.0).unwrap().is_finite());
    }

    #[test]
    fn tiny_periods_add_up_without_decay() {
        assert_eq!(peak_after_pulses(1.0, 1e-320, 1.0, 1000).unwrap(), 1000.0);
        assert_eq!(peak_after_pulses(2.0, 1e-30, 1.0, 7).unwrap(), 14.0);
        assert_eq!(summation_ratio(1e-320, 1.0).unwrap_err().name(), "result");
    }

    #[test]
    fn spike_peak_time_covers_extreme_ratios() {
        assert_eq!(spike_peak_time(1e300, 1e-10).unwrap(), 1e-10);
        let t = spike_peak_time(1e-300, 1e300).unwrap();
        assert!((t - 1e-300 * 600.0 * 10f64.ln()).abs() < 1e-312);
    }
}
