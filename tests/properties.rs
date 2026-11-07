// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! Property tests over many inputs, for what single golden vectors cannot show:
//! NaN and infinity (which JSON cannot hold), identities and invariants.

use dynamics_core::*;

/// Ratios from 1e-3 to 1e3, spread evenly on a log scale.
fn ratios() -> impl DoubleEndedIterator<Item = f64> + Clone {
    (-30..=30).map(|i| 10f64.powf(f64::from(i) / 10.0))
}

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * b.abs()
}

/// Whether a function rejects its arguments.
type Rejects = fn(&[f64]) -> bool;

#[test]
fn rejects_nan_and_infinity_everywhere() {
    // Each function with valid arguments; every position is made bad in turn.
    let functions: &[(&str, &[f64], Rejects)] = &[
        ("decay_factor", &[1.0, 1.0], |a| {
            decay_factor(a[0], a[1]).is_err()
        }),
        ("relax", &[1.0, 1.0, 1.0], |a| {
            relax(a[0], a[1], a[2]).is_err()
        }),
        ("half_life", &[1.0], |a| half_life(a[0]).is_err()),
        ("peak_after_pulses", &[1.0, 1.0, 1.0, 3.0], |a| {
            pulse_count(a[3])
                .and_then(|n| peak_after_pulses(a[0], a[1], a[2], n))
                .is_err()
        }),
        ("steady_state_peak", &[1.0, 1.0, 1.0], |a| {
            steady_state_peak(a[0], a[1], a[2]).is_err()
        }),
        ("steady_state_trough", &[1.0, 1.0, 1.0], |a| {
            steady_state_trough(a[0], a[1], a[2]).is_err()
        }),
        ("summation_ratio", &[1.0, 1.0], |a| {
            summation_ratio(a[0], a[1]).is_err()
        }),
        ("summation_ratio_at_hz", &[1.0, 1.0], |a| {
            summation_ratio_at_hz(a[0], a[1]).is_err()
        }),
        ("pulses_to_fraction", &[0.5, 1.0, 1.0], |a| {
            pulses_to_fraction(a[0], a[1], a[2]).is_err()
        }),
        ("spike_shape", &[1.0, 1.0, 1.0], |a| {
            spike_shape(a[0], a[1], a[2]).is_err()
        }),
        ("spike_peak_time", &[1.0, 1.0], |a| {
            spike_peak_time(a[0], a[1]).is_err()
        }),
    ];
    for &(name, good, rejects) in functions {
        assert!(!rejects(good), "{name} accepts {good:?}");
        for i in 0..good.len() {
            for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut args = good.to_vec();
                args[i] = bad;
                assert!(rejects(&args), "{name} rejects {args:?}");
            }
        }
    }
}

#[test]
fn relaxing_twice_is_relaxing_for_the_sum() {
    for (dt1, dt2) in ratios().zip(ratios().rev()) {
        let twice = relax(relax(2.5, dt1, 3.0).unwrap(), dt2, 3.0).unwrap();
        let once = relax(2.5, dt1 + dt2, 3.0).unwrap();
        // e^(-x) amplifies rounding of x by x, here up to about 300.
        assert!(close(twice, once, 1e-13), "{dt1} + {dt2}");
    }
    assert_eq!(decay_factor(half_life(75.0).unwrap(), 75.0).unwrap(), 0.5);
}

#[test]
fn peaks_follow_the_recurrence_and_approach_the_steady_state() {
    for x in ratios() {
        let (mut a, steady) = (0.0, steady_state_peak(1.0, x, 1.0).unwrap());
        let r = decay_factor(x, 1.0).unwrap();
        for n in 1..=200 {
            a = a * r + 1.0;
            let peak = peak_after_pulses(1.0, x, 1.0, n).unwrap();
            assert!(close(peak, a, 1e-12), "T/tau = {x}, n = {n}");
            assert!(peak <= steady * (1.0 + 1e-15));
        }
        let trough = steady_state_trough(1.0, x, 1.0).unwrap();
        assert!(
            close(trough + 1.0, steady, 1e-14),
            "trough + h = peak at {x}"
        );
    }
}

#[test]
fn pulses_to_fraction_is_the_first_peak_past_the_fraction() {
    for x in ratios() {
        let steady = steady_state_peak(1.0, x, 1.0).unwrap();
        for p in [1e-9, 0.01, 0.1, 0.5, 0.632, 0.9, 0.99, 0.999_999] {
            let n = pulses_to_fraction(p, x, 1.0).unwrap();
            let peak = |n| peak_after_pulses(1.0, x, 1.0, n).unwrap() / steady;
            // Floating point decides only away from the boundary, hence the slack.
            assert!(
                peak(n) >= p * (1.0 - 1e-12),
                "p = {p}, T/tau = {x}: n = {n}"
            );
            assert!(
                n == 1 || peak(n - 1) < p * (1.0 + 1e-12),
                "p = {p}, T/tau = {x}: n = {n}"
            );
        }
    }
}

#[test]
fn summation_ratio_falls_towards_one() {
    let mut last = f64::INFINITY;
    for x in ratios() {
        let ratio = summation_ratio(x, 1.0).unwrap();
        assert!(1.0 <= ratio && (ratio < last || ratio == 1.0), "{x}"); // 1 in f64 past T = 37τ
        last = ratio;
    }
    assert_eq!(summation_ratio(800.0, 1.0).unwrap(), 1.0);
    assert_eq!(
        summation_ratio_at_hz(1.0 / 600.0, 75.0).unwrap(),
        summation_ratio(600.0, 75.0).unwrap()
    );
}

#[test]
fn spike_peak_time_is_the_maximum_of_the_shape() {
    for rise in ratios().step_by(6) {
        for decay in ratios().step_by(6) {
            let t = spike_peak_time(rise, decay).unwrap();
            assert!(0.0 < t && t <= decay, "rise {rise}, decay {decay}");
            let top = spike_shape(t, rise, decay).unwrap();
            for i in 0..=2000 {
                let s = t * f64::from(i) / 1000.0; // from 0 to twice the peak time
                assert!(spike_shape(s, rise, decay).unwrap() <= top * (1.0 + 1e-15));
            }
        }
    }
}
