// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! Python bindings for `dynamics-core`. Every function keeps the Rust name, so the
//! golden vectors can call it by name; invalid input raises `ValueError`.

use dynamics_core_rs as imp;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Declares one Python function per Rust function, then a `register` that adds them all.
macro_rules! bind {
    ($($doc:literal $name:ident($($arg:ident: $ty:ty),+) -> $ret:ty;)+) => {
        $(
            #[doc = $doc]
            #[pyfunction]
            fn $name($($arg: $ty),+) -> PyResult<$ret> {
                imp::$name($($arg),+).map_err(|e| PyValueError::new_err(e.to_string()))
            }
        )+

        fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
            $(m.add_function(wrap_pyfunction!($name, m)?)?;)+
            Ok(())
        }
    };
}

bind! {
    "The fraction of the state left after `dt` with time constant `tau`: exp(-dt / tau)."
        decay_factor(dt: f64, tau: f64) -> f64;
    "The state `a` after relaxing for `dt`: a * exp(-dt / tau)." relax(a: f64, dt: f64, tau: f64) -> f64;
    "The time in which the state halves: tau * ln 2." half_life(tau: f64) -> f64;
    "The peak right after the `n`-th pulse (n >= 1), starting from 0: h * (1 - r**n) / (1 - r)."
        peak_after_pulses(h: f64, period: f64, tau: f64, n: i64) -> f64;
    "The peak that periodic pulses approach: h / (1 - r)." steady_state_peak(h: f64, period: f64, tau: f64) -> f64;
    "The steady state just before each pulse: h * r / (1 - r)."
        steady_state_trough(h: f64, period: f64, tau: f64) -> f64;
    "How much repeated pulses add up: 1 / (1 - r), with r = exp(-period / tau)."
        summation_ratio(period: f64, tau: f64) -> f64;
    "summation_ratio for pulses at `hz` pulses per unit of time." summation_ratio_at_hz(hz: f64, tau: f64) -> f64;
    "The smallest n >= 1 whose peak reaches the fraction `p` of the steady peak, exact at the boundaries."
        pulses_to_fraction(p: f64, period: f64, tau: f64) -> i64;
    "The spike shape `s` after it starts: (1 - exp(-s / rise)) * exp(-s / decay)."
        spike_shape(s: f64, rise: f64, decay: f64) -> f64;
    "When spike_shape peaks: rise * ln(1 + decay / rise)." spike_peak_time(rise: f64, decay: f64) -> f64;
}

#[pymodule]
fn dynamics_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    register(m)
}
