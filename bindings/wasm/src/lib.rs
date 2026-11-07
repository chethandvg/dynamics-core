// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! WebAssembly bindings for `dynamics-core`. Every function keeps the Rust name, so
//! the golden vectors can call it by name; invalid input throws a JS `Error`.
//! Pulse counts cross the boundary as `f64`, so JavaScript never sees a `BigInt`,
//! and a count that is not a whole number from 1 to 2⁵³ throws instead of being
//! truncated.

use dynamics_core_rs as imp;
use wasm_bindgen::prelude::*;

/// The crate version.
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_owned()
}

/// Declares one exported JS function per Rust function of `f64`s.
macro_rules! bind {
    ($($doc:literal $name:ident($($arg:ident),+);)+) => {
        $(
            #[doc = $doc]
            #[wasm_bindgen]
            pub fn $name($($arg: f64),+) -> Result<f64, JsError> {
                Ok(imp::$name($($arg),+)?)
            }
        )+
    };
}

bind! {
    "The fraction of the state left after `dt` with time constant `tau`: e^(−dt/τ)." decay_factor(dt, tau);
    "The state `a` after relaxing for `dt`: a·e^(−dt/τ)." relax(a, dt, tau);
    "The time in which the state halves: τ·ln 2." half_life(tau);
    "The peak that periodic pulses approach: h/(1 − r)." steady_state_peak(h, period, tau);
    "The steady state just before each pulse: h·r/(1 − r)." steady_state_trough(h, period, tau);
    "How much repeated pulses add up: 1/(1 − r), with r = e^(−T/τ)." summation_ratio(period, tau);
    "`summation_ratio` for pulses at `hz` pulses per unit of time." summation_ratio_at_hz(hz, tau);
    "The spike shape `s` after it starts: (1 − e^(−s/rise))·e^(−s/decay)." spike_shape(s, rise, decay);
    "When `spike_shape` peaks: rise·ln(1 + decay/rise)." spike_peak_time(rise, decay);
}

/// The peak right after the `n`-th pulse (n ≥ 1), starting from 0: h·(1 − rⁿ)/(1 − r).
#[wasm_bindgen]
pub fn peak_after_pulses(h: f64, period: f64, tau: f64, n: f64) -> Result<f64, JsError> {
    Ok(imp::peak_after_pulses(
        h,
        period,
        tau,
        imp::pulse_count(n)?,
    )?)
}

/// The smallest n ≥ 1 whose peak reaches the fraction `p` of the steady peak,
/// exact at the boundaries. At most 2⁵³, so it is exact as a JS number.
#[wasm_bindgen]
pub fn pulses_to_fraction(p: f64, period: f64, tau: f64) -> Result<f64, JsError> {
    Ok(imp::pulses_to_fraction(p, period, tau)? as f64)
}
