// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! Golden-vector tests. The Python and WebAssembly bindings check the same
//! `vectors/*.json` cases, so all three implementations give the same answers.

use dynamics_core as d;
use golden_vectors::{Outcome, f64_arg, i64_arg};
use serde_json::Value;

/// Maps each vector's `fn` name to this crate's implementation.
fn dispatch(name: &str, args: &[Value]) -> Outcome {
    let x = |i| f64_arg(args, i);
    match name {
        "decay_factor" => d::decay_factor(x(0), x(1)).into(),
        "relax" => d::relax(x(0), x(1), x(2)).into(),
        "half_life" => d::half_life(x(0)).into(),
        "peak_after_pulses" => d::peak_after_pulses(x(0), x(1), x(2), i64_arg(args, 3)).into(),
        "steady_state_peak" => d::steady_state_peak(x(0), x(1), x(2)).into(),
        "steady_state_trough" => d::steady_state_trough(x(0), x(1), x(2)).into(),
        "summation_ratio" => d::summation_ratio(x(0), x(1)).into(),
        "summation_ratio_at_hz" => d::summation_ratio_at_hz(x(0), x(1)).into(),
        "pulses_to_fraction" => d::pulses_to_fraction(x(0), x(1), x(2)).into(),
        "spike_shape" => d::spike_shape(x(0), x(1), x(2)).into(),
        "spike_peak_time" => d::spike_peak_time(x(0), x(1)).into(),
        _ => Outcome::unknown(name, args),
    }
}

#[test]
fn golden_vectors() {
    golden_vectors::check_dir(concat!(env!("CARGO_MANIFEST_DIR"), "/vectors"), dispatch);
}
