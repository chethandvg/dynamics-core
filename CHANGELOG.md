# Changelog

All notable changes are listed here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-11-07

### Added

- Relaxation: `decay_factor`, `relax` and `half_life` for a linear leaky integrator. `relax` (and `steady_state_trough`) never form e^(−Δt/τ) on its own, so a large state stays accurate after more than 708 time constants.
- Periodic pulses: `peak_after_pulses` (closed form with `exp_m1`), `steady_state_peak`, `steady_state_trough`, `summation_ratio` and `summation_ratio_at_hz` (built on `rhythm_units::period_from_hz`).
- `pulses_to_fraction`, exact at the boundaries: a first guess is checked against the definition in double-double arithmetic.
- Spikes: `spike_shape` and `spike_peak_time`.
- `DynamicsError` for invalid input and out-of-range results, and `pulse_count`, which turns a JavaScript number into a pulse count or rejects it.
- Python (PyO3) and WebAssembly bindings with the same function names. In JavaScript the pulse count is a number that must be a whole number from 1 to 2⁵³.
- Lean 4 proofs of ten P-level claims, including the closed form and limit of the peaks, `dynamics:SUMMATION-MONOTONE`, `dynamics:SPIKE-PEAK` and `dynamics:HUMAN-SPIKES-BARELY-ADD-UP`, and `peakAfterPulses_eq_peaks`, which links the closed form in T and τ that `peak_after_pulses` computes to the recurrence.
- Golden vectors (124 cases) from an independent 400-digit decimal oracle, shared by Rust, Python and WebAssembly, and property tests.
- The model assumption `dynamics:TAU-75S` (τ = 75 s), with its source.
