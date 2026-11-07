# Dynamics core

**Leaky integrators and periodic pulses, with Lean 4 proofs of when repeated pulses add up and when they are forgotten.**

[![CI](https://github.com/chethandvg/dynamics-core/actions/workflows/ci.yml/badge.svg)](https://github.com/chethandvg/dynamics-core/actions/workflows/ci.yml)
[![Proofs](https://github.com/chethandvg/dynamics-core/actions/workflows/proofs.yml/badge.svg)](https://github.com/chethandvg/dynamics-core/actions/workflows/proofs.yml)
[![Licence](https://img.shields.io/badge/licence-Apache--2.0-blue)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24038470.svg)](https://doi.org/10.5281/zenodo.24038470)

Part of **[Frequencies of Life](https://chethandvg.github.io/frequencies-of-life/)**, open code for the rhythms of life in which every claim says how it is known: proved, verified, tested, referenced or assumed.

## Claims and evidence

<!-- claims-table:start -->
| Claim | Level | Statement | Evidence |
|---|---|---|---|
| `dynamics:RELAX-COMPOSE` | P | Relaxing a leaky integrator for dt₁ and then for dt₂ is the same as relaxing it for dt₁ + dt₂: a·e^(−dt₁/τ)·e^(−dt₂/τ) = a·e^(−(dt₁+dt₂)/τ). | [`DynamicsCore.relax_compose`](proofs/DynamicsCore/Basic.lean) |
| `dynamics:HALF-LIFE` | P | After one half-life, τ·ln 2, exactly half of the state is left. | [`DynamicsCore.decayFactor_halfLife`](proofs/DynamicsCore/Basic.lean) |
| `dynamics:PEAK-AFTER-N` | P | The peaks of pulses of size h with residual factor r per period, A(n+1) = A(n)·r + h with A(1) = h, have the closed form A(n) = h·(1 − rⁿ)/(1 − r) for every r ≠ 1. | [`DynamicsCore.peaks_closed_form`](proofs/DynamicsCore/Basic.lean) |
| `dynamics:STEADY-STATE` | P | For 0 ≤ r < 1 the peaks tend to the steady peak h/(1 − r) as the number of pulses grows. | [`DynamicsCore.peaks_tendsto`](proofs/DynamicsCore/Basic.lean) |
| `dynamics:STEADY-TROUGH` | P | The steady peak h/(1 − r) is the fixed point of 'relax for one period, then add h', and relaxing it for one period gives the steady trough h·r/(1 − r). | [`DynamicsCore.steady_state_fixed_point`](proofs/DynamicsCore/Basic.lean) |
| `dynamics:PEAKS-INCREASE` | P | For h > 0 and 0 < r < 1 the peaks strictly increase with every pulse and always stay below the steady peak h/(1 − r). | [`DynamicsCore.peaks_increase`](proofs/DynamicsCore/Basic.lean) |
| `dynamics:FRACTION-REACHED` | P | For h > 0 and r < 1, the n-th peak reaches the fraction p of the steady peak exactly when rⁿ ≤ 1 − p. | [`DynamicsCore.fraction_reached`](proofs/DynamicsCore/Basic.lean) |
| `dynamics:SUMMATION-MONOTONE` | P | For τ > 0 the summation ratio 1/(1 − e^(−T/τ)) strictly decreases as the period T > 0 grows and tends to 1 as T → ∞: pulses far apart stop adding up. | [`DynamicsCore.summation_monotone`](proofs/DynamicsCore/Summation.lean) |
| `dynamics:HUMAN-SPIKES-BARELY-ADD-UP` | P | With τ = 75 s, spikes 10 minutes (600 s) apart leave r = e^(−8), and for every period of at least 600 s the summation ratio is between 1 and 1.0004: what is left of all earlier spikes adds at most 0.04 % to each peak. | [`DynamicsCore.human_spikes_barely_add_up`](proofs/DynamicsCore/Summation.lean) |
| `dynamics:SPIKE-PEAK` | P | For rise, decay > 0 the spike shape (1 − e^(−s/rise))·e^(−s/decay) is largest at s* = rise·ln(1 + decay/rise) > 0, over every real s. | [`DynamicsCore.spike_peak`](proofs/DynamicsCore/Spike.lean) |
| `dynamics:PULSES-TO-FRACTION-EXACT` | T | pulses_to_fraction(p, T, τ) returns the smallest n ≥ 1 with rⁿ ≤ 1 − p, deciding each candidate in double-double arithmetic. It matches a 400-digit reference on every case, including p one double either side of 1 − rⁿ, a subnormal period and an exact tie, where the floating-point formula ⌈ln(1 − p)/ln r⌉ can be off by one. | [`vectors/thresholds.json`](vectors/thresholds.json), [`src/exact.rs`](src/exact.rs), [`tests/properties.rs`](tests/properties.rs) |
| `dynamics:SUMMATION-AT-HZ` | T | summation_ratio_at_hz(f, τ) is summation_ratio(1/f, τ), with the period computed by rhythm-units' period_from_hz. | [`vectors/pulses.json`](vectors/pulses.json), [`tests/properties.rs`](tests/properties.rs) |
| `dynamics:IMPLEMENTATIONS-AGREE` | T | The Rust crate, the Python bindings and the WebAssembly bindings all match golden vectors computed by an independent 400-digit decimal reference: within 1e-13 relative, and pulse counts exactly. | [`vectors/relaxation.json`](vectors/relaxation.json), [`vectors/pulses.json`](vectors/pulses.json), [`vectors/thresholds.json`](vectors/thresholds.json), [`vectors/spikes.json`](vectors/spikes.json), [`tests/vectors.rs`](tests/vectors.rs), [`bindings/python/tests/test_vectors.py`](bindings/python/tests/test_vectors.py), [`bindings/wasm/tests/vectors.test.mjs`](bindings/wasm/tests/vectors.test.mjs) |
| `dynamics:INPUT-VALIDATION` | T | Every function rejects time constants, periods, pulse sizes and frequencies that are not finite and positive, elapsed times and states that are negative or not finite, p outside (0, 1), pulse counts that are not whole numbers of at least 1 (and, in JavaScript, counts above 2^53, which are never wrapped or truncated) and results that overflow or exceed 2^53 pulses; Python raises ValueError (TypeError for a count that is not an int, OverflowError for an int count outside the 64-bit range) and JavaScript throws an Error. | [`vectors/errors.json`](vectors/errors.json), [`tests/properties.rs`](tests/properties.rs), [`src/lib.rs`](src/lib.rs), [`bindings/python/tests/test_vectors.py`](bindings/python/tests/test_vectors.py), [`bindings/wasm/tests/vectors.test.mjs`](bindings/wasm/tests/vectors.test.mjs) |
| `dynamics:TAU-75S` | M | The project's egg-activation model reads calcium spikes with a CaMKII-like leaky integrator whose time constant is 75 s. | assumption: CaMKII activity in the egg relaxes like a linear leaky integrator with τ = 75 s (half-life 52 s), chosen so that it returns to baseline within about 5 minutes of a calcium rise, as reported for CaMKII in mouse eggs; e^(−300/75) leaves 1.8 % after 5 minutes. The 75 s is a chosen model parameter, not a fitted or measured time constant, and using it for human eggs assumes the mouse behaviour carries over. |
<!-- claims-table:end -->

Levels: **P** machine-checked proof (Lean 4) · **V** verified numerics · **T** tested · **R** referenced (DOI or PMID) · **M** model assumption.

## Use

The model is a linear leaky integrator: the state decays as e^(−Δt/τ) between pulses and each pulse adds h. With a pulse every T, one period leaves r = e^(−T/τ), the n-th peak is h·(1 − rⁿ)/(1 − r), and the summation ratio 1/(1 − r) says how much pulses add up (1 means not at all). Times can be in any one unit. Time constants, periods and pulse sizes must be finite and greater than 0; anything else is an error.

**Rust**

```toml
[dependencies]
dynamics-core = { git = "https://github.com/chethandvg/dynamics-core", tag = "v0.1.0" }
```

```rust
use dynamics_core::{pulses_to_fraction, spike_peak_time, steady_state_peak, summation_ratio};

// Spikes 10 minutes apart, read by an integrator with τ = 75 s: r = e^-8.
let ratio = summation_ratio(600.0, 75.0)?;          // 1.0003355752008412
// Pulses every τ/2 add up to 2.54 pulses, and pass 90 % of that at the fifth.
let top = steady_state_peak(1.0, 37.5, 75.0)?;      // 2.541494082536798
assert_eq!(pulses_to_fraction(0.9, 37.5, 75.0)?, 5);
// A spike rising with τ = 15 s and decaying with τ = 60 s peaks at 15·ln 5 s.
let t = spike_peak_time(15.0, 60.0)?;               // 24.141568686511505
```

**Python** (builds the Rust core with maturin)

```bash
pip install "git+https://github.com/chethandvg/dynamics-core@v0.1.0#subdirectory=bindings/python"
```

```python
import dynamics_core as dc

dc.summation_ratio(600, 75)  # 1.0003355752008412
# The same, with the period taken from rhythm-units:
dc.summation_ratio_at_hz(1 / 600, 75)  # 1.0003355752008412
dc.peak_after_pulses(1, 37.5, 75, 5)  # 2.332875544269118
dc.pulses_to_fraction(0.9, 37.5, 75)  # 5
dc.pulses_to_fraction(1, 600, 75)
# ValueError: p must be a number strictly between 0 and 1, got 1
```

**JavaScript / WebAssembly** (build steps in [bindings/wasm](bindings/wasm))

```js
import init, { summation_ratio, pulses_to_fraction } from "./pkg/web/dynamics_core_wasm.js";

await init();
summation_ratio(600, 75);             // 1.0003355752008412
pulses_to_fraction(0.9, 37.5, 75);    // 5
```

The functions are `decay_factor`, `relax`, `half_life`, `peak_after_pulses`, `steady_state_peak`, `steady_state_trough`, `summation_ratio`, `summation_ratio_at_hz`, `pulses_to_fraction`, `spike_shape` and `spike_peak_time`, with the same names in all three languages. Errors raise `ValueError` in Python (a pulse count that is not an `int` raises `TypeError`, and an `int` outside the 64-bit range `OverflowError`) and throw `Error` in JavaScript. JavaScript passes the pulse count `n` as a number, which must be a whole number from 1 to 2⁵³; anything else throws rather than being truncated.

## How it is verified

- **Proofs:** the Lean 4 files in [`proofs/`](proofs) prove the ten P-level claims about the definitions `decayFactor`, `relax`, `peaks`, `steadyStatePeak`, `summationRatio`, `spikeShape` and `spikePeakTime`, which mirror the Rust functions: the closed form of the peaks (and, in `peakAfterPulses_eq_peaks`, that the form in T and τ which `peak_after_pulses` computes is the n-th peak with r = e^(−T/τ)), their limit, the steady trough, the summation ratio falling to 1, the human-spike bound (from e⁸ > 2501, using twelve Taylor terms) and the spike peak (from one use of the convexity of exp). CI rebuilds them with Mathlib, rejects `sorry`, and allows only the standard axioms (`propext`, `Classical.choice`, `Quot.sound`).
- **Implementation:** one Rust implementation, compiled to Python (PyO3) and WebAssembly. The same 124 golden vectors in [`vectors/`](vectors) are checked in all three: Rust and Python on Linux and Windows, WebAssembly on Linux. Their oracle, [`vectors/generate.py`](vectors/generate.py), shares no code with the crate: it works in 400-digit decimal arithmetic from the exact input doubles, and computes peaks by iterating the recurrence rather than with the closed form. Property tests in [`tests/properties.rs`](tests/properties.rs) add NaN and infinity, the recurrence over 200 pulses, monotonicity and a dense check of the spike peak.
- **Exact thresholds:** `pulses_to_fraction` checks its first guess against the definition in double-double arithmetic (about 106 bits), so it is exact at the boundaries where ⌈ln(1 − p)/ln r⌉ in floating point is off by one. Its vectors include p one double either side of 1 − rⁿ, a subnormal period and an exact tie, with zero tolerance.
- **What the proofs don't cover:** they are about exact real-number mathematics. The Rust code uses 64-bit floating point and `exp`, `exp_m1` and `ln_1p`, whose rounding is checked by the golden-vector tolerances (1e-13 relative), not by the proofs; the double-double decision in `pulses_to_fraction` is tested, not proved. Input validation and `summation_ratio_at_hz` are tested. The proofs say nothing about biology: τ = 75 s is a model parameter chosen to match mouse data (`dynamics:TAU-75S`, level M), and the 10–35 min spike interval is a reference, not a theorem.

## Cite

If you use this work, please cite it. GitHub's **Cite this repository** button gives APA and BibTeX from [CITATION.cff](CITATION.cff).

## Licence

Code: Apache-2.0 ([LICENSE](LICENSE)). Copyright notices are in [NOTICE](NOTICE).
