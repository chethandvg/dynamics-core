/-
SPDX-License-Identifier: Apache-2.0
SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama
-/
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecificLimits.Basic

/-!
# Leaky integrators and periodic pulses

Definitions mirror the Rust functions of the same names (`decay_factor` is
`decayFactor`, and so on). A pulse train is described by its size `h` and the
residual factor `r` that one period leaves, which for a period `T` and time
constant `τ` is `r = decayFactor T τ = e^(−T/τ)`.

The theorems are about exact real numbers. The Rust code computes the same
formulas in floating point, and the golden vectors check its rounding.
-/

namespace DynamicsCore

open Filter Topology

/-- The fraction of the state left after `dt`, with time constant `τ`: `e^(−dt/τ)`. -/
noncomputable def decayFactor (dt τ : ℝ) : ℝ := Real.exp (-dt / τ)

/-- The state `a` after relaxing for `dt`. -/
noncomputable def relax (a dt τ : ℝ) : ℝ := a * decayFactor dt τ

/-- The time in which the state halves. -/
noncomputable def halfLife (τ : ℝ) : ℝ := τ * Real.log 2

/-- The peak right after the `n`-th pulse, starting from 0: relax for one period
(factor `r`), then add `h`. -/
def peaks (h r : ℝ) : ℕ → ℝ
  | 0 => 0
  | n + 1 => peaks h r n * r + h

/-- The peak that periodic pulses approach. -/
noncomputable def steadyStatePeak (h r : ℝ) : ℝ := h / (1 - r)

/-- The steady state just before each pulse. -/
noncomputable def steadyStateTrough (h r : ℝ) : ℝ := h * r / (1 - r)

/-- `dynamics:RELAX-COMPOSE`: relaxing for `dt₁` and then for `dt₂` is relaxing
for `dt₁ + dt₂`. -/
theorem relax_compose (a dt₁ dt₂ τ : ℝ) :
    relax (relax a dt₁ τ) dt₂ τ = relax a (dt₁ + dt₂) τ := by
  unfold relax decayFactor
  rw [mul_assoc, ← Real.exp_add]
  congr 2
  ring

/-- `dynamics:HALF-LIFE`: after one half-life `τ · ln 2`, half the state is left. -/
theorem decayFactor_halfLife {τ : ℝ} (hτ : τ ≠ 0) : decayFactor (halfLife τ) τ = 1 / 2 := by
  unfold decayFactor halfLife
  rw [neg_div, mul_div_cancel_left₀ _ hτ, Real.exp_neg, Real.exp_log two_pos]
  norm_num

/-- The first peak is one pulse. -/
theorem peaks_one (h r : ℝ) : peaks h r 1 = h := by
  simp [peaks]

/-- `dynamics:PEAK-AFTER-N`: the recurrence `A(n+1) = A(n)·r + h`, `A(1) = h`, has
the closed form `h·(1 − rⁿ)/(1 − r)` for every `r ≠ 1`. -/
theorem peaks_closed_form (h : ℝ) {r : ℝ} (hr : r ≠ 1) (n : ℕ) :
    peaks h r n = h * (1 - r ^ n) / (1 - r) := by
  have h1 : 1 - r ≠ 0 := sub_ne_zero.mpr hr.symm
  induction n with
  | zero => simp [peaks]
  | succ n ih =>
    rw [peaks, ih, eq_div_iff h1]
    field_simp
    ring

/-- The closed form in `T` and `τ` that the Rust `peak_after_pulses` computes:
`h·(1 − e^(−nT/τ))/(1 − e^(−T/τ))`. -/
noncomputable def peakAfterPulses (h T τ : ℝ) (n : ℕ) : ℝ :=
  h * (1 - Real.exp (-(n * T) / τ)) / (1 - Real.exp (-T / τ))

/-- `peakAfterPulses` is the `n`-th peak for the residual factor `r = e^(−T/τ)`,
because `e^(−nT/τ) = rⁿ`. -/
theorem peakAfterPulses_eq_peaks (h : ℝ) {T τ : ℝ} (hT : T ≠ 0) (hτ : τ ≠ 0) (n : ℕ) :
    peakAfterPulses h T τ n = peaks h (decayFactor T τ) n := by
  have hr : decayFactor T τ ≠ 1 := by
    rw [decayFactor, Ne, Real.exp_eq_one_iff]
    exact div_ne_zero (neg_ne_zero.mpr hT) hτ
  rw [peaks_closed_form h hr, peakAfterPulses, decayFactor, ← Real.exp_nat_mul]
  congr 4
  ring

/-- `dynamics:STEADY-STATE`: for `0 ≤ r < 1` the peaks tend to `h / (1 − r)`. -/
theorem peaks_tendsto (h : ℝ) {r : ℝ} (h0 : 0 ≤ r) (h1 : r < 1) :
    Tendsto (peaks h r) atTop (𝓝 (steadyStatePeak h r)) := by
  rw [funext (peaks_closed_form h h1.ne), steadyStatePeak]
  have hlim : Tendsto (fun n : ℕ => h * (1 - r ^ n) / (1 - r)) atTop (𝓝 (h * (1 - 0) / (1 - r))) :=
    ((tendsto_const_nhds.sub (tendsto_pow_atTop_nhds_zero_of_lt_one h0 h1)).const_mul h).div_const _
  simpa using hlim

/-- `dynamics:STEADY-TROUGH`: the steady peak is the fixed point of "relax for one
period, then add `h`", and relaxing it for one period gives the trough `h·r/(1 − r)`. -/
theorem steady_state_fixed_point (h : ℝ) {r : ℝ} (hr : r ≠ 1) :
    steadyStatePeak h r * r + h = steadyStatePeak h r ∧
      steadyStatePeak h r * r = steadyStateTrough h r := by
  have h1 : 1 - r ≠ 0 := sub_ne_zero.mpr hr.symm
  unfold steadyStatePeak steadyStateTrough
  refine ⟨?_, div_mul_eq_mul_div h (1 - r) r⟩
  field_simp
  ring

/-- Each pulse raises the next peak by `h·rⁿ`. -/
theorem peaks_succ_sub (h r : ℝ) (n : ℕ) : peaks h r (n + 1) - peaks h r n = h * r ^ n := by
  induction n with
  | zero => simp [peaks]
  | succ n ih =>
    have step : peaks h r (n + 2) - peaks h r (n + 1) =
        (peaks h r (n + 1) - peaks h r n) * r := by
      simp only [peaks]
      ring
    rw [step, ih, pow_succ]
    ring

/-- `dynamics:PEAKS-INCREASE`: for `h > 0` and `0 < r < 1` the peaks strictly
increase and stay below the steady peak `h / (1 − r)`. -/
theorem peaks_increase {h r : ℝ} (hh : 0 < h) (h0 : 0 < r) (h1 : r < 1) :
    StrictMono (peaks h r) ∧ ∀ n, peaks h r n < steadyStatePeak h r := by
  refine ⟨strictMono_nat_of_lt_succ fun n => ?_, fun n => ?_⟩
  · have := peaks_succ_sub h r n
    have := mul_pos hh (pow_pos h0 n)
    linarith
  · rw [peaks_closed_form h h1.ne, steadyStatePeak, div_lt_div_iff_of_pos_right (sub_pos.mpr h1)]
    nlinarith [mul_pos hh (pow_pos h0 n)]

/-- `dynamics:FRACTION-REACHED`: for `h > 0` and `r < 1`, the `n`-th peak reaches
the fraction `p` of the steady peak exactly when `rⁿ ≤ 1 − p`. This is the
condition that `pulses_to_fraction` decides. -/
theorem fraction_reached {h r : ℝ} (p : ℝ) (hh : 0 < h) (h1 : r < 1) (n : ℕ) :
    p * steadyStatePeak h r ≤ peaks h r n ↔ r ^ n ≤ 1 - p := by
  have hs : 0 < steadyStatePeak h r := div_pos hh (sub_pos.mpr h1)
  have : peaks h r n = (1 - r ^ n) * steadyStatePeak h r := by
    rw [peaks_closed_form h h1.ne, steadyStatePeak]
    ring
  rw [this, mul_le_mul_iff_of_pos_right hs]
  constructor <;> intro <;> linarith

end DynamicsCore
