/-
SPDX-License-Identifier: Apache-2.0
SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama
-/
import Mathlib.Analysis.Convex.SpecificFunctions.Basic

/-!
# The peak of a calcium spike

A spike that rises with time constant `rise` and decays with time constant
`decay` has the shape `(1 − e^(−s/rise))·e^(−s/decay)`. Its peak is at
`s* = rise·ln(1 + decay/rise)`. The proof needs no calculus: one application of
the convexity of `exp` bounds the shape at every `s` by its value at `s*`.
-/

namespace DynamicsCore

/-- The spike shape `s` after it starts. -/
noncomputable def spikeShape (s rise decay : ℝ) : ℝ :=
  (1 - Real.exp (-s / rise)) * Real.exp (-s / decay)

/-- When the spike peaks. -/
noncomputable def spikePeakTime (rise decay : ℝ) : ℝ := rise * Real.log (1 + decay / rise)

/-- `dynamics:SPIKE-PEAK`: for `rise, decay > 0` the peak time is positive and the
spike shape is largest there: no `s`, before or after, gives a higher value. -/
theorem spike_peak {rise decay : ℝ} (hr : 0 < rise) (hd : 0 < decay) :
    0 < spikePeakTime rise decay ∧
      ∀ s, spikeShape s rise decay ≤ spikeShape (spikePeakTime rise decay) rise decay := by
  have hq : 0 < decay / rise := div_pos hd hr
  refine ⟨mul_pos hr (Real.log_pos (by linarith)), fun s => ?_⟩
  set t := spikePeakTime rise decay with ht
  obtain ⟨δ, rfl⟩ : ∃ δ, s = t + δ := ⟨s - t, by ring⟩
  have hrd : 0 < rise + decay := by linarith
  have hr' := hr.ne'
  have hd' := hd.ne'
  have hrd' := hrd.ne'
  -- At the peak, e^(−t/rise) = rise/(rise + decay).
  have hT : Real.exp (-t / rise) = rise / (rise + decay) := by
    rw [ht, spikePeakTime, neg_div, mul_div_cancel_left₀ _ hr', Real.exp_neg,
      Real.exp_log (by linarith)]
    field_simp
  -- Convexity of exp between 0 and −δ·(rise + decay)/(rise·decay).
  have hconv := convexOn_exp.2 (Set.mem_univ 0)
    (Set.mem_univ (-δ * (rise + decay) / (rise * decay)))
    (div_nonneg hd.le hrd.le) (div_nonneg hr.le hrd.le) (by field_simp; ring)
  simp only [smul_eq_mul, mul_zero, zero_add, Real.exp_zero, mul_one] at hconv
  have e₁ : rise / (rise + decay) * (-δ * (rise + decay) / (rise * decay)) = -δ / decay := by
    field_simp
  have e₂ : -δ * (rise + decay) / (rise * decay) = -δ / rise + -δ / decay := by
    field_simp
    ring
  rw [e₁, e₂, Real.exp_add] at hconv
  unfold spikeShape
  rw [show -(t + δ) / rise = -t / rise + -δ / rise by ring,
    show -(t + δ) / decay = -t / decay + -δ / decay by ring, Real.exp_add, Real.exp_add, hT,
    show 1 - rise / (rise + decay) = decay / (rise + decay) by field_simp; ring]
  nlinarith [mul_le_mul_of_nonneg_left hconv (Real.exp_pos (-t / decay)).le]

end DynamicsCore
