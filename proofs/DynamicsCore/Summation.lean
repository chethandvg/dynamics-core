/-
SPDX-License-Identifier: Apache-2.0
SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama
-/
import DynamicsCore.Basic
import Mathlib.Analysis.SpecialFunctions.Exp

/-!
# How much repeated pulses add up

The summation ratio `1 / (1 − e^(−T/τ))` is the steady peak in units of one
pulse. It falls towards 1 as the period `T` grows, and for the human egg's
calcium spikes, read by a 75 s integrator, it is already below 1.0004.
-/

namespace DynamicsCore

open Filter Topology

/-- How much repeated pulses add up: `1 / (1 − e^(−T/τ))`. -/
noncomputable def summationRatio (T τ : ℝ) : ℝ := 1 / (1 - Real.exp (-T / τ))

/-- The summation ratio is the steady peak of unit pulses. -/
theorem summationRatio_eq (T τ : ℝ) :
    summationRatio T τ = steadyStatePeak 1 (decayFactor T τ) :=
  rfl

/-- For `T, τ > 0`, repeated pulses always add up a little: the ratio exceeds 1. -/
theorem one_lt_summationRatio {T τ : ℝ} (hT : 0 < T) (hτ : 0 < τ) : 1 < summationRatio T τ := by
  have hneg : -T / τ < 0 := div_neg_of_neg_of_pos (neg_neg_of_pos hT) hτ
  exact one_lt_one_div (by linarith [Real.exp_lt_one_iff.mpr hneg])
    (by linarith [Real.exp_pos (-T / τ)])

/-- `dynamics:SUMMATION-MONOTONE`: for `τ > 0`, the summation ratio strictly
decreases as the period `T > 0` grows, and tends to 1 as `T → ∞`: pulses far
apart stop adding up. -/
theorem summation_monotone {τ : ℝ} (hτ : 0 < τ) :
    StrictAntiOn (fun T => summationRatio T τ) (Set.Ioi 0) ∧
      Tendsto (fun T => summationRatio T τ) atTop (𝓝 1) := by
  constructor
  · intro T₁ h₁ T₂ _ hlt
    rw [Set.mem_Ioi] at h₁
    have e₂ : Real.exp (-T₂ / τ) < Real.exp (-T₁ / τ) :=
      Real.exp_lt_exp.mpr (div_lt_div_of_pos_right (neg_lt_neg hlt) hτ)
    have e₁ : Real.exp (-T₁ / τ) < 1 :=
      Real.exp_lt_one_iff.mpr (div_neg_of_neg_of_pos (neg_neg_of_pos h₁) hτ)
    exact one_div_lt_one_div_of_lt (by linarith) (by linarith)
  · have h0 : Tendsto (fun T => Real.exp (-T / τ)) atTop (𝓝 0) :=
      Real.tendsto_exp_atBot.comp (tendsto_neg_atTop_atBot.atBot_div_const hτ)
    have hlim : Tendsto (fun T => 1 / (1 - Real.exp (-T / τ))) atTop (𝓝 (1 / (1 - 0))) :=
      tendsto_const_nhds.div (tendsto_const_nhds.sub h0) (by norm_num)
    simpa [summationRatio] using hlim

/-- `e⁸ > 2501`, from the first twelve terms of its Taylor series. -/
theorem exp_eight_gt : (2501 : ℝ) < Real.exp 8 :=
  lt_of_lt_of_le (by norm_num [Finset.sum_range_succ, Nat.factorial])
    (Real.sum_le_exp_of_nonneg (by norm_num) 12)

/-- `dynamics:HUMAN-SPIKES-BARELY-ADD-UP`: with `τ = 75` s, spikes 10 minutes
(600 s) apart leave `r = e^(−8)`, and for any period of at least 600 s the
summation ratio is between 1 and 1.0004: what is left of all earlier spikes
adds at most 0.04 % to each peak. -/
theorem human_spikes_barely_add_up {T : ℝ} (hT : 600 ≤ T) :
    decayFactor 600 75 = Real.exp (-8) ∧ 1 < summationRatio T 75 ∧
      summationRatio T 75 < 1.0004 := by
  refine ⟨by norm_num [decayFactor], one_lt_summationRatio (by linarith) (by norm_num), ?_⟩
  refine lt_of_le_of_lt
    ((summation_monotone (by norm_num)).1.antitoneOn (Set.mem_Ioi.mpr (by norm_num))
      (Set.mem_Ioi.mpr (by linarith)) hT) ?_
  have hinv : (Real.exp 8)⁻¹ < 1 / 2501 := by
    rw [inv_eq_one_div]
    exact one_div_lt_one_div_of_lt (by norm_num) exp_eight_gt
  show 1 / (1 - Real.exp (-600 / 75)) < 1.0004
  rw [show (-600 : ℝ) / 75 = -8 by norm_num, Real.exp_neg, div_lt_iff₀ (by linarith),
    show (1.0004 : ℝ) = 2501 / 2500 by norm_num]
  linarith

end DynamicsCore
