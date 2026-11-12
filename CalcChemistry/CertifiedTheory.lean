import Mathlib.Data.Complex.Exponential
import Mathlib.Data.Real.Sqrt
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

namespace CalcChemistry

/-- Bookkeeping on real, not integer, energy intervals. -/
theorem reaction_interval
    {Lts Uts Lr Ur ets er : ℝ}
    (h_ts : Lts ≤ ets ∧ ets ≤ Uts)
    (h_r : Lr ≤ er ∧ er ≤ Ur) :
    Lts - Ur ≤ ets - er ∧ ets - er ≤ Uts - Lr := by
  constructor <;> linarith

/-- Each independently bounded additive model correction widens the interval. -/
theorem add_error_interval {L U E correction radius : ℝ}
    (hE : L ≤ E ∧ E ≤ U) (hc : |correction| ≤ radius) :
    L - radius ≤ E + correction ∧ E + correction ≤ U + radius := by
  rw [abs_le] at hc
  constructor <;> linarith

/-- A two-variable positive-semidefinite form, including a zero diagonal. -/
theorem psd_quadratic_nonneg (u v b x y : ℝ)
    (hu : 0 ≤ u) (hv : 0 ≤ v) (hdet : b ^ 2 ≤ u * v) :
    0 ≤ u * x ^ 2 - 2 * b * x * y + v * y ^ 2 := by
  by_cases hz : u = 0
  · have hb : b = 0 := by nlinarith [sq_nonneg b]
    subst hb
    subst hz
    have : 0 ≤ v * y ^ 2 := mul_nonneg hv (sq_nonneg y)
    nlinarith
  · have hpos : 0 < u := lt_of_le_of_ne hu (Ne.symm hz)
    have hnum : 0 ≤ (u * x - b * y) ^ 2 + (u * v - b ^ 2) * y ^ 2 :=
      add_nonneg (sq_nonneg _) (mul_nonneg (sub_nonneg.mpr hdet) (sq_nonneg y))
    have heq : u * x ^ 2 - 2 * b * x * y + v * y ^ 2 =
        ((u * x - b * y) ^ 2 + (u * v - b ^ 2) * y ^ 2) / u := by
      field_simp [hz]
      ring
    rw [heq]
    exact div_nonneg hnum (le_of_lt hpos)

/-- Sharp two-block scalar lower bound. The independent matrix proof reduces
norms of block vectors to x,y and then uses this real-algebra result. -/
theorem sharp_block_scalar_lower (a c b x y : ℝ) :
    let L := (a + c - Real.sqrt ((a - c) ^ 2 + 4 * b ^ 2)) / 2
    L * (x ^ 2 + y ^ 2) ≤ a * x ^ 2 - 2 * b * x * y + c * y ^ 2 := by
  let s := Real.sqrt ((a - c) ^ 2 + 4 * b ^ 2)
  have hs0 : 0 ≤ s := Real.sqrt_nonneg _
  have hs2 : s ^ 2 = (a - c) ^ 2 + 4 * b ^ 2 :=
    Real.sq_sqrt (add_nonneg (sq_nonneg _) (mul_nonneg (by norm_num) (sq_nonneg _)))
  have hsa : a - c ≤ s := by
    by_contra h
    have hh : s < a - c := lt_of_not_ge h
    nlinarith [sq_nonneg b]
  have hsc : c - a ≤ s := by
    by_contra h
    have hh : s < c - a := lt_of_not_ge h
    nlinarith [sq_nonneg b]
  have h := psd_quadratic_nonneg (a - (a + c - s) / 2)
    (c - (a + c - s) / 2) b x y (by linarith) (by linarith) (by nlinarith)
  dsimp only
  change (a + c - s) / 2 * (x ^ 2 + y ^ 2) ≤ _
  nlinarith

/-- The scalar completion-of-squares inequality behind the gapped block bound. -/
theorem quadratic_gap_lower (a c b x y : ℝ) (hgap : 0 < c - a) :
    (a - b ^ 2 / (c - a)) * (x ^ 2 + y ^ 2) ≤
      a * x ^ 2 - 2 * b * x * y + c * y ^ 2 := by
  have hn : 0 ≤ (((c - a) * y - b * x) ^ 2 + (b * y) ^ 2) / (c - a) :=
    div_nonneg (add_nonneg (sq_nonneg _) (sq_nonneg _)) (le_of_lt hgap)
  calc
    (a - b ^ 2 / (c - a)) * (x ^ 2 + y ^ 2) ≤
        (a - b ^ 2 / (c - a)) * (x ^ 2 + y ^ 2) +
        (((c - a) * y - b * x) ^ 2 + (b * y) ^ 2) / (c - a) :=
      le_add_of_nonneg_right hn
    _ = a * x ^ 2 - 2 * b * x * y + c * y ^ 2 := by
      field_simp [ne_of_gt hgap]
      ring

/-- Under a block quadratic-form hypothesis, any normalized state obeys the
energy lower bound. Obtaining the hypothesis for molecular integrals remains
an explicit obligation; it is not asserted as an axiom. -/
theorem normalized_block_lower {a c b x y energy : ℝ}
    (hgap : 0 < c - a) (hnorm : x ^ 2 + y ^ 2 = 1)
    (hform : a * x ^ 2 - 2 * b * x * y + c * y ^ 2 ≤ energy) :
    a - b ^ 2 / (c - a) ≤ energy := by
  have h := quadratic_gap_lower a c b x y hgap
  rw [hnorm, mul_one] at h
  exact le_trans h hform

/-- Transfer a uniform Rayleigh lower bound and one trial upper bound to an
attained minimum. The finite-dimensional spectral theorem supplies such a
minimum, but that matrix-to-Rayleigh step is not formalized here. -/
theorem ground_energy_interval {S : Type*} {rayleigh : S → ℝ}
    {e lower upper : ℝ}
    (hmin : (∀ s, e ≤ rayleigh s) ∧ ∃ s, rayleigh s = e)
    (h_lower : ∀ s, lower ≤ rayleigh s)
    (h_trial : ∃ s, rayleigh s ≤ upper) : lower ≤ e ∧ e ≤ upper := by
  obtain ⟨smin, hsmin⟩ := hmin.2
  obtain ⟨strial, hstrial⟩ := h_trial
  constructor
  · simpa [hsmin] using h_lower smin
  · exact le_trans (hmin.1 strial) hstrial

/-- Monotonic propagation through the real exponential; theta is k_B*T for
per-particle energy, or R*T for molar energy. -/
theorem exponential_barrier_interval {theta L B U : ℝ}
    (htheta : 0 < theta) (hL : L ≤ B) (hU : B ≤ U) :
    Real.exp (-U / theta) ≤ Real.exp (-B / theta) ∧
    Real.exp (-B / theta) ≤ Real.exp (-L / theta) := by
  constructor
  · exact Real.exp_le_exp.mpr
      (div_le_div_of_nonneg_right (neg_le_neg hU) (le_of_lt htheta))
  · exact Real.exp_le_exp.mpr
      (div_le_div_of_nonneg_right (neg_le_neg hL) (le_of_lt htheta))

/-- Rate interval allowing uncertainty in both prefactor and barrier.
This proves propagation conditional on the rate law, not TST itself. -/
theorem rate_interval {theta L B U Alo A Ahi : ℝ}
    (htheta : 0 < theta) (hL : L ≤ B) (hU : B ≤ U)
    (hAlo : 0 ≤ Alo) (hA₁ : Alo ≤ A) (hA₂ : A ≤ Ahi) :
    Alo * Real.exp (-U / theta) ≤ A * Real.exp (-B / theta) ∧
    A * Real.exp (-B / theta) ≤ Ahi * Real.exp (-L / theta) := by
  obtain ⟨hl, hu⟩ := exponential_barrier_interval htheta hL hU
  constructor
  · exact mul_le_mul hA₁ hl (le_of_lt (Real.exp_pos _)) (le_trans hAlo hA₁)
  · exact mul_le_mul hA₂ hu (le_of_lt (Real.exp_pos _))
      (le_trans hAlo (le_trans hA₁ hA₂))

/-- Restricting a minimization domain cannot lower its minimum. This is the
order argument behind outer-relaxation lower bounds; feasibility and global
optimality are explicit hypotheses. -/
theorem relaxation_lower_bound {S : Type*} {physical relaxed : Set S}
    {energy : S → ℝ} {lower exactEnergy : ℝ}
    (hincl : physical ⊆ relaxed)
    (hrelax : ∀ s ∈ relaxed, lower ≤ energy s)
    (hground : ∃ s ∈ physical, energy s = exactEnergy) : lower ≤ exactEnergy := by
  obtain ⟨s, hs, hE⟩ := hground
  simpa [hE] using hrelax s (hincl hs)

/-- Disjoint sound intervals certify a strict ranking. -/
theorem certified_ranking {La Ua Lb Ub a b : ℝ}
    (ha : La ≤ a ∧ a ≤ Ua) (hb : Lb ≤ b ∧ b ≤ Ub) (hsep : Ua < Lb) :
    a < b := lt_of_le_of_lt ha.2 (lt_of_lt_of_le hsep hb.1)

/-- Correct characteristic polynomial of the paired block in H_z convention. -/
theorem pairing_characteristic (d g E : ℝ) :
    (d - g - E) * (-d - g - E) - g ^ 2 =
      (E + g) ^ 2 - (d ^ 2 + g ^ 2) := by ring

/-- The lower analytic root satisfies the characteristic equation. This is
not a formalization of the entire fermionic/molecular Hamiltonian. -/
theorem pairing_lower_root (d g : ℝ) :
    let E := -g - Real.sqrt (d ^ 2 + g ^ 2)
    (d - g - E) * (-d - g - E) - g ^ 2 = 0 := by
  dsimp
  rw [pairing_characteristic]
  have hs := Real.sq_sqrt (add_nonneg (sq_nonneg d) (sq_nonneg g))
  nlinarith

end CalcChemistry
