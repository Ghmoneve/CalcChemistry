import Mathlib.Data.Real.Basic
import Mathlib.Data.Finset.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Certified energy intervals for surface reactions

This file isolates the real-arithmetic core used to turn computed energies with
error bars into certified statements about a surface reaction. It deliberately
contains no physical assumption: every physical input enters as an explicit
hypothesis, and the theorems below only propagate intervals.

The central honest limitation is proved as a theorem
(`sampling_gives_no_upper_bound`): sampling energies at a finite set of
configurations can certify a *lower* bound on a barrier but can never certify an
*upper* bound, because the sampled set is not the path. Any claim of an upper
bound requires a connected path as an explicit input.
-/

namespace CalcChemistry

/-! ## Adsorption energy -/

/-- Adsorption (or reaction) energy built from three separately bounded total
energies. This is the arithmetic behind `E_ads = E_complex - E_slab - E_gas`. -/
theorem adsorption_energy_interval
    {Ec EcL EcU Es EsL EsU Em EmL EmU : ℝ}
    (hc : EcL ≤ Ec ∧ Ec ≤ EcU)
    (hs : EsL ≤ Es ∧ Es ≤ EsU)
    (hm : EmL ≤ Em ∧ Em ≤ EmU) :
    EcL - EsU - EmU ≤ Ec - Es - Em ∧ Ec - Es - Em ≤ EcU - EsL - EmL := by
  constructor <;> linarith

/-- A bounded adsorption energy is bounded after adding one independently
bounded model correction. -/
theorem adsorption_with_correction
    {L U E correction radius : ℝ}
    (hE : L ≤ E ∧ E ≤ U) (hc : |correction| ≤ radius) :
    L - radius ≤ E + correction ∧ E + correction ≤ U + radius := by
  rw [abs_le] at hc
  constructor <;> linarith

/-! ## Growth per cycle and surface coverage -/

/-- A measured growth per cycle together with the layer spacing pins the
fractional surface coverage. With `gpc` in Å/cycle and `spacing` the thickness of
one complete layer, `coverage` is the fraction of a layer deposited per cycle. -/
theorem gpc_coverage_interval {gpc g1 g2 spacing coverage : ℝ}
    (hspacing : 0 < spacing) (h1 : g1 ≤ gpc) (h2 : gpc ≤ g2)
    (heq : gpc = spacing * coverage) :
    g1 / spacing ≤ coverage ∧ coverage ≤ g2 / spacing := by
  have hcov : coverage = gpc / spacing := by
    rw [heq]
    field_simp
  rw [hcov]
  exact ⟨div_le_div_of_nonneg_right h1 (le_of_lt hspacing),
         div_le_div_of_nonneg_right h2 (le_of_lt hspacing)⟩

/-- The converse direction: a coverage fraction and layer spacing give the
growth per cycle, with the same interval structure. -/
theorem coverage_gpc_interval {coverage c1 c2 spacing gpc : ℝ}
    (hspacing : 0 < spacing) (h1 : c1 ≤ coverage) (h2 : coverage ≤ c2)
    (heq : gpc = spacing * coverage) :
    spacing * c1 ≤ gpc ∧ gpc ≤ spacing * c2 := by
  rw [heq]
  exact ⟨mul_le_mul_of_nonneg_left h1 (le_of_lt hspacing),
         mul_le_mul_of_nonneg_left h2 (le_of_lt hspacing)⟩

/-! ## Barriers from path sampling -/

/-- At a single configuration on the path, a certified lower bound on its energy
minus an upper bound on the reactant energy is a certified lower bound on that
configuration's relative energy. Because the barrier is the maximum of the
relative energy along the path, this is a valid lower bound on the barrier
provided the configuration lies on the path. That proviso is a physical input,
not a consequence of the arithmetic. -/
theorem barrier_lower_from_samples
    {ι : Type*} (L E : ι → ℝ) (i₀ : ι) (Er Ur : ℝ)
    (hE : L i₀ ≤ E i₀) (hr : Er ≤ Ur) :
    L i₀ - Ur ≤ E i₀ - Er := by
  linarith

/-- The corresponding upper bound at a single configuration, which requires the
same path proviso plus the configuration being the actual maximum. -/
theorem barrier_upper_at_max
    {ι : Type*} (U E : ι → ℝ) (i₀ : ι) (Er Lr : ℝ)
    (hE : E i₀ ≤ U i₀) (hr : Lr ≤ Er) :
    E i₀ - Er ≤ U i₀ - Lr := by
  linarith

/-- Combining the two directions: if a configuration is known to be the maximum
of the relative energy along a connected path, its energy interval and the
reactant energy interval give a two-sided certified barrier interval. -/
theorem certified_barrier_interval
    {E L U : ℝ → ℝ} {i₀ : ℝ} {Er Lr Ur : ℝ}
    (hlo : L i₀ ≤ E i₀) (hhi : E i₀ ≤ U i₀) (hr : Lr ≤ Er) (hr' : Er ≤ Ur) :
    L i₀ - Ur ≤ E i₀ - Er ∧ E i₀ - Er ≤ U i₀ - Lr := by
  constructor <;> linarith

/-- **The honest limitation.** Sampling a finite set of configurations cannot
certify an upper bound on a barrier. Here the sampled set has all energies zero
while an unsampled configuration carries energy 10, so the sampled maximum is 0
and the true maximum is 10.

This is not a defect of the method: a lower bound is what sampling can give, and
any upper bound must come from a connected path supplied as an explicit input. -/
theorem sampling_gives_no_upper_bound :
    ∃ (E : Fin 3 → ℝ) (S : Finset (Fin 3)),
      S.Nonempty ∧ (∀ i ∈ S, E i = 0) ∧ (∃ i : Fin 3, E i = 10) := by
  refine ⟨fun i => if i = 1 then 10 else 0, ({0, 2} : Finset (Fin 3)), ?_, ?_, ⟨1, by simp⟩⟩
  · simp
  · intro i hi
    simp only [Finset.mem_insert, Finset.mem_singleton] at hi
    rcases hi with h | h <;> simp [h]

/-- The lower-bound direction does hold for every sampled configuration
simultaneously, so the certified lower bound is the maximum over samples. -/
theorem lower_bound_holds_at_every_sample
    {ι : Type*} (L E : ι → ℝ) (Er Ur : ℝ) (i : ι)
    (hE : L i ≤ E i) (hr : Er ≤ Ur) :
    L i - Ur ≤ E i - Er := by
  linarith

/-! ## Error budget -/

/-- Repeated independent corrections widen an interval additively. Two
corrections: the radius is the sum of radii. -/
theorem two_corrections_widen
    {L U E c1 r1 c2 r2 : ℝ}
    (hE : L ≤ E ∧ E ≤ U)
    (h1 : |c1| ≤ r1) (h2 : |c2| ≤ r2) :
    L - (r1 + r2) ≤ E + c1 + c2 ∧ E + c1 + c2 ≤ U + (r1 + r2) := by
  rw [abs_le] at h1 h2
  constructor <;> linarith

/-- A positive quantity bracketed away from zero stays positive under an
interval perturbation smaller than its lower bound. Used to certify that a
computed gap or barrier cannot vanish within the error bars. -/
theorem stays_positive {L x : ℝ} (h : L ≤ x) (hpos : 0 < L) : 0 < x :=
  lt_of_lt_of_le hpos h

/-- Disjoint intervals certify a sign: if the upper end of a quantity's interval
is below zero, the quantity is negative regardless of its exact value. -/
theorem certified_negative {U x : ℝ} (h : x ≤ U) (hU : U < 0) : x < 0 :=
  lt_of_le_of_lt h hU

/-! ## Steric footprint and a precursor design rule -/

/-- Non-overlap of adsorbates: if `sigma` precursors deposit per unit area and
each occupies at least area `A`, the occupied fractions cannot exceed one. -/
theorem packing_bound {sigma A : ℝ} (hsigma : 0 < sigma) (hpack : sigma * A ≤ 1) :
    A ≤ 1 / sigma := by
  rw [le_div_iff₀ hsigma]
  linarith

/-- **Precursor design rule.** A ligand of footprint at least `A0` caps the
deposited areal density at `1 / A0`, no matter how fast the surface chemistry is.
This is a geometric statement, not a kinetic one: it bounds what a cycle *can*
deposit, and says nothing about whether the cycle reaches that bound. -/
theorem max_density_from_ligand_size {sigma A A0 : ℝ}
    (hA0 : 0 < A0) (hlig : A0 ≤ A) (hpack : sigma * A ≤ 1) :
    sigma ≤ 1 / A0 := by
  have h1 : sigma ≤ 1 / A := by
    rw [le_div_iff₀ (lt_of_lt_of_le hA0 hlig)]
    linarith
  have h2 : 1 / A ≤ 1 / A0 := by
    apply div_le_div_of_nonneg_left (by linarith) hA0 hlig
  linarith

/-- The same design rule expressed as a bound on growth per cycle. With `sigma`
the areal density of deposited metal, `Asite` the area per surface site, and
`spacing` the thickness of one complete layer, `gpc = spacing * sigma * Asite`. -/
theorem gpc_upper_from_ligand_footprint
    {gpc spacing sigma A A0 Asite : ℝ}
    (hspacing : 0 < spacing) (hA0 : 0 < A0) (hAsite : 0 < Asite)
    (hlig : A0 ≤ A) (hpack : sigma * A ≤ 1)
    (hdef : gpc = spacing * sigma * Asite) :
    gpc ≤ spacing * Asite / A0 := by
  have hs : sigma ≤ 1 / A0 := max_density_from_ligand_size hA0 hlig hpack
  rw [hdef]
  have hstep : spacing * sigma * Asite ≤ spacing * (1 / A0) * Asite := by
    have hpos : 0 < spacing * Asite := mul_pos hspacing hAsite
    nlinarith [hs, hpos, hspacing, hAsite]
  calc spacing * sigma * Asite ≤ spacing * (1 / A0) * Asite := hstep
    _ = spacing * Asite / A0 := by field_simp

/-- **Falsification test.** The measured areal density `sigma` and an
independently measured footprint `A` are jointly consistent only if
`sigma * A ≤ 1`. When the product exceeds one, the hypothesis that steric
occupancy by that footprint controls the cycle is refuted: no arrangement of
these molecules can cover the surface at that density without overlap.

This is the formal content of the falsifiable statement in the worked example:
the observed `sigma` and the crystallographic footprint agree, but a future
measurement placing the footprint below `1 / sigma` would exclude steric
self-limitation by the intact dimer. -/
theorem steric_falsification {sigma A : ℝ} (hsigma : 0 < sigma) (hA : 0 < A)
    (hincons : 1 < sigma * A) : ¬ (sigma ≤ 1 / A) := by
  intro h
  rw [le_div_iff₀ hA] at h
  linarith

/-- The consistency direction, stated so that a numerical check has a named
target: a footprint `A` is admissible for a measured density `sigma` exactly
when `sigma * A ≤ 1`. -/
theorem steric_admissible {sigma A : ℝ} (hsigma : 0 < sigma) (hA : 0 < A)
    (hpack : sigma * A ≤ 1) : sigma ≤ 1 / A := by
  rw [le_div_iff₀ hA]
  linarith

end CalcChemistry
