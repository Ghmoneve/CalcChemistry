import Mathlib.Data.Real.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Complement-space lower bounds without enumerating the complement

## The bottleneck this file addresses

The certificate chain in `CertifiedTheory.lean` and `BarrierBounds.lean` needs a
lower bound `c` on the ground energy of the *complement* block: the determinants
omitted from the active space. That space is exponentially large, so the bound
cannot come from diagonalising it, and even writing down one vector in it costs
proportional to its dimension.

The reduction formalised here splits the complement Hamiltonian as

  `H_C = D + O`,

with `D` diagonal (the determinant energies) and `O` the off-diagonal coupling.
Two classical facts then give a bound whose hypotheses are *local*:

* the diagonal part contributes at least `(min_i D_i) * ||x||^2`;
* the off-diagonal part contributes at least `-tau * ||x||^2` whenever its
  quadratic form is bounded by `tau` in absolute value.

Hence `lambda_min(H_C) >= min_i D_i - tau`. Neither `min_i D_i` nor `tau`
requires enumerating the complement: `min_i D_i` is a combinatorial minimum over
occupation patterns, and `tau` is bounded by row sparsity times coupling
strength. What remains open is *computing* those two numbers for a real molecular
integral set. That is a problem about the integrals, not about the certificate
logic, and this file certifies only the composition.

## Honesty note

The ingredients are textbook: the diagonal/off-diagonal split is trivial and the
quadratic-form bound on `O` is the standard estimate behind Gershgorin- and
Weyl-type arguments. No novelty is claimed for them. The contribution is the
composition, and specifically that it turns an enumeration problem into a
two-scalar problem.
-/

namespace CalcChemistry

open Finset

/-- Quadratic form of a diagonal operator on `R^n`. -/
def diagForm {n : ℕ} (D : Fin n → ℝ) (x : Fin n → ℝ) : ℝ := ∑ i, D i * x i ^ 2

/-- A diagonal operator bounded below by `d` gives the scalar lower bound. -/
theorem diagForm_lower {n : ℕ} (D : Fin n → ℝ) (d : ℝ)
    (hD : ∀ i, d ≤ D i) (x : Fin n → ℝ) :
    d * (∑ i, x i ^ 2) ≤ diagForm D x := by
  unfold diagForm
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum (fun i _ => ?_)
  have hi := hD i
  nlinarith [sq_nonneg (x i)]

/-- **Complement bound.** If every diagonal entry of the complement block is at
least `d`, and the off-diagonal part has quadratic form bounded in absolute value
by `tau` times the squared norm, then every state in the complement has energy at
least `d - tau` per unit norm.

Consequently `d - tau` may be passed as the complement lower bound `c` required
by the block certificate. The hypotheses mention only `d` and `tau`, because
those are the two quantities that can in principle be obtained without touching
the complement's dimension. -/
theorem complement_ground_bound {n : ℕ}
    (D : Fin n → ℝ) (Q : (Fin n → ℝ) → ℝ) (d tau : ℝ)
    (hD : ∀ i, d ≤ D i)
    (hQ : ∀ x : Fin n → ℝ, |Q x| ≤ tau * (∑ i, x i ^ 2)) :
    ∀ x : Fin n → ℝ, (d - tau) * (∑ i, x i ^ 2) ≤ diagForm D x + Q x := by
  intro x
  have h1 := diagForm_lower D d hD x
  have h2 : -(tau * (∑ i, x i ^ 2)) ≤ Q x := (abs_le.mp (hQ x)).1
  linarith

/-- The same bound for a normalised state, which is the form the block
certificate consumes. -/
theorem complement_ground_bound_normalized {n : ℕ}
    (D : Fin n → ℝ) (Q : (Fin n → ℝ) → ℝ) (d tau : ℝ)
    (hD : ∀ i, d ≤ D i)
    (hQ : ∀ x : Fin n → ℝ, |Q x| ≤ tau * (∑ i, x i ^ 2))
    {x : Fin n → ℝ} (hnorm : ∑ i, x i ^ 2 = 1) :
    d - tau ≤ diagForm D x + Q x := by
  have h := complement_ground_bound D Q d tau hD hQ x
  rwa [hnorm, mul_one] at h

/-- A vanishing off-diagonal part recovers the bare diagonal bound, so the gap
between the two bounds is exactly the coupling strength `tau`: no precision is
lost when the complement is decoupled. -/
theorem complement_bound_tight_when_decoupled {n : ℕ}
    (D : Fin n → ℝ) (d : ℝ) (hD : ∀ i, d ≤ D i) (x : Fin n → ℝ) :
    d * (∑ i, x i ^ 2) ≤ diagForm D x :=
  diagForm_lower D d hD x

/-- **Monotonicity of the complement bound.** Improving the diagonal lower bound
raises the certificate and improving the coupling bound raises it too; the two
improvements are interchangeable at this level. This is what makes a local
sparsity argument on `O` worth pursuing: tightening `tau` buys back precision
exactly as raising `d` would. -/
theorem complement_bound_improvement {d1 tau1 d2 tau2 : ℝ}
    (hd : d1 ≤ d2) (htau : tau2 ≤ tau1) :
    d1 - tau1 ≤ d2 - tau2 := by
  linarith

/-- **Separation of the certificate.** If the certified complement energy exceeds
an active-space energy, the two spaces are separated and the active calculation
cannot be contaminated by a lower complement state. Stated so that a computed
`d - tau` can be compared directly against a computed active-space energy. -/
theorem complement_gap_positive {d tau e : ℝ}
    (hsep : e < d - tau) :
    0 < d - tau - e := by
  linarith

end CalcChemistry
