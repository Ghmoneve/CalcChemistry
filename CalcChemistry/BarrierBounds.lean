import Mathlib.Data.Real.Sqrt
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import CalcChemistry.CertifiedTheory

/-!
# Certified reaction barriers from finite-dimensional block data

Setting. A real finite-dimensional Hilbert space `R^n (+) R^m` (blocks A and C,
the complement). `H` enters only through three pieces:

* `A : Fin n -> Fin n -> R`,  the A-block of `H`;
* `C : Fin m -> Fin m -> R`,  the complement block, restricted to the
  "transition region";
* `B : Fin n -> Fin m -> R`,  the coupling.

The quadratic form of `H` on a *pair* `(u, v)` is `blockForm A B C u v`. The
block data available in practice are

* lower bounds `a`, `c` on the two diagonal block quadratic forms
  (`ha`, `hc`), i.e. Rayleigh-quotient lower bounds on the block ground
  energies `lambda_A`, `lambda_C`;
* an upper bound `b` on the coupling quadratic form (`hb`).

This file proves what that data does and does not certify. Nothing here is
specialised to a molecule: every physical input is an explicit hypothesis.
-/

namespace CalcChemistry

open Finset

/-! ## 1. The block quadratic form -/

/-- Quadratic form of a self-adjoint block operator with blocks `A`, `C` and
coupling `B`, evaluated at a state `(u, v)`. Hermiticity makes the two
off-diagonal cross terms equal, which is why a single factor `2` appears. -/
def blockForm {n m : ℕ} (A : Fin n → Fin n → ℝ) (B : Fin n → Fin m → ℝ)
    (C : Fin m → Fin m → ℝ) (u : Fin n → ℝ) (v : Fin m → ℝ) : ℝ :=
  (∑ i, ∑ j, u i * A i j * u j) + 2 * (∑ i, ∑ j, u i * B i j * v j)
    + (∑ i, ∑ j, v i * C i j * v j)

/-- The squared norm of the pair `(u,v)` in the orthogonal direct sum. -/
def blockNormSq {n m : ℕ} (u : Fin n → ℝ) (v : Fin m → ℝ) : ℝ :=
  (∑ i, u i ^ 2) + ∑ j, v j ^ 2

/-- **Block certificate (sharp form).** If `a`, `c` are Rayleigh-quotient lower
bounds on the two diagonal blocks and `b` bounds the coupling form, then every
state of the full space satisfies

`((a + c - sqrt((a-c)^2 + 4 b^2))/2) * ||(u,v)||^2 <= <(u,v), H (u,v)>`.

The constant is the lower eigenvalue of the 2x2 matrix `[[a,b],[b,c]]`, so it is
sharp for the two-block *scalar* problem. -/
theorem block_energy_certificate {n m : ℕ}
    {A : Fin n → Fin n → ℝ} {B : Fin n → Fin m → ℝ} {C : Fin m → Fin m → ℝ}
    {a c b : ℝ}
    (ha : ∀ u : Fin n → ℝ, a * (∑ i, u i ^ 2) ≤ ∑ i, ∑ j, u i * A i j * u j)
    (hc : ∀ v : Fin m → ℝ, c * (∑ i, v i ^ 2) ≤ ∑ i, ∑ j, v i * C i j * v j)
    (hb : ∀ (u : Fin n → ℝ) (v : Fin m → ℝ),
        |∑ i, ∑ j, u i * B i j * v j| ≤
          b * Real.sqrt (∑ i, u i ^ 2) * Real.sqrt (∑ j, v j ^ 2))
    (u : Fin n → ℝ) (v : Fin m → ℝ) :
    ((a + c - Real.sqrt ((a - c) ^ 2 + 4 * b ^ 2)) / 2) * blockNormSq u v ≤
      blockForm A B C u v := by
  set U := ∑ i, u i ^ 2 with hU
  set V := ∑ j, v j ^ 2 with hV
  have hU0 : 0 ≤ U := by rw [hU]; exact Finset.sum_nonneg (fun i _ => sq_nonneg (u i))
  have hV0 : 0 ≤ V := by rw [hV]; exact Finset.sum_nonneg (fun j _ => sq_nonneg (v j))
  have hx2 : Real.sqrt U ^ 2 = U := Real.sq_sqrt hU0
  have hy2 : Real.sqrt V ^ 2 = V := Real.sq_sqrt hV0
  have hQA : a * U ≤ ∑ i, ∑ j, u i * A i j * u j := by
    rw [hU]; exact ha u
  have hQC : c * V ≤ ∑ i, ∑ j, v i * C i j * v j := by
    rw [hV]; exact hc v
  have hW : |∑ i, ∑ j, u i * B i j * v j| ≤ b * Real.sqrt U * Real.sqrt V := by
    have := hb u v
    rwa [← hU, ← hV] at this
  have hW' : -(b * Real.sqrt U * Real.sqrt V) ≤
      ∑ i, ∑ j, u i * B i j * v j := (abs_le.mp hW).1
  have hsc := sharp_block_scalar_lower a c b (Real.sqrt U) (Real.sqrt V)
  dsimp only at hsc
  rw [hx2, hy2] at hsc
  have hrw : blockForm A B C u v =
      (∑ i, ∑ j, u i * A i j * u j) + 2 * (∑ i, ∑ j, u i * B i j * v j)
        + (∑ i, ∑ j, v i * C i j * v j) := rfl
  have hns : blockNormSq u v = U + V := by rw [blockNormSq, ← hU, ← hV]
  rw [hrw, hns]
  nlinarith [hsc, hW', hQA, hQC]

/-- **Block certificate (gap / Loewdin-Temple form).** When the complement gap
`g = c - a` is positive the sharper root is bounded below by the familiar
`a - b^2/g`, the second-order downfolding correction. This is the form quoted in
practice because it needs no square root. -/
theorem block_energy_certificate_gap {n m : ℕ}
    {A : Fin n → Fin n → ℝ} {B : Fin n → Fin m → ℝ} {C : Fin m → Fin m → ℝ}
    {a c b : ℝ} (hgap : 0 < c - a)
    (ha : ∀ u : Fin n → ℝ, a * (∑ i, u i ^ 2) ≤ ∑ i, ∑ j, u i * A i j * u j)
    (hc : ∀ v : Fin m → ℝ, c * (∑ i, v i ^ 2) ≤ ∑ i, ∑ j, v i * C i j * v j)
    (hb : ∀ (u : Fin n → ℝ) (v : Fin m → ℝ),
        |∑ i, ∑ j, u i * B i j * v j| ≤
          b * Real.sqrt (∑ i, u i ^ 2) * Real.sqrt (∑ j, v j ^ 2))
    (u : Fin n → ℝ) (v : Fin m → ℝ) :
    (a - b ^ 2 / (c - a)) * blockNormSq u v ≤ blockForm A B C u v := by
  set U := ∑ i, u i ^ 2 with hU
  set V := ∑ j, v j ^ 2 with hV
  have hU0 : 0 ≤ U := by rw [hU]; exact Finset.sum_nonneg (fun i _ => sq_nonneg (u i))
  have hV0 : 0 ≤ V := by rw [hV]; exact Finset.sum_nonneg (fun j _ => sq_nonneg (v j))
  have hx2 : Real.sqrt U ^ 2 = U := Real.sq_sqrt hU0
  have hy2 : Real.sqrt V ^ 2 = V := Real.sq_sqrt hV0
  have hQA : a * U ≤ ∑ i, ∑ j, u i * A i j * u j := by
    rw [hU]; exact ha u
  have hQC : c * V ≤ ∑ i, ∑ j, v i * C i j * v j := by
    rw [hV]; exact hc v
  have hW : |∑ i, ∑ j, u i * B i j * v j| ≤ b * Real.sqrt U * Real.sqrt V := by
    have := hb u v
    rwa [← hU, ← hV] at this
  have hW' : -(b * Real.sqrt U * Real.sqrt V) ≤
      ∑ i, ∑ j, u i * B i j * v j := (abs_le.mp hW).1
  have hsc := quadratic_gap_lower a c b (Real.sqrt U) (Real.sqrt V) hgap
  rw [hx2, hy2] at hsc
  have hrw : blockForm A B C u v =
      (∑ i, ∑ j, u i * A i j * u j) + 2 * (∑ i, ∑ j, u i * B i j * v j)
        + (∑ i, ∑ j, v i * C i j * v j) := rfl
  have hns : blockNormSq u v = U + V := by rw [blockNormSq, ← hU, ← hV]
  rw [hrw, hns]
  nlinarith [hsc, hW', hQA, hQC]

/-- Normalised form: a unit state in the full space has energy at least the
sharp two-block constant. -/
theorem normalized_energy_lower {n m : ℕ}
    {A : Fin n → Fin n → ℝ} {B : Fin n → Fin m → ℝ} {C : Fin m → Fin m → ℝ}
    {a c b : ℝ}
    (ha : ∀ u : Fin n → ℝ, a * (∑ i, u i ^ 2) ≤ ∑ i, ∑ j, u i * A i j * u j)
    (hc : ∀ v : Fin m → ℝ, c * (∑ i, v i ^ 2) ≤ ∑ i, ∑ j, v i * C i j * v j)
    (hb : ∀ (u : Fin n → ℝ) (v : Fin m → ℝ),
        |∑ i, ∑ j, u i * B i j * v j| ≤
          b * Real.sqrt (∑ i, u i ^ 2) * Real.sqrt (∑ j, v j ^ 2))
    {u : Fin n → ℝ} {v : Fin m → ℝ} (hnorm : blockNormSq u v = 1) :
    ((a + c - Real.sqrt ((a - c) ^ 2 + 4 * b ^ 2)) / 2) ≤ blockForm A B C u v := by
  have h := block_energy_certificate ha hc hb u v
  rwa [hnorm, mul_one] at h

end CalcChemistry
