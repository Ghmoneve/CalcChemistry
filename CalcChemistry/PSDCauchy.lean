import Mathlib.LinearAlgebra.Matrix.PosDef
import Mathlib.Tactic

/-!
# Cauchy--Schwarz in the inner product of a positive semidefinite matrix

## Why this file exists

`SchurStep.schur_block_bound_shifted` proves that the block quadratic form
`d x^2 + 2 x <b,w> + <A w,w>` is bounded below by `(lambda - s)(x^2 + ||w||^2)`
**given** the scalar hypothesis `<b,w>^2 <= s (<A w,w> - lambda ||w||^2)` for every
`w`. The step from `A > 0` to that hypothesis was proved on paper (Proposition 1 of
the manuscript) and checked numerically in `verified_cholesky.py`, but was outside
the kernel. This file closes it.

## The mathematics

For `A` positive semidefinite the bilinear form `<u,v>_A = u^T A v` satisfies the
Cauchy--Schwarz inequality

  `<x,y>_A^2 <= <x,x>_A <y,y>_A`,

proved by the discriminant argument: the quadratic `t |-> <x+ty, x+ty>_A` is
nonnegative for every real `t`, so its discriminant is at most zero.

Taking `y = A^{-1} b` turns this into `<b,w>^2 <= (b^T A^{-1} b)(w^T A w)`, which is
the hypothesis of the Schur step with `s = b^T A^{-1} b`.

## What is proved here

* `disc_le` --- the scalar discriminant lemma.
* `posSemidef_cauchy` --- the Cauchy--Schwarz inequality for a positive semidefinite
  matrix, in the real case.

## What is not proved here

The identification `s = b^T A^{-1} b` and the substitution `y = A^{-1}b` are left to
the caller: this file supplies the inequality, not the choice of `s`.
-/

namespace CalcChemistry

open Matrix Finset

/-- **Scalar discriminant lemma.** If `a + 2 t b + t^2 c` is nonnegative for every
real `t` and `c >= 0`, then `b^2 <= a c`. No hypothesis on `a` is needed: for `c > 0`
the vertex argument gives `b^2/c <= a` directly, and for `c = 0` nonnegativity for all
`t` forces `b = 0`.

For `c > 0` the minimum sits at `t = -b/c`; for `c = 0` the expression is affine and
nonnegativity for all `t` forces `b = 0`. -/
theorem disc_le {a b c : ℝ} (hc : 0 ≤ c)
    (h : ∀ t : ℝ, 0 ≤ a + 2 * t * b + t ^ 2 * c) : b ^ 2 ≤ a * c := by
  rcases eq_or_lt_of_le hc with hc0 | hcpos
  · -- c = 0: the expression is affine in t, so b must vanish
    subst hc0
    have hb : b = 0 := by
      by_contra hb
      have hinst := h (-(a + 1) / (2 * b))
      have hval : a + 2 * (-(a + 1) / (2 * b)) * b + (-(a + 1) / (2 * b)) ^ 2 * (0:ℝ)
          = -1 := by
        field_simp
        ring
      rw [hval] at hinst
      linarith
    rw [hb]
    simp
  · -- c > 0: evaluate at the vertex
    have hinst := h (-b / c)
    have hval : a + 2 * (-b / c) * b + (-b / c) ^ 2 * c = a - b ^ 2 / c := by
      field_simp
      ring
    rw [hval] at hinst
    have h1 : b ^ 2 / c ≤ a := by linarith
    rw [div_le_iff₀ hcpos] at h1
    linarith


/-- Entrywise symmetry of a positive semidefinite real matrix. -/
theorem posSemidef_symm_apply {n : ℕ} {A : Matrix (Fin n) (Fin n) ℝ}
    (hA : A.PosSemidef) (i j : Fin n) : A i j = A j i := by
  have h := hA.isHermitian.apply i j
  simpa using h.symm

/-- **Cauchy--Schwarz for a positive semidefinite matrix.** For real `A ⪰ 0` and any
`x`, `y`, the bilinear form `u ⬝ᵥ A.mulVec v` satisfies `<x,y>_A^2 <= <x,x>_A <y,y>_A`.

This is the inequality that turns `A ≻ 0` into the scalar hypothesis consumed by the
Schur step in `SchurStep.lean`. -/
theorem posSemidef_cauchy {n : ℕ} (A : Matrix (Fin n) (Fin n) ℝ) (hA : A.PosSemidef)
    (x y : Fin n → ℝ) :
    (x ⬝ᵥ A.mulVec y) ^ 2 ≤ (x ⬝ᵥ A.mulVec x) * (y ⬝ᵥ A.mulVec y) := by
  have hnn : ∀ z : Fin n → ℝ, 0 ≤ z ⬝ᵥ A.mulVec z := by
    intro z
    have h := hA.re_dotProduct_nonneg z
    simpa using h
  -- the form is symmetric in its two arguments
  have hsym : ∀ u v : Fin n → ℝ, u ⬝ᵥ A.mulVec v = v ⬝ᵥ A.mulVec u := by
    intro u v
    simp only [dotProduct, mulVec, Finset.mul_sum, Finset.sum_mul, mul_assoc]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
    rw [posSemidef_symm_apply hA j i]
    ring
  -- expanding the quadratic form along the line x + t y
  have hexp : ∀ t : ℝ,
      (x + t • y) ⬝ᵥ A.mulVec (x + t • y)
        = (x ⬝ᵥ A.mulVec x) + 2 * t * (x ⬝ᵥ A.mulVec y) + t ^ 2 * (y ⬝ᵥ A.mulVec y) := by
    intro t
    rw [mulVec_add, mulVec_smul]
    simp only [dotProduct_add, add_dotProduct, dotProduct_smul, smul_dotProduct,
               smul_eq_mul]
    rw [hsym y x]
    ring
  exact disc_le (hnn y) (fun t => by rw [← hexp t]; exact hnn (x + t • y))

/-- **The hypothesis consumed by the Schur step, derived from `A` positive definite.**

For a positive definite real `A` and any `b`, the scalar `s = b^T A^{-1} b` makes the
per-`w` inequality `<b,w>^2 <= s <Aw,w>` true. This is `posSemidef_cauchy` applied
with `x = A^{-1} b` and `y = w`: the left factor collapses to `<b,w>` and the middle
factor to `s`, using `A^T = A` and `A A^{-1} = 1`.

This is exactly the hypothesis that `SchurStep.schur_block_bound_shifted` assumes. -/
theorem schur_scalar_hypothesis {n : ℕ} {A : Matrix (Fin n) (Fin n) ℝ}
    (hA : A.PosDef) (b w : Fin n → ℝ) :
    (b ⬝ᵥ w) ^ 2
      ≤ (b ⬝ᵥ A⁻¹.mulVec b) * (w ⬝ᵥ A.mulVec w) := by
  have hdet : IsUnit A.det := isUnit_iff_ne_zero.mpr (ne_of_gt hA.det_pos)
  have ht : Aᵀ = A := by
    simpa using hA.isHermitian.eq
  have hleft : Aᵀ.mulVec (A⁻¹.mulVec b) = b := by
    rw [Matrix.mulVec_mulVec, ht, Matrix.mul_nonsing_inv A hdet, Matrix.one_mulVec]
  -- <A^{-1}b, Aw> = <b, w>
  have h1 : A⁻¹.mulVec b ⬝ᵥ A.mulVec w = b ⬝ᵥ w := by
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hleft]
  -- <A^{-1}b, A A^{-1}b> = <b, A^{-1}b>
  have h2 : A⁻¹.mulVec b ⬝ᵥ A.mulVec (A⁻¹.mulVec b) = b ⬝ᵥ A⁻¹.mulVec b := by
    rw [Matrix.dotProduct_mulVec, ← Matrix.mulVec_transpose, hleft]
  have hcs := posSemidef_cauchy A hA.posSemidef (A⁻¹.mulVec b) w
  rw [h1, h2] at hcs
  exact hcs

end CalcChemistry
