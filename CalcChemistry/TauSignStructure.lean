import Mathlib.Data.Real.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum
import Mathlib.Algebra.Order.Field.Basic

/-!
# Why the coupling bound `tau` cannot be improved by changing basis

Companion to `notes/tau_sign_structure.md`.

## The question this file answers

`ComplementBound.complement_ground_bound` needs a constant `tau` with

  `|Q_O(x)| <= tau * ||x||^2`   for every `x`,

where `O` is the off-diagonal part of the complement Hamiltonian. Measured for
water in STO-3G, the smallest such constant is `-lambda_min(O) = 1.750272` Ha
one-sided and `||O||_2 = 2.026458` Ha two-sided, while the correlation energy the
certificate would need to resolve is only `0.049305` Ha: the certificate pays 35
to 41 times what the physics requires and is therefore vacuous.

The natural hope is that a change of basis, or a diagonal similarity
`O |-> S O S^{-1}` chosen to exploit the signs of `O`, reduces the constant. This
file isolates what can and cannot be formalised about that hope.

## What is proved here

1. `abs_quadform_le_absO` --- signs are invisible to a magnitude bound: the
   quadratic form of `O` is dominated by that of `|O|`. Every bound of the form
   `tau = f(|O|)` loses all sign cancellation at this step.

2. `two_mul_abs_le_weighted` and `quadform_le_diag_scaled` --- the
   diagonal-similarity-optimised (weighted Collatz--Wielandt) off-diagonal bound.
   For symmetric `O` and any positive weight `v` with weighted row sums dominated
   by `tau`, the quadratic form is dominated by `tau * ||x||^2`. Plain row sums are
   the case `v = 1` (`quadform_le_rowsum`).

3. `posOrthant_gain` --- the one place signs genuinely help: on the positive
   orthant the negative entries can be discarded. Converting this into a
   ground-energy bound needs the hypothesis that the ground vector is
   sign-definite, which is *not* proved here.

4. `two_by_two_scaling_cannot_help` --- a positive diagonal similarity cannot lower
   the 2-norm of a symmetric matrix, witnessed in dimension 2.

## What is NOT proved here, and why it matters

The strongest statement in the note is that for symmetric `O`,

  `inf over positive diagonal S of ||S O S^{-1}||_2  =  ||O||_2`,

because `S O S^{-1}` is similar to `O` and a symmetric matrix has 2-norm equal to
its spectral radius. That argument needs spectral theory of normal matrices
(eigenvalues of a similar matrix, and the equality of 2-norm and spectral radius
for symmetric matrices), which the available mathlib cache in this project does not
provide; it is a five-line textbook argument and is checked numerically instead
(20000 random scalings, minimum 59.33 against a floor of 2.026: every attempt is
strictly worse, none lower). The dimension-2 witness above is the part that is
kernel-checked. Do not read this file as a formal proof that no basis change helps.
-/

namespace CalcChemistry

open Finset

variable {n : ℕ}

/-- **Signs are invisible to the quadratic form.** The absolute value of the
quadratic form of `O` is bounded by the quadratic form of the entrywise absolute
value `|O|`. Every magnitude-only bound passes through this inequality, and this is
exactly where sign cancellation is discarded. -/
theorem abs_quadform_le_absO (O : Fin n → Fin n → ℝ) (x : Fin n → ℝ) :
    |∑ i, ∑ j, O i j * x i * x j| ≤ ∑ i, ∑ j, |O i j| * |x i| * |x j| := by
  calc |∑ i, ∑ j, O i j * x i * x j|
      ≤ ∑ i, |∑ j, O i j * x i * x j| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ i, ∑ j, |O i j * x i * x j| :=
        Finset.sum_le_sum (fun i _ => Finset.abs_sum_le_sum_abs _ _)
    _ = ∑ i, ∑ j, |O i j| * |x i| * |x j| := by
        refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
        rw [abs_mul, abs_mul]

/-- The pointwise AM–GM step behind the weighted bound: for positive weights `v`,
`2 |x i| |x j| <= (v j / v i) (x i)^2 + (v i / v j) (x j)^2`. -/
theorem two_mul_abs_le_weighted (vi vj : ℝ) (hvi : 0 < vi) (hvj : 0 < vj) (xi xj : ℝ) :
    2 * |xi| * |xj| ≤ (vj / vi) * xi ^ 2 + (vi / vj) * xj ^ 2 := by
  have hpos : 0 < vi * vj := mul_pos hvi hvj
  have key : 2 * |xi| * |xj| * (vi * vj) ≤ vj ^ 2 * xi ^ 2 + vi ^ 2 * xj ^ 2 := by
    have hsq : 0 ≤ (vj * |xi| - vi * |xj|) ^ 2 := sq_nonneg _
    have hxi : |xi| ^ 2 = xi ^ 2 := sq_abs xi
    have hxj : |xj| ^ 2 = xj ^ 2 := sq_abs xj
    nlinarith [hsq, hxi, hxj]
  have e : (vj / vi) * xi ^ 2 + (vi / vj) * xj ^ 2
      = (vj ^ 2 * xi ^ 2 + vi ^ 2 * xj ^ 2) / (vi * vj) := by
    field_simp
    ring
  rw [e, le_div_iff₀ hpos]
  nlinarith [key]

/-- **Diagonal-similarity-optimised (weighted Collatz--Wielandt) off-diagonal bound.**
For symmetric `O` and any positive weight vector `v`, if `tau` dominates the
`v`-weighted row sums of `|O|`, then `tau` bounds the quadratic form.

The optimal weights are the Perron vector of `|O|`, which makes the weighted row sum
equal `rho(|O|)`; that optimality is Perron--Frobenius and is not formalised here.
Note this bound is stated for `|O|`, so the sign pattern of `O` still does not enter.
-/
theorem quadform_le_diag_scaled (O : Fin n → Fin n → ℝ)
    (hO : ∀ i j, O j i = O i j) (v : Fin n → ℝ) (hv : ∀ i, 0 < v i) (τ : ℝ)
    (hτ : ∀ i, ∑ j, |O i j| * v j ≤ τ * v i) (x : Fin n → ℝ) :
    |∑ i, ∑ j, O i j * x i * x j| ≤ τ * ∑ i, x i ^ 2 := by
  have h1 : |∑ i, ∑ j, O i j * x i * x j| ≤ ∑ i, ∑ j, |O i j| * |x i| * |x j| :=
    abs_quadform_le_absO O x
  have hpair : ∀ i j : Fin n,
      (2:ℝ) * |O i j| * |x i| * |x j| ≤
        (|O i j| * v j / v i) * x i ^ 2 + (|O i j| * v i / v j) * x j ^ 2 := by
    intro i j
    have h := two_mul_abs_le_weighted (v i) (v j) (hv i) (hv j) (x i) (x j)
    have hOij : 0 ≤ |O i j| := abs_nonneg _
    calc (2:ℝ) * |O i j| * |x i| * |x j| = |O i j| * (2 * |x i| * |x j|) := by ring
      _ ≤ |O i j| * ((v j / v i) * x i ^ 2 + (v i / v j) * x j ^ 2) :=
          mul_le_mul_of_nonneg_left h hOij
      _ = (|O i j| * v j / v i) * x i ^ 2 + (|O i j| * v i / v j) * x j ^ 2 := by ring
  have hswap : (∑ i, ∑ j, (|O i j| * v i / v j) * x j ^ 2)
      = ∑ i, ∑ j, (|O i j| * v j / v i) * x i ^ 2 := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl (fun i _ => Finset.sum_congr rfl (fun j _ => ?_))
    rw [hO j i]
  have hdiv : ∀ i : Fin n, (∑ j, |O i j| * v j) / v i = ∑ j, (|O i j| * v j / v i) := by
    intro i
    rw [div_eq_mul_inv, Finset.sum_mul]
    exact Finset.sum_congr rfl (fun j _ => by rw [div_eq_mul_inv])
  have hR : (∑ i, ∑ j, (|O i j| * v j / v i) * x i ^ 2)
      = ∑ i, x i ^ 2 * ((∑ j, |O i j| * v j) / v i) := by
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [hdiv i, ← Finset.sum_mul]
    ring
  have hsum : 2 * (∑ i, ∑ j, |O i j| * |x i| * |x j|)
      ≤ 2 * (∑ i, ∑ j, (|O i j| * v j / v i) * x i ^ 2) := by
    have h2S : 2 * (∑ i, ∑ j, |O i j| * |x i| * |x j|)
        = ∑ i, ∑ j, (2:ℝ) * |O i j| * |x i| * |x j| := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl (fun i _ => ?_)
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl (fun j _ => by ring)
    have hsplit : 2 * (∑ i, ∑ j, (|O i j| * v j / v i) * x i ^ 2)
        = (∑ i, ∑ j, (|O i j| * v j / v i) * x i ^ 2)
          + ∑ i, ∑ j, (|O i j| * v i / v j) * x j ^ 2 := by
      rw [hswap]; ring
    rw [h2S, hsplit, ← Finset.sum_add_distrib]
    refine Finset.sum_le_sum (fun i _ => ?_)
    rw [← Finset.sum_add_distrib]
    exact Finset.sum_le_sum (fun j _ => hpair i j)
  have hfinal : (∑ i, x i ^ 2 * ((∑ j, |O i j| * v j) / v i)) ≤ τ * ∑ i, x i ^ 2 := by
    rw [Finset.mul_sum]
    refine Finset.sum_le_sum (fun i _ => ?_)
    have hdivτ : (∑ j, |O i j| * v j) / v i ≤ τ := by
      rw [div_le_iff₀ (hv i)]
      simpa using hτ i
    have := mul_le_mul_of_nonneg_left hdivτ (sq_nonneg (x i))
    nlinarith [this]
  rw [hR] at hsum
  linarith

/-- **Row-sum (Gershgorin/Schur) bound**: the `v = 1` case, i.e. the bound actually
measured at `3.758707` Ha for water, against a true norm of `2.026458` Ha. -/
theorem quadform_le_rowsum (O : Fin n → Fin n → ℝ) (hO : ∀ i j, O j i = O i j)
    (τ : ℝ) (hτ : ∀ i, ∑ j, |O i j| ≤ τ) (x : Fin n → ℝ) :
    |∑ i, ∑ j, O i j * x i * x j| ≤ τ * ∑ i, x i ^ 2 := by
  refine quadform_le_diag_scaled O hO (fun _ => 1) (fun _ => one_pos) τ ?_ x
  intro i
  simpa using hτ i

/-- **The one genuine sign gain, on the positive orthant.** If `x i * x j ≥ 0` for
every pair (in particular if `x ≥ 0`), the negative entries of `O` can be discarded
entirely: only the entrywise positive part survives.

Turning this into a ground-energy bound requires knowing the ground vector is
sign-definite, i.e. a Perron--Frobenius hypothesis on `O`'s sign pattern. That
hypothesis is not established here, and for a determinant-basis Hamiltonian with
mixed-sign off-diagonals it is generally false. -/
theorem posOrthant_gain (O : Fin n → Fin n → ℝ) (x : Fin n → ℝ)
    (hx : ∀ i j, 0 ≤ x i * x j) :
    ∑ i, ∑ j, O i j * x i * x j ≤ ∑ i, ∑ j, max (O i j) 0 * x i * x j := by
  refine Finset.sum_le_sum (fun i _ => Finset.sum_le_sum (fun j _ => ?_))
  have hxij := hx i j
  rcases le_or_lt 0 (O i j) with h | h
  · rw [max_eq_left h]
  · rw [max_eq_right (le_of_lt h)]
    nlinarith

/-- **A positive diagonal similarity cannot lower the 2-norm of a symmetric matrix**,
witnessed in dimension 2. For `O = [[0,1],[1,0]]` and `S = diag(a,b)` the Rayleigh
quotient of `S O S^{-1}` at `x = (1,1)` is `(a/b + b/a)/2`, never below the unscaled
value `1`.

This is the kernel-checked core of the general no-go statement, whose full form
(`inf_S ||S O S^{-1}||_2 = ||O||_2` for symmetric `O`) is a spectral argument not
available in this project's mathlib cache; see the module docstring. -/
theorem two_by_two_scaling_cannot_help (a b : ℝ) (ha : 0 < a) (hb : 0 < b) :
    1 ≤ (a / b + b / a) / 2 := by
  have hab : 0 < a * b := mul_pos ha hb
  rw [le_div_iff₀ (by norm_num : (0:ℝ) < 2), one_mul]
  have e : a / b + b / a = (a ^ 2 + b ^ 2) / (a * b) := by
    field_simp
    ring
  rw [e, le_div_iff₀ hab]
  nlinarith [sq_nonneg (a - b)]

/-- **The certificate's price is the constant, not the split.** Whatever `tau` is
chosen, the certified complement energy is `d - tau`, so the shortfall against the
true ground energy is the correlation energy plus `tau`. This records why a `tau`
thirty-five times the correlation energy cannot produce a useful bound,
independently of how `d` was obtained. -/
theorem certificate_shortfall {d tau e_true : ℝ} :
    e_true - (d - tau) = tau - (d - e_true) := by
  ring

end CalcChemistry
