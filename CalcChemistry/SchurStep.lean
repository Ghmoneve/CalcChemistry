import Mathlib.Data.Real.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum

/-!
# The Schur step: one scalar inequality buys an n-dimensional bound

## The statement

Split a symmetric operator so that a single reference coordinate `x` is separated
from the remaining `n` coordinates `w`:

  `H(x, w) = d x^2 + 2 x <b, w> + <A w, w>`,

with `d` the reference diagonal, `b` the coupling and `A` the remaining block.

The theorem below says: if the coupling satisfies ONE scalar inequality,
`<b,w>^2 <= s * (<A w,w> - lambda <w,w>)` for every `w` and some `s > 0`, then

  `H(x, w) >= (lambda - s) (x^2 + <w,w>)`   whenever `lambda <= d`.

So an n-dimensional positivity question collapses to a single scalar bound. This is
the step that separates the Schur certificate from the block certificate: the block
certificate needs the operator norm of the whole off-diagonal part, which is fixed by
the spectrum and cannot be lowered by a change of basis, while this needs only the
scalar `<b, A^{-1} b>`, which is one row of the problem.

## What is proved

* `two_abs_mul_le` --- the scaled AM-GM that does all the work.
* `schur_block_bound_shifted` --- the theorem above, for any admissible shift.
* `schur_block_bound` --- the case `lambda = d`, which is the form used in the
  measured certificate.

## What is NOT proved, and is not implied

Nothing here establishes that any particular scalar `s` exists. For `A - lambda I`
positive definite the sharp choice is `s = lambda_max(B^T (A - lambda I)^{-1} B)`,
and for the water / STO-3G system that sharp value is `0.050675`, giving a certified
lower bound of `-75.01338150 Ha` against `E_FCI = -75.01201163 Ha`. Producing the
scalar is a numerical obligation, discharged elsewhere; this file only proves that
once it is produced, the bound follows.

In particular the exponential wall is untouched: bounding the scalar without forming
`A` is the open problem, and it is not addressed by anything here.
-/

namespace CalcChemistry

open Finset

/-- **Scaled AM-GM.** For `s > 0`, `2 a b <= s a^2 + b^2 / s`.

Obtained from `(s a - b)^2 >= 0` divided by `s`. This is the only analytic input the
Schur step needs. -/
theorem two_abs_mul_le {s a b : ℝ} (hs : 0 < s) :
    2 * a * b ≤ s * a ^ 2 + b ^ 2 / s := by
  have hsq : 0 ≤ (s * a - b) ^ 2 := sq_nonneg _
  have hdiv : 0 ≤ (s * a - b) ^ 2 / s := div_nonneg hsq (le_of_lt hs)
  have heq : (s * a - b) ^ 2 / s = s * a ^ 2 - 2 * a * b + b ^ 2 / s := by
    field_simp
    ring
  linarith [heq ▸ hdiv]

/-- **The Schur step, with a shift.**

Let `H(x,w) = d x^2 + 2 x <b,w> + <A w,w>`. Suppose a scalar `s > 0` dominates the
coupling in the sense

  `<b,w>^2 <= s * (<A w,w> - lambda <w,w>)`   for every `w`,

and suppose `lambda <= d`. Then `H(x,w) >= (lambda - s) (x^2 + <w,w>)`.

The proof is the AM-GM bound applied to `2 x <b,w>`: it costs `s x^2 + <b,w>^2/s`, and
the hypothesis converts the second term into `<A w,w> - lambda <w,w>`, which cancels
against the block term exactly, leaving `(d - lambda) x^2 + s <w,w>`.

The hypothesis is where all the numerical work lives: it asks for one scalar, not for
a norm of the whole off-diagonal part. -/
theorem schur_block_bound_shifted {n : ℕ}
    (A : Fin n → Fin n → ℝ) (b : Fin n → ℝ) (d lam s : ℝ)
    (hs : 0 < s)
    (hlam : lam ≤ d)
    (hb : ∀ w : Fin n → ℝ,
      (∑ i, b i * w i) ^ 2
        ≤ s * ((∑ i, ∑ j, w i * A i j * w j) - lam * ∑ i, w i ^ 2))
    (x : ℝ) (w : Fin n → ℝ) :
    (lam - s) * (x ^ 2 + ∑ i, w i ^ 2)
      ≤ d * x ^ 2 + 2 * x * (∑ i, b i * w i) + (∑ i, ∑ j, w i * A i j * w j) := by
  set B := ∑ i, b i * w i with hB
  set Q := ∑ i, ∑ j, w i * A i j * w j with hQ
  set N := ∑ i, w i ^ 2 with hNdef
  have hbw : B ^ 2 ≤ s * (Q - lam * N) := hb w
  have ham := two_abs_mul_le (s := s) (a := |x|) (b := |B|) hs
  have hB2 : |B| ^ 2 / s ≤ Q - lam * N := by
    rw [div_le_iff₀ hs, sq_abs B]
    linarith [hbw]
  have h1 : 2 * |x| * |B| ≤ s * x ^ 2 + (Q - lam * N) := by
    rw [sq_abs x] at ham
    linarith [ham, hB2]
  have habs : |2 * x * B| ≤ s * x ^ 2 + (Q - lam * N) := by
    have htwo : |(2 : ℝ)| = 2 := by norm_num
    rw [abs_mul, abs_mul]
    calc |2| * |x| * |B| = 2 * |x| * |B| := by rw [htwo]
      _ ≤ s * x ^ 2 + (Q - lam * N) := h1
  have h2 : 2 * x * B ≥ -(s * x ^ 2 + (Q - lam * N)) := by
    linarith [neg_abs_le (2 * x * B), habs]
  have hN : 0 ≤ N := by
    rw [hNdef]
    exact Finset.sum_nonneg (fun i _ => sq_nonneg _)
  have hdl : 0 ≤ (d - lam) * x ^ 2 := mul_nonneg (sub_nonneg.mpr hlam) (sq_nonneg x)
  have hsN : 0 ≤ s * N := mul_nonneg (le_of_lt hs) hN
  nlinarith [h2, hdl, hsN]

/-- **The Schur step at the reference energy.** The case `lambda = d` of the shifted
form, which is what the measured certificate uses: the reference diagonal is fixed by
Rayleigh--Ritz, so no shift is needed and none is wasted. -/
theorem schur_block_bound {n : ℕ}
    (A : Fin n → Fin n → ℝ) (b : Fin n → ℝ) (d s : ℝ)
    (hs : 0 < s)
    (hb : ∀ w : Fin n → ℝ,
      (∑ i, b i * w i) ^ 2
        ≤ s * ((∑ i, ∑ j, w i * A i j * w j) - d * ∑ i, w i ^ 2))
    (x : ℝ) (w : Fin n → ℝ) :
    (d - s) * (x ^ 2 + ∑ i, w i ^ 2)
      ≤ d * x ^ 2 + 2 * x * (∑ i, b i * w i) + (∑ i, ∑ j, w i * A i j * w j) :=
  schur_block_bound_shifted A b d d s hs le_rfl hb x w

/-- **No coupling means no loss.** If the reference is decoupled (`b = 0`), the only
admissible constraint is `0 ≤ s (Q - lambda N)`, and the bound loses exactly `s`. As
`s -> 0` the bound collapses onto `lambda`, which is the statement that a decoupled
reference is an exact eigenvector. Recorded as the limiting case, since the measured
`0.86 kcal/mol` is exactly this loss. -/
theorem schur_zero_coupling_exact {n : ℕ}
    (A : Fin n → Fin n → ℝ) (d s : ℝ) (hs : 0 ≤ s)
    (hA : ∀ w : Fin n → ℝ,
      0 ≤ (∑ i, ∑ j, w i * A i j * w j) - d * ∑ i, w i ^ 2)
    (x : ℝ) (w : Fin n → ℝ) :
    (d - s) * (x ^ 2 + ∑ i, w i ^ 2)
      ≤ d * x ^ 2 + (∑ i, ∑ j, w i * A i j * w j) := by
  have hq := hA w
  have hN : 0 ≤ ∑ i, w i ^ 2 := Finset.sum_nonneg (fun i _ => sq_nonneg _)
  have hsx : 0 ≤ s * x ^ 2 := mul_nonneg hs (sq_nonneg x)
  have hsN : 0 ≤ s * (∑ i, w i ^ 2) := mul_nonneg hs hN
  nlinarith [hq, hsx, hsN]

end CalcChemistry
