import Mathlib.Data.Real.Sqrt
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

open Finset

namespace Probe

def blockform {n m : ℕ} (A : Fin n → Fin n → ℝ) (B : Fin n → Fin m → ℝ)
    (Cc : Fin m → Fin m → ℝ) (u : Fin n → ℝ) (v : Fin m → ℝ) : ℝ :=
  (∑ i, ∑ j, u i * A i j * u j) + 2 * (∑ i, ∑ j, u i * B i j * v j)
    + (∑ i, ∑ j, v i * Cc i j * v j)

example {n : ℕ} (u : Fin n → ℝ) : 0 ≤ ∑ i, u i ^ 2 :=
  Finset.sum_nonneg (fun i _ => sq_nonneg _)

example {n : ℕ} (u : Fin n → ℝ) : (Real.sqrt (∑ i, u i ^ 2)) ^ 2 = ∑ i, u i ^ 2 :=
  Real.sq_sqrt (Finset.sum_nonneg (fun i _ => sq_nonneg _))

end Probe
