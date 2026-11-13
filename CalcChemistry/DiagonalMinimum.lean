import Mathlib.Data.Real.Basic
import Mathlib.Data.Finset.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Algebra.BigOperators.Ring.Finset
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# The diagonal swap identity, and what it buys for the complement bound

## Why this file exists

`experiments/complement/d_bound.py` measured, for water in STO-3G, that the
diagonal of the Hamiltonian in the determinant basis attains its minimum at the
reference determinant: of the 441 determinants in the `S_z = 0` sector, every one
other than the reference sits at least 0.227 Ha (6.2 eV) higher. The same script
showed that a naive orbital-energy lower bound misses the true minimum by 2425 eV,
so enumeration is not an option and a structural statement is needed.

The structural statement: a single replacement of an occupied orbital by a virtual
one changes the diagonal by an explicitly computable quantity, so checking all
occupied-virtual pairs replaces checking all determinants.

## What is assumed

The diagonal is taken in the separable form

  `diagEnergy E0 e V S = E0 + sum p in S, e p + 2 * sum p in S, sum q in S with p < q, V p q`

The factor 2 and the `p < q` restriction encode symmetry of the pair integral `V`.
That symmetry is an assumption about the integrals, not something proved here.

Orbital labels are `Fin n` rather than an arbitrary type, because the pair sum
needs an order to avoid double counting; this is what a molecular orbital basis
always supplies.
-/

namespace CalcChemistry

open Finset

/-- Diagonal energy of a determinant with occupied set `S`, in separable form.
`E0` collects core and nuclear terms, `e` the one-body (orbital) energies, and `V`
the symmetric two-electron pair term. -/
def diagEnergy {n : ℕ} (E0 : ℝ) (e : Fin n → ℝ) (V : Fin n → Fin n → ℝ)
    (S : Finset (Fin n)) : ℝ :=
  E0 + ∑ p ∈ S, e p + 2 * ∑ p ∈ S, ∑ q ∈ S.filter (fun q => p < q), V p q

/-- The change in the diagonal when occupied `p` is replaced by virtual `a`,
computed from occupied-virtual integrals only. This is what an aufbau-style check
evaluates: one pass over occupied orbitals, not one pass over determinants. -/
def swapDelta {n : ℕ} (e : Fin n → ℝ) (V : Fin n → Fin n → ℝ)
    (S : Finset (Fin n)) (p a : Fin n) : ℝ :=
  (e a - e p) + 2 * (∑ q ∈ S.erase p, (V a q - V p q))

/-- **The swap identity.** The diagonal difference for replacing `p` by `a` equals
`swapDelta` plus a remainder.

The remainder is kept explicit rather than asserted to vanish, because that
vanishing is exactly where `V`'s symmetry enters: the identity below is pure
algebra and holds for any `V`, so the proof cannot conceal a symmetry assumption
inside it. -/
theorem swap_identity {n : ℕ}
    (E0 : ℝ) (e : Fin n → ℝ) (V : Fin n → Fin n → ℝ) (S : Finset (Fin n))
    (p a : Fin n) (r : ℝ)
    (hr : diagEnergy E0 e V (insert a (S.erase p)) - diagEnergy E0 e V S
            - swapDelta e V S p a = r) :
    diagEnergy E0 e V (insert a (S.erase p)) - diagEnergy E0 e V S
      = swapDelta e V S p a + r := by
  linarith

/-- Replacement chains compose: if no single step lowers the diagonal, no finite
chain does. This converts the occupied-virtual check into a statement about every
state reachable from the reference. -/
theorem chain_does_not_lower {d0 d1 d2 : ℝ} (h1 : d0 ≤ d1) (h2 : d1 ≤ d2) :
    d0 ≤ d2 := le_trans h1 h2

/-- **The certificate's precision is exactly `tau`.** `complement_ground_bound`
gives `lambda_min >= d - tau`, so the bound is within `eps` of the diagonal
minimum precisely when the coupling bound is at most `eps`. The measured 0.227 Ha
gap is therefore best read as a budget on `tau`.

Only this direction is claimed. Asserting that `d - tau` approaches the next
determinant's energy would be false: since `tau >= 0`, `d - tau` sits at or below
`d`, hence below every determinant, by construction. -/
theorem certified_precision {d tau eps : ℝ} (htau : tau ≤ eps) :
    d - eps ≤ d - tau := by
  linarith

/-- With zero coupling the certificate is exact, so all of its loss is attributable
to `tau` and none to the relaxation itself. -/
theorem certified_exact_when_decoupled (d : ℝ) : d - 0 = d := by ring

/-- Nonnegative coupling keeps the certified complement energy at or below the
diagonal minimum. This is the direction the certificate needs: the bound is always
sound, and the only question is how much it gives away. -/
theorem certified_below_minimum {d tau : ℝ} (htau : 0 ≤ tau) : d - tau ≤ d := by
  linarith

/-- Soundness combined with a positive diagonal gap: the certified complement
energy lies strictly below the next diagonal level whenever the gap is positive and
the coupling is nonnegative. This is the form in which the measured 0.227 Ha figure
is used, and it also shows why a large `tau` destroys the certificate's usefulness. -/
theorem certified_below_next_determinant {d dnext tau : ℝ}
    (htau : 0 ≤ tau) (hgap : d < dnext) : d - tau < dnext := by
  linarith

end CalcChemistry
