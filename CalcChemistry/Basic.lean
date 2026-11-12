import Mathlib.Data.Real.Basic

namespace QuantumChem

/-- A finite-model energy certificate over the reals. Soundness is an explicit
input, not an assertion that an uncomputed molecular energy has been bounded. -/
structure EnergyCertificate where
  lower : ℝ
  upper : ℝ
  target : ℝ
  ordered : lower ≤ upper
  lower_sound : lower ≤ target
  upper_sound : target ≤ upper

def EnergyCertificate.Certifies (c : EnergyCertificate) (x : ℝ) : Prop :=
  c.lower ≤ x ∧ x ≤ c.upper

end QuantumChem
