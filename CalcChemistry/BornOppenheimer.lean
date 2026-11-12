import Mathlib.Analysis.SpecialFunctions.Pow.Real

namespace BornOppenheimer

noncomputable def parameter (electronMass nuclearMass : ℝ) : ℝ :=
  (electronMass / nuclearMass) ^ (1 / 4 : ℝ)

/-- Only positivity of the mass parameter is asserted here. No molecular
spectral gap or Born--Oppenheimer error estimate is postulated. -/
theorem parameter_nonneg (electronMass nuclearMass : ℝ)
    (he : 0 ≤ electronMass) (hn : 0 < nuclearMass) :
    0 ≤ parameter electronMass nuclearMass := by
  exact Real.rpow_nonneg (div_nonneg he (le_of_lt hn)) _

/-- A bound supplied by an independent analysis, not a BO theorem. The units
of constant and error must match. No instance for CaF or CaF2 is supplied. -/
structure FifthOrderCertificate where
  epsilon : ℝ
  constant : ℝ
  error : ℝ
  epsilon_nonneg : 0 ≤ epsilon
  constant_nonneg : 0 ≤ constant
  bound : |error| ≤ constant * epsilon ^ 5

end BornOppenheimer
