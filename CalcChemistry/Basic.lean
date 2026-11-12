/-
# CalcChemistry.Basic

Foundational definitions for quantum chemistry formalization:
- Hilbert spaces for quantum states
- Molecular Hamiltonian structure
- Born-Oppenheimer approximation
-/

import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.LinearAlgebra.Matrix.Spectrum
import Mathlib.Data.Complex.Exponential

namespace QuantumChem

variable {𝕜 : Type*} [RCLike 𝕜]

/-- Hilbert space for N-particle quantum system -/
structure QuantumState (N : ℕ) where
  space : Type*
  [innerProd : InnerProductSpace 𝕜 space]
  [complete : CompleteSpace space]

/-- Molecular Hamiltonian structure: kinetic + potential -/
structure MolecularHamiltonian (N M : ℕ) where
  /-- Number of electrons -/
  n_electrons : N > 0
  /-- Number of nuclei -/
  n_nuclei : M > 0
  /-- Kinetic energy operator -/
  kinetic : LinearMap 𝕜 (QuantumState N).space (QuantumState N).space
  /-- Coulomb potential -/
  coulomb : LinearMap 𝕜 (QuantumState N).space (QuantumState N).space
  /-- Self-adjoint property -/
  selfAdjoint : IsSelfAdjoint (kinetic + coulomb)

/-- Born-Oppenheimer approximation: separate nuclear and electronic motion -/
axiom BornOppenheimer_separation {N M : ℕ} (H : MolecularHamiltonian N M) :
  ∃ (ε : ℝ), ε = (MassElectron / MassNucleus : ℝ) ∧ ε < 1/1000

/-- Placeholder for electron/nucleus mass ratio -/
axiom MassElectron : ℝ
axiom MassNucleus : ℝ
axiom mass_ratio_small : MassElectron / MassNucleus < 1/1000

end QuantumChem
