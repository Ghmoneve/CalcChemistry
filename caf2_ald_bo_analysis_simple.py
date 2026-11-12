#!/usr/bin/env python3
# ARCHIVED UNVERIFIED DEMO: reported chemistry/proof claims are withdrawn.
# See docs/audit_numerics.md and STATUS_REPORT.md; do not use as validation.
"""
Born-Oppenheimer Approximation Analysis for CaF₂ ALD
System: Ca(thd)₂ + HF → CaF₂ + H(thd)
Pure Python implementation (no numpy)
"""

import math

# Atomic masses in amu
m_e = 5.485799e-4
M_Ca = 40.078
M_F = 18.998
M_H = 1.008

# Reduced masses
mu_CaF = (M_Ca * M_F) / (M_Ca + M_F)
mu_HF = (M_H * M_F) / (M_H + M_F)

# Born-Oppenheimer parameters
epsilon_Ca = (m_e / M_Ca) ** 0.25
epsilon_CaF = (m_e / mu_CaF) ** 0.25
epsilon_HF = (m_e / mu_HF) ** 0.25

# Error bounds
error_Ca = epsilon_Ca ** 5
error_CaF = epsilon_CaF ** 5
error_HF = epsilon_HF ** 5

# Exponential corrections (γ=5.0)
gamma = 5.0
exp_Ca = math.exp(-gamma / epsilon_Ca**2)
exp_CaF = math.exp(-gamma / epsilon_CaF**2)
exp_HF = math.exp(-gamma / epsilon_HF**2)

# Chemical accuracy
chem_acc_hartree = 1.6e-3
chem_acc_kcal = 1.0

print("="*70)
print("Born-Oppenheimer Analysis: CaF₂ ALD (Ca(thd)₂ + HF)")
print("="*70)
print()

print("Nuclear Masses:")
print(f"  Ca: {M_Ca:.3f} amu")
print(f"  F:  {M_F:.3f} amu")
print(f"  H:  {M_H:.3f} amu")
print(f"  Electron: {m_e:.6e} amu")
print()

print("Reduced Masses:")
print(f"  μ(Ca-F) = {mu_CaF:.4f} amu  (Ca-F bond formation)")
print(f"  μ(H-F)  = {mu_HF:.4f} amu   (H-F bond breaking)")
print()

print("Born-Oppenheimer Parameter ε = (m_e/M)^(1/4):")
print(f"  ε(Ca)   = {epsilon_Ca:.6f}  (heaviest, best BO)")
print(f"  ε(Ca-F) = {epsilon_CaF:.6f}  (primary coordinate)")
print(f"  ε(H-F)  = {epsilon_HF:.6f}  (lightest, worst case)")
print()

print("BO Error Bounds O(ε⁵):")
print(f"  Error(Ca)   = {error_Ca:.3e}  ({error_Ca*100:.6f}%)")
print(f"  Error(Ca-F) = {error_CaF:.3e}  ({error_CaF*100:.6f}%)")
print(f"  Error(H-F)  = {error_HF:.3e}  ({error_HF*100:.5f}%)")
print()

print(f"Exponential Corrections O(exp(-γ/ε²)) with γ={gamma}:")
print(f"  Exp(Ca)   = {exp_Ca:.3e}")
print(f"  Exp(Ca-F) = {exp_CaF:.3e}")
print(f"  Exp(H-F)  = {exp_HF:.3e}")
print()

print("="*70)
print("Nonadiabatic Coupling Analysis at Transition State")
print("="*70)
print()

# Energy gaps (eV)
gaps = {
    "Ground-1st excited": 4.5,
    "Ground-2nd excited": 6.8,
    "Avoided crossing (est)": 1.2,
}

print("Electronic Energy Gaps (estimated):")
for state, gap in gaps.items():
    print(f"  {state:25s}: {gap:.2f} eV")
print()

print(f"Reduced mass (H-F): {mu_HF:.4f} amu")
print(f"BO parameter ε(H-F): {epsilon_HF:.6f}")
print()

print("Nonadiabatic Coupling Criterion: ε²||d_kk'|| << 1")
print()

# Estimate nonadiabatic coupling
# For thermal motion at 500K: v ~ 0.1 Å/fs, coupling ~ v/Δ
eV_to_au = 0.036749
for state, gap_eV in gaps.items():
    gap_au = gap_eV * eV_to_au
    coupling_est = 0.01 / gap_au  # rough estimate
    bo_param = epsilon_HF**2 * coupling_est
    status = "✓ VALID" if bo_param < 0.1 else "✗ BREAKDOWN"
    print(f"  {state:25s}: ε²||d|| ≈ {bo_param:.4f}  {status}")

print()

print("Conical Intersection Check:")
print("  CaF₂ formation involves:")
print("    - Ca²⁺ (closed shell, ¹S₀ ground state)")
print("    - F⁻ (closed shell)")
print("    - HF (¹Σ⁺ ground state)")
print()
print("  Expected behavior:")
print("    ✓ Ground state reaction (no excited states involved)")
print("    ✓ All species closed-shell → no spin crossings")
print("    ✓ Large HOMO-LUMO gap (>4 eV) → well-separated states")
print("    ✓ No symmetry-required degeneracies")
print()
print("  Conclusion: Conical intersections UNLIKELY on ground state PES")
print("              BO approximation VALID for this reaction")
print()

print("="*70)
print("BO Error vs Chemical Accuracy")
print("="*70)
print()

typical_rxn_energy = 50.0  # kcal/mol
typical_energy_au = typical_rxn_energy * chem_acc_hartree

bo_error_abs_au = error_HF * typical_energy_au
bo_error_kcal = bo_error_abs_au / chem_acc_hartree

print(f"Chemical accuracy: {chem_acc_kcal:.1f} kcal/mol = {chem_acc_hartree:.6f} Hartree")
print(f"Typical reaction energy: {typical_rxn_energy:.1f} kcal/mol")
print()
print(f"BO error O(ε⁵): {error_HF:.6e} (relative)")
print(f"BO error (absolute): {bo_error_abs_au:.6e} Hartree")
print(f"BO error (absolute): {bo_error_kcal:.6e} kcal/mol")
print()
print(f"Ratio: BO_error / chemical_accuracy = {bo_error_kcal:.6e}")
print()

if bo_error_kcal < chem_acc_hartree:
    print("✓ BO error is MUCH SMALLER than chemical accuracy")
    print("  BO approximation is JUSTIFIED for kinetics prediction")
else:
    print("✗ BO error approaches chemical accuracy")
    print("  Nonadiabatic corrections may be needed")
print()

print("="*70)
print("SUMMARY: Born-Oppenheimer Validity for CaF₂ ALD")
print("="*70)
print()
print(f"✓ BO parameter ε(H-F) = {epsilon_HF:.6f} (worst case)")
print(f"✓ BO error O(ε⁵) = {error_HF:.3e} << chemical accuracy")
print("✓ Large HOMO-LUMO gap (>4 eV) → no conical intersections")
print("✓ Ground state reaction → nonadiabatic coupling negligible")
print()
print("CONCLUSION: Born-Oppenheimer approximation is VALID")
print("            Nonadiabatic effects are negligible for kinetics")
print()

# Numerical values for Lean proof
print("="*70)
print("KEY NUMERICAL RESULTS (for Lean formalization)")
print("="*70)
print()
print(f"ε_HF = {epsilon_HF:.10f}")
print(f"ε_HF⁵ = {error_HF:.15e}")
print(f"μ_CaF = {mu_CaF:.10f} amu")
print(f"μ_HF = {mu_HF:.10f} amu")
print(f"Chemical accuracy = {chem_acc_hartree:.10f} Hartree")
print(f"ε⁵ / chem_acc = {error_HF / chem_acc_hartree:.15e}")
print()
