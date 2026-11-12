# Yangian Y(sl(2)) Symmetry in Polyacetylene Model: Proof and Refutation

**Date**: November 12, 2026  
**Status**: REFUTED — Conjecture 4.3 is FALSE

---

## Executive Summary

**RESULT**: The polyacetylene tight-binding Hamiltonian does **NOT** possess Yangian Y(sl(2)) symmetry in the general case.

**Key Findings**:
1. ✓ Level-0 generators (standard SU(2) spin symmetry) commute with H when U=0
2. ✗ Level-1 Yangian generators do **NOT** commute with H in general
3. Failure occurs even for U=0 (non-interacting case) with open boundary conditions
4. Periodic boundary conditions with specific generator definitions can preserve symmetry only in special limits

---

## 1. The Polyacetylene Model

### Hamiltonian

```
H = H_hop + H_int

H_hop = -t Σ_{j=1}^{L-1} Σ_{σ} (c†_{j,σ} c_{j+1,σ} + h.c.)

H_int = U Σ_{j=1}^{L} n_{j,↑} n_{j,↓}
```

Where:
- L = number of carbon sites
- t = nearest-neighbor hopping (π-electron delocalization)
- U = on-site Coulomb repulsion
- c†_{j,σ}, c_{j,σ} = fermionic creation/annihilation operators
- σ ∈ {↑,↓} = electron spin

### Physical Context

This models conjugated π-electrons in polyacetylene chains (-CH=CH-CH=CH-).

---

## 2. Yangian Y(sl(2)) Structure

### Definition

A Yangian Y(g) for Lie algebra g is a Hopf algebra with generators:
```
{J^a, J^a_{(1)}, J^a_{(2)}, ...}  for a ∈ {+, -, z}
```

Level-0: Standard sl(2) ≅ su(2)
```
J^+ = Σ_j S^+_j = Σ_j c†_{j,↑} c_{j,↓}
J^- = Σ_j S^-_j = Σ_j c†_{j,↓} c_{j,↑}
J^z = Σ_j S^z_j = (1/2)Σ_j (n_{j,↑} - n_{j,↓})
```

Level-1: Yangian extension (multiple conventions)
```
Convention A (Position-weighted):
J^z_{(1)} = Σ_j j·S^z_j

Convention B (Nested commutator):
J^z_{(1)} = Σ_{j<k} [S^z_j, h_k]  where h_k is local Hamiltonian
```

### Symmetry Criterion

For Yangian symmetry:
```
[H, J^a] = 0         (Level-0)
[H, J^a_{(1)}] = 0   (Level-1)
```

---

## 3. Analytical Verification

### 3.1 Level-0 Symmetry (Standard SU(2))

**Claim**: [H, J^a] = 0 for U=0 and arbitrary boundary conditions.

**Proof**:

For J^z:
```
[H_hop, J^z] = -t Σ_j Σ_σ [c†_j c_{j+1,σ} + h.c., (1/2)Σ_k (n_{k,↑} - n_{k,↓})]

Using [c†_j c_{j+1,σ}, n_{k,σ'}] = δ_{jk} c†_j c_{j+1,σ} δ_{σσ'} - δ_{j+1,k} c†_j c_{j+1,σ} δ_{σσ'}
```

For spin-conserving hopping, each term with σ=↑ contributes equally and oppositely to σ=↓:
```
[H_hop, J^z] = 0  ✓
```

Similarly for J^±:
```
[H_hop, J^+] = 0  ✓
[H_hop, J^-] = 0  ✓
```

**Conclusion**: Level-0 SU(2) spin symmetry is exact for U=0.

When U≠0:
```
[H_int, J^z] = 0   (number operator commutes with spin)
[H_int, J^±] ≠ 0   (spin flip changes occupancy)
```

So full SU(2) only when U=0.

---

### 3.2 Level-1 Symmetry (Yangian Extension)

**Claim**: [H, J^z_{(1)}] ≠ 0 in general.

**Proof** (Convention A: position-weighted):

```
J^z_{(1)} = Σ_j j·(1/2)(n_{j,↑} - n_{j,↓})

[H_hop, J^z_{(1)}] = -t Σ_j Σ_σ [c†_j c_{j+1,σ} + h.c., Σ_k k·(1/2)(n_{k,↑} - n_{k,↓})]
```

Focus on j=1, σ=↑:
```
[c†_1 c_{2,↑}, Σ_k k·n_{k,↑}] = [c†_1 c_{2,↑}, 1·n_{1,↑} + 2·n_{2,↑} + ...]
                                = c†_1 c_{2,↑}·1 - 2·c†_1 c_{2,↑}
                                = -c†_1 c_{2,↑}
```

This gives non-zero contributions:
```
[H_hop, J^z_{(1)}] = t Σ_j (j - (j+1))(c†_j c_{j+1,↑} - c†_{j+1} c_j,↑ + ...)
                   = -t Σ_j (c†_j c_{j+1,↑} + ...) ≠ 0
```

**Result**: [H, J^z_{(1)}] ≠ 0 for position-weighted generators.

---

**Proof** (Convention B: nested commutator):

For the Hubbard chain with local Hamiltonian:
```
h_j = -t(c†_j c_{j+1} + c†_{j+1} c_j) + U n_{j,↑}n_{j,↓}
```

The proper Yangian generator involves:
```
J^z_{(1)} = Σ_{j<k} [S^z_j, h_k]
```

This construction works for integrable models where:
1. Yang-Baxter equation is satisfied
2. Transfer matrix exists
3. Local Hamiltonians form a compatible family

**For polyacetylene (single-band tight-binding)**:

The model is **NOT** Yang-Baxter integrable because:
- Scattering matrix S(u) does not satisfy Yang-Baxter equation for single-band case
- Transfer matrix formalism requires graded structure not present here
- Local h_j do not form a commuting family: [h_j, h_k] ≠ 0 for |j-k|=1

**Result**: No proper Yangian structure exists.

---

## 4. Failure Conditions

### 4.1 Why Yangian Symmetry Fails

**Fundamental obstruction**: The polyacetylene Hamiltonian is a **free fermion** problem (U=0) or **density-density interaction** (U≠0), but:

1. **Boundary breaks boost symmetry**: Open boundaries break the "boost" generator J^z_{(1)} which weights by position

2. **No R-matrix**: The model lacks an R-matrix satisfying Yang-Baxter equation
   - XXZ spin chain: Has R-matrix, is Yangian-symmetric
   - Hubbard model (1D, infinite): Has Yangian in specific limits
   - Polyacetylene: No R-matrix structure

3. **Integrability vs Symmetry**: 
   - U=0 polyacetylene is trivially solvable (free fermions) but NOT integrable in Bethe sense
   - Integrability requires infinite conserved charges
   - Yangian symmetry generates these charges

### 4.2 Conditions Where Symmetry Could Hold

**Special case: Infinite periodic chain, U=0, proper generator definition**

If we take L→∞, periodic BC, and define:
```
J^z_{(1)} = lim_{ε→0} Σ_j e^{iεj} S^z_j
```

Then in momentum space:
```
H = Σ_k ε_k c†_k c_k,  ε_k = -2t cos(k)
J^z = Σ_k S^z_k
J^z_{(1)} = Σ_k (∂/∂k) S^z_k
```

In this limit, [H, J^z_{(1)}] → 0 by translational invariance.

**But**: This is trivial — it's just momentum conservation in disguise, not genuine Yangian symmetry.

---

## 5. Correct Statement of Conjecture 4.3

### Original Conjecture (WRONG)

> "Polyacetylene tight-binding chains under single-band approximation possess Yangian Y(sl(2)) symmetry."

### Corrected Statement

> "The 1D Hubbard model with infinite chain length and specific parameter fine-tuning possesses Yangian Y(su(2)|su(2)) symmetry (charge × spin sectors), enabling Bethe ansatz solution. The simpler polyacetylene tight-binding model (U=0) does NOT have Yangian symmetry, but is exactly solvable as free fermions."

---

## 6. Exact Energy Spectrum (Alternative Method)

Since Yangian symmetry fails, we use the **correct** method for polyacetylene:

### 6.1 Free Fermion Solution (U=0)

**Bloch waves**:
```
c_k = (1/√L) Σ_j e^{ikj} c_j
```

**Dispersion relation**:
```
ε_k = -2t cos(k),  k ∈ {2πn/L : n=0,1,...,L-1}
```

**Ground state energy** (N electrons, filling lowest states):
```
E_0 = 2 Σ_{k in filled} ε_k  (factor 2 for spin degeneracy)
```

For half-filling (N=L):
```
E_0 = -(4t/π) L  (thermodynamic limit)
```

**Excited states**: Particle-hole excitations across Fermi level.

---

### 6.2 Interacting Case (U≠0)

**No Yangian symmetry** → **No Bethe ansatz** for single-band.

Instead, use:
1. **Numerical**: Density Matrix Renormalization Group (DMRG)
2. **Perturbative**: U/t ≪ 1 or U/t ≫ 1 limits
3. **Variational**: Gutzwiller wavefunctions

**Known result**: For 1D Hubbard at half-filling, U>0 opens Mott gap.

---

## 7. Minimal Correction to Conjecture

### Option A: Restrict to Integrable Models

**Replace**: "Polyacetylene chains"  
**With**: "1D Hubbard model on infinite lattice"

**Justification**: Hubbard model has proven Yangian Y(sl(2))⊕Y(sl(2)) symmetry [Uglov & Korepin, 1993].

---

### Option B: Clarify Symmetry Type

**Replace**: "Yangian Y(sl(2)) symmetry"  
**With**: "Yangian Y(gl(1|1)) graded symmetry for spinless fermions, or standard SU(2) for spinful"

**Justification**: Graded Yangians exist for supersymmetric chains, not generic tight-binding.

---

### Option C: Correct Method Statement

**Original**: "Exact energy spectrum can be derived using Yangian algebra"

**Corrected**: "Exact energy spectrum for U=0 derives from free fermion analysis; for U≠0, Bethe ansatz applies only to the 1D Hubbard model, not general polyacetylene"

---

## 8. Comparison: Models with/without Yangian Symmetry

| Model                    | Yangian Y(sl(2))? | Exact Solution Method       |
|--------------------------|-------------------|-----------------------------|
| Polyacetylene (U=0)      | ✗ No              | Free fermions (Bloch waves) |
| Polyacetylene (U≠0)      | ✗ No              | Numerical (DMRG)            |
| XXZ spin chain           | ✓ Yes             | Bethe ansatz                |
| 1D Hubbard (infinite)    | ✓ Yes*            | Bethe ansatz                |
| Richardson-Gaudin pairing| ✓ Yes (su(2))     | Richardson equations        |

*With dual Yangian structure Y(su(2))_charge ⊕ Y(su(2))_spin

---

## 9. Conclusions

### Verdict on Conjecture 4.3

**REFUTED**: Polyacetylene tight-binding chains do NOT possess Yangian Y(sl(2)) symmetry.

### Key Reasons

1. **No Yang-Baxter R-matrix** for single-band tight-binding
2. **Boundary conditions break boost generator** J^z_{(1)}
3. **Free fermion structure** is trivially solvable but not Yangian-symmetric
4. **Confusion with Hubbard model**: The conjecture conflates polyacetylene (U typically small) with 1D Hubbard (U dominant)

### Correct Exact Solution

For polyacetylene (U=0 or small):
- **Method**: Solve as free fermions in momentum space
- **Spectrum**: ε_k = -2t cos(k)
- **Ground state**: Fill Fermi sea
- **No Yangian needed**

For general conjugated chains with correlations:
- **Method**: Numerical (DMRG, exact diagonalization)
- **Analytical**: Pariser-Parr-Pople model, not Yangian

### Implications for Molecular Physics

- **Conjugated molecules**: Use Hückel theory, not Yangian
- **Strongly correlated systems**: Use Hubbard model (which HAS Yangian)
- **Exact solutions**: Reserve for integrable models (XXZ, Richardson-Gaudin)

---

## References

1. **Uglov, D. & Korepin, V. E.** (1993). "The Yangian symmetry of the Hubbard model." *Phys. Lett. A* **190**, 238.
2. **Essler, F. H. L. et al.** (2005). *The One-Dimensional Hubbard Model*. Cambridge University Press.
3. **Drinfeld, V. G.** (1985). "Hopf algebras and the quantum Yang-Baxter equation." *Sov. Math. Dokl.* **32**, 254.
4. **Gaudin, M.** (1976). "Diagonalisation d'une classe d'hamiltoniens de spin." *J. Phys.* **37**, 1087.

---

**Document Status**: Complete  
**Verification**: Analytical proof by explicit commutator calculation  
**Next Steps**: Update molecular_hamiltonian_algebraic_structure.md to correct Conjecture 4.3
