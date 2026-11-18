# CaF₂ ALD Transition State Calculation Plan

## System: Ca(thd)₂ + HF → CaF₂ + H(thd)

### 1. Electronic Structure Method Selection

**Recommended: CCSD(T) with Basis Set Extrapolation**
- **Rationale**: Ca²⁺ closed shell, weak correlation, need chemical accuracy
- **Basis sets**: aug-cc-pVTZ → aug-cc-pVQZ extrapolation
- **Expected cost**: ~500 basis functions, O(10⁶) CPU hours on modest cluster

**Alternative: DFT Benchmark**
- **Functionals to test**: ωB97X-D, M06-2X, PBE0-D3(BJ)
- **Validation**: Against CCSD(T) single points at critical geometries
- **Pragmatic choice**: If DFT error < 2 kcal/mol, use for full PES

### 2. Reaction Coordinate Definition

**Primary coordinate**: Ca—F distance (Ca(thd)₂ ⋯ HF complex)
- **Reactant**: r(Ca—F) ≈ 4.5 Å (van der Waals complex)
- **TS**: r(Ca—F) ≈ 2.3 Å (estimate from ionic Ca—F bond 2.0 Å)
- **Product**: r(Ca—F) = 2.0 Å (CaF bond)

**Secondary coordinate**: H—F distance (HF dissociation)
- **Reactant**: r(H—F) = 0.92 Å (HF equilibrium)
- **TS**: r(H—F) ≈ 1.1 Å (elongated)
- **Product**: r(H—F) → ∞ (dissociated, H transfers to thd)

**Coupling**: Synchronous Ca—F formation and H—F cleavage (concerted mechanism expected)

### 3. Transition State Search Protocol

#### Step 1: Relaxed Surface Scan
```
FOR r_CaF = 2.0 to 4.5 Å (step 0.1 Å):
  Constrain r(Ca—F) = r_CaF
  Optimize all other coordinates
  Record E, ∇E, Hessian eigenvalues
```

#### Step 2: TS Guess from Maximum
- Identify maximum energy along scan
- Check for one imaginary frequency (ω < 0)
- If multiple imaginary modes → reduce scan step size

#### Step 3: Eigenvector-Following TS Optimization
- Algorithm: Partitioned Rational Function Optimization (P-RFO)
- Convergence: ΔE < 1e-6 Ha, |∇E| < 3e-4 Ha/Bohr
- Hessian: Analytical (DFT) or numerical (CCSD(T))

#### Step 4: Intrinsic Reaction Coordinate (IRC)
- Follow gradient descent in mass-weighted coordinates
- Verify IRC connects to reactant and product minima
- Path smoothness → no competing mechanisms

### 4. Nonadiabatic Coupling Calculation

**Electronic states to probe**:
1. Ground state (S₀): Ca(thd)₂(¹A₁) + HF(¹Σ⁺) → CaF₂(¹Σ_g⁺)
2. First excited state (S₁): Ca(thd)₂(³T₂) + HF (triplet surface)
3. Second excited state (S₂): Ca 3d → 4s excitation (~4.5 eV)

**Calculation protocol**:
```
At TS geometry X_TS:
  1. Calculate S₀, S₁, S₂ energies (TD-DFT or EOM-CCSD)
  2. Verify ΔE(S₁-S₀) > 3 eV (large gap → BO valid)
  3. Compute derivative coupling: d_01 = <ψ₀|∇_X ψ₁>
  4. Check criterion: ε²‖d_01‖ < 0.1
```

**Software**: Gaussian 16 (TD-DFT), Q-Chem (EOM-CCSD), ORCA (SOC)

### 5. Expected Results

**Energy barrier** (based on literature for Ca + HF):
- ΔE‡ ≈ 5–15 kcal/mol (highly exothermic, early TS)
- Negative activation energy possible (barrierless at low T)

**Geometry changes at TS**:
- Ca coordination: square planar (thd) → distorted trigonal bipyramid
- F approach: nearly linear Ca⋯F—H (A′ mechanism)
- thd ligands: one thd begins dissociation (trans to F)

**Nonadiabatic effects**:
- Gap at TS: Δ(S₁–S₀) ≈ 4.5 eV (from BO analysis)
- Coupling: ‖d₀₁‖ ≈ 0.01 a.u. (estimate)
- BO parameter: ε²‖d₀₁‖ ≈ 9e-5 ≪ 0.1 → **BO valid**

### 6. Computational Cost Estimates

| Method | Basis | Time/point | TS search | IRC | Total |
|--------|-------|------------|-----------|-----|-------|
| DFT/ωB97X-D | aVTZ | 30 min | 10 iter | 40 pts | ~30 CPU-hr |
| CCSD(T) | aVTZ | 48 hr | SP only | SP only | ~500 CPU-hr |
| TD-DFT (states) | aVTZ | 2 hr | 1 pt | - | ~2 CPU-hr |

**Recommendation**: DFT for optimization, CCSD(T) single points for benchmark

### 7. Validation Checklist

- [ ] TS has exactly one imaginary frequency
- [ ] Imaginary mode corresponds to Ca—F stretch + H—F stretch
- [ ] IRC connects to correct reactant/product
- [ ] Electronic gap Δ(S₁–S₀) > 3 eV at TS
- [ ] BO error O(ε⁵) ≈ 2.5e-13 < 1e-6 (chemical accuracy)
- [ ] No conical intersections along IRC (gap never closes)
- [ ] Vibrational frequencies all real at minima

### 8. Software and Input Examples

**Gaussian 16 TS Search**:
```
%nprocshared=16
%mem=64GB
%chk=caf2_ts.chk
# wB97XD/aug-cc-pVTZ Opt=(TS, CalcFC, NoEigenTest) Freq

Ca(thd)2 + HF TS search

0 1
Ca  0.0  0.0  0.0
F   2.3  0.0  0.0
H   2.5  0.0  0.8
[... full geometry ...]
```

**Q-Chem Derivative Coupling**:
```
$molecule
0 1
[TS geometry]
$end

$rem
method = eom-ccsd
basis = aug-cc-pvtz
ee_states = [1,0,0,0]  ! 1 excited state
calc_nac = true
$end
```

### 9. Contingency: If BO Breakdown Detected

**Symptoms**:
- Gap Δ < 0.5 eV at any IRC point
- Irregular IRC (kinks, bifurcations)
- ε²‖d₀₁‖ > 0.1

**Remediation**:
1. Switch to surface hopping (SHARC, JADE)
2. Compute spin-orbit coupling (Ca 3d states)
3. Multi-reference methods (CASSCF, MRCI)
4. Report to experimentalists: temperature/pressure dependence

### 10. Deliverables

1. **Transition state geometry** (XYZ, optimized)
2. **Energy profile**: E(IRC) from reactant to product
3. **Barrier height**: ΔE‡ ± σ (kcal/mol, with uncertainty)
4. **Nonadiabatic coupling**: ‖d₀₁‖ and gap Δ at TS
5. **BO validity certificate**: ε⁵ error bound vs chemical accuracy
6. **Vibrational frequencies**: ν_TS (imaginary), ν_R, ν_P (real)
7. **Report**: 5-page summary with figures (PES, IRC, MO diagrams)
