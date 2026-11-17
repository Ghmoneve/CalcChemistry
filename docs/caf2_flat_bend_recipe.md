# Converging CaF2 on its flat bend — a numerical recipe

Scope: why the existing damped-BFGS loop stalls on CaF2, what to change, whether the
basis set / grid / symmetry constraints are responsible for the bent minimum, and the
minimum defensible treatment for a reaction-energy paper. Everything numerical below was
measured in this repository with PySCF 2.15.0 (`/usr/bin/python3`, PBE/def2-svp unless
stated) using `experiments/caf2/flat_mode_study.py`; raw numbers are in
`experiments/caf2/flat_mode.json`. A working optimiser is `experiments/caf2/caf2_relax.py`.

---

## 0. Bottom line

1. **The optimiser, not the physics, is the problem.** The bend is an ordinary
   low-frequency mode (ω ≈ 103 cm⁻¹, curvature 4.1e-4 Ha bohr⁻² amu⁻¹, only ~35× softer
   than the Ca–F stretch). BFGS handles κ≈35 in a dozen iterations. Your loop fails for
   three specific, fixable reasons (§1).
2. **Your 1e-5 convergence target is unreachable in Cartesian coordinates** with
   atom-centred grids, independent of the optimiser: the analytic gradient carries a
   spurious net force of 1.9e-4 Ha Å⁻¹. Converge on the projected (internal) gradient, or
   set `grad.grid_response = True` (§2).
3. **The bent minimum is real, not a basis-set artifact** — CCSD(T) and CASPT2 also give
   a bent, quasilinear CaF2. But the barrier is basis-sensitive by a factor of 6
   (2.20 → 0.35 kcal mol⁻¹), so def2-svp overstates it (§3).
4. **The grid matters for the Hessian, not for the geometry.** The bend eigenvalue at the
   linear geometry flips sign at levels 0–3 (+414 cm⁻¹ at level 2!); levels ≥4 give
   −103.5, −102.4, −102.8 cm⁻¹, and level 6 agrees with an unpruned (150,770) reference to
   0.01 cm⁻¹ (§4).
5. **Do not constrain to C∞v** unless you are prepared to call the result a saddle. For a
   reaction energy, report the bent minimum and quote the linear reference as the
   systematic uncertainty (§5, §6).

---

## 1. Why the damped BFGS stalls, and what to change

### 1.1 Measured replay

Replaying `reoptimize.py:relax()` exactly (H₀ = 0.3·I, cap 0.08 Å, no line search,
`sy > 1e-12` guard), started from θ = 140°, r = 1.9925 Å:

| iter | E / Ha | ‖g‖ / Ha Å⁻¹ | θ / ° | eigenvalues of B_k (internal space) |
|---:|---|---:|---:|---|
| 0 | −876.76991703 | 1.76e−2 | 140.00 | 0.30, 0.30, 0.30 |
| 4 | −876.77021484 | 7.13e−3 | 137.12 | 0.30, 0.67, 13.6 |
| 6 | −876.77027571 | 1.33e−3 | 134.19 | 0.30, 0.89, 24.0 |
| 8 | −876.77027565 | **1.10e−4** | 134.25 | 0.30, 1.05, 22.3 |
| 14 | −876.77026939 | 2.40e−3 | 133.59 | 0.29, 1.88, 16.3 |
| 20 | −876.77012907 | 1.35e−2 | 131.07 | 0.29, 9.3, 125 |
| 30 | −876.76990655 | 2.22e−2 | 129.52 | 0.29, 8.1, 325 |
| 35 | −876.76982978 | **2.46e−2** | 129.15 | 0.29, 8.0, 504 |

True mass-weighted curvature: bend 4.06e−4, stretches 1.07e−2 and 1.42e−2
(ω = 103.5, 531.9, 612.4 cm⁻¹). This reproduces your symptom — ≈2.5e−2 Ha Å⁻¹ and a
drifting angle — after ~35 iterations; the exact value 3.0e−2 at 128.2° presumably came
from a slightly different start/round count. The *mechanism* is unambiguous:

**(a) H₀ = 0.3 is 740× too stiff along the bend, and the bend never learns.**
The smallest eigenvalue of B_k stays at exactly 0.3000 for 30+ iterations: no secant pair
with bend character is ever generated, so the rank-2 update never touches that direction.
The bend is therefore minimised by steepest descent with step 0.3·g, i.e. with an
effective relaxation factor 1 − 0.3·4.06e−4 = 0.99988 per iteration — ~10⁴ iterations to
converge. This is the "tiny, poorly-scaled steps" you observed, and it is a *scaling*
failure, not a flatness failure.

**(b) The step cap destroys the secant relation, and B then diverges.**
Capping ‖−B_k g‖ at 0.08 Å rescales the step, so s ≠ −B_k g and the pair (s, y) no longer
estimates curvature along the step. Once the geometry reaches the noise floor, s and y are
both ~1e−7, so s·y ~ 1e−14 passes the `sy > 1e-12` test, ρ = 1/s·y explodes, and B_k
inflates without bound: 0.30 → 22 → 125 → 504 (a factor 3.5e4 too *soft*). The next steps
overshoot, the energy *rises* by 4e−4 Ha, and the angle walks from 134° to 129°. This is
the "stall": a runaway, not a plateau.

**(c) No line search means nothing enforces descent.** `relax()` accepts every step. Any
quasi-Newton method needs either a Wolfe/Armijo line search (Nocedal & Wright, *Numerical
Optimization*, 2nd ed., Ch. 3) or a trust region (Ch. 4); a hard cap is neither.

### 1.2 What to change, in order of payoff

**Fix 1 — seed the inverse Hessian with the analytic Hessian (largest single win).**
Compute the RKS Hessian once at the start, project out translations/rotations, and use
H₀⁻¹ = Σ_k v_k v_kᵀ/λ_k over the 3N−6 physical modes, with the null space set to exactly
zero. Measured (PBE/def2-svp, grid 4): θ 140.00 → 134.25°, ‖g‖ = 1.1e−4, **5 iterations**;
with `grid_response=True` it reaches ‖g‖ = 5.6e−6 in 8. Cost: one 8 s Hessian. This is the
standard "exact initial Hessian" option offered by production optimisers (Schlegel 1982;
Baker/Kessi/Delley 1996).

**Fix 2 — a relative curvature guard on the update.**
Replace `if sy > 1e-12:` with

```python
if sy > 1e-8 * np.linalg.norm(s) * np.linalg.norm(y):
```

With this single change the same H₀ = 0.3·I loop stops diverging and reaches ‖g‖ = 1.09e−4
(‖Pg‖ = 9.4e−6) at iteration 8, B_max stabilised at 22 instead of 500. The same idea is
exposed as `min_curvature` in `scipy.optimize.BFGS`; the classical alternative is Powell
damping of the update.

**Fix 3 — a backtracking line search (or a trust radius accepted on ΔE).**
Halve the step until E decreases. With Fixes 1–3 the whole relaxation is monotone,
6 iterations, ‖Pg‖ = 7.7e−7 (table below, §6). Cost: 1–2 extra SCFs per iteration.

### 1.3 Trust region vs line search vs Newton/RFO — honest comparison

* **Line search (with Wolfe conditions) + scaled BFGS**: cheapest, sufficient here once
  B₀ is right. Without B₀ it cannot fix a direction that is wrong by 740×; a line search
  only rescales the step, it does not rotate it.
* **Trust region**: better when the quadratic model is unreliable, because the step is
  bounded by a radius and accepted on the ratio of predicted to actual ΔE. This is what
  makes RFO work (Banerjee/Adams/Simons/Shepard 1985; ASE's `ase.optimize.RFO`, installed
  here as ASE 3.29, defaults: `damping = 1.0 Å⁻¹`, `maxstep = 0.2 Å`, `alpha = 70 eV Å⁻²`).
* **Plain Newton with the analytic Hessian and a *global Cartesian* trust radius is not
  enough**: measured, the first step overshot by 0.0014 Ha (0.9 kcal mol⁻¹) because the
  bend component of the step is unbounded (the harmonic model is useless where λ ≈ 4e−4)
  while the stretch component is stiff; a single radius on the Cartesian 9-vector then
  spends the whole budget on the bends and truncates the stretches. If you use Newton,
  either work in internal coordinates (r₁, r₂, θ) or use per-mode (eigenvector-following)
  scaling. geomeTRIC (Wang & Song 2016) does exactly this with translation–rotation
  internal coordinates; `pip install geometric` is the one-line fix, and PySCF's
  `pyscf.geomopt.geometric_solver` will then drive it with analytic gradients.
* **RFO** is the safest general-purpose choice, but with Fix 3 you do not need it for a
  triatomic.

**Do not** simply lower the step cap (the usual reflex). A smaller cap makes (b) worse:
every step is then cap-dominated, no curvature information is accumulated, and the run
becomes pure steepest descent with a 740×-wrong metric.

---

## 2. What actually limits the convergence: the gradient is not the derivative of the energy

PySCF's RKS gradient defaults to `grid_response = False`
(`pyscf/grad/rks.py:750` → `grad_rks_Gradients_grid_response`, default `False`), and the
RKS Hessian refuses grid response outright
(`pyscf/hessian/rks.py:47`: `NotImplementedError('CPU RKS Hessian grid response is not implemented')`).
With an atom-centred quadrature the analytic gradient is then *not* the derivative of the
energy being minimised.

At the converged bent geometry (θ = 134.254°, r = 1.97039 Å, PBE/def2-svp, grid 4):

| quantity | value |
|---|---|
| ‖g‖ (as used by your loop) | **1.10e−4** Ha Å⁻¹ |
| ‖P g‖, projected on the internal space | **5.2e−6** Ha Å⁻¹ |
| ‖Σ_atoms F‖ (net force, must be 0) | **1.9e−4** Ha Å⁻¹ |
| ‖g_FD − g_analytic‖, same grid | 8.6e−4 |
| same, grid level 6 | 7.9e−5 |
| same, `grid_response=True`, grid 4 | **1.3e−5** |
| net force with `grid_response=True` | 1e−14 |

So at the geometry your optimiser is sitting on, the molecule is already stationary to
5e−6 in every internal direction; the entire 1.1e−4 is a spurious rigid-body force from
the omitted grid-response term. **No optimiser can drive that to 1e−5**, which is why the
seeded run froze at 1.099e−4 for 30 iterations while E stayed constant to 1e−9 Ha.

Actionable:

```python
g = mf.nuc_grad_method()
g.grid_response = True          # restores translational invariance exactly
grad = g.kernel()
```

and judge convergence on the projected gradient:

```python
P, _ = internal_projector(mol, coords)      # rank 3N-6 (3N-5 if linear)
gnorm = np.linalg.norm(P @ grad)
```

Note the flip side: `grid_response=True` gives a *different* (correct) surface derivative,
and its stationary point sits at θ = 134.67°, r = 1.9715 Å — 0.4° and 0.001 Å away from the
`grid_response=False` point, essentially isoenergetic with it (ΔE = 1.5e−9 Ha). The 0.4° is
the size of the inconsistency; it is only visible because the bend is so soft.

---

## 3. Basis set: are diffuse functions needed, and is a bent minimum a basis artifact?

### 3.1 Measured basis dependence (relaxed r at every angle)

Energy of the linear geometry above the bent minimum, kcal mol⁻¹:

| level | barrier to linearity | θ_min / ° | r(Ca–F) at θ_min / Å |
|---|---:|---:|---:|
| PBE/def2-svp | **2.20** | ~134 | 1.9716 |
| PBE/def2-tzvp | 1.14 | 140 | 2.0006 |
| PBE/def2-svpd (diffuse) | **0.50** | 144.0 | 2.0074 |
| PBE/ma-def2-svp (diffuse) | 0.65 | 140 | 2.0070 |
| PBE/def2-tzvpd (diffuse) | **0.88** | 140 | 2.0055 |
| B3LYP/def2-svp | 1.53 | 140 | 1.9765 |
| B3LYP/def2-tzvpd | **0.35** | 150 | 2.0152 |

(θ_min from a 10° relaxed scan, except svp and svpd which were fully relaxed with
`caf2_relax.py`; the scan brackets are ±5°.)

Harmonic bend frequency at the linear geometry (analytic Hessian, projected):

| level | bend / cm⁻¹ |
|---|---:|
| PBE/def2-svp | −103.9 |
| PBE/def2-svpd | −58.5 |
| PBE/ma-def2-svp | −63.6 |
| PBE/def2-tzvp | −127.2 |

Two clear effects, both in the same direction:

* **def2-svp over-binds the bend** (barrier 2.20 vs 0.50–0.88 with diffuse functions, i.e.
  2.5–4×) and **contracts Ca–F by 0.036 Å** (1.969 vs 2.005 Å). This is the expected
  behaviour of a basis without diffuse functions on a fluoride: F⁻ is too compact.
* **The bend survives every basis.** Diffuse functions soften it substantially but never
  make the linear geometry the minimum: the imaginary bend stays between −59 and −127 cm⁻¹.

A caution when reading these two tables together: the harmonic bend frequency at the
*linear* geometry is not a reliable proxy for the barrier, because the well is strongly
anharmonic away from the top. def2-tzvp has the largest |ω_bend| (−127 cm⁻¹) but a smaller
barrier than def2-svp (1.14 vs 2.20 kcal). Use the barrier, not the frequency, when
comparing levels.

### 3.2 Is a bent minimum in an ionic triatomic a sign of BSSE or missing diffuse functions?

**Not BSSE.** Basis-set superposition error is defined for the *interaction* of two or more
fragments computed in a common basis (Boys & Bernardi 1970); it has no meaning for the
geometry of a single molecule. The correct concern is plain basis-set incompleteness, and
the answer for CaF2 is: it changes the barrier by 2–4× but not the topology.

**And the "conventional linear CaF2" premise is out of date.** Gas-phase CaF2 is bent:

* Electric-deflection measurements placed CaF2 among the bent alkaline-earth dihalides
  (Klemperer et al., as summarised in Kaupp et al. 1991).
* Prasad, Wittmaack & Donald (2018) state plainly that "the Ba dihalides, SrF2, SrCl2 and
  CaF2 are bent", and that CaF2 and SrCl2 have surfaces "so flat that they are described as
  floppy or quasi-linear".
* The most recent high-level study (Romeu et al., JACS 2026, CCSD(T)/awQ(-PP) and
  CASPT2/awQ(-PP), all-electron and PP) finds CaF2 bent and adds the key sentence: "the
  bending harmonic frequencies of CaF2 and SrCl2 are comparable to the inversion barriers,
  so the molecules are predicted to be quasilinear at most temperatures and highly
  fluxional."
* The NIST CCCBDB entry lists CaF2 as C2v with fundamentals 490, 120, 575 cm⁻¹
  (Gurvich) — i.e. a 120 cm⁻¹ bend, essentially the mode you are computing (103 cm⁻¹
  harmonic at PBE/def2-svp).

So a bent PBE minimum near 135° with a 0.5–2 kcal mol⁻¹ barrier is **not** a numerical
artifact and should not be described as one. It is the same physics as the CCSD(T)
result, exaggerated by the functional and by the missing diffuse functions. Where DFT
does err systematically is the angle and the barrier magnitude (the JACS 2026 paper notes
that DFT-type methods "tend to predict smaller bond angles ... and more bent minimum
energy cases"): PBE/def2-svp gives ~134°, the correlated methods give a shallower,
wider well.

One caveat worth stating in the paper: because the barrier is comparable to the bend
frequency, the harmonic-oscillator/rigid-rotor thermochemistry that would normally be
applied to this molecule is not valid; a large-amplitude (2-D bend–rotation) treatment is
what the spectroscopy literature uses. I did not implement one, so any ZPE/entropy
correction you add for CaF2 from a normal-mode analysis should be flagged as an
approximation.

---

## 4. Grid convergence

PBE/def2-svp, energies at fixed r = 1.9925 Å, Hessian at θ = 180° (projected):

| level | (nrad, nang) F / Ca | ΔE(140°−180°) / kcal | bend at 180° / cm⁻¹ | largest of the 5 null modes / cm⁻¹ |
|---:|---|---:|---:|---:|
| 0 | (15,86) / (30,110) | −1.419 | −1350 | 4860 |
| 1 | (40,194) / (60,194) | −2.316 | **+184.7** | n.m. |
| 2 | (60,302) / (75,302) | −2.044 | **+414.3** | n.m. |
| 3 (PySCF default) | (75,302) / (90,434) | −2.040 | −172.9 | 142 |
| 4 | (90,590) / (105,590) | −2.056 | −103.5 | 16.7 |
| 5 | (105,770) / (120,770) | −2.058 | −102.4 | n.m. |
| 6 | (120,974) / (135,974) | −2.057 | −102.8 | 7.9 |
| unpruned (150,770) | | — | −102.79 | — |

(Values from `pyscf.dft.gen_grid.RAD_GRIDS` / `ANG_ORDER` / `LEBEDEV_ORDER`; the docstring
of `Grids` states that level 3 is (75,302) for second-row and (80–105,434) for heavier
atoms. "n.m." = not measured.)

**Verdict.** For the *energy* of the bend, anything from level 3 up is fine: it changes by
0.02 kcal between levels 3 and 5, and only levels 0–1 are visibly wrong (0.6 kcal). For the
*Hessian* of this mode the grid is decisive: levels 0–3 give −1350, +185, +414 and
−173 cm⁻¹, i.e. the sign of the curvature at the linear geometry flips; from level 4 it is
stable and converged by level 5 to better than 1 cm⁻¹. **Use `grids.level >= 4` for the
bend Hessian (level 5 if you want sub-cm⁻¹ stability).** Their current `GRID = 4` is
acceptable for energies and marginal-but-usable for frequencies; `GRID = 3` would not be.

Two further grid facts that matter here:

* **The rotational-symmetry breaking shrinks with the grid but is large at the default.**
  At θ = 180° the five null mass-weighted eigenvalues (3 translations + 2 rotations; the
  rotation about the axis is identically null) are spread up to 142 cm⁻¹ at level 3,
  17 cm⁻¹ at level 4 and 8 cm⁻¹ at level 6. In Ha bohr⁻² (un-mass-weighted) that is
  ~1.5e−2 at level 3, 2e−4 at level 4 and 5e−5 at level 6 — so the "~1e−2 Ha/bohr²
  spurious rotation" you observed is the *default grid* value, not a universal one.
  At level 3 the contamination (142 cm⁻¹) is larger than the bend itself (103 cm⁻¹), so
  projection is mandatory there; at level 4+ it is 2–6 % of the bend and projection is
  merely good practice. A fixed |λ| threshold cannot replace the projection in either
  case, because the bend eigenvalue is the smallest physical eigenvalue.
* **The grid sets the achievable gradient.** ‖g_FD − g_analytic‖ = 8.6e−4 (level 4) and
  7.9e−5 (level 6). Since the bend force constant is only 4e−4 in mass-weighted units, a
  1e−4 gradient error corresponds to a ~0.4° uncertainty in the bend angle. This is the
  quantitative reason a "3e−2" gradient says nothing about the bend being unconverged.

---

## 5. Should you constrain to C∞v?

**No — not as the primary geometry.** At PBE/def2-svp the linear point has a doubly
degenerate imaginary bend (−103.9 cm⁻¹); constraining C∞v and optimising r puts you on a
saddle, and the "bend frequency" you report is the imaginary eigenvalue of that saddle. It
is not wrong to *do* the calculation, but calling the result "the relaxed geometry" is.

The honest trade-offs:

* **Constraining is defensible as a defined reference.** "Linear-reference energy at fixed
  C∞v, with the perpendicular bend treated separately" is a legitimate convention, and it
  is what the older literature effectively does when it assumes VSEPR linearity. It
  becomes misleading only if the constraint is presented as a minimum.
* **It costs 0.35–2.2 kcal mol⁻¹** depending on level (§3.1). That is between 0.35× and 2×
  chemical accuracy, on a reaction where CaF2 carries coefficient +1 — so it is *not* a
  negligible convention. On your reaction, replacing the linear product with the bent
  minimum shifts ΔE by −2.3 kcal mol⁻¹ at PBE/def2-svp (E_lin = −876.766619 at r = 1.9925,
  E_bent = −876.770245), but by only −0.5 kcal mol⁻¹ with diffuse functions on F.
* **Starting a Hessian-seeded optimiser at exactly 180° will "converge" immediately** at
  the saddle, because the bend gradient is zero by symmetry. Measured with the recipe
  script: from θ = 180° it reaches ‖P g‖ = 5.8e−7 in 3 iterations, r = 1.99196 Å,
  frequencies −103.2 / −103.2 / 517.4 / 637.9 cm⁻¹. Every gradient criterion is satisfied
  and the geometry is a saddle. This is the trap that produced the original saddle-point
  geometries in `reaction2.json` (CaF2 there is exactly linear, r = 1.9925 Å, and is a
  saddle with two imaginary bend modes). Always displace along the lowest Hessian
  eigenvector before re-optimising, or start bent.
* **A genuine standard practice exists, but it is different from constraining:** report
  the molecule as quasilinear, give the inversion barrier and the bend frequency, and
  state that the geometry is fluxional (Prasad et al. 2018; Romeu et al. 2026). If you
  need a single number for a thermodynamic cycle, the linear point is the *upper* bound
  and the bent minimum the lower bound, and the difference is your uncertainty.

---

## 6. Minimum defensible treatment for a reaction energy

If the goal is a defensible ΔE for `Ca(famd)2 + 2 HF → CaF2 + 2 formamidine`, and not a
spectroscopic study of CaF2:

1. **Basis.** Put diffuse functions on F: `def2-svpd` is the cheap fix (barrier 0.50 vs
   2.20 kcal), `def2-tzvpd` if you can afford it (0.88 kcal, and the Ca–F distance is
   0.036 Å longer). Keep the same basis for the other species so the reaction energy stays
   balanced. This is the single change that most affects the CaF2 term.
2. **Optimise without symmetry constraints, from a bent start, with the analytic Hessian
   as the seed and a line search.** `experiments/caf2/caf2_relax.py` does this; measured,
   6 iterations, monotone:

   | iter | E / Ha | ‖P g‖ / Ha Å⁻¹ | θ / ° |
   |---:|---|---:|---:|
   | 0 | −876.769915781 | 1.75e−2 | 140.00 |
   | 1 | −876.770206567 | 7.39e−3 | 136.55 |
   | 2 | −876.770237947 | 2.52e−3 | 134.16 |
   | 3 | −876.770244959 | 4.37e−4 | 134.92 |
   | 4 | −876.770245070 | 1.01e−4 | 134.95 |
   | 5 | −876.770245072 | 6.37e−5 | 134.95 |
   | 6 | −876.770245073 | **7.7e−7** | 134.95 |

   Final: r(Ca–F) = 1.9716 Å, θ = 134.95°, frequencies 108.4 / 560.3 / 638.1 cm⁻¹ (no
   imaginary mode). With diffuse functions on F the same script gives
   PBE/def2-svpd: θ = 144.0°, r = 2.0074 Å, E = −876.803599 Ha,
   frequencies 66.4 / 512.6 / 594.5 cm⁻¹ — a much shallower well (barrier 0.50 kcal)
   and a bend that is only 66 cm⁻¹.
3. **Verify.** Frequencies with translations/rotations projected out, counting 3N − n_proj
   eigenvalues of largest |λ| — your existing `freqs` code already does this correctly.
4. **Report the spread, not a single angle.** State ΔE(linear → bent) at your production
   level and note that it is the systematic uncertainty of the CaF2 term (0.5–2.3 kcal
   mol⁻¹ depending on basis). That is smaller than your existing 3.0 kcal mol⁻¹ spread
   across the four method levels, so it does not overturn the paper's conclusion, but it
   must replace the claim that 128° is a verified minimum.
5. **Do not report ZPE/thermal corrections for CaF2 from the harmonic analysis without a
   caveat** — the barrier is comparable to the bend frequency (§3.2).

Is there precedent for "fix linear and acknowledge it"? Yes, in the sense that essentially
all the older CaF2 work assumes linearity from VSEPR, and the modern literature
(Prasad et al. 2018; Romeu et al. 2026) reports both the linear point and the shallow
minimum because the two are within the same kilocalorie. So a paper may legitimately say
"CaF2 was treated as linear (C∞v); the computed inversion barrier is 0.5–2 kcal mol⁻¹, which
bounds the effect on ΔE". What is *not* defensible is presenting a bent geometry at 128.2°
with ‖g‖ = 3e−2 as a converged minimum, or presenting a linear geometry as the minimum
when the Hessian at that geometry has an imaginary bend.

---

## 7. What I am not sure about

* **The exact CCSD(T) angle and barrier for CaF2.** The JACS 2026 paper is paywalled; I
  have only the abstract (bent, barrier comparable to the bend frequency). I did not
  extract a number, so I do not quote one.
* **The experimental r(Ca–F).** CCCBDB has no geometry for CaF2 in the fields I could
  read; I therefore compare basis sets to each other, not to experiment.
* **`grid_response=True` robustness.** It works and is FD-consistent here, but
  `pyscf/grad/rks.py` contains a commented-out warning that grid pruning affects gradient
  accuracy; I did not test pruned vs unpruned converged geometries systematically
  (the numbers above show pruning is irrelevant at level 4, 1.87e−4 vs 1.91e−4 net force).
* **The "128.2° / 3.0e−2" numbers from your run.** My replay gives 129.15° / 2.46e−2 at
  iteration 35; the difference is presumably the start geometry or the round count. The
  failure mode is identical.
* **RRHO invalidity for quasilinear molecules** is standard knowledge in spectroscopy but
  I have not cited a specific source for it here.

---

## 8. References

1. M. Kaupp, P. v. R. Schleyer, H. Stoll, H. Preuss, "The question of bending of the
   alkaline earth dihalides MX2 … An ab initio pseudopotential study",
   *J. Am. Chem. Soc.* **113**, 6012–6020 (1991). doi:10.1021/ja00016a013
2. S. Prasad, B. K. Wittmaack, K. J. Donald, "Bending Ternary Dihalides: A Single
   Functional Form for Linearization Energies", *J. Phys. Chem. A* **122**, 9065–9072
   (2018). doi:10.1021/acs.jpca.8b09068
3. J. G. F. Romeu et al., "When Core Orbitals Act as Valence Orbitals: Linear and Bent
   Geometries of MX2 for X a Halogen and M an Alkaline Earth", *J. Am. Chem. Soc.* **148**,
   828–836 (2026). doi:10.1021/jacs.5c16443
4. NIST Computational Chemistry Comparison and Benchmark Database, CaF2 entry
   (CAS 7789-75-5), vibrational fundamentals from Gurvich.
5. S. F. Boys, F. Bernardi, "The calculation of small molecular interactions by the
   differences of separate total energies", *Mol. Phys.* **19**, 553 (1970).
   doi:10.1080/00268977000101561
6. F. Weigend, R. Ahlrichs, "Balanced basis sets of split valence, triple zeta valence and
   quadruple zeta valence quality for H to Rn", *Phys. Chem. Chem. Phys.* **7**, 3297
   (2005). doi:10.1039/b508541a
7. D. Rappoport, F. Furche, "Property-optimized Gaussian basis sets for molecular response
   calculations", *J. Chem. Phys.* **133**, 134105 (2010). doi:10.1063/1.3484283
   (def2-SVPD / def2-TZVPD, the diffuse sets used above)
8. Diffuse functions for anions: T. Clark, J. Chandrasekhar, G. W. Spitznagel,
   P. v. R. Schleyer, "Efficient diffuse function-augmented basis sets for anion
   calculations. III. The 3-21+G basis set for first-row elements, Li–F",
   *J. Comput. Chem.* **4**, 294–301 (1983); and R. A. Kendall, T. H. Dunning Jr.,
   R. J. Harrison, "Electron affinities of the first-row atoms revisited. Systematic basis
   sets and wave functions", *J. Chem. Phys.* **96**, 6796 (1992).
9. J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed., Springer (2006):
   Ch. 3 (line search, Wolfe conditions), Ch. 4 (trust region), Ch. 6 (BFGS, initial
   scaling).
10. D. F. Shanno, K.-H. Phua, "Matrix conditioning and nonlinear optimization",
    *Math. Program.* **14**, 149–160 (1978). (initial-Hessian scaling γ = s·y/‖y‖²)
11. H. B. Schlegel, "Optimization of equilibrium geometries and transition structures",
    *J. Comput. Chem.* **3**, 214–218 (1982). doi:10.1002/jcc.540030212
12. A. Banerjee, N. Adams, J. Simons, R. Shepard, "Search for stationary points on
    surfaces", *J. Phys. Chem.* **89**, 52–57 (1985). doi:10.1021/j100247a015 (RFO)
13. P. Pulay, G. Fogarasi, "Geometry optimization in redundant internal coordinates",
    *J. Chem. Phys.* **96**, 2856 (1992). doi:10.1063/1.462844
14. J. Baker, A. Kessi, B. Delley, "The generation and use of delocalized internal
    coordinates in geometry optimization", *J. Chem. Phys.* **105**, 192 (1996).
    doi:10.1063/1.471864
15. L.-P. Wang, C. Song, "Geometry optimization made simple with translation and rotation
    coordinates", *J. Chem. Phys.* **144**, 214108 (2016). doi:10.1063/1.4952956
    (geomeTRIC; also documents why Cartesian coordinates are a poor choice for flat modes)
16. Q. Sun et al., "Recent developments in the PySCF program package",
    *J. Chem. Phys.* **153**, 024109 (2020). doi:10.1063/5.0006074
17. PySCF source, this installation: `pyscf/grad/rks.py:750` (`grid_response` default
    `False`), `pyscf/hessian/rks.py:47` (grid response not implemented for RKS Hessians),
    `pyscf/dft/gen_grid.py` (`Grids` docstring and `RAD_GRIDS`/`ANG_ORDER` tables).
