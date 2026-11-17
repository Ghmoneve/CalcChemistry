# CaF₂ bending potential: what level of theory is actually required?

Assessment for the model reaction Ca(formamidinate)₂ + 2 HF → CaF₂ + 2 formamidine.
Date: 2026-11-17. Literature and methodology only — no calculations were run.

---

## 0. Executive answer

- **CaF₂ is not linear.** It is bent/quasilinear with a *small* barrier to linearity.
  The "conventionally linear" premise in `sections/caf2_worked_example.tex` is factually wrong
  and should be corrected.
- The observed imaginary bend at D∞h is therefore **the correct physical result, not a symptom
  of a DFT artefact.** Every *correlated* treatment on record (MP2, QCISD, B3LYP, CCSD(T), MRCI,
  CASPT2) finds CaF₂ bent or quasilinear; HF is the only level ever reported as giving a strictly
  linear geometry, and Seijo et al. still classify CaF₂ as quasilinear at that level. The sentence
  "the imaginary bend survives at all four levels, so it is not a basis artefact" reasons in the
  wrong direction: survival of the imaginary mode is expected and carries no diagnostic
  information. What would have been diagnostic is agreement in the *depth* of the well.
- What *is* wrong is the **magnitude**: our 1.9–2.1 kcal mol⁻¹ well depth and 128° minimum are
  roughly 2× too deep and ~15–35° too bent against every converged number in the literature,
  which span 0.01–0.05 eV (0.2–1.2 kcal mol⁻¹) at 142–168°.
- The dominant error is **the Ca d/f basis**, not the functional. `def2-SVP`/`def2-TZVP` for Ca
  carry a 2- or 3-function d manifold with the tightest primitive at 5.5 and **no f functions**;
  the converged literature requires 4–6 largely uncontracted d functions with exponents reaching
  ~7–10.5, plus f. This is precisely the manifold that Wright, Lee & Dyke showed controls the
  bent/linear energy difference.
- **Recommendation: (c) bracketing, upgraded to (a)** — report ΔE for CaF₂ at D∞h and at the
  unconstrained minimum as an interval, and replace CaF₂ with a CCSD(T)/CBS focal-point value
  (a 3-atom molecule: hours, not days). Do not report the 128° structure at all.
- The most important consequence for the paper: the ~1–1.5 kcal mol⁻¹ CaF₂ error is
  **correlated across all four method/basis combinations**, so it does *not* appear in the
  3.41 kcal mol⁻¹ spread. The spread is a lower bound on the uncertainty, not an estimate of it.

---

## 1. How do PBE/B3LYP-class functionals do?

**Direction of the error: over-bending, not under-bending.** Nobody in this literature reports
DFT making CaF₂ too linear.

The one directly citable statement is Donald & Hoffmann (JACS 2006, 128, 11236), who compare
B3LYP against MP2 and QCISD for all 20 MX₂:

> "the B3LYP method tends to exaggerate the bending when compared to those two methods but shows
> the best overall agreement with the experimental geometries for the dihalides"

Their B3LYP numbers for CaF₂ (10-valence-electron ECP, 6s6p**5d** with *uncontracted* d on Ca,
cc-pVTZ on F): **θ = 142.8°, E_lin = 0.04 eV (0.92 kcal mol⁻¹), k_lin = −2.09 eV rad⁻².**
Against experiment they tabulate 142(1)°.

So B3LYP, *with an adequate Ca d set*, lands on the right angle and a ~1 kcal mol⁻¹ barrier.
Our four calculations land 14° further bent with a 2× larger barrier — from a basis that is
much smaller on Ca, at both the SVP and TZVP level. The functional is not the variable that
changed; the basis is.

Two honest qualifications:

- **I found no dedicated multi-functional benchmark** (a 10–20 functional study) of group-2
  dihalide bending, and none that isolates PBE or PBE0 for CaF₂. The DFT evidence base is
  essentially "B3LYP works acceptably here" plus the general observation that quasilinear
  geometries are model-chemistry-sensitive.
- **The quantity is intrinsically hard for any semi-local functional.** The barrier is
  2–3 % of the Ca–F bond energy and is a *correlation* effect: HF gives linear or nearly linear
  (Dyke & Wright: 163.5° or 180° depending on basis; Seijo et al. classify CaF₂ as quasilinear
  at the HF/AIMP level). Recovering a 1 kcal mol⁻¹ correlation-driven barrier to better than
  ~50 % is beyond what a GGA/hybrid can be trusted to do, independent of basis.

**Basis sensitivity is the dominant, well-documented effect, and it is large:**

| Evidence | Number |
|---|---|
| Salzner & Schleyer, 16 basis sets, HF and MP2 (as summarised by Sugarman et al. 1999 and Hassett & Marsden 1995) | **θ between 126° and 180°** |
| Wright, Lee & Dyke, MP2, same method, three Ca d sets | E(linear) − E(bent,150°) = **−101, −4, +19 cm⁻¹** (sign flips) |
| Wright, Lee & Dyke, counterpoise on the bend | BSSE alone *reduces the barrier* from ~106 to ~20 cm⁻¹, i.e. uncorrected BSSE **over-bends by ~86 cm⁻¹ (0.25 kcal mol⁻¹)** |
| Dyke & Wright, HF | ν₂ = 46 cm⁻¹ vs 115 cm⁻¹ experimental |

The "126° to 180°" figure is the single most useful number in this assessment: **our 128° sits
at the extreme edge of the range that basis-set incompleteness alone can produce.** It is not
a converged prediction.

The mechanism matters, and it points at the same basis functions. Romeu & Dixon (JACS 2026)
show that the bending is driven by mixing of the alkaline-earth **outer-core (n−1)p** orbitals
into the halide valence space, and Wright et al. showed that changing the Ca *d* space changes
the composition of the "essentially Ca 2s and 3s" MOs — i.e. the tight d functions act as a
proxy for core-polarisation / 3p participation. A basis that under-describes that manifold (no
tight uncontracted d, no f) distorts the bending coordinate specifically.

**Our Ca basis, concretely** (Basis Set Exchange, def2-SVP / def2-TZVP, all-electron, Z=20):

- `def2-SVP`: d primitives 5.498, 1.318, 0.3219 contracted into **one** function, plus one
  diffuse 0.0705 → **[2d]**, no f.
- `def2-TZVP`: (5.498, 1.318)→1, 0.3219→1, 0.0705→1 → **[3d]**, no f.
- Literature requirement: **4–6 d functions, largely uncontracted, largest exponent 7–10.5,
  plus at least one f.**

The Ca–F distance we get (1.9915/1.9927 Å) agrees with the converged B3LYP value (1.996 Å),
so the *stretching* coordinate is fine. Only the bending coordinate — the one that needs the
flexible d manifold — is wrong. That is a clean, self-consistent diagnosis.

## 2. What is benchmark quality, and with what basis?

**Method ladder for this specific problem:**

- **HF** — qualitatively wrong direction (linear / quasilinear). Not usable, but a useful
  reference point: the bend is entirely a correlation effect.
- **MP2** — the workhorse of the 1990s literature (Dyke & Wright: 151.8°/153°; Hassett &
  Marsden: five d and seven f functions, exponents optimised *for the molecule*). MP2 gets the
  angle roughly right but the barrier is ±0.1 kcal mol⁻¹ noise against basis choice. Adequate
  for the angle, inadequate for the barrier.
- **CCSD(T)** — the current standard. Romeu & Dixon's CCSD(T)/aug-cc-pVQZ(-PP) survey of all
  24 MX₂ is the reference dataset. For CaF₂ they report a bent minimum whose **barrier is
  comparable to the bending harmonic frequency** — the operational definition of quasilinear.
- **CASPT2** — used by the same authors as an independent check. Note the calibration of how
  much this matters: **SrCl₂ is linear at CASPT2 but slightly bent at CCSD(T)**. Two
  benchmark-quality methods disagree on the sign of the distortion for a neighbouring molecule.
  That is the size of the irreducible method uncertainty at the quasilinear threshold, and it is
  the same order as CaF₂'s barrier.
- **MRCI(+Q)** — used by Sardar et al. for a full 6-D Ca₂F₂ surface. Their CaF₂:
  **θ = 168.1°**, versus 157.5° and 154.9° from two earlier theory references they tabulate.
  Note the ~13° spread between high-level methods on the same molecule.

**Basis set requirements:**

- **Ca d manifold: essential and large.** Uncontracted or nearly so, 4–6 functions, with tight
  exponents (~7–10) as well as diffuse ones. This is the single requirement the literature is
  unanimous and emphatic about. Wright et al.: a d function of exponent 7.75 added to a set that
  already had four contracted d functions still moved the linear/bent gap by ~100 cm⁻¹.
- **f functions on Ca: yes.** Hassett & Marsden: adding f changed ν₂ by 12–20 cm⁻¹ at correlated
  level; for CO₂ (isoelectronic) f functions change the bend by ~50 cm⁻¹.
- **Diffuse functions: helpful, second-order.** Diffuse p on F changes the bending force constant
  by ~10 % at DZ, only 1.2 % at TZ (Hassett & Marsden). The Ca d/f manifold dominates.
  Modern practice (Romeu & Dixon) uses augmented sets anyway: `aug-cc-pVQZ(-PP)`.
- **Semi-core:** with an ECP on Ca, use the 10-valence-electron one (3s3p in valence), which is
  what both `def2` (all-electron) and Peterson's `cc-pwCVXZ-PP` provide. The mechanism requires
  the 3p to be present and *correlated*, so a large-core ECP is disqualified.

**Numbers for CaF₂ (and the closest analogues):**

| Method / basis | θ (deg) | E_lin | Source |
|---|---|---|---|
| HF, extended spd AIMP | "quasilinear" | small | Seijo, Barandiarán & Huzinaga 1991 |
| HF | 163.5 **or** 180 (basis-dependent) | — | Dyke & Wright 1990 |
| MP2 | 151.8 / 153 | — | Dyke & Wright 1990 |
| MP2, 3 Ca d sets | — | −101 / −4 / +19 cm⁻¹ | Wright, Lee & Dyke 1991 |
| B3LYP, 5 uncontracted d | 142.8 | 0.04 eV (0.92 kcal mol⁻¹) | Donald & Hoffmann 2006 |
| CCSD(T)/awQ(-PP) | bent, quasilinear | ≈ bend frequency | Romeu & Dixon 2026 |
| MRCI+Q/cc-pwCVTZ-PP | 168.1 | — | Sardar et al. 2023 |
| MRCI (literature refs in Sardar) | 157.5 / 154.9 | — | via Sardar et al. 2023 |
| "experiment" (matrix IR, Ar) | 142(1) | — | tabulated in Donald & Hoffmann 2006; Ramondo et al. 1989 |
| matrix IR, Kr | 131 | — | Donald & Hoffmann 2006, footnote |
| closest analogue, SrF₂ | 128.8 | 0.14 eV | Donald & Hoffmann 2006 |

The experimental 142° is an **effective (r₀-type) angle from matrix-isolated IR, and it moves to
131° in a krypton matrix.** It is not a gas-phase r_e and should not be quoted as one; matrix
polarisation of a highly polar molecule is the suspected cause (Sugarman, Wilson & Madden 1999).

## 3. Linear, bent, or quasilinear — and is the question well posed?

**The correct picture is "quasilinear with a flat potential", and the sharp question "linear or
bent?" is not well posed for CaF₂.** Both of those statements are quantitative, not rhetorical.

The three energy scales:

- barrier to linearity: **0.03–0.05 eV ≈ 250–400 cm⁻¹** (0.7–1.2 kcal mol⁻¹)
- bending harmonic frequency: **~110–120 cm⁻¹**
- bending zero-point energy: **~55–60 cm⁻¹**

The barrier supports only about **three bending quanta**. When the barrier is a handful of
vibrational quanta, the bending wavefunction is delocalised over tens of degrees, the bend–rotation
coupling is strong (the A rotational constant of a near-linear triatomic is enormous), and there
is no sharp "equilibrium angle" — only a vibrationally averaged ⟨θ⟩ that depends on the
rotational state and on the temperature. Romeu & Dixon phrase exactly this: the bending harmonic
frequencies of CaF₂ and SrCl₂ are *comparable to the inversion barriers*, so these molecules are
"quasilinear at most temperatures and highly fluxional". Donald & Hoffmann note that the
linear/quasilinear/bent label depends on an **arbitrary cutoff** on the linearisation energy.

So the operationally meaningful questions are:
1. What is the barrier height V_lin (a well-defined number)?
2. What is the vibrationally averaged angle (a well-defined, J-dependent observable)?

and **not** "what is the equilibrium angle", which for this molecule is a property of a
flat-bottomed potential well rather than a sharp structural parameter. A defensible treatment of
the *shape* needs a variational 1-D or 2-D (r, θ) bend calculation, not a single optimisation.

## 4. Feasibility and cost of CCSD(T)

**This is cheap, and the fear of it is misplaced.**

- CaF₂ has 3 atoms. Optimisation in C₂ᵥ has 2 geometry parameters; at D∞h, 1. With a ~200–300
  function basis (aug-cc-pVQZ(-PP) on Ca, aug-cc-pVQZ on F) a CCSD(T) energy is seconds to a
  couple of minutes, and CCSD(T) **analytic gradients** are standard in CFOUR, PSI4 and Molpro.
- A numerical-gradient optimisation costs ~4–6 single points per step and converges in 5–10
  steps from a good start — **a few core-hours, one workstation, well under a day.**
- Romeu & Dixon did exactly this for **24 molecules** including geometry optimisation; Sardar
  et al. built a **full six-dimensional MRCI surface for the 4-atom Ca₂F₂**. A single CaF₂
  CCSD(T) optimisation is a small fraction of either.
- Add a CCSD(T)/aug-cc-pV5Z(-PP) single point plus an MP2 CBS extrapolation and you have a
  focal-point-quality CaF₂ for the price of an afternoon.

**But**: CCSD(T) is not the end of the argument, because at the quasilinear threshold the
method uncertainty (CCSD(T) vs CASPT2, SrCl₂) is comparable to the barrier itself. A CCSD(T)
number buys you the right order of magnitude and a defensible central value; it does not buy
0.1 kcal mol⁻¹. That is worth saying explicitly rather than hiding.

**Verdict: yes, it is a sensible target** — but as a *focal-point correction for one species in
an error budget*, not as the centrepiece of the paper. The paper's contribution is the
certificate machinery; CaF₂ should be handled in two paragraphs, not a methods section.

## 5. Recommendation

**Take (c), implemented with (a). Reject (b) alone.**

Reasoning:

- **(b) is not defensible on its own.** Constraining CaF₂ to D∞h and "acknowledging" it produces
  a reaction energy to a *saddle point* of the product, i.e. a number that is not the reaction
  energy, with a one-sided error of known sign (too endothermic by ~0.9–1.2 kcal mol⁻¹). An
  acknowledged systematic error of known sign is still an error; fixing it costs nothing.
- **(a) alone is over-engineering for this paper.** A CCSD(T)/CBS CaF₂ is affordable, but it
  does not remove the residual method ambiguity (CCSD(T) vs CASPT2 at the quasilinear limit),
  so it cannot be presented as "the answer".
- **(c) is the honest frame** — but the *width* of the bracket should come from the literature
  barrier (0.03–0.05 eV), not from our own 2 kcal mol⁻¹ value, which is itself the artefact.

**Concretely, in order of cost:**

1. **Drop the 128.2° structure entirely.** Gradient norm 3.0×10⁻² Ha Å⁻¹, mutually inconsistent
   with its own Hessian (see §6), and outside the 126–180° envelope that basis incompleteness
   alone spans. Do not report it, even as a caveat.
2. **Report ΔE twice**: once with CaF₂ constrained to D∞h (well defined, reproducible), once with
   CaF₂ at its own minimum (the ~140° angle-scan value is more defensible than the 128° point).
   Give both, and give the difference (~1.9–2.1 kcal mol⁻¹) as the bracketing width.
3. **Add one CCSD(T)/aug-cc-pVQZ(-PP) optimisation of CaF₂** and splice it in as a documented
   focal-point correction Δ = E_CCSD(T)(CaF₂) − E_DFT(CaF₂). Report Δ explicitly, with its sign.
4. **State the direction and size of the residual**: DFT/basis of this class over-bends CaF₂ by
   roughly a factor 2 in the barrier, i.e. ΔE is ~1.0–1.5 kcal mol⁻¹ too exothermic, common to
   all four combinations.

**The single most important sentence to add to the paper:** because all four method/basis pairs
share the same Ca basis deficiency, the CaF₂ error is *correlated*, so the 3.41 kcal mol⁻¹
spread does not contain it. **The spread should be presented as a lower bound on the total
uncertainty**, with the CaF₂ bending contribution added separately (~1.0–1.5 kcal mol⁻¹,
one-sided). This is a stronger and more honest claim than "the method spread is 3.4 and a flat
bending potential adds comparable uncertainty" — and it is a claim the certificate framing of
the paper actually needs, since a bound that omits a correlated systematic is not a bound.

Also correct: "Gas-phase CaF₂ is conventionally linear" → "CaF₂ is quasilinear: bent with a
barrier to linearity of a few tenths of a kcal mol⁻¹, a saddle at D∞h, and a barrier comparable
to its bending frequency."

## 6. Internal consistency check on our own numbers (arithmetic only, no QC run)

Using the bend G-matrix element for F–Ca–F, g_θθ = 2/(m_F r²) with r = 2.05 Å:
**μ_bend = 39.9 amu Å² = 6.63×10⁻⁴⁶ kg m².**

- ω_imag = 150 cm⁻¹ ⟺ |V″| at 180° = **0.53 aJ rad⁻² = 3.3 eV rad⁻²**.
- Cross-check against the converged B3LYP value k_lin = 2.09 eV rad⁻² ⟹ ω_imag ≈ **119 cm⁻¹**.
  Our −103.5 to −154.6 cm⁻¹ is *consistent* with that. **The curvature at linear is fine.**
- Now the depth. For a quartic well V = a(θ−θ₀)²(θ−180°)² with that curvature:
  - depth 2.06 kcal mol⁻¹ ⟹ θ_min ≈ **161°**, not 140°, and certainly not 128°;
  - θ_min = 140° ⟹ depth ≈ **9 kcal mol⁻¹**;
  - θ_min = 128° ⟹ depth ≈ **16 kcal mol⁻¹**;
  - θ_min = 140° with depth 2.06 kcal mol⁻¹ ⟹ ω_imag ≈ **71 cm⁻¹** (observed: 103–155);
  - θ_min = 128° with depth 2.06 kcal mol⁻¹ ⟹ ω_imag ≈ **54 cm⁻¹**.
- A cubic+quartic potential fixed to (θ_min = 128°, depth 2 kcal mol⁻¹, V″(180°) = −0.53 aJ rad⁻²)
  has **V″(θ_min) < 0**: the 128° point is not a minimum of any such potential.

**Interpretation:** the angle, the depth and the imaginary frequency in hand cannot all come from
one smooth bending potential. The imaginary frequency is right (it matches B3LYP's k_lin), and the
depth/angle are not. Note also that the 140° angle-scan value and the 128.2° re-optimised value
disagree with each other by 12° — consistent with the reported gradient norm of 3×10⁻² Ha Å⁻¹,
i.e. the geometry is sliding on a nearly flat, poorly described surface rather than converging
to a minimum.

This is a *model-consistency* argument only (a 1-D quartic fit to three cited numbers); it does
not identify which of the three numbers is wrong. It does establish that quoting "128°, 2 kcal
mol⁻¹ below linear" as a physical well depth is not supported even internally.

---

## 7. What I could not find (explicit)

1. **No numerical CCSD(T)/awQ(-PP) angle or barrier for CaF₂ from Romeu & Dixon 2026.** The
   article PDF is paywalled and returned only the abstract; the Supporting Information URL
   returned HTTP 403; the local browser could not start (missing `libnspr4.so`), so I could not
   try the HTML full text. I have only their abstract-level statements (quasilinear; barrier
   comparable to the bending frequency). I therefore quote **no** CCSD(T) number for CaF₂.
2. **No dedicated multi-functional DFT benchmark** for group-2 dihalide bending angles, and no
   study isolating PBE or PBE0 on CaF₂. The "DFT over-bends" conclusion rests on B3LYP in
   Donald & Hoffmann plus the general model-chemistry sensitivity both they and Kaupp report.
3. **Primary records not retrieved** (seen only as summaries inside retrieved papers, and
   therefore cited only that way): von Szentpály & Schwerdtfeger, Chem. Phys. Lett. 1990, 170,
   555 (SCF 159°, MP2 154.8°); the Salzner & Schleyer 16-basis-set study (via Sugarman et al.
   1999 and Hassett & Marsden 1995). I did not locate their DOIs and do not cite them directly.
4. **No gas-phase rotationally resolved r_e angle** for CaF₂. The 142° figure is an effective
   angle from matrix-isolated IR (Ar), shifting to 131° in Kr. The matrix-free gas-phase
   geometry appears never to have been determined to that precision — the electric-deflection
   experiment establishes only that a dipole exists.
5. **No published CCSD(T) or MRCI barrier height in cm⁻¹ or kcal mol⁻¹ for CaF₂** that I could
   retrieve. All quantitative barrier values I can cite trace to B3LYP (0.04 eV) or to MP2
   (±0.1 kcal mol⁻¹, Wright et al.).
6. **No anharmonic / vibrationally averaged angle** for CaF₂ at any level. The quasilinear
   discussion in the literature is qualitative (barrier vs frequency), not a computed ⟨θ⟩.
7. **No counterpoise-corrected DFT results** for CaF₂. The BSSE study (Wright, Lee & Dyke) is
   MP2 only, and its conclusion — BSSE is *not* the root cause, the Ca sp space is — has not
   to my knowledge been repeated at DFT level.

---

## 8. Sources (all retrieved during this task; metadata verified against Crossref)

1. J. G. F. Romeu, D. A. Dixon, "When Core Orbitals Act as Valence Orbitals: Linear and Bent
   Geometries of MX₂ for X a Halogen and M an Alkaline Earth", *J. Am. Chem. Soc.* **2026**,
   148(1), 828–836. DOI: 10.1021/jacs.5c16443.
2. K. J. Donald, R. Hoffmann, "Solid Memory: Structural Preferences in Group 2 Dihalide
   Monomers, Dimers, and Solids", *J. Am. Chem. Soc.* **2006**, 128, 11236–11249.
   DOI: 10.1021/ja062817j.
3. L. Seijo, Z. Barandiarán, S. Huzinaga, "Ab initio model potential study of the equilibrium
   geometry of alkaline earth dihalides: MX₂ (M=Mg, Ca, Sr, Ba; X=F, Cl, Br, I)",
   *J. Chem. Phys.* **1991**, 94, 3762–3773. DOI: 10.1063/1.459748.
4. M. Kaupp, P. v. R. Schleyer, H. Stoll, H. Preuss, "The question of bending of the alkaline
   earth dihalides MX₂ … An ab initio pseudopotential study", *J. Am. Chem. Soc.* **1991**,
   113, 6012–6020. DOI: 10.1021/ja00016a013.
5. M. Kaupp, P. v. R. Schleyer, H. Stoll, H. Preuss, "Pseudopotential approaches to Ca, Sr, and
   Ba hydrides. Why are some alkaline earth MX₂ compounds bent?", *J. Chem. Phys.* **1991**,
   94, 1360–1366. DOI: 10.1063/1.459993.
6. J. M. Dyke, T. G. Wright, "An ab initio study of the equilibrium geometries and vibrational
   frequencies of the group IIa difluorides MF₂ (M=Be, Mg, Ca)", *Chem. Phys. Lett.* **1990**,
   169, 138–144. DOI: 10.1016/0009-2614(90)85178-F.
7. T. G. Wright, E. P. F. Lee, J. M. Dyke, "Effects of BSSE and d-orbital space in the
   calculation of the equilibrium geometry of CaF₂", *Mol. Phys.* **1991**, 73, 941–951.
   DOI: 10.1080/00268979100101671.
8. D. M. Hassett, C. J. Marsden, "The shapes and bending potentials of monomeric MgF₂ and CaF₂",
   *J. Mol. Struct.* **1995**, 346, 249–263. DOI: 10.1016/0022-2860(94)09020-P.
9. R. Sugarman, M. Wilson, P. A. Madden, "Matrix effects on the geometry of CaF₂ — a consequence
   of polarization effects?", *Chem. Phys. Lett.* **1999**, 308, 509–515.
   DOI: 10.1016/S0009-2614(99)00610-7.
10. F. Ramondo, L. Bencivenni, S. Nunziante Cesaro, K. Hilpert, "FTIR matrix isolation study on
    gaseous calcium dihalide molecules", *J. Mol. Struct.* **1989**, 192, 83–94.
    DOI: 10.1016/0022-2860(89)87008-5.
11. M. Hargittai, "Molecular Structure of Metal Halides", *Chem. Rev.* **2000**, 100, 2233–2302.
    DOI: 10.1021/cr970115u.
12. L. von Szentpály, "Hard Bends Soft: Bond Angle and Bending Force Constant Predictions for
    Dihalides, Dihydrides, and Dilithides of Groups 2 and 12", *J. Phys. Chem. A* **2002**,
    106, 11945–11949. DOI: 10.1021/jp026658b.
13. D. Sardar, A. Christianen, H. Li, J. L. Bohn, "Four-body singlet potential-energy surface for
    reactions of calcium monofluoride", *Phys. Rev. A* **2023**, 107, 032822.
    DOI: 10.1103/PhysRevA.107.032822.
14. L. Wharton, R. A. Berg, W. Klemperer, "Geometry of the Alkaline-Earth Dihalides",
    *J. Chem. Phys.* **1963**, 39, 2023–2031. DOI: 10.1063/1.1734577.
15. J. G. Hill, K. A. Peterson, "Gaussian basis sets for use in correlated molecular
    calculations. XI. Pseudopotential-based and all-electron relativistic basis sets for alkali
    metal (K–Fr) and alkaline earth (Ca–Ra) elements", *J. Chem. Phys.* **2017**, 147, 244106.
    DOI: 10.1063/1.5010587.
16. R. L. De Kock, M. A. Peterson, L. K. Timmer, E. J. Baerends, P. Vernooijs, "A theoretical
    study of the linear versus bent geometry for several MX₂ molecules: MgF₂, CaH₂, CaF₂, CeO₂
    and YbCl₂", *Polyhedron* **1990**, 9, 1919–1934. DOI: 10.1016/S0277-5387(00)84004-8.
