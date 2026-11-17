# Gas-phase CaF2 and the alkaline-earth dihalides: linear, bent, or quasilinear?
Literature assessment compiled 2026-11-17. Only records whose metadata was retrieved
(Crossref DOI records and/or publisher abstract pages) are cited. Items not found are listed at the end.

## 1. Verdict for CaF2

CaF2 is **quasilinear**: a bent equilibrium geometry sitting behind a very small barrier at
the linear configuration, with a large-amplitude bending vibration. It is not "linear" in the
rigid-rotor sense, but at ordinary temperatures it behaves as an effectively linear molecule.

High-level numbers (all electronic-structure, no gas-phase experiment):

| Quantity | Value | Source |
|---|---|---|
| theta_e (F-Ca-F) | 154.3 deg | Koput & Roszczak 2004, CCSD(T)/aug-cc-pVQZ+5Z |
| barrier to linearity | 54 cm-1 (0.15 kcal/mol) | Koput & Roszczak 2004 |
| nu_2 fundamental (variational) | 9 cm-1 | Koput & Roszczak 2004 |
| theta_e | 157 deg | Vasiliu et al. 2010, CCSD(T)/awCVTZ |
| theta_e | 157.5 deg | Kaupp et al. 1991, CISD+Q, 6s6p5d1f + ECP |
| classification "quasilinear at most temperatures, highly fluxional" | — | Romeu & Dixon 2026, CCSD(T)/awQ(-PP), CASPT2/awQ(-PP) |

Experimental constraints on the angle are matrix-isolation numbers only:
Kr matrix, 20 K: mean apex angle 140 deg (Calder et al. 1969); quoted as 135 +/- 7 deg
(Vasiliu et al. 2018, citing Calder); 139-156 deg across Ar/Kr/N2/CO matrices
(Beattie, as summarised by Sugarman et al. 1999). Ar matrix 142 deg, Kr matrix 131 deg
(Ramondo, as tabulated by Donald & Hoffmann 2006, footnote to Table 1).
The frequently quoted "about 145 deg" is a matrix value, not a gas-phase value
(Hassett & Marsden 1995 call it an "apparent" angle).

**No gas-phase determination of the CaF2 angle exists.** There is no microwave/rotational
spectrum, no gas-phase high-resolution IR, and the only electron-diffraction study cited in
this literature (Akishin & Spiridonov 1957) is described by Sugarman et al. (1999) as
indicating "a linear ground-state configuration".

## 2. Experimental evidence, by technique

- **Electric deflection (molecular beam):** CaF2 deflects -> permanent dipole -> not
  centro-symmetric. Wharton, Berg & Klemperer 1963 list CaF2, SrF2, SrCl2, BaF2, BaCl2,
  BaBr2, BaI2 as polar, and CaCl2, CaBr2, SrBr2, SrI2 as non-polar. Buchler, Stauffer &
  Klemperer 1964 confirm CaF2 "deflects in the manner expected for a polar molecule".
  This proves bending on the timescale of the experiment but yields no angle.
- **Matrix-isolation IR (Ar, Kr, N2, CO; 10-30 K):** isotopic shifts of nu_2/nu_3 give the
  apex angle; the angle moves by >20 deg between matrices, so the method measures
  molecule-plus-cage, not the free molecule (Calder et al. 1969; Snelson 1966;
  Ramondo, Rossi & Bencivenni 1988; Ramondo et al. 1989; Sugarman et al. 1999).
- **Gas-phase electron diffraction:** low-precision, linear-consistent (Akishin & Spiridonov
  1957, as cited by Vasiliu et al. 2018 and Sugarman et al. 1999). For the analogue SrCl2 a
  modern high-temperature GED study gives 142.4 +/- 4.0 deg (Varga et al. 2006).
- **Microwave / high-resolution gas-phase IR:** none found for CaF2.

## 3. Family table (M = Be, Mg, Ca, Sr, Ba, Ra; X = F, Cl, Br, I)

Classification follows Romeu & Dixon 2026 (CCSD(T)/awQ(-PP) and CASPT2/awQ(-PP)); experimental
notes from the sources named. "QL" = quasilinear/floppy (flat potential, barrier comparable to
or below the bending frequency).

| | F | Cl | Br | I |
|---|---|---|---|---|
| **Be** | linear | linear | linear | linear |
| **Mg** | linear (matrix IR "158 deg" is a matrix artefact; GED + deflection linear) | linear | linear | linear |
| **Ca** | **bent / QL**, theta_e 154-157 | linear | linear | linear |
| **Sr** | **bent** (matrix 108 deg; calc ~130-139) | **QL**: CASPT2 linear, CCSD(T) slightly bent; GED 142.4 +/- 4.0 deg | linear at CCSD(T) but very floppy (barrier ~0.001 eV) -> QL | linear, listed QL by Kaupp's <4 kJ criterion |
| **Ba** | **bent** (calc 118-123; matrix 100 deg) | **bent** (calc 127; exp 120 +/- 10 deg) | **bent** (calc 130; exp 137.1 +/- 4.9 deg) | **bent** (calc 134; exp 137.6 +/- 0.9 deg) |
| **Ra** | **bent** | **bent** | **bent** | **bent** |

Trend: bending increases down the group and from I to F. Only CaF2 is bent in the CaX2 row;
in the Sr row bending starts at SrF2/SrCl2; all Ba and Ra dihalides are bent.

Definitions and classification conventions:
- Working definition of quasilinearity (barrier at linear <= bending frequency): Varga et al. 2006.
- Quantitative quasilinearity parameter gamma_0: Yamada & Winnewisser 1976.
- Kaupp's operational cut-off (quasilinear if <4 kJ/mol to deviate >20 deg from linear):
  summarised in Donald & Hoffmann 2006, which puts all Ca dihalides, SrCl2, SrBr2, SrI2 and
  BaI2 in that class.
- Three-tier scheme (linear / quasilinear / bent) reviewed in Hargittai & Hargittai 1993 and
  Hargittai 2000.

## 4. Bending frequency of CaF2

- Matrix-isolated: bending fundamental reported near 162 cm-1 (value quoted by Sugarman et al.
  1999 as the observed matrix frequency, vs <100 cm-1 from high-quality calculations).
- Evaluated gas-phase compilation: NIST CCCBDB (SRD 101) lists fundamentals 490 (a1),
  120 (a1), 575 (b2) cm-1 for CaF2, attributed to Gurvich. Vasiliu et al. 2018 note that such
  compiled values derive from Kr-matrix IR + normal-coordinate analysis, not from gas-phase
  measurement.
- Best variational prediction on a CCSD(T) surface: 9 cm-1 (Koput & Roszczak 2004).
- **No direct gas-phase experimental nu_2 found.** The spread 9 / 120 / 162 cm-1 is itself the
  signature of a quasilinear double-minimum potential plus strong matrix perturbation.

## 5. Methods reported to reproduce the CaF2 bending potential

- CCSD(T) with quadruple/quintuple-zeta basis sets and explicit treatment of Ca semicore:
  theta_e 154.3 deg, barrier 54 cm-1 (Koput & Roszczak 2004); 157 deg at CCSD(T)/awCVTZ
  (Vasiliu et al. 2010, 2018).
- CASPT2/awQ(-PP) agrees with CCSD(T) on the quasilinear classification (Romeu & Dixon 2026).
- CISD+Q with 6s6p5d1f valence basis + Ca pseudopotential: 157.5 deg (Kaupp et al. 1991).
- SA-CASSCF natural orbitals explain the physics: mixing of the Ca outer-core (n-1)p orbitals
  with the halide valence space drives bending; core d and valence d contribute little
  (Romeu & Dixon 2026). This is why the Ca 3p/core-valence treatment is decisive.
- Basis-set requirements: large, uncontracted d sets on Ca are needed even at HF to obtain a
  bent structure (Hassett & Marsden 1990); BSSE can artificially favour bending, and 16
  different basis sets at HF/MP2 spanned 126-180 deg (Salzner & von Rague Schleyer 1990;
  BSSE study, Mol. Phys. 1991).
- DFT warning: B3LYP "tends to exaggerate the bending" relative to MP2/QCISD
  (Donald & Hoffmann 2006). PBE/B3LYP with def2 bases giving a ~2 kcal/mol well here is
  therefore ~13x deeper than the 0.15 kcal/mol CCSD(T) barrier - qualitatively right
  (linear is a saddle point) but quantitatively wrong.

## 6. Not found / not retrievable

- No gas-phase microwave, rotational, or high-resolution IR spectrum of CaF2.
- No modern gas-phase electron-diffraction structure of CaF2; the 1957 Akishin & Spiridonov
  paper (Kristallografiya 2, 475) could not be retrieved, only its title and its use as a
  citation by others.
- Full texts of Koput & Roszczak 2004, Romeu & Dixon 2026 and Hargittai 2000 were paywalled;
  their numbers here come from retrieved abstracts and from tabulations in Donald & Hoffmann
  2006 and Vasiliu et al. 2018.
- No numeric CCSD(T)/CASPT2 angles from Romeu & Dixon 2026 (abstract gives classification
  only).
- No experimental gas-phase structure for SrBr2 or SrI2 beyond "quasilinear" assignments.

## References (metadata verified)

1. Wharton, L.; Berg, R. A.; Klemperer, W. "Geometry of the Alkaline-Earth Dihalides."
   J. Chem. Phys. 1963, 39 (8), 2023-2031. DOI 10.1063/1.1734577
2. Buchler, A.; Stauffer, J. L.; Klemperer, W. "Geometry of the Transition-Metal Dihalides:
   The Fluorides of Manganese, Cobalt, Nickel, Copper, and Zinc." J. Chem. Phys. 1964,
   40, 3471-3474. DOI 10.1063/1.1725038
3. Snelson, A. "Infrared Spectra of Some Alkaline Earth Halides by the Matrix Isolation
   Technique." J. Phys. Chem. 1966, 70 (10), 3208-3217. DOI 10.1021/j100882a031
4. Mann, D. E.; Calder, V.; Seshadri, K. S.; White, D. "Geometry and Vibrational Spectra of
   the Alkaline-Earth Dihalides. I. MgF2." J. Chem. Phys. 1967, 46, 1138-1143.
   DOI 10.1063/1.1840780
5. Calder, V.; Mann, D. E.; Seshadri, K. S.; Allavena, M.; White, D. "Geometry and Vibrational
   Spectra of Alkaline-Earth Dihalides. II. CaF2, SrF2, and BaF2." J. Chem. Phys. 1969,
   51 (5), 2093-2099. DOI 10.1063/1.1672304
6. White, D.; Calder, G. V.; Hemple, S.; Mann, D. E. "Geometry and vibrational spectra of the
   alkaline earth dihalides. III. MgCl2, CaCl2, SrCl2, and BaCl2." J. Chem. Phys. 1973,
   59 (12), 6645-6651. DOI 10.1063/1.1680046
7. Yamada, K.; Winnewisser, M. "A Parameter to Quantify Molecular Quasilinearity."
   Z. Naturforsch. A 1976, 31, 139-144. DOI 10.1515/zna-1976-0206
8. Ramondo, F.; Rossi, V.; Bencivenni, L. "Structure and vibrational frequencies of CaF2
   aggregates isolated in a low temperature matrix." Mol. Phys. 1988, 64, 513-520.
   DOI 10.1080/00268978800100363
9. Ramondo, F.; Bencivenni, L.; Nunziante Cesaro, S.; Hilpert, K. "FTIR matrix isolation study
   on gaseous calcium dihalide molecules." J. Mol. Struct. 1989, 192, 83-94.
   DOI 10.1016/0022-2860(89)87008-5
10. Hassett, D. M.; Marsden, C. J. "The influence of d orbitals on the shape of monomeric
    CaF2." J. Chem. Soc., Chem. Commun. 1990, 667-669. DOI 10.1039/C39900000667
11. Salzner, U.; von Rague Schleyer, P. "Ab initio studies of the geometries and electronic
    structures of CaF2 and CaCl2." Chem. Phys. Lett. 1990, 172 (6), 461-470.
    DOI 10.1016/0009-2614(90)80140-9
12. von Szentpaly, L.; Schwerdtfeger, P. "Which double-octet ABC molecules are bent? CI
    calculations on CaF2, and a softness criterion to predict bending." Chem. Phys. Lett.
    1990, 170 (5-6), 555-560. DOI 10.1016/S0009-2614(90)87101-V
13. Kaupp, M.; Schleyer, P. v. R.; Stoll, H.; Preuss, H. "The question of bending of the
    alkaline earth dihalides MX2 (M = Be, Mg, Ca, Sr, Ba; X = F, Cl, Br, I). An ab initio
    pseudopotential study." J. Am. Chem. Soc. 1991, 113 (16), 6012-6020.
    DOI 10.1021/ja00016a013
14. Seijo, L.; Barandiaran, Z.; Huzinaga, S. "Ab initio model potential study of the
    equilibrium geometry of alkaline earth dihalides: MX2 (M=Mg, Ca, Sr, Ba; X=F, Cl, Br, I)."
    J. Chem. Phys. 1991, 94 (5), 3762-3773. DOI 10.1063/1.459748
15. Kaupp, M.; Schleyer, P. v. R.; Stoll, H. "Model calcium hydride CaH2(L) and calcium
    fluoride CaF2(L) complexes (L = Ne, Ar, Kr, Xe, CO, N2)." J. Phys. Chem. 1992, 96,
    9801-9805. DOI 10.1021/j100203a041
16. Hargittai, M.; Hargittai, I. "Linear, Bent, and Quasilinear Molecules." In Structures and
    Conformations of Non-Rigid Molecules; 1993; pp 465-489.
    DOI 10.1007/978-94-011-2074-6_23
17. Hassett, D. M.; Marsden, C. J. "The shapes and bending potentials of monomeric MgF2 and
    CaF2." J. Mol. Struct. 1995, 346, 249-263. DOI 10.1016/0022-2860(94)09020-P
18. Sugarman, R.; Wilson, M.; Madden, P. A. "Matrix effects on the geometry of CaF2 - a
    consequence of polarization effects?" Chem. Phys. Lett. 1999, 308 (5-6), 509-515.
    DOI 10.1016/S0009-2614(99)00610-7
19. Hargittai, M. "Molecular Structure of Metal Halides." Chem. Rev. 2000, 100 (6),
    2233-2302. DOI 10.1021/cr970115u
20. von Szentpaly, L. "Hard Bends Soft: Bond Angle and Bending Force Constant Predictions for
    Dihalides, Dihydrides, and Dilithides of Groups 2 and 12." J. Phys. Chem. A 2002,
    106 (49), 11945-11949. DOI 10.1021/jp026658b
21. Koput, J.; Roszczak, A. "CaF2 As a Quasilinear Molecule: the Vibrational-Rotational Energy
    Levels Predicted by Ab Initio Quantum Chemistry Approach." J. Phys. Chem. A 2004,
    108 (42), 9267-9273. DOI 10.1021/jp046682k
22. Varga, Z.; Lanza, G.; Minichino, C.; Hargittai, M. "Quasilinear Molecule par Excellence,
    SrCl2: Structure from High-Temperature Gas-Phase Electron Diffraction and
    Quantum-Chemical Calculations." Chem. Eur. J. 2006, 12 (32), 8345-8357.
    DOI 10.1002/chem.200600328
23. Donald, K. J.; Hoffmann, R. "Solid Memory: Structural Preferences in Group 2 Dihalide
    Monomers, Dimers, and Solids." J. Am. Chem. Soc. 2006, 128 (34), 11236-11249.
    DOI 10.1021/ja062817j
24. Garcia-Fernandez, P.; Bersuker, I. B.; Boggs, J. E. "Why Are Some ML2 Molecules (M = Ca,
    Sr, Ba; L = H, F, Cl, Br) Bent while Others are Linear? Implications of the Pseudo
    Jahn-Teller Effect." J. Phys. Chem. A 2007, 111, 10409-10415. DOI 10.1021/jp073207o
25. Vasiliu, M.; Li, S.; Feller, D.; Gole, J. L.; Dixon, D. A. "Structures and Heats of
    Formation of Simple Alkaline Earth Metal Compounds: Fluorides, Chlorides, Oxides, and
    Hydroxides for Be, Mg, and Ca." J. Phys. Chem. A 2010, 114, 9349-9358.
    DOI 10.1021/jp1050657
26. Vasiliu, M.; Hill, J. G.; Peterson, K. A.; Dixon, D. A. "Structures and heats of formation
    of simple alkaline earth metal compounds II: fluorides, chlorides, oxides, and hydroxides
    for Ba, Sr, and Ra." J. Phys. Chem. A 2018, 122 (1), 316-327. DOI 10.1021/acs.jpca.7b09056
27. Romeu, J. G. F.; Dixon, D. A. "When Core Orbitals Act as Valence Orbitals: Linear and Bent
    Geometries of MX2 for X a Halogen and M an Alkaline Earth." J. Am. Chem. Soc. 2026,
    148 (1), 828-836. DOI 10.1021/jacs.5c16443
