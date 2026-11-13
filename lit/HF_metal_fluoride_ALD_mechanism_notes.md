# HF-based metal-fluoride ALD: surface mechanism, HF/CaF2(111), amidine byproduct, activation energies

Compiled 2026-11-13. Every DOI below was resolved through `https://api.crossref.org/works/<DOI>`
on that date; bibliographic fields are transcribed from the Crossref record.
Labels: **[S]** = stated in the retrieved source; **[I]** = my inference, not in any source.

---

## 1. Accepted mechanism when HF (or HF-pyridine) is the fluorine co-reactant

### 1.1 What the George-group data actually show

**AlF3 ALD from Al(CH3)3 (TMA) + HF-pyridine** — Lee, DuMont, Cavanagh, George,
*J. Phys. Chem. C* 2015, 119, 14185–14194, DOI 10.1021/acs.jpcc.5b02625.

- **[S]** In situ QCM at 150 °C: MGPC = 31 ng cm^-2 cycle^-1, split into ΔM(TMA) = 22 and
  ΔM(HF) = +9 ng cm^-2 cycle^-1; ΔM(TMA)/MGPC = 0.71 vs 0.67 predicted for x = 1
  (one adsorbed HF per AlF3 deposited). Derived x = 0.78 at 150 °C, 0.99 at 100 °C,
  where the maximum growth rate 1.43 Å/cycle occurs.
- **[S]** Both half-reactions saturate (self-limiting) within one minidose.
- **[S]** In situ FTIR: after HF exposure, positive features at 880 and 920 cm^-1 assigned to
  **adsorbed HF molecules**, 900 cm^-1 to an out-of-plane librational ("reciprocating") HF mode,
  and 3000–3675 cm^-1 to isolated and hydrogen-bonded Al–(HF)* stretching vibrations.
  After TMA exposure those HF features go negative while AlF(CH3)2 features appear
  (725 cm^-1 CH3 rock, 600–650 cm^-1 Al–F, 2800–3000 cm^-1 C–H).
- **[S]** QMS: CH4 (m/z 16) appears during **both** the TMA and the HF exposures;
  no pyridine signal (m/z 52, 79) is seen, i.e. the HF-pyridine azeotrope delivers
  essentially pure HF vapour (measured static HF vapour pressure 90–100 Torr at RT).
- **[S]** Proposed scheme: (A) TMA reacts with HF adsorbed on the AlF3 surface to give
  adsorbed AlF(CH3)2 + CH4; (B) HF converts adsorbed AlF(CH3)2 to AlF3, again releasing CH4,
  with additional HF left adsorbed. The paper also invokes a Lewis acid–base picture
  (AlF3 = F^- acceptor, F in HF = Lewis base) and cites DFT for strong HF coordination
  on α- and β-AlF3.

**ZrF4, HfF4, MnF2, MgF2, ZnF2, LiF ALD with HF-pyridine** — Lee, Sun, Young, George,
*Chem. Mater.* 2016, 28, 2022–2032, DOI 10.1021/acs.chemmater.5b04360.

- **[S]** Mechanism written for ZrF4/TEMAZ: Zr(NEtMe)4 reacts with HF adsorbed on the ZrF4
  surface to give adsorbed ZrF2(NEtMe)2 + **HNEtMe** ("a volatile reaction product"); a second
  HF exposure converts ZrF2(NEtMe)2 to ZrF4, again releasing HNEtMe, with HF left adsorbed.
  MnF2/Mn(EtCp)2 is written the same way with **HCpEt** as the volatile product.
- **[S]** QCM mass balance at 150 °C (ng cm^-2 cycle^-1; MP = metal precursor):

  | film | precursor | MGPC | ΔM(MP) | ΔM(HF) | HF per MFy (x) |
  |---|---|---|---|---|---|
  | ZrF4 | TEMAZ | 35.5 | 39.2 | **−3.7** | 2.4 |
  | HfF4 | TDMAH | 57.2 | 58.4 | **−1.2** | 2.2 |
  | ZrF4 | Zr(OtBu)4 | 29.8 | 42.0 | **−12.2** | 2.0 |
  | AlF3 | TMA | 31.4 | 22.3 | **+9.1** | 0.8 |
  | MnF2 | Mn(EtCp)2 | 15.6 | 40.8 | **−25.2** | 0 |
  | MgF2 | Mg(EtCp)2 | 12.3 | 39.5 | **−27.2** | 0.1 |
  | ZnF2 | DEZ | 33.5 | 41.4 | **−7.9** | 0 |
  | LiF | LiHMDS | 12.2 | 67.6 | **−55.4** | 0.1 |

- **[S]** The sign of ΔM(HF) is the diagnostic: for strong Lewis-acid fluorides (ZrF4, HfF4)
  the HF exposure is a net **mass gain** because HF adsorbs; for the weak Lewis-acid
  difluorides the HF exposure is a large net **mass loss** even though little/no HF is
  retained, which the authors attribute to loss of the organic ligand as a volatile
  hydride/hydrocarbon. x tracks Lewis acidity: 2.0–2.4 for ZrF4/HfF4 vs ≈0 for
  MnF2/MgF2/ZnF2, 0.8 for AlF3, 0.1 for LiF.

### 1.2 Verdict on the three candidate mechanisms

- **[S]** The HF half-reaction is written as **protonation of the metal–ligand bond**
  (HF + M–CH3 → M–F + CH4; HF + M–NR2 → M–F + HNR2), i.e. a Brønsted acid–base /
  σ-bond-metathesis step in which HF is the proton donor and F^- ends up on the metal.
- **[S]** HF additionally acts as an **adsorbed fluorine reservoir**: HF* is present on the
  surface after the HF exposure, quantified by x, and its coverage (not gas-phase HF) sets
  the stoichiometry and correlates with Lewis acidity and with growth rate.
- **[S]** "Ligand exchange" in this literature is reserved for the *other* half-reaction —
  a metal precursor accepting F from a fluorinated surface while transferring its own
  ligand (e.g. AlMe3 + AlF3 → AlF(3−x)Mex + Al2F2Me4; Hu & Schuster, *J. Phys. Chem. C* 2022,
  126, 7410–7420, DOI 10.1021/acs.jpcc.2c00158). It is not the HF step.
- **[I]** So the honest statement is: *protonation-driven ligand removal plus a molecular-HF
  surface reservoir*, and "ligand exchange" is a mislabel if applied to the HF pulse.

### 1.3 What QCM can and cannot establish

- **[S]** QCM measures only total mass; x is **derived** from the measured ΔM ratios under an
  assumed stoichiometry (eqs 3–5 of the AlF3 paper, eqs 3–4 of the Chem. Mater. paper).
- **[I]** QCM therefore does not identify the byproduct. The identity of the leaving species
  (CH4, HNEtMe, HCpEt) comes from QMS/FTIR (AlF3 system) or from assumed chemistry plus mass
  balance (ZrF4/HfF4/MnF2/MgF2/ZnF2 systems — no QMS was reported there).
- **[S]** No in situ QCM study of a CaF2 process was found. The CaF2 ALD paper's abstract
  (Hennessy et al. 2023, below) reports only linear growth, 0.4→0.3 Å/cycle over 175–250 °C,
  optical properties and SF6-plasma etch selectivity — no mass-balance or mechanism data.

---

## 2. Does HF adsorb molecularly or dissociatively on CaF2(111)?

Only quantitative source found: Kaawar & Paulus, *J. Fluorine Chem.* 2019, 224, 67–72,
DOI 10.1016/j.jfluchem.2018.09.006 (periodic DFT, CRYSTAL09, PBE/PBE0/HF/PBE-D3, CaF2, SrF2,
BaF2 (111), (110), (100) at 100/50/25 % coverage).

- **[S]** Starting geometry: HF placed atop a surface cation, F pointing down, 2.5 Å.
- **[S]** Adsorption is **molecular on CaF2(111) at all three coverages** by the tabulated
  intramolecular H–F distances: 0.98 Å (100 %), 1.00 Å (50 %), 1.03 Å (25 %) — only slightly
  stretched from free HF — with the H hydrogen-bonded to surface F (H···F 1.56 / 1.47 / 1.33 Å)
  and the adsorbed F 2.34–3.02 Å from the nearest surface F.
- **[S]** E_ads(CaF2(111)) = −0.40 (100 %), −0.48 (50 %), −0.62 eV (25 %); overall range across
  all surfaces/materials −0.4 to −1.7 eV, "chemisorption in most of the cases", interaction
  "mainly governed by electrostatics" plus H-bonding to surface fluorines.
- **[S]** Clear **dissociation to an adsorbed bifluoride FHF^-** is reported for the **(100)**
  surface (H–F stretched to 1.16 Å, F···F ≈ 2.3 Å, confirmed by Mulliken charges); the paper
  also states generally that "at half coverage … the HF molecule is dissociated and the
  fluorine of the adsorbate is strongly bound to the surface", lowering E_ads by ≈0.5 eV.
- **[I]** That general sentence is in tension with the (111)/(110) table values, which stay at
  1.00–1.07 Å. Read literally, only (100) shows FHF^- formation. Flag this ambiguity if the
  model needs a single number for CaF2(111): molecular, weakly chemisorbed, no dissociation.
- **[S]** No barrier for HF dissociation on CaF2(111) is reported; only adsorption energies.
- **Not found:** any experimental determination (TPD, LEED, IRAS, XPS) of molecular vs
  dissociative HF on CaF2(111). Searches returned only DFT and CaF2 surface-structure papers.

---

## 3. Byproduct of HF + a metal amidinate, and its volatility

- **[S]** Protonating an amidinate ligand gives the free **amidine** (amidines are the
  conjugate acids of amidinates). Zaera/Gordon frame the ALD second half-cycle exactly this
  way: the amidinate ligands should "survive the first half of the cycle so they can be
  hydrogenated to the corresponding amidine in the second half" — Ma, Guo, Gordon, Zaera,
  *Chem. Mater.* 2011, 23, 3325–3334, DOI 10.1021/cm200432t.
- **[I]** For Ca(amd)2 = bis(N,N'-di-i-propylformamidinato)calcium(II), the byproduct is
  therefore **N,N'-di-i-propylformamidine** (C7H16N2, MW 128.22; PubChem CID 437419,
  IUPAC N,N'-di(propan-2-yl)methanimidamide). The stoichiometry follows the George-group
  pattern: Ca(amd)2 + 2 HF → CaF2 + 2 amidine.
- **[S]** Volatility evidence (indirect): in TPD from Ni(110), **intact
  N-sec-butylacetamidine (MW 114) desorbs molecularly as a sharp peak at 300 K**, while the
  copper amidinate itself desorbs at ~215 K (Ma et al. 2011). The same paper notes that
  released acetamidine can **re-adsorb** on SiO2 and is a suspected carbon-contamination
  source — i.e. the amidine is volatile but "sticky" enough to readsorb.
- **[I]** A 300 K desorption temperature on a metal surface implies the amidine is well
  below its own sublimation limit at 175–250 °C, so thermodynamic desorption is not the
  limiting step in CaF2 ALD. Readsorption/purge efficiency, not vapour pressure, is the
  plausible risk.
- **Not found:** any measured vapour pressure or boiling point for
  N,N'-di-i-propylformamidine (PubChem has no experimental BP/VP entry for CID 437419);
  any mass-spectrometric detection of the amidine in a CaF2 (or any fluoride) ALD reactor
  exhaust; any TPD/IR study of an amidine on a fluoride surface.

---

## 4. Activation energies (measured or computed) for HF + metal-organic surface reactions

What exists:

| Reaction / system | Value | Type | Source |
|---|---|---|---|
| HF fluorination on amorphous Al2O3 (releasing H2O and CH4), 250 °C | **E_a = 0.18–1.14 eV** | DFT (AIMD-generated slabs) | Khumaini, Cho, Kim, Lee, *Surf. Interfaces* 2025, 72, 107114, DOI 10.1016/j.surfin.2025.107114 |
| AlMe3 removal step (AlF3 → AlF(3−x)Mex, Al2F2Me4 release), 250 °C | **E_a < 1.08 eV** | same DFT study | ibid. |
| H2O formation from HF on θ/amorphous-alumina-type slabs | **barrier ≥ +0.9 eV** (endoergic); H2O desorption as cheap as +0.2 eV; AlF3 desorption unfavourable | DFT | Kondati Natarajan & Elliott, *Chem. Mater.* 2018, 30, 5912–5922, DOI 10.1021/acs.chemmater.8b01930 |
| HF on CaF2(111) | E_ads = −0.40 / −0.48 / −0.62 eV at 100/50/25 % coverage; **no barrier reported** | periodic DFT | Kaawar & Paulus 2019 (above) |
| HF on β-AlF3(100) | record verified; abstract not retrievable | DFT | Bailey, Wander, Mukhopadhyay, Searle, Harrison, *Phys. Chem. Chem. Phys.* 2008, 10, 2918, DOI 10.1039/b718733e |

What does **not** exist (searched, not found):

- No measured activation energy for any HF + metal-organic (amidinate, amido, alkyl)
  surface reaction — no Arrhenius analysis of an HF half-cycle in any ALD paper retrieved.
  The George-group kinetics work reports pressure/time dependence (Langmuir + parabolic law)
  but no E_a: Cano, Marquardt, DuMont, George, *J. Phys. Chem. C* 2019, 123, 10346–10355,
  DOI 10.1021/acs.jpcc.9b00124.
- No computed barrier for HF + Ca(amd)2, or for HF + any amidinate/amidinate-derived surface.
- No computed barrier for HF dissociation on CaF2(111)/(110)/(100).

**[I]** The 0.18–1.14 eV HF-fluorination window and the ≥0.9 eV H2O-formation barrier are the
only defensible numbers to carry into a CaF2 model, and both are for alumina, not CaF2.
Any barrier for HF + Ca(amd)2 in the model must be labelled an assumption, not a literature value.

---

## 5. CaF2 ALD with an amidinate + HF — what the record says

- **[S]** Hennessy, Rodríguez, West, Brandon, *J. Vac. Sci. Technol. A* 2023, 41 (4),
  DOI 10.1116/6.0002563: CaF2 ALD from **bis(N,N-di-i-propylformamidinato)calcium(II)**
  (= Ca(amd)2) and **anhydrous HF**; linear growth 0.4 → 0.3 Å/cycle over 175–250 °C;
  nanocrystalline films transparent to 193 nm; Si:CaF2 etch ratio > 2000:1 in SF6 plasma.
  *Full text is paywalled and was not retrieved; the abstract states no mechanism,
  no byproduct identification and no QCM/QMS data.*
- **[S]** Kim, Yang, Powers, Davis, Lou, Gordon, *Angew. Chem. Int. Ed.* 2016, 55, 10228–10233,
  DOI 10.1002/anie.201602406: synthesis and structures of bis(N,N'-di-i-propylformamidinato)-
  and bis(N,N'-di-i-propylacetamidinato)calcium(II) (dimeric in solution and solid state,
  [(η2-L)Ca(μ-η2:η2-L)(μ-η2:η1-L)Ca(η2-L)]); "unprecedentedly high volatility" compared with
  [Ca3(tmhd)6]; CaS ALD demonstrated with H2S.
- **[S]** Pilvi, Arstila, Leskelä, Ritala, *Chem. Mater.* 2007, 19, 3387–3392,
  DOI 10.1021/cm0629412: CaF2 ALD using Ca(thd)2 + TiF4; states that only one earlier CaF2 ALD
  study existed, which used **Ca(thd)2 and HF**.
- **[S]** Patent record (not peer reviewed, no DOI): US 11,111,578, Elam, Mane, Gebhard
  (UChicago Argonne), "Atomic layer deposition of fluoride thin films" — CaF2 by ALD half-cycle
  exposures of an alkaline-metal amidinate and a fluorinated compound (OSTI record 1840416).

---

## 6. Verified bibliographic records (Crossref-resolved 2026-11-13)

1. 10.1021/acs.jpcc.5b02625 — Lee, Y.; DuMont, J. W.; Cavanagh, A. S.; George, S. M.
   "Atomic Layer Deposition of AlF3 Using Trimethylaluminum and Hydrogen Fluoride."
   *J. Phys. Chem. C* **2015**, 119 (25), 14185–14194.
2. 10.1021/acs.chemmater.5b04360 — Lee, Y.; Sun, H.; Young, M. J.; George, S. M.
   "Atomic Layer Deposition of Metal Fluorides Using HF–Pyridine as the Fluorine Precursor."
   *Chem. Mater.* **2016**, 28 (7), 2022–2032.
3. 10.1116/6.0002563 — Hennessy, J.; Rodríguez, R. E.; West, W. C.; Brandon, E. J.
   "Atomic layer deposition of calcium fluoride for barrier coating applications."
   *J. Vac. Sci. Technol. A* **2023**, 41 (4) (article number not returned by Crossref/OpenAlex).
4. 10.1016/j.jfluchem.2018.09.006 — Kaawar, Z.; Paulus, B. "Adsorption of hydrogen fluoride on
   alkaline earth fluoride surfaces: A first-principles study." *J. Fluorine Chem.* **2019**,
   224, 67–72.
5. 10.1021/cm200432t — Ma, Q.; Guo, H.; Gordon, R. G.; Zaera, F. "Surface Chemistry of
   Copper(I) Acetamidinates in Connection with Atomic Layer Deposition (ALD) Processes."
   *Chem. Mater.* **2011**, 23 (14), 3325–3334.
6. 10.1021/acs.chemmater.8b01930 — Kondati Natarajan, S.; Elliott, S. D. "Modeling the Chemical
   Mechanism of the Thermal Atomic Layer Etch of Aluminum Oxide: A Density Functional Theory
   Study of Reactions during HF Exposure." *Chem. Mater.* **2018**, 30 (17), 5912–5922.
7. 10.1016/j.surfin.2025.107114 — Khumaini, K.; Cho, G.; Kim, H.-L.; Lee, W.-J. "Thermal atomic
   layer etching mechanism of amorphous aluminum oxide by hydrogen fluoride and
   trimethylaluminum: First-principles study." *Surf. Interfaces* **2025**, 72, 107114.
8. 10.1021/acs.jpcc.2c00158 — Hu, X.; Schuster, J. "Chemical Mechanism of AlF3 Etching during
   AlMe3 Exposure: A Thermodynamic and DFT Study." *J. Phys. Chem. C* **2022**, 126, 7410–7420.
9. 10.1021/acs.jpcc.9b00124 — Cano, A. M.; Marquardt, A. E.; DuMont, J. W.; George, S. M.
   "Effect of HF Pressure on Thermal Al2O3 Atomic Layer Etch Rates and Al2O3 Fluorination."
   *J. Phys. Chem. C* **2019**, 123 (16), 10346–10355.
10. 10.1002/anie.201602406 — Kim, S. B.; Yang, C.; Powers, T.; Davis, L. M.; Lou, X.;
    Gordon, R. G. "Synthesis of Calcium(II) Amidinate Precursors for Atomic Layer Deposition
    through a Redox Reaction between Calcium and Amidines." *Angew. Chem. Int. Ed.* **2016**,
    55 (35), 10228–10233.
11. 10.1021/cm0629412 — Pilvi, T.; Arstila, K.; Leskelä, M.; Ritala, M. "Novel ALD Process for
    Depositing CaF2 Thin Films." *Chem. Mater.* **2007**, 19 (14), 3387–3392.
12. 10.1039/b718733e — Bailey, C. L.; Wander, A.; Mukhopadhyay, S.; Searle, B. G.;
    Harrison, N. M. "Adsorption of HF and HCl on the β-AlF3(100) surface."
    *Phys. Chem. Chem. Phys.* **2008**, 10 (20), 2918.
13. 10.1021/acs.chemmater.1c01950 — Lii-Rosales, A.; Cavanagh, A. S.; Fischer, A.; Lill, T.;
    George, S. M. "Spontaneous Etching of Metal Fluorides Using Ligand-Exchange Reactions:
    Landscape Revealed by Mass Spectrometry." *Chem. Mater.* **2021**, 33 (19), 7719–7730.
    (Record verified; full text not retrieved — content not used here.)
14. 10.1002/admi.202400372 — Atosuo, E.; Mäntymäki, M.; Ritala, M. "ALD of Metal Fluorides–
    Potential Applications and Current State." *Adv. Mater. Interfaces* **2024**, 12 (4).
    (Record verified; only abstract/early sections retrieved.)

---

## 7. Explicit unknowns for the CaF2 model

1. No in situ QCM/QMS/IR study of Ca(amd)2 + HF exists in the retrieved literature; the
   QCM "x" (adsorbed HF per CaF2) for CaF2 is unmeasured.
2. No experimental molecular-vs-dissociative determination for HF on CaF2(111).
3. No measured E_a for HF + Ca(amd)2 (or any HF + amidinate/amido/alkyl surface step).
4. No experimental vapour pressure for N,N'-di-i-propylformamidine; volatility inferred only
   from a 300 K TPD desorption of a different (sec-butyl) amidine on Ni(110).
5. No mass-spectrometric identification of the amidine byproduct in any fluoride ALD reactor.
6. The CaF2 ALD paper's full text (possible mechanism discussion) is paywalled and was not
   read; nothing in its abstract should be quoted as a mechanism.
