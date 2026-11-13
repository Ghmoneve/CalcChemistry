# CaF2(111) 表面结构数据汇编（用于 ALD 半反应 DFT slab 建模）

所有数值均来自实际检索到的原文（正文、摘要或 ICSD 条目），来源列于 [S#]。凡未能检索到的，明确标注「未找到」。
自算值一律标注「derived」，并给出所用公式，不冒充文献值。

---

## 1. 体相 CaF2 晶格常数（fluorite, Fm-3m, No. 225）

| 量 | 值 | 说明 |
|---|---|---|
| a | **5.4626 Å** | 室温 XRD 粉末数据；Z = 4；V = 163.001 Å³ [S1] |
| 空间群 | Fm-3m (225)，Wyckoff: Ca 4a (0,0,0)，F 8c (¼,¼,¼) | [S1] |
| Ca–F 键长 | 2.365 Å (derived: a√3/4) | 与文献常引的 2.39 Å 差异源于取整 |

[S1] H. E. Swanson, E. Tatge, *Standard X-ray Diffraction Powder Patterns*, National Bureau of Standards (U.S.) Circular 539, 1953, pp. 69–70.（数字化记录：ICSD 52754；CIF 见 https://topas.webspace.durham.ac.uk/wp-content/uploads/sites/261/2021/08/caf2_52754.cif，其中 `_cell_length_a 5.4626`、`_symmetry_space_group_name_H-M 'F m -3 m'`）

注：同一 5.4626 Å 值亦被矿物学手册 (handbookofmineralogy.org, fluorite) 与 Wikipedia fluorite 条目重复引用，属同一原始数据。

---

## 2. 表面终止与 F–Ca–F triple layer 间距

**结论：确认为 F 终止的 (1×1) 表面，无重构。**

依据（均为检索到的原文表述）：
- 「The perfect (111) cleavage plane is F⁻ terminated and the topmost F⁻ ions form a hexagonal array.」[S2]
- 「The (111) surface is thus terminated with fluoride ions and seven-coordinated metal ions occur.」[S3]
- 「the F-M-F stoichiometric termination」；文中定义 ideal surface 为 bulk-truncated、无重构、仅有微小原子位移 [S4]
- 「charged planes are arranged into neutral repeating unit (F-Ca-F) terminated with fluoride ion F⁻ on the surface」[S5]
- LEED 分析基于「F–M–F layer structure」，拟合量即最外层 fluoride 层位移 [S6]

**triple layer 间距：**

| 量 | 值 | 来源 |
|---|---|---|
| F–Ca–F trilayer (TL) 高度 | **3.15 Å** | 「The basic unit of the fluorite structure of CaF₂ along the (111) direction is a F-Ca-F trilayer (TL), 3.15 Å high.」[S7] |
| 台阶高度 = 1 TL | **0.32 nm** (≈3.2 Å) | [S2] |
| d(111) = a/√3 | **3.1538 Å** (derived, a = 5.4626 Å) | 自算，与上述实验值一致 |

子层结构（derived，理想 fluorite 位置，a = 5.4626 Å，沿 [111] 投影）：
- 每个 3.154 Å 周期内含 1 个 Ca 面 + 2 个 F 面；
- F(0) – **0.7885 Å** – Ca – **0.7885 Å** – F(1.577 Å)：即 F–Ca–F 板厚 1.577 Å（= a/(2√3)）；
- 相邻 triple layer 之间空隙 **1.577 Å**（= a/(2√3)），即未占据的 octahedral 位所在处；
- 面内 F–F 最近邻 = a/√2 = **3.8626 Å**，与 AFM 文献所给 386.2 pm 一致 [S8]。

[S2] V. E. Puchin, A. V. Puchina, M. Huisinga, M. Reichling, "Theoretical modelling of steps on the CaF2(111) surface", *J. Phys.: Condens. Matter* **13** (2001) 2081–2094. DOI: 10.1088/0953-8984/13/10/304
[S3] Z. Kaawar, S. Mahn, E. Kemnitz, B. Paulus, "On the Morphology of Group II Metal Fluoride Nanocrystals at Finite Temperature and Partial Pressure of HF", *Molecules* **22** (2017) 663. DOI: 10.3390/molecules22040663
[S4] A. S. Foster, T. Trevethan, A. L. Shluger, "Structure and diffusion of intrinsic defects, adsorbed hydrogen, and water molecules at the surface of alkali-earth fluorides calculated using density functional theory", *Phys. Rev. B* **80** (2009) 115421. DOI: 10.1103/PhysRevB.80.115421
[S5] "Type B epitaxy of Ge on CaF2(111) surface", *Surface Science* (2010), DOI: 10.1016/j.susc.2010.06.008（仅摘要检索到，作者名未获取，见第 7 节）
[S6] J. Vogt, J. Henning, H. Weiss, "The structure of CaF2(111) and BaF2(111) single crystal surfaces: A tensor low energy electron diffraction study", *Surface Science* **578** (2005) 57–70. DOI: 10.1016/j.susc.2005.01.011
[S7] F. Calleja, J. J. Hinarejos, A. L. Vázquez de Parga, S. M. Suturin, N. S. Sokolov, R. Miranda, "Epitaxial growth of CaF2(111) on Cu(111) visualized by STM", *Surface Science* **582** (2005) 14–20. DOI: 10.1016/j.susc.2005.03.003
[S8] F. J. Giessibl, M. Reichling, "Investigating Atomic Details of the CaF2(111) Surface with a qPlus Sensor", arXiv:cond-mat/0412430 (2004)；期刊版 DOI 10.1088/0957-4484/16/3/022（原文：「the corresponding surface layers are trigonal arrangements of F⁻ ions spaced by a0/√2 = 386.2 pm… the CaF2(111) surface must be terminated by a complete triple layer F⁻–Ca²⁺–F⁻ with a F⁻ layer at the surface」）

---

## 3. 弛豫 / 重构数据（LEED、MD、DFT）

| 方法 | 结果 | 文献 |
|---|---|---|
| Tensor LEED, 20 K, 80–360 eV | 与 truncated bulk 几乎无差别；最大位移为**最外层 F⁻ 向内 0.03 ± 0.06 Å**；Ca 与第二 F 子层、以及更深层均在体相位置 | [S6] |
| Tensor LEED + MD, 300 K | 结构接近 bulk-terminated，仅最外离子小幅弛豫；最外 F–Ca–F triple layer 的**振动幅度增强** | [S9] |
| DFT (VASP, GGA-PAW) | ΔF¹ = +0.01 Å, ΔM = +0.01 Å, ΔF² = −0.02 Å（相对理想终止） | [S4, Table I] |
| Hybrid B3PW (CRYSTAL) | 顶层 upper F 向内弛豫 **0.05 Å** | [S10] |
| Ab initio HF (CRYSTAL95) | 位移 < 0.01 Å；6/9/12 层弛豫量 −1.0/−0.7/−0.2 pm | [S2, Table 1] |
| 全局结构搜索 (DFT + minima hopping) | 重构只出现在 **(100)** 表面；(111) 表面稳定、部分疏水，未报告 (111) 重构 | [S11] |

结论：CaF2(111) 在真空中是 **(1×1)、几乎无弛豫**的表面，位移量级 ≤0.05 Å，且小于 LEED 的误差棒。

[S9] A. Gotte, M. Baudin, A. G. Cabello-Cartagena, J. Vogt, H. Weiss, "Theoretical and experimental studies of the structure and dynamics of the CaF2(111) surface", *Surface Science* **601** (2007) 411–418. DOI: 10.1016/j.susc.2006.10.023
[S10] H. Shi, R. I. Eglitis, G. Borstel, "Ab initio calculations of the CaF2 electronic structure and F centers", *Phys. Rev. B* **72** (2005) 045109. DOI: 10.1103/PhysRevB.72.045109
[S11] G. Fisicaro, M. Sicher, M. Amsler, S. Saha, L. Genovese, S. Goedecker, "Surface reconstruction of fluorites in vacuum and aqueous environment", *Phys. Rev. Materials* **1** (2017) 033609. DOI: 10.1103/PhysRevMaterials.1.033609

---

## 4. 表面能 γ(111)

| 方法 | γ(111) | γ(110) | γ(100) | 文献 |
|---|---|---|---|---|
| Ab initio HF（6/9/12 原子层） | **0.466–0.471 J/m²** | 0.82–0.86 | 1.12–1.19 | [S2, Table 1] |
| Hybrid B3PW（7 层 / 9 层） | **0.437 / 0.438 J/m²** | 0.717 / 0.719 | 0.957 / 0.979 | [S10, Table VIII] |
| PBE 周期 slab | **0.47 J/m²** | 0.71 | 0.95 | [S3] |
| 实验（cleavage 估计） | **0.45 J/m²** | — | — | 作为二手引用见于 [S3] 与 [S10]；原始文献未能检索到（见第 7 节） |

台阶能（(111) 上两类 [110] 取向台阶）：0.25 与 0.35 nJ/m；[110] 双高台阶在能量上有利 [S2]。

一致性：三种独立方法给出 γ(111) ≈ 0.44–0.47 J/m²，且一致认为 (111) 是最稳定的低指数面 → 建模时用 (111) 终止是合理的。

---

## 5. 已发表的 CaF2(111) DFT / 第一性原理 slab 模型

| 文献 | 方法 / 泛函 | Slab 厚度 | 真空层 | 其他设置 |
|---|---|---|---|---|
| [S2] Puchin 2001 | Ab initio **Hartree–Fock**, CRYSTAL95, Gaussian 基组 (8-6-5-11G* for Ca²⁺, 7-311G for F⁻) | 6 / 9 / 12 原子层（≈480 / 800 / 1120 pm） | 正文未给出数值 | a 固定为 HF 值 0.5538 nm；平移矢量固定 |
| [S10] Shi 2005 | **B3PW** hybrid（另有 HF、LDA、PW91、PBE 用于体相带隙对比） | 7 层与 9 层（(110)/(100) 用 15 层） | 正文未给出数值 | 局域 Gaussian 基组；表面能对层数已收敛 |
| [S4] Foster 2009 | **VASP, PAW, GGA**；500 eV cutoff | **3 个 F–Ca–F triple layer**（2×2 表面胞） | **1.7 nm (17 Å)** | 5×5×5 MP k-mesh（体相 2×2×2 超胞）；monopole/dipole 校正，ε = 6.76；NEB 迁移路径 |
| [S12] Barzilai 2008 | DFT, **APW+lo (WIEN2k)**, 全势 | 正文未给出数值 | 有 vacuum separation，数值未获取 | 3 个吸附位点；金属层距表面 1–4 Å；结论：仅在 terminal F 之上吸附才放热 |
| [S3] Kaawar 2017 | 周期 slab DFT，**PBE**（体相亦对比 B3LYP） | 正文未获取 | 未获取 | 表面胞固定为体相值，仅弛豫原子 |

**建模建议（基于上述数据的推论，非文献结论）**：以 [S4] 的设置（3 TL + ≥15 Å 真空 + 偶极校正）为下限做收敛测试；由于 (111) 面弛豫极小，可直接用 bulk-truncated 结构起算而误差在 0.05 Å 量级。

[S12] S. Barzilai, N. Argaman, N. Froumin, D. Fuks, N. Frage, "First-principles modeling of metal layer adsorption on CaF2(111)", *Surface Science* **602** (2008) 1517–1524. DOI: 10.1016/j.susc.2008.02.022

---

## 6. 带隙与 Si 晶格失配

### 带隙
| 量 | 值 | 文献 |
|---|---|---|
| 直接带隙 (Γ→Γ)，实验 | **12.1 eV** | [S10]（引言）；[S13] 亦引此值 |
| 间接带隙 (X→Γ)，实验估计 | **~11.8 eV** | [S13] |
| DFT-LDA | 8.72 eV | [S10] |
| DFT-PW91 / PBE | 8.51 / 8.45 eV | [S10] |
| DFT-GGA (VASP-PAW) | 8.0 eV | [S4] |
| B3PW hybrid | 10.96 eV | [S10] |
| G0W0@LDA | 11.38 eV | 作为引用值见 [S13] |
| HF (ab initio) | 20.77 eV（严重高估） | [S10] |

**直接后果**：GGA 低估 ~4 eV，不能用于判断 ALD 半反应的电荷转移/带边对齐；电子结构必须用 hybrid 或 GW，几何与吸附能可用 GGA（[S4] 自评误差 ~0.1 eV）。

[S13] F. Matusalem, M. Marques, L. K. Teles, A. Filippetti, G. Cappellini, "Electronic properties of fluorides by efficient approximated quasiparticle DFT-1/2 and PSIC methods: BaF2, CaF2 and CdF2 as test cases", *J. Phys.: Condens. Matter* **30** (2018) 365501. DOI: 10.1088/1361-648X/aad654

（另：Y. Ma, M. Rohlfing, "Quasiparticle band structure and optical spectrum of CaF2", *Phys. Rev. B* **75** (2007) 205114, DOI 10.1103/PhysRevB.75.205114 — 仅检索到题录与摘要片段「direct band gap is 12.1 eV and the indirect gap is estimated to be 11.8 eV」，全文受订阅墙限制，故数值以上述 [S13] 为准。）

### 与 Si 的晶格失配
| 量 | 值 | 文献 |
|---|---|---|
| a(Si), 300 K | **5.43095 Å** | [S14]（参考表，300 K） |
| 失配 (a_CaF2 − a_Si)/a_Si | **+0.583 %**（derived） | 自算：5.4626 vs 5.43095 |
| 文献值 | **+0.6 %**，室温 | [S15]（同时讨论 Si(100) 与 Si(111)） |

- 对 **Si(100)**：cube-on-cube，(100) 面内 a_CaF2 ∥ a_Si，失配同为 +0.58%。
- 对 Si(111)：CaF2 与 Si 的面内周期均为 a/√2（3.8626 vs 3.8401 Å），失配同样是 +0.58% → 两种取向的失配数值相同，来源都是体相 a 的 0.58%。
- 后果：外延 CaF2/Si 呈压缩应变（CaF2 较大），~0.58% 失配在几 nm 厚度内可靠弹性容纳；若 slab 用 Si 衬底建模，需注意此应变。

[S14] Advanced Photon Source, Sector 7-ID, "Crystal Lattice Parameters"（300 K 表，Si = 5.43095 Å）, https://7id.xray.aps.anl.gov/calculators/crystal_lattice_parameters.html
[S15] M. Maeda, N. Matsudo, S. Watanabe, K. Tsutsui, "Crystalline structure of epitaxial CaxMg1−xF2 alloys on Si(100) and (111) substrates", *Thin Solid Films* **515** (2006) 448–451. DOI: 10.1016/j.tsf.2005.12.254（原文：「their lattice mismatches with Si at room temperature (RT) are +0.6%…」）

### 附：ALD 相关背景（供建模选体系参考）
T. Pilvi, K. Arstila, M. Leskelä, M. Ritala, "Novel ALD Process for Depositing CaF2 Thin Films", *Chem. Mater.* **19** (2007) 3387–3392. DOI: 10.1021/cm0629412 — Ca(thd)₂ + TiF₄，300–450 °C，生长速率 1.6 Å/cycle，薄膜为多晶（即实验 ALD 并未外延在 (111) 单晶面上）。

---

## 7. 明确未找到 / 未核实的内容

1. **实验表面能 0.45 J/m² 的原始文献**：仅作为二手引用见于 [S3]（其 ref. 27）与 [S10]（其 ref. 40），原始出处未检索、未核实。
2. **Puchin [S2]、Shi [S10]、Barzilai [S12]、Kaawar [S3] 的真空层厚度具体数值**：可获取的正文/摘要中均未给出（仅 [S4] 明确 1.7 nm）。
3. **针对 CaF2(111) 上 ALD 半反应（如 Ca(thd)₂ / TiF₄ / amidinate + fluoride 前驱体）的已发表 DFT 研究**：检索不到。检索到的 CaF2 表面 DFT 工作只涉及 H₂O/H 吸附 [S4]、金属吸附 [S12]、表面能与形貌 [S2, S3, S10]、重构搜索 [S11]。这意味着本项目的 ALD 半反应建模没有可直接对照的文献基准。
4. **CaF2(111) 上绝对层间距的实验值（Å）**：LEED 只报道相对体相的位移（±0.03–0.06 Å），未检索到绝对 z 坐标表；本报告的 0.7885/1.577 Å 子层值为理想几何推算。
5. **专论 Si(100) 上 CaF2 失配的独立测量文献**：只获取到 [S15] 的通用「+0.6%」表述（该文同时涵盖 (100) 与 (111)），未检索到单独针对 (100) 的失配测定。
6. **[S5]（Type B epitaxy of Ge on CaF2(111)）的作者名单**：仅检索到摘要与 DOI，Crossref 未查询该条，引用时作者待补。
7. **CaF2 热膨胀系数与高温晶格常数**：未检索，故本报告所有数值均按室温处理；若 slab 建模在 ALD 温度（300–450 °C）下比较晶格，需另行核实热膨胀。

---

## 8. 关键数字速查（供 slab 构建）

```
a (bulk)            = 5.4626 Å          [S1]
space group         = Fm-3m (225), Ca 4a, F 8c
surface             = (111), F-terminated, (1×1), no reconstruction   [S2,S3,S4,S6]
d(111)              = 3.1538 Å  (= a/√3, derived)
TL period / step    = 3.15 Å (exp 3.15 / 0.32 nm)   [S7,S2]
F–Ca interplanar    = 0.7885 Å  (= a/(4√3), derived)
inter-TL gap        = 1.577 Å   (= a/(2√3), derived)
in-plane F–F        = 3.8626 Å  (= a/√2; lit. 386.2 pm)  [S8]
relaxation (top F)  = −0.03 ± 0.06 Å (LEED) ; +0.01 Å (DFT)   [S6,S4]
γ(111)              = 0.44–0.47 J/m² (HF / B3PW / PBE); exp 0.45 (二手)  [S2,S10,S3]
E_g direct / indirect = 12.1 / ~11.8 eV; GGA ≈ 8.0–8.7 eV   [S13,S4,S10]
mismatch vs Si      = +0.58 % (derived) ; lit. +0.6 %   [S14,S15]
recommended slab    = ≥3 TL, ≥15 Å vacuum, dipole correction  [S4]
```
