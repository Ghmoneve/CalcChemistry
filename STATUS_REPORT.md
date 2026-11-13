# CalcChemistry：当前状态

旧版"突破已证明"的结论已撤回，证据见 `docs/audit_numerics.md`。
本文档只记录**本次实际执行并核验**的结果。

---

## 一、论文

主稿 `paper_main.tex` + 6 个 `sections/*.tex`，**25 页**。

- Tectonic 0.15.0 实际编译：主文档与单文件版各 25 页，**均无 undefined reference、无致命 TeX 错误**。
- `references.bib` **26 条**，其中 **24 条经正文引用**（另 2 条未引用：Atosuo2024、Chan2010）。
- 全部新文献由 `scripts/add_refs.py` **从 Crossref 元数据程序化生成**，无一手写条目；原始响应存于 `build/added-refs-readback.json`。
- 一条命令复现两版 PDF：`python3 scripts/build_paper.py`（改 tex 后必须运行）。

### 新增章节 `sections/caf2_worked_example.tex`

给出一个完整走通的工作实例：CaF₂ 晶体沉积（`[Ca(amd)₂]₂` + 无水 HF 在 CaF₂(111) 上）。

**可核验结果：**

1. **每循环覆盖度**：由已发表 GPC（0.3–0.4 Å/cycle，Hennessy 2023）与三层板高度
   `a₀/√3 = 3.1538 Å`（由实测晶格常数 `a₀ = 5.4626 Å` 导出）得
   **9.51%–12.68%**，即每循环完成约十分之一层。
2. **每个沉积 Ca 占据 101.9–135.8 Å²**，是裸 Ca 位点面积 `12.9211 Å²` 的 **7.9–10.5 倍**。
3. **配体空间位阻设计判据**（内核验证）：足迹面积 `A₀` 把生长率封顶在 `p·A_site/A₀`。
4. **两个独立已发表数据集互相印证**：晶体结构给出每个二聚体 `938.5 Å³`，
   隐含密度 `1.042 g/cm³`（有机金属的合理值），等效球截面 **115.9 Å²**——
   正落在 GPC 反推所需的窗口内。两者都约为裸位点面积的 9 倍。
5. **采样无法给出势垒上界**（内核验证的反例定理），只能给下界。

**明确未主张**：本体系没有任何势垒数值。文中明确写出三条"一致性不等于证明"的理由
（动力学限制未被排除、100 K 晶体结构 vs 175–250 °C 工艺、模型中配体被截断），
并给出证伪条件。

---

## 二、Lean 形式化（内核检查通过）

- Lean 4.19.0 + mathlib v4.19.0（commit `c44e0c8e`）。
- `lake build`：**Build completed successfully**。
- **33 个声明**，全部经 `#print axioms` 检查，**只依赖 `propext`、`Classical.choice`、`Quot.sound`**。
- 全仓 Lean 源码**无 `sorry`、无 `admit`、无自定义公理**。

文件与内容：

| 文件 | 内容 |
|---|---|
| `CalcChemistry/CertifiedTheory.lean` | 区间减除、误差传播、尖锐分块标量界、速率区间、松弛下界、配对多项式 |
| `CalcChemistry/BarrierBounds.lean` | **分块能量证书推广到任意有限维 n×m**；无平方根下折形式 `a - b²/(c-a)`；单位向量形式 |
| `CalcChemistry/SurfaceReaction.lean` | 覆盖度区间、吸附能区间、**配体足迹判据**、**证伪条件**、两层势垒区间、**采样无上界反例** |
| `CalcChemistry/BornOppenheimer.lean` | 质量参数非负（仅算术，不含 BO 定理） |

论文 §Evidence and formal verification status 逐条映射全部 33 个声明，并明确写出：
**数量是声明数而非科学内容**，除分块证书与配体不等式外均为实算术管线，无一条认证分子能量。

---

## 三、数值计算

### 已完成

| 项目 | 结果 | 性质 |
|---|---|---|
| 每循环覆盖度换算 | 9.51%–12.68% | 已发表数据的算术，可复算 |
| Ca 位点几何 | 12.9211 Å²/位点 | 由实测晶格常数导出 |
| 立体足迹交叉核验 | 115.9 Å² vs 需求 ≤135.8 Å² | 两个独立数据集互证 |
| pyscf 环境探测 | pyscf 2.15.0 可用；HF/def2-svp 基态 `-99.9325 Ha` | 自检 |

### 基准测试（诊断结果，非科学结论）

| 设置 | 结果 |
|---|---|
| mesh30 + 1×1×1 k + density fitting | **收敛，49 s** |
| mesh40 + 2×2×2 k + density fitting | 546 s，**未收敛** |
| mesh40 + 2×2×2 k **无 DF** | 532 s，**E = 20009 Ha，完全发散** |

结论：周期性 DFT 必须开启 density fitting；k 点数是主要成本。FFT 网格若不显式设置会
自动选到 1301³（约 49 GB），直接 OOM。

### 进行中

- 气相模型反应 `Ca(formamidinate)₂ + 2 HF → CaF₂ + 2 formamidine`
  （PBE/B3LYP × def2-svp/def2-tzvp，BFGS 几何优化）。
  电子数已校验配平（86 = 86）。**结果尚未产生。**

---

## 四、已撤回的旧声明

- SOS "0.52% gap / O(K⁸)"：来自硬编码表格，非实际求解。
- 拓扑 "0% DFT 误差"：参考值预先写死为 3；且一般 `β₀` 不等于容量（已有反例）。
- CaF₂ 的能隙、BO 常数、无锥形交叉：无计算支持。
- "85% 分子多项式可解"、流形维数定理、普遍 Yangian 证伪：未证明。
- MoCl₄ 文档：研究构想，非已验证算法或复杂度证明。
- 根目录 Richardson/Yangian 旧 Lean 草稿：含 `sorry`/错误公理，不参与构建。
- `drafts/BarrierBounds.lean` 初版含三处证明错误，修正并自行复核后才入库。

---

## 五、化学范围

用户已确认：目标为 `Al2O3/MoCl5` 与 **CaF₂ 晶体沉积**
（`MoCl4` 与 `CaF` 均为用户本人笔误，已按其确认更正）。

已确定的建模选择（依据与 DOI 见 `docs/chemistry_decisions.md`）：

- **表面**：CaF₂(111)，F 终止，(1×1) 无重构。
- **前驱体**：`[Ca(amd)₂]₂`（二聚体，Ca···Ca 3.323 Å）+ 无水 HF。
- **HF 吸附态**：分子态，`E_ads = −0.40/−0.48/−0.62 eV`（100/50/25% 覆盖度）。
- **机理框架**：HF 脉冲既是质子化剂，又留下分子态 HF 作为表面氟储库。

**仍未知**：CaF₂ 工艺的实测势垒、QCM 质量平衡、副产物检测；
`Al2O3/MoCl5` 的晶面与羟基覆盖度。

---

## 六、复现

```sh
python3 scripts/build_paper.py           # 论文双版本（改 tex 后必跑）
python3 scripts/add_refs.py references.bib <DOI> ...   # 程序化加文献
python3 scripts/check_manuscript.py      # 静态检查
python3 scripts/verify.py                # 六门统一验证（旧门）
export PATH="$PWD/.tools/lean-4.19.0/bin:$PATH"
lake build                               # Lean 全量构建
lake env lean scripts/CheckAxioms.lean   # 33 个声明的公理审计
python3 experiments/caf2/gpc_coverage.py        # 覆盖度换算
python3 experiments/caf2/steric_crosscheck.py   # 立体足迹交叉核验
```

Lean/mathlib 固定 v4.19.0；Git 直连失败时用 `scripts/bootstrap_archives.py`
从官方固定提交归档恢复。tectonic 必须加 `--keep-intermediates`，否则不落盘
`.bbl`，单文件版会静默套用旧文献。
