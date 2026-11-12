# CalcChemistry：审计后状态

旧版“突破已证明”的结论已撤回，证据见 `docs/audit_numerics.md`。

## 已完成并验证

### 论文
- `paper_main.tex` + 5 个 `sections/*.tex`：19 个定理类陈述，均含完整证明。
- Tectonic 0.15.0 实际编译：19 页，无 undefined reference、无 LaTeX error、无 overfull。
- `references.bib` 8 篇与正文引用双向完全匹配（无缺失、无未引用）。
- 8 篇 DOI 元数据经 Crossref 逐条回读，字段不一致数 0；原始响应存于 `build/reference-readback.json`。

### Lean 形式化（内核检查通过）
- Lean 4.19.0 + mathlib v4.19.0（commit `c44e0c8e`）。
- `lake build`：**Build completed successfully**。
- `lake env lean scripts/CheckAxioms.lean`：14 个定理，全部只依赖 `propext`、`Classical.choice`、`Quot.sound`。
- 无 `sorry`、无 `admit`、无项目自定义公理。
- 定理与论文逐条对应（14↔14↔14）。

### 数值
- 30 项测试全部通过（警告视为错误、NumPy 浮点异常抛出）。
- 360 个随机分块、9 个 2×2 尖锐案例、20 个 CTMC 案例，重新运行并解析核对计数。
- 最大包围界违约 `4.44e-16`；最大尖锐性误差 `8.88e-16`；CTMC 无超容差违约。
- 人工有限模型浮点验证，非严格区间算术，非化学数据。

### 配对模型
- 补入 `p=q` 自配对项；独立 Jordan–Wigner/CAR 构造交叉验证。
- 参考例基态由 `-0.583095…` 更正为 `-0.883095…`。
- 明确为两轨道约化配对玩具模型，非从头算 H₂。

## 已撤回的旧声明

- SOS “0.52% gap / O(K⁸)”：来自硬编码表格，非实际求解。
- 拓扑 “0% DFT 误差”：参考值预先写死为 3；且一般 `β₀` 不等于容量（已有反例）。
- CaF₂ 的能隙、BO 常数、无锥形交叉：无计算支持。
- “85% 分子多项式可解”、流形维数定理、普遍 Yangian 证伪：未证明。
- MoCl₄ 文档：研究构想，非已验证算法或复杂度证明。
- 根目录 Richardson/Yangian 旧 Lean 草稿：含 `sorry`/错误公理，已标记归档且不参与构建。

## 化学范围

保留用户原记号 `Al2O3/MoCl4` 与 `CaF`（未擅自改写为 CaF₂）。表面、化学计量、前驱体与共反应物仍需指定。当前没有任何真实 ALD 势垒、吸附能或 GPC 数据。

## 复现

```
python -B -m unittest discover -s tests -v
python -B richardson_h2_verification.py
python -B -m experiments.certificates
python scripts/check_manuscript.py
python scripts/verify_references.py
.tools/tectonic --outdir build paper_main.tex
```

Lean/mathlib 固定 v4.19.0；Git 直连失败时用 `scripts/bootstrap_archives.py` 从官方固定提交归档恢复（已校验字节并记录 SHA256）。`proofwidgets` 缺失的 Git 标签元数据已从固定源码快照重建，以启用官方缓存流程。
