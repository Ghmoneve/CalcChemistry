# 下次会话从哪里开始

写于 2026-11-17 会话结束时。HEAD `52ea5d8`，工作区干净，六个验证门全过。

**触发约定：用户在对话中发"开始"，就执行下面"待办"里指定的那一项，不要重新讨论，不要再问一遍。**

---

## 已完成：A-内积 Cauchy–Schwarz（2026-11-18）

方案 A 做成了。`Mathlib.Tactic` 编译完成，`CalcChemistry/PSDCauchy.lean` 四条定理全部进内核：
`disc_le`（标量判别式）、`posSemidef_symm_apply`（对称性）、`posSemidef_cauchy`（A-内积 C-S）、
`schur_scalar_hypothesis`（实例化到 s = bᵀA⁻¹b）。60→61 条声明，全部仅三标准公理，无 sorry。
等式在 w = A⁻¹b 处精确取到（比值 1.000000000000），确认陈述 sharp。

---

## 原待办记录（已解决，保留供追溯）

### 缺口是什么

`CalcChemistry/SchurStep.lean` 已经证明：**给定**逐 $w$ 的标量不等式
$(b^\top w)^2 \le s(\langle Aw,w\rangle-\lambda\|w\|^2)$，则块二次型有下界
$(\lambda-s)(x^2+\|w\|^2)$。

但**如何从 $A\succ0$ 得到那个不等式**这一步没有形式化。数学上它就是 $A$-内积下的
Cauchy–Schwarz：

$$(b^\top w)^2 \le (b^\top A^{-1}b)\,(w^\top A w)$$

- 论文里是 Proposition 1 的证明（三行，已写）
- 数值上 `experiments/complement/verified_cholesky.py` 已核验
- **内核看不到**

### 三个方案（用户要选一个）

**A. 做完整**——预期 1.5–3 小时，可能失败
1. 编译 `Mathlib.Tactic`（当前无 olean；`PosDef` 那次拉了 2238 模块，这次还要一批）
2. 证 PSD 二次型的判别式论证：`∀t, a + 2tb + t²c ≥ 0 ⟹ b² ≤ ac`
3. 接到 `SchurStep.schur_block_bound_shifted` 的假设上

**B. 不做，维持现状的写法**——0 成本
现在的论文措辞已在学术上成立：命题有证明、有数值核验、并明确写出"内核看不到这一步"。
位置：`paper_main.tex` 的 Lean mapping 一节，"One step of the chain is *not*
formalized..." 那段。

**C. 只编译 `Mathlib.Tactic`**——20–40 分钟
低成本前置投资，之后无论 A 做不做都用得上。做完不改变任何论断。

### 已确认的技术事实（省得明天重查）

- `Matrix.PosDef` / `Matrix.Spectrum` 的 olean **已经编好**：
  `.lake/packages/mathlib/.lake/build/lib/lean/Mathlib/LinearAlgebra/Matrix/`
- `Mathlib.Tactic` 的 olean **仍然缺失**
- 可用的相关引理：
  - `Matrix.PosSemidef.re_dotProduct_nonneg : 0 ≤ re ⟨x, M.mulVec x⟩`
  - `Matrix.PosDef.re_dotProduct_pos`
  - `Matrix.PosDef.posSemidef`
  - `Matrix.PosSemidef.isHermitian`
- **没有现成的 A-内积 Cauchy–Schwarz**，必须自己证
- 上次失败点：`Finset` 求和展开二阶型 `Q(u,v) = ∑ i ∑ j, u i * A i j * v j` 时卡住

---

## 当前项目状态（供参考，不需要重做）

| 项 | 状态 |
|---|---|
| 论文 | 34 页，两版源码 0 undefined / 0 fatal |
| Lean | 57 条声明，9 模块，仅 `propext, Classical.choice, Quot.sound`，无 sorry |
| 参考文献 | 34 条，全部经 Crossref 核验，全部被正文引用 |
| 验证门 | python-tests / pairing / lean-build / lean-axioms / manuscript-static / latex-build 全 exit 0 |
| PDF | `/mnt/e/CalcChemistry/CalcChemistry_paper.pdf` |
| 单文件源 | `CalcChemistry_paper_standalone.tex` |

### 本会话已完成的三件事

1. **理论**（`84531dd`）：Schur 证书此前**没有任何定理**。补了 Prop 1（Schur 下界）、
   Lemma 2（验证式 Cholesky）、Lemma 3（标量夹逼）+ Remark（何时优于块界），全部带证明，
   全部数值核验。
2. **CaF₂**（`c66166d`）：CCSD(T) 实测势垒 0.00–0.11 kcal/mol（文献 0.15），
   角度 155°（文献 154.3°）；DFT 给 1.1–2.2，**高估 20 倍**。因冻结核敏感性与势垒同量级，
   论文只声称**上界 ≤0.15**，不假装知道精确值。
3. **代码**（`52ea5d8`）：`reaction3.py` 的 `linear=True` 分支 5 条 Pyright 告警清零，
   行为验证不变（|dE| = 6.8e-13 Ha）。

### 更早的重要结论（不要推翻）

- **块证书路线已封闭**：对称 $O$ 的相似变换保谱，$\inf_S\|SOS^{-1}\|_2=\|O\|_2$；
  符号、低秩、对角加权都逃不掉。实测 20000 次随机缩放，最小值 59.33，无一低于 2.026。
- **Schur 路线严格认证到 0.860 kcal/mol**（对 STO-3G 水，E_FCI 差 0.860），代价 $O(n^3)$。
  质量提升，不是复杂度提升；指数墙未动。
- **标量 $b^\top A^{-1}b$ 无法便宜地界定** —— 这是真正的开放问题。
  两个旋钮（shift $t$、参考块大小 $k$）都朝错误方向转，已测。
- **CaF₂ 是 quasilinear**，不是线性。之前论文写"惯例为线性"是错的，已撤回。
- 反应能是**区间**不是数：八级别 $\times$ 两种 CaF₂ 处理，$[-11.95,-5.69]$，
  可接受级别上 $-6$ 到 $-8$。spread 是**不确定度下界**，不是估计值。

---

## 环境备忘

```bash
export PATH="$PWD/.tools/lean-4.19.0/bin:$PATH"     # Lean 4.19.0
/usr/bin/python3 scripts/build_paper.py              # 论文（两遍编译）
/usr/bin/python3 scripts/verify.py                   # 六门
lake env lean scripts/CheckAxioms.lean               # 57 条公理审计
```

- pyscf 2.15.0 只在 `/usr/bin/python3`，默认 `python` 没有
- 无 geomeTRIC；`reaction3.py` 自带 Hessian 播种的优化器
- `cc.kernel()` 返回**相关能**，总能是 `mycc.e_tot`
- pyscf 的 `frozen` 是**冻结轨道数**不是电子数
- Tectonic 在 `.tools/tectonic`

## 已知的未验证项

- UNVERIFIED：A-内积 Cauchy–Schwarz 的形式化（就是上面的待办）
- UNVERIFIED：CCSD(T) 的 CBS 外推没做，只用了两个基组
- UNVERIFIED：CCSD(T) 几何优化没做，只做了固定键长的点扫描
- UNVERIFIED：F–Ca–F 键长对势垒影响明显（def2-tzvpd 在 r=2.0256 时势垒塌到 0.005）
- UNVERIFIED：Schur 标量能否在 $O(k\cdot\mathrm{nnz})$ 内界定（Gauss–Radau 需要谱下界，
  又回到同一个要求，降了一层）
