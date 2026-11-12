# SOS、拓扑与 BO 数值证据审计

审计基线：`89cb9248bb514cf23bb1503a42364c6a7fbf6b46`。下列 `文件:行号` 均相对此源码基线。已通读两个指定 SOS 文件、`ald_topology/*.py` 的四个文件、两个 BO 脚本，并交叉检查 `sos_n_representability.py`、保存结果及相关 Lean/LaTeX。仅新增本报告；没有运行 DFT、提交代码或修改已有文件。

**结论：目前能复现的是常量表格的算术和人工点云上的计算，不能据此声称 SOS 已求解并达到化学精度、拓扑已通过独立 DFT 验证，或 CaF₂ ALD 的 BO 误差已有严格上界。** “β₀ 等于最大非重叠位点数”有明确反例。

## 1. SOS：硬编码数值与未成立的数学推论

| 发现 | 源码证据与含义 |
|---|---|
| 能量未由求解器生成 | `sos_relaxation_core.py:197–235` 直接构造 HF/PQG/SOS/3-pos/FCI 数值；`sos_tight_construction.py:134–174` 同样返回常量，注释明确写有 “literature data + interpolation”。该处没有逐项文献定位、积分输入或求解日志。 |
| “构造约束”主要返回文字 | `sos_relaxation_core.py:59–162`、`sos_tight_construction.py:45–117` 返回字符串列表；没有据此组装、求解并验收 SDP。补查的 `sos_n_representability.py:73–77` 将标准正性留为 `pass`，`:174–180` 返回零矩阵，`:301–353` 仍使用硬编码能量。 |
| “SOS 验证”不是证书 | `sos_n_representability.py:224–248` 只在随机点测试正性；有限采样既不能证明全局非负，也不能证明存在 SOS Gram 分解。 |
| `<1%` 与化学精度是不同指标 | tight 表计算的是相对相关能、相对 **3-pos** 的 gap（`sos_tight_construction.py:23–24`），不是相对 FCI 的绝对能量误差。 |

实际调用两个 `numerical_validation*()`，得到：

| 常量表中的体系 | basic → 3-pos gap | tight → 3-pos gap | tight 与 FCI 的绝对差 / Hartree |
|---|---:|---:|---:|
| H₄ | 1.712963% | 0.555556% | 0.0026 |
| N₂ | 3.512015% | 0.739372% | 0.0060 |
| 苯 | 1.462194% | 0.258034% | 0.0124 |
| 算术平均 | 2.229057% | 0.517654% | — |

这些是**输入常量的复算**，不是独立量子化学结果。三项 tight–FCI 差均大于仓库采用的化学精度阈值 `0.0016 Hartree`（`caf2_ald_bo_analysis_simple.py:36–38`）。运行 core 主程序还同时输出 `Target achieved: False (< 1%)` 与 `✓ Gap to 3-positivity: 2.23% < 1% target`；前者按条件计算，后者是无条件模板（`sos_relaxation_core.py:494–496,513–519`）。

还有以下数学问题，不能靠重算百分比修复：

- **变分松弛能量方向不符。** 对同一 Hamiltonian、同一基组，若集合满足 \(\mathcal D_N\subseteq R_{\rm tight}\subseteq R_{\rm loose}\)，最小化必有
  \[
  E_{\rm loose}\le E_{\rm tight}\le E_{\rm FCI}.
  \]
  两张表却把各“松弛最优能量”放在 FCI 之上，且加约束后能量降低。除非这些只是尚未收敛的试探值或采用了另一问题定义，否则与所称外松弛最优值不相容。
- **维数与复杂度未建立。** `sos_relaxation_core.py:144–151` 将 \(\binom{K^4+3}{3}\) 误写成 \(O(K^6)\)；按其声明的 \(K^4\) 个交换变量，完整三阶单项式基规模应为 \(\Theta(K^{12})\)。对 \(K\) 个自旋轨道，3-RDM 的行索引是三元组，反对称基维数为 \(\binom K3\)，矩阵元素数才是 \(O(K^6)\)，不是 `:257` 所写的 \(K^6\times K^6\)。`:164–178` 的复杂度与加速比只是格式化公式，不能证明求解耗时；若采用稀疏或费米子算符基，须另行定义并分析。
- **2-RDM 不能普遍重建 3-RDM。** `sos_relaxation_core.py:259–269` 声称后者是前者的多项式。取固定三电子态
  \[
  |\psi_\pm\rangle=(|123\rangle\pm|456\rangle)/\sqrt2.
  \]
  任意二体算符不能连接相差三个占据轨道的这两个行列式，因此两态的全部 2-RDM 元素相同；但按 \({}^3D_{123,456}=\langle a_1^\dagger a_2^\dagger a_3^\dagger a_6a_5a_4\rangle\) 定义，值分别为 \(1/2,-1/2\)。已用标准库 `Fraction` 和费米子位串符号规则穷举核验全部 2-RDM 元素。
- **增强不等式会排除物理态。** 令 \(|\psi\rangle=(|12\rangle+|34\rangle)/\sqrt2\)，采用 \(D_{ij,kl}=\langle a_i^\dagger a_j^\dagger a_l a_k\rangle\)。有 \(n_i=1/2\)、\(D_{12,12}=1/2\)、\(D_{13,13}=D_{12,13}=0\)，且 1-RDM 非对角元为零。`sos_relaxation_core.py:89` 要求 \(D_{12,12}^2\le n_1n_2D_{12,12}\)，即 \(1/4\le1/8\)，为假；`sos_tight_construction.py:78–79` 取两对指标 \((1,2),(1,3)\)，得到 \(0\le(1/4)(-1/4)=-1/16\)，也为假。已以有理数精确核验。这些式子不能作为通用 N-可表示性的必要条件。
- **所谓证明存在逻辑缺口。** `sos_tight_construction.py:191` 要构造 \(D\in SOS_6\setminus R_3\)，`:215` 却得出 \(SOS_6\subsetneq R_3\)，方向自相矛盾。`:205–213` 没有给出把 Motzkin 多项式转为满足全部 RDM 条件的实例；`:253–267` 直接声称两个 `0.005` 上界及“90–95%”信息覆盖率，缺少推导。`sos_relaxation_core.py:336–437` 的 Lean 内容是带 `True`、`sorry` 等占位的返回字符串，不能提供核验过的证明。

## 2. 拓扑：“DFT reference”是预设值

四个文件均未通过 DFT 计算产生参考位点数：

| 文件与行号 | 实际参考来源 |
|---|---|
| `ald_topology/persistent_homology_ald.py:390–412` | `dft_sites = 3`；只有概略作者、期刊和年份说明。 |
| `ald_topology/ald_topology_predictor.py:362–381` | 直接返回 `1, 3, 4, 6` 字典；`:151–169` 的预测权重同样直接指定。 |
| `ald_topology/ald_site_predictor_v2.py:366–367` | `dft_sites = 3`；`:178–200` 手定权重，`:381–384` 在同一参考上选“最佳方法”。 |
| `ald_topology/ald_persistent_homology_final.py:221–275,297–312` | 手写坐标预设“两簇加一孤立点”，再令 `dft_sites = 3` 比较并打印 `PROVEN`。 |

最终版确实在这些坐标上算出 3 个连通分量；参考值也确实是 3。但这只能核验点云算法与预设结果一致，不能验证真实表面的吸附能力。`ald_topology/validation_final.json:8–22` 的 `0.0% / PROVEN` 来自上述过程。当前代码及 `references.bib` 没有给出可把该 **特定表面结构与三个位点** 对应到原始 DFT 输出的证据；不能据此断言相关文献不存在。

实测数值函数：初版 `persistent_homology_ald.py` 预测 20；`ald_topology_predictor.py` 预测 20；v2 的 steric/topology/combined 都为 10；final 为 3。后面三项与保存的 `results.json:14–24`、`validation_results.json:12–27`、`validation_final.json:14–22` 一致。不同版本改变了点云、规则和手定参数，这不是同一独立测试集上已证明的精度改进。

### 严格反例：连通性不等于不重叠容量

定义冲突图 \(G_r\)：候选位点为顶点，距离 \(\le r\) 时连边，表示两点不能同时占据。Rips 复形的 \(\beta_0\) 是该图的连通分量数；在**仅有这些两两排斥约束**的模型下，最大同时占据数是独立数 \(\alpha(G_r)\)。

取 \(p_0=(0,0,0),p_1=(3,0,0),p_2=(6,0,0)\)，单位 Å，\(r=3.5\) Å。边恰为 \(\{0,1\},\{1,2\}\)，形成三顶点路径。整图连通，故 \(\beta_0=1\)；两端距离为 6，可同时选取，三个点又不能全部选取，故 \(\alpha=2\)。严格 `<` 或 `≤` 阈值均得相同反例。

一般地，非空连通分量 \(C_i\) 满足
\[
\alpha(G)=\sum_i\alpha(C_i)\ge\#\{C_i\}=\beta_0(G).
\]
等号成立当且仅当每个连通分量都是完全图：完全图至多选一个顶点；非完全的分量至少有两个不相邻顶点。因此，最终版“每个连通簇只能容纳一个 TMA”的规则（`ald_persistent_homology_final.py:91–95,364–368`）需要额外假设。上述反例否定的是规则的一般性，并未宣称真实表面的容量完全由此图决定。

最小标准库复现，在仓库根目录执行：

```sh
python -B - <<'PY'
from itertools import combinations
from math import dist
p = [(0, 0, 0), (3, 0, 0), (6, 0, 0)]
edges = [(i, j) for i, j in combinations(range(3), 2)
         if dist(p[i], p[j]) <= 3.5]
parts = [{i} for i in range(3)]
for i, j in edges:
    a = next(s for s in parts if i in s)
    b = next(s for s in parts if j in s)
    if a is not b:
        a.update(b)
        parts.remove(b)
alpha = max(len(s) for n in range(4) for s in combinations(range(3), n)
            if all(not (i in s and j in s) for i, j in edges))
assert (len(parts), alpha) == (1, 2)
print(edges, len(parts), alpha)
PY
```

已实测输出：`[(0, 1), (1, 2)] 1 2`。调用 final 的 `compute_betti_at_scale` 也返回 1。另有两个可复现实现问题：

- v2 的 `apply_steric_constraint` 声称求最大值，却按输入下标贪心（`ald_site_predictor_v2.py:110–139`）。上述点按 `[左端,中点,右端]` 输入返回 2，按 `[中点,左端,右端]` 输入返回 1；它只能保证得到极大独立集，不能保证最大。
- final 的同一函数在上述点、`scale=1.0` 时返回 **0**，而应有 **3** 个分量。GUDHI 3.13.0 默认不输出最高维持久同调；无边时复形实际维数为 0，`:68–75` 的默认调用漏掉全部 H₀。对相同复形调用 `persistence(persistence_dim_max=True)` 实测得到三个 `[0,∞)` 区间。

此外，初版 `persistent_homology_ald.py:99–127`、中间版 `ald_topology_predictor.py:95–135`、v2 `ald_site_predictor_v2.py:91–108` 统计的是超过寿命阈值的条形数。一般的尺度 Betti 数应统计 \(b\le r<d\) 的存活区间；不能把所有维度的“长条形数”直接解释为同一尺度的 Betti 数。H₀ 在全部出生于零等特定条件下可以对应某阈值处的分量数。

## 3. BO：质量尺度可复算，反应误差界未证明

- **输入与输出的性质。** `caf2_ald_bo_analysis.py:28–34` 将 `epsilon**5` 命名为误差界；`:99` 指定 `gamma=5.0`；`:164–172` 指定能隙 `4.5,6.8,1.2 eV` 和速度 `0.1 Å/fs`；`:224–229` 再指定 `50 kcal/mol` 能标。简化版对应 `:25–34,80–85,99–105,132–136`。没有电子态、过渡态或非绝热耦合的实际计算输入。
- **大 O 不是数值上界。** 即便适用的渐近定理给出 \(|\delta E|\le C\varepsilon^5\)，仍须说明近似阶数、适用域、谱隙及正则性条件，并界定常数 \(C\) 和单位。直接取系数 1，再把它当相对误差乘任意反应能标，不能推出具体反应的误差保证。
- **耦合估计缺少电子态导数。** 对非简并本征态，微分本征方程得到
  \[
  \langle\phi_n|\nabla_R\phi_m\rangle
  =\frac{\langle\phi_n|\nabla_R H_e|\phi_m\rangle}{E_m-E_n},\quad n\ne m.
  \]
  代码 `caf2_ald_bo_analysis.py:139–148` 用速度除能隙替代导数耦合，没有该矩阵元及一致的量纲推导。它还将速度乘 `0.01`，而简化版 `:104` 直接用 `0.01/gap`。同取 `4.5 eV` 时，实测两者的所谓耦合参数分别为 `0.000144763085286` 和 `0.001447630852864`，相差十倍。
- **单位比较错误与无条件结论。** `caf2_ald_bo_analysis.py:241`、简化版 `:148` 比较 kcal/mol 数值与 Hartree 数值。简化版实际输出 `bo_error_kcal = 0.004433632`，应与 `1.0 kcal/mol` 同单位比较，却与 `0.0016` 比较，导致打印“approaches chemical accuracy”；`:160–166` 又无条件打印 BO 有效。修复这个比较也不能补足实际误差界。
- **无交叉的结论没有计算支撑。** `caf2_ald_bo_analysis.py:194–208,259–265` 和简化版 `:111–124,160–166` 仅打印闭壳层、假定能隙及结论。单个轨道能隙数值和闭壳层标签不能证明整条反应路径上的电子谱隙，更不能排除所有锥形交叉。

运行 `python -B caf2_ald_bo_analysis_simple.py` 确认：

```text
μ(H-F) = 0.9572120364 amu
ε(H-F) = 0.1547240893
ε(H-F)^5 = 8.867264742147110e-05
50 * ε(H-F)^5 = 0.004433632 kcal/mol  # 采用脚本假设能标
```

这里最后一项是模型算术，不能标注为测得或严格界定的 BO 误差。`CalcChemistry/BornOppenheimer.lean:175–183` 却声明 `ε_HF ∈ (0.00302,0.00304)`、五次方在 `(2.4e-13,2.6e-13)`，与其自身质量定义及实际脚本复算不相容，证明处为 `sorry`；`:193–199` 直接把能隙和无交叉设为公理。`:92–98` 甚至要求一个固定有限常数约束**任意两个实数能量**之差，该命题按当前量词就不成立。未运行 Lean 编译，不把源码中的 theorem 名称当作已验证证据。

## 4. 验证范围与可保留的结论

默认 `python` 为 3.14.9，缺少 NumPy/SciPy/Matplotlib，导入 GUDHI 报缺 NumPy；改用已有 `/usr/bin/python3` 3.10.12，实测 NumPy 2.2.6、SciPy 1.15.3、GUDHI 3.13.0。拓扑测试在内存加载源码，仅去掉未使用的 `matplotlib.pyplot` 导入，再调用数值函数；未调用会覆盖 JSON、报告或图片的主程序。标准库反例、精确费米子算符检查、SOS 常量复算及 BO 简化版均已实际执行。

可保留：常量表的 gap 算术、人工点云的拓扑描述、给定质量下的无量纲参数计算。应撤回或改标为未验证：SOS 求解精度及加速保证、独立 DFT 验证、一般位点容量等式、具体 CaF₂ ALD 的严格 BO 误差和无交叉结论，以及“首次”新颖性判断。`STATUS_REPORT.md:14–33`、`sections/rdm_theory.tex:19–30`、`sections/approximations.tex:3–13` 已传播这些结论，本次仅标记问题，没有改动它们。

建立后续证据至少需要：固定 Hamiltonian/基组并提交真实 SDP/FCI 求解日志与残差；为具体表面提供可追溯结构和独立 DFT 位点标签；区分拓扑描述与排斥图优化；为 BO 指明物理量、适用定理条件、常数和真实电子态数据。把手定值转写成 Lean 或 LaTeX 不会增加它们的物理可信度。
