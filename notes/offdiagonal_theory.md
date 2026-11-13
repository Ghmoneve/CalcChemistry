# 非对角耦合的可约化性：数学理论综述

**任务**：为"complement-space certificate 可以被做成 non-vacuous"这一命题提供理论依据。
**实测参数**（水分子 / STO-3G，441 determinants）：diagonal gap `d = 0.447404 Ha`，off-diagonal spectral norm `‖E‖ = 2.026458 Ha`，`β/d = 4.529`，一阶 certificate `d − τ` 为 vacuous。

**本文件的核心判定（先说结论）**：

| 路线 | 是否有定理支持 | 在本参数下能否救活 certificate |
|---|---|---|
| (1) Lieb–Robinson / quasi-adiabatic continuation | 是，强定理 | **不能**。全部定理的假设是 `β ≪ Δ`，此处 `β/Δ = 4.53` |
| (1') Bravyi–DiVincenzo–Loss 弱耦合算法 | 是，强定理 | **不能**。阈值常数小到 `2^{-17}` 量级，差约 6 个数量级 |
| (2) numerical range / field of values | 是，**反向**定理 | **不能，且可证明不能**。Hermitian 情形 `w(A) = ‖A‖`，恒等，无任何改进 |
| (3) Schrieffer–Wolff 严格误差界 | 是，certified（非渐近） | **不能**。收敛半径 `β/Δ < 1` 已被违反 |
| (4) 二阶泛函：Schur/Feshbach、Kato–Temple | 是，精确判据 | **能，但机制不同**——靠的是残差平方，不是 `‖E‖` 变小 |
| (5) 对角相似（gauge scaling）压低 `‖E‖` | 是，有可达下界 `ρ(|E|)` | **未知**，需实测；但被 `ρ(|E|)` 硬性限制 |
| (6) 酉变换压低 `‖E‖` | 是，**反向**定理 | **不能**。`‖UEU*‖ = ‖E‖`，酉不变 |

> **范围说明**：441×441 的 Hermitian 矩阵本可直接精确对角化（`O(n³) ≈ 10^8` flop）。因此本文件的"certificate"只在其**作为可扩展的、不依赖对角化的结构性论断**时才有意义；对单个 441 维体系，精确答案本身就是可计算的。

---

## 0. 文献更正（重要）

任务描述中的两条引用**在 Crossref 中不存在**，需要更正：

1. **"Bravyi–DiVincenzo–Loss, *Polynomial-time classical simulation of quantum ferromagnets*"** —— 该标题的论文作者是 **Bravyi & Gosset**（PRL **119**, 100503, 2017）。Bravyi–DiVincenzo–Loss 三人在 2008 年发表的、与本任务相关的论文标题是 *Polynomial-Time Algorithm for Simulation of Weakly Interacting Quantum Spin Systems*（CMP **284**, 481–507）。两者常被混引。本文件分别引用。
2. **"Hastings, quasi-adiabatic continuation, PRB 69, 104431 (2004)"** —— PRB **69**, 104431 的实际标题是 *Lieb-Schultz-Mattis in higher dimensions*。quasi-adiabatic continuation 的原始文献是 **Hastings & Wen, PRB 72, 045141 (2005)**；严格化版本是 **Bachmann–Michalakis–Nachtergaele–Sims, CMP 309, 835 (2011)**。

---

## 1. 近似 / 几乎对角化的严格理论

### 1.1 设定

设 `H = H₀ + εV`，`H₀` 在某个 product basis / determinant basis 下对角，`V` 是局域（local）或有限程相互作用。记 `P₀` 为 `H₀` 的低能谱投影，`Δ` 为 `H₀` 在 `P₀` 与其补空间之间的谱隙。"almost diagonalisation" 问的是：**是否存在一个（拟局域的）酉变换 `U`，使 `U H U*` 在同一个基下近似分块对角？** 以及**有效耦合能压到多小？**

### 1.2 Lieb–Robinson 界与指数聚类（定理）

**定理 1.1（Lieb–Robinson, 1972）**〔[LR72]〕
设 `H = Σ_X Φ(X)` 是图 `Γ` 上的有限程相互作用，`‖Φ‖_μ := sup_x Σ_{X∋x} ‖Φ(X)‖ e^{μ|X|} < ∞`。则对支撑在 `X`、`Y` 上的算子 `A, B`、任意 `a > 0` 与 `t ∈ R`：

```
‖[τ_t(A), B]‖  ≤  2‖A‖‖B‖ C_a · e^{ a( v|t| − d(X,Y) ) }
```

其中 `τ_t(A) = e^{itH} A e^{-itH}`，`v`（Lieb–Robinson velocity）与 `C_a` 只依赖 `‖Φ‖_μ` 与图的度数。**这是一个定理，不是估计**；但它给出的是"信息传播速度有限"，本身并不直接说 off-diagonal 矩阵元衰减。

**定理 1.2（Exponential Clustering, Nachtergaele–Sims 2006）**〔[NS06]〕
若 `H` 有唯一基态且基态之上有谱隙 `γ > 0`，则对支撑分离距离为 `d` 的 `A, B`：

```
|⟨A B⟩ − ⟨A⟩⟨B⟩|  ≤  C ‖A‖‖B‖ · e^{ −γ d / (2v) }   ( + 边界项 )
```

即 **谱隙 + Lieb–Robinson ⇒ 指数聚类**。这是定理；常数 `C, v` 非显式最优。Hastings–Koma〔[HK06]〕独立给出同类结果（并覆盖 1D 任意温度下的小 `β`）。

**对 determinant-basis off-diagonal 的直接推论（定理）**：若单粒子轨道指数局域、且体系有谱隙，则两个 determinant 之间的 Hamiltonian 矩阵元随它们"差异区域"的直径指数衰减。**注意**：这衰减的是**单个矩阵元**，不是 `‖E‖` 这个算子范数。`‖E‖` 由最大奇异值控制，可以被大量中等矩阵元"堆"起来——这正是水分子 `‖E‖ = 2.03 Ha` 的机制。**"矩阵元指数衰减"⇒"`‖E‖` 小"这一步是无效推理，本文件中标记为 heuristic，且在此处是错的。**

### 1.3 Hastings 的 quasi-adiabatic continuation / spectral flow（定理）

**定理 1.3（Hastings–Wen 2005；严格化：AMNS 2011, Theorem 3.4 & 5.5）**〔[HW05], [AMNS11]〕
设 `H(s) = H(0) + Φ(s)`，`s ∈ [0,1]`，`‖H′(s)‖` 一致有界；并设（AMNS Assumption 2.1）谱可分解为 `Σ(s) = Σ₁(s) ∪ Σ₂(s)`，两者距离 `≥ γ > 0` **一致于 `s`**。令 `P(s)` 为 `Σ₁(s)` 上的谱投影。则存在酉族 `U(s)`，`U(0)=1`，使

```
U(s) P(0) U(s)* = P(s) ,
```

且 `U(s)` 由**almost-local generator** 生成：

```
D(s) = i ∫_R dt  W_γ(t) · e^{itH(s)} Φ′(s) e^{−itH(s)} ,
```

`W_γ(t) = ∫_t^∞ w_γ`，`ŵ_γ` 支撑在 `[−γ, γ]`。关键结论：`D(s)` 的对易子 `[D(s), A_X]` 的范数随 `dist(X, Y)` **亚指数（almost exponentially）衰减**。AMNS 的表述是 `G(I)(γR/2v)`，`G` 为 subexponential function。

**推论（Local Perturbations Perturb Locally, AMNS Theorem 3.4）**：对任意 `R>0` 存在支撑在 `X_R` 内的酉 `V_R(s)`，使
`‖U(s) − V_R(s)‖ ≤ C G(I)(γR/2v)`，从而 `‖P(1) − V_R(1)P(0)V_R(1)*‖ ≤ 2C G(I)(γR/2v)`。

**这是本任务最重要的严格结果**：它说明"用算子代数变换把 off-diagonal 耦合赶出去"**确实有定理**，而且误差随距离衰减。但——

> **假设是 `γ > 0` 一致于 `s`，且 `‖Φ′‖` 有界。** 它**不**给出 `‖E‖` 的衰减率与 `‖E‖`、`Δ` 的定量关系。定理保证"能分块对角化"，但**不保证有效耦合小到 `< d`**。把"指数局域 ⇒ certificate non-vacuous"当作推论，是 **NOT KNOWN**（见 §6）。

### 1.4 Bravyi–DiVincenzo–Loss：弱耦合下的定量定理

**定理 1.4（BDL 2008, Theorem 1）**〔[BDL08]〕
`H(ε) = H₀ + εV`，`H₀` 为单粒子项之和、谱隙 `Δ`，`V` 为图（最大度 `d`）上的 2-local 项、最大范数 `J`。若

```
|ε|  ≤  2 ε₀ ,     ε₀ = 2^{−18} · Δ / (d J)
```

则最小本征值**重数为 1**，且它与第二本征值的间隙 **≥ Δ/2**。

这是一个真正的"certified gap"定理（承自 Osborne / Yarotsky 的路线）。同一文中 Theorem 2 给出微扰级数的收敛半径。算法本身基于 generalized Kirkwood–Thomas ansatz（与 coupled cluster 同源）。

**定理 1.5（BDL 2011, Lemma 4.1）**〔[BDL11]〕
设 `V = Σ_A V_A` 是 `k`-local 微扰，strength 为 `J`，且**每一项 `V_A` 都 block-diagonal（保持 `P₀`）**。若

```
2^{k+2} J  <  Δ ,
```

则 `P₀` 至少包含 `H₀ + V` 的一个基态。

**对本参数的判定**：取 `k=2`，阈值是 `J < Δ/16 = 0.02796 Ha`；实测 `β = 2.026 Ha`，**超出 72 倍**。BDL 2008 的阈值 `2^{−17} Δ/(dJ)` 更严，超出约 `6×10⁵` 倍。

### 1.5 小结：定理 vs 启发式

| 陈述 | 地位 |
|---|---|
| Lieb–Robinson 界 | **定理**（LR72），常数非最优 |
| 谱隙 + LR ⇒ 指数聚类 | **定理**（NS06, HK06） |
| 存在 almost-local unitary 使 `U P₀ U* = P₁`，generator 亚指数局域 | **定理**（HW05, AMNS11 Thm 3.4） |
| 有效相互作用在 renormalisation 后指数衰减 | **定理**（在 LR + gap 假设下；见 AMNS Thm 3.4/5.5, NSY22） |
| 上述衰减 ⇒ FCI determinant-basis 的 `‖E‖` 变小 | **启发式，且在此处为假**（`‖E‖` 是算子范数，不是矩阵元） |
| 谱隙稳定性在 `β/Δ = 4.5` 时仍成立 | **NOT KNOWN**；已知定理只覆盖 `β/Δ ≲ 10^{−5}` |

---

## 2. Numerical range / field of values 能否比谱范数更紧？

**这是本任务中唯一可以给出"可证明的否定答案"的问题。**

### 2.1 核心定理：Hermitian 情形恒等，无改进

**定理 2.1（Toeplitz–Hausdorff 1918/1919；normality）**〔[Toe18], [Hau19], [GR97]〕
`A ∈ M_n(C)`。若 `A` 正规（normal，含 Hermitian），则

```
W(A) := { ⟨x|A|x⟩ : ‖x‖ = 1 }  =  conv σ(A) .
```

**推论 2.2**：若 `A = A*`（Hermitian），则 `W(A) = [λ_min(A), λ_max(A)] ⊂ R`，于是

```
w(A) := max{ |z| : z ∈ W(A) }  =  max( |λ_min|, |λ_max| )  =  ‖A‖₂ .
```

**因此对任何 Hermitian off-diagonal 块 `E`，numerical radius 恒等于 spectral norm，改进量为零。** 这不是常数松，是**恒等**。

> 对水分子：`E = E*`，故 `w(E) = ‖E‖₂ = 2.026458 Ha`，一字不差。用 numerical range 做 certificate 不会得到任何好处。

数值验证（`numpy`，`n = 20, 60, 200` 随机 Hermitian 零对角矩阵，用 support function `w(A) = max_θ λ_max( Herm(e^{−iθ}A) )` 精确计算）：

```
n= 20  ‖B‖₂=5.601921122  max|λ|=5.601921122  w(B)=5.601921122  |diff|=8.9e-15
n= 60  ‖B‖₂=10.494584260 max|λ|=10.494584260 w(B)=10.494584260 |diff|=8.9e-15
n=200  ‖B‖₂=19.800230833 max|λ|=19.800230833 w(B)=19.800230833 |diff|=7.1e-15
```

### 2.2 非正规情形：最好只能赢 2 倍

**定理 2.3（Goldberg–Tadmor 1982；Kittaneh 2005）**〔[GT82], [Kit05], [GR97]〕
对任意 `A ∈ M_n(C)`：

```
(1/2)‖A‖₂  ≤  w(A)  ≤  ‖A‖₂ ,
```

且常数 `1/2` **不可改进**。

数值验证（精确 support-function 计算）：

```
nilpotent J = [[0,1],[0,0]] : ‖J‖=1.00  w(J)=0.50   -> 比值恰为 1/2
[[0,1],[1e-03,0]]           : ‖K‖=1.0000 w(K)=0.5005 -> 逼近 1/2
随机非正规 n=10             : w/‖A‖ = 0.8401
随机非正规 n=30             : w/‖A‖ = 0.7056
```

**推论**：即使 `E` 非正规（例如经过非酉相似变换之后），numerical radius 相对 spectral norm 的**最大**收益是因子 2。对 `β/d = 4.53` 的缺口，因子 2 远远不够。

### 2.3 numerical range 真正有用的地方（但用不上）

**定理 2.4（field-of-values resolvent bound）**〔[TE05]〕
若 `z ∉ W(A)`，则

```
‖(zI − A)^{−1}‖₂  ≤  1 / dist( z, W(A) ) .
```

证明：若 `‖(zI−A)^{−1}‖ = 1/δ`，取单位向量 `x` 使 `‖(zI−A)^{−1}x‖ = 1/δ`，令 `y = (zI−A)^{−1}x / ‖·‖`，则 `‖(zI−A)y‖ = δ`，故 `|z − ⟨y|Ay⟩| ≤ δ`，即 `dist(z, W(A)) ≤ δ`。∎

这是 numerical range 的**正确用途**：给**非正规**算子的**预解式**界。但对 Hermitian `H`，`W(H) = [λ_min, λ_max]`，该界退化为
`‖(zI−H)^{−1}‖ ≤ 1/dist(z, [λ_min,λ_max])`，与直接用谱的界**完全一致**——依然零收益。

### 2.4 结论

> **对 Hermitian off-diagonal 矩阵，numerical range / field-of-values 界不可能严格优于 spectral norm。这是定理，不是经验。** 任何声称"用 field of values 收紧 Hermitian off-diagonal 界"的方案在此处必然无效。
>
> 唯一例外：若把问题搬到**非 Hermitian** 的相似变换后的算子（例如 §5.2 的对角 scaling），则最多可赢 2 倍（Thm 2.3），且此时必须同时处理非正交基带来的度规问题。

### 2.5 但存在另一类"更紧的泛函"（不是 numerical range）

必须区分：**"更紧的泛函"存在，但它不是 numerical range。** 见 §4：
- **Schur / Feshbach 判据**：耦合以 `β²` 进入（二阶），不是 `β`；
- **Kato–Temple**：以**残差平方** `ε²` 进入，与 `‖E‖` 无关。

这两者都是**精确判据（iff）**或**定理**，而且在本参数下确实可以非空。

---

## 3. Schrieffer–Wolff：certified 误差界，还是渐近展开？

### 3.1 严格表述

SW 变换是一个酉 `U = e^S`（`S` anti-Hermitian），使 `U(H₀+εV)U*` **分块对角**（保持 `P₀`）。严格构造用 **direct rotation**：`U = (R_{P₀}R_P)^{1/2}`，其中 `R_P = 2P − 1`。

**引理 3.1（BDL 2011, Lemma 2.2 / 2.3）**〔[BDL11]〕
若 `‖P − P₀‖ < 1`，则 `R_{P₀}R_P` 无负实轴本征值，`U` 唯一确定且 `U P U* = P₀`。进一步存在唯一 anti-Hermitian `S`，`U = e^S`，`S` 关于 `P₀` 块对角外（`P₀SP₀ = 0`），且 **`‖S‖ < π/2`**。

### 3.2 误差界是定理，不是渐近

**定理 3.2（BDL 2011, Theorem 1）**〔[BDL11]〕
存在阈值 `ε_c > 0`（只依赖阶数 `n`、谱隙 `Δ`、相互作用图最大度），使得对**所有** `|ε| < ε_c`，第 `n` 阶有效 Hamiltonian `H_eff^{⟨n⟩}` 的基态能量与 `H₀ + εV` 的基态能量之差满足

```
δ_n  ≤  c_m N |ε|^{n+1}          ( 显式： δ_n = O(N|ε|^{n+1}) )
```

`c_m, ε_m` 只依赖图度、自旋相互作用强度、`Δ` 与截断阶 `n`。

**定理 3.3（BDL 2011, Theorem 3，local SW）**〔[BDL11]〕
用 local SW 构造的 `H_eff^{⟨n⟩}`，存在阈值

```
ε_c = Ω( Δ · n^{−2} )
```

使 `|ε| < ε_c` 时 `δ_n = O(N|ε|^{n+1})`。

**关键点（回答任务问题）**：**这是 certified 误差界，不是渐近展开。** BDL 明确写道，该界"even in the region `‖εV‖ ≫ Δ` where the perturbative series for `H_eff` are expected to be divergent"——即误差界在微扰级数本身发散的区域仍然有效。但注意：

1. 误差界 `O(N|ε|^{n+1})` 是**误差**的界，不是**有效耦合**的界；
2. `ε_c` 的**显式常数极小**（BDL 2008 的 `2^{−18}`，BDL 2011 的 `α, β` 未优化）；
3. 定理**不**保证 `‖H_eff − P₀H P₀‖ < d`。

### 3.3 收敛半径

微扰级数在 `|ε|` 的收敛半径 `ε_c ~ Δ/(N²)`（BDL 2011 §4.4 提到 `ε_c ∼ 1/N²`）。用无量纲比 `β/Δ` 表述：**级数只在 `β/Δ < 1` 时收敛**。实测 `β/Δ = 4.529 > 1`——**在收敛半径之外**。

因此 SW 在**任意阶**都不给出改进：第 `m` 阶残余是 `~ (β/Δ)^m · β`，随 `m` **发散**。

### 3.4 判定

| 陈述 | 地位 |
|---|---|
| SW 变换存在、唯一、`‖S‖ < π/2` | **定理**（BDL11 Lem 2.2/2.3） |
| 截断误差 `δ_n ≤ c N |ε|^{n+1}` | **定理**（BDL11 Thm 1/3）——certified，非渐近 |
| 在 `‖εV‖ ≫ Δ` 时误差界仍有效 | **定理**（BDL11 明言） |
| 在 `β/Δ = 4.53` 时存在某阶 SW 使有效耦合 `< d` | **NOT KNOWN，且与收敛半径矛盾** |
| 把 SW 用在 FCI determinant space（无格点结构） | **不适用**：BDL 的定理要求格点/图上的 `k`-local 结构，FCI 的 determinant 图不是这种结构 |

> **重要限制**：BDL 的 SW 定理建立在**格点自旋系统**（`H₀ = Σ_u H_{0,u}` 单格点项 + `k`-local `V`）之上。分子电子 Hamiltonian 在 determinant basis 下**不是**这种结构（Jordan–Wigner 后是非局域的，或说存在长程 string）。因此 BDL 定理**不能直接搬到**本题。这是一个**结构性不匹配**，不是常数问题。

---

## 4. 真正可行的替代：二阶泛函

**这是本文件对"如何让 certificate non-vacuous"的正面回答。**

### 4.1 Schur complement / Feshbach 判据（精确，iff）

**定理 4.1（Feshbach 1962；Löwdin 1962；Kato 1995）**〔[Fes62], [Low62], [Kat95]〕
设 `H` 自伴，`H = P𝔥 ⊕ Q𝔥`，`P + Q = 1`，分块为

```
H = [[ H_PP , H_PQ ],
     [ H_QP , H_QQ ]]
```

设 `E < λ_min(H_QQ)`（即 `H_QQ − E ≻ 0`）。则

```
E ≤ λ_min(H)     ⟺     H_PP − E − H_PQ (H_QQ − E)^{−1} H_QP  ⪰  0   (在 P𝔥 上)
```

**这是充要条件，不是近似。** 推导来自恒等式

```
H − E = [[I, H_PQ(H_QQ−E)^{−1}],[0, I]] · [[ H_PP−E−H_PQ(H_QQ−E)^{−1}H_QP , 0],[0, H_QQ−E]] · [[I,0],[(H_QQ−E)^{−1}H_QP, I]]
```

的同余变换（congruence），同余保持半正定性。

**推论 4.2（可计算的一阶充分条件）**
设 `γ := λ_min(H_QQ)`，`β := ‖H_PQ‖`。则

```
λ_min(H)  ≥  d − β² / (γ − d)          ( d := λ_min(H_PP) , 要求 γ > d )
```

**耦合以 `β²` 进入**——这是与 `d − β` 的本质区别。

**数值验证（`numpy`，随机 Hermitian，`p=2, q=15`，`H_QQ ≥ 3`）**：Schur 判据的最小本征值在 `E = λ_min(H)` 处**精确变号**：

```
E = λ_min         min eig Schur = -3.8e-17
E = λ_min + 1e-6  min eig Schur = -1.25e-06
E = λ_min − 1e-6  min eig Schur = +1.25e-06
```

（3 组随机试验全部如此。）⇒ 判据确实是 **iff**。

### 4.2 Schur vs 一阶（Weyl）的交叉点

模型（单参考态 + `q=60` 维补空间，`d = 0.447404`）：

```
β=0.10 | true λ₀=-0.0079 | Weyl=-0.1000 | Schur=-0.0224 | Schur 更紧
β=0.30 | true λ₀=-0.0711 | Weyl=-0.3000 | Schur=-0.2012 | Schur 更紧
β=0.44 | true λ₀=-0.1290 | Weyl=-0.4400 | Schur=-0.4327 | Schur 更紧
β=0.60 | true λ₀=-0.2298 | Weyl=-0.6000 | Schur=-0.8046 | Schur 更松
β=2.03 | true λ₀=-1.5064 | Weyl=-2.0300 | Schur=-9.2107 | Schur 更松 4.5 倍
```

**交叉点恰在 `β ≈ d`**，与推论 4.2 一致。**在实测 `β = 2.026` 处，Schur 界比一阶界更差 4.5 倍。** 所以 Schur 路线**不能**救活本参数下的 certificate。

### 4.3 Kato–Temple 界（残差平方，与 `‖E‖` 无关）

**定理 4.3（Temple 1928；Kato 1949）**〔[Tem28], [Kat49], [Kat95]〕
`H` 自伴，`ψ` 归一化，`ρ := ⟨ψ|H|ψ⟩`，`ε := ‖(H − ρ)ψ‖`。设存在 `γ` 满足 `ρ < γ ≤ λ₁(H)`。则

```
λ₀(H)  ≥  ρ − ε² / (γ − ρ) .
```

**这是二阶的**：误差是 `ε²`，不是 `ε`（对比 Weinstein 界 `λ₀ ≥ ρ − ε`）。

**自洽选取 `γ`（关键实操点）**：取 `Q = 1 − |ψ⟩⟨ψ|`，令 `γ := λ_min(QHQ)`。由 Cauchy interlacing，`γ ≤ λ₁(H)`；且 `γ ≥ λ₀(H)`。因此 `γ` 是合法且**完全可计算**的（一个 `n×n` 对称本征问题）。

**数值验证**（随机 Hermitian `n=300`，`λ₀ = −23.774413077`，`λ₁ = −23.413098952`，gap `= 0.361314`）：

```
mix     ρ              ε          γ            KT 下界        true λ₀        certified error
1e-3   -23.77441272   3.61e-04   -23.413099   -23.77441308   -23.77441308   3.7e-13 Ha
1e-2   -23.77437695   3.61e-03   -23.413135   -23.77441308   -23.77441308   3.6e-09 Ha
3e-2   -23.77408819   1.08e-02   -23.413424   -23.77441337   -23.77441308   2.9e-07 Ha
1e-1   -23.77083571   3.58e-02   -23.416676   -23.77444921   -23.77441308   3.6e-05 Ha
3e-1   -23.74457980   9.94e-02   -23.442932   -23.77736362   -23.77441308   3.0e-03 Ha
5e-1   -23.70215025   1.45e-01   -23.485362   -23.79850069   -23.77441308   2.4e-02 Ha
```

全部满足 `KT 下界 ≤ λ₀`（**无 violation**），且在残差 `3.6×10⁻²` 时仍给出 `3.6×10⁻⁵ Ha` 的 certified 误差。

**结论**：**Kato–Temple 是唯一在本参数下可非空的一类严格界**，因为它的收敛量是**试探向量的残差**，而非 `‖E‖`。对一个 CISD 质量的 `ψ`，残差通常在 `10⁻²–10⁻¹ Ha`，对应 certified 误差 `10⁻⁵–10⁻³ Ha`。

### 4.4 代价

Kato–Temple 需要：
1. 一个归一化试探向量 `ψ`（CISD / selected CI 即可）；
2. `ρ = ⟨ψ|H|ψ⟩`、`ε = ‖Hψ − ρψ‖`（两次 matvec）；
3. `γ = λ_min(QHQ)`——**一个 `n×n` 对称本征问题**（`n = 441` 时约 `10⁷` flop）。

即：**不需要任何 off-diagonal 范数界**。这正面回答了"certificate 能否 non-vacuous"：**能，但必须换成残差型泛函，而不是范数型。**

---

## 5. 算子代数变换能做什么、不能做什么

### 5.1 酉变换：`‖E‖` 是酉不变量（定理，否定）

**定理 5.1**：对任意酉 `U`，`‖U E U*‖₂ = ‖E‖₂`（奇异值不变）。因此**任何**酉基变换都不改变 off-diagonal 块的 spectral norm。

**但**：这不等于"基选择无用"。基选择改变的是**分裂** `H = D + E` 本身（`D` 依赖基），不是共轭 `UEU*`。所以：
- 在**固定** determinant 基下，不存在酉变换能压低 `‖E‖`——**定理**；
- 改变**单粒子轨道**（⇒ 改变整个 determinant 基）会改变 `‖E‖`——**这是真的，但无定量定理**（见 §6）。

### 5.2 对角相似（gauge scaling）：有精确的可达下界

设 `D` 为可逆对角矩阵（正对角）。`D H D^{−1}` 与 `H` **谱相同**，但一般不再 Hermitian——它是关于内积 `⟨x,y⟩_D := ⟨Dx, Dy⟩` 的自伴算子。这是一个"算子代数式"的变换（在 `D`-度规下保持谱）。

**定理 5.2（Bauer 1963；Osborne 1960；Perron–Frobenius / Collatz–Wielandt）**〔[Bau63], [Osb60], [BP79]〕
设 `M ≥ 0` 不可约。则

```
min_{D > 0 diagonal}  ‖D M D^{−1}‖_∞  =  ρ(M) ,
```

且最优 `D` 为 `D = diag(1/x)`，`x` 为 `M` 的 Perron 向量。

**数值验证**（`numpy`，3 组随机非负不可约矩阵，`n=10`）：

```
ρ(|M|)=8.721563006   ‖V M V^{-1}‖_∞ = 8.721563006   rel.diff=1.6e-15
ρ(|M|)=8.162715679   ‖V M V^{-1}‖_∞ = 8.162715679   rel.diff=6.5e-16
ρ(|M|)=8.319883752   ‖V M V^{-1}‖_∞ = 8.319883752   rel.diff=1.1e-15
```

**推论 5.3（对 off-diagonal 耦合的硬下界）**：对任意正对角 `D` 与任意 absolute norm（含 `‖·‖₂` 与 `‖·‖_∞`），

```
‖D E D^{−1}‖  ≥  ρ(|E|) .
```

证明：`|D E D^{−1}| = D |E| D^{−1}`（entrywise，`D` 正对角），故 `‖D E D^{−1}‖ = ‖D|E|D^{−1}‖ ≥ ρ(D|E|D^{−1}) = ρ(|E|)`。∎

即**对角 scaling 能压到的最好程度不低于 `ρ(|E|)`**（entrywise 绝对值矩阵的 Perron 根）。这是一个 `O(n³)` 可计算的下界：如果 `ρ(|E|) ≥ d`，则任何对角 gauge 都**不可能**让 certificate 非空——路线被证明堵死。

**注意 `∞`-范数与 `2`-范数的差别**：`min_D ‖D M D^{−1}‖_∞ = ρ(M)` 是**定理**（Thm 5.2，可达）；`min_D ‖D M D^{−1}‖₂` 的**闭式未知**（§6 第 3 条）。因此推论 5.3 的 `2`-范数版本只给出**下界**，不给出可达性。

> **对本项目**：`ρ(|E|)` 的实测值**尚未计算**（需要实际的 441×441 `E`）。这是一个**可立即执行、成本 `O(n³)`、结论二分**的判定实验。见 §7。

**已知的缺口**：`min_D ‖D H D^{−1}‖₂`（2-范数版本）**没有闭式**——这是 **NOT KNOWN**（见 §6）。Thm 5.2 只覆盖 `∞`-范数。

### 5.3 stoquastic / 符号问题：一般情形下被证明堵死

若能把 `H` 经**对角酉**（gauge）变成所有 off-diagonal ≤ 0（**stoquastic**），则 Perron–Frobenius 给出基态为正，且 Collatz–Wielandt 商给出远比范数界紧的严格下界。

**定理 5.4（Troyer–Wiese 2005）**〔[TW05]〕
判定一个费米子 Hamiltonian 是否 sign-problem-free（等价于：是否存在使 off-diagonal 全同号的基）**是 NP-hard 的**。因此不存在一般的 gauge 变换能统一符号。

**对分子**：determinant basis 下 fermionic antisymmetry 给出交错符号，signed graph 一般**不平衡（frustrated）**。故 §5.3 路线在一般情况下**被定理堵死**——不是"还没找到"，是"一般不存在"。

（具体到某个分子，问题是可判定的，但没有先验保证。）

---

## 6. 明确 NOT KNOWN

以下各项**没有**已知定理，本文件不把它们当作依据：

1. **"指数局域 + 谱隙 ⇒ FCI determinant-basis 的 `‖E‖` 小于 diagonal gap"** —— **NOT KNOWN**，且如 §1.2 所述，这一步推理本身是无效的（矩阵元衰减 ⇏ 算子范数小）。
2. **"存在单粒子轨道基使 `β/d < 1`"** —— **NOT KNOWN**。Maximally localised Wannier functions〔[MV97]〕最小化的是 spread functional，**没有**任何定理把 spread 与 determinant-basis `‖E‖` 联系起来。这是纯粹的启发式。
3. **`min_D ‖D H D^{−1}‖₂` 的闭式** —— **NOT KNOWN**（`∞`-范数版本是定理，见 Thm 5.2）。
4. **"`β/Δ = 4.53` 是常数松，还是真的在收敛半径外"** —— 对 SW 的**渐近级数**，`β/Δ > 1` 确实在收敛半径外（定理级结论）。但 BDL 的**误差界**在 `‖εV‖ ≫ Δ` 时仍有效（定理），所以"BDL 定理在此处失效"仅指其**阈值条件**不满足，不是定理被推翻。二者必须区分。**是否存在更好的（未发现的）certified 界覆盖此区域，NOT KNOWN。**
5. **把 BDL/AMNS 的格点定理搬到 FCI determinant space** —— **不成立**：FCI determinant 图不是 `k`-local 格点结构（Jordan–Wigner string 非局域）。**NOT KNOWN** 是否有对应的、针对 determinant graph 的定理。
6. **`ρ(|E|)` 对水/STO-3G 的实际值** —— **未计算**（本文档不伪造）。这是 §7 的第一个实验。
7. **`E` 的 row norm `‖P H Q‖`（单参考态）的实际值** —— **未提供**（给定数据只有 `‖E‖` 的谱范数）。这是 §7 的第二个实验，也是判定 Schur 路线是否可行的关键量。

---

## 7. 可执行判据（给主 agent 的具体建议）

按成本从低到高，三个**二分判定**实验。全部只需已有的 441×441 `H`。

### 实验 A：`ρ(|E|)` —— 对角 gauge 的理论极限

```python
import numpy as np
E = H - np.diag(np.diag(H))          # 441x441, Hermitian, zero diagonal
rho_abs = max(abs(np.linalg.eigvals(np.abs(E))))   # Perron root of |E|
```
- 若 `rho_abs ≥ d = 0.447404` → **对角 scaling 路线被 Thm 5.2/5.3 证明堵死**，放弃。
- 若 `rho_abs < d` → 该路线**在原理上可行**，需进一步构造 `D`（`D = diag(1/x)`，`x` 为 `|E|` 的 Perron 向量）并检验 `‖D E D^{−1}‖`。

### 实验 B：`β_row = ‖P H Q‖` —— Schur 路线是否可行

```python
P = np.zeros((n,1)); P[ref_index,0] = 1.0     # 单参考 determinant
Q = np.eye(n) - P @ P.T
HPQ = P.T @ H @ Q                              # 1 x 440
beta_row = np.linalg.norm(HPQ, 2)              # = row norm
```
判定式（推论 4.2）：Schur 界优于一阶界 ⟺ `beta_row < d`。
- 若 `beta_row < 0.447404` → **Schur certificate 非空**，直接用 `λ₀ ≥ d − beta_row²/(γ−d)`。
- 若 `beta_row ≥ 0.447404` → Schur 也不行，转实验 C。

> 注：`‖E‖ = 2.026` 是**整个** off-diagonal 的谱范数（440×440 块的奇异值），与 `β_row`（单行范数）是**不同的量**，后者通常小得多（`β_row ~ ‖E‖/√n` 量级，`n=441` 时约 `0.1 Ha`）。**这是一个有实质希望的缺口，但数值未知。**

### 实验 C：Kato–Temple —— 唯一不依赖 `‖E‖` 的路线

```python
psi = cisd_vector(); psi /= np.linalg.norm(psi)
rho = psi @ H @ psi
eps = np.linalg.norm(H @ psi - rho*psi)
Q   = np.eye(n) - np.outer(psi, psi)
gamma = np.linalg.eigvalsh(Q @ H @ Q)[0]       # <= lambda_1 by interlacing
lower = rho - eps**2/(gamma - rho)             # 定理 4.3，certified
```
- 已在随机体系上验证：残差 `3.6e-2` → certified 误差 `3.6e-5 Ha`（§4.3）。
- **此路线与 `β/d = 4.53` 完全无关**，是本参数下唯一有理论保证且数值上有希望的方案。

---

## 8. 参考文献（全部 DOI 已通过 `https://api.crossref.org/works/<DOI>` 验证解析，33/33）

**[LR72]** E. H. Lieb, D. W. Robinson, "The finite group velocity of quantum spin systems", *Commun. Math. Phys.* **28**(3), 251–257 (1972). DOI: 10.1007/BF01645779

**[NS06]** B. Nachtergaele, R. Sims, "Lieb-Robinson Bounds and the Exponential Clustering Theorem", *Commun. Math. Phys.* **265**(1), 119–130 (2006). DOI: 10.1007/s00220-006-1556-1

**[HK06]** M. B. Hastings, T. Koma, "Spectral Gap and Exponential Decay of Correlations", *Commun. Math. Phys.* **265**(3), 781–804 (2006). DOI: 10.1007/s00220-006-0030-4

**[BHV06]** S. Bravyi, M. B. Hastings, F. Verstraete, "Lieb-Robinson Bounds and the Generation of Correlations and Topological Quantum Order", *Phys. Rev. Lett.* **97**(5), 050401 (2006). DOI: 10.1103/PhysRevLett.97.050401

**[NS10]** B. Nachtergaele, R. Sims, "Lieb-Robinson bounds in quantum many-body physics", *Contemp. Math.* **529**, 141–176 (2010). DOI: 10.1090/conm/529/10429

**[HW05]** M. B. Hastings, X.-G. Wen, "Quasiadiabatic continuation of quantum states: The stability of topological ground-state degeneracy and emergent gauge invariance", *Phys. Rev. B* **72**(4), 045141 (2005). DOI: 10.1103/PhysRevB.72.045141

**[Has07a]** M. B. Hastings, "Quasi-adiabatic continuation in gapped spin and fermion systems: Goldstone's theorem and flux periodicity", *J. Stat. Mech.* **2007**(05), P05010 (2007). DOI: 10.1088/1742-5468/2007/05/P05010

**[Has07b]** M. B. Hastings, "An area law for one-dimensional quantum systems", *J. Stat. Mech.* **2007**(08), P08024 (2007). DOI: 10.1088/1742-5468/2007/08/P08024

**[AMNS11]** S. Bachmann, S. Michalakis, B. Nachtergaele, R. Sims, "Automorphic Equivalence within Gapped Phases of Quantum Lattice Systems", *Commun. Math. Phys.* **309**, 835–871 (2011). DOI: 10.1007/s00220-011-1380-0

**[NSY21]** B. Nachtergaele, R. Sims, A. Young, "Quasi-Locality Bounds for Quantum Lattice Systems. Part II. Perturbations of Frustration-Free Spin Models with Gapped Ground States", *Ann. Henri Poincaré* **23**(2), 393–511 (2021). DOI: 10.1007/s00023-021-01086-5

**[BHM10]** S. Bravyi, M. B. Hastings, S. Michalakis, "Topological quantum order: Stability under local perturbations", *J. Math. Phys.* **51**(9), 093512 (2010). DOI: 10.1063/1.3490195

**[BTW26]** L. Becker, S. Teufel, M. Wesle, "Automorphic equivalence within gapped phases of infinitely extended fermion systems", *J. Funct. Anal.* **290**, 111421 (2026). DOI: 10.1016/j.jfa.2026.111421

**[BDL08]** S. Bravyi, D. DiVincenzo, D. Loss, "Polynomial-Time Algorithm for Simulation of Weakly Interacting Quantum Spin Systems", *Commun. Math. Phys.* **284**(2), 481–507 (2008). DOI: 10.1007/s00220-008-0574-6

**[BDL11]** S. Bravyi, D. DiVincenzo, D. Loss, "Schrieffer–Wolff transformation for quantum many-body systems", *Ann. Phys.* **326**(10), 2793–2826 (2011). DOI: 10.1016/j.aop.2011.06.004

**[BG17]** S. Bravyi, D. Gosset, "Polynomial-Time Classical Simulation of Quantum Ferromagnets", *Phys. Rev. Lett.* **119**(10), 100503 (2017). DOI: 10.1103/PhysRevLett.119.100503

**[Weg94]** F. Wegner, "Flow-equations for Hamiltonians", *Ann. Phys. (Leipzig)* **506**(2), 77–91 (1994). DOI: 10.1002/andp.19945060203

**[Fes62]** H. Feshbach, "A unified theory of nuclear reactions. II", *Ann. Phys.* **19**(2), 287–313 (1962). DOI: 10.1016/0003-4916(62)90221-X

**[Low62]** P.-O. Löwdin, "Studies in Perturbation Theory. IV. Solution of Eigenvalue Problem by Projection Operator Formalism", *J. Math. Phys.* **3**(5), 969–982 (1962). DOI: 10.1063/1.1724312

**[Kat95]** T. Kato, *Perturbation Theory for Linear Operators*, Classics in Mathematics, Springer (1995). DOI: 10.1007/978-3-642-66282-9

**[Tem28]** G. Temple, "The theory of Rayleigh's principle as applied to continuous systems", *Proc. R. Soc. Lond. A* **119**, 276–293 (1928). DOI: 10.1098/rspa.1928.0098

**[Kat49]** T. Kato, "On the Upper and Lower Bounds of Eigenvalues", *J. Phys. Soc. Jpn.* **4**, 334–339 (1949). DOI: 10.1143/JPSJ.4.334

**[Toe18]** O. Toeplitz, "Das algebraische Analogon zu einem Satze von Fejér", *Math. Z.* **2**, 187–197 (1918). DOI: 10.1007/BF01212904

**[Hau19]** F. Hausdorff, "Der Wertvorrat einer Bilinearform", *Math. Z.* **3**, 314–316 (1919). DOI: 10.1007/BF01292610

**[GR97]** K. E. Gustafson, D. K. M. Rao, *Numerical Range: The Field of Values of Linear Operators and Matrices*, Universitext, Springer (1997). DOI: 10.1007/978-1-4613-8498-4

**[GT82]** M. Goldberg, E. Tadmor, "On the numerical radius and its applications", *Linear Algebra Appl.* **42**, 263–284 (1982). DOI: 10.1016/0024-3795(82)90155-0

**[Kit05]** F. Kittaneh, "Numerical radius inequalities for Hilbert space operators", *Studia Math.* **168**, 73–80 (2005). DOI: 10.4064/sm168-1-5

**[BF60]** F. L. Bauer, C. T. Fike, "Norms and exclusion theorems", *Numer. Math.* **2**, 137–141 (1960). DOI: 10.1007/BF01386217

**[Bau63]** F. L. Bauer, "Optimally scaled matrices", *Numer. Math.* **5**, 73–87 (1963). DOI: 10.1007/BF01385880

**[Osb60]** E. E. Osborne, "On Pre-Conditioning of Matrices", *J. ACM* **7**, 338–345 (1960). DOI: 10.1145/321043.321048

**[BP79]** A. Berman, R. J. Plemmons, *Nonnegative Matrices in the Mathematical Sciences*, Academic Press (1979), Ch. 2 "Nonnegative Matrices", pp. 26–62. DOI: 10.1016/B978-0-12-092250-5.50009-6

**[TE05]** L. N. Trefethen, M. Embree, *Spectra and Pseudospectra: The Behavior of Nonnormal Matrices and Operators*, Princeton University Press (2005). DOI: 10.1515/9780691213101

**[MV97]** N. Marzari, D. Vanderbilt, "Maximally localized generalized Wannier functions for composite energy bands", *Phys. Rev. B* **56**, 12847–12865 (1997). DOI: 10.1103/PhysRevB.56.12847

**[TW05]** M. Troyer, U.-J. Wiese, "Computational Complexity and Fundamental Limitations to Fermionic Quantum Monte Carlo Simulations", *Phys. Rev. Lett.* **94**(17), 170201 (2005). DOI: 10.1103/PhysRevLett.94.170201

---

## 附录 A：数值验证清单

全部在本机 `numpy 2.5.4` 上运行（脚本：`~/.hermes/cache/scratch/offdiag/checks.py`, `kt.py`, `kt2.py`）。

| 断言 | 验证方式 | 结果 |
|---|---|---|
| Hermitian `w(A) = ‖A‖₂` | `n=20/60/200` 随机零对角 Hermitian，support-function 精确算 `w` | 差 `≤ 8.9e-15`，恒等 |
| `(1/2)‖A‖ ≤ w(A) ≤ ‖A‖`，`1/2` 最佳 | nilpotent `J`；`[[0,1],[ε,0]]` 族 | `w(J)/‖J‖ = 0.500`；`ε=1e-3` 时 `0.5005` |
| `min_D ‖DMD^{−1}‖_∞ = ρ(M)` | 3 组随机非负不可约，`D = diag(1/x)` | 相对差 `≤ 1.6e-15` |
| Feshbach/Schur 判据是 iff | 3 组随机 `2+15` 分块，扫描 `E` | 在 `λ_min` 处精确变号（`±1e-6`） |
| Schur vs Weyl 交叉点 = `β ≈ d` | `d=0.447404`，`β ∈ {0.10…2.03}` | `β ≤ 0.44` 时 Schur 更紧；`β ≥ 0.60` 时更松 |
| Kato–Temple 自洽 `γ = λ_min(QHQ)` | `n=300` 随机 Hermitian，6 个混合比 | 全部 `≤ λ₀`，残差 `3.6e-2` → 误差 `3.6e-5` |

## 附录 B：符号表

| 符号 | 含义 |
|---|---|
| `d` | diagonal gap，`λ_min(H_PP)` 或对角部分的基态-第一激发间隔（实测 `0.447404 Ha`） |
| `β` | `‖H_PQ‖`（off-diagonal 耦合的算子范数；实测谱范数 `2.026458 Ha`） |
| `β_row` | `‖P H Q‖`，单参考态的行范数（**未知**，见 §7 实验 B） |
| `γ` | `λ_min(H_QQ)`，补空间最低本征值（**未知**） |
| `Δ` | `H₀` 在 `P₀` 与补空间之间的谱隙 |
| `E` | off-diagonal 块 `H − diag(H)`（零对角，Hermitian） |
| `w(A)` | numerical radius，`max{|z| : z ∈ W(A)}` |
| `W(A)` | numerical range / field of values，`{⟨x|A|x⟩ : ‖x‖=1}` |
| `ρ(|E|)` | `E` 的 entrywise 绝对值矩阵的 Perron 根（对角 scaling 的理论极限） |
