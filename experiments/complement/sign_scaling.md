# 对角相似变换能否让 complement certificate 变为非空？

water/STO-3G，441 个 determinant（$S_z=0$），$d$ = diag 最小值 = $E_{\mathrm{RHF}}$。
本文件由 `sign_scaling.py` 直接从 `sign_scaling.json` 生成，所有数字无手工转抄。

## 结论

**不能。** 在任何对角相似变换 $D\,O\,D^{-1}$ 下，

$$\tau \;\ge\; \rho(D O D^{-1}) \;=\; \rho(O) \;=\; 2.026458\ \mathrm{Ha},$$

而 diagonal gap 只有 $0.447404$ Ha。最优变换是 $D=I$（$s=0$），数值优化从 6 个起点
出发得到的**最好**值就是 baseline 本身，改进量 $-8.9e-16$ Ha（机器零）。
Gershgorin 代理确实从 $3.7587$ 降到 $2.7660$（缩短 26%），
但它被钉在 $\rho(|O|) = 2.6270$ Ha 之上，仍比 gap 大 5.87 倍。

因此该证书在这个体系上：**sound but vacuous**，且这个 vacuity 不是数值精度问题，
而是 certificate 的假设本身是 basis-free 的谱条件（见 §5），换基无法改变其真假。

## 1. 设定与校验

Hamiltonian 由 pyscf 的 `contract_2e` 逐列构造，并与 `fci.FCI` 对照：

- $E_{\mathrm{RHF}} = -74.962706609$ Ha
- $E_{\mathrm{FCI}} = -75.012011634$ Ha（与 `fci.FCI` 相差 $-1.2\times10^{-12}$ Ha）
- $d = \min_i H_{ii} = -74.962706609$ Ha，$d = E_{\mathrm{RHF}}$（True）
- diagonal gap $= 0.447404347$ Ha
- $O = H - \mathrm{diag}(H)$：对称 = True，对角严格为 0
- $n = 441$ 个 determinant

## 2. 基线谱量

| 量 | 值 (Ha) |
|---|---|
| $\lVert O\rVert_2$ | 2.026457687 |
| $\rho(O)$ | 2.026457687 |
| $\lambda_{\max}(O)$ | 2.026457687 |
| $\lambda_{\min}(O)$ | -1.750272167 |
| $\lVert O\rVert_F$ | 11.534234817 |
| max row 2-norm | 1.090010848 |
| Gershgorin $\max_i\sum_j\lvert O_{ij}\rvert$ | 3.758707204 |

$O$ 对称 $\Rightarrow$ normal $\Rightarrow \|O\|_2 = \rho(O)$（数值上差 $-8.9e-16$）。
这一点很重要：**谱范数在此不是谱半径的上界，而是相等**，所以它是一个 similarity invariant。

注意证书只需要单侧不等式 $x^T O x \ge -\tau \|x\|^2$，因此真正 sharp 的常数是
$-\lambda_{\min}(O) = 1.750272$ Ha（$\tau/\mathrm{gap} = 3.91$），
而不是双侧的 $\|O\|_2 = 2.026458$（$\tau/\mathrm{gap} = 4.53$）。
`ComplementBound.lean` 里的假设是双侧的 $|Q x| \le \tau\|x\|^2$，故论文用的是后者。
**两种取法都远大于 gap，结论不变。**

## 3. 最小化 $\sigma_{\max}(D O D^{-1})$

### 3.1 定理（这就是答案）

对任意可逆对角 $D$：相似变换保谱，故 $\rho(DOD^{-1}) = \rho(O)$；又对任意矩阵
$\|\cdot\|_2 \ge \rho(\cdot)$，于是

$$\|D O D^{-1}\|_2 \;\ge\; \rho(O) \;=\; 2.026457687\ \mathrm{Ha},$$

且 $D=I$ 取等号。**下确界就是 baseline，无法改进。** 数值搜索若能找到低于
$\rho(O)$ 的值，那只能是 bug —— 脚本里就是这样断言的（`if best < nrm2 - 1e-9: abort`）。

### 3.2 解析梯度（为什么梯度法一开始就"收敛"）

令 $M(s) = \mathrm{diag}(e^s)\,O\,\mathrm{diag}(e^{-s})$，则
$\partial M/\partial s_i = E_{ii}M - M E_{ii}$。取 $Mv = \sigma u$、$M^T u = \sigma v$：

$$\frac{\partial \sigma_{\max}}{\partial s_i} = u^T(E_{ii}M - ME_{ii})v = \sigma\, v_i (u_i - v_i).$$

$M$ 正规（$u=\pm v$）时该梯度恒为 0。$s=0$ 时 $O$ 对称，所以 $\nabla_s\sigma_{\max} = 0$ ——
实测范数 $7.1e-16$。这不是优化器的缺陷，是定理穿着数值的外衣。

### 3.3 曲率：$s=0$ 是严格局部极小

沿 24 个随机的 sum-zero 方向做二阶差分：

| 统计 | 值 |
|---|---|
| min | +0.0114 |
| median | +0.0238 |
| max | +0.0613 |
| 全正？ | True |

曲率全正，说明 $s=0$ 附近单调下降**只能**朝 $s=0$ 走，不可能朝更小的 $\sigma$ 走。

### 3.4 多维 L-BFGS-B（如实报告）

解析梯度 + `scipy.optimize.minimize`，441 维，6 个起点：

| # | $\sigma_0$ | $\sigma_{\text{final}}$ | iters |
|---|---|---|---|
| 0 | 2.0461 | 2.026457731 | 400 |
| 1 | 2.2837 | 2.026459064 | 400 |
| 2 | 4.5696 | 3.323500571 | 9 |
| 3 | 26.9711 | 8.137508617 | 12 |
| 4 | 2075.9724 | 144.903023551 | 19 |
| 5 | 2.0265 | 2.026457687 | 0 |

| 汇总 | 值 |
|---|---|
| 所有起点中最好值 | 2.026457687 Ha |
| baseline | 2.026457687 Ha |
| 改进 | -8.88e-16 Ha |
| 停在 baseline 之上 >1e-4 的起点数 | 3 |

**如实说明：3 个起点没有收敛到 floor，停在 3.3, 8.1, 144.9 Ha。** 原因是
$\sigma_{\max}$ 在 $s$ 上非凸（正如题设所说），且高振幅处 landscape 尺度极差
（$\sigma_0 = 2076$ 时梯度范数 ~$10^2$），L-BFGS-B 靠 `ftol` 提前退出。这不是理论失败，
是优化器失败；而且**失败方向是停在 baseline 之上**，从未低于它。
收敛的那 3 个起点（$\sigma_0 = 2.05, 2.28, 2.03$）全部精确落在 baseline 上。

### 3.5 去掉尺度病态的定向扫描

为把上面的"优化器失败"与"定理"分开，令 $s = t\,u$，$\|u\|=1$，$\sum u = 0$（去掉 gauge），
对每个方向在 $t\in[-2,2]$ 上做一维极小化（16 个随机方向）：

| 统计 | 值 |
|---|---|
| $\max\|t^\star\|$ | 2.66e-07 |
| 各方向最小值与 baseline 的最大偏差 | 3.55e-15 Ha |
| 无一低于 baseline | True |

每个方向的最优点都在 $t=0$，值就是 baseline，到机器精度。这是比多维多起点更干净的证据。

## 4. Gershgorin 与符号缩放

### 4.1 合法上界层级与各自的地板

| 上界 | 值 (Ha) | 地板 | 说明 |
|---|---|---|---|
| $\max_i\sum_j\lvert O_{ij}\rvert$（raw） | 3.758707 | — | 行和，sign-blind |
| Osborne 平衡后 Gershgorin | 2.765962 | $\rho(\lvert O\rvert)=2.627031$ | 凸（log 尺度），真实改进 26% |
| $\rho(\lvert O\rvert) = \lambda_{\max}(\lvert O\rvert)$ | 2.627031 | — | 绝对值守恒的地板 |
| 谱范数（sharp） | 2.026458 | 自身 | similarity invariant |

平衡幅度：$\max p - \min p = 2.1319$（1118 次 L-BFGS 迭代）。

$\rho(|O|)$ 是合法的 $\tau$，因为
$|x^T O x| \le \sum_{ij}|x_i||O_{ij}||x_j| = |x|^T|O||x| \le \lambda_{\max}(|O|)\|x\|^2$。
它比 $\rho(O)$ 大 $+0.600574$ Ha —— **这个差就是 certificate 丢掉的符号抵消**。

任何绝对值型上界都无法低于 $\rho(|O|)$：$A_s = DAD^{-1}$ 与 $A=|O|$ 相似，Gershgorin 圆盘包含谱，
故行和 $\ge \rho(A_s) = \rho(A)$。所以连"完美平衡的 Gershgorin"也到不了 $\rho(O)$，更到不了 gap。

### 4.2 符号缩放：恒等无效

若 $D = \mathrm{diag}(\varepsilon)$，$\varepsilon\in\{\pm1\}^n$，则

$$|(DOD^{-1})_{ij}| = |\varepsilon_i \varepsilon_j O_{ij}| = |O_{ij}|,$$

**任何绝对值型界都完全不变。** 实测值 $3.758707204$ Ha，与 raw 只差
$+0.0e+00$ Ha（`sign_scaling_is_vacuous` = True）。
所以"用符号缩放改进 Gershgorin"在数学上是不可能的：Gershgorin 已经取了绝对值。

### 4.3 一个必须记录的死胡同（易犯错误）

直观上会想去平衡**符号**以让项相互抵消，即最小化
$\max_i\big|\sum_j \varepsilon_i\varepsilon_j |O_{ij}|\big|$。该目标确实能压到
$0.7476$ Ha —— **低于 $\rho(O) = 2.0265$**。

任何真正的谱上界都不可能低于谱半径，所以这不是上界。更直接的证伪：把它自称要界住的那个矩阵
建出来，$B = \mathrm{diag}(\varepsilon)\,|O|\,\mathrm{diag}(\varepsilon)$，它的
$\lambda_{\max} = 2.627031$ Ha，
比该"上界"$0.7476$ Ha 高
3.51 倍 —— 而 $|B|$ 的行和仍是
$3.758707$（因为 $|\varepsilon_i\varepsilon_j|=1$）。也就是说：符号抵消只发生在
**没有谱内容**的量里，Gershgorin 求和的是 $|S_{ij}|$，完全不受符号影响。
脚本把它保存在 `invalid_signed_sum` 里并标记 `is_a_valid_bound: false`。

## 5. 为什么任何对角变换都不可能奏效（证明，非观测）

证书的假设是 $\forall x:\; x^T O x \ge -\tau \|x\|^2$，等价于**谱条件**

$$O + \tau I \succeq 0.$$

- **相似变换** $DOD^{-1}$：不改变二次型 $x^T O x$，改的是另一个矩阵的二次型；它能影响的只有范数型上界，
  而该上界 $\ge \rho(O)$，保谱。实测：20 个随机 $D$ 下谱扰动
  $2.8e-14$（虚部 $0.0e+00$）。
- **合同变换** $D(O+\tau I)D$：这正是把 $x = D^{-1}y$ 代入二次型后唯一自然的变换。
  由 Sylvester 惯性定律，$O+\tau I \succeq 0 \iff D(O+\tau I)D \succeq 0$。
  实测：20 个随机 $D$ 下惯性改变 $0$ 次（从未改变）。

于是假设的真假是 **basis-free** 的。既然它是假的（$\lambda_{\min}(O) = -1.7503$，
且 $O$ 有 **87/441** 个特征值低于 $-\mathrm{gap}$），
任何换基都不把它变成真的。**这不是需要更聪明的 $D$，而是这个假设本身太弱。**

## 6. 松弛预算：问题出在哪

| 量 | 值 (Ha) | 相对 gap |
|---|---|---|
| diagonal gap | 0.447404 | 1.00 |
| 证书实际付出的 off-diagonal 风险 $\lVert O\rVert_2$ | 2.026458 | 4.53 |
| 同上（单侧 sharp） | 1.750272 | 3.91 |
| **真实需要的** $d - \lambda_{\min}(H)$ | 0.049305 | 0.11 |

证书为覆盖真实的 0.0493 Ha 支付了 1.7503 Ha（单侧）或 2.0265 Ha（双侧），
即 **35× / 41× 的过度保险**。
$\lambda_{\min}(O)$ 与 $\lambda_{\min}(H)$ 之间的全部差距就是能量的符号抵消 ——
而 $|x^T O x| \le \rho(O)\|x\|^2$ 这类界按定义把抵消全部丢掉。

## 7. 证书数值汇总

| $\tau$ 来源 | $\tau$ (Ha) | $d-\tau$ (Ha) | 与 $E_{\mathrm{FCI}}$ 之差 (Ha) | $\tau/\mathrm{gap}$ | $\tau < \mathrm{gap}$? |
|---|---|---|---|---|---|
| $\lVert O\rVert_2$ (two-sided) | 2.026458 | -76.989164 | -1.977153 | 4.53 | 否 |
| $-\lambda_{\min}(O)$ (one-sided, sharp) | 1.750272 | -76.712979 | -1.700967 | 3.91 | 否 |
| Gershgorin raw (row sums) | 3.758707 | -78.721414 | -3.709402 | 8.40 | 否 |
| Gershgorin balanced (Osborne) | 2.765962 | -77.728669 | -2.716657 | 6.18 | 否 |
| $\rho(\lvert O\rvert)$ absolute-value floor | 2.627031 | -77.589738 | -2.577726 | 5.87 | 否 |

没有任何一行是非空的。最小可行 $\tau$ 为 1.750272 Ha，是 gap 的 3.91 倍。

## 8. 复现

```
/usr/bin/python3 experiments/complement/sign_scaling.py   # ~40 s（H 已缓存）
```

输出 `sign_scaling.json`（全部中间量）、`sign_scaling.log`（完整运行日志）、
`sign_scaling.partial.json`（各 stage 后的检查点）。Hamiltonian 缓存在
`H_water_sto3g.npy`，仍会与 `fci.FCI` 重新对账。

内置的一致性断言（全部通过，且抓到过两个真实 bug）：
1. $\lVert O\rVert_2 \ge \rho(O)$；
2. 任何搜索得到的 $\sigma_{\max}$ 不得低于 $\rho(O)$；
3. 平衡后 Gershgorin 不得低于 $\rho(|O|)$；
4. 相似变换下谱不变（用 `eigvals`，**不能用** `eigvalsh` ——$DOD^{-1}$ 不对此称，
   这个错误一度给出 $1.6\times10^3$ 的假偏差）；
5. 符号和尝试必须低于 $\rho(O)$（否则该前提本身要重新检查）。

### 独立核验

`audit_sign_scaling.py` 不引用本文件任何函数，从缓存的 Hamiltonian 重建 $O$、用另一条
代码路径重算全部头条数字，并用**暴力搜索**（20000 个随机对角相似变换 + 2000 个随机
幅度缩放）直接检验中心论断，而不是相信 L-BFGS-B。结果写入
`sign_scaling_audit.json`。这一步存在的理由：本次工作中两个 bug 都是被
一致性断言抓到的，而不是被"看起来对"抓到的。

## 9. 局限

- 单分子、单基组、单 sector（441 determinants）。这是"真值可知"的模型测量，
  不推断任何 scaling law。
- 结论 $\tau \ge \rho(O)$ 是 exact 的、与体系无关的；但"$\tau > \mathrm{gap}$"
  是这一体系的经验事实，不是定理。
- 未测试非对角（非 similarity、非 congruence）的其它界，例如利用 excitations 分块的
  block Gershgorin、或对 $x^T O x$ 做 sign-aware 的 resolvent 估计。这类界的成立与否
  不违反 §5 的论证（因为 §5 只排除"换基 + 范数上界"这一族），但不在本次范围内，
  也没有做任何验证。**不要把本文件的负结果读作"该证书在原则上无法收紧"。**
