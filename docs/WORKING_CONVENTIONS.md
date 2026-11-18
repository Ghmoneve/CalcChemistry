# 工作纪律（每次会话开始先读这份）

写于 2026-11-18。因多次工具调用格式损坏导致流程中断，把规则固化在这里。

---

## 一、工具调用：唯一最重要的是这条

**每次只发一个干净的调用，参数正确闭合，后面不留任何多余标签。**

失败模式（我反复犯的）：
- 写了 `</invoke>`、`</calls>` 这类不存在的闭合标签
- `</parameter>` 写了两遍
- 把参数内容写在思考块里，然后下面跟一个空的 invoke
- 一个消息里塞多个调用时括号嵌套错乱

正确做法：

```
<invoke>
  <parameter name="command">
cd /mnt/e/CalcChemistry && ls
  </parameter>
</invoke>
```

要点：
1. 一个 invoke 内有且仅有它需要的 parameter
2. 闭合标签数必须与开标签数一致
3. 调用之后**不要**再写任何标签
4. 不确定时，宁可分多次单调用，也不要在一条消息里堆多个

**自查动作**：发出前在心里数一遍开闭标签对。

---

## 二、硬约束：这些文件位置不能动

移动文件前必须确认没破坏下列依赖。

| 必须留在仓库根目录 | 原因 |
|---|---|
| `lakefile.lean` | `lean_lib «CalcChemistry»` 相对于它定位库 |
| `lean-toolchain` | elan 定位工具链 |
| `lake-manifest.json` | lake 依赖锁 |
| `CalcChemistry/` | Lean 库源码目录，名字固定 |
| `CalcChemistry.lean` | 库根导入文件 |

改动论文位置后必须同步修改：
- `scripts/build_paper.py`：`ROOT / "paper_main.tex"`、`CalcChemistry_paper_standalone.tex`、`references.bib`、`build/` 输出
- `scripts/check_manuscript.py`：`ROOT / 'paper_main.tex'`、`references.bib`、`\input{}` 子文件的相对解析
- `scripts/verify.py`：调用上面两个脚本，路径为 `scripts/xxx.py`

**每次移动后立刻跑 `scripts/verify.py`，六个门必须全 exit 0。**

---

## 三、验证门（提交前必跑）

```bash
export PATH="$PWD/.tools/lean-4.19.0/bin:$PATH"
/usr/bin/python3 scripts/verify.py
```

六项：`python-tests`、`pairing`、`lean-build`、`lean-axioms`、`manuscript-static`、`latex-build`。

论文：

```bash
/usr/bin/python3 scripts/build_paper.py    # 两遍编译，处理引用
```

---

## 四、环境事实（省得反复查）

- Lean 4.19.0：`.tools/lean-4.19.0/bin/`，需 `export PATH`
- Tectonic：`.tools/tectonic`
- pyscf 2.15.0 **只在 `/usr/bin/python3`**，默认 `python` 没有
- 无 geomeTRIC；几何优化用 `experiments/caf2/reaction3.py` 自带的 Hessian 播种优化器
- 内存 15 GB 总，约 14 GB 可用；20 核
- `cc.kernel()` 返回**相关能**，总能是 `mycc.e_tot`
- pyscf 的 `frozen` 是**冻结轨道数**，不是电子数
- `mol.atom_coord()` 返回**玻尔**，`gto.M(atom=...)` 默认**埃**
- 频率分析必须投影掉平动转动；计数取 |w| 最大的 `3N−n_proj` 个

---

## 五、学术纪律

- 引用只用 Crossref 实际取回的记录，绝不手写 DOI 或页码
- 已知定理 / 猜测 / 未知必须明确区分
- 数值断言必须能由仓库内脚本复现，JSON 落盘
- 未验证的写进 UNVERIFIED，不写成结论
- 陈述过时立刻撤回并说明，不留模糊措辞

---

## 六、当前状态入口

- `docs/NEXT_SESSION.md`：断点与待办
- `STATUS_REPORT.md`：分层状态总览
- `PAPER_STATUS.md`（若存在）：投稿进度

---

## 七、中断的三类真实原因（2026-11-18 实测）

1. **工具调用格式损坏**（主因）。多余闭合标签、空调用、只写意图不发调用。
   对策：一次一个调用，短命令，发完不留任何标签。
2. **前台命令约 420 秒硬上限**。`sleep 540` 配 `timeout 590` 仍被截断在 420 s。
   对策：长作业一律 `background=true`，前台只做 `sleep <= 300` 的轮询。
3. **OOM**。`free -g` 显示的是机器总量，进程实际可用远小于此。
   实测 4 GB 左右会被杀。对策见下。

## 八、大矩阵计算的内存纪律

14400×14400 float64 就是 1.66 GB，任何一次多余拷贝都会翻倍。

禁止：
- `O = H - np.diag(D)`（`np.diag` 先造一个完整矩阵，再减，两份 1.66 GB）
- `np.abs(L)` 后做 `absL @ absL.T`（又两份）
- 同时持有 H 和 A 和 Cholesky 结果

正确做法：
- 需要 `O` 的作用时用 `LinearOperator(lambda v: H@v - D*v)`，不物化
- 舍入界用 `‖ |L||Lᵀ| ‖_F ≤ ‖L‖_F² = trace(A)`，只要一个标量，
  在 `dpotrf` **之前**取，因为分解会覆盖 A
- `lapack.dpotrf(..., overwrite_a=1)` 原地分解；`dtrtri(..., overwrite_c=1)` 原地求逆
- `dpotrf(lower=1)` 只写下三角，严格上三角仍是原 A 的残留 ——
  **分解后必须 `L[np.triu_indices_from(L,1)] = 0.0`**，否则残差和范数全错
- 求 `‖L⁻¹‖_∞` 时逐行循环取 max，不要造 `np.abs(Linv)` 大数组

分两趟：先建 H 取极值、参考行、对角，释放；再直接按列建 A。
峰值从 5 GB 降到 1.7 GB。

## 九、图必须用统一能量基准

`flat_mode.json` 的 `dE_kcal` 相对**各自极小点**，`cct_barrier.json` 的
`rel_kcal` 相对**线性构型**。混画会把结论说反。
统一到 `E(θ) − E(180°)` 后再画。

## 十、落盘的 JSON 也要验

`route_comparison.json` 的 `temple_true_gap` 条目 slack 为 **正值**
（在 E_FCI 之上），即根本不是下界。画图前必须逐个断言 `value < E_FCI`。

## 十一、图的尺寸必须等于印刷宽度

`\includegraphics[width=\columnwidth]` 会把图拉伸到 15.94 cm。
若 `figsize` 只有 3.3 英寸（8.4 cm），就是**放大 1.9 倍**，
图内 7.5 pt 的字在页面上变成约 14 pt，比 11 pt 正文还大。

规则：`figsize` 宽度一律设成 `6.3` 英寸（= 15.94 cm），保持纵横比只改宽度。
这样显示尺寸不变，字回到标称字号。

还有一个字体因素：DejaVu Serif 的 x-height 比 LaTeX 的 Computer Modern 大，
同样 pt 数看起来更壮。所以图内基准字号取 8 pt（标题 8.5、刻度 7.5）。

核验方法（不要靠肉眼）：

    import pymupdf
    doc = pymupdf.open('build/paper_main.pdf')
    for b in doc[p].get_text('dict')['blocks']:
        for l in b.get('lines', []):
            for sp in l['spans']:
                print(round(sp['size'],1), sp['text'][:40])

正文应变 10.9 pt，图内所有 span 都必须小于它。

## 十二、`tabular` 不能跨页

12 行的长表塞进一页时，超出部分**不会报错**，而是直接画到纸张外面。
本项目 `tab:lean` 就这样把 `ComplementBound` 那一行溢出页底 77 pt（2.7 cm）。

logs 里查不到 "Float too large" —— 因为它不是浮动体问题。
必须用坐标检查：

    M, H = 70.87, 841.9          # 2.5cm 边距
    for i in range(doc.page_count):
        for b in doc[i].get_text('blocks'):
            if b[3] - (H - M) > 36:   # 36 pt 以上才算，32 pt 是页码页脚
                print(i+1, b[3]-(H-M), b[4][:50])

长表一律用 `longtable`（带 `\endfirsthead` / `\endhead` / `\endfoot`），
且不能放在 `table` 浮动体内。

## 十三、窄 `p{}` 列里的长标识符

两端对齐的 `p{}` 列里放 `\texttt{sampling\_gives\_no\_upper\_bound}`
这类标识符时，`\_` 处不能断行，TeX 只能拉伸词距，badness 冲到 10000。
对策：列定义加 `>{\raggedright\arraybackslash}`（需 `\usepackage{array}`），
或者干脆不要在正文里罗列十几个标识符，改用表格。
