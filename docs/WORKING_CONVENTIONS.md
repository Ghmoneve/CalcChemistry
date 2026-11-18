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
