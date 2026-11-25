# 投稿指南

面向：上海交通大学硕士研究生，本人第一作者，上海交通大学为第一单位。
论文：`paper/paper_main.tex`（34 页，12 图，3 表，31 条参考文献）。

所有期刊指标标注来源。**凡标「需你确认」的，我无法从公开渠道核实，不要当作事实使用。**

---

## 一、毕业硬性条件（已核实）

依据《上海交通大学关于申请授予硕士学位（学术型）的规定》沪交研〔2017〕93 号
（https://gk.sjtu.edu.cn/Data/View/651），以及各学院公开发布的细则：

| 条件 | 内容 | 状态 |
|---|---|---|
| 数量 | 答辩前至少发表**或录用** 1 篇与学位论文主要内容相关的论文 | 已核实 |
| 第一单位 | **上海交通大学为第一单位**（不是"包含"，是排第一） | 已核实 |
| 作者位次 | 本人第一作者计 1 篇；第二作者且第一作者为导师计 1/2 篇；第三作者及以后不计 | 已核实 |
| 质量 | 多数学院细则为 **SCI 或 EI 或中文核心期刊**收录 | 已核实（以工科各学院为代表） |
| 刊物目录 | **由各学位评定分委员会认定**，学院可制定更高标准 | 已核实 |
| 化学化工学院硕士具体目录 | 公开渠道只检索到该院博士标准（要求中科院分区 Q1/Q2）与论文评阅办法，**未找到硕士刊物目录** | **需你确认** |

**行动项：** 向化学化工学院研究生教务或学位秘书索取你入学年份适用的
《硕士学位授予标准》与刊物目录认定表。这是唯一能决定"投哪里算数"的文件。

---

## 二、论文现状对投稿的影响（审稿人视角）

**优势**

- 有严格下界且给出可复算的认证程序（verified Cholesky），不是"近似估计"
- Lean 4 内核验证，并明确列出哪些假设仍由调用者提供
- 报告了负面结果（块证书在两体系均失效），可信度加分
- 第二体系（N₂，14 400 行列式）证明不是单一矩阵的偶然
- 图表全部由落盘 JSON 生成，可复算

**风险，按严重程度排序**

1. **化学内容偏薄。** CaF₂ 部分全是气相模型反应，没有任何表面/周期性计算。
   投纯计算化学刊（如 Comput. Theor. Chem.）会被审稿人直接指出
   "no surface calculation"。→ 优先投数学化学或理论化学刊。
2. **摘要自述"不主张新颖性"。** 原摘要写着 "no novelty is claimed for the inequalities"。
   对 JOMC 这类要求 "original, chemically important mathematical results" 的刊物是自我削弱。
   **已处理**：现摘要改为陈述链条本身与量化失效边界，不再主动弃权，但也没有夸大。
3. **未做 CBS 外推、未做几何优化方法学说明、未用 ECP。** 这是计算化学刊的常见质询点，
   数学化学刊则不问。
4. **参考文献 31 条偏少**，但数学化学类文章可接受。

---

## 三、推荐期刊（全部 SCIE，已逐条核实）

### 第一选：Journal of Mathematical Chemistry（Springer）

| 项目 | 数值 | 来源 |
|---|---|---|
| 索引 | **SCIE**、Scopus、Math. Reviews、zbMATH | Springer 期刊页索引列表 |
| 影响因子 | 2.0（2025），5 年 1.5 | Springer 期刊页 |
| 首次决定 | **中位数 12 天** | Springer 期刊页 |
| 出版模式 | Hybrid（走订阅路线**不付 APC**） | 同上 |
| 评审 | 单盲，通常 2–3 位审稿人，可推荐审稿人 | 投稿指南 |
| 篇幅 | **无页数上限**；摘要 150–250 词；关键词 4–6 | 投稿指南 |
| 栏目 | original papers / review / brief communications / letters | 投稿指南 |

定位原文："publishes original, chemically important mathematical results which use
non-routine mathematical methodologies often unfamiliar to the usual audience of
mainstream experimental and theoretical chemistry journals."

**为什么最匹配**：本论文的贡献形态正是"用主流化学读者不熟悉的数学工具处理化学问题"。
34 页无上限、12 天首决、不强制付费，对毕业时间压力最友好。
近期确有量子化学方法学论文（如 Feniou & Toulouse 的局域轨道工作）。

**风险**：要求数学上的原创性。本论文的原创性应表述为
"证书链条的组合 + 经验证的下界认证程序 + 量化的失效边界"，而非单个不等式。

### 第二选：Theoretical Chemistry Accounts（Springer）

SCIE、Scopus，IF 约 1.5–1.8、Q3，h-index 124（来源：researchcloud / journalindexes / resurchand）。
接受理论方法学稿件，Springer 体系投稿流程相同。第一选被拒后的自然顺位。

### 第三选：International Journal of Quantum Chemistry（Wiley）

SCIE，IF 2.7，Web of Science Q3（来源：wos-journal.info、journalindexes）。
量子化学对口，但 Wiley 侧近年稿件量下降，审稿周期不确定性较高。

### 第四选：Computational and Theoretical Chemistry（Elsevier）

SCIE，Scopus CiteScore 4.8、Q2（来源：journalindexes.com、ScienceDirect 期刊页）。
范围含电子结构、机理、催化。**但侧重"方法学改进 + 应用"**，纯数学链条可能被判化学内容不足。

### 备选（快，但必须付费）

**AIP Advances**（AIP Publishing）：SCIE，Gold OA，**APC 1800 USD**，IF 1.7、Q4，
首次决定平均 30 天、录用 70 天、见刊 94 天（来源：pubs.aip.org 官方指标页与 DOAJ）。
范围宽、接受率高，但**没有经费就不能投**。不要作为首选。

---

## 四、投递顺序建议

```
JOMC（12 天首决，免费）
  └─ 若拒 → 按审稿意见修改 → Theoretical Chemistry Accounts
        └─ 若拒 → International Journal of Quantum Chemistry
              └─ 若拒 → Computational and Theoretical Chemistry
                    └─ 时间告急且有经费 → AIP Advances
```

- **不要一稿多投**，SJTU 学术规范明令禁止，且会被刊物拉黑。
- JOMC 首决仅 12 天，即使被拒损失也小，值得先试。
- 每次被拒后**必须**按意见修改再投，不要原样转发。

---

## 五、投 JOMC 前必须补齐的材料

| 项目 | 现状 | 待办 |
|---|---|---|
| 摘要 150–250 词 | **249 词，已合规** | 完成 |
| 关键词 4–6 个 | 无 | **待加** |
| 作者与单位 | `\author{}` 为空 | **待填**（见第六节） |
| 通讯作者邮箱 | 无 | **待填** |
| ORCID | 无 | 建议注册并填写 |
| Statements & Declarations | 无 | **待加**：Funding、Competing interests、Data availability、Author contributions |
| 图件单独上传 | 现为 PDF 内嵌 | Springer 要求图单独提交，300 dpi 以上；现有 `paper/figures/*.pdf` 为矢量，可直接用 |
| 参考文献格式 | natbib + unsrtnat | Springer 用 Basic 风格，投稿时按其模板转换或交由系统处理 |
| 代码与数据可获得性 | 有 GitHub 仓库 | 声明中给出仓库链接（审稿人可复算，是本论文的强项） |
| Cover letter | 无 | **待写**，见第七节 |

---

## 六、作者与单位署名

JOMC 要求标题页含：作者姓名、单位、通讯作者邮箱、ORCID。

**英文单位写法**（上海交通大学官方英文名）：

```
Shanghai Jiao Tong University, Shanghai 200240, China
```

需精确到学院的写法（**请替换为你实际所属学院**）：

```
School of Chemistry and Chemical Engineering,
Shanghai Jiao Tong University, Shanghai 200240, China
```

**必须遵守的两条**（来自沪交研〔2017〕93 号）：

1. 上海交通大学必须是**第一单位**，即第一作者名下的第一个 affiliation。
2. 你必须是**第一作者**。若导师需挂通讯作者，不影响你的第一作者计 1 篇。

---

## 七、Cover letter 要点（据本论文实际情况写）

- 一句话说明问题：有限基组下的能量证书，前半段（区间算术）成立，后半段（耦合常数）在实测中失效
- 一句话说明贡献：给出失效的精确边界 + 一个避免该常数的单参考 Schur 步 + 经验证的 Cholesky 认证
- 一句话说明为何投 JOMC：数学工具与化学问题交叉，且不依赖表面计算即可给出可核验结论
- 主动声明：不含表面计算、不含实验数据、Lean 验证的边界已如实标注
  （主动说明局限比被审稿人发现更有利）

---

## 八、我无法核实、必须你自己确认的四件事

1. **化学化工学院硕士的刊物目录** —— 决定 JOMC/TCA 是否被认可的唯一依据
2. **你入学年份适用的培养方案** —— 各学院明确提示"不同入学年份要求可能有差异"
3. **导师是否要求通讯作者、是否另有目标刊** —— 导师可制定高于学校的标准
4. **是否有 APC 经费** —— 决定备选刊物（AIP Advances 等纯 OA 刊）是否可用

前两项直接决定投稿策略，建议明天就向学院教务确认。

---

## 九、当前论文与投稿要求的差距汇总

已完成：摘要压缩至 249 词；PDF 与 tex 内容级一致性核验通过；
图件尺寸与正文匹配；表格跨页溢出已修复；作者历史已改为本人。

待完成：关键词、作者与单位、通讯邮箱、ORCID、Statements 四项声明、
Cover letter、图件单独打包、参考文献风格转换。
