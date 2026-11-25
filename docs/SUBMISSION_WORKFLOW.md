# 首次投稿完整流程

对象：董京瑞（Jingrui Dong），上海交通大学机械与动力工程学院，硕士毕业用。
目标刊：Journal of Mathematical Chemistry（Springer），备选 Theoretical Chemistry Accounts、Int. J. Quantum Chem.
当前稿件：34 页主文件 / 35 页含声明，12 图 3 表，31 条参考文献，摘要 249 词。

---

## 零、时间表（先看这个，决定你什么时候开始）

| 阶段 | 乐观 | 现实 | 保守 |
|---|---|---|---|
| 材料准备 | 3 天 | 1 周 | 2 周 |
| 首决（JOMC 中位数 12 天） | 12 天 | 3 周 | 6 周 |
| 小修/大修 + 返修 | 3 周 | 6 周 | 10 周 |
| 接收后在线发表 | 2 周 | 4 周 | 8 周 |
| **合计** | **约 2 个月** | **约 3.5 个月** | **约 6 个月** |

**如果被拒**，重投下一家要从头走一遍，加 2–3 个月。
**所以最迟在答辩前 6 个月必须投出去。** 你现在就该确认自己是几月答辩。

SJTU 规定是"发表**或录用**"，拿到录用通知即可，不必等见刊。这能省 1–2 个月。

---

## 一、投前必须先做的三件事（顺序不能反）

### 1. 确认刊物在机动学院的认可目录内

这是唯一能决定"投中算不算数"的事。学校的原话是"刊物目录由各学位评定分委员会认定"。

- 找谁：机动学院研究生教务办，或你的学位秘书
- 问什么：「我入学年份适用的硕士学位授予标准里，Journal of Mathematical Chemistry
  是否在认可目录内？SCI 收录是否自动认可？」
- 一并确认：是否需要**录用通知**即可，还是要见刊

**没有确认之前不要投。** 投中但不算数，是最坏的结果。

### 2. 告知导师并取得同意

你是唯一作者兼通讯作者，这是合法的。但如果工作用了课题组的资源、经费或设备，
投前必须告知导师。部分单位要求导师签署投稿同意书。

建议这样说：「我准备投 Journal of Mathematical Chemistry，Springer 的，
内容是我做的证书理论和 CaF₂ 的模型反应计算。想先跟您说一声，看您有没有意见。
如果您希望署名或挂通讯，我按您的意思改。」

### 3. 填完稿件里剩余的两个红色占位符

打开 `paper/paper_main.tex`，搜索 `fillme`，还剩两处：

**（a）Funding（第 549 行附近）**

如果这项工作来自课题组基金，必须写明基金名称和编号，例如：

```latex
\paragraph{Funding.}
This work was supported by the National Natural Science Foundation of China
(grant number XXXXXXXXX).
```

如果没有专项基金：

```latex
\paragraph{Funding.}
The author received no specific funding for this work.
```

**注意**：用了课题组经费却写"no specific funding"是不实陈述，被查出会影响学术诚信记录。
拿不准就问导师。

**（b）Use of generative artificial intelligence（第 567 行附近）**

这条必须如实写。本项目在开发过程中使用了 AI 辅助（代码编写、绘图、文字编辑）。
建议措辞：

```latex
\paragraph{Use of generative artificial intelligence.}
Generative artificial intelligence tools were used to assist with software
development and language editing. The author independently verified every
numerical result, every proof step and every bibliographic entry, and takes
full responsibility for the content of the manuscript.
```

如果你实际使用范围与此不同，按实际情况改。**Springer 明确要求披露，隐瞒的风险远大于披露。**

---

## 二、文件准备清单

### 需要上传的文件

| 文件 | 在仓库的位置 | 说明 |
|---|---|---|
| 主 tex | `paper/paper_main.tex` | |
| 章节 tex | `paper/sections/*.tex`（8 个） | 全部要传 |
| 参考文献 | `build/paper_main.bbl` | **传 .bbl，不用传 .bib** |
| 图 | `paper/figures/*.pdf`（12 个） | 矢量，已符合 300 dpi 要求 |
| 编译好的 PDF | `build/paper_main.pdf` | 供编辑快速浏览 |
| Cover letter | 见第五节 | 单独一个 PDF |

### 上传前的自检

- [ ] 搜索 `fillme`，确认 0 处
- [ ] 搜索 `TODO`、`XXX`，确认 0 处
- [ ] 编译无 error，页数正确
- [ ] 所有图在 PDF 里正常显示（翻一遍）
- [ ] 摘要 150–250 词（现在是 249）
- [ ] 关键词 4–6 个（现在 6 个）
- [ ] 通讯邮箱正确：jingruidong@sjtu.edu.cn
- [ ] 单位顺序：Shanghai Jiao Tong University 是**第一个**
- [ ] 参考文献 31 条，全部能在正文找到引用位置

### 建议加行号（审稿人友好）

Springer 建议投稿稿件带行号，方便审稿人指出位置。在 `paper_main.tex` 导言区加：

```latex
\usepackage{lineno}
\linenumbers
```

接收后终稿再去掉。这个不是必须，但能减少"审稿意见指不清位置"的麻烦。

---

## 三、投稿系统的实际操作（逐屏）

JOMC 用的是 Springer Nature 的投稿系统。

### 步骤 1：注册账号

打开 https://submission.nature.com/new-submission/10910/3

用 **jingruidong@sjtu.edu.cn** 注册。不要用 QQ 邮箱——部分期刊系统会过滤。
注册后先验证邮箱，否则收不到通知。

### 步骤 2：选择文章类型（Article Type）

选 **Original Paper**。不要选 Review（那是综述）、Brief Communication（那是短讯）。

### 步骤 3：填元数据

- **Title**：Energy Certificates and Conditional Kinetic Bounds for Finite-Basis
  Surface-Reaction Models
- **Abstract**：从 PDF 里复制粘贴。系统按纯文本存，公式会变形，**提交前在系统里读一遍**，
  确认 441、1.75–2.03、1241、0.860 这些数字没被吃掉
- **Keywords**：certified energy bounds; interval arithmetic; Schur complement;
  certified Cholesky; Lean 4; atomic layer deposition
- **Authors**：Jingrui Dong，机构填 School of Mechanical Engineering,
  Shanghai Jiao Tong University, Shanghai 200240, China
  - 标记为 **Corresponding Author**
  - ORCID 建议填（没有就先去 orcid.org 免费注册，10 分钟）

### 步骤 4：上传文件

系统会让你在下拉框里给每个文件指定类型：

| 你的文件 | 选这个类型 |
|---|---|
| paper_main.tex + sections/*.tex | Manuscript |
| paper_main.bbl | Manuscript（或 LaTeX Support File） |
| figures/*.pdf | Figure |
| cover_letter.pdf | Cover Letter |

**顺序很重要**：把主 tex 放在第一位。系统最后会合成一个 PDF，顺序会影响编辑的第一印象。

### 步骤 5：推荐审稿人（Suggested Reviewers）

JOMC 允许推荐。建议填 3 位，**必须是与你不存在利益关系的**：
- 不要填你的导师、同门、近期合作者
- 从你参考文献里找引用过的、不同国家/机构的作者最合适
- 从 31 条参考文献的作者里挑，这是最省力也最合理的做法

如果实在找不到合适的，可以留空，不影响送审。

### 步骤 6：声明（Declarations）

系统会逐条问你：
- 是否一稿多投 → **必须选 No**
- 是否已发表或部分发表 → 如实回答
- 是否有利益冲突 → No
- 是否使用 AI → **如实填写，与稿件里的声明一致**

### 步骤 7：审核生成的 PDF（最容易出错的一步）

系统会把你上传的文件合成一个 PDF，**必须逐页看完**。重点看：

- 图有没有丢失或错位
- 公式有没有变成乱码
- 表格有没有跨页截断
- 参考文献编号对不对

发现问题就返回修改重新生成，不要抱着"应该没事"的心态点提交。
**这一步出错是最常见的退稿原因，而且是自己造成的。**

### 步骤 8：提交

点 Submit 后会立刻收到确认邮件和稿件编号（形如 JOMC-D-26-00XXX）。
**把这个编号保存好**，之后所有通信都要引用它。

---

## 四、投稿后的状态与应对

### 状态含义

| 状态 | 含义 | 你要做什么 |
|---|---|---|
| Submitted | 已提交 | 等 |
| With Editor / Editor Assigned | 编辑在处理 | 等 |
| Reviewers Invited | 在找审稿人 | 等，**这步可能卡 2–4 周** |
| Under Review | 审稿中 | 等 |
| Required Reviews Complete | 意见已回齐 | 等编辑决定 |
| Decision in Process | 编辑在写决定 | 等 |

**不要催稿。** 首决中位数 12 天是 JOMC 的官方数字，但"找审稿人"阶段常超时。
超过 6 周没动静，可以礼貌发一封询问信给编辑。

### 四种结果

**1. Accept（直接接收）** —— 极少见，别期待

**2. Minor Revision（小修）** —— 最好的常见结果
- 通常给 2–4 周
- 逐条回复每个意见，说明改了什么、在哪一页
- 不同意的意见可以反驳，但要给理由，不要沉默略过

**3. Major Revision（大修）**
- 通常给 4–8 周
- 需要补计算或补论证。如果要求补表面计算，那正是本论文的弱点，要如实说明
  范围限制，或补充能做的部分

**4. Reject（拒稿）**
- 拿到意见后**先放 3 天**，不要立刻改投
- 逐条看意见：哪些是真问题（必须改），哪些是审稿人偏好（可以不改）
- 改完再投下一家。**原样转发会被下一家识别**

### 回复信的写法（Response to Reviewers）

每一条按这个格式写：

```
Reviewer 1, Comment 3:
"<原文照抄审稿人的话>"

Response: We thank the reviewer for this point. <你的回应>
We have revised Section 4.2 (page 12, lines 340-352) accordingly.
```

三条铁律：
1. 每条都回，不许跳过
2. 改了就要给出具体的页/行位置
3. 不同意就说明理由，态度要客气但立场要清楚

---

## 五、Cover Letter 草稿

单独存成 PDF 上传。以下可直接用，按需改：

```
Dear Editor,

I am pleased to submit my manuscript "Energy Certificates and Conditional
Kinetic Bounds for Finite-Basis Surface-Reaction Models" for consideration
as an Original Paper in the Journal of Mathematical Chemistry.

The manuscript asks what happens when an energy-certificate chain is measured
rather than assumed. On systems where every determinant can be diagonalised
exactly, the block certificate fails: its coupling constant exceeds the
correlation energy it is meant to resolve by a factor of forty, and I prove
that no change of basis can repair this, since a similarity transform preserves
the spectrum. I then show that a single-reference Schur step avoids the
constant entirely, provided one positive-definiteness hypothesis is certified
by a verified Cholesky factorisation with an explicit rounding bound. On
N2 with 14,400 determinants the same chain is repeated, and the slack grows,
which I report rather than suppress.

I believe this fits the journal's scope because the contribution is a
mathematical result of chemical relevance obtained with tools that are not
standard in mainstream theoretical-chemistry venues, and because its negative
half is stated as precisely as its positive half.

I should be explicit about the limits of the work: it contains no surface
calculation and no experimental data. The CaF2 system is treated as a gas-phase
model reaction, and the deposition arithmetic is applied to published growth
data. The Lean 4 development is reported with an explicit statement of which
hypotheses remain supplied by the caller. The repository at
https://github.com/Ghmoneve/CalcChemistry regenerates every number and figure.

This manuscript is not under consideration elsewhere, has not been published
previously, and I am the sole author. I have no competing interests to declare.

Thank you for your consideration.

Sincerely,
Jingrui Dong
School of Mechanical Engineering
Shanghai Jiao Tong University, Shanghai 200240, China
jingruidong@sjtu.edu.cn
```

**为什么这么写**：主动说明"没有表面计算"看起来像自曝其短，但审稿人一定会发现，
由你先说、并说明为什么这不影响结论的价值，比被审稿人指出要好得多。
这是有经验的作者的标准做法。

---

## 六、录用之后

1. **保存录用通知**（PDF），这是申请学位的凭证
2. 填出版协议（Licence to Publish），选订阅路线则**不需要付费**
3. 校对 proof：Springer 会发排好版的校样，**必须在 48 小时内逐字核对**
   - 重点核对：你的名字拼写、单位名称、公式编号、图是否错位
   - 这一步之后再改要发勘误，很麻烦
4. 提交最终版，等待在线发表（DOI 会先出来，DOI 出来就算正式发表）
5. 把 DOI 和录用通知交给学院，办学位申请

---

## 七、新手最容易踩的十个坑

1. **没确认刊物是否在学院目录内就投** —— 投中也不算数，最严重的错误
2. **一稿多投** —— 会被所有相关期刊拒稿并通报学校，后果远超延期
3. **图表只放在 PDF 里，没单独上传** —— 编辑部会退回让你重传
4. **不检查系统生成的 PDF** —— 图丢失、公式乱码，多数是这一步没看
5. **摘要超字数** —— JOMC 上限 250 词，超了系统会拒收
6. **单位顺序写错** —— SJTU 必须是第一个 affiliation
7. **AI 使用声明不实或不写** —— Springer 强制要求披露
8. **被拒后原样改投** —— 同一批审稿人可能再抽到，会被认为不认真
9. **审稿意见挑着回** —— 漏回一条，编辑会退回让你补全
10. **校样不逐字核对** —— 名字或单位拼错，见刊后只能发勘误

---

## 八、你现在的具体待办

**本周内：**

- [ ] 找研究生教务确认 JOMC 是否在认可目录内（最重要）
- [ ] 告知导师并取得同意
- [ ] 决定 Funding 怎么写（问导师是否有基金）
- [ ] 确认 AI 使用声明的措辞与实际一致
- [ ] 注册 ORCID

**确认之后：**

- [ ] 填完两处 `fillme`
- [ ] 加行号（可选）
- [ ] 编译，逐页看一遍 PDF
- [ ] 生成 cover letter PDF
- [ ] 注册投稿系统账号并验证邮箱
- [ ] 按第三节逐屏提交

**提交后：**

- [ ] 保存稿件编号
- [ ] 每 3–5 天看一次状态，不要天天刷
- [ ] 6 周无进展再发询问信

---

## 九、我能帮你做的和不能替你做的

**能帮：** 改稿、补计算、生成图、检查格式、写回复信、按审稿意见修改论文、
核对参考文献、重新编译、检查 PDF 与源码一致性（`scripts/check_pdf_sync.py`）。

**不能替：** 确认学院目录、决定基金表述、确认 AI 使用范围、与导师沟通、
判断某条审稿意见该不该接受。这些涉及你本人和单位的责任，必须你自己做。
