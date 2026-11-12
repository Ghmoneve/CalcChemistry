# 指定八篇文献核验审计

本审计只覆盖 `references/verified_sources.bib` 中的八篇文献，使用指定 BibTeX 键及原始 DOI。元数据来自本次实际 Crossref REST 请求及出版商页面提取，没有将现有 `references.bib` 作为证据。检索日期为 2026-11-12（实际系统时钟记录：2026-11-12T14:30:45Z）。八篇均已取得题名、作者、期刊、年、卷期、页码或文章号；内容证据包括六篇摘要，以及 Richardson、Mazziotti 的研究说明／引言。摘要证据不等同于逐项核查全文定理假设；下文逐项列明可支持的范围与限制。

## 核验结果

| BibTeX 键 | 已核验作者 | 已核验题名 | 期刊、年、卷(期)、页码或文章号 |
|---|---|---|---|
| OGorman2022 | Bryan O’Gorman; Sandy Irani; James Whitfield; Bill Fefferman | Intractability of Electronic Structure in a Fixed Basis | PRX Quantum, 2022, 3(2), 020322（文章号） |
| Liu2007 | Yi-Kai Liu; Matthias Christandl; F. Verstraete | Quantum Computational Complexity of the N-Representability Problem: QMA Complete | Physical Review Letters, 2007, 98(11), 110503（文章号） |
| Mazziotti2012 | David A. Mazziotti | Two-Electron Reduced Density Matrix as the Basic Variable in Many-Electron Quantum Chemistry and Physics | Chemical Reviews, 2012, 112(1), 244–262 |
| Dukelsky2004 | J. Dukelsky; S. Pittel; G. Sierra | Colloquium: Exactly solvable Richardson-Gaudin models for many-body quantum systems | Reviews of Modern Physics, 2004, 76(3), 643–662 |
| Richardson1963 | R. W. Richardson | A restricted class of exact eigenstates of the pairing-force Hamiltonian | Physics Letters, 1963, 3(6), 277–279 |
| Eyring1935 | Henry Eyring | The Activated Complex in Chemical Reactions | The Journal of Chemical Physics, 1935, 3(2), 107–115 |
| Hagedorn1988 | George A. Hagedorn | High order corrections to the time-independent Born-Oppenheimer approximation II: Diatomic Coulomb systems | Communications in Mathematical Physics, 1988, 116(1), 23–44 |
| Pieck2025 | Fabian Pieck; Ralf Tonner-Zech | Computational Ab Initio Approaches for Area-Selective Atomic Layer Deposition: Methods, Status, and Perspectives | Chemistry of Materials, 2025, 37(9), 2979–3021 |

以上字段分别由下列同名键的 Crossref URL 实际返回，原始字段摘录嵌入本文件，以便只写入用户指定的两个文件。完整 JSON 可通过相应 URL 再取；这里不冒称另有保存的完整响应文件。

## Crossref 证据与摘要核验

### OGorman2022

- DOI：`10.1103/PRXQuantum.3.020322`。
- 元数据：[Crossref](https://api.crossref.org/works/10.1103/PRXQuantum.3.020322)。返回字段摘录：`"volume":"3", "issue":"2", "article-number":"020322", "published":{"date-parts":[[2022,5,2]]}`。
- 已提取：[APS 摘要及正文](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.3.020322)。摘要原文：“We prove that the electronic-structure problem, when restricted to a fixed single-particle basis and a fixed number of electrons, is QMA-complete.”
- 摘要还陈述固定基下最低 Slater determinant 能量估计的 NP-complete 结果。证明的复杂性输入编码使用所选空间轨道；“固定”表示给定实例的基与电子数，不应误解成整个问题族的维数恒定。
- 可支持：一般问题族的最坏情形复杂性限制。不能据此证明某个实际 ALD 实例必然困难；QMA-complete 本身也不是已证明 BQP 与 QMA 分离，更不表示所有特例或误差证书均不可高效处理。
- 提取缓存：`/home/moneve/.hermes/cache/web/journals.aps.org-2a110315f3.md`。提取正文有数学符号缺失，不能用该纯文本代替原版公式。

### Liu2007

- DOI：`10.1103/PhysRevLett.98.110503`。
- 元数据：[Crossref](https://api.crossref.org/works/10.1103/PhysRevLett.98.110503)。返回字段摘录：`"volume":"98", "issue":"11", "article-number":"110503", "published":{"date-parts":[[2007,3,16]]}`。
- 已提取：[DOI / APS 页面](https://doi.org/10.1103/PhysRevLett.98.110503)，含作者、刊期与摘要。
- 摘要原文：“We show that this problem is quantum Merlin-Arthur complete, which is the quantum generalization of nondeterministic polynomial time complete.”
- 可支持：N-representability 的复杂性背景。不可据此声称每个实际化学实例都困难，或任何固定层级的必要约束自动成为充分条件。

### Mazziotti2012

- DOI：`10.1021/cr2000493`。
- 元数据：[Crossref](https://api.crossref.org/works/10.1021/cr2000493)。返回字段摘录：`"volume":"112", "issue":"1", "page":"244-262", "published-print":{"date-parts":[[2012,1,11]]}, "published-online":{"date-parts":[[2011,8,25]]}`。
- 已提取：[DOI / ACS 页面](https://doi.org/10.1021/cr2000493)，明确显示 “Cite this: Chem. Rev. 2012, 112, 1, 244–262”。BibTeX 使用卷期年 2012，不能把在线年 2011 当作该卷期年份。
- 摘要状态：[ACS PDF 入口实际返回的网页](https://pubs.acs.org/doi/pdf/10.1021/cr2000493) 明确说明：“In lieu of an abstract, this is the article's first page.” 因而本审计不编造正式摘要。Crossref 也未返回摘要；首次 DOI 提取器标记 Abstract 的内容实际是网页导航。
- 内容证据：[PubMed 入口的提取结果](https://pubmed.ncbi.nlm.nih.gov/21863900/) 返回带有 `10.1021/cr2000493 | Chem. Rev. 2012, 112, 244-262` 页脚的论文文本。其引言原文：“In this work, we review recent advances in one of these new approaches to the study of electron correlation, the direct calculation of the two-electron reduced density matrix (2-RDM).”
- 引言同时区分 contracted Schrödinger equation 方法和 2-RDM 能量变分最小化，后者又区分必要 N-representability positivity 条件与参数化近似。可支持：2-RDM 方法及约束的综述背景。不可支持：有限必要约束等于全部充分条件，或未经数值误差控制的优化输出自动是严格证书。
- 来源限制：上述工具返回的内容超出通常 PubMed 记录，可能来自提取服务关联的全文；本审计只声明读取了该入口返回的论文文本，不声称 NLM 托管或认证了完整正文。引文可在缓存 `/home/moneve/.hermes/cache/web/pubmed.ncbi.nlm.nih.gov-a46b5bda93.md` 查找；ACS“无摘要、以首页代替”的说明缓存于 `/home/moneve/.hermes/cache/web/pubs.acs.org-8da7cffa70.md`。首页图片另行访问返回 403，未能视觉核验。内容证据的来源强度低于直接取得出版商原版全文。

### Dukelsky2004

- DOI：`10.1103/RevModPhys.76.643`。
- 元数据：[Crossref](https://api.crossref.org/works/10.1103/RevModPhys.76.643)。返回字段摘录：`"volume":"76", "issue":"3", "page":"643-662", "published":{"date-parts":[[2004,8,6]]}`。
- [APS](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.76.643) 提取报告付费限制；已实际读取对应 [arXiv 预印本摘要](https://arxiv.org/abs/nucl-th/0405011)。直接 HTTP 读取其 HTML 元数据返回 `citation_doi=10.1103/RevModPhys.76.643`，作者依次为 J. Dukelsky、S. Pittel、G. Sierra，投稿历史只有 v1（2004-05-05）；固定版本入口为 [nucl-th/0405011v1](https://arxiv.org/abs/nucl-th/0405011v1)。期刊元数据仍使用 Crossref，而不是预印本题名或发布日期。
- 摘要原文：“The use of exactly-solvable Richardson-Gaudin (R-G) models to describe the physics of systems with strong pair correlations is reviewed.”
- 可支持：Richardson–Gaudin 配对可积模型、强配对关联及其应用的综述；摘要还讨论纯配对 Hamiltonian 的大 N 极限与 BCS 理论。不可支持：任意分子电子 Hamiltonian 与该模型严格等价，或未经构造与误差估计就把可积解当作真实表面能量。
- 题名中的 Colloquium 属于出版记录，予以保留；Crossref 的 HTML 斜体标签只作格式清理。

### Richardson1963

- DOI：`10.1016/0031-9163(63)90259-2`。
- 元数据：[Crossref](https://api.crossref.org/works/10.1016%2F0031-9163%2863%2990259-2)。返回字段摘录：`"volume":"3", "issue":"6", "page":"277-279", "published":{"date-parts":[[1963,2]]}`。
- 已提取：[DOI 返回的扫描正文文本](https://doi.org/10.1016/0031-9163(63)90259-2)。Crossref 没有摘要；以下使用论文开头的研究说明，不伪称正式摘要。
- 原文：“In this note we will show that there exists a restricted class of eigenstates of the pair-force Hamiltonian which contain N pairs of particles and which are natural generalisations of those states which contain one pair.”
- 可支持：特定 pairing-force Hamiltonian 的受限精确本征态构造。不可支持：任意分子电子 Hamiltonian 都是该可积模型。
- 提取质量：扫描文本开头混入前一篇文章、公式存在 OCR 损坏；仅选取辨识明确的本篇文字，公式必须另查清晰原版。

### Eyring1935

- DOI：`10.1063/1.1749604`。
- 元数据及摘要：[Crossref](https://api.crossref.org/works/10.1063/1.1749604)。返回字段摘录：`"volume":"3", "issue":"2", "page":"107-115", "published":{"date-parts":[[1935,2,1]]}`。
- Crossref 摘要原文：“The probability of the activated state is calculated using ordinary statistical mechanics. This probability multiplied by the rate of decomposition gives the specific rate of reaction.”
- 已提取：[DOI 对应论文正文](https://doi.org/10.1063/1.1749604)。正文明确讨论配分函数与反复越垒修正：“For some reactions it will happen that the same activated complex may cross the barrier and return without decomposing. This fact reduces the actual reaction rate.”
- 可支持：活化复合物与统计力学速率框架。不能仅凭两个电子能量之差声称得到绝对真实反应速率；温度、配分函数及动力学校正须另行处理。
- 提取器缓存：`/home/moneve/.hermes/cache/web/doi.org-f5d79249b7.md`（提取工具产生的完整文本缓存，非本次新增仓库文件）。该扫描提取尾部混入下一篇文章，不能不加区分地引用。

### Hagedorn1988

- DOI：`10.1007/BF01239023`。
- 元数据：[Crossref](https://api.crossref.org/works/10.1007%2FBF01239023)。返回字段摘录：`"volume":"116", "issue":"1", "page":"23-44", "published":{"date-parts":[[1988,3]]}`。
- 已提取：[Springer 摘要](https://link.springer.com/article/10.1007/BF01239023)。直接 DOI 提取最初报告付费限制，文章入口随后成功返回摘要（返回 URL 带 `error=cookies_not_supported` 查询参数，仍有实际题名、摘要及参考文献；此处引用稳定入口）。
- 摘要原文：“We study the bound states of diatomic molecular systems.” 随后陈述：“We prove that if the nuclear masses are proportional to ε−4 then certain eigenvalues and eigenvectors of the Hamiltonian have asymptotic expansions to arbitrarily high order in powers of ε, as ε→0.” 这里 ε−4 是提取器对上标 ε 的负四次方的纯文本呈现。
- 可支持：在该文双原子体系条件下，某些束缚态本征值与本征矢随核质量缩放的高阶渐近展开；摘要还说明本征值展开零至四阶为 Born–Oppenheimer 项，第五阶为零。摘要不提供任意 ALD 表面体系可直接套用的显式统一数值误差常数。
- 范围约束：题名明确限定 time-independent、II、diatomic Coulomb systems；不能将其冒用为适用于任意表面反应的无条件 Born–Oppenheimer 误差定理。

### Pieck2025

- DOI：`10.1021/acs.chemmater.4c03477`。
- 元数据：[Crossref](https://api.crossref.org/works/10.1021/acs.chemmater.4c03477)。返回字段摘录：`"volume":"37", "issue":"9", "page":"2979-3021", "published-print":{"date-parts":[[2025,5,13]]}, "published-online":{"date-parts":[[2025,4,4]]}`。
- 已提取：[ACS PDF](https://pubs.acs.org/doi/pdf/10.1021/acs.chemmater.4c03477)，包括题名、两位作者、刊期、摘要和部分正文。摘要原文：“Special emphasis is placed on the use of ab initio methods to explore surface chemistry, optimize precursor and inhibitor properties, and improve selectivity.”
- 可支持：AS-ALD 的计算方法、表面相互作用、前驱体／抑制剂设计及选择性建模背景。正文讨论以量子化学反应势垒和能量作为 kinetic Monte Carlo 的输入，并说明反应器流动与扩散可能需要更大尺度模型；这提供多尺度建模的背景，不是本项目的数值验证。
- 文献类型：领域方法综述；不表示本仓库已完成任何 ALD 反应计算，也不提供任意电子能量近似到真实 ALD 速率的无条件严格证书。
- 来源限制：DOI 页面最初返回 403；PDF 提取随后成功，尽管返回结果的 URL 字段被提取服务归并为 `https://pubs.acs.org/`，正文题名、作者、页脚 DOI 与指定记录一致。缓存为 `/home/moneve/.hermes/cache/web/pubs.acs.org-9221f2eee4.md`；它是服务返回的文本，长度截到 50,000 字符，不能声称已取得完整 43 页全文。PDF 内已有的历史下载时间不是本次核验时间。

## 文件验收

已实际运行只读 Python 校验，得到以下结果：

- BibTeX 恰有 8 个互不重复的指定键；DOI 原始拼写和 DOI URL 均与任务输入逐字一致；作者、题名、期刊、年、卷、期、页码／文章号字段齐全，花括号嵌套平衡。
- 对落盘 BibTeX 的八篇记录逐条重新查询 Crossref，比对 `author/title/journal/year/volume/issue/pages-or-article-number/DOI`，全部 `mismatches: {}`。首次重查在 Richardson 处遇到 429，改用先前成功的编码 URL 后，余下四篇比对也全部通过。
- 比对时仅规范空白、HTML／MathML 标签、BibTeX 大小写保护括号、页码双连字符和弯／直撇号；DOI 匹配 Crossref 时忽略大小写，但文件内另行严格核对并保留用户原始拼写。作者缩写没有根据记忆擅自展开；`R.W.` 与 `R. W.` 只作空白排版处理。
- 中文审计有对应的 8 个条目，无“待补查”占位；保留 `Al2O3/MoCl4` 与 `CaF`。
- 从五段有本地提取缓存的引文逐字核对原文（仅折叠空白），全部通过；其余摘要／开篇引文来自本次实际工具响应，已在各条目给出入口和摘录。
- 校验是字段、结构及证据一致性检查，不冒称已执行 LaTeX／BibTeX 排版编译，也不冒称已验证每篇论文证明。

本次仓库写入仅为 `references/verified_sources.bib` 与 `docs/reference_audit.md`；既有工作区修改不属于此次文献审计。

## 使用边界与保留事项

1. 原始用户术语 `Al2O3/MoCl4`、`CaF` 保留原样。这里不将 `MoCl4` 改成别的氯化物，也不将 `CaF` 擅改为 `CaF2`；它们的物种、相态与角色需独立定义。
2. 本次只核验文献，不宣称任何具体 ALD 反应、过渡态、势垒、选择性或生长速率已经计算。
3. 从电子能量证书到 ALD 速率需要另外明确 Hamiltonian、离散化与模型误差、几何和过渡态、自由能修正、动力学模型适用条件与参数界限；这些要求是本审计的使用限制，不冒称由上述某一篇论文提供了完整端到端证明。
4. `020322`、`110503` 是文章号，BibTeX 为常见样式兼容放入 `pages`，不是页数。
5. 初次对 Crossref 的编码 DOI 批量请求有六条返回 429；后续实际请求取得全部八条元数据。检索遇到的访问错误和缺失摘要如实列出，不以猜测补齐。
