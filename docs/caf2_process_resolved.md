# CaF₂ ALD：前驱体问题已解决（文献核验）

日期 2026-11-13。所有 DOI 已经 Crossref 核验，摘要经原文页面取得。

## 原先的卡点

用户对 Ca 前驱体不确定，并提出质疑：「对于 Ca 前驱体是 β-二酮盐还是酰胺盐我不清楚，
这个不应该是含 F 的吗？」

**这个质疑是对的，而且指出了两条不同的真实工艺。** 两者都存在、都已发表，
区别在于氟从哪来。

## 已核验的四条工艺

### 路线 A：酰胺盐 + 外部 HF（本项目 `reaction2.py` 建模的就是这条）

Hennessy J., Rodríguez R., West W., Brandon E.,
*Atomic layer deposition of calcium fluoride for barrier coating applications*,
J. Vac. Sci. Technol. A **41**, 042405 (2023). DOI `10.1116/6.0002563`

- 前驱体：**bis(N,N-di-i-propylformamidinato)calcium(II)**，即 `Ca(amd)₂`
- 共反应物：**无水 HF**
- 温度：**175–250 °C**
- 生长速率：**0.4 Å/cycle（175 °C）→ 0.3 Å/cycle（250 °C）**，随循环数线性
- 产物：纳米晶 CaF₂，193 nm 以上透明
- 衬底：硅；在 SF₆ 等离子体中 Si:CaF₂ 刻蚀比 > 2000:1（300 W）

这条路线的前驱体与我们在 `experiments/caf2/reaction2.py` 里用的
`Ca(formamidinate)₂ + 2 HF → CaF₂ + 2 formamidine` **完全一致**。
所以先前的模型选择是对的，现在有了明确引用。

### 路线 B1：β-二酮盐自带氟（用户直觉指向的这条）

Putkonen M., Szeghalmi A., Pippel E., Knez M.,
*Atomic layer deposition of metal fluorides through oxide chemistry*,
J. Mater. Chem. **21**, 14461 (2011). DOI `10.1039/c1jm11825k`

- 前驱体：**Ca(hfac)₂**，六氟乙酰丙酮钙 —— **配体自带 6 个 F**
- 共反应物：**O₃**
- 温度：300 °C
- 生长速率：**0.3 Å/cycle**
- 产物：RBS 测得计量比 **CaF₂.₁₇**，氧含量 < 5 at%

这条不需要外部 HF，氟由配体提供。所谓「含 F 的前驱体」即指此类。

### 路线 B2：先成氧化物再氟化

同上文献。用传统非氟 β-二酮盐 + 氧源先生长金属氧化物，再用
**Hhfac + O₃** 氟化。作者判断此路线更可行。氧含量低于 RBS 检测限 2 at%，
折射率 MgF₂ 1.429、CaF₂ 1.472、LaF₃ 1.687。

### 最早的工作：Ca(thd)₂ + HF

Pilvi T., Arstila K., Leskelä M., Ritala M.,
*Novel ALD Process for Depositing CaF₂ Thin Films*,
Chem. Mater. **19**, 3387–3392 (2007). DOI `10.1021/cm0629412`

这是 CaF₂ ALD 的第一项研究（该文引言指出此前仅有 Ca-(thd)₂ + HF 一例）。

**对衬底问题直接有用的数据（出自该文引言）：**
- CaF₂ 与 Si 的晶格失配：**室温 0.6%，600 °C 时 1.6%**
- 因此 CaF₂ 可作为 Si 的晶体生长模板 / 缓冲层
- 该文列举的外延体系绝大多数是 **Si(111)**：
  Fe₃Si/CaF₂/Si(111)、CeO₂/CaF₂/Si(111)、ZnO/CaF₂/Si(111) 等
- 也有 CaF₂-(110)/Si(001) 的报道

**注意**：用户提到「衬底实验用 100Si」。文献里的外延体系以 Si(111) 为主，
Si(001) 上报道的是 CaF₂-(110)。所以「Si(100) 上继续生长 CaF₂」需要单独确认
取向关系，不能直接沿用 (111) 的结果。

## 综述

Atosuo E., Mäntymäki M., Ritala M.,
*ALD of Metal Fluorides–Potential Applications and Current State*,
Adv. Mater. Interfaces **12** (2024). DOI `10.1002/admi.202400372`

Hennessy J., Jewell A., Nikzad S.,
*(Invited) ALD Metal Fluorides for Ultraviolet Filter and Reflective Coating Applications*,
ECS Trans. **80**, 107–118 (2017). DOI `10.1149/08003.0107ecst`

## 对本项目的直接影响

1. **`reaction2.py` 的模型反应已获文献支撑**：Ca(amd)₂ + 2 HF 就是路线 A 的
   半反应基础，不再是凭空的模型。
2. **GPC 有了实测锚点**：0.3–0.4 Å/cycle（175–250 °C）。论文里此前用
   0.3–0.4 Å/cycle 做覆盖度换算，现在这个区间有直接文献来源。
3. **衬底取向需要重新确认**：文献以 Si(111) 为主，用户提的是 Si(100)。
4. **氟源有两类**，后续若要比较机理，应当分别处理：
   - 外部 HF 质子化配体（路线 A）
   - 配体内部 C–F 键断裂（路线 B1/B2，需要更强的氧化条件）

## 仍未确定

- Si(100) 上 CaF₂ 的具体取向关系与界面结构
- 两种氟源路线在原子尺度的基元步骤与势垒（均无计算数据）
- 表面终止：CaF₂(111) 是 F 终止这一点在文献中常见，但本项目尚未取得
  可引用的直接证据
