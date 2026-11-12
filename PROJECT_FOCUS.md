# 项目聚焦：ALD 前驱体设计与动力学预测

> 审计更正：以下是早期方案，不是已完成结果。用户原记号为 Al2O3/MoCl4 与 CaF；下文 CaF2、Ca(thd)2/HF 由先前助手未经确认引入。O(K⁶)、固定精度及GPC目标均未证明或验证。当前可信范围以 STATUS_REPORT.md 与 docs/audit_numerics.md 为准。

## 目标体系
1. **Al₂O₃/MoCl₄**：Mo 前驱体在羟基化氧化铝表面的反应机理
2. **CaF₂**：Ca(thd)₂ + HF 的 ALD 生长动力学

## 核心问题
- 前驱体吸附能计算（替代传统 DFT）
- 配体交换活化能预测
- 生长速率（GPC）理论估计

## 方法创新
- 2-RDM（3-positivity SDP，O(K⁶)）替代 CCSD(T)（O(N⁷)）
- 张量网络 + 纠缠感知轨道优化降低活性空间
- Lean 4 形式化验证核心定理

## 验证标准
- 与已知 DFT 计算对比（TMA/H₂O 为基准）
- 活化能误差 <0.1 eV（~2 kcal/mol）
- GPC 预测误差 <10%

## 文献对标
- ALD 工艺：JVST-A, Chem. Mater.
- 计算方法：JCTC, J. Phys. Chem.
- 交叉创新：APL, Adv. Mater.

## 数值实验计划
1. **H₂ Richardson 模型**（概念验证）
2. **TMA/Al₂O₃ 反应能**（方法校准）
3. **MoCl₄/Al₂O₃ 新体系**（原创预测）
4. **CaF₂ GPC 估计**（应用验证）
