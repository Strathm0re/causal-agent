# CausalQC-Agent

**LLM agent that adjudicates competing mechanism hypotheses in molecular materials using quantum-chemistry calculations as interventional data (do-operators).**

> 状态：规划阶段（设计树已闭合，2026-09-08）· 主文档见 [PLAN.md](PLAN.md)

## 一句话

让 agent 把"构效关系的竞争性解释"翻译成因果图，主动选择信息增益最大的量子化学干预（结构编辑 / CDFT·外场钳制 / 协议干预），用计算结果裁决哪个机制存活——把小样本材料数据上不稳定的 SHAP 式相关性解释，升格为因果裁决。

## 核心设计（详见 PLAN.md）

- **因果域**：分子/电子结构域（域 A）；干预分类学 = 结构硬干预 / 模拟器内钳制 / 协议干预
- **主案例**：D-A 二聚体 S1 电荷转移，裁决 η（轨道重叠）vs 驱动力（IP−EA）vs 几何（堆积）三假设
- **引擎**：RTX 4090 + GPU4PySCF（干预主力）→ HPC Gaussian（终审仲裁）；xTB/CREST 粗筛
- **架构**：单智能体 ReAct + 显式假设登记簿（Hypothesis Registry）
- **干预选择**：两层保真，acquisition = EIG / expected CPU-hour（成本感知 BOED）
- **验证**：L1 Hückel 重发现 → L2 Hammett σ 复现 → L3c 光热对偶设计（反事实判别性预测）；无湿实验
- **目标**：Digital Discovery（EiC: Alán Aspuru-Guzik）· arXiv 2027-02 · 投稿 2027-03

## 仓库结构（规划）

```
PLAN.md            # 方案备忘录（决定清单 Q1–Q20、理论、架构、实验矩阵、时间线）
docs/              # 文献笔记、设计讨论记录
src/
  engine/          # PySCF/xTB/Gaussian 封装：JSON in → 指标 out
  registry/        # 假设登记簿（因果图 × 后验 × 证据日志）
  selector/        # EIG/cost 干预选择器（多保真 GP）
  agent/           # ReAct 循环、干预菜单、LLM 配置
experiments/
  l1_huckel/       # Hückel 玩具闭环
  l2_hammett/      # σ 常数重推导
  main_da/         # D-A 二聚体主裁决战役
  shap_control/    # TADF 档案 SHAP 不稳定性对照组（图 1）
  l3c_photothermal/# 光热对偶应用
```

## 时间线锚点

2026-09 infra 骨架 → 2026-11 主裁决战役完成 → 2027-02 中 arXiv → 2027-03 底投稿

## 相关项目

- MFBO / D-A 共晶窄带隙（分子层因果结论为其输送结构先验）
- 光热知识图谱（论文 2 种子：因果规则驱动分子生成 + 湿实验闭环）
