# 来源、版本与采用范围

## S1：ALPHA-g 已发表读数

E. K. Anderson et al. (ALPHA Collaboration), *Observation of the effect of gravity on the motion of antimatter*, Nature **621**, 716–722 (2023), version of record 2023-09-27。

- [DOI / 原刊](https://doi.org/10.1038/s41586-023-06527-1)
- [作者 Purdue PDF](https://www.physics.purdue.edu/~robichf/papers/nature621_716.pdf)
- [PMC 同文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10533407/)
- [实验组出版页](https://alpha.web.cern.ch/observation-effect-gravity-motion-antimatter/)

本次直接读取作者 PDF 的文本，不取图上点，不存大 PDF 或图像。网页直接打开有访问失败，但搜索缓存和作者 PDF 可核；没有绕过访问控制。该文为 CC BY 4.0，本件为有标注的中文摘要及独立符号整理，并非原文逐字转录。

| 精确位置 | 用于本件 |
|---|---|
| pp.717–719，“Antihydrogen and ALPHA-g”“The effect of gravity”“The release experiment” | 磁能俘获、20 s 放出、读口及偏置定义；偏置是实验配置标签而非全阱均匀重力替代 |
| p.718 梯度段、p.719 未编号偏置公式 | 原文 $1.77\times10^{-3}\,T/m$、25.6 cm、$4.53\times10^{-4}\,T$；只核舍入单位运输 |
| p.721，“Classification of uncertainties”、Table 2/3、Conclusion | 0.75、0.13、0.16 的分类与发表结论；不制造新置信度 |
| PDF pp.10–11，“Simulations of the dynamics of trapped antihydrogen” | 磁场标定及三维轨迹、初态/离轴场模型依赖 |
| PDF p.12，“Analysis for escape curve and gravitational acceleration” | 校准效率/背景、Poisson 计数和模拟曲线到加速度的似然；不把裸二项计数当正式分析 |
| 同页 Data/Code availability | 原作者按合理请求提供，本件没有原始数据重分析 |

核验中只转录报告值并作单位乘法；不将 Table 3 各项当成已知独立 Gaussian 变量。原文的模拟模型误差不是本项目已重做的误差认证。

## S2：同一度规上的粒反粒关系

U. D. Jentschura, *Antimatter Gravity: Second Quantization and Lagrangian Formalism*, arXiv:**2003.08733v3**，版本记录 2020-07-18（PDF 首页另有 July 21, 2020）；Physics **2** (2020), 397–411。

- [版本摘要](https://arxiv.org/abs/2003.08733v3)
- [实际已读 PDF v3](https://arxiv.org/pdf/2003.08733v3)
- [期刊 DOI](https://doi.org/10.3390/physics2030022)

定位：§1 的 CPT/反地球区别；§3 式 (14)、(16)、(33) 的共同度规和电荷共轭。式 (33) 后正文也区分粒反粒关系与惯性/引力质量匹配。本件的弱场式 (3) 从明确采用的点粒子作用展开，不误标为该文 Eq.(35)；后者属于其它力参数讨论。

只采用标准电磁—度规支，不继承文中的更广统一模型、第五力解释、天体反中微子或全局正能量旁论。本件两模式 CAR 核验是自由场正规序恒等式的有限检查，不复现曲率场论证明。

## 本地复用与排重

- [P981](../../../archive_956_989/research_note_981.md)：采用的物质父作用，未由认知原则推出。
- [971](../../../archive_956_989/research_note_971.md)：复合材料内能与几何来源的既有接口。
- [1027](../../../archive_1009_1043/research_note_1027.md)、[1028](../../../archive_1009_1043/research_note_1028.md)、[1029](../../../archive_1009_1043/research_note_1029.md)：同一采用作用及来源的条件运输。
- [普通物质等效原理](../equivalence_principle_adoption.md)、[核子来源](../nucleon_source_adoption.md)：总质量与特定标量响应分开，不重复计实验。
- [钟初始化读口](../clock_initialization_readout_adoption.md)：被动落体与内部态分辨的钟任务不同，不从本件自动领取。

本次净覆盖为反氢实际磁放出记录；成熟采用计 0，不追加科学编号或累计数。作者文件为相邻 [正文](../antihydrogen_gravity_adoption.md)、本表、[check.py](check.py)、[results.json](results.json)。
