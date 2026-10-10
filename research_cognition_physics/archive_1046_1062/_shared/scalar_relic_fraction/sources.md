# 标量热遗迹共同字典：原始来源与采用边界

2026-10-09。配套[采用正文](../scalar_relic_fraction_adoption.md)。固定版本、成熟采用、计 0；未复跑任何似然或热史求解器。

## 一手来源定位

|来源|本次实际核查位置|用途与不领取的内容|
|---|---|---|
|Cline、Kainulainen、Scott、Weniger，[1306.4710v3](https://arxiv.org/html/1306.4710v3)，*Update on scalar singlet dark matter*|§I 式 (1)/(2)；§II 式 (3)；§III 式 (4)—(6)；§V；附录 A/B|势和门户归一、不可见宽度、湮灭与热平均、丰度分数重标度。旧实验边界与未来投影不当成新观测。|
|GAMBIT Collaboration，[1705.07931v3](https://arxiv.org/html/1705.07931v3)，*Status of the scalar singlet dark matter model*|§2.1 式 (1)/(2)；§2.2 式 (3)—(6) 与共振末段；§2.3 式 (7)—(11)；§3.1 表 2；§3.3—3.7；§4.2 表 4/5|所采用的固定版共同分析。数值来自表格／正文，不读图估值；不将其全参数图运输成 B993 子域。|
|Binder、Bringmann、Gustafsson、Hryczuk，[1706.07433v2](https://arxiv.org/html/1706.07433v2)，*Early kinetic decoupling of dark matter: when the standard way of calculating the thermal relic density fails*|§II.1—II.3；§III；§V|局部动能平衡失效时，单数密度方程不足；采用其所示量级限制，不计算新热遗迹。HTML 自动展示日期不作为原论文发表日期。|
|Hoferichter、Klos、Menéndez、Schwenk，[1708.02245](https://arxiv.org/html/1708.02245)，*Improved limits for Higgs-portal dark matter from LHC searches*|式 (9)—(15)，特别 §III 一／二体约定及式 (15) 标量分支|核靶与自由核子字典不能混用；衰变—散射比例沿用本项目旧采用，不计新结果。|

## 固定版资料与理论角色

GAMBIT §3.3 的 Planck 2015 总暗丰度为 $0.1188\pm0.0010$，在其主分析中只作上限；理论误差采用 5%，共振限制另见正文。§3.4/3.5 使用当时 Higgs 与直接探测似然；表 4 还含间接搜索和共同干扰参数，不能删去这些项后仍把表 5 称作重算的“三资料拟合”。这里没有提取新的置信区间。

同一核矩、晕密度及 SM 输入被共同使用；§2.3 采用截断 Maxwell 速度分布，固定速度参数并保留局部密度干扰参数，并未穷尽晕模型自由。本文式 (6) 明确把局部质量分数等于宇宙分数列为额外输入。欠丰度的另一个暗成分不属于 B993 自动已有的字段。

## 与旧合同逐项连接

- [B993 v1](../../../archive_990_1008/993/common_candidate_v1.md)：原 $\kappa/2$ 门户、$\mu_D^2>0$ 及不宣称现实身份／丰度的范围保持。
- [核子来源采用](../nucleon_source_adoption.md)：复用同一 $f_N$ 与零动量比例；[有限形状审计](../../_admission/nucleon_form_factors/public_data_audit.md)的统一误差停止线保持。
- [995](../../../archive_990_1008/research_note_995.md)、[996](../../../archive_990_1008/research_note_996.md)：额外相互作用或冷占据可改变碰撞项和初态，不能无证据归入纯门户热分支。
- [早晚膨胀采用](../early_late_expansion_audit.md)：本次旧 Planck 输入不与另一版本／共享天空资料作独立证据相乘。
- [共同范围校正](../joint_scope_after1062.md)：参数、历史和核子匹配在实际重叠处一致，不要求同一装置承载所有任务。

## 最小校核记录

2026-10-09 使用现有 Python 的 `decimal.Decimal`，不导入扫描／外部模型代码，实际计算打印点的 $m_D^2-\kappa\,246^2/2$：

$$
62.51^2-0.00065\,246^2/2=3887.83240,
\qquad
132.5^2-9.9\,246^2/2=-281997.95.
$$

单位为 $\mathrm{GeV}^2$。这仅核正裸质量子域，既不复核原文最佳拟合，也不认证热／高阶总误差。无需为两个初等代数恒等式另保存程序与结果副本。正文及本表的本地链接和最终 SHA 在交付时另作只读核验。
