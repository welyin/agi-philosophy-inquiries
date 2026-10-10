# 双脉冲星共同平衡采用：来源与提取角色

2026-10-09。对应[采用正文](../double_pulsar_balance_adoption.md)，计0。固定下列版本；不是截至当前的全球最新合并分析。原文公式／表格通过文本核读，没有提图、重跑计时或重新拟合后验。

## 1. 2021实际观测与分析

M. Kramer et al., *Strong-Field Gravity Tests with the Double Pulsar*, **Physical Review X 11, 041050 (2021)**, published 13 December 2021.

- [最终期刊PDF](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.11.041050)
- [arXiv:2112.06795v1，2021-12-13](https://arxiv.org/abs/2112.06795v1)
- [同版可检索HTML](https://arxiv.org/html/2112.06795v1)

HTML重新生成时显示的日期不作为论文版本年；本件使用2021发表版本、以下节号和式号。

|原文位置|本件采用内容／角色|
|---|---|
|§V.1.1—2，式(9)—(17)|1PN、2PN、A的LT进动及同一 $I_A$；近自旋对齐和所忽略项的精度范围|
|§V.1.4，式(18)—(19)|Shapiro形状与质量函数的1PN修正，不只用Newton质量函数|
|§V.1.5，式(20)—(21)|自转质量损失与运动学修正；A的 $I_A$ 共用|
|§VI.1，表4及表后说明、脚注20—22|计时参数、末位误差、TDB约定；观测与派生质量分列。自转频率导数的红噪声误差讨论不能自动删除|
|§VI.2.1，式(35)—(38)与脚注23|采用 $k,s$ 及外部MoI求质量；共同Doppler尺度，质量首先是 $Gm$ 的提取|
|§VI.2.2，式(39)—(42)|Shklovskii、银河加速度及其相关性；不得用两个误差独立RSS代替作者联合运输|
|§VI.2.2，式(43)—(48)|质量损失修正、GW项、2.5PN＋3.5PN预测及联合MC比值；本件表值均从这里直接采用|
|式(48)后的两段|共同随机化观测参数与 $I_A$；宽MoI对照不是新资料|
|§VI.2.3|另一方向的GR条件MoI推断会使用 $\dot P_b$，不能反过来冒充未拟合辐射验证|
|§IV、§VIII.1及Appendix A|到达时刻、计时／VLBI天体测量与距离合同；本文不重做原管线|

这里“先提取、再检验”描述物理参数的使用角色，不表示各拟合参数统计独立。原文全部关键比值按其已发表分析采用；没有取得并重建完整原始到达时刻似然。公布中心数舍入后相加／相除不保证复现联合MC的最后一位，故不将中心比值抄算称为独立复算。

## 2. 同一物质的质量损失与辐射账

H. Hu, M. Kramer, N. Wex, D. J. Champion, M. S. Kehl, *Constraining the dense matter equation-of-state with radio pulsars*, **MNRAS 497, 3118–3130 (2020)**， [arXiv:2007.07725](https://arxiv.org/abs/2007.07725)，[已核HTML](https://arxiv.org/html/2007.07725)。

- §3、式(2)—(6)：进动、质量、转动惯量相连。
- §4、式(9)—(10)：同一 $\dot P_b$ 的内部／外部项。
- §4.1、式(11)—(13)：领先2.5PN及3.5PN辐射阻尼；本文只列领先显式式，最终数值沿2021原分析。
- §4.2、式(17)：自行对应的Shklovskii项。
- §4.3、式(18)—(19)：$\dot m c^2\simeq I\Omega\dot\Omega$ 与质量损失引起的正轨道周期贡献，和进动共用MoI。

这篇文献含预测未来精度的模拟；那些未来精度不是2021实测，更不是本件成果。其表2使用旧质量／距离作量级示例，本件不与2021表4混合取值。

领先辐射的原始理论来源亦见[P. C. Peters, Phys. Rev. 136, B1224 (1964)](https://doi.org/10.1103/PhysRev.136.B1224)。本次只核其原刊入口；显式计时约定及本件式(4)直接对照已读Hu等式(11)，不声称重读重证Peters全文。

## 3. 结构先验的实际资料依赖

Kramer原文参考101为T. Dietrich et al., *Multi-messenger constraints on the neutron-star equation of state and the Hubble constant*, **Science 370, 1450–1453 (2020)**，[arXiv:2002.11355v3](https://arxiv.org/abs/2002.11355v3)，修订于2020-12-21。

本次核该原始摘要的资料组成：GW170817及其电磁对应、GW190425、既有脉冲星X射线／射电测量和核理论。它不是仅由双脉冲星当前 $\dot P_b$ 得到的输入，但与本项目已有GW170817采用共享资料，不能称全部独立。半径分布到MoI还采用Kramer参考102的半径—MoI关系；本件按Kramer的已发表运输角色接受，未重建这一关系或EOS后验。

## 4. 项目链接及完成范围

- [1047](../../research_note_1047.md)：完整束缚源与辐射匹配，不能代替双星计时。
- [近期TT审计](../../_admission/after1059_dynamic_source/selection.md)：未准入的量子连续TT／完整记录桥保持原停止判断。
- [标准汽笛采用](../standard_siren_amplitude_adoption.md)：源与传播共同标定，未采用这里的轨道周期检验。
- [中子星潮汐采用](../neutron_star_tidal_adoption.md)：EOS及共享GW170817边界保持；不能独立重复累计显著性。

本件仅新增这一成熟共同映射的项目覆盖。没有计算新的系数、原始经验组或量子来源定理，也未降低ROADMAP的原验收要求。
