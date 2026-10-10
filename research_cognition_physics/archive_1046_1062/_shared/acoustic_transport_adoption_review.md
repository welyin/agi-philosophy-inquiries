# 光子—重子声学采用：独立审阅

2026-10-09。审阅者 `/root/phase1044_goal_scope_review`。**结论：通过；成熟采用计0，无实质修订请求。** 未修改作者正文、程序、结果或历史冻结资产。

## 1. 对象与核验记录

核读[正式采用](acoustic_transport_adoption.md)、[准入](../_admission/acoustic_transport_after1055/selection.md)、[1005](../../archive_990_1008/research_note_1005.md)、[程序](acoustic_transport_adoption/check.py)和[结果](acoustic_transport_adoption/results.json)。直接核对下列四份一手原文的相关方程／节，没有图像检查、拟合或额外求解器。

使用现有Python、`-B -X utf8`运行作者程序默认只读模式，退出0；保存结果完全一致，13条本地链接与三份历史输入哈希通过。复核固定资产SHA256：

|资产|SHA256|
|---|---|
|`acoustic_transport_adoption.md`|`58ff0fac768c98fe415976d4c7307eb37f0ddd9df727d1eb0558e009e5d5f0e6`|
|`acoustic_transport_adoption/check.py`|`924c172e007eadc6a6151d3d19310ee165064131218dc9b8e5b9097134b69a47`|
|`acoustic_transport_adoption/results.json`|`5918c7bd906201238e4638b7d9e34e63eec0d1b30868f8f55f139de04b0f725e`|

## 2. 独立核对意见

### 2.1 正散射率、惯性与共同几何源

采用稿的 $\Gamma_T=a n_e\sigma_T>0$ 与 $R_b=3\rho_b/(4\rho_\gamma)$ 自洽；Ma–Bertschinger紧耦合节的 $R$ 确为其倒数。碰撞矩阵

$$
L=\Gamma_T
\begin{pmatrix}-1&1\\R_b^{-1}&-R_b^{-1}\end{pmatrix}
$$

有惯性加权左零向量、共同速度右零向量，另一特征值为 $-\Gamma_T(1+R_b^{-1})$。这验证两端相消与滑移衰减，不将碰撞抵消冒称整个含膨胀系统的普通总动量不变。

正文式(4)也正确：原文 $\delta T^0{}_0=-\delta\rho$ 给负的密度源，剪切约定给 $k^2(\phi-\psi)$ 的正 $12\pi G$ 系数；不能删掉其它组分应力后仍领取原GR关系。[Ma–Bertschinger，式(23)、(25)、(59)—(74)](https://arxiv.org/html/astro-ph/9506072)

### 2.2 偏振与扩散

从MB的偏振碰撞矩直接解领先代数，得到 $G_0=5F_2/4$、$G_2=F_2/4$，故有效四极阻尼 $f_2=3/4$。Hu–Sugiyama的 $4/(5f_2)$ 因而为 $16/15$；式(8)的两项、分母和振幅阻尼约定均相符。声速还依赖 $R_b$，不能由一项散射率独自决定。

原文推导采用紧耦合及相对慢变的声学展开；正文没有把该近似外推成最后散射附近全谱或统一误差定理。[Hu–Sugiyama，附录A式(A7)—(A14)，刊页555—556](https://background.uchicago.edu/~whu/Papers/small.pdf)

独立算术核对诊断点：$c_{\gamma b}^2=5/24$，滑移／剪切扩散分别为 $3/128000$、$1/9000$，总和 $31/230400$；两端碰撞加速度为 $-50,250/3$，误将重子一端加倍产生残差50。与结果一致。这些数字是无量纲代数检查，非宇宙数据。

### 2.3 电离史、材料温度与交换阶次

式(5)对应HyRec式(12)，用固有时而非式(2)的共形时；$x_e/(1+f_{\rm He}+x_e)$ 和冷却项正确。HyRec的实际联立对象还包括能级人口和辐射分布，不能只以Saha公式代替。领先弹性Thomson动量核不单独给完整热交换，正文已保留Compton、原子放能及反向能源来源的相应阶次与热库近似边界。[HyRec，式(12)、§V.4—5及§VI](https://arxiv.org/html/1011.3758v2)

### 2.4 拖曳与最后散射

光学深度向末时积分时，$g_\gamma=\Gamma_T e^{-\tau_\gamma}$ 为正；有限起点的未散射项未被遗忘。重子解的积分因子涉及 $a v_b$，所以正文式(10)的 $a(\Gamma_T/R_b)e^{-\tau_{\rm drag}}$ 及再归一化符合原文C8—C9。它不是单纯把光子权重中的 $\Gamma_T$ 换成 $\Gamma_T/R_b$。

采用稿区分 $r_s(\eta_*)$、$r_s(\eta_{\rm drag})$，保留可见分布宽度、引力驱动与初始项，未把一个时刻当成全传递核。[Hu–Sugiyama，附录C式(C8)—(C12)，刊页562](https://background.uchicago.edu/~whu/Papers/small.pdf)

### 2.5 公开数据的角色

Planck §3.1式(5)—(7)中的TT、TE、EE与原联合 $100\theta_*$ 四行转录正确，均按原68%范围；不是 $\theta_{\rm MC}$。§5.1使用拖曳声学尺，§7.7的复合史无显著偏离只在原模型与分析范围内成立。正文明确共享天空、lowE、参数校准及其它已采用Planck／BAO资料，没有新造独立Gaussian组合或经验显著性。[Planck2018 VI v4，§3.1、§5.1、§7.7](https://arxiv.org/html/1807.06209v4)

## 3. 历史净增加与整合范围

1005已经有成对四力与灰、各向同性、弹性示例。它不包含此次所需的Thomson角／偏振系数、电离史和不同可见权重。将这些成熟关系落到同一宇宙分支有实际解释增量，但不是新守恒定理或新的科学试验组；作者计0正确。

另核[1055后整合](integration_adoptions_after1055.md)，本次读取SHA为 `9bca614da07422445db0a840f3da5d5771f1e14bac3967d2250bb282cc1d1896`。其声学与[π⁰采用](pion_anomaly_adoption.md)总结未扩大范围：保留原QCD／QED分支、物理质量修正、2009预测与PrimEx-II差异；未将π⁰宽度扩成完整相位验证或异荷族预测；未把声学方程与既有BBN摘要相加成已重算的完整热史。计0、未全路线结项及经验共享的表述均一致。

本签审支持将两项成熟关系纳入共同输入／匹配账。它不认证全部B993暗态、任意未知量子输入与仪器后态、所有来源反馈或整个ROADMAP；也不要求先完成全探测器或无截断UV理论才能使用现有有限有效描述。到此停止，没有追加技术门槛。
