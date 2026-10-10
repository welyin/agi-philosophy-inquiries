# 光子视界接口：同一吸收系数、出射占据与能源账

2026-10-09。**成熟物理采用，科学新增0、认知公理新增0，不立1060。** 属P11，并连接P10的来源／反馈；[准入与历史去重](../_admission/after1059_horizon_scattering/selection.md)已保存。本项不重新推导温度，不计算径向数表。

## 1. 与旧1004的区别及共同对象

[1004](../../archive_990_1008/research_note_1004.md)已给4D标量HH／Unruh有限响应差、1+1校准和Hawking机制；其§6也已说明散射改变远处谱、净出射必须有能源。[1008的P11](../../archive_990_1008/1008/overall_operation_hypothesis_v2_2.md)保留这些范围。

这里把该机制接成**实际光子同一模式的入射吸收—量子出射字典**：不能为两个出口独立选“灰度”。采用[P981](../../archive_956_989/981/drafts/common_parent_contract_v1.md)低能电磁部门的领先自由Maxwell项、渐近平直Schwarzschild外域、坍缩后的晚时Unruh分支及其渐近时间归一。电荷／旋转、外部介质、来流、QED和引力高阶修正不在本领先公式中隐含消失；若保留它们须重新匹配。

用 $G=c=\hbar=k_B=1$，$M$ 为该单位下质量，$T_H=(8\pi M)^{-1}$。Maxwell辐射有 $\ell\ge1$、$m=-\ell,\ldots,\ell$ 和两个物理偏振 $p$；静电单极不是第三种或 $\ell=0$ 光子辐射模式。

## 2. 同一吸收系数连接两种任务

令 $\Gamma_{\omega\ell mp}$ 是纯向内视界边界下，来自无穷远的该经典模式被吸收的能量比例。Schwarzschild没有旋转超辐射，故

$$
0\le\Gamma_{\omega\ell mp}\le1,\qquad
|{\cal R}_{\omega\ell mp}|^2=1-\Gamma_{\omega\ell mp}.
\tag{1}
$$

这里用能流归一的反射振幅，不把任意径向函数振幅直接当概率。[Hawking1975，式(2.29)后及电磁场段落](https://astrofrelat.fcaglp.unlp.edu.ar/agujeros_negros/media/Papers/Hawking_1975-Particle_creation_by_black_holes.pdf)明确同一吸收份额进入晚时出射占据。[Page本人2004，§2式(5)](https://arxiv.org/html/hep-th/0409024v3)给带物种和偏振标签的写法。对无入射激发的参考分支，

$$
N_{\omega\ell mp}
=\frac{\Gamma_{\omega\ell mp}}{e^{8\pi M\omega}-1}.
\tag{2}
$$

**有限能量入射比较。** 固定一个角模／偏振，取归一振幅 $\int |h(\omega)|^2d\omega=1$，支撑在 $0<\omega_1\le\omega\le\omega_2<\infty$。在同一零均值Unruh参考上加一份弱的相干入射 $\alpha(\omega)=a h(\omega)$；这是带实际能源的另一次准备，非免费控制。在线性测试场合同内，新增平均应力是该经典波的应力，故增量分账为

$$
\begin{aligned}
\delta E_{\rm in}&=|a|^2\int \omega |h|^2\,d\omega,\\
\delta E_{\rm refl}&=|a|^2\int\omega|h|^2(1-\Gamma)\,d\omega,\\
\delta E_H&=|a|^2\int\omega|h|^2\Gamma\,d\omega,
\qquad
\delta E_{\rm in}=\delta E_{\rm refl}+\delta E_H .
\end{aligned}
\tag{3}
$$

式(3)是相对于未加探测波的同一背景的**增量**，不是把自发Hawking能流也算成反射光。其依据是自由场线性和参考零均值；不使用未知应力仅由迹距控制的错误替代。

**同一系数的量子出射。** 在晚时平稳测试场近似下，出射模按 $[b_\omega,b^\dagger_{\omega'}]=\delta(\omega-\omega')$ 归一，$b_h=\int h^*(\omega)b_\omega d\omega$，则未加探测波时

$$
\langle b_h^\dagger b_h\rangle_U
=\int |h(\omega)|^2
\frac{\Gamma_{\omega\ell mp}}{e^{8\pi M\omega}-1}\,d\omega .
\tag{4}
$$

式(3)—(4)是成熟模式关系在同一频带任务中的运输，不是本项目的新定理。一般不能把式(4)写成平均吸收率乘中心频率的Planck因子。它也不是1004局部探测器的响应概率：实际探测还需轨迹、开关、耦合与效率。

有限频带／有限能量并不同时意味着紧支撑时间访问。这里没有实现远处波包计数器，没有认证严格有限时间的全仪器后态或任意参考输入。$\Gamma$ 本身由Maxwell方程、同一几何和边界求得；本项目本次未重新数值求解它。

## 3. 光子功率、模式计数及有限范围反馈

在同一渐近时间 $u$ 的准平稳通量描述中，光子部门给

$$
P_\gamma(M)=\frac1{2\pi}
\sum_{\ell=1}^{\infty}\sum_{m=-\ell}^{\ell}\sum_{p=1}^{2}
\int_0^\infty \omega
\frac{\Gamma_{\omega\ell mp}}{e^{8\pi M\omega}-1}\,d\omega.
\tag{5}
$$

归一波包式(4)没有额外 $1/(2\pi)$；式(5)的该因子来自单位渐近时间的频率模式密度。[Page2004式(13)、(19)](https://arxiv.org/html/hep-th/0409024v3)报告原计算的

$$
P_\gamma=3.3638\times10^{-5}M^{-2}.
\tag{6}
$$

这个光子数值**已含两个偏振，不再乘2**；它是公开理论数值，本项目未复算灰体积分，也不是天体Hawking辐射实测。1976原刊的APS全文接口本次401，另一原刊扫描副本未提供可读文本；因此本项明确引用Page本人2004对原计算的报告，不冒称直接核过1976表格。不把其中历史假定的无质量中微子总和当作现实全部SM辐射率。

孤立出射时，式(5)进入总质量账的光子项 $dM/du|_\gamma=-P_\gamma$。若有式(3)的探测波，还须记其正的吸收增量；其它物种、来流及外部供给继续进入总账。**$P_\gamma$ 不是全部物种的总失重率。**

[Hawking1975§4，式(4.10)及其后](https://astrofrelat.fcaglp.unlp.edu.ar/agujeros_negros/media/Papers/Hawking_1975-Particle_creation_by_black_holes.pdf)以守恒应力联系远处通量与视界能源交换，并说明质量远大于Planck质量时采用缓慢变化的近静态背景。这里沿用这项领先半经典能源反馈：在所用窗口内需检查全部来源导致的相对质量变化足够小，且已越过所关心的坍缩瞬变。小变化本身不是严格总余项界。

远处功率不唯一确定全部局域 $\langle T_{ab}\rangle_{\rm ren}$；固定Schwarzschild测试场不是带净失重的精确时空解。本项不交付完整反作用解、任意有限窗误差或全部蒸发历史。

## 4. 本次接入与停止

新增的是可共同核对的物理字典：同一Maxwell传播与边界的 $\Gamma$ 同时控制入射能量吸收、晚时光子占据和相应光子质量账；若分别任意更换三个出口的系数，就违反这一采用合同。

认知原则、几何来源和初态选择没有因此被重新证明。未新增独立经验资料，未关闭P11或M4全部任务。到此停止径向扫描、Page小数精修和视界器件构造；完整蒸发、微观熵及UV问题不被提升为此有限适用描述的先决门槛。
