# 1283 证明：有限相位校准反向限制规范链路动能

2026-10-17。审阅状态以 [independent_review.md](independent_review.md) 为准。[报告](../research_note_1283.md)、[复算](check.py)、[结果](results.json)、[后继](NEXT.md)。本轮没有新增公理；一组作者诊断支持下面的解析命题。

## 1. 查重、对象与新增条件

[557 §§2、5—6](../../archive_554_584/research_note_557.md)已经给同一 $SU(2)$ 链路、正常物质源、规范不变边势与正向测量注能；[566 §7](../../archive_554_584/research_note_566.md)、[585 §5](../../archive_585_628/research_note_585.md)、[586 §3](../../archive_585_628/research_note_586.md)也已有给定二阶动能的相位或仪器能量恒等式。本轮不重报这些一般共轭公式。新增是：在 [1282](../1282/proof.md) 同一正常 Gauss 态上，把有限幅度相位调整的源能量响应作为**额外待验证校准条件**，反向限制一个公开扩大的链路动能类别。

沿用 $\mathscr H=L^2(\mathbb R^{10}\times SU(2))$、归一 Haar 测度、原两端局部作用及物理 Gauss 部门。记 $\mathcal L=-\Delta_G$，其谱为 $j(j+1)$，基本四元数生成元归一为 $J_a^2=-I/4$。设

$$
H_{b,c}=H_{\rm mat}+F+b\mathcal L+c\mathcal L^2,
\qquad F=\frac\kappa2|X_1-R(U)X_2|^2,
\qquad b>0,\ c\ge0,\ \kappa>0. \tag{1}
$$

$H_{\rm mat}=-a_{\rm kin}\sum_i(\Delta_{X_i}+\partial_{s_i}^2)+W$，$a_{\rm kin}>0$；$W=W_1+W_2$ 是557原正四次势，系数保留，不因本轮校准重新拟合。这里 $\kappa$ 已包括557的体积因子 $v$。式(1)由正闭形式的 Friedrichs 扩张定义；$c\|\mathcal L\psi\|^2$ 是新增项，不能假称557早已采用了它。各项保持物理 Gauss 部门。

新增操作合同是同一源上允许 $U_\theta=e^{i\theta F}$，其中 $\theta F$ 无量纲。它是保 Gauss 的有效相位控制原语，并非已由旧指针自动制造的自治控制器。每个比较从同一准备开始；$0,\theta,2\theta$ 是三个任务的幅度，不能把两次源增能直接相加。新成本条件只针对这份**链路调整菜单**，不推广到所有物理相位操作或相对论粒子能源。

## 2. 同一有限正常准备和完整源能源

令

$$
y=\frac{X_1+R(U)X_2}{\sqrt2},\quad z=\frac{X_1-R(U)X_2}{\sqrt2},\quad
\psi_{s,t}=\frac{e^{-|y|^2/(4s)-|z|^2/(4t)}}{(2\pi s)(2\pi t)},\quad s,t>0. \tag{2}
$$

再乘两个概率方差为1的实归一 singlet Gaussian。$|\psi|^2$ 下 $y,z$ 分别为四维独立中心 Gaussian，协方差 $sI,tI$。这是1282原正常物质—链路相关态，不是固定链路或零测支撑态。它及 $U_\theta\psi$ 的每个固定阶微分都是多项式乘 Gaussian，因而属于本轮所有固定多项式能源所需的形式域，且能源有限。

为避免与 $c$ 混淆，记

$$
\zeta=X_1\cdot R(U)X_2=\frac{|y|^2-|z|^2}{2},\quad
D=|y|^2|z|^2-(y\cdot z)^2,\quad q=\frac{s-t}{4st}.
$$

基本表示恒等式给 $\Delta_G\zeta=-3\zeta/4$、$|\nabla_G\zeta|^2=D/4$。直接 Gaussian 积分得到

$$
\begin{aligned}
\langle D\rangle&=12st,&\langle D^2\rangle&=360s^2t^2,\\
\langle\zeta^2\rangle&=6(s^2+t^2)-8st,&\langle\zeta D\rangle&=36st(s-t).
\end{aligned}\tag{3}
$$

例如固定 $y$ 时，$D=|y|^2|z_\perp|^2$，故 $\langle D^2\rangle=(24s^2)(15t^2)$。同理所有正整数 $n$ 有 $\langle D^n\rangle=(4st)^n(2)_n(3/2)_n>0$。

不遗漏源的原基线：

$$
E_0=a_{\rm kin}(s^{-1}+t^{-1}+1/2)+\langle W\rangle+4\kappa t
+\frac{3b(s-t)^2}{16st}
+c\left[\frac{9(s-t)^2}{64st}+\frac{45(s-t)^4}{512s^2t^2}\right].\tag{4}
$$

原 $W_i=\frac v4(|X_i|^2-u_h,s_i^2-u_s)^T L_{\rm pot}(|X_i|^2-u_h,s_i^2-u_s)$ 的均值为，记 $h=s+t$，

$$
\langle W\rangle=\frac v2\left[L_{11}(6h^2-4u_hh+u_h^2)
+2L_{12}(2h-u_h)(1-u_s)+L_{22}(3-2u_s+u_s^2)\right].\tag{5}
$$

它有限且非负，后续相位保持其值。代码保留这一公共项的解析身份，不用任意新势替换旧势。旁挂参考的操作为 $U_\theta\otimes I_R$；能量数值仅对式(2)准备作证，不扩大为所有未知输入的相同响应。

## 3. 精确四次响应及有限判据

定义 $G_F=|\nabla_GF|^2=\kappa^2D/4$、$\mathcal D_F=2\nabla_GF\cdot\nabla_G+\Delta_GF$。在共同光滑域上

$$
U_\theta^\dagger\mathcal L U_\theta
=\mathcal L-i\theta\mathcal D_F+\theta^2G_F,\qquad
\frac{\mathcal L\psi}{\psi}=\ell=\frac{3q\zeta}{4}-\frac{q^2D}{4},\quad
\frac{\mathcal D_F\psi}{\psi}=d_F=\kappa\left(\frac{3\zeta}{4}-\frac{qD}{2}\right).\tag{6}
$$

实准备使实部和虚部正交，因此

$$
\|\mathcal L U_\theta\psi\|^2
=\left\langle(\ell+\theta^2G_F)^2+\theta^2d_F^2\right\rangle.
$$

坐标动能的增量为 $16a_{\rm kin}\kappa^2t\theta^2$，二阶群动能增量为 $3b\kappa^2st\theta^2$，乘法势不变。代入式(3)，得到完整源能源的精确身份

$$
E(\theta)-E_0=A_2\theta^2+cC_4\theta^4,
\quad C_4=\frac{45}{2}\kappa^4s^2t^2>0,
\quad A_2=16a_{\rm kin}\kappa^2t+3b\kappa^2st
+\frac{9c\kappa^2}{16}[5(s^2+t^2)-6st].\tag{7}
$$

于是任何固定非零 $\theta$ 都给

$$
\mathcal R_\theta:=E(2\theta)-4E(\theta)+3E_0=12cC_4\theta^4.\tag{8}
$$

在式(1)类别内，精确倍幅残差为零当且仅当 $c=0$。只要求 $\theta\to0$ 时领先二次 Taylor 响应不能排除 $c>0$：整个候选族已经满足 $A_2\theta^2+O(\theta^4)$。这是一项必要的删弱区分。

设三份源能源报告满足独立登记的绝对误差 $\eta_0,\eta_1,\eta_2$，包含各自准备、切换、测读和统计误差，则

$$
|\widehat{\mathcal R}_\theta-\mathcal R_\theta|\le
\eta_2+4\eta_1+3\eta_0=:\eta_R,
\qquad
0\le c\le\frac{|\widehat{\mathcal R}_\theta|+\eta_R}{12C_4\theta^4}.\tag{9}
$$

参数及幅度须独立校准；若有参数误差，应使用分母的已证正下界。本轮没有制造这些测读合同，也不把程序舍入误差当实验误差。有限精度不推出 $c=0$。

## 4. 一次覆盖有限多项式类别

进一步允许实有限多项式 $p(\mathcal L)=\sum_{j=0}^n a_j\mathcal L^j$，在 Casimir 谱上有下界；非恒定时最高非零系数 $a_n>0$。可加常数使谱下界非负，常数在响应中抵消。其余源和式(2)保持。

式(6)的 $n$ 次有序乘积是 $\theta$ 的有限算符多项式，最高次数为 $2n$，唯一最高项为 $\theta^{2n}G_F^n$。实系数、实准备及复共轭把 $\theta$ 送为 $-\theta$，所以能源均值是偶多项式。对 $n\ge2$，其最高系数

$$
a_n\langle G_F^n\rangle
=a_n\left(\frac{\kappa^2}{4}\right)^n(4st)^n(2)_n(3/2)_n>0.\tag{10}
$$

故在任意非空有限开放幅度窗上，**精确**二次响应恒等式已经迫使 $n\le1$；这是多项式恒等定理，不需要执行无限大幅度或去掉截止。它不选择 $b$、$W$、规范群、图或其它许可动力学。有限数据若无次数及系数控制，不能把式(9)推广为任意多项式或非局域算符的统一逼近结论。

## 5. 固定诊断菜单、控制账与结论

在首次计算前固定以下代数诊断参数，不扫描参数：$a_{\rm kin}=1/2,b=3/10,\kappa=4/5,s=2,t=1,\theta=1/2$，比较 $c=0$ 与 $c=1/10$。旧 $W$ 只作同一公共有限项，不改其实际系数。于是

$$
C_4=\frac{4608}{125},\quad12C_4\theta^4=\frac{3456}{125},\quad
\mathcal R_{1/2}\big|_{c=1/10}=\frac{1728}{625}.
$$

若每份报告误差为 $10^{-4}$ 且观察残差为零，式(9)给 $c\le1/34560$。这是条件误差例，不是本轮取得的测量数据。

代码以八个 Gaussian 坐标的精确有理矩展开核式(3)、(6)—(9)，并用固定群路径的二阶差分独立诊断群导数归一。只记一组作者核验；一般多项式结论由解析证明承担，不逐阶另算轮次。

$\Delta E$ 是源的增能。闭合实现必须把相反能量变化与控制、切换、交互、指针及记录的余额共同登记；本轮没有证明控制资源为零、最小功等于此增能或实现了自治新设备。$U_{2\theta}=U_\theta^2$ 本身不蕴涵式(8)为零，完整能量守恒也不蕴涵新成本角色。该角色若要采用，仍须认知含义、删弱及共同实现审查。候选的排除不构成全部认知公理或强生成猜想的反证。
