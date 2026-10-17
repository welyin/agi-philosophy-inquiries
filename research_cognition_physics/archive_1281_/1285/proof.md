# 1285：同一共形规范—标量有效作用的量子等待、完整来源与 FLRW 响应

2026-10-18。本轮是明确有效模型类别中的共同实例，新增认知公理为0。作者核验1组；状态以 [独审](independent_review.md) 为准。[报告](../research_note_1285.md)、[复算](check.py)、[结果](results.json)、[后继](NEXT.md)。

## 1. 对象与查重

[1031 §§3—5](../../archive_1009_1043/research_note_1031.md)已有同一经典规范父作用的非零来源、辐射型 FLRW 解与有限规范不变量预测；没有给本轮的正常量子态及 Born 概率。[1221 §7](../../archive_1219_1233/1221/proof.md)的32维压缩不保持原 Hamiltonian，不能直接充当真实传播。[625 §§3—8](../../archive_585_628/research_note_625.md)已有同过程、记录、来源联合压缩，但未给本例所需的具体截止率。

本轮直接量子化一个事前声明的齐次有效模型，不使用不变性未证的有限维压缩。它与3+1局部父作用的**经典约化**精确；把该约化量子模型作为有限有效类是公开建模选择，不是全4D量子场的精确不变扇区。未保留非齐次模式、完整场真空、重整化反常及高阶修正的匹配误差没有被赋0。

本例交付：同一正常量子初态、带相互作用的真实酉等待过程、由同一有效作用独立变分所得的全部均匀平均来源、推导后 Einstein 方程的解，以及有严格余项的非平凡有限量子概率。末端数学效果与实际读出装置资格分开；完整来源证书只覆盖未测等待过程。

## 2. 同一父作用、边界项与有效量子对象

取符号 $(-+++)$、$\hbar=c=1$。局部经典父为质量零 $U(1)$ 规范场与带荷复标量，写 $\phi=(\varphi_1+i\varphi_2)/\sqrt2$：

$$
S_m=\int\sqrt{-g}\,d^4x\left[-\frac14F_{\mu\nu}F^{\mu\nu}
-\frac12(D_\mu\varphi)\cdot(D^\mu\varphi)
-\frac{R}{12}|\varphi|^2-\frac\lambda4|\varphi|^4\right],\qquad
D_\mu\varphi=\partial_\mu\varphi+gA_\mu J\varphi,\quad
J=\begin{pmatrix}0&-1\\1&0\end{pmatrix}. \tag{1}
$$

即非最小项为 $-R|\phi|^2/6$。$g,\lambda$、共形耦合和字段类别均为建模输入，没有由本轮认知公理选出。空间取固定坐标周期、坐标体积 $V$；只保留齐次代表，使用其剩余时间依赖 $U(1)$ 作用。不把此有限模式类宣称为全部局域或大规范变换下的量子截断定理。

令 $ds^2=-N(t)^2dt^2+a(t)^2d\mathbf x^2$，$\chi=a\varphi$。直接展开曲率项及标量动能，得到

$$
S_m=S_{\rm red}-\left[\frac{V\dot a}{2N}|\chi|^2\right]_{t_i}^{t_f},\qquad
S_{\rm red}=V\int dt\left[\frac{a}{2N}\left(|\dot{\mathbf A}|^2+|D_t\chi|^2\right)
-\frac Na\left(\frac{g^2}{2}|\mathbf A|^2|\chi|^2+\frac\lambda4|\chi|^4\right)\right]. \tag{2}
$$

边界项不能悄悄删除后又在不同固定变量下变分。以下使用式(2)对应的固定端点/正则边界约定；在原场变量下的体变分与在 $\chi$ 下的体变分相差场方程项及已显示的边界项。它们在同一物质解上给相同均匀改进 Hilbert 来源。

固定 $V=1$。令 $q=(A_1,A_2,A_3,\chi_1,\chi_2)\in\mathbb R^5$，有效量子母体为 $L^2(\mathbb R^5,dq)$ 中的剩余 Gauss 不变子空间，

$$
H_m=\frac Na H_0-gA_0G,\qquad
H_0=-\frac12\Delta_q+\frac{g^2}{2}|\mathbf A|^2|\chi|^2+\frac\lambda4|\chi|^4,\qquad
G=\chi_1p_{\chi_2}-\chi_2p_{\chi_1}. \tag{3}
$$

$D_t=\partial_t+gA_0J$ 与上述 $G$ 约定共同给出约束项 $-gA_0G$；下文取时间规范 $A_0=0$，且始终保持 $G\psi=0$，不省略规范耦合归一。非负闭二次形式定义 $H_0$ 的 Friedrichs 自伴实现；本轮 Gaussian 及其每个固定次 $H_0$ 作用都是多项式乘 Gaussian，均在相应算符域中。$H_0$ 与 $U(1)$ 及对三个 $A_i$ 的 $SO(3)$ 作用对易。这里后者保证均匀应力期望各向同性，不要求整个周期空间存在任意连续旋转的全局几何对称。

事前固定

$$
g=1,\quad\lambda=\frac14,\quad
\psi_0(q)=(2\pi)^{-5/4}e^{-|q|^2/4},\qquad
\psi_\eta=e^{-i\eta H_0}\psi_0,\quad d\eta=\frac Na dt. \tag{4}
$$

这是正常 Gauss 态；$|\psi_0|^2$ 的五个实坐标方差均为1。有限矩及不变对称性由同一酉传播保持。没有把分布态、固定经典规范坐标或热平衡态替代此准备。

## 3. 独立来源变分与同一几何

量子有效作用使用同一正则字典，物质部分为

$$
S_{m,\rm eff}=\int dt\left[\frac i2\langle\psi\mid\overleftrightarrow{\partial_t}\psi\rangle
-\frac Na\langle H_0\rangle+gA_0\langle G\rangle\right]. \tag{5}
$$

规范约束为 $G\psi=0$；量子方程为式(4)。这定义了本次有限模式量子有效来源，**不等于声称全量子场的重整化有效作用没有反常**。先对 $N,a$ 变分、再使用物质方程，可得

$$
\rho=-\frac1{Va^3}\frac{\delta S_{m,\rm eff}}{\delta N}
=\frac{E}{Va^4},\qquad
P=\frac1{3NV a^2}\frac{\delta S_{m,\rm eff}}{\delta a}
=\frac{E}{3Va^4},\qquad E=\langle\psi_0,H_0\psi_0\rangle. \tag{6}
$$

归一化时间积分密度按式(5)理解。旋转及 Gauss 不变性使平均动量流为0、空间应力各向同性。$E$ 沿真实量子演化守恒，因此完整均匀 Ward 恒等式为

$$
\frac{d\rho}{d\tau}+3\frac1a\frac{da}{d\tau}(\rho+P)=0,\qquad d\tau=Ndt. \tag{7}
$$

这里只签均匀有效域的完整来源及守恒；不由式(7)反推任意非齐次局域 Ward 定理。

在已完成的有效 Einstein 理论内验证实例，取固定 $M_P^2=11/2$、$\Lambda_{\rm eff}=0$、$a(0)=1$。完整齐次引力作用为 $-3M_P^2V\int a\dot a^2/N\,dt$（相应端点项已计）；其独立 lapse 约束与 $a$ 方程给

$$
3M_P^2\left(\frac{a'}{a^2}\right)^2=\rho,
\qquad a''=0,
\qquad a(\eta)=1+\frac\eta2. \tag{8}
$$

下面的精确矩给 $E=33/8$，所以式(8)满足两个独立 FLRW 方程。没有用 Einstein 张量定义来源。固定同一 $g,\lambda,V,M_P$，从式(3)、(5)、(8)得到量子过程、来源与几何，未并列拼接1221的量子模型和1031的另一几何。

全部原能源包括

$$
\langle T_A\rangle=\frac38,\quad
\langle T_\chi\rangle=\frac14,\quad
\langle V_{A\chi}\rangle=3,\quad
\langle V_\chi\rangle=\frac12,\quad E=\frac{33}{8}. \tag{9}
$$

这些是初片分账；各分项不被假定分别守恒。没有减掉量子动能基线或规范—标量相互作用存能。没有调用读出、重置或时变外部控制，因此等待过程中不另欠一个免费开关能源。把 $\eta$ 换成固有时仅是同一度规的规范选择，不是已构造了材料钟仪器。

## 4. 精确矩与严格有限量子预测

置 $X=|\mathbf A|^2,Y=|\chi|^2$。原准备下二者独立，$X\sim\chi^2_3,Y\sim\chi^2_2$。对任意多项式 $P(X,Y)$，直接微分而非截断传播给

$$
\begin{aligned}
\frac{H_0(P\psi_0)}{\psi_0}
={}&-2X P_{XX}+(X-3)P_X-2Y P_{YY}+(Y-2)P_Y\\
&+\left[\frac54-\frac{X+Y}{8}+\frac{XY}{2}+\frac{Y^2}{16}\right]P.\end{aligned}\tag{10}
$$

由 $P_0=1,P_{j+1}=\mathcal DP_j$ 和
$\mathbb E[X^iY^j]=\prod_{k=0}^{i-1}(3+2k)\prod_{k=0}^{j-1}(2+2k)$，独立精确有理运算得到

$$
(m_1,m_2,m_3,m_4)=
\left(\frac{33}{8},\frac{2683}{64},\frac{495575}{512},\frac{161114457}{4096}\right),\quad
m_j=\langle H_0^j\rangle,\qquad
\operatorname{Var}(H_0)=\frac{797}{32}. \tag{11}
$$

$P_0^{\rm surv}=|\psi_0\rangle\langle\psi_0|$ 是有界、Gauss 及旋转不变数学效果。其同一真实过程概率为
$p(\eta)=|\langle\psi_0,e^{-i\eta H_0}\psi_0\rangle|^2$。以 $\mu$ 为此准备的 $H_0$ 谱测度，则

$$
1-p(\eta)=\iint[1-\cos(\eta(x-y))]\,\mu(dx)\mu(dy).
$$

对所有实数 $z$ 有 $z^2/2-z^4/24\le1-\cos z\le z^2/2$，所以这是无需谱数值、无时间幂级数收敛假设的全实参数界：

$$
\operatorname{Var}(H_0)\eta^2-\frac{D_4\eta^4}{24}
\le1-p(\eta)\le\operatorname{Var}(H_0)\eta^2,
\quad D_4=2m_4-8m_3m_1+6m_2^2=\frac{14661753}{256}. \tag{12}
$$

事前固定 $T=1/100$，在同一窗口同时得到

$$
\boxed{\frac{505192749}{204800000000}\le1-p(T)\le\frac{797}{320000}},\qquad
\boxed{\frac3{1250}<1-p(T)<\frac1{400}},\qquad
\boxed{a(T)=\frac{201}{200},\quad\tau(T)=\int_0^Ta\,d\eta=\frac{401}{40000}}. \tag{13}
$$

两份独立复制谱测度只是证明式(12)的积分工具，不是新增物理副本或实施两次测量。几何 $a$ 是本例独立度规/尺度量；其数值不是用 survival 概率重新标定得到的。

## 5. 后态、伙伴与实际权限边界

对任意旁挂 $R$，等待通道是 $\operatorname{Ad}(e^{-i\eta H_0})\otimes\mathrm{id}_R$。数学终末二结果仪器可写

$$
\mathcal I_0(\rho)=P_0^{\rm surv}\rho P_0^{\rm surv},\qquad
\mathcal I_1(\rho)=(I-P_0^{\rm surv})\rho(I-P_0^{\rm surv}), \tag{14}
$$

保留两结果及其全部未归一后态，并张量 $\mathrm{id}_R$；没有后选。式(14)的 CPTP 完整性不自动授予实际空间访问、制备或读出权限。主准备为纯 $\psi_0$，其任何 $R$ 扩展必为乘积；不能把本例几何误差推广成任意未知输入/纠缠初态的统一半经典背景，也不能给测后分支瞬时分配一套 Einstein 几何。

数学仪器的非选择后态能源仍有限：令 $Q=I-P_0^{\rm surv}$，则
$\|H_0^{1/2}Q\psi_T\|\le\sqrt E+|\langle\psi_0,\psi_T\rangle|\sqrt E\le2\sqrt E$，故

$$
0\le \operatorname{Tr}H_0[\mathcal I_0+\mathcal I_1](|\psi_T\rangle\langle\psi_T|)\le5E=\frac{165}{8}. \tag{15}
$$

它一般不等于 $E$。若实施此读出，器件/控制/记录必须承担相应能量和压力，式(6)—(8)不能未经追加源合同穿过测量窗。有限成本和正常后态并不等于装置已被构造。本轮不以另开一串仪器工程为已证数学共同实例的前置条件。

## 6. 误差承接与真正未证项

本模型内部式(6)—(13)为精确或显式有限余项结论。若以后把它用作另一系统的有效描述，须分别给同过程概率误差和完整来源匹配误差，不能由小迹距离推出无界能源接近。[07的场桥](../../archive_1159_1160/_shared/field_bridge_after_1160.md)已有此区分；[10论文](../../10_认知操作公理与有效广义相对论_阶段论文.md)允许公开宏观有效仪器输入，但其旧数值预算不自动迁入本例。

例如若匹配后的均匀方程已独立认证为
$\widetilde a'=\sqrt{\widetilde E/(3M_P^2V)}+r_a$、
$|\widetilde E-E|\le\epsilon_E<E$，则直接稳定性界为

$$
|\widetilde a(T)-a(T)|\le|\widetilde a(0)-a(0)|+
\frac{T\epsilon_E}{\sqrt{3M_P^2V}\,[\sqrt{E-\epsilon_E}+\sqrt E]}
+\int_0^T|r_a|d\eta. \tag{16}
$$

这里的 $\epsilon_E,r_a$ 是待认证的独立输入，未设为0而声称全4D匹配完成。概率预测若有独立总误差 $\epsilon_p$，仅将式(13)两端各放宽 $\epsilon_p$。本轮不选择规范群、共形耦合、$g,\lambda,M_P$ 或初态；不签热参考、实际温标、全部局域量子场、全空间约束或量子引力。净新增是一个同一声明有效对象内的**正常量子过程—同作用完整均匀源—几何—有限概率**接线，不是重新推导已验收的 Einstein 方程。
