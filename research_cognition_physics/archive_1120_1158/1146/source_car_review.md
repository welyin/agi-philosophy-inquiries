# 1146 独立来源与 CAR 资格审查

2026-10-12。审阅者：independent_record_parent_screen。绑定[事前准入](../_shared/1120_1156/round1146_admission.md)，复用[交互资格](../_shared/1120_1156/after1145_interacting_ward_screen.md)、1128—1129 的完整 CAR/正常来源工具。此文件只验来源及未知码的接合资格；作者的移动包误差和主线完整证明冻结后另作最终独审。

## 1. 固定范围及意见

采用 $T=1$、$\delta=10^{-3}$、$m=10^{10}$、$p_0=3m e_z/4$、$\epsilon=10^{-16}$、$A=10^{100}$；交互系数沿 1129 的 $\lambda=\pi/2$，$F$、$\eta$、$\rho$ 与来源准备规则沿旧定义。材料改为两异质量味，旋量耦合为 $\beta\sigma_x$。实际读口是声明有限区域内两单位 Cauchy 模式的完整味矩阵代数，不领取原 288 维 $D/M_r$ 的身份。

**资格通过，但完整轮次仍待波包及读数反差证明。** 以下结论保实际量子来源、完整粒子/空穴、全 CAR 余部及任意未知 Q/R。没有采用期望源代替实际源，也没有把局部 Cauchy 模式先截成纯正频。

本材料须明确选取新自由真空 $S=\operatorname{diag}(S_m,S_{m+\delta})$ 的正常表示；不能声称已经与旧 $m=1$ 的无限空间真空建立酉等价。旧工具按相同结构迁移，新背景和新正常准备属于本轮明示材料合同。

## 2. $\beta$ 核保留有界 H 与完整来源

以真实 Cauchy 两谱

$$
h=\boldsymbol\alpha\cdot p+\beta(mI+\delta P_1)
$$

定义正能源粒子/空穴 Fock Hamiltonian $H_F$。$J^\beta_\epsilon(x)$ 是完整有限 CAR 双线性，内部矩阵 $\beta\sigma_x$ Hermitian、范数 1。两味四旋量的精确粗界为 $4q_\epsilon$；本轮保守使用旧 $8q_\epsilon$：

$$
q_\epsilon\le8\epsilon^{-3},\quad
\|J^\beta_\epsilon(x)\|\le8q_\epsilon,\quad
C_*:=8q_\epsilon\|\eta\|_1\le33536\epsilon^{-3}
=3.3536\times10^{52}.
$$

因此

$$
H=H_F+H_\phi+
\lambda\int\eta(x)F(\phi_\epsilon(x)/A)J^\beta_\epsilon(x)\,d^3x
$$

是 $D(H_F+H_\phi)$ 上的自伴半有界算符，且 $\|V\|\le16\lambda C_*<9\times10^{53}$。此结论在完整正常 Fock 上成立，不靠一粒子矩阵或点态正规序密度的未经证明定义。

1129 的正常源比较仅使用 $F$ 的 Lipschitz 常数 $L\le512$、上述 CAR 范数和自由来源涨落 $v_\epsilon<3/\epsilon$。异质量与 $\beta$ 不改变该部分，因此对全部材料输入及任意参考，固定独立来源准备时，

$$
\epsilon_\phi(t)
\le\lambda C_*L\,\frac{3t}{\epsilon A},
\qquad
\epsilon_\phi(1)<8.1\times10^{-29}.
$$

比较过程的材料演化必须是新的
$h+\lambda F(t)K_\epsilon\eta K_\epsilon\beta\sigma_x$；
不可以沿用等质量矢量源的精确旋转结论。

来源态在实际 H 下仍为量子态并允许与材料纠缠。完整等距界比较其真实后态与明确参照，不是先迹掉来源后补一份经典能源账。

## 3. 真实正常码，不能再用 $1/4$ 空概率

令 $f=e^{ip_0\cdot x}b(x)\xi$，$b$ 为旧归一紧支实包络，$G=\|\nabla b\|^2<400$。$\xi$ 为 $h_m(p_0)$ 的正能旋量，$E_0=5m/4$；$f_i$ 为该模式的两个味标签。

异质量真空跨味协方差为零，$s_i=\langle f_i,Sf_i\rangle$ 一般不同。对 $m_i=m+\delta_i$，$\delta_0=0,\delta_1=\delta$，

$$
s_i
\le\frac{G+\delta_i^2}{(m_i+E_0)^2}
<\frac8{10^{19}}.
$$

证明由
$P_i^-(p)(h_i(p)-E_0)=-(\omega_i+E_0)P_i^-(p)$
作用于 $\xi$，再用 Clifford 反对易及 Plancherel 得到。该估计保负谱小尾，而非删除它。$f$ 的紧支性与真实负谱部分须同时保留。

置 $e_0=(1-n_0)(1-n_1)$，则

$$
Z=\|e_0\Omega\|^2=(1-s_0)(1-s_1)>1-2\times10^{-18}>0,
\quad
\Omega_0=e_0\Omega/\sqrt Z,
\quad
J_i=c^\dagger(f_i)\Omega_0 .
$$

直接 CAR 身份给 $c(f_i)\Omega_0=0$、$J^\dagger J=I_2$。这是精确正常输入码；未知逻辑态和 R 可任意纠缠，余部不依赖该未知态。准备权限仍明示，不将 $Z$ 大解释为已经完成自动重置。

若以后改成经过混味的初始模式，一般应使用 $\det(I-C)$，不能无条件继续两个对角因子相乘。本轮未作该改变。

## 4. 新输入能源与实际后态域

令 $a_i=\|h_i f\|$。由相同波包的中心谱式，

$$
a_i\le E_0+\sqrt{G+\delta_i^2}<5m/4+21 .
$$

正常真空上 $H_F\Omega=0$，有限 CAR 平移导子给

$$
\|H_F\Omega_0\|
\le\frac{2(a_0+a_1)}{\sqrt Z},\qquad
\|H_FJ_i\|
\le a_i+\frac{2(a_0+a_1)}{\sqrt Z}.
$$

使用 $1/\sqrt Z<1+2\times10^{-18}$ 及两列 Cauchy–Schwarz，

$$
\sup_{\|z\|=1}\|H_FJz\|<9\times10^{10}.
$$

这是新材料的实际 Fock 图范数界，不是旧输入数 268。参考扩张不改变该算符范数。

正常 KG 源仍满足 $\langle H_\phi\rangle\le1073A^2=1.073\times10^{203}$，且属于 $D(H_\phi)$；光滑紧支源准备保证其加权二阶矩有限。输入因而在 $D(H_F+H_\phi)=D(H)$，自治传播保该域，真实输出材料和来源各自均能有限。独立能源账为

$$
\langle H_F+H_\phi\rangle_t
\le 9\times10^{10}+1073A^2+2\|V\|,
\qquad
\Delta\langle H_F\rangle+
\Delta\langle H_\phi\rangle+
\Delta\langle V\rangle=0 .
$$

这里没有以小态距控制无界能源；巨大但有限的来源资源必须公开，不能作为效率或可制造性结论。

## 5. $\beta$ 核下的实际源反作用

源 CCR 与材料因子仍分立，1129 的测试场有界 Weyl 微分证明直接适用。只是把原 $J_\epsilon$ 全部一致替为 $J^\beta_\epsilon$，实际 retarded 项保持：

$$
\left\|
\frac{\phi_{\epsilon,H}(t,x)-\Phi_0(t,x)}A
\right\|
\le \frac{\lambda L C_*q_\epsilon t^2}{2A^2},
\qquad
b_{\rm back}(1)<1.1\times10^{-96}.
$$

同 H 的测试动量力为

$$
\frac d{dt}\Pi_H(g)
=\phi_H(\Delta g)
-\frac\lambda A\int\eta(x)(\rho_\epsilon*g)(x)
\big[F'(\phi_\epsilon(x)/A)J^\beta_\epsilon(x)\big]_H\,d^3x .
$$

来源真实后态、反作用算符与关联均保留。某些输入的平均力为零，不意味着交互被取消。

因 $F(0)=F(1)=0$，端点实际交互均能独立满足

$$
|\langle V(0)\rangle|,\ |\langle V(1)\rangle|
\le\lambda C_*L
\left[\frac3{\epsilon A}+b_{\rm back}(1)\right]
<8.1\times10^{-29}.
$$

端点小交互能不代表中途无能源交换，也不让静态部署成为免费来源。共同空间平移仍给固定 $\eta$ 的外力责任式；动态支撑未实现，不能宣布全装置四 P 闭合。

## 6. 两异质量目标列的完整 CAR 校正

若作者证明同一经典比较传播的两正交实际列 $v_i$ 接近两正交目标列 $w_i$，$\|v_i-w_i\|\le d_i$，可以逐列作完整 CAR 校正：

1. 第一项至多两模式 SU(2) 把 $v_0$ 送到 $w_0$，完整 CAR 酉距恒等不超过 $d_0$。
2. 经第一步后第二列仍正交于 $w_0$，距 $w_1$ 不超过 $d_0+d_1$。
3. 在 $w_0^\perp$ 内作第二次 SU(2)，固定第一目标列。

因此

$$
\|\Gamma(R)-I\|\le2d_0+d_1\le3d\quad(d_i\le d).
$$

该界覆盖全部占据与正常余部，不是仅一粒子范数；复杂重叠和纯相位也由保行列式的两模式实施处理。原独立 rail 的 $2d$ 不应不加说明地照搬。

保实际 CAR 实施相位，置
$\Xi=\Gamma(R)\widehat U\Omega_0$，
则 $c(w_i)\Xi=0$ 且
$\Gamma(R)\widehat UJ_i=c^\dagger(w_i)\Xi$。
所以完整未知输入及任意 R 的联合等距误差可取

$$
\epsilon_{\rm full}
\le\epsilon_\phi+2d_0+d_1 .
$$

$\Xi$ 一般不是旧真空。固定的合法 Cauchy 模式读数或仪器可收缩此误差；来源/材料的真实后态没有被删除。目标列、源相位、初始资料和最后读口必须全部采用同一字典。

## 7. 因果范围及正式验收前剩余检查

新的自由 Dirac 质量不改变主部有限传播速度，$\beta\sigma_x$ 不扩大有限空间核的跳跃半径。若仍保 $K_\epsilon\eta K_\epsilon$、$\lambda\int_0^tF\le\lambda$ 及 1129 原材料偶干预/占据菜单，则其 Dyson 尾方法可以迁移；源干预、完整无限局部代数或逐点材料微因果不能自动加入。

还须由完整证明核验：

- 作者在新完整两谱上的移动包界、两列正交性及内部 $\delta P_1/\sigma_x$ 演化；
- 同一个实际区域读口的有限差超过全部误差；错误类型对照保持相同准备和来源；
- 因果菜单是否实际采用，全部使用中的常数与时限是否一致；
- 新有限 Cauchy 读口的材料权限，不能改称已构造原 $D/M_r$；
- 不以局部成功领取共同 boosts 全域、自治支撑、总四 P 或整个 1120 阶段。

本资格的固定参数预算已用现有 Python 精确分数只读核算，采用 $\lambda<11/7$：
$\epsilon_\phi<8.094633\times10^{-29}$、
$b_{\rm back}<1.079285\times10^{-96}$、
端点交互均能 $<8.094633\times10^{-29}$。
没有新实验文件或独立校准计数；波包部分尚未据此代验收。
