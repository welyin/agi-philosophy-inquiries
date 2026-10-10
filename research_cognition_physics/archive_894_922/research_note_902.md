# 第902轮：同一物质—几何背景的实际短时演化与能量交换

日期：2026-10-06。接[901](research_note_901.md)与[902冻结入口](902/drafts/STATUS.md)。

配套：[完整背景演化](902/common_background_evolution.py) · [独立作用/传播检查](902/joint_evolution_checks.py) · [全部诊断数据](902/joint_evolution_checks_results.json) · [范围核验代码](902/validate_joint_evolution.py) · [正式结果](902/joint_background_validation_results.json) · [历史与文档核验](902/research_round_902_checks.json) · [复算入口](902/verify_round902.py)。

## 1. 本轮取得什么

**已把原859经典初值推进成一份同时演化完整度规、颜色/弱/圆规范场、五实标量及已有探针的三维周期计算。** 同一作用量的canonical力、度规应力和物质—几何做功通过独立检查；短时演化的约束残差随空间细化明显下降。

这填补了原模型只有初始约束解与局部符号演化、缺少完整未来背景计算实现的一部分缺口。结果是浮点PDE求解和细化证据，**尚无有限时间区间误差证书**；不能将残差小直接解释成同样大小的真实解误差。

|层次|本轮地位|
|---|---|
|认知动机|同一物质同时提供参考和来源，几何与物质的能量交换必须相容|
|继承输入|四维Einstein—规范—H5—探针作用、群和参数、原周期初值、ε=.02|
|复用|572/573约束和局部发展；753颜色资料；730初始时间数据；859/860背景；901初片误差方法|
|新增实现与证据|完整三维耦合未来计算、全非Abelian力与完整应力、原背景的短时空间/时间细化|
|未完成|连续时间误差包络、参考支持区统一证书、实际量子模式、872总来源及899的M_T|

预置Einstein作用后检验共同演化不是生成引力。这里经典费米和接收器字段为零，与原分支一致；量子涨落并未被证明可以忽略，也未被本程序求解。860和869的新顶点在这份零费米经典背景不改原经典方程，其量子作用仍需另接。

## 2. 一份共同作用，而非分别指定各场的力

记Φ=(φ¹,…,φ⁵,p)，G=diag(𝒦,1)，U_total=U_H5+p²/2。沿原Einstein单位与Hermitian规范约定，经典核心为

$$
S_{\rm cl}=\int\sqrt{-g}\left[{R\over2}
-{1\over2}G_{AB}D_\mu\Phi^A D^\mu\Phi^B-U_{\rm total}
-{1\over4}\sum_a K_a F^a_{\mu\nu}F^{a\mu\nu}\right]d^4x .
\tag{1}
$$

K_a在每个简单部门恒定；颜色、弱、圆三块分别为1/g_c²、1/g_w²、36/g_Y²，与原K=1/(2b)相同。规范连接吸收耦合，圆荷用Q=6Y。五标量保完整角向，不能只演化h与s。探针中性；经典颜色物质源为零，但颜色电磁应力保留。

在lapse α、shift β、空间度规γ下，沿原canonical密度π_A、E^{ia}定义

$$
\begin{gathered}
v^A={G^{AB}\pi_B\over\sqrt\gamma},\qquad
e_i^a={\gamma_{ij}E^{ja}\over K_a\sqrt\gamma},\qquad
M_i=\pi_A D_i\Phi^A+E^{ja}F^a_{ij},\\
\rho={1\over2}G_{AB}(v^Av^B+\gamma^{ij}D_i\Phi^A D_j\Phi^B)
+U_{\rm total}+{1\over2}\sum_aK_a\gamma^{ij}e_i^ae_j^a
+{1\over4}\sum_aK_aF^a_{ij}F^{aij} .
\end{gathered}
\tag{2}
$$

e_i=F_{ni}是物理法向电场，不是canonical E。内部时间规范A₀=0时，整份物质Hamiltonian为

$$
H_m=\int_{T^3}\left(\alpha\sqrt\gamma\rho+\beta^iM_i\right)d^3x,\qquad
\dot\Phi={\delta H_m\over\delta\pi},\quad
\dot\pi=-{\delta H_m\over\delta\Phi},\quad
\dot A_i={\delta H_m\over\delta E^i},\quad
\dot E^i=-{\delta H_m\over\delta A_i} .
\tag{3}
$$

尤其Φ̇=αv+βⁱD_iΦ，Ȧ_i=αe_i+βʲF_ji；动量方程包含全部目标度量导数、势梯度、非Abelian交换子和shift荷项。代码的独立Hamiltonian差分特意加入非零shift，防止仅在β=0时碰巧正确。

应力来自同一(1)：标量部分为G_AB D_iΦ^A D_jΦ^B+γ_ij[(vGv−|DΦ|²_G)/2−U_total]；规范部分为Σ_a K_a[−e_i^ae_j^a+γ^{kl}F^a_ikF^a_jl+γ_ij(e²_a/2−F²_a/4)]。按T_ni=M_i/√γ、T_nn=ρ与空间应力S_ij拼回完整T_μν。

## 3. 完整度规的谐和演化

[730](../archive_702_741/730/joint_dynamic_continuum_reference.py)已给原ADM符号K_ij=−γ̇_ij/2及初始jet；805明确其局部符号插值不是原未来解。本轮用谐和坐标C_μ=g_μν g^{αβ}Γ^ν_αβ=0，解

$$
R_{\mu\nu}-\nabla_{(\mu}C_{\nu)}
=T_{\mu\nu}-{1\over2}g_{\mu\nu}T,\qquad
R_{\mu\nu}-\nabla_{(\mu}C_{\nu)}
=-{1\over2}g^{\alpha\beta}\partial_\alpha\partial_\beta g_{\mu\nu}
+Q_{\mu\nu}(g,\partial g).
\tag{4}
$$

Q由Christoffel表达的二阶导数零部分直接构造，源使用(2)—(3)同一物质。将20组任意对称二阶metric jets代入完整Ricci和谐和项，与右端独立比较，最大差2.78×10⁻¹⁶；不靠对角度规或单一波矢校准。

方法对接[Sarbach—Tiglio，2012](https://arxiv.org/abs/1203.6443)的双曲约化、约束传播和离散化框架。旧573的局部存在方法继续承担连续存在论；内部时间规范可由光滑局部规范变换取得。数值时间规范实现与谐和Einstein方程的联合全误差稳定性，尚未由此自动证明。

初始α=1、β=0、γ=ψ⁴δ，采用原K=ψ⁻²Ã+(τ/3)γ。同一Cauchy资料的谐和初始时间导数为

$$
\partial_tg_{ij}=-2K_{ij},\qquad
\partial_tg_{00}=2\tau,\qquad
\partial_tg_{0i}=\gamma_{ij}\,{}^{(3)}\Gamma^j
=-2\psi^{-1}\partial_i\psi .
\tag{5}
$$

代入C_μ得到零。与730的Gaussian时间不同的是坐标规范，物理初始γ、K和canonical物质资料相同。未来不限制空间度规继续共形平坦，不固定lapse/shift为初值，也不把全部原场简化成独立振子。

连续层面同一物质方程给∇^μT_μν=0，Bianchi使C满足齐次约束传播；初始Einstein约束及C=0给相容法向资料。规范Gauss由同一作用的内部对称性传播。这里复用成熟身份，不将数值小残差当成它们的证明。

## 4. 与作用独立交叉检查的能量交换

保持canonical资料，在初始α=1、β=0时对空间度规作变分：

$$
\delta H_m=-{1\over2}\int\sqrt\gamma\,S^{ij}\delta\gamma_{ij}\,d^3x,
\qquad
{d\over dt}\int\sqrt\gamma\rho\,d^3x
=\int\sqrt\gamma\,K_{ij}S^{ij}\,d^3x .
\tag{6}
$$

第二式在此初始切片成立：canonical场演化的Hamiltonian配对相消，空间散度在环面积分为零，γ̇=−2K留下右端。它不是一般shift/lapse下的同样写法，也不是物质能量单独守恒。

四类canonical变分在最后差分步的误差均小于2.8×10⁻¹¹；度规应力变分误差约1.43×10⁻¹¹。原N=12背景的右端约154.110700612846（模型单位），沿完整耦合向量场计算左端，步长1e−4、5e−5、2.5e−5的误差分别为2.50e−6、6.24e−7、1.54e−7，呈二阶差分下降。

这核对同一原作用中的物质—几何能源交换；该功是内部几何耦合，不是加入外部能源。数值检查不是该数值等于连续真实值的区间证明。

## 5. 实际短时计算

空间用周期Fourier配点导数，时间用四阶Runge–Kutta；终点T=.01为原模型坐标单位。N=12、16、24、32各用8步，另在N=24用16步。所有部门同时演化；每个N复用原859约束求解器准备初值，没有冻结非Abelian电场、标量角向或shift。

独立监测

$$
\begin{gathered}
\mathcal H={}^{(3)}R+K^2-K_{ij}K^{ij}-2\rho,\qquad
\mathcal M_i=D_jK^j{}_i-\partial_iK+{M_i\over\sqrt\gamma},\\
\mathcal G_a=D_iE^{ia}+\pi_A(R_a\Phi)^A,
\qquad C_\mu=g^{\alpha\beta}\partial_\alpha g_{\beta\mu}
-\tfrac12g^{\alpha\beta}\partial_\mu g_{\alpha\beta}.
\end{gathered}
\tag{7}
$$

|N³，8时间步|终点𝓗最大绝对值|终点𝓜最大绝对值|终点𝓖最大绝对值|终点C最大绝对值|
|---:|---:|---:|---:|---:|
|12³|2.348e−3|1.621e−4|7.763e−8|2.302e−5|
|16³|1.064e−4|1.055e−5|3.982e−9|7.232e−7|
|24³|2.153e−7|2.684e−8|8.406e−12|1.456e−9|
|32³|3.278e−10|6.515e−11|2.049e−14|2.292e−12|

表头中的max表示网格最大值。N=24时间步减半，终点全状态的最大分量差2.82×10⁻¹²，度规部分4.81×10⁻¹³；没有借此推断未解析空间误差也同样小。

初始𝓗独立重建残差同样从2.291e−3降到3.301e−10。它高于原共形方程内部配点残差，是因为独立重建先形成非线性度规、再计算空间曲率，存在不同的频率混叠。这些数据均保留，不能写成“初始约束在所有离散尺度严格为零”。

终点网格上的γ最小本征值> .730、F>1.862，模型在这组数值轨迹中未离开正支。仅是采样轨迹证据，不是全时空严格下界。

## 6. 材料参考也随同一演化改变

读取仍使用859的完整内部磁标量

$$
M_h={1\over2}P^{\mu\rho}P^{\nu\sigma}
\sum_{a\in\mathrm{color}}F^a_{\mu\nu}F^a_{\rho\sigma},\qquad
P^{\mu\nu}=g^{\mu\nu}+n_h^\mu n_h^\nu,\qquad
n_h^\mu={-g^{\mu\nu}\partial_\nu h\over\sqrt{-(dh)^2}}.
\tag{8}
$$

未来保实际h法向、非零shift和颜色电场，不能继续使用特殊初片上的S_cψ⁻⁸表达。在原坐标点(0,π/2,π/4)，N=32的诊断为：

|量|初始|T=.01|
|---|---:|---:|
|h|.698715113865|.698788307696|
|s|.545863690227|.543611737219|
|M_h|.007099480533|.040069411207|
|−(dh)²|5.38579e−5|4.34426e−5|

N=24和32终点M_h相差约1.33e−8。这是同一坐标点的场与参考诊断，不是已实现的量子仪器输出。没有证明未来参考Jacobian处处满秩。N=12索引实际对应z=π/6，正式结果另列坐标元数据；该行点值不参与原π/4点的跨网格比较。

## 7. 证据强度与下一步

本轮填入真实三维背景发展，且演化、应力和能量交换共同使用原作用。与730局部64维符号或899有限振子校准不同，它实际计算了原周期空间中的全部上述经典场。

**尚不能把本轮当成902范围的连续严格误差证书。** 具体还差：

1. 把901已证明的初始近似及其误差运输到本谐和/时间规范的全部Cauchy变量；当前各N仍各自求初值，未宣称901包络已自动覆盖所有时间变量和所需高范数。
2. 将整段近似轨迹的网格间、时间插值及频率尾部残差接入适用的稳定性估计。守恒约束小不能单独替代解误差。
3. 在原860/869固定材料规则及支持上，推进实际费米/玻色模式与872共同源；当前没有算量子M_T或真实有限耦合的记录概率。

这些是当前有限预测所需的连接，不是无截断UV完成门槛。旧全局表示等障碍不重新加回。代码保存紧凑参数与结果，不保存大型场数组；全部场可按同一初值与参数复算。

接[903入口](903/drafts/STATUS.md)，优先把当前实际轨迹与同一原量子传播的输入合同接通，并为必要的背景误差选定可执行证书；不再仅用越来越细网格的漂亮残差代替证明。正式902，新增1组，累计3687；完整统一目标仍未完成。
