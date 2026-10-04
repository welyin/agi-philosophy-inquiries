# 第583轮：共同真空中的联合涨落与冻结背景

日期：2026-10-01。接[582共同有效项](../../research_note_582.md)、[580冻结变量审计](../../research_note_580.md)及[合并清单](../../582/joint_condition_compression_update_582.md)。链接按archive_231_解析；完成以独立终审和冻结核验为准。

## 1. 本轮解决哪个连接

553在固定Jordan度规上处理物质量子修正，574与580—582使用固定Einstein度规。580已经指出这两个冻结方向不同；581—582的经典／一阶EFT场重定义没有解决部分量子化的匹配。

本轮在**原模型的同一个平直、零势、驻定真空**给出明确结果：两个固定度规的标量二次算符不同，统一动能标准化后仍有不同的径向谱；保留度规—标量混合的联合二次型则按同一变量变换对应。直接代入原非线性作用进一步给出可行接法：固定条件本身也必须拉回，固定Einstein度规对应随标量变化的Jordan度规。

这排除了“共同在壳背景足以让两种固定度规的标量行列式相同”的直接拼接。它不是完整量子框架不等价，更不是量子引力构造完成。

## 2. 文献映射与范围

[Kamenshchik与Steinwachs，1408.5769v3](https://arxiv.org/pdf/1408.5769)在第II节将多场结果约化到单标量，第V节显式检验一个宇宙学背景的在壳一致性；附录D计算含度规、标量、规范固定及ghost。本项目不能把该完整变量比较当成不同冻结子块相同。文中Einstein框架还要求标量canonical；本项目保留非平直五场目标，不套用单场全局展平公式。

[Finn、Karamitsos与Pilaftsis，1910.06661v4，第IX节](https://arxiv.org/pdf/1910.06661)把度规与标量纳入共同配置空间；式IX.12—18还涉及测度、规范固定、FP因子及模型函数ell。本轮不计算这些完整对象，六维切向块的行列式关系不替代其量子构造。

|层次|本轮地位|
|---|---|
|认知动机|同一个物质—几何系统的不同表示应共同转换操作与背景条件|
|输入|原五实场K、U与领先四维标量—张量作用，F>0分支，周期边界及固定共同单位|
|范围|零规范背景、共同常真空、标量与共形度规方向；规范场和其它物种未作内部环计算|
|解析增量|两个冻结子块及其统一标准化、原模型联合混合、实际非线性冻结条件拉回|
|未完成|完整规范固定Hessian、ghost与测度、量子约束、非零共同来源和原图连续匹配|

Schur补、二次型变换及Weyl变换是成熟工具。新增的是把它们映射到原共同模型并定位旧两套计算不能直接合并的具体部分。

## 3. 同一作用与同一真空

按既有Euclidean、κ=1约定，保留原参数和势；在当前单位中gE=F gJ。零规范背景的领先作用为

$$
I_J=\int\sqrt{g_J}\left[-\frac12FR_J+\frac12\delta_{ab}\partial\phi^a\cdot\partial\phi^b+V\right],\qquad
I_E=\int\sqrt{g_E}\left[-\frac12R_E+\frac12K_{ab}\partial\phi^a\cdot\partial\phi^b+U\right],
\quad U=\frac{V}{F^2},\quad K=\frac{I}{F}+\frac{3}{2F^2}dF\otimes dF.
\tag{1}
$$

F=M−φ²/6、M=2，V=dᵀLd/4，d=(h²−u_h,s²−u_s)，L及u直接复用旧参数，h²=Σ前四个φ分量平方。本轮不另选质量。取

$$
\phi_0=(0,\sqrt{u_h},0,0,\sqrt{u_s}),\quad
V_0=U_0=0,\quad dU|_0=0,\quad
g_{E0}=\delta,\quad g_{J0}=\frac{\delta}{F_0},\qquad F_0=1.8765364152074429.
\tag{2}
$$

两度规不能同时取δ再以同一坐标动量比较。该背景满足所列领先标量与度规方程；不是572非均匀有规范应力来源的替代解。

记Cartesian势Hessian为H_U。它只有h、s径向块非零，因为V及dV在真空消失：

$$
(H_U)_{\{h,s\}}=\frac{2}{F_0^2}
\begin{pmatrix}
L_{hh}u_h&L_{hs}\sqrt{u_hu_s}\\
L_{hs}\sqrt{u_hu_s}&L_{ss}u_s
\end{pmatrix},\qquad
\operatorname{rank}H_U=2.
\tag{3}
$$

原L正定，径向块正；三个角向零方向为当前冻结规范背景的标量方向，不将它们直接计作完整规范理论的物理粒子。

## 4. 冻结不同变量时丢失了什么

在共同坐标中取gE=e^(2σE)δ、gJ=e^(2σJ)δ，背景σE0=0、σJ0=−logF0/2。周期边界分部积分后，原Einstein作用沿共形方向为

$$
I_E[\sigma_E,\phi]=\int\left[-3e^{2\sigma_E}(\partial\sigma_E)^2
+\frac12e^{2\sigma_E}K_{ab}\partial\phi^a\partial\phi^b+e^{4\sigma_E}U\right],
\qquad \sigma_E=\sigma_J+\frac12\log F(\phi).
\tag{4}
$$

因此在线性阶

$$
\delta\sigma_E=\delta\sigma_J+c^\top\delta\phi,\qquad
c=\frac{dF|_0}{2F_0}=-\frac{\phi_0}{6F_0},\qquad
K_0=\frac I{F_0}+6cc^\top.
\tag{5}
$$

一般变量变换的普通Hessian满足

$$
\frac{\delta^2 I_J}{\delta q^i\delta q^j}
=\frac{\partial Q^A}{\partial q^i}\frac{\delta^2 I_E}{\delta Q^A\delta Q^B}\frac{\partial Q^B}{\partial q^j}
+\frac{\delta I_E}{\delta Q^A}\frac{\partial^2Q^A}{\partial q^i\partial q^j}.
\tag{6}
$$

当前共同在壳背景使第二项消失，联合Hessian可作合同变换。取非零Fourier模p²>0，在“共形幅度＋五标量”六个切向坐标中：

$$
H_E=\begin{pmatrix}-6p^2&0\\0&D_E\end{pmatrix},\quad D_E=p^2K_0+H_U,\quad
J=\begin{pmatrix}1&c^\top\\0&I_5\end{pmatrix},\quad
H_J=J^\top H_EJ=\begin{pmatrix}-6p^2&-6p^2c^\top\\-6p^2c&C_J\end{pmatrix},\quad
C_J=\frac{p^2}{F_0}I_5+H_U.
\tag{7}
$$

固定σE得到D_E，固定σJ得到C_J。它们相差6p²ccᵀ，不能因背景在壳而忽略。保留混合后

$$
C_J-H_{J,\phi\sigma}H_{J,\sigma\sigma}^{-1}H_{J,\sigma\phi}=D_E,
\qquad \det H_J=\det H_E,\qquad (H_J^{-1})_{\phi\phi}=D_E^{-1}.
\tag{8}
$$

式(8)只是在p²>0的该受限块上求Schur补；不是已经完成度规路径积分。两联合块均有一个Euclidean共形负方向，不能当成正Gaussian协方差或健康独立模。没有处理完整张量方向、规范固定、ghost或约束，p=0的共形零模也不能在此求逆。

## 5. 不是单纯的动能归一或单位差

同一坐标p²下，由矩阵行列式引理及H_U径向正性，

$$
\frac{\det D_E}{\det C_J}=1+6p^2c^\top C_J^{-1}c,
\qquad 1<\frac{\det D_E}{\det C_J}<\frac{M}{F_0},
\qquad \lim_{p^2\to\infty}\frac{\det D_E}{\det C_J}=\frac M{F_0}.
\tag{9}
$$

不能把这个高动量常数直接当物理效应；先分别标准化两个标量动能，定义

$$
\mathsf M_E=K_0^{-1/2}H_UK_0^{-1/2},\qquad
\mathsf M_J=F_0H_U,\qquad
\frac{\det(p^2I+\mathsf M_E)}{\det(p^2I+\mathsf M_J)}
=\frac{F_0}{M}\frac{\det D_E}{\det C_J}
\ \xrightarrow[p^2\to\infty]{}1.
\tag{10}
$$

低动量比值趋F0/M，因为c位于正径向子空间。三个共同角向零方向的p²因子相消；本式以p²>0取极限，不直接对零模求行列式逆。

实际旧参数给两个标准化径向质量平方，均按共同Einstein背景坐标单位表示：

$$
\operatorname{spec}_{+}\mathsf M_E=(0.0468385104353,\ 0.138056726486),\qquad
\operatorname{spec}_{+}\mathsf M_J=(0.0468422679361,\ 0.147128135649).
\tag{11}
$$

例如p²=.1时标准化行列式比=.963268040290，仍非1。这是两个不同冻结处方的部分标量谱差，不是两个完整物理框架预测不一致；也不是实验极点质量或质量拟合。动能标准化排除一项平凡常数差，不等于已选择完整配置空间测度。

## 6. 原非线性作用给出的明确接法

固定gJ=δ/F0时，σE=log(F/F0)/2。式(4)中共形梯度项恰好抵消K的额外径向部分，所以在周期边界或相应边界项处理下

$$
I_E\!\left[\frac12\log\frac F{F_0},\phi\right]
=I_J\!\left[\frac{\delta}{F_0},\phi\right]
=\int\left[\frac{(\partial\phi)^2}{2F_0}+\frac V{F_0^2}\right].
\tag{12}
$$

反过来，固定gE=δ在Jordan变量中必须使用gJ(φ)=δ/F(φ)，于是

$$
I_J\!\left[\frac{\delta}{F(\phi)},\phi\right]
=I_E[\delta,\phi]
=\int\left[\frac12K_{ab}\partial\phi^a\partial\phi^b+U\right].
\tag{13}
$$

这是同一受限配置族的精确作用映射，无需先宣称引力已被量子化。它明确给出“保留固定Einstein物质处方”时另一表示必须如何改变；传统固定Jordan背景并不满足此条件。是否能将其提升为相同量子处方，还要共同给测度、正则化、源、边界和读取，不仅代入同一个经典作用。

## 7. 实际非线性检查

[代码](../joint_frame_hessian_matching.py)和[保存结果](../joint_frame_hessian_matching_results.json)只用原Python与NumPy。实际动作积分直接使用原U、K和Jordan的F、V，并分别计算原曲率表达式与分部积分表达式，没有用目标二次型生成非线性结果。

在周期x∈[0,2π)、横向坐标体积1上，取φ=φ0+ηv cos(nx)、σJ=σJ0+ηa cos(nx)，n=2、v=(.13,.8,−.11,.07,.6)。原作用的中心二阶变分应满足

$$
\lim_{\eta\to0}\frac{I_J(\eta)+I_J(-\eta)-2I_J(0)}{\eta^2}
=\pi\begin{pmatrix}a\\v\end{pmatrix}^{\!\top}H_J(n^2)\begin{pmatrix}a\\v\end{pmatrix},
\qquad
\delta^2 I_{\text{fixed }J}-\delta^2 I_{\text{fixed }E}=-6\pi n^2(c^\top v)^2.
\tag{14}
$$

|五组复算|结果|
|---|---|
|原真空与五维势Hessian|原U直接中心差分，步长.002到.0005，最大误差1.07×10⁻⁶降至6.69×10⁻⁸；径向秩2|
|联合二次型与混合|五个p²，合同／Schur误差≤8.89×10⁻¹⁶，行列式比误差≤1.89×10⁻¹⁵；联合逆的标量块误差≤4.55×10⁻¹³；全部惯性为1负、5正|
|统一动能标准化|七个p²核精确行列式关系；高动量归一比趋1，而有限p²仍不同；径向谱见式(11)|
|两种固定背景的非线性作用|二阶极限分别7.609647822、7.169890755，解析差−.4397570672；η减半误差约缩四倍，非线性框架及冻结条件拉回的积分差均小于10⁻¹⁶|
|联合度规—标量路径|a=.17、−.21时，二阶极限6.948671108、1.426383833；直接原Jordan与Einstein作用积分一致，有限η二阶余差按η²收敛|

这些数值检验支持前述解析等式，不以有限样本替代证明。浮点谱中的10⁻¹⁸量级负数是三个解析零方向的舍入误差。

## 8. 条件账与下一步

[583条件账](../unified_physics_condition_ledger_583.md)将“量子处方相同”拆成可核义务：相同配置族、相同受限积分变量、正确混合和标准化、共同测度／正则化／源。这些是比较既有处方的条件，不自动成为新的认知公理。

本轮已给原作用层面的接通，并排除即使在共同真空也直接交换两个冻结子块的做法。尚未证明完整量子匹配或553全部物种运行与574固定图量子化相同；578—579的固定ℏ细化能源问题也仍在。

下一步先检验能否在**同一有限截止**下，把固定条件的拉回、量子测度及实际算符一起保持；这比直接宣称全部量子引力等价更接近当前可审计缺口。若需要另选测度、排序或源，明确标为新增输入并检验效果。随后将共同应力与几何连接，而非继续扫描本轮两个质量数值。
