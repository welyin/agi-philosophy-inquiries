# 第549轮：同一物质参考与动态几何——局部约束相容构造

日期：2026-09-30。继承[548](research_note_548.md)、[条件压缩表](548/joint_condition_compression_table.md)；[代码](549/joint_reference_gravity_constraints.py)、[结果](549/joint_reference_gravity_constraints_results.json)、[核验](549/research_round_549_checks.json)、[条件账](549/unified_physics_condition_ledger_549.md)。

## 1. 本轮合并了什么

548在固定Minkowski几何中证明：已有Higgs径向场h和实中性场s的值与导数，可以给四个局部关系参考。现在证明一个更强但仍有边界的命题：**在明确给定的四维、二导数标量—张量作用中，存在同一h、s和动态几何的局部解析解，既满足引力初始约束，又使这四个参考量满秩。**

因此，在这个经典局部实现层面，“物质提供参考”和“参考参与几何反作用”不必再由两个不相容的背景模型分别承担。不需要新添四个基本参考场，也不减去势的常数项。它是同模型联合存在，既不是条件等价，也不是从认知原则生成Einstein作用。

|层次|内容|
|---|---|
|认知动机|参考、物质和几何源应在同一模型内核算|
|继承输入|四维Lorentz流形、h＋s场内容、非最小曲率耦合和给定势；s仍是543显式添加的候选|
|本轮附加范围|二导数经典分支；F、V实解析，取h>0、F>0局部片；选择解析非均匀初值|
|解析增量|用协变场动量消掉动量约束源；解出局部Hamiltonian约束；构造真正耦合解上的满秩关系图|
|复算作用|检查归一化、完整场方程二阶导数、Ricci定义、混合动量和参考雅可比；有限导数数据不冒充完整PDE数值解|
|仍独立|维数与作用来源、谱截断的受控性、曲率系数运行、状态制备、全局延拓和实际量子读取|

## 2. 对接成熟工具及历史去重

使用多标量Einstein框架，取[Kaiser等，1210.7487，式(1)—(9)](https://arxiv.org/pdf/1210.7487)中的f=F/2、目标Planck单位为1。该变换保留非对角场空间度量，不能把两场当成独立自由标量。其后文的宇宙学专门解不是本轮定理的前提或证明。

约束符号、三维初始片共形指数可核[Choquet-Bruhat、Isenberg、Pollack，gr-qc/0610045，§2式(4)—(8)](https://arxiv.org/pdf/gr-qc/0610045)。这里只继承约束和共形工具，不借用其紧致单标量存在分类，亦不把另一篇渐近Euclidean存在定理强加到本地两场问题。

解析约束存在用标准Cauchy–Kowalevskaya（CK）定理；非特征解析初值的表述见[Zheng-Chao Han，PDE讲义，§9.4—9.5](https://sites.math.rutgers.edu/~zchan/current-preprint/PDEIntro25.pdf)。下文明确写出法向最高导数可解的形式。

303、319—325已有正F、变分和守恒工具；533已有谱系数；548已有导数参考方法。新结果是它们在同一局部受约束状态中的连接。没有把成熟定理、重新整理或场变量替换计作新发现。

## 3. 同一作用、曲率归一和常势项

采用号差(−+++)，在h>0的H=(0,h/√2)实方向片，规范场、费米子取零；Higgs规范电流为零，故继承548的一致经典子部门。作用定义为

$$
I_J=\int\sqrt{-g_J}\left[\frac{F(h,s)}2R_J-\frac12\delta_{IJ}\partial\phi^I\cdot\partial\phi^J-V(h,s)\right]d^4x,
\qquad \phi=(h,s).
\tag{1}
$$

这已经输入了引力作用的形式。一般存在定理只需F、V解析且F>0；接入533、543的同尺度裸归一时，取f₀T=2π²、T₂=2T(h²+s²)，得到

$$
F=2K=M_0^2-\frac{h^2+s^2}{6},\quad
C_0=\frac{4f_2\Lambda^2}{f_0},\quad M_0^2=\frac{nC_0}{24T},\quad
V=V_0-\frac{C_0(h^2+s^2)}2+\frac{\lambda_Hh^4+2p h^2s^2+\lambda_s s^4}{4},\quad
V_0=\frac{nf_4\Lambda^4}{2\pi^2}.
\tag{2}
$$

n=Tr Z必须包括所有代和零Yukawa态；在536的N=3、n_g代、夸克权1与轻子权r的约定下，n=8n_g(3+r)。不能把一代的n套到三代中性谱。常数V₀作为几何源保留。

**有效阶数单列：** 533完整截断另有Weyl平方等项，完整谱作用还有更高阶。式(1)只是明确选出的二导数分支。本轮不证明高阶修正小，不称其解为完整谱作用解；也没有把548的低能参数直接与尚未运行的裸F相配。下面数例使用同一裸尺度，属于形式相容见证。

## 4. 混合场空间度量与约束

令g_E=Fg_J，目标Planck单位固定为1。除局部全导数外，式(1)变为

$$
I_E=\int\sqrt{-g_E}\left[\frac{R_E}2-\frac12G_{IJ}\partial\phi^I\cdot\partial\phi^J-U\right]d^4x,
\quad G_{IJ}=\frac{\delta_{IJ}}F+\frac{3F_I F_J}{2F^2},\quad U=\frac V{F^2}.
\tag{3}
$$

F>0使G正定，但通常G_hs不为0。记初始片法向速度πᴵ=nᵘ∂ᵤφᴵ，取空间度规γ、外曲率K。Einstein约束是

$$
R(\gamma)-|K|^2+(\operatorname{tr}K)^2
=G_{IJ}\pi^I\pi^J+G_{IJ}\gamma^{ij}\partial_i\phi^I\partial_j\phi^J+2U,
\qquad \nabla^jK_{ij}-\partial_i\operatorname{tr}K=-G_{IJ}\pi^I\partial_i\phi^J.
\tag{4}
$$

把548的“s速度为零”直接搬过来一般失败：h有速度、s有空间梯度时，右侧会出现G_hs交叉源。正面构造要选择的是协变动量，而非凭直觉清零某个速度。

## 5. 局部解析初值的构造

固定任意(h★,s★)，h★>0、F★>0，不要求它是任何势的驻点。用固定单位，取ε>0，在三维初始片原点附近设置

$$
\gamma_{ij}=\psi^4\delta_{ij},\quad K_{ij}=0,\quad
h=h_\star,\quad s=s_\star+\varepsilon(x+xz),\quad
\mathfrak p_I:=\psi^6G_{IJ}\pi^J=P\delta_{Ih},\quad
P=\frac{\varepsilon(1+y)}{(G^{hh})_\star}.
\tag{5}
$$

于是πᴵ=ψ⁻⁶GᴵʰP，πˢ一般非零。整个局部初始片上的动量源严格为

$$
j_i=-\psi^{-6}P\partial_i h=0,
\qquad
\Delta_\delta\psi=-\frac18\left[\psi^{-7}G^{hh}P^2+\psi G_{ss}\varepsilon^2\big((1+z)^2+x^2\big)+2U\psi^5\right].
\tag{6}
$$

第二式由R(ψ⁴δ)=−8ψ⁻⁵Δψ直接得到，包括全部势。接着在空间平面z=0给解析数据，并将最高z导数解出：

$$
\psi(x,y,0)=1,\qquad \partial_z\psi(x,y,0)=0,\qquad
\psi_{zz}=-\psi_{xx}-\psi_{yy}-\frac18\left[\psi^{-7}G^{hh}P^2+\psi G_{ss}\varepsilon^2((1+z)^2+x^2)+2U\psi^5\right].
\tag{7}
$$

右侧在ψ=1、F>0附近对所有自变量和所列低阶导数解析；ψ_zz的系数为1。CK定理给原点附近的实解析解。缩小邻域可保持ψ>0、h>0、F>0。由式(5)—(7)，全部初始约束在这个开邻域内成立，而不只是某点的代数残差为零。

这里是椭圆方程的**局部解析存在构造**。它不意味着对一般光滑、带噪声平面数据适定，不保证全初始片、紧致／渐近平直边界或有限总能量。尤其不能像548的固定背景问题那样直接乘紧支撑cutoff：这种做法通常破坏引力约束。

## 6. 由初始约束到同一动态解

将Einstein方程作谐波约化，并联立sigma模型方程

$$
\Box_E\phi^I+\Gamma^I{}_{JK}(G)\,\partial\phi^J\cdot\partial\phi^K-G^{IJ}U_J=0.
\tag{8}
$$

主部是Lorentz度规的波算子；初始片为空间片，G正定可逆，系数解析。选与初始几何兼容的谐波规范数据后，解析非特征局部存在给约化方程的解。物质方程给∇T=0，Bianchi恒等式使规范约束满足齐次传播方程；初始Hamiltonian和动量约束保证初始规范约束及其法向导数为零，故所得解满足未约化Einstein方程。

这使用成熟局部发展方法；不需要将单标量、特定边界或全局存在结论误推广到这里。之后可换为初始片附近的Gaussian时间坐标，取初始 lapse=1、shift=0。在原点ψ=1、∂ᵢψ=0、K=0，因此度规一阶导数及全部时空Christoffel在该点为0。

## 7. 动态几何没有消掉参考的满秩

记c=(Gˢʰ/Gʰʰ)★=−(G_hs/G_ss)★。由初值和上节坐标选择，在原点o有

$$
dh=(\varepsilon,0,0,0),\quad ds=(c\varepsilon,\varepsilon,0,0),\quad
h_{0y}=\varepsilon,\quad s_{0y}=c\varepsilon,\quad h_{0z}=s_{0z}=0,\quad s_{xz}=\varepsilon.
\tag{9}
$$

h的纯空间Hessian为0，s的纯空间对角二阶导数也为0。h₀ₓ、s₀ₓ由G⁻¹的s导数决定；h₀₀、s₀₀由式(8)决定，均不能擅自清零。**即使V_I=0，U_I=−2VF_I/F³也通常非零。** 保留真空常数时尤其如此。

对Einstein框架菜单X_E=(h,s,A_E,B_E)，A_E=(∂h)²_E、B_E=(∂s)²_E，直接得到

$$
J_E\big|_o=\begin{pmatrix}
\varepsilon&0&0&0\\c\varepsilon&\varepsilon&0&0\\
*&*&-2\varepsilon^2&0\\ *&*&-2c^2\varepsilon^2&2\varepsilon^2
\end{pmatrix},\qquad \det J_E\big|_o=-4\varepsilon^6\ne0.
\tag{10}
$$

星号包含完整势与场空间联络，影响前两列，不影响行列式。原Jordan菜单为X_J=(h,s,FA_E,FB_E)；该菜单目标空间变换的雅可比是块三角矩阵，行列式为F²，所以

$$
\det J_J\big|_o=F_\star^2\det J_E\big|_o=-4F_\star^2\varepsilon^6\ne0,
\qquad (dh)^2_E=-\varepsilon^2<0,\qquad(dh)^2_J=-F_\star\varepsilon^2<0.
\tag{11}
$$

这是完整链式法则，已包括∂F项。局部逆函数定理和连续性给非零邻域中的关系坐标，并让h保持类时标签。参考字段仍使用度规，不是生成度规的独立输入。它们本身也不是固定抽象点上的Dirac可观测量；仍应如548那样组合为f∘X⁻¹和度规拉回。没有保证全局单射、所有状态、长期参考或实际仪器。

## 8. 同一裸尺度的复算见证

选N=3、n_g=3、r=7/3、T=1、q=1/4、S=3/28、C₀=1/4，以及正单调截止f(u)=f₀exp(−u²)。于是f₂=f₄=f₀/2、Λ²=1/8；543的同尺度关系和本轮曲率归一给

$$
n=128,\quad M_0^2=\frac43,\quad(\lambda_H,p,\lambda_s)=\left(\frac3{14},\frac3{28},\frac37\right),\quad
(h_\star^2,s_\star^2)=\left(1,\frac13\right),\quad F_\star=\frac{10}9,\quad
V_0=1,\quad V_\star=\frac{11}{12},\quad U_\star=\frac{297}{400}.
\tag{12}
$$

这是人为有理参数见证，不是546的人造低能目标，更不是现实拟合。它的常势项不小，**不提供谱导数展开受控的证据**。局部定理可以在这种形式模型上成立，而将其解释成高能谱模型的低能解仍缺尺度论证。

8组复算覆盖：全代身份迹及常势；G、G⁻¹及其独立差分；协变动量和错误零速度反例；任意前列的精确有理行列式；完整sigma方程与从Christoffel定义计算的空间Ricci；两个框架的菜单差分；改变常势会改变曲率和加速度但不改变秩；ε=0或F≤0的边界。

动量残差最大5.60×10⁻¹⁹，所核场方程残差为浮点0，菜单差分误差分别4.42×10⁻¹²和1.04×10⁻¹¹。差分用的是二阶Taylor代表和原点切空间度规：真实度规一阶导数为0，足以核对菜单的一阶导数；并未求解完整PDE或给数值收敛证明。[独立审查](549/drafts/independent_review.txt)、[终审](549/drafts/final_review.txt)。

## 9. 合并后的真正卡点与接续

本轮关闭的是二导数经典局部模型内“参考与引力约束不能同时满足”的疑问。这里的Einstein动力学已输入，不能签收认知到GR、三维生成或完整谱理论。

下一项优先核**同一裸参数、保留真空常数后，耦合解能否处于可控的曲率尺度**。533的静态系数界要继承；新问题是它是否经完整非最小场方程变成动态曲率约束，以及参考振幅、场空间取向或尺度重标是否真能消除该限制。不先优化参考精度，也不把一个任意减常数的势当作旧模型。若裸分支无受控极限，就记录准确范围，并转向声明清楚的有效匹配分支；不把局部分支失败扩大成统一计划失败。
