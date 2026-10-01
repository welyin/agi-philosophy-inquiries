# 第650轮：原引力约束、面积匹配与物质边界通量的联合条件

日期：2026-10-02。接[649](research_note_649.md)及[起步审计](round650_drafts/boundary_flux_entry.md)，回查[573](research_note_573.md)、[600](research_note_600.md)、[618](research_note_618.md)、[619](research_note_619.md)、[647](research_note_647.md)、[648](research_note_648.md)。[代码](joint_gravity_boundary_flux.py)、[结果](joint_gravity_boundary_flux_results.json)、[核验](research_round_650_checks.json)、[条件账](unified_physics_condition_ledger_650.md)。三组复算、十六式；主代理审查，无新增独立代理审查。

## 1. 联合结果及其范围

**原领先Einstein—五标量作用中，面积匹配不能替代物质—几何的完整边界配对。** 本轮给原径向物质模的具体见证，并证明它们可成为满足全部原初始约束的渐近平直小数据解族的切向。两侧取同一解的限制，面积匹配精确成立；在真空点，两条切向的面积一阶变化均为零，单侧边界辛通量仍非零。

在不增加边界相空间、不限制这些变化的指定标准区域辛形式上，这使时间平移的候选Hamiltonian一形式不闭合。它排除的是“面积匹配就足以获得独立封闭区域演化”的接法，不排除有通量的开放区域、边界扩展、其它边界条件或整体演化。

正面连接也明确：618的完整共同边界数据与法向动量匹配，使两侧辛通量相消。区域组合、物质来源和几何边界因而须接受同一套变分条件，不能分别指定。

|层次|本轮采用|
|---|---|
|认知动机|组合后的同一整体应同时保持状态变化、物质来源与边界交换|
|继承输入|原正F目标、势及规范表示；给定3+1维领先Einstein作用|
|检验分支|规范场和费米场为零、Higgs固定实径向方向与singlet；固定光滑类时切口|
|新见证输入|渐近平直R³初片、时间对称小数据、紧支撑剖面；不替换572的周期来源|
|解析增量|原模的非线性约束完成、面积不可见的非零通量及指定生成元的局部不可积性|
|数值角色|原势Hessian、两框架完整边界配对、有理收缩界和紧支撑通量积分|
|保留缺口|关系区域、全量子引力、完整手征连续过程、全部EFT阶及独立认知实现|

成熟协变相空间方法复用[Iyer—Wald，式(78)—(82)](https://arxiv.org/pdf/gr-qc/9403028)。边界势和角点项不能默认省略；[Harlow—Wu，式(1)—(6)、§3.2](https://arxiv.org/html/1906.08616v3)将允许变分、边界作用及修正辛形式共同处理。本轮采用明确的标准未扩展区域结构检验，不把这种指定结构的障碍扩大成所有边界Hamiltonian均不存在。

## 2. 原作用和原真空的物理径向模式

沿618的Lorentz约定，记目标度量为K，空间外曲率另记k_ij。当前一致截断的作用为

$$
S_E=\int\sqrt{-g_E}\left[\frac12R_E-\frac12K_{AB}(\phi)\partial\phi^A\cdot\partial\phi^B-U(\phi)\right]+\int_{\partial M}\sqrt{|h_E|}\,K_E^{\rm ext},\quad
F=2-\frac{|\phi|^2}{6},\quad K_{AB}=\frac{\delta_{AB}}F+\frac{\phi_A\phi_B}{6F^2},\quad U=\frac V{F^2}.
\tag{1}
$$

V保留原全部常数。原参数函数给正矩阵L及u=(h★²,s★²)，V=(|H|²−u₁,s²−u₂)L(|H|²−u₁,s²−u₂)ᵀ/4。在φ★=(0,h★,0,0,s★)处U=dU=0。固定实径向Higgs的规范电流为零：规范生成元作用沿角向，K保持径向—角向正交；势和目标联络也保持径向子空间。因此零规范场与零费米场构成此处的真实经典截断，不是漏去非零来源。

平直g_E与常φ★是同一作用的真空。令K_r为(h,s)块，原势径向Hessian及归一化模满足

$$
S_r=\frac{2}{F_\star^2}\operatorname{diag}(h_\star,s_\star)L\operatorname{diag}(h_\star,s_\star),\qquad
S_ra=m^2K_ra,\quad a^\top K_ra=1,\qquad
m^2\simeq0.04683851044,\ 0.13805672649.
\tag{2}
$$

将a嵌入五分量的第2、5分量。因原真空梯度及dU均为零，标量应力从扰动二阶开始，δg_E=0满足线性Einstein方程。局部平面波可取

$$
\delta_1\phi=a\cos(\omega t)\cos(kx),\qquad
\delta_2\phi=a\cos(\omega t)\sin(kx),\qquad
\omega^2=k^2+m^2,\qquad \delta_1g_E=\delta_2g_E=0.
\tag{3}
$$

常背景的径向标量变化不是坐标变换：Lie导数作用于常φ★为零。这些波只作为局部线性诊断；下面另外构造紧支撑初值，不把无限平面波称为有限能量渐近平直解。

## 3. 面积读取不到的边界配对

618在类时边界、向外空间法向n上的一形式及其场空间外微分为

$$
\vartheta_B=\Pi^{ab}\delta h_{ab}+\pi_A\delta\phi^A,\qquad
\omega_B(\delta_1,\delta_2)=\delta_1\Pi^{ab}\delta_2h_{ab}+\delta_1\pi_A\delta_2\phi^A-(1\leftrightarrow2),\qquad
\pi_{E,A}=-\sqrt{|h_E|}K_{AB}n\cdot\partial\phi^B.
\tag{4}
$$

在平面x=0，式(3)给每单位切向面积ω_B=k cos²(ωt)，而c=√q_E=F√q_J的两个一阶变化均为零。局部值非零本身还不足以证明全约束解族上的生成元不可积；以下关闭这一缺口。

## 4. 满足原非线性初始约束的二参数族

在R³选半径R的光滑紧支撑函数χ=exp[−r²/(R²−r²)]（r<R），r≥R时为零。φ与初始数据指定为

$$
\phi_{\epsilon,\zeta}(\mathbf x)=\phi_\star+a\chi(r)\,[\epsilon\cos(kx)+\zeta\sin(kx)],\qquad
\gamma_{ij}=\psi^4\delta_{ij},\qquad k_{ij}=0,\qquad p_A=0,\qquad A_i=E^i=0.
\tag{5}
$$

规范和动量约束均严格为零。设B=K_AB(φ)∂_iφ^A∂_iφ^B；在G_μν=T_μν单位下，R(γ)=2ρ及R(ψ⁴δ)=−8ψ⁻⁵Δψ给唯一剩余约束

$$
-8\Delta\psi=B\psi+2U\psi^5,\qquad \psi\longrightarrow1\quad(r\longrightarrow\infty),\qquad B,U\ge0.
\tag{6}
$$

这是先改变物质、再求实际几何约束，并非把旧解乘cutoff后宣称约束仍成立。B,U光滑紧支撑，均为O((|ε|+|ζ|)²)。在1≤ψ≤3/2上定义Newton映射

$$
(T\psi)(\mathbf x)=1+\frac1{32\pi}\int_{B_R}\frac{B(\mathbf y)\psi(\mathbf y)+2U(\mathbf y)\psi(\mathbf y)^5}{|\mathbf x-\mathbf y|}\,d^3y,
\qquad \sup_{\mathbf x}\frac1{32\pi}\int_{B_R}\frac{d^3y}{|\mathbf x-\mathbf y|}=\frac{R^2}{16}=C_R.
\tag{7}
$$

球的Newton势直接积分可核最后一式。若B_+、U_+为全域上界，则

$$
\|T\psi-1\|_\infty\le C_R\left(\tfrac32B_++2(\tfrac32)^5U_+\right)<\tfrac12,\qquad
\|T\psi-T\tilde\psi\|_\infty\le\ell\|\psi-\tilde\psi\|_\infty,\quad
\ell=C_R\left(B_++10(\tfrac32)^4U_+\right)<1.
\tag{8}
$$

于是C₀(R³)中的闭正区间有唯一不动点，位于该小数据分支。Newton正核给ψ≥1，椭圆正则性逐阶给光滑，外部调和给ψ−1=O(1/r)。这不是所有幅度解的全局唯一性定理。

可复算的宽界如下：原输入满足|φ★|<.9、|a|<1.5、0<L<I；χ≤1且|∇χ|≤8/(eR)<3/R。取R=1、k=.7、η=|ε|+|ζ|≤.01，令d=1.5η、F₋=2−(.9+d)²/6。由K≤2I/F₋²及|(|H|²−h★²,s²−s★²)|≤1.8d+d²得到

$$
F\ge F_-=1.8604625,\qquad
B_+=\frac{2d^2(k+3/R)^2}{F_-^2}<0.001780,\quad
U_+=\frac{(1.8d+d^2)^2}{4F_-^2}<0.00005354,\quad
\|T\psi-1\|_\infty<0.000218,\quad \ell<0.000281.
\tag{9}
$$

结果保存精确有理分数；原有限参数的宽界另作数值核对。点采样仅检查实现，没有以采样最大值代替连续上界。此处不用网格求全PDE。

映射在参数和ψ上光滑，I−D_ψT由一致收缩可逆。因此小邻域中的不动点对两参数可微至所需阶，且

$$
\psi_{\epsilon,\zeta}=1+O(\eta^2),\qquad
\partial_\epsilon\psi|_0=\partial_\zeta\psi|_0=0,\qquad
\delta_1\gamma=\delta_2\gamma=0,\qquad \delta_1k_{ij}=\delta_2k_{ij}=0.
\tag{10}
$$

沿573已核的正F拟线性波主部及Bianchi约束传播，光滑AF小数据在适当规范中有共同短时间发展；可取更高初始正则性保证参数切向的线性发展。只需局部存在，不借用任何小数据全局稳定定理。原线性应力为零，故可选δg_E=0的线性发展。此构造避开紧T³有Killing场时的线性化稳定性障碍，没有把AF结论移用于旧周期初片。

## 5. 非零积分通量与面积匹配的准确量词

初始切面x=0上χ的法向导数为零，δ₁φ=aχ、δ₂φ=0、∂_xδ₂φ=akχ。因此单侧通量及积分为

$$
\omega_B(\delta_1,\delta_2)|_{t=0}=k\chi(\sqrt{y^2+z^2})^2,\qquad
\mathcal F=2\pi k\int_0^R r\exp\!\left[-\frac{2r^2}{R^2-r^2}\right]dr>0,
\qquad \mathcal F\simeq0.609908597826\quad(R=1,k=.7).
\tag{11}
$$

半空间在无穷远的本组一阶扰动为零，不产生抵消该值的外表面辛通量。切口取同一光滑全局解两侧的限制，诱导几何及面积对全部小参数精确相同。每侧的面积一般在二阶变化；本轮没有证明一个“各侧面积固定不变”的非线性解族。排除的是**两侧面积匹配蕴含单侧无通量**，不是“所有定面积边界条件均失败”。

紧支撑AF解的初始切向是截断波形，其后不等于式(3)的全空间平面波。式(11)在t=0已经足以检验局部可积性，无须借用错误的未来波形。

## 6. 同一结果必须保持两框架的完整混合项

在真空点保持δg_E=0，Jordan侧不能同时冻结法向外曲率。记f=−φ★/3，v_E=n_E·∂δφ，则

$$
\delta h_J=-\frac{f\cdot\delta\phi}{F_\star}h_J,\qquad
v_J=\sqrt{F_\star}v_E,\qquad
\delta K^{\rm ext}_{J,ab}=-\frac{h_{J,ab}}{2F_\star}f\cdot v_J,\qquad
\delta\Pi_J=0,\quad \delta\pi_J=-K(\phi_\star)v_E .
\tag{12}
$$

这里h_E=diag(−1,1,1)、√|h_E|=1。最后两式由618原P=F(Kh⁻¹−K^up)+h⁻¹f·v、j=fK−v直接得到。故原Jordan与Einstein完整配对相同。

对较重原模，k=.7、t=0，正确局部通量为.7。只冻结Jordan外曲率、继续把该配对当作同一Einstein固定几何变化，会得到.656845674249，漏.043154325751。这个比较检验的是不一致冻结的错误，不是两个正确框架的物理差异。

## 7. 对区域生成元的限定结论

令Θ为原作用标准协变辛势、ω=δΘ，ξ为字段无关的边界切向时间流，S为固定空间切片的切口。对在壳场及线性化解，标准未扩展区域的候选生成元一形式及其curl为

$$
\alpha_\xi[\delta]=\int_S(\delta Q_\xi-i_\xi\Theta),\qquad
(\delta\alpha_\xi)(\delta_1,\delta_2)=-\int_S i_\xi\omega(\delta_1,\delta_2),\qquad
\big|\delta\alpha_{\partial_t}(\delta_1,\delta_2)\big|=\mathcal F>0.
\tag{13}
$$

两参数的混合变化可交换，δ₁δ₂Q相消；由此直接得到第二式。整体符号依S与向外法向的取向，非零结论不依符号。Gibbons—Hawking等纯场空间全微分不改变本组δΘ；纯几何角点修正也因两条δg_E=0而不抵消物质项。对切向紧支撑和本固定切口，这与式(4)的辛配对相同。

因此不存在局部函数H_ξ使δH_ξ=α_ξ并同时保留上述全部变化。给H增加一个普通函数的精确微分不能消去curl。这并不是禁止更换边界辛结构、加入额外边界自由度、改用有通量的荷，或限制允许变化。

例如边界理论的完整替代资料一般应满足

$$
(\Theta+\delta\ell_B)|_B=dC_B,\qquad
\Omega_C=\int_C\omega-\int_{\partial C}\delta C_B,
\tag{14}
$$

其中允许变分、边界作用ℓ_B、角点势C_B必须一起说明；这引用成熟方法，不声称本轮已选出唯一方案。只给面积而保留任意原物质变化，并没有给出这些资料。关系嵌入的字段依赖还会改变ξ和变分，式(13)不能原样跳用于647。

## 8. 完整透射拼接保留整体演化

继承618同一正F、光滑类时无源界面：两侧h、φ识别，相对各自向外法向要求

$$
\Pi_+^{ab}+\Pi_-^{ab}=0,\qquad \pi_{+,A}+\pi_{-,A}=0
\quad\Longrightarrow\quad
\vartheta_{B,+}+\vartheta_{B,-}=0,\qquad
\omega_{B,+}+\omega_{B,-}=0
\tag{15}
$$

最后一步是在共同匹配子空间上取外微分，包含动量匹配的切向变化；不能只令某个基点上的动量相消。原F+3|f|²/2=2使混合块可逆，原法向一阶数据与这些动量相互确定。因此这套拼接不是仅有面积的重复条件。

沿同一无源光滑解、固定切口和匹配的变分，区域辛形式可相加且内部通量取消：

$$
\Omega_{C_+\cup C_-}=\Omega_{C_+}+\Omega_{C_-},\qquad
\Delta\Omega_{C_+}=-\int_B\omega_+,\qquad
\Delta\Omega_{C_-}=-\int_B\omega_-,\qquad
\Delta(\Omega_{C_+}+\Omega_{C_-})=0
\tag{16}
$$

此处只展示内部界面贡献，外部通量须另按边界条件处理。整体可保持辛结构，单个区域仍可交换物质和信息；既没有推出独立区域Hilbert张量因子，也没有量子化完整引力约束。619的费米透射和617的固定图规范切分保留各自范围，不能凭同名“拼接”直接签收其与本式的全量子等价。

## 9. 核验、减少的自由与下一步

三组复算分别核原势模式及两框架通量；连续AF收缩证书与原非线性系数；紧支撑积分及一般非恒定边界一形式的独立差分。64、128、192阶积分给.609908597826；边界一形式在两框架的curl差小于2×10⁻¹⁰。一般边界数据差分是离壳运动学检查，真正的在壳见证由第4节构造提供。

减少的自由：不能分别任取面积匹配、允许物质变化和区域Hamiltonian，而不检查同一边界辛通量。完整透射识别使物质来源、几何动量及区域组合共同受约束。仍有竞争路线：开放区域及通量账、限制允许变化的边界条件、扩展边界相空间。尚未由认知原则唯一选择。

[651入口](round651_drafts/STATUS.md)优先接647的字段依赖关系区域：同一原物质参考定义移动切口时，体积分与边界配对必须如何共同变分；检查完整匹配是否自动保留，及退化参考是否留下真实缺口。不得只重述一般Lie导数，也不继续做控制器或波形优化。总目标保持，认知系统设计后置。
