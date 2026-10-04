# 第664轮：原物质、联络选择与共同几何来源的匹配条件

日期：2026-10-02。接[663](research_note_663.md)、[664入口](664/drafts/STATUS.md)。[代码](664/joint_connection_matter_matching.py)、[结果](664/joint_connection_matter_matching_results.json)、[核验](664/research_round_664_checks.json)、[条件账](664/unified_physics_condition_ledger_664.md)。三组复算、十六式；主代理推导与审查，无独立代理审查。

## 1. 本轮整合的条件

**原来的无挠度规模型，与在原Jordan作用中把spin联络独立化的最小模型，不能直接认作同一个物质理论。** 对项目原五标量和整代费米内容，后一操作同时改变标量目标度规并生成总轴流的接触作用。只增加自旋来源而继续原样使用旧标量几何和Gaussian费米权重，会混合不同分支。

这不迫使选择Einstein–Cartan理论。若目的只是以独立联络重写原模型，可以给出成对的匹配项，在经典体内消元后恢复原作用及其共同来源。该等价尚不包含量子测度、边界或独立挠率观测。

上一条用户问答仅确认研究顺序，没有新增结果，按无进展计。本次先读导航、663报告／结果／事后核验及664入口，实际进程查询未找到Python进程，不假设旧计算仍在运行。相关历史回查如下：

|历史|继承的结果与范围|
|---|---|
|[358](../archive_342_369/research_note_358.md)|无挠独立仿射联络的真空路线；原报告明确排除未处理的自旋／挠率分支|
|[375](../archive_370_428/research_note_375.md)|固定Levi-Civita的连接分类，以及有挠传播的不可辨识核；不重算分类|
|[580](../archive_554_584/research_note_580.md)|原五标量目标为H⁵，截面曲率−1/6；已有目标几何及圈项依赖这一输入|
|[598](../archive_585_628/research_note_598.md)—[600](../archive_585_628/research_note_600.md)|实际32CAR、完整质量及度规／标量来源；连续作用使用无挠spin连接|
|[619](../archive_585_628/research_note_619.md)|原左右手表示、Majorana计数、Weyl权重和区域透射|
|[643](../archive_629_652/research_note_643.md)、[663](research_note_663.md)|条件Gaussian费米权重；给定无挠几何下的旋量输送与完整字典导数|

本轮区分认知动机与物理输入：认知动机是同一物质与几何不能采用互不相容的更新规则；额外输入是四维Lorentz／spin背景、指定Einstein领先作用、Jordan最小独立联络规则及有限单元的正常排序。解析结论负责一般局部关系；数值只核原模型系数与明示条件扇区。三维选择、引力作用选择、连续量子极限均未解决。

## 2. 固定原模型及比较约定

原五个实场由复Higgs和实singlet组成，保留原势V、规范作用及所有Yukawa／Majorana参数。用M_P²=1单位，原F、Einstein目标和质量为

$$
F(\phi)=M-\frac{|\phi|^2}{6}>0,\quad M=2,\qquad
K^{\rm met}_{ab}=\frac{\delta_{ab}}F+\frac{\phi_a\phi_b}{6F^2},\qquad
\mathcal M_E=\frac{\mathcal N(\phi)}{\sqrt F}.
\tag{1}
$$

本轮Lorentz符号η=diag(−,+,+,+)，标量作用写为−K(Dφ)²/2；600的Euclidean正动能约定不直接照搬符号。连续M分支的Jordan体内作用写为

$$
S_{\rm met,J}=\int\sqrt{-g_J}\left[
\frac F2 R(g_J)-\frac12\delta_{ab}D\phi^a\cdot D\phi^b-V
\right]+S_{\rm gauge}+S_{\rm ferm}[\omega_{\rm LC},\phi].
\tag{2}
$$

这是原模型的给定连续表示，不是从图H新推出的引力。比较分支C_J把式(2)曲率和Hermitian费米主部中的同一spin联络改为独立、度规相容联络；不加Holst项、非最小挠率项或新的标量流耦合。规范曲率仍用外微分和内部规范联络，不改为带挠率的时空协变反对称导数。

成熟接口采用[Karananas等，*Matter matters in Einstein–Cartan gravity*，式(21)、(28)、(34)、(46)—(50)](https://arxiv.org/html/2106.13811)：挠率分解为向量、轴向量和无迹张量，可在该阶代数消元。这里只取明示最小分支，未接受其全部非最小参数。下面代入原F及原物种，独立完成匹配。

## 3. 原F使标量与费米部门一起改变

记向量挠率v、轴挠率a、无迹张量挠率τ。令A为文献的轴流约定；它与下文J₅只差整体负号。先处理紧支撑变分或配套边界项，将F的散度项分部积分后，联络依赖密度为

$$
\mathcal L_{\rm tor}
=-\frac F3v_\mu v^\mu+\frac F{48}a_\mu a^\mu
 +\frac F4\tau_{\mu\nu\rho}\tau^{\mu\nu\rho}
-v^\mu\partial_\mu F-\frac18a_\mu A^\mu.
\tag{3}
$$

直接变分得到v=−3dF/(2F)、a=3A/F、τ=0。代回即为

$$
\Delta\mathcal L_J
=\frac{3}{4F}(\partial F)^2-\frac{3}{16F}J_{5,J}^2,
\qquad
g_E=Fg_J,\quad \psi_E=F^{-3/4}\psi_J,\quad
J_{5,J}^{\hat a}=F^{3/2}J_{5,E}^{\hat a}.
\tag{4}
$$

帽指标为正交标架指标。原Majorana质量不含联络，既不改变消元方程，也不使Weyl物种加倍。原规范及Dirac质量仍按同一619字典输送。

式(4)第一项在Einstein表示恰抵消原度规Weyl变换产生的径向动能。轴流项的F次幂同时抵消：体积F⁻²、系数F⁻¹、轴流平方F³。于是

$$
K^{C_J}_{ab}=\frac{\delta_{ab}}F,\qquad
K^{\rm met}-K^{C_J}=\frac{3}{2F^2}\,dF\otimes dF,
\qquad
\Delta\mathcal L_{E,4\psi}=-\frac3{16}J_{5,E}^2.
\tag{5}
$$

这不是修改Y后可以忽略的差异。若改在已经得到的Einstein作用中才实施“最小独立联络”，则标量K保持原值、仍增加轴流接触项，得到另一个C_E分支。“最小化”必须注明在哪个表示和哪份作用中实施；独立联络本身不能唯一决定所有系数。

## 4. 两个标量目标不是单纯换坐标

在相同五维场域、相同正曲率球的Riemann约定下，原目标标量曲率为−10/3。对K^{C_J}=e^{2u}δ取u=−lnF/2，用共形曲率公式直接求得

$$
\mathscr R[K^{\rm met}]= -\frac{10}{3},\qquad
\mathscr R[K^{C_J}]= -\frac{20}{3}-\frac{7|\phi|^2}{9F}.
\tag{6}
$$

一个是常数，另一个变化且始终更小。因而在保留Einstein度规表示的普通标量坐标重定义下，两目标不局部等距。此处不排除连时空度规及其它耦合一起改变的更广义场重定义；那需要重新匹配完整作用。

径向r从0向√(6M)移动，直接积分径向线元还给

$$
\ell_{\rm met}(r)=\sqrt6\,\operatorname{artanh}\frac r{\sqrt{6M}},\qquad
\ell_{C_J}(r)=\sqrt6\,\arcsin\frac r{\sqrt{6M}}.
\tag{7}
$$

原F=0边界距离无穷，C_J边界距离为π√6/2≈3.84765。因此依赖原完整双曲目标、原测度及其Laplace–Beltrami算符的旧证明，不能自动搬到C_J。有限距离不单独证明量子不稳定或没有自伴实现；势、域与边界仍须审查。

独立微分复算从两份度规本身生成Christoffel与Ricci，没有把式(6)作为数值算法。三点、两步长误差由至多6.87×10⁻⁸降至1.72×10⁻⁸；原点得到−3.333333336和−6.666666672。

## 5. 想保留原模型，需要哪些匹配项

若独立联络只是辅助表示，可以在C_J加入

$$
\mathcal L_{\rm match,J}
=-\frac{3}{4F}(\partial F)^2+\frac{3}{16F}J_{5,J}^2,
\qquad
\left(S_{C_J}+\int\sqrt{-g_J}\mathcal L_{\rm match,J}\right)_{\omega=\omega_*}
=S_{\rm met,J}.
\tag{8}
$$

这是在指定体内变量和消元阶数下的等式，匹配项不是新认知公理，也不是唯一物理选择。一般非最小联络作用会给别的接触矩阵；本轮未排除它们。由于独立联络驻定，任意原物质或几何参数λ满足

$$
\frac{\delta S_{\rm eff}}{\delta\lambda}
=\left.\frac{\delta S}{\delta\lambda}\right|_{\omega_*}
 +\left.\frac{\delta S}{\delta\omega}\right|_{\omega_*}
\frac{\delta\omega_*}{\delta\lambda}
=\left.\frac{\delta S}{\delta\lambda}\right|_{\omega_*}.
\tag{9}
$$

故标量力、度规来源不能分别拟合；必须对同一个完整作用求导。比较独立挠率读出时，式(8)不保证等价。独立联络的量子积分还可能有行列式、测度和局部匹配项，本轮不把经典驻定消元当作量子积分相等。

## 6. 共同联络耦合的是原全部物种的总流

保留598物理右手约定：Q、L手性χ=−1，u、d、e、ν为+1。每个内部二分量Weyl场的数流空间部分为χσ，因此总轴流的共同局部表达是

$$
J_5^0=\sum_A\chi_A\psi_A^\dagger\psi_A,\qquad
J_5^i=\sum_A\psi_A^\dagger\sigma_i\psi_A,\qquad
J_5^2=-(J_5^0)^2+\sum_{i=1}^3(J_5^i)^2.
\tag{10}
$$

原内部重数6、3、3、2、1、1不改，旋量指标与SU(2)_L指标分开。直接用原32×32商群表示检查四个流系数，均与原内部群对易。共同联络消元必然保留跨物种作用：

$$
\left(\sum_A J_{5,A}\right)^2
=\sum_A J_{5,A}^2+2\sum_{A<B}J_{5,A}\cdot J_{5,B}.
\tag{11}
$$

不能让每个部门独立消去一份联络后再相加；那相当于改变了共同耦合结构。也不能用各流的平均值乘积替代量子二点。

## 7. 实际CAR的接触项及原热权重

有限体积单元w中取ψ=c/√w，c为原canonical CAR。指定相对于同一空真空的正常排序，仅作为有限单元调节约定，不声称它已给连续重整化。记流系数C₀=diagχ、C_i=⊕σ_i，Ĵ_a=c†C_ac。CAR恒等式给

$$
\widehat Q=-(\widehat J_0)^2+\sum_i(\widehat J_i)^2,
\qquad \widehat Q=:\widehat Q:+c^\dagger\left(-C_0^2+\sum_iC_i^2\right)c
=:\widehat Q:+2\widehat N.
\tag{12}
$$

所以直接平方与正常排序相差确定的双线性项。按式(5)的Lorentz作用和无导数相互作用，Einstein单元Hamiltonian项为

$$
H_{\rm tor}(w)=\frac{3}{16w}:\widehat Q:,
\qquad w=\epsilon^3e^{6\theta},\qquad
\partial_\theta H_{\rm tor}=-6H_{\rm tor}.
\tag{13}
$$

诊断取ε=1、κ=1自然单位，不是现实引力强度拟合。该项在固定w有限CAR中有界，但不是原二次费米Hamiltonian；不据此证明动态w→0或连续极限存在。

代码保留全部八个轻子模式及原Dirac／Majorana矩阵。夸克空真空是条件单节点作用的一个不变扇区；它不是完整Gauss热积分。实际e_R与ν_R各单占据的旋转singlet，以及原ν双占据的对角四点差给

$$
\langle s_{e\nu}|:\widehat Q:|s_{e\nu}\rangle=-8,
\qquad
\left\langle s_{e\nu}\left|\sum_{A=L,e,\nu}:\widehat Q_A:\right|s_{e\nu}\right\rangle=0,
\qquad Q_{\uparrow\downarrow}-Q_\uparrow-Q_\downarrow+Q_0=-8.
\tag{14}
$$

末式对任何二次CAR Hamiltonian的相应对角组合都为零，故不能仅通过重新选择32×32质量／传播矩阵吸收此作用。带电见证不被称为整个封闭系统的Gauss物理态；算符自身的内部规范不变性另由完整表示检查。

## 8. 同一个条件热过程必须给同一个几何来源

在原φ=(0.45,0.31,−0.17,0.26,0.62)、θ=0.13、β=0.7处，使用原八轻子模式256维Fock矩阵，保留复Dirac及Majorana质量。完整条件热迹满足

$$
Z=\operatorname{Tr}e^{-\beta(H_0+H_{\rm tor})},\qquad
\partial_\theta\log Z=-\beta\langle\partial_\theta H_{\rm tor}\rangle
=6\beta\langle H_{\rm tor}\rangle.
\tag{15}
$$

不要求H₀与H_tor对易；等式由迹的循环性得到。直接谱计算：

|条件单元方案|log Z|
|---|---:|
|原二次质量H₀|5.559136407732|
|总流正常排序接触项|5.885124337122|
|使用未排序的流平方|5.388980254270|
|删除跨模块流作用|6.305643884585|

同一几何导数为−2.443719917196；独立中央差分为−2.443719929812，误差1.27×10⁻⁸。冻结接触项的体积系数会漏掉全部这一来源。加上同处方的相反接触项则在算符层精确恢复H₀，而不是仅拟合某个温度的Z。

## 9. 条件压缩、保留分支与下一步

本轮条件图可概括为

$$
\begin{array}{c}
\text{原Jordan作用＋最小独立共同spin联络}\\
\Downarrow\\
\text{新的标量目标＋总轴流接触项＋同源几何变分}\\
\text{再加匹配项并经典消元}\ \Downarrow\\
\text{原无挠度规／物质作用}
\end{array}
\qquad\not\Longrightarrow\quad
\text{量子测度或全部物理理论已经相同}.
\tag{16}
$$

三组核验分别覆盖原F的消元／几何不变量、实际物种流与CAR排序、原质量的条件热权重及几何来源。原模型没有被改写，旧冻结结果保持有效于其原分支。并未新增一套任意参数来拟合差异：指定最小分支后，径向差、共同接触系数及其体积来源已经相互固定。

仍有四个未完全统一的主分支：原固定图量子过程、声明的连续物质、给定作用的经典动态几何、原手征辅助测度。这里的联络选择是连续几何／物质接口中的子分支。保留原无挠方案作为已验证基线；C_J与C_E不能静默取代它。维数、群、谱、Y参数、完整区域条件、连续极限、量子引力仍独立开放。

接[665共同量子权重与联络辅助表示](665/drafts/STATUS.md)：检查原完整CAR中的接触作用能否接到同一有序热历史与来源，明确辅助积分、归一化和原Gauss投影的关系。先回查643—655的实际权重，不能复述一般Gaussian消元凑轮次，也不转入认知装置设计。
