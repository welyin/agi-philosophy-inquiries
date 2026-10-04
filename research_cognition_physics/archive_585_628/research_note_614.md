# 第614轮：原物质、全局商群与成熟手征表示的共同字典

日期：2026-10-01。接[613](research_note_613.md)及[614入口](614/drafts/STATUS.md)。[代码](614/joint_spinor_subgroup_mass.py)、[结果](614/joint_spinor_subgroup_mass_results.json)、[核验](614/research_round_614_checks.json)、[条件账](614/unified_physics_condition_ledger_614.md)。三组复算、十二式；主代理审查，无新增独立代理审查。

## 1. 正面连接与需要分开的两种用途

613排除了用非幂等手征筛选直接实现原规范群。本轮回到一个合法的有限维群表示，给原物质到成熟手征候选的**实际字典**：

- 原一代16个内部Weyl分量（含ν）、原整数超荷、原Z₆商，都可在同一个Spin(10)旋量表示的子群限制中保持。
- 原32模式CAR质量，通过固定的粒子—空穴变换与内部基变换，保全部原复Yukawa、Majorana、BdG谱及共同标量来源。没有补成64个物理CAR模式，也没有为嵌入重拟合原Y。
- **采用表示容器，不等于规范化整个容器。** 若真的扩大为完整Spin(10)规范群，同时仍把原实singlet作为全群不变场，则原非零ν Majorana项不合法。这条额外要求有明确Ward身份反例；只保原子群时没有此障碍。

因此可以对接相应成熟候选而暂不增加完整大统一规范群与Higgs多重态。连接仍只到内部表示与有限质量过程，没有完成真实格点Weyl测度、局域性、反射正性或GR。

## 2. 文献中本来就有子群路线，不能遗漏

[Kikukawa，1710.11618，§5.2—5.3](https://arxiv.org/html/1710.11618v3)明确讨论把链路限制到子群，以及标准模型加右手中微子与Higgs／Yukawa的接入。其§3.6保留测度对规范场的光滑／局域性问题。该候选使用右手辅助部门的饱和构造，不是613朴素投影，也不因“规范不变”便自动提供全部Hamiltonian性质。

本轮复用这条成熟方向，补本项目原数据的映射。Spin(10)旋量分解与Yukawa表示限制是成熟群论；[Babu、Bajc、Saad，Yukawa sector of minimal SO(10) unification，§1](https://inspirehep.net/files/007e59be20072a6ca30b0d3835a12f58)给出相应费米双线性分解。本轮不把这些文献结果算作新物种推导。

辅助积分变量仍需在测度与重建中交代，不能从本内部字典推出其无物理代价。2026报告大PDF未取得可信正文，本轮不依赖它。

## 3. 一份表示同时保荷、颜色、弱同位旋和Z₆商

使用全左手约定。原598物理右手u、d、e、ν转为其电荷共轭；Q、L保持。原群与表示为

$$
G=\frac{SU(3)\times SU(2)\times U(1)}{\mathbb Z_6},\qquad
\mathcal R_L=(3,2)_1\oplus(\bar3,1)_{-4}
\oplus(\bar3,1)_2\oplus(1,2)_{-3}
\oplus(1,1)_6\oplus(1,1)_0 .
\tag{1}
$$

下标是原整数6Y。定义到SU(5)的映射

$$
\iota(C,W,z)=\operatorname{diag}(z^{-2}C,z^3W)\in SU(5),\qquad
\ker\iota=
\left\langle(e^{2\pi i/3}I_3,-I_2,e^{i\pi/3})\right\rangle
=\mathbb Z_6 .
\tag{2}
$$

证明：若像为I₅，则C=z²I₃、W=z⁻³I₂，行列式条件迫使z⁶=1。反之六个中心元素均在核内，故保留的是原全局商，而不只是相同Lie代数。

在内部空间V=C⁵上取偶外代数

$$
S=\Lambda^{\rm even}V=\Lambda^0V\oplus\Lambda^2V\oplus\Lambda^4V,
\qquad \dim_\mathbb C S=1+10+5=16 .
\tag{3}
$$

这是Spin(10)一个半旋量的SU(5)限制。它是**内部表示维数**，不是时空维数。V中的前三个基向量荷−2，后两个荷3；外积荷相加，得到如下显式基字典。令a,b,c∈{1,2,3}为颜色，α∈{1,2}为弱指标，★为按固定体积形式定义的补指标基，ε=iσ₂：

$$
\begin{aligned}
Q_{a\alpha}&\mapsto e_a\wedge e_{3+\alpha},&
u^c_a&\mapsto\tfrac12\epsilon_{abc}e_b\wedge e_c,&
e^c&\mapsto e_4\wedge e_5,\\
d^c_a&\mapsto\star e_a,&
L_\alpha&\mapsto\sum_\beta(\star e_{3+\beta})\varepsilon_{\beta\alpha},&
\nu^c&\mapsto1 .
\end{aligned}
\tag{4}
$$

这里★仅在指定实正交基上定义补指标的线性对应；不是误把一般复向量上的反线性Hodge算符当成复线性等距。SU(2)的ε将共轭双重态与原双重态交织。按代码固定符号，形成酉矩阵J，满足

$$
\Lambda^{\rm even}\!\iota(g)\,J=J\,R_L(g).
\tag{5}
$$

各外代数模块直接验证该身份；数值对18个非Abelian变换独立用子矩阵行列式计算，最大误差2.21×10⁻¹⁵。由于Λ⁴包含忠实的SU(5)共轭基本表示，(3)不会再引入额外群核。原Higgs也满足(0,H)↦(0,z³WH)，其四实分量子空间在原子群下保持；原实s继续不变。

这减少的是分别指定荷表、商群和另一份表示字典的自由，不证明为何自然界选择该群、该表示或三代。

## 4. 原有限CAR也能使用同一全左手字典

(1)的电荷共轭不能在原Hamiltonian上只改标签。对原右手两分量模式定义b_R=εa_R†，左手b_L=a_L。在列向量约定中

$$
b=Ua+Va^\dagger,\qquad
U=P_{\rm left},\quad
V=\bigoplus_{A\in{\rm right}}I_{R_A}\otimes\varepsilon,\qquad
\mathcal W=
\begin{pmatrix}U&V\\V^*&U^*\end{pmatrix},
\qquad \mathcal W\mathcal W^\dagger=I .
\tag{6}
$$

未写的块为零。左右支持正交，且UU†+VV†=I、UVᵀ+VUᵀ=0，故保CAR。有限模式的Bogoliubov变换可由Fock酉实现；此处没有改变模式数。全左手内部群作用正是(1)，因为ε只作用spin，右手内部表示变成复共轭。

令原598 BdG矩阵为B(φ)，再令A=J⊗I₂、K=diag(A,A*)W，则

$$
B(\phi)=
\begin{pmatrix}h(\phi)&\Delta(\phi)\\
\Delta^\dagger(\phi)&-h^T(\phi)\end{pmatrix},
\qquad
B_{\rm carrier}(\phi)=\mathcal K B(\phi)\mathcal K^\dagger,\qquad
\operatorname{spec}B_{\rm carrier}=\operatorname{spec}B .
\tag{7}
$$

原正常质量h仅连接左右模块；变换后它们成为全左手的配对质量。原ν项也同步变换。以配对产生项的内部对称系数M记号，有

$$
\mathcal W B\mathcal W^\dagger=
\begin{pmatrix}0&M\otimes\varepsilon\\
(M\otimes\varepsilon)^\dagger&0\end{pmatrix},
\qquad M=M^T,\qquad
M_{L,R^c}=-y_{LR}(\phi),\quad
M_{\nu^c\nu^c}=-Y_s^*\frac{s}{\sqrt{F(\phi)}} .
\tag{8}
$$

符号随已写的ε约定固定，并非新的Y参数。原h和新正常块的迹均为零，因此这个质量部门的正规序常数没有遗漏；若将来变换含对角跳跃的完整H，须重新保留其常数与几何来源，不能照搬“零”。

原32×32正常／配对矩阵、64×64 BdG谱和原全部复Y均进入复算。诊断点的中性系数为−.1142157996+.0331594257i，与(8)相符。随机12个原标量背景中谱误差≤1.12×10⁻¹⁵，块身份误差为零。

## 5. 质量、规范与标量来源一起保留

原质量的规范协变通过(5)—(7)继承到容器表示，代码也独立变换原Higgs再调用原mass_matrices验证。J、W均不依赖φ或几何参数，所以

$$
\partial_{\phi_a}B_{\rm carrier}
=\mathcal K(\partial_{\phi_a}B)\mathcal K^\dagger,\qquad
\partial_\lambda B_{\rm carrier}
=\mathcal K(\partial_\lambda B)\mathcal K^\dagger
\quad\text{若}\quad\partial_\lambda\mathcal K=0 .
\tag{9}
$$

来源不能另外重选；五个原标量方向的独立差分核验通过。这与611在变化规范背景中的移动手征纤维不同：本轮只做固定内部／粒子空穴字典，没有把它冒称为一般GW场依赖投影的完整过程交织。

到此证明的是可用的**原群—原物种—原有限质量与来源接口**。没有证明此有限CAR就是文献Euclidean Weyl测度所重建的Hilbert空间，也没有在本轮添加真实手征空间动能。

## 6. 若额外要求整个Spin(10)都成为物理规范群，会发生什么

在五个内部振子的偶子空间中，令N为占据数。Spin(10)含额外Cartan生成元X=2N−5；它在Λ⁰、Λ²、Λ⁴上分别为−5、−1、3。把ν^c对应的配对系数方向记为B₀=|0〉〈0|，则

$$
X B_0+B_0X^T=-10B_0\ne0 .
\tag{10}
$$

因此非零sν^cν^c项若仍使用全Spin(10)不变的实s，不满足额外生成元的Ward身份。原子群却对ν^c作用平凡，故没有该障碍。这里的X是规范生成元候选，不是额外时空坐标。

还可以确定缺了多少生成方向。Spin(10) Lie代数可用数守恒u(5)与十个复配对方向表示。令T B₀+B₀Tᵀ=0，等价于T|0〉=0：配对产生部分必须为零，u(5)迹也必须为零。因此

$$
\dim_{\mathbb R}\ker
\{T\mapsto TB_0+B_0T^T\}=24,\qquad
45-24=21 .
\tag{11}
$$

保留下来的是su(5) Lie代数；没有据此声称已分类全部离散稳定子。代码从10个Clifford矩阵生成45个Hermitian生成元，独立核实实线性约束秩21和24个su(5)方向的零残差。

若坚持未破缺的全Spin(10)、线性基本标量和重整化型双费米—单标量耦合，ν^cν^c极端权的轨道生成一个126维复表示。用D₅的ρ=(4,3,2,1,0)、λ=(1,1,1,1,−1)，成熟Weyl维数公式给

$$
\dim V_\lambda=
\prod_{1\le i<j\le5}
\frac{(\lambda_i+\rho_i)^2-(\lambda_j+\rho_j)^2}
{\rho_i^2-\rho_j^2}
=126 .
\tag{12}
$$

本例是半旋量极端权的对称平方，其最高权为两倍旋量最高权；共轭约定改变的是126或其共轭名称，不改变维数。公式由有理数精确计算。对应标量需有对偶多重态，而非仍只有原一个实s。复合算符、先破缺再取有效理论、其它中微子机制均是不同分支，本限制不排除它们。

这一额外多重态要求**只属于“物理规范群扩大”分支**。采用(2)的子群限制去利用成熟手征候选，不要求加入126个分量，也不要求新增大统一规范玻色子。

## 7. 联合条件得到怎样的压缩

三组检查分别覆盖：实际子群交织／Z₆／Higgs；原完整有限CAR质量／BdG／来源；完整容器规范化的Majorana障碍与精确维数。它们全部通过，但“检查通过”在第三组意味着确认指定扩大接法失败。

由此可以共用一份字典处理C14群、C15物种、C17质量、C18原参数和C22来源。我们无需为对接该候选先增加物理Spin(10)假设，也无需将原Y改成未经论证的大统一质量关系。这是条件整合的实际减少输入机会，尚非唯一候选或完整模型存在证明。

剩下的连接是实质性的：同一允许配置域上的Weyl测度及其规范变化、局域性、真实时间重建、与原玻色作用／热态／来源的匹配。文献的辅助部门饱和、真正物理右手ν，以及本轮用于字典的内部五振子，属于三个不同对象，不能混叫新增或已删除粒子。

[615入口](615/drafts/STATUS.md)继续核**限制到原子群的测度与共同来源**。先取成熟定理或候选中已明确的条件，检查哪些在限制映射下严格继承、哪些需要额外积分或反项；不重新做反常荷表计数，不宣布2017的开放问题在本项目中已被解决。完整量子动力学、时空及GR仍开放，认知系统设计后置。

