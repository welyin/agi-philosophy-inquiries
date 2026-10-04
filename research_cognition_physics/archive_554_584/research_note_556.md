# 第556轮：Higgs不变参考的合法读取——径向域反例、共同能量界与后态约束

日期：2026-09-30。接续[555](research_note_555.md)、[554](research_note_554.md)和[548](../archive_531_553/research_note_548.md)。[代码](556/joint_higgs_invariant_reference_readout.py)、[结果](556/joint_higgs_invariant_reference_readout_results.json)、[核验](556/research_round_556_checks.json)、[条件账](556/unified_physics_condition_ledger_556.md)。早期[候选稿](556/drafts/research_note_556_draft.md)保留。

## 1. 共同接口与实际增量

554—555量化了全实轴上的两实场h、s，未证明它们构成完整Higgs规范理论的量子子部门。本轮返回**单胞元Higgs双重态的四个实分量及一个实singlet**，采用全局SU(2)不变部门，核对同一势、读取与状态更新能否相容。

得到三项相连的结论：

1. 一个同势、有限完整二阶能量矩的源，在径向表示中却不属于拆开的平直动能与逆平方乘法各自的算符域。旧两实场的图估计及Weyl读口不能直接搬过去。
2. 保留完整Cartesian动能，并将相空间输出按规范不变量合并，可以构造合法的不变POVM。一个完整H的二阶矩预算继续控制四个实际读数，允许非Gaussian不变源。
3. 同样的不变测量概率，不自动保证测量后态留在不变子空间：明确的heterodyne测量后制备仪器，在一个正规不变Gaussian输入上有**11/12**的后态权重离开该子空间。有限分箱的平方根仪器可抽象地避免这项泄漏，但其物质实现仍需另行交代。

|层次|本轮地位|
|---|---|
|认知动机|共同参考必须属于同一物质结构，读取过程也要符合声明的约束|
|继承|544的L、C、u及正四次势；555的图范数方法；554的相空间POVM|
|新增模型输入|单胞元Cartesian量子化、指定Higgs表示、全局SU(2) singlet部门；无空间梯度和规范链路|
|新增仪器输入|各向同性Gaussian种子、不变后处理；明确区分两种后态更新|
|解析结论|域反例、实际不变读数及风险界、特定仪器精确泄漏与抽象保部门替代|
|没有得到|局域Gauss理论、完整电弱场论、548时空坐标、探针来源、连续极限或引力生成|

四个Higgs实分量是给定内部表示，不能据此推断四维时空。此处“Higgs”指沿用表示与势的有限量子力学候选，并非已证明完整标准模型的受控截断。

## 2. 去重与成熟接口

[Paz，§I—II](https://arxiv.org/pdf/math-ph/0009016)给径向动量的单位变换及非自伴边界；[Dereziński–Richard，§1—2](https://www.fuw.edu.pl/~derezins/inverse-square_revised.pdf)给半直线逆平方算符的实现。本轮从Cartesian不变部门确定算符，不任意选择新边界条件。径向公式本身不是新发现。

POVM继承[Werner，§3.1式(23)—(26)](https://arxiv.org/pdf/quant-ph/0405184)。新增的是它与本势、正确不变菜单及仪器后态的共同连接，不借用文献中另一种误差定义的常数。旧[360](../archive_342_369/research_note_360.md)、[365](../archive_342_369/research_note_365.md)的区域拼接以及[421](../archive_370_428/research_note_421.md)的控制权限直接继承；本轮不重复一般Gauss或WAY障碍。

## 3. 完整对象与径向表示

取X∈R⁴、s∈R，r=|X|，v>0，a=ℏ²/(2v)。沿用正定L=[[λ_h,p],[p,λ_s]]、C=Lu及u_h,u_s>0：

$$
\mathcal H=L^2(\mathbb R^4_X\times\mathbb R_s),\qquad
H=-a(\Delta_X+\partial_s^2)+W,\quad
W=\frac v4(r^2-u_h,s^2-u_s)L(r^2-u_h,s^2-u_s)^{\mathsf T}\ge0,
\qquad \mathcal H_0=P_0\mathcal H,\quad
P_0=\int_{SU(2)}U_g\,dg.
\tag{1}
$$

SU(2)以双重态作用在X上，在每个非零r的S³上传递。因此不变波函数只依赖r、s。H与U_g强对易，故H₀=H|𝓗₀是自伴的约化部门；全局对称性与这个部门的选择均已声明。吸收球面积后，原内积测度为r³drds，酉映射u_rad=r^(3/2)ψ给

$$
H_{\rm rad}=-a\partial_r^2-a\partial_s^2+\frac{3a}{4r^2}+W(r,s),
\qquad
T_R=-a\partial_r^2+\frac{3a}{4r^2}
=\left.(-a\Delta_X)\right|_{\mathcal H_0}.
\tag{2}
$$

右侧等号包含相应酉映射。这里的T_R应整体保留；它不是单独的平直半直线动能。形式径向动量变换为−iℏ∂r，不能直接继承全实轴的自伴平移群及全部Weyl位移权限。

## 4. 同一有限能量源的域反例

取归一不变Gaussian ψ=C exp[−ω(r²+s²)/2]，ω>0。它是参考振子的真空，不是H的相互作用真空。Cartesian Hψ为多项式乘Gaussian，因此E₂=‖Hψ‖²有限。然而

$$
\frac{\partial_r^2u_{\rm rad}}{u_{\rm rad}}
=\frac3{4r^2}-4\omega+\omega^2r^2,\qquad
T_Ru_{\rm rad}=a(4\omega-\omega^2r^2)u_{\rm rad},\qquad
\int_{r>\delta}\!\left|\frac{3a}{4r^2}u_{\rm rad}\right|^2dr\,ds
=\frac{9a^2\omega^2}{8}\log\frac1\delta+O(1).
\tag{3}
$$

‖−a∂r²u_rad‖²有相同发散系数，两项组合的平方范数却为6a²ω²。故u_rad∈D(H_rad)，但不在平直Dirichlet动能及逆平方乘法的各自算符域。其平直动能第一矩可有限，第二矩却无限。

同一冻结势、a=1/2、ω=1.3下，E₂=3.81028584427322，发散系数0.4753125，完整T_R的第二矩2.535。对数坐标独立积分残差6.22×10⁻¹⁵。δ只用于诊断原源的范数，不改变边界、源或H。

这限制直接搬用555，而不反驳555：奇异势3a/(4r²)+W既非光滑全空间势，其二阶导数也不能由有限AW+B界住。若把平直径向平方动量当作独立可测正算符，其未裁剪第二矩在本源上无限；不能据此否定完整T_R或所有径向测量。

## 5. 同一Cartesian能量预算仍可用

令ℓ=λ_min(L)>0，β=max(0,6λ_h+p,4p+3λ_s)，η>0。五变量Laplacian及其逐点界是

$$
\Delta_{X,s}W
=v[(6\lambda_h+p)r^2+(4p+3\lambda_s)s^2-4C_h-C_s]
\le A W+B,\qquad
A=\frac{2\beta}{\ell\eta},\quad
B=v\beta(u_h+u_s+\eta)+v(4|C_h|+|C_s|).
\tag{4}
$$

系数不同于两实场的3λ_h+p，不能直接照抄555。用r²+s²≤u_h+u_s+η+[(r²−u_h)²+(s²−u_s)²]/(2η)，复用555的全空间自伴／图域证明。对任意正常源ρ，E₂=TrρH²有限时，

$$
\begin{aligned}
E_1&=\operatorname{Tr}\rho H,\qquad
\mathcal G=E_2+aAE_1+aB,\qquad
\langle T^2+W^2\rangle\le\mathcal G,\\
\langle r^2+s^2\rangle&\le S:=u_h+u_s+\eta+\frac{2E_1}{v\ell\eta},\\
\langle r^4\rangle&\le2u_h^2+\frac{8E_1}{v\ell},\qquad
\langle(P_X^2)^2\rangle,\ \langle p_s^4\rangle\le4v^2\mathcal G .
\end{aligned}
\tag{5}
$$

T=−a(Δ_X+∂s²)，P_X²=−ℏ²Δ_X。最后一步使用T_X、T_s的强对易及非负性，不使用非对易平方序。所有平方期望仍按二次型范数理解。将源限制到𝓗₀后界不变；其正确径向像是完整T_R。一个E₂≤𝓔²预算足够，因为E₁≤𝓔。

## 6. 不变POVM与实际矩

在完整五模空间选择独立零均值Gaussian种子τ，协方差T_τ=diag(t_qI₅,t_pI₅)，t_qt_p≥ℏ²/4。它不选择Higgs内部方向。Cartesian相空间POVM记M(dz)=D(z)τD(z)†dz/(2πℏ)⁵。定义输出映射

$$
f(z)=(|X_{\rm out}|^2,s_{\rm out},|P_{X,{\rm out}}|^2,p_{s,{\rm out}}^2),\qquad
E(B)=\int_{f^{-1}(B)}M(dz),\qquad
[E(B),U_g]=0,\qquad
E_0(B)=P_0E(B)P_0,\quad E_0(\mathbb R^4)=I_{\mathcal H_0}.
\tag{6}
$$

证明：U_gD(z)U_g†=D(R_gz)，τ不变，f(R_gz)=f(z)，相空间测度保持；换元即得对易性，压缩保正性与归一性。这构造的是约束子空间上的概率POVM，不需要半直线位移。种子的相位／位置动量尺度及仪器本身仍是资源输入。

只扣除已知种子偏置，记Y=(|X_out|²−4t_q,s_out,|P_X,out|²−4t_p,p_s,out²−t_p)。其均值恰为(r²,s,P_X²,p_s²)的量子期望，源偏置没有被删除。对任意源有

$$
\begin{aligned}
\mathbb EY_1^2&=\langle r^4\rangle+4t_q\langle r^2\rangle+8t_q^2,&
\mathbb EY_2^2&=\langle s^2\rangle+t_q,\\
\mathbb EY_3^2&=\langle(P_X^2)^2\rangle+4t_p\langle P_X^2\rangle+8t_p^2,&
\mathbb EY_4^2&=\langle p_s^4\rangle+4t_p\langle p_s^2\rangle+2t_p^2 .
\end{aligned}
\tag{7}
$$

这是555任意源卷积矩公式在纯位置／纯动量块的应用，相关Moyal常数为零；不把非Gaussian源Wigner函数当正概率。核上多项式恒等式再由式(5)的图范数控制延拓到E₂有限的混态。得到

$$
\begin{aligned}
B_1&=2u_h^2+\frac{8E_1}{v\ell}+4t_qS+8t_q^2,&
B_2&=S+t_q,\\
B_3&=4v^2\mathcal G+8vt_pE_1+8t_p^2,&
B_4&=4v^2\mathcal G+8vt_pE_1+2t_p^2,\qquad
\mathbb EY_i^2\le B_i .
\end{aligned}
\tag{8}
$$

本轮菜单与548的完整时空导数不变量不是同一个对象：这里只有单胞元，未含空间梯度；P_X²包括完整内部动能，也不等于不加修正的径向canonical动量平方。因此签收的是物质—概率—资源的单胞元接口，不是548坐标图已量子实现。

## 7. 同一概率并不确定守约束后态

取纯最小Gaussian种子，以同频参考振子coherent态|α〉记五模POVM。一个具体兼容仪器是测量后制备：

$$
\mathcal I(d^{10}\alpha)(\rho)
=|\alpha\rangle\langle\alpha|\,
\langle\alpha|\rho|\alpha\rangle\,\frac{d^{10}\alpha}{\pi^5}.
\tag{9}
$$

其效应正是该coherent POVM，粗粒化记录后也正是式(6)。但输入参考Gaussian真空属于𝓗₀，非选择Higgs四模输出却为均值占据数1的四份乘积热态。这里“热”仅指该仪器输出的精确振子分布，不是相互作用H的Gibbs态。

对总振子数n，每个Fock基态权重为2^(−n−4)。SU(2)实四维表示的复化为两个基本双重态，Symⁿ(2⊕2)=⊕_{j=0}ⁿSymʲ(2)⊗Symⁿ⁻ʲ(2)。两自旋分别j/2和(n−j)/2，只有n偶且j=n/2时各给一个singlet。等价地，径向振子每个偶数激发恰有一个径向态。于是

$$
\rho_{{\rm out},X}=2^{-(N_X+4)},\qquad
\operatorname{Tr}(P_0\rho_{\rm out})
=\sum_{k=0}^{\infty}2^{-(2k+4)}=\frac1{12},\qquad
\operatorname{Tr}[(I-P_0)\rho_{\rm out}]=\frac{11}{12}.
\tag{10}
$$

s部门对投影无影响。丢弃方向记录不消除这项非选择泄漏。该输出密度仍与SU(2)对易；**密度的对称性不等于支持在singlet子空间**，所以没有违反仪器协变性。

代码以n=0至8的完整Fock部门独立构造三个生成元与Casimir，核到重数1,0,1,0,…,1，所加singlet权重341/4096；未计singlet尾严格为1/12288，两者和为1/12。该尾不是全热态的总截断误差。中心−I的奇偶检验另给40/81泄漏下界，与精确值相容。

这个反例只排除把原测量后制备仪器直接称为保𝓗₀操作，不是规范不变测量普遍无法实现，也不是物理SU(2)被破坏。

## 8. 抽象保部门仪器与未消除的实现条件

把输出分成有限个箱及溢出类别，E_b=E(B_b)。因为E_b有界、非负且与U_g对易，其平方根亦然。可定义另一种具有相同分箱概率的仪器与中性记忆等距：

$$
K_b=\sqrt{E_b},\qquad
\mathcal J_b(\rho)=K_b\rho K_b,\qquad
V_{\rm rec}\psi=\sum_bK_b\psi\otimes|b\rangle,\qquad
V_{\rm rec}^\dagger V_{\rm rec}=I,\qquad
V_{\rm rec}P_0=(P_0\otimes I)V_{\rm rec}.
\tag{11}
$$

所以输入及任意被动参考在𝓗₀时，所有分支都保留这一部门。这证明抽象相容仪器存在；并不说明同一四次势、给定物质和有限资源已能产生K_b、完成开关或保存记录。它也不意味着未知输入态不受扰动，或测后E₂预算继续保持。全局对称的抽象等距不能直接改称局域Gauss守恒装置。

在声明的读数单位及预设常数基准c_i、非零标定j_i下，可沿555给一次有限记录保证：

$$
R_i=\frac{(\sqrt{B_i}+|c_i|)^2}{j_i^2},\qquad
r_0=\sqrt{\frac{\sum_iR_i}{\delta}},\qquad
\Pr\!\left\{\max_i\left|\frac{Y_i-c_i}{j_i}\right|>r_0\right\}\le\delta .
\tag{12}
$$

分箱可继续加有限量化误差。此保证只涉及一次末读的输入统计；固定预算下未测H演化保持E₁、E₂，不能据此签收重复测量后的来源预算或持续坐标精度。

## 9. 复算结果

源采用SU(2)不变的径向Laguerre激发k与singlet振子激发n_s，ω=1.3。前者不是逐分量独立Fock源；k>0时是相关非Gaussian源。直接在五个Cartesian变量上计算H及菜单，另以十个相空间变量上的正Bargmann输出密度积分独立核实际记录矩。

|k，n_s，v|E₁|E₂|实际P_X²读口二阶矩|共同上界B₃|
|---|---:|---:|---:|---:|
|0，0，1|1.790042|3.810286|20.280000|123.598711|
|1，1，1|4.503825|21.532717|57.460000|343.476958|
|2，0，0.7|7.211630|60.007265|114.920000|401.425116|

六组检查覆盖径向域反例、完整五变量势与能量界、独立正POVM积分、SU(2)作用与不变菜单、Fock Casimir及精确泄漏尾、统一记录预算。独立实际矩最大差2.93×10⁻¹²。一般源、全部激发重数及域结论由解析证明承担，数值未替代全族论证。

## 10. 条件账与下一项

本轮在同一势下接通C17物质、C19来源、C22能量及C03／C21概率读取，并把“保约束状态更新”从效应对称中明确分离。原canonical径向读口的直接映射已被具体反例排除，完整不变菜单分支则保留。

后继应回到共同模型的空间连接：单胞元的不变POVM尚不能替代区域间的参考读取，需要明确有限图的规范链路、相邻物质交互和实际仪器。先回查360—365与548，列出哪些现有结构可直接合并，再检验同一链路是否同时承担协变梯度、通信和边界拼接；不继续只优化单胞元的泄漏常数或分箱精度。维数、连续时空、完整引力及整体统一仍开放。
