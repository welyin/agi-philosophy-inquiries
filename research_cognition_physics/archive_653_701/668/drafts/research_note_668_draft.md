# 第668轮：原非零质量、物理观测与辅助测度的共同反射结构

日期：2026-10-02。接[667](../../research_note_667.md)、[实际入口](mass_observation_entry.md)。[代码](../joint_mass_auxiliary_reflection.py)、[结果](../joint_mass_auxiliary_reflection_results.json)、[核验](../research_round_668_checks.json)、[条件账](../unified_physics_condition_ledger_668.md)。两组新核验、十六式；主代理审查，无独立代理审查。

## 1. 问题、继承与范围

上一目标轮完成667并执行668入口，属于进展。本次核对三份导航、最新报告、结果和发布后核验，无运行中Python。660—661已有自由物理Weyl观测与原S⁹辅助测度的共同正泛函，但尚未把原非零质量加入该对象。本轮关闭这个**固定原标量背景、单位规范链路、有限盒**的接口，并核对质量、观测字典及来源的相容性。

|层次|本轮地位|
|---|---|
|认知动机|同一物质的质量、观测、概率归一和响应不可来自不同过程|
|继承输入|原598质量，614全左手及外幂载体，660变量变换，661物理代数，658辅助严格归一|
|新增声明|在物理Weyl代数中加入下述Euclidean质量二次型；原φ冻结且空间时间常数|
|几何输入|给定有限周期空间、偶数反周期时间、m₀=1的自由overlap核；维数没有由此推出|
|解析结果|准确质量拉回、全多项式反射正性、任意实质量强度的严格归一、共同来源公式|
|数值验证|全部原内部质量的有限Pfaffian、Fourier积及来源；不是全部玻色／Gauss积分|
|仍开放|原动态标量、完整Hamiltonian识别、一般规范背景、实际连续极限及量子GR|

直接复用604主符号与倍增边界，611／614静止质量；不重复自由传播实验。382—386、425、522—523按667继承表继续有效，不恢复384已消去的正则性，也不叠加386与425的替代输入。

本轮收到侧聊的历史回顾线索，按当前接口回读了[221](../../../archive_217_222/research_note_221.md)、[222](../../../archive_217_222/research_note_222.md)、[230](../../../archive_223_230/research_note_230.md)、[404§5](../../../archive_370_428/research_note_404.md)、[553](../../../archive_531_553/research_note_553.md)。具体复用：权重相同须继续核实际观测／过程；有限仪器误差标准不能自动覆盖无界场；局部相容而无原Fock正常极限已有概念反例，667是原CAR对象上的具体连接；本轮有限格恒等式不能升级为553所限制的全尺度谱匹配。369模包含仍需同一表示与参考等实际前提，本轮不调用它生成几何。其它线索留在原继承账中，不重开回顾轮次。

## 2. 原质量在同一物理观测中的定义

记614变换后的原配对块Δ为32×32反对称矩阵。在实际外幂载体中Δ=(J⊗I₂)Δ_old(Jᵀ⊗I₂)，重排为spin⊗internal后有Δ=ε⊗M，M=Mᵀ。所有原复Y和右手中微子Majorana均保留。物理观测Y=(w,bar-wᵀ)的新增项定义为

$$
P_\phi=\begin{pmatrix}\Delta_\phi^*&0\\0&-\Delta_\phi\end{pmatrix},
\qquad V_{\lambda,\phi}=\frac\lambda2\sum_xY_x^{\mathsf T}P_\phi Y_x,
\qquad \lambda\in\mathbb R.
\tag{1}
$$

这里使用指数+二次型的原Grassmann约定；配对符号与exp(−质量Hamiltonian)对应的ww／bar-bar排序一致。它是明确的有限Euclidean相互作用选择；静止质量字典不独自证明全时传播、态和记录已与原CAR过程相同。

660—661的变量为ξ=Sχ、η=Tχ、Y=Rη，χ=(b,c,d,e)，η=(b,c,d′,e′)，w=(J_-†v)c。因此

$$
\Phi=RTS^\dagger,\qquad Y=\Phi\xi,\qquad
N_\lambda=N+\lambda\Phi^{\mathsf T}P\Phi.
\tag{2}
$$

在同一index-zero、K_ℓ可逆帧片上，det T=1，S酉，准确得到

$$
T^{-\mathsf T}S^{\mathsf T}N_\lambda ST^{-1}
=N_0+\lambda R^{\mathsf T}PR.
\tag{3}
$$

R只取c、e′，故新项不混入η中的辅助b、d′。不要求J_-†v可逆，也不对可能退化的物理Wick核取逆。定义L=diag(J_-†v,I)，则(c,e′)块为

$$
N_{{\rm W},\lambda}=
\begin{pmatrix}0&-K_\ell^{\mathsf T}\\K_\ell&0\end{pmatrix}
+\lambda L^{\mathsf T}PL.
\tag{4}
$$

原辅助A(E)、bar-B(E)不变，全部原相位和det S保留。固定次序的非零自由分母给

$$
\operatorname{Pf}N_\lambda(E)=
R_{\rm W}(\lambda,\phi;D)\operatorname{Pf}N(E),\qquad
R_{\rm W}=\frac{\operatorname{Pf}N_{{\rm W},\lambda}}
{\operatorname{Pf}N_{{\rm W},0}}.
\tag{5}
$$

这是有限多项式身份，即使A(E)有零模仍成立；只是不能再除以为零的Pf N(E)。R_W与E无关。因此完整S⁹积分、所有E观测及同一物理Weyl观测按该质量因子共同连接，不由少量E采样猜测球面积分。

## 3. 全多项式反射正性

成熟输入采用[Kikukawa–Usui，§IV式(46)—(50)及§V](https://arxiv.org/html/1005.3751v3#S4)：自由overlap的指定Weyl分量代数反射正；相容的逐点反射配对可在该证明结构中加入。原文的具体Yukawa模型不是本轮模型。本轮使用661已映射的原联合泛函Ω₀，并逐项检验式(1)。

当前spin基中J_+†γ₄J_-=I₂，反射交换w与bar-w，反线性且反转乘积次序。常数复M不必取实；式(1)中的共轭及负号恰给

$$
P_\phi+\Theta_{\rm lin}^{\mathsf T}P_\phi^*\Theta_{\rm lin}=0,
\qquad V_{\lambda,\phi}=V_++\Theta V_+.
\tag{6}
$$

这里第二项使用真实空间时间反射，第一式表示其线性字段矩阵；负号来自Grassmann次序反转。偶数时间的link反射没有固定时间片，V_+只含正半区逐点质量。V_+偶，故对所有正半区物理多项式及E观测F，

$$
\Omega_0\!\left(e^{V}\Theta F F\right)
=\Omega_0\!\left(\Theta(e^{V_+}F)(e^{V_+}F)\right)\ge0.
\tag{7}
$$

这是全多项式证明，不是从两点Gram数值推断。有限格Grassmann指数是有限多项式，冻结的质量系数有限。归一还要严格非零，不能只凭式(7)的非负性省掉下一节。

## 4. 任意实强度的严格归一

单位规范链路、常数φ下内部动能为身份。用M的Takagi奇异值μ_a≥0分解物理Weyl块；这是计算质量因子的内部基变换，不声称辅助作用有任意U(16)对称性。给定四维自由核的动量记号为

$$
b(p)=-1+\sum_\nu(1-\cos p_\nu),\quad
s^2(p)=\sum_\nu\sin^2p_\nu,\quad
\omega(p)=\sqrt{b^2+s^2},\quad
r^2(p)=\frac{s^2(p)}{(\omega(p)+b(p))^2}.
\tag{8}
$$

偶数反周期时间使sin p₄≠0，因此分母正，无自由零模。两分量动能逆块S(p)=−iτ_νsin p_ν/(ω+b)，满足SS†=r²I₂和εS(−p)ᵀε†=S(p)†。每个Takagi通道的四分量Nambu行列式比为(1+λ²μ_a²r²)²。p↔−p配对、Pfaffian连续性和λ=0时比值1固定符号，得到

$$
R_{\rm W}(\lambda,\phi;D)
=\prod_{p}\prod_{a=1}^{16}\left(1+\lambda^2\mu_a^2r^2(p)\right)>0
\quad(\lambda\in\mathbb R).
\tag{9}
$$

零质量通道因子为1，不需质量下界。658已证原自由有限盒辅助积分严格正，661自由归一非零；所以在共同归一约定中

$$
\Omega_\lambda(F)=\frac{\Omega_0(e^VF)}{R_{\rm W}},\qquad
\Omega_\lambda(\Theta F F)\ge0,\qquad
\Omega_\lambda(1)=1.
\tag{10}
$$

这是原物理观测代数的同一有限正泛函；不要求660的非局部Φ保候选字段的正时间支撑。661已排除该支撑身份，本轮没有恢复它。一般时变标量、动态标量积分或规范场不能直接使用常数质量Fourier积，另须证明其收敛和归一。

## 5. 原质量、原辅助和共同来源

式(5)给归一Weyl质量部门与原辅助部门的乘积结构，但D是共同输入。改变D时两部门都响应；不能把乘积误读为独立外加的环境。固定D、φ只变λ时有

$$
\partial_\lambda\log Z_\lambda
=\sum_{p,a}\frac{2\lambda\mu_a^2r^2(p)}{1+\lambda^2\mu_a^2r^2(p)}
=\tfrac12\operatorname{tr}(N_\lambda^{-1}\Phi^{\mathsf T}P\Phi).
\tag{11}
$$

右式在非零固定E片上亦成立，因为原辅助权重不依赖λ。一般背景θ同时改变D、φ和实际观测Φ时，完整导数是

$$
\partial_\theta N_\lambda=
N'+\lambda\left(\Phi'^{\mathsf T}P\Phi+
\Phi^{\mathsf T}P'\Phi+\Phi^{\mathsf T}P\Phi'\right),
\tag{12}
$$

这里λ固定；若λ也变化另加λ′ΦᵀPΦ。原帧来源已含于矩阵身份，换帧计算时仍须保det S，不能将其重置。完整积分存在并可求导时，

$$
\partial_\theta\log Z_\lambda
=\partial_\theta\log Z_0+\partial_\theta\log R_{\rm W}.
\tag{13}
$$

这继承661的共同背景原则；新增的是实际非零质量及Φ′补项的定量检验。不是首次发现归一来源或一般链式法则。

## 6. 原矩阵复算及排除的接法

入口已核1×2、2×2盒、每盒两组完整E，λ=.37。准确质量拉回保权重；若将Φ错换为裸局部字段B=diag(J_-†,J_+ᵀ)，权重比变成：

$$
\begin{array}{c|c|cc}
\text{盒}&R_{\rm W}&R_{\rm bare}(E_1)&R_{\rm bare}(E_2)\\
1\times2&1.89769069134&4.48368635366&4.39597683166\\
2\times2&1.96789224651&18.0786164968&5.86207261769
\end{array}
\tag{14}
$$

所以“原质量＋原辅助权重”不能靠往局部镜像字段中直接加同名裸质量完成。省掉Φ改变了模型。这只排除这个具体接法，不排除其它局部实现。

正式第一组另在1×2、2×4、3×4盒保全部16内部通道，分别λ=.37、.37、1.1。Pfaffian比与式(9)相对误差<9.8×10⁻¹⁵；λ来源误差<7.2×10⁻¹⁵。最后一盒比值约4.01526×10¹¹，其虚部绝对残差.00354只占约8.9×10⁻¹⁵，按相对误差验收，不宣称它严格为实的浮点数。反射Gram最小值的−1.3×10⁻¹⁵为舍入诊断，解析式(7)承担正性。

第二组在2×2盒使用原空间Wilson系数κ=1+.2θ及φ=PHI[0]+θ(.1,−.07,.05,.03,−.12)，λ=.37。此θ是声明的共同来源族，不是原全部引力应力。固定代码所列E时，

$$
\partial_\theta\log\operatorname{Pf}N_\lambda=0.680888812660,
\quad \partial_\theta\log\operatorname{Pf}N=0.626961492873,
\quad \partial_\theta\log R_{\rm W}=0.053927319785.
\tag{15}
$$

省掉Φ′会遗漏−0.015011020923−0.000182040206i；原准确全部来源虚部仅约4.3×10⁻¹⁴。矩阵导数残差<2.7×10⁻¹²、Pfaffian差分误差<1.1×10⁻¹⁰、两部门分解误差<2.4×10⁻¹²。非物理虚部来自错误的部分求导，不据此解释为反常。

入口的矩形分块代码错误、首稿与修正结果全部保留；已复算但不重复计入两组新检查。没有重跑旧维数、自由倍增或一般参考表示实验。

## 7. 本轮压缩及下一项

在明确有限范围内，C17原质量现与C01物理观测、C16辅助正泛函、C19共同归一和C22来源一起成立。需要的质量定义受到真实观测映射约束；错误裸质量候选被定量排除。共同对象仍满足

$$
\text{固定背景的共同质量／反射泛函}
\ \not\Rightarrow\
\text{原动态标量、完整Gibbs过程或全尺度量子几何已识别}.
\tag{16}
$$

接[669](../../669/drafts/STATUS.md)：将原h、s的真实热演化及624记录接到本物理质量部门，优先核原标量传递核、非紧目标、反射配对和积分范围；不另拟独立Gaussian标量替代它，也不把Euclidean历史权重直接叫作实际仪器。保留653荷谱、612全时谱、553运行关系及既有空间接口限制。认知设计后置，总目标不改。
