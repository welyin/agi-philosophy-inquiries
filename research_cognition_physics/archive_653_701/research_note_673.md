# 第673轮：共同规范边界、变秩手征字典与有限动态泛函

日期：2026-10-02。接[672](research_note_672.md)、[实际入口](673/drafts/boundary_holonomy_entry.md)。[代码](673/joint_gauss_boundary_functional.py)、[结果](673/joint_gauss_boundary_functional_results.json)、[核验](673/research_round_673_checks.json)、[条件账](673/unified_physics_condition_ledger_673.md)。两组新核验、十八式；主代理审查，无独立代理审查。

## 1. 接续目标与增量

上轮完成672并执行673入口，属于有效进展。本轮先核README、direction、state、最新报告与发布后核验；无运行中的Python研究。回查598、603、605—610、643—644、653、669—672。旧空间接口照667 §1.1复用；384、386／425及523已完成的条件不重新列成缺口。

本轮不继续扫描672的热时。新增两项连接：**把两个规范拼接边界与原标量／辅助／物理来源一起输送；把670的等秩字典扩成矩形字典，定义包含完整边界平均的有限未归一候选。**

[Kikukawa原文§3及结论](https://arxiv.org/html/1710.11618v3)已经讨论全拓扑部门的测度定义，同时区分可容许域保障的光滑局域性与一般背景下尚未证明的非零辅助积分。本轮不把“允许不同拓扑部门”称作新理论；新增的是本项目原非零质量、正则物理观测与完整边界拼接的具体字典。原文不是本候选反射正性或原Hamiltonian身份的证明。

|层次|本轮地位|
|---|---|
|认知动机|表示、边界、物质来源和共同时间应属于同一过程|
|继承输入|原商群、16内部通道、原质量／辅助处方、原H_b及配置热矩|
|声明的候选|有限周期空间盒、两片反周期费米时间、m₀=1；用原H_b标量热核拼接指定overlap物质权重|
|解析结果|任意非零Wilson谱下的矩形来源身份；双边界全群积分、绝对可积及规范不变性|
|数值检验|不等秩代数样本；原全群链路、变化标量及辅助场上的实际边界与来源身份|
|仍缺|双边界积分的正性与严格归一、同一原有序CAR过程、共同连续极限和量子GR|

“不等秩”指某个背景的正负谱子空间维数不同，不表示系统通过连续酉演化随意改变Hilbert空间维数。

## 2. 670字典不必把K固定成方阵

设原单粒子维数n=2r，r为32的倍数。J±是固定γ₅手征基，u、v分别为Wilson H的负、正谱基；H在此点无零模。允许

$$
u:\mathbb C^{r_u}\to\mathbb C^n,\quad
v:\mathbb C^{r_v}\to\mathbb C^n,\quad
r_u+r_v=n,\qquad
K_\ell=J_+^\dagger Dv=J_+^\dagger v\in\mathbb C^{r\times r_v}.
\tag{1}
$$

原D和配对M仍满足670式(2)—(3)。令A=uᵀMu、w=J₋†v、bar-B=J₋†bar-M J₋*；bar-B为原可逆单位配对、Pf=1。670的三角变换按(b,c,d,e)维数(rᵤ,rᵥ,r,r)原样有意义：

$$
d'=d+\bar B^{-1}Kb,\qquad
e'=e-J_+^{\mathsf T}Mu\,b-\tfrac12J_+^{\mathsf T}Mv\,c,\qquad
\det T=1.
\tag{2}
$$

证明只用D v=Q₊v与M Q₊=Q₊ᵀM，并没有用Kₗ为方阵。变换后原二次型分成A、bar-B和

$$
N_{W,0}=
\begin{pmatrix}0_{r_v}&-K_\ell^{\mathsf T}\\K_\ell&0_r\end{pmatrix},
\quad
L=\begin{pmatrix}w&0\\0&I_r\end{pmatrix}
:\mathbb C^{r_v+r}\longrightarrow\mathbb C^n,\quad
N_{W,\lambda}=N_{W,0}+\lambda L^{\mathsf T}P_\phi L.
\tag{3}
$$

r为偶数，交换c、d块的方向符号为(−1)^(rrᵥ)=1。故原S的方向约定和Pf(bar-B)=1都不变。物理观测仍是670的固定n行映射：

$$
\Phi=
\begin{pmatrix}
J_-^\dagger P_v&0\\
-J_+^{\mathsf T}M(P_u+\tfrac12P_v)&J_+^{\mathsf T}
\end{pmatrix},
\qquad
\mathcal N_\lambda=\mathcal N+\lambda\Phi^{\mathsf T}P_\phi\Phi.
\tag{4}
$$

对任意物理Grassmann来源j，

$$
\int d\xi\,e^{\xi^{\mathsf T}\mathcal N_\lambda\xi/2+j^{\mathsf T}\Phi\xi}
=
\frac{\operatorname{Pf}A}{\det S}
\int d\eta\,e^{\eta^{\mathsf T}N_{W,\lambda}\eta/2+j^{\mathsf T}L\eta}.
\tag{5}
$$

这包括矩形L以及不同维数的物理积分变量。若rᵤ为奇数，左侧所有纯物理来源系数均为零：变换后奇数个b只能由二次A配对，物理来源不含b。右侧以奇数维反对称A的Berezin积分为零理解。该奇数代数情况不宣称是原商群实际存在的拓扑部门。

所以变化的是坐标积分的分块尺寸，而非任意改变原可观测对象。归一逆传播子仍仅在相应非零权重、可逆核处使用；本轮只需未归一多项式。

## 3. 两时间片上Wilson零集不必成为硬积分边界

对固定两片空间配置，先把所有时间链路设为单位。反周期两片时间移位T₀满足T₀†=−T₀、T₀†T₀=I，因而S₀=0。任意空间规范链路下

$$
X=W_s+C_s+\gamma_0T_0,\qquad
W_s=\sum_j\left[I-\tfrac12(T_j+T_j^\dagger)\right]\succeq0,\quad
C_s^\dagger=-C_s.
\tag{6}
$$

若Xψ=0，则取实部有〈ψ,Wₛψ〉=0。每项为(1/2)‖(I−Tⱼ)ψ‖²，故所有Tⱼψ=ψ、Cₛψ=0；于是γ₀T₀ψ=0，推出ψ=0。因此**对于每份空间配置，单位时间边界处H一定可逆**。此论证适用任意有限周期空间盒，不要求空间链路静态或交换；仅使用本项两片时间及m₀=1。

H关于两个边界群元a、b的矩阵元实解析，原𝒢=G^V连通，det H在(a,b)=(I,I)不为零。所以

$$
\mu_{\rm Haar}\{(a,b):\det H(a,b)=0\}=0.
\tag{7}
$$

这没有证明所有背景都有统一谱隙。零集上可选sign(0)=0定义有界Borel矩阵；任何不同的零集赋值对下面的Haar积分无影响。H非零时式(5)适用；不同谱秩由同一个固定尺寸左侧表示承接。

继承670的范数估计，cΦ=(1+√5)/2，有

$$
\|D\|\le1,\quad \|\Phi\|\le c_\Phi,\quad
\|\mathcal N_\lambda\|\le 2+|\lambda|c_\Phi^2\|P_\phi\|,
\qquad
|\operatorname{Pf}\mathcal N_\lambda|
\le(2+|\lambda|c_\Phi^2\|P_\phi\|)^n .
\tag{8}
$$

因此有限未归一权重不会因H最小谱值趋零而自动发散。**这里去掉的是有限积分的统一硬谱隙要求，不是光滑性、局域性、背景导数或连续极限所需的谱控制。** 605的硬可容许域／原热态限制继续保留；未把它改名消除。

## 4. 原物质与两个边界必须一起输送

单片配置记x=(q,E)，q=(φ,U)，E为每节点的原S⁹变量。半区代表为xᵢ、xⱼ。a、b∈𝒢分别从j侧向i侧、从i侧向j侧输送；费米反周期符号仍单独保留。记完整带相位候选为Fλ(xᵢ,xⱼ;a,b)。

局部规范变换gᵢ、gⱼ作用为

$$
(x_i,x_j;a,b)\longmapsto
(g_ix_i,g_jx_j;\ g_iag_j^{-1},\ g_jbg_i^{-1}),
\qquad
F_\lambda\ \text{不变}.
\tag{9}
$$

这复用670协变性，现将两处边界都放入同一等式。取gᵢ=I、gⱼ=a，设Ω=ab：

$$
F_\lambda(x_i,x_j;a,b)
=F_\lambda(x_i,a x_j;I,\Omega).
\tag{10}
$$

不能只把a置I，却保持j侧的φ、E、U不变。实际数值中漏掉这些输送，使权重相对误差达到632.5。完整输送残差约3.66×10⁻¹⁴。

若同时有物理来源，令Y′=G_Y(a)Y，则

$$
j'=G_Y(a)^{-\mathsf T}j,\qquad
j'^{\mathsf T}Y'=j^{\mathsf T}Y.
\tag{11}
$$

实际二来源Pfaffian系数的输送误差约3.97×10⁻¹⁴；原质量合同、Φ身份误差约10⁻¹⁵。源坐标的变换不是已经证明所有带电插入是Gauss物理观测；物理来源还需规范不变的收缩和边界输送。

## 5. 包含原H_b的完整双边界候选

用原H_b（保留全部曲目标、电动能、边势与磁势）的正标量热核kτ，固定τ>0。定义

$$
\mathcal K_\lambda(x_i,x_j)=
\int_{\mathcal G}da\int_{\mathcal G}db\
k_\tau(q_i,a q_j)\,
k_\tau(q_j,b q_i)\,
F_\lambda(x_i,x_j;a,b).
\tag{12}
$$

这里不是672的独立Casimir分步近似；kτ指已定义的完整原H_b半群。声明的模型输入是“把这个标量热核与当前overlap物质权重配在一起”，不是声称该配对已等于原H_F。

Haar不变性、原kτ协变性及式(10)给

$$
\mathcal K_\lambda(x_i,x_j)=
\int da\,d\Omega\
k_\tau(q_i,a q_j)\,
k_\tau(a q_j,\Omega q_i)\,
F_\lambda(x_i,a x_j;I,\Omega).
\tag{13}
$$

这正好暴露两个不同角色：a输送中间半区的整个配置，Ω保留最终闭合。643使用一个Gauss闭合元，但还对中间配置积分；其并未许可在固定半区代表时同时删除a。

式(12)在xᵢ、xⱼ各自独立的规范变换下不变，因此下降到规范轨道。反射交换两个半区并使(a,b)变成(a⁻¹,b⁻¹)；若完整来源泛函满足相应共轭身份，结合标量核对称性即可得到𝒦(xⱼ,xᵢ)=𝒦(xᵢ,xⱼ)*。本轮在所测可逆样本核验该反射身份，不把采样扩大成跨全部谱部门的一般证明；Hermitian和正性均须继续验收。

## 6. 全配置与辅助积分确实可积

不是仅凭Haar体积有限宣称完整积分有限。原质量界给常数C及单片w(q)=1+Σ_v F(φ_v)^(-1/2)，使‖Pφ‖≤C[w(qᵢ)+w(qⱼ)]。由式(8)，每个固定物理来源系数被这两个w的有限次多项式控制。

标量正半群的Cauchy–Schwarz和规范协变性给

$$
|k_\tau(q_i,a q_j)k_\tau(q_j,b q_i)|
\le k_\tau(q_i,q_i)k_\tau(q_j,q_j).
\tag{14}
$$

669已经证明原H_b的全部配置热矩，w的每个固定整数次幂受原W的某个有限次幂控制。因此

$$
M_p(\tau):=\int_Q w(q)^p k_\tau(q,q)\,d\mu(q)<\infty
\quad(p<\infty,\ \tau>0).
\tag{15}
$$

结合归一Haar、归一S⁹及(a+b)^n≤2^(n−1)(a^n+b^n)，得到某个有限Cλ：

$$
\int d\mu(q_i)d\mu(q_j)d\nu(E_i)d\nu(E_j)da\,db\
|k_\tau k_\tau F_\lambda|
\le C_\lambda[M_n(\tau)M_0(\tau)+M_0(\tau)M_n(\tau)]<\infty.
\tag{16}
$$

同理全部固定Grassmann来源系数可积。Fubini换序和式(13)因而有支配依据；不是把振荡路径权当作概率。积分后的未归一总权重是λ的有限多项式，但**未证明它在所有参数下为实、非负、非零或可以归一成物理态**。对λ或不影响Wilson核的有限质量参数，可逐系数求导；不能据此对穿越Wilson零集的任意规范／几何来源交换微分和积分。

## 7. 可复算证据及当前验收边界

第一组用原64维单粒子配对和质量，分别取(rᵤ,rᵥ)=(30,34)、(31,33)、(34,30)。这些是检查矩形恒等式的代数样本，**不是人为指定拓扑荷的Wilson物理配置**。实际积分维数为66、65、62，原物理观测仍为64行。三角及Φ误差小于4.5×10⁻¹⁶；偶数配对的权重／二来源／四来源相对误差小于1.4×10⁻¹⁴；奇数辅助块权重约7×10⁻²⁹，解析值为零。

第二组使用原2空间×2时间、全部16通道、四份不同原φ及不同E，真实颜色／弱／超荷边界a、b；同时核式(10)、式(11)、质量合同、反射和原H⁵距离。另给三份较强规范场，最小Wilson谱隙为

$$
0.0279867,\quad 0.0178884,\quad 0.0297422,
\tag{17}
$$

权重与来源字典正常定义。旧探针的“.5谱隙门槛”只是其样本检查，现不将其用作物理积分域。其中两份样本的完整费米矩阵各有四个约10⁻¹⁵的数值零模，直接Pfaffian约10⁻⁵⁰及10⁻⁵³；它们的浮点相位与相对比值不可靠，不能作为物理复权重或反射失败的证据。另一份可逆样本反射相对误差约2.9×10⁻¹⁴。

[首版代码](673/drafts/first_boundary_attempt.py)和[首版结果](673/drafts/first_boundary_results.json)保留；后补核共同物理来源，并保存[来源核验版](673/drafts/source_boundary_attempt.py)、[结果](673/drafts/source_boundary_results.json)、[相对反射检查失败诊断](673/drafts/complex_reflection_diagnostic.json)及[审查前稿](673/drafts/research_note_673_before_review.md)。最终版明确屏蔽近零权重的相对数值判据，没有删除这些样本；初稿未经充分证明的一般Hermitian／实性表述已收回。没有执行完整Haar积分或声称计算了其配分函数。解析部分给的是有限候选的存在性与绝对可积：

$$
\text{全边界、全配置、原质量及来源有共同有限泛函}
\ \not\Longrightarrow\
\text{该泛函正且等于原Gauss量子过程}.
\tag{18}
$$

## 8. 下一项

接[674入口](674/drafts/STATUS.md)：直接研究式(12)—(13)在物理半区观测上的正性及其与643原有序影响的身份。数值如使用分步热核，须明确与完整kτ的区别及收敛条件；不能把有限抽样均值当作完整投影。优先核同一过程的缺口，不继续调整672热时或设计认知组织。原空间条件及已消去假设保留，统一目标不变。
