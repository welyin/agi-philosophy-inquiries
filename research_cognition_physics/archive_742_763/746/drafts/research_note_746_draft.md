# 第746轮：原连通图中的占据读出与共同资源

日期：2026-10-04。接[745](../../research_note_745.md)。[代码](../joint_connected_population_readout.py)、[结果](../joint_connected_population_readout_results.json)、[核验](../research_round_746_checks.json)、[条件账](../unified_physics_condition_ledger_746.md)。两组复算、十八式；主代理审查，无新增独立代理审查。

## 1. 增量与研究方向调整

745认证的是一节点sin s联合读取。继续推向连通图时，不应把所有边单独调弱，却保留与这些边不相容的节点几何。本轮回查后使用**项目已经存在的Higgs幅度T=|X|²/2读口**，发现其占据响应在四阶就出现，且该阶精确不含其他节点的标量势与原跳跃。因此可以直接处理原固定边强度。

**结论：在下述原固定有限图、固定正几何、原保持物种模块的规范协变跳跃模型中，存在正规、严格Gauss、有限能源准备，使本地原T二元读口在充分小正时间区分未知偶编码的空态／配对态。所有原相互作用都在演化中保留。**

输入与读口位于同一个节点。本轮没有证明信息已传到其他节点，没有证明无扰占据测量、自治准备或自治终端。744的单次sin s盲性、745的sin s联合信号均保留；本轮换用旧菜单中的T读口，不声称已将原sin s协议推广到全部连通图。

|层次|本轮地位|
|---|---|
|认知动机|局部可读能力应与加入更大共同系统相容，资源与后态同时交代|
|继承模型|574／589固定正几何、598原完整CAR与Gauss、623共同域|
|旧操作输入|624／723已有T字段平方根读口；没有新增指针、耦合或物种|
|准备输入|同一未知偶编码；明确的反射对称、正常有限能源玻色准备|
|解析增量|四阶占据响应的固定图局部性、严格符号、完整未知态仪器及同态能源|
|计算|完整96模原质量／非平凡链路的算符词校准；严格有理数差密度与区间积分|
|未消除|原图、几何、物种与耦合的输入地位；终端装置、连续、动态几何及统一目标|

## 2. 为什么不直接弱化全部边

在原共形几何，节点体积和测地边系数共用ε与ψ，直接得到

$$
w_v=\epsilon^3\psi_v^6,\qquad
k_{vw}=\epsilon\left(\frac{\psi_v+\psi_w}{2}\right)^2
=\frac14\left(w_v^{1/6}+w_w^{1/6}\right)^2 .
\tag{1}
$$

所以均匀w=1时k=1。固定w同时令所有k趋零不属于同一共形几何路径。598的一般有限图系数类可以另行声明这种比较，但它不自动满足当前共同几何合同。本轮不采用该捷径。589的剪切也必须同步改变其他几何系数，不能把三个方向全都独立缩小。

回查723的四参考顺序反例、717的质量力与718的完整历史，均不等于本轮固定图的未知占据读出；其通用域、CP与资源工具直接复用。382—386、425、522—523及604、649／699不重证。

## 3. 原完整图、编码与读口

固定任意有限连通图Γ，固定正外部几何系数；不含动态量子几何。按598记

$$
H=\mathcal T+V_{\rm sc}+B,\qquad
\mathcal T=\sum_v-\frac{\hbar^2}{2w_v}\Delta_{\mathcal K,v}+\mathcal T_{\rm el},\quad
B=\sum_v(D_v+M_v)+C_{\rm hop}.
\tag{2}
$$

V_sc包括原位势、原全部测地／Log边势和磁势，是Fock恒等上的标量乘法。D_v是原Dirac质量，M_v是原sterile Majorana。C_hop保各原规范物种模块，可含原spin矩阵、复系数和非平凡规范链路；不依赖节点标量。不同物种之间的新跨节点配对或跳跃不在本轮声明范围。

在选定v编码空态与sterile双占据，其余CAR初始为空。玻色准备在每个节点规范不变，链路取归一常数Haar波函数：

$$
|0_L\rangle=\Omega_\Gamma,\qquad
|1_L\rangle=e^{i\arg Y_s}a_v^\dagger b_v^\dagger\Omega_\Gamma,\qquad
V_\Gamma|j\rangle=\psi_\Gamma\otimes|j_L\rangle .
\tag{3}
$$

以下先用局部Gaussian计算四阶系数，再在第7节给全图光滑紧支撑准备的严格实现，不将一节点的Gaussian高阶域未经检验直接搬到图上。

沿原T菜单，设

$$
T_v=\frac{|X_v|^2}{2},\quad
E_{T,\pm}=\frac12\pm\frac14\sin T_v,\quad
L_{T,\pm}=\sqrt{E_{T,\pm}},\qquad
\mathcal I_{\pm,t}(\rho)=L_{T,\pm}U(t)V_\Gamma\rho
V_\Gamma^\dagger U(t)^\dagger L_{T,\pm}.
\tag{4}
$$

ρ为未知编码态，可带被动参考。完整物理输出保留，求和完全正且保迹；未压回编码子空间。准备与终端仪器仍是输入。

## 4. 旧对称性对新读口的准确作用

744的全局singlet反射与CAR四分之一相位联合对称Θ保原H。反射对称准备给ΘVΓ=VΓZ；T为反射偶。因此

$$
[\Theta,H]=0,\qquad \Theta E_{T,+}\Theta^\dagger=E_{T,+},
\qquad
ZM_{T,+}(t)Z=M_{T,+}(t),\quad
M_{T,+}(t)=V_\Gamma^\dagger U^\dagger E_{T,+}UV_\Gamma .
\tag{5}
$$

与sin s效果的反射奇性不同，T效果在编码中严格为

$$
M_{T,+}(t)=\alpha(t)I+\beta(t)Z,\qquad
\frac14I\le M_{T,+}(t)\le\frac34I.
\tag{6}
$$

这里只证明效果对占据取向对角，并不证明实际后态保留占据，也不等于锐投影仪器。

## 5. 四阶算符词为什么不含外部边

设A=sin(T_v)/4；常数1/2在差中消去。对任意含微分的CAR算符词，令ℓ取两个固定Fock输入的对角期望之差。它保留一个作用在玻色波函数上的算符；不能把ℓ误作全系统迹。

用x=φ/√F、a=|x_H|²、b=x₅²、u=1+(a+b)/6，则T=a/u及原逆度量给

$$
K_x^{-1}=I+\frac{xx^T}{6},\qquad
\operatorname{grad}_{K}T=(2x_H/u,0),\qquad
[B,[\mathcal T,A]]
=\frac{\hbar^2}{w_v}\frac{2f'(T)}u D_v,\quad f(T)=\sin T/4 .
\tag{7}
$$

这直接应用717的完整质量力，不丢角微分。特别是首次质量插入只给局部Dirac块；跳跃不依赖x_v，所以此处不出现。

在ad_H⁴ A中按B插入次数分类。纯玻色项的ℓ为零；一次B及其任意玻色导数仍是Dirac异物种、Majorana变粒子数或跨节点双线性，均无所需对角差。两次B时，唯一可能非零的词为

$$
\mathcal W_2(A)=
\operatorname{ad}_B^2\operatorname{ad}_{\mathcal T}^2A+
\operatorname{ad}_B\operatorname{ad}_{\mathcal T}
\operatorname{ad}_B\operatorname{ad}_{\mathcal T}A+
\operatorname{ad}_{\mathcal T}\operatorname{ad}_B^2
\operatorname{ad}_{\mathcal T}A .
\tag{8}
$$

若将式(8)中额外的一次𝒯换成V_sc，则词逐一为零：V与矩阵乘法对易，或先产生一个标量乘法再与B对易。一次B词即使含V的微分，其CAR对角差也仍为零。这证明不是“势小所以忽略”。

式(8)中最先遇到的B必须是v节点质量；其他节点的偶质量与它及v微分对易。若另一个B来自跨节点跳跃，词在远端留下奇CAR次数，远端真空对角期望为零。链路／其他节点动能即使微分这些系数也不改变该奇偶判断。因而

$$
\ell(\mathcal W_2(A))=
\ell\!\left(
\operatorname{ad}_{B_v}^2\operatorname{ad}_{\mathcal T_v}^2A+
\operatorname{ad}_{B_v}\operatorname{ad}_{\mathcal T_v}
\operatorname{ad}_{B_v}\operatorname{ad}_{\mathcal T_v}A+
\operatorname{ad}_{\mathcal T_v}\operatorname{ad}_{B_v}^2
\operatorname{ad}_{\mathcal T_v}A\right).
\tag{9}
$$

三次B只剩ad_B²作用于式(7)的局部D_v。三个局部质量因子中，奇数Majorana不能给定粒子数对角；偶数Majorana时有奇数Dirac，不能恢复原左右物种占据。一次跳跃在远端留奇数，亦为零。两次跳跃可回到v，但保物种，剩下的一次Dirac仍改变物种，不能给对角。其他节点质量同样服从此计数。因此

$$
\ell(\operatorname{ad}_B^3\operatorname{ad}_{\mathcal T}A)=0,\qquad
\ell(\operatorname{ad}_H^j A)=0\ (j=0,1,2,3),\qquad
\ell(\operatorname{ad}_H^4A)=\ell(\mathcal W_2(A)).
\tag{10}
$$

第三阶中的局部两质量词使用式(7)：[D_v+M_v,D_v]=[M_v,D_v]是异常配对块，对角差为零。由此明确关闭所有较低阶，而不只依赖数值抵消。

式(9)(10)是固定图局部系数定理，允许任意已声明的固定边强度；没有说后续高阶或有限时间演化不依赖其他节点。

## 6. 同一原波包的四阶严格符号

取v上的比较Gaussian exp(−|x|²/2)，在真实H⁵测度下归一化。745已经保存其精确四阶差密度，本轮利用式(9)说明该局部密度在全图中的用途。令g=|Yν|²=109/625：

$$
P_4(a,b)=g\left[-\frac{76}3-8b-\frac{17}3a+
\frac{31}9a(a+b)+\frac23a(a+b)^2\right].
\tag{11}
$$

它来自完整五方向微分与原32CAR，不是先删去其他物种后的拟合。两次𝒯_v使一般节点体积只贡献w_v⁻²。原位势、其他节点及跳跃对这一阶无贡献已由算符词证明，不能把它误读成这些作用整体可删除。

用745的真实归一化Z，定义T读口的未归一化系数

$$
J_T=\int_0^\infty\frac{R^4e^{-R^2}}{\sqrt{1+R^2/6}}
\int_{-1}^1(1-z^2)P_4(R^2(1-z^2),R^2z^2)
\frac{\sin\!\left(\frac{6R^2}{6+R^2}(1-z^2)\right)}4\,dz\,dR .
\tag{12}
$$

本轮只复用745的区间算法和Cauchy误差方法，重新计算当前角矩、sin余项和界。取sin的21个非零项至41次、500段20次精确Newton—Cotes。复半径1/2中|6ζ²/(6+ζ²)|<13；P₄的径向次数至多6。解析求积、级数与真实尾部误差分别小于8.171×10⁻¹⁵、1.215×10⁻¹⁷、1.757×10⁻²¹。得到[严格证书](higgs_readout_certificate_results.json)：

$$
-0.151458842046344<J_T<-0.151458842046326<0,\qquad
c_T=J_T/Z<0,\qquad
c_T\simeq-0.20107945876584 .
\tag{13}
$$

只有负号依赖严格证书；归一化小数是辅助校准。此次四阶系数只含原诊断Yν；不把它说成已从认知原则或RG精确算出的物理耦合。

## 7. 全图有限能源准备及真实正时间

为严格落实原全图高阶域，在v取反射偶、径向光滑截止χ_R，|x|≤R时为1、|x|≥2R时为0。其他节点任取正常反射偶径向紧支撑波函数；群链路常数。v的波函数为

$$
\psi_{v,R}=\frac{\chi_R(x)e^{-|x|^2/2}}
{\|\chi_R e^{-|x|^2/2}\|},\qquad
\psi_{\Gamma,R}\otimes|j_L\rangle\in
C_c^\infty((H^5)^V\times G^E)\otimes\mathcal F
\ \cap\ \mathcal H_{\rm phys}.
\tag{14}
$$

原光滑全H把该核映入自身，因此全部有限H矩存在。它不是隐含的低能准备：截止和准备本身是输入，能源须计入。没有把Gauss物理空间作一般区域张量分解；这些是在运动学空间明确构造后已经不变的特殊态。

式(9)中的局部有限阶微分算符系数至多多项式增长；原f(T)及其微分为有界三角函数乘有理系数，分母1+|x|²/6不为零。截止导数支撑在Gaussian尾部，故通过带多项式权的L²支配收敛得

$$
c_{T,R}\longrightarrow c_T<0,\qquad
\exists R_0:\ R\ge R_0\Rightarrow c_{T,R}<0,\qquad
p_{1,+}(t)-p_{0,+}(t)
=\frac{c_{T,R}}{24w_v^2}t^4+o_{\Gamma,R}(t^4).
\tag{15}
$$

因此固定Γ、R≥R₀后，所有充分小正t都有严格负差。不必调弱原边。此处没有计算R₀、可用时间窗或最小可观测概率差；余项依赖图及准备，不能因局部首项不依赖图就宣布跨图统一速度或资源界。

## 8. 未知输入、后态与代价属于同一过程

选择的玻色态在singlet反射下为偶，原Majorana期望的x₅均值为零；Dirac与跳跃在编码中无对角／编码间矩阵元。因此

$$
V_\Gamma^\dagger H V_\Gamma=E_{\Gamma,R}I,\qquad
E_{\Gamma,R}<\infty .
\tag{16}
$$

所有未知编码态有同一初始总平均能源，但不声称同一H分布或同一能源噪声。资源不是通过给某个输入额外初始平均能量藏进去。

读取后所有分支仍使用式(4)。598的共同注能身份对既有T菜单直接适用。M=2时

$$
K^{-1}(dT,dT)=Fh^2(1-h^2/12)
\le24q(1-q)^2\le32/9,\quad q=h^2/12,
\qquad
0\le D_T:=\mathcal I_T^*(H)-H
\le\frac{\hbar^2}{9w_v}.
\tag{17}
$$

这里用Σ|L′_r(T)|²≤1/16。演化保原H，非选择读取的总注能有界，实际条件分支也有有限移位能源。有限次数重复读取按原过程逐次计费，不免费重置为原准备。

任何已有原来源Jγ及实际后态仍按同一物理输出计算：

$$
\langle J_\gamma\rangle_{\rm out}
=\sum_r\operatorname{Tr}
(J_\gamma L_{T,r}U_\gamma V_\Gamma\rho V_\Gamma^\dagger
U_\gamma^\dagger L_{T,r}),\qquad
\mathcal I_T^*(\partial_\gamma H)-\partial_\gamma H
=\partial_\gamma D_T .
\tag{18}
$$

最后的算符／形式身份在L与所取外部几何参数无关时成立；变化准备或动力的来源导数仍须包括它们自身响应，复用623—625及704，不能只对D_T求导替代全部后态导数。这里没有构造新的Einstein解。

## 9. 核验、适用边界及下一步

第一组保全部96模、三节点三条非平凡规范／spin边，逐项校准一次质量、边—质量及三质量／跳跃词的对角抵消，并核原梯度、体积—边系数及注能界。第二组复算式(11)的精确系数和式(13)全部误差项。数值没有模拟原全图的玻色时间演化；图结论来自第5、7节的解析论证。

保存的[各向异性Gaussian诊断](local_population_jet_probe_results.json)显示，仅改变singlet宽度并未产生单节点sin s偶报告的四阶差；它不是全准备的不可能定理，不重新计轮次。它促使本轮回到已经存在的完整读口菜单。

**本轮减少的独立条件：** 在原固定图联合过程内，局部占据信息可读取不再需要独立新增探针或“边足够弱”的假设。其代价与原几何、原能源及同一来源共同核算。

下一项[747](../../747/drafts/STATUS.md)：把本地已可读的信息连接到其他节点的实际记录，回查577的原传播通道并保未知态与真实后态；输入、控制、消息和来源必须属于同一过程。若只得到本地读取或既有相位通信的复述，不计新轮次。不重新优化本轮区间积分或添加无成本经典通信。
