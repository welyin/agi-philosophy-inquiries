# 1123证明：固定局部CAR通道的完整四动量变化误差

2026-10-11。接[1121](../1121/proof.md)。本轮给一项有限强度的条件证书：若近似操作的全部Kraus确实支撑在原两味模式内，逐Kraus保该模式的CAR占据数，则普通通道误差控制**真实完整场四动量的变化之差**。额外预算是实际生成元作用于原模式的漏出范数。无需另给全场统一二阶矩，但模式支撑和输入／输出合同不能省略。

新增认知公理0、科学校准0。这是材料操作的必要来源账与充分误差界，未构造实际控制器，未把通道误差提升为完整局部应力、来源Hessian或全部能源统计误差。

## 1. 旧材料和净增量

[592§6—7](../../archive_585_628/research_note_592.md)、[593§5](../../archive_585_628/research_note_593.md)已给同一物理能源二阶矩、迹距离与来源均值的尾界，也已有稀薄高能尾反例。[623—625](../../archive_585_628/research_note_625.md)给同一源族的一二阶图范数逼近及来源插入后的尾控制。[879](../../archive_873_893/research_note_879.md)和[881§4—7](../../archive_873_893/research_note_881.md)给有限区域CAR和共同系数逼近；881的辅助数权重不能代替真实物理能源。[1110§6](../../archive_1106_1119/1110/proof.md)明确指出完整输出的diamond误差不自动控制无界四动量。

本轮不重复这些一般尾定理。新连接是：在1121**同一**平直双味Dirac分支，物理生成元相对于固定模式代数有一个精确有限的交换子代表；整个无限余部严格抵消，因此得到线性且可饱和的误差常数。有限码若只检查原两味末端，仍不足以套用它。

## 2. 同一真实正能表示与操作合同

完全沿用[1121§1—3](../1121/proof.md)的正质量双味Dirac真空表示，记

$$
G_0=h_D=-i\alpha^j\partial_j+\beta m,\qquad
G_j=p_j=-i\partial_j,\qquad
\mathbb P_0=\mathbb H\ge0,\qquad \mathbb P_j=\text{物理空间动量}.
\tag{1}
$$

这里$G_\nu$作用于单味Cauchy空间，物理表示用$S=1_{h_D<0}$填充负支。$G_0$不是正算子；$\mathbb H$是粒子和空穴的正能源。四个真实平移分别满足

$$
e^{ia\mathbb P_\nu}c(v)e^{-ia\mathbb P_\nu}
=c(e^{iaG_\nu}v).
\tag{2}
$$

$c(v)$按$v$反线性。紧支Cauchy模式不是纯正能粒子态。本轮没有改变Dirac海、重选能源零点或用有限矩阵替换$\mathbb P_\nu$。

固定归一$f\in C_c^\infty(B_R,\mathbb C^4)$，令

$$
F=f\otimes\mathbb C^2,\qquad
\mathcal A_F=\operatorname{CAR}(F)\simeq M_4,\qquad
N_F=\sum_{a=1}^2c^\dagger(f\otimes e_a)c(f\otimes e_a).
\tag{3}
$$

$N_F$是这两份Cauchy CAR模式的占据数，不是物理粒子／反粒子总数。逐Kraus保$N_F$保证这里的局部电荷旋转不变性，也保证每个Kraus为偶元。

考察有限Kraus通道

$$
\Phi^*(A)=\sum_{\ell=1}^L K_\ell^\dagger A K_\ell,\qquad
K_\ell\in\mathcal A_F,\quad
\sum_\ell K_\ell^\dagger K_\ell=I,\quad [K_\ell,N_F]=0.
\tag{4}
$$

这是明确的结构合同。仅“通道整体偶”或“平均数守恒”不能替代逐Kraus条件；不在此合同的过程仍可用其它预算研究。$M_4$上总可用有限Kraus表达，但(4)要求存在一份**全部**落在所声明代数且逐个保$N_F$的表达。

目标为同一个有限角度混味门

$$
U_\theta=\exp[-i\theta(c_1^\dagger c_2+c_2^\dagger c_1)],
\qquad \Psi(\rho)=U_\theta\rho U_\theta^\dagger.
\tag{5}
$$

这些是局部代数操作合同；双点模式核不等于已建成的点状Hamiltonian或任意短时实施权限。候选$\Phi$的实际控制、时长、反冲材料仍须独立认证。

## 3. 无界外侧如何精确抵消

逐$\nu=0,1,2,3$定义真实模式数据

$$
\mu_\nu=\langle f,G_\nu f\rangle\in\mathbb R,\qquad
r_\nu=(G_\nu-\mu_\nu)f,\qquad
\sigma_\nu=\|r_\nu\|.
\tag{6}
$$

令有限秩自伴算子和相应CAR元为

$$
C_\nu=(|r_\nu\rangle\langle f|+|f\rangle\langle r_\nu|)\otimes I_2,
\qquad Q_\nu=q(C_\nu)=\sum_{ab}(C_\nu)_{ab}c_a^\dagger c_b.
\tag{7}
$$

若$\sigma_\nu>0$，$C_\nu$只在
$E_\nu=F\oplus(r_\nu/\sigma_\nu)\otimes\mathbb C^2$上非零，矩阵为

$$
C_\nu=\sigma_\nu
\begin{pmatrix}0&I_2\\I_2&0\end{pmatrix},\qquad
\operatorname{spec}C_\nu=(\sigma_\nu,\sigma_\nu,-\sigma_\nu,-\sigma_\nu),
\qquad \boxed{\|Q_\nu\|=2\sigma_\nu}.
\tag{8}
$$

最后一个范数是完整四模式CAR的范数，来自占据本征值求和；不是一粒子范数。$\sigma_\nu=0$时$C_\nu=Q_\nu=0$。有限CAR在真实Dirac表示中的含幺$*$表示忠实，故(8)保持。

**定义域与强对易。** 对$v\in F$，$v\in D(G_\nu)$。由(2)，$a\mapsto\alpha_a(c(v))$范数可微，所以$c(v)$及其伴随保持$D(\mathbb P_\nu)$，并在此域上满足

$$
[\mathbb P_\nu,c(v)]=-c(G_\nu v)
=[\mu_\nu N_F+Q_\nu,c(v)].
\tag{9}
$$

有限多项式、从而$\mathcal A_F$所有元素同样保持该域。定义

$$
R_\nu=\mathbb P_\nu-\mu_\nu N_F-Q_\nu,
\qquad D(R_\nu)=D(\mathbb P_\nu).
\tag{10}
$$

有界自伴扰动定理保证$R_\nu$自伴。(9)及其伴随式给
$R_\nu A\psi=A R_\nu\psi$，对所有$A\in\mathcal A_F$、$\psi\in D(R_\nu)$成立。由此$A$与$(R_\nu-z)^{-1}$对易，再与$R_\nu$全部谱投影对易。这里的“强对易”指$R_\nu$的谱投影与该有限代数对易；无需假设四个$R_\nu$彼此强对易。

因此整个无限余部$R_\nu$在每个Kraus下严格不变。对$\psi\in D(\mathbb P_\nu)$，

$$
\begin{aligned}
\sum_\ell K_\ell^\dagger\mathbb P_\nu K_\ell-\mathbb P_\nu
&=\sum_\ell K_\ell^\dagger Q_\nu K_\ell-Q_\nu\\
&=:\mathcal D_\nu(\Phi).
\end{aligned}
\tag{11}
$$

第一步同时用了$R_\nu$强对易、$\sum K_\ell^\dagger K_\ell=I$以及逐Kraus保$N_F$。右端是有界算符，故给左端唯一有界延拓。未假设$G_\nu$保持$E_\nu$；$G_\nu r_\nu$及更远模式归入$R_\nu$，不截掉真实外侧传播。

在具有相应有限矩且定义域合法的状态上，(11)就是完整物理四动量均值的变化。即使某个输入没有可相减的绝对均值，右端仍定义有界变化证书；不能把两个发散均值形式相减。所有有限个$D(\mathbb P_\nu)$中的向量输入均可直接使用，推广到其它形式域须另核该输入合同。

若采用正常序，$Q_\nu=:Q_\nu:_S+\operatorname{Tr}(SC_\nu)I$。在$\Phi^*(Q_\nu)-Q_\nu$中同一标量因保幺性抵消；这不允许从1121的最终$q(B_\nu)$中擅删其正常序补偿。二者给的是同一真实平移导数。

## 4. 从有限通道误差到完整四动量变化

令$\Phi_F,\Psi_F$为完整$M_4$输入、完整$M_4$输出的通道，取半diamond距离

$$
\delta=\tfrac12\|\Phi_F-\Psi_F\|_\diamond.
\tag{12}
$$

它包括任意被动参考和全部模式占据部门，不只包括$N_F=1$的逻辑码。把$F$排在CAR基的最前面，则补空间奇元带$F$的宇称串；(4)的偶Kraus在分解中为$K_\ell\otimes I$，所以诱导到$E_\nu$及任意重数空间的过程是普通张量恒等扩展。没有把原外侧态重置成独立准备。

由diamond范数的对偶完全有界界与(8)、(11)，

$$
\boxed{\ 
\|\mathcal D_\nu(\Phi)-\mathcal D_\nu(\Psi)\|
\le 2\delta\,\|Q_\nu\|
=4\sigma_\nu\delta,\quad \nu=0,1,2,3.\ }
\tag{13}
$$

证明只需
$\mathcal D_\nu(\Phi)-\mathcal D_\nu(\Psi)
=(\Phi^*-\Psi^*)(Q_\nu)$，再在四模式矩阵因子上应用通道对偶界。任意被动参考、原场外侧关联和物理表示重数均不改变常数。对任意参考有界效果$B$，相应
$(\mathcal D_\nu(\Phi)-\mathcal D_\nu(\Psi))\otimes B$的范数至多$4\sigma_\nu\delta\|B\|$；因此原输入与参考的联合矩也保同一证书。

四个$\sigma_\nu$来自同一$f$和同一真实$G_\nu$。它们可由原模式的空间梯度／Dirac作用矩计算或独立估计，是本合同下足够的额外物理预算；不宣称它们是在所有可能操作类别中必要且最小的预算。尤其$\sigma_0$是Cauchy $h_D$的模式漏出量，不是将$f$当纯正能粒子后的能源方差。

目标门的变化由[1121§3—4](../1121/proof.md)的同一计算给

$$
\|\mathcal D_\nu(\Psi)\|=4\sigma_\nu|\sin(\theta/2)|,
\qquad
\|\mathcal D_\nu(\Phi)\|\ge
4\sigma_\nu\bigl(|\sin(\theta/2)|-\delta\bigr)_+.
\tag{14}
$$

时间方向的相同步骤用$G_0=h_D$及真实正$\mathbb H$的(2)，而非负能Fock Hamiltonian。

## 5. 常数可饱和，且原一粒子码误差可以完全看不见

固定一个非平凡有限角度，例如$\theta=\pi/3$。令局部宇称
$Z_F=(-1)^{N_F}$，并取$0\le\eta\le1$，

$$
\Phi_\eta=(1-\eta)\operatorname{Ad}_{U_\theta}
+\eta\operatorname{Ad}_{Z_FU_\theta}.
\tag{15}
$$

两份Kraus$\sqrt{1-\eta}U_\theta$、$\sqrt\eta Z_FU_\theta$满足(4)。因为$Z_FQ_\nu Z_F=-Q_\nu$，

$$
\mathcal D_\nu(\Phi_\eta)-\mathcal D_\nu(\Psi)
=-2\eta U_\theta^\dagger Q_\nu U_\theta,\qquad
\|\mathcal D_\nu(\Phi_\eta)-\mathcal D_\nu(\Psi)\|=4\sigma_\nu\eta.
\tag{16}
$$

$Z_F$同时有$+1,-1$本征值，所以$\tfrac12\|\operatorname{Ad}_{Z_F}-\mathrm{id}\|_\diamond=1$。上界来自两个通道的距离至多1；下界由两个宇称部门的等幅相干输入直接取得。因此$\delta=\eta$，(13)在有限角度例中饱和。若要求总宇称为偶，可把相干输入与一个参考费米模式配对为
$(|0_F0_R\rangle+|1_F1_R\rangle)/\sqrt2$，仍取得该距离。

但在原$N_F=1$码$\operatorname{span}\{c_1^\dagger|0_F\rangle,c_2^\dagger|0_F\rangle\}$上，$Z_F=-I$，故$\Phi_\eta$和$\Psi$对全部未知码输入及任意参考的通道**完全一致**。该码通道误差为0，而(16)非零。这只说明码误差不能领取完整场变化范数；若码与余部采用固定数部门的乘积准备，(16)的某些压缩也可能为0，不能反过来声称每个码输入都有非零平均误差。

这里的逻辑码是有限CAR内部的两维码，没有把局部CAR真空记号$|0_F\rangle$等同于全场Dirac真空。

## 6. 同一场中的隐藏模式反例

固定同一$B_R$和$f$。由于$f$紧支于球内部，可在$B_R$内、与$\operatorname{supp}f$不交的一个小球选$x_0$。取归一实标量光滑紧支函数$h$乘固定旋量，满足$\sigma_j(h)>0$，并设

$$
g_\lambda(x)=\lambda^{3/2}h(\lambda(x-x_0)),\qquad
\lambda\ge\lambda_0.
\tag{17}
$$

所有$g_\lambda$仍是原双味Dirac场的紧支Cauchy模式，支撑留在该小球，并与$f$及其各阶导数正交。它们不是新增物种或任意离散动量模型。空间动量数据满足
$\mu_j(g_\lambda)=0$、$\sigma_j(g_\lambda)=\lambda\sigma_j(h)$。

令$V_\lambda$在$g_\lambda$的两味上作同一个有限角度$\theta=\pi/3$的混味。比较完整场过程$\operatorname{Ad}_{U_\theta}$与$\operatorname{Ad}_{V_\lambda U_\theta}$。由于$V_\lambda$是外侧偶元，两者在原$F$末端及任意被动参考的约化通道精确相同；所有只访问原$F$的码实验误差均为0。另一方面，支撑不交使两份有限变化算符也相互对易，所以

$$
\left\|\mathcal D_j(\operatorname{Ad}_{V_\lambda U_\theta})
-\mathcal D_j(\operatorname{Ad}_{U_\theta})\right\|
=4\lambda\sigma_j(h)|\sin(\theta/2)|
\longrightarrow\infty.
\tag{18}
$$

此处$\mathcal D$直接按完整场门的真实平移导数定义；$V_\lambda U_\theta$不满足固定$F$支撑合同(4)，所以不违反(13)。每个有限$\lambda$的模式均光滑、门有有界交换子并保持物理域。范数端点可在真实有限CAR谱投影内取得；该投影保持光滑物理核，选非零投影向量可得有限物理矩的见证。没有宣称所有$\lambda$共用一个有界能源准备预算或同一输入向量。

该反例说明必须核“全Kraus支撑不泄漏”，不能只核原端口通道。它尚未构造这些局部门的实际有限时控制；反例和目标沿用同一种局部操作合同。

## 7. 可以怎样帮助1121的有限码任务

若已认证(4)，并有实际固定有限输入等距$J$及原准备，未必需要验证全$M_4$的(12)。可把$E_\nu$的**完整四模式输出**保留，检查同一输入码（含任意参考）诱导的两份输出通道$\Phi_{E_\nu,J},\Psi_{E_\nu,J}$，记半diamond误差$\delta_{\nu,J}$。由于无限外侧已在(11)精确抵消，

$$
\|J^\dagger[\mathcal D_\nu(\Phi)-\mathcal D_\nu(\Psi)]J\|
\le4\sigma_\nu\delta_{\nu,J}.
\tag{19}
$$

这里的输出限制是对真实$\operatorname{CAR}(E_\nu)\simeq M_{16}$因子取态，不重置原关联。对自伴压缩，范数由所有码输入期望的上确界得到；每个期望差是同一四模式输出中有界$Q_\nu$的差，故(19)成立。可分别检查四个$E_\nu$，也可保它们有限并集。若$J$是更一般固定输入CP准备，可将左侧改为其保幺对偶作用，不必免费假设纯准备。

这给1121有限码一个明确可交付的合同：真实$f,G_\nu$预算、逐Kraus支撑和占据数、实际输入$J$、扩展模式输出误差。仅原接收者末端或原两维码的误差不足；也没有证明旧准备／读出权限已允许对所有扩展模式作实际过程层析。

完整控制器若与材料共同满足总四动量守恒，其在同一输入合同上的平均变化须补偿(11)，由(13)或(19)获得同一容差。被动参考稳定不等于已给出承担反冲的物理参考。控制器来源、有限时局域实现和共同D-P时间生成元仍须提供。

## 8. 验证与严格范围

[check.py](check.py)在固定角度$\pi/3$、固定混合概率$1/5$下核四组：

1. 六模式CAR的完整Kraus共轭、外侧混合及有限$Q_\nu$闭合；外侧矩阵故意连接漏出模式与更远模式，避免假设四模式不变。
2. 四模式源算符范数、误差常数饱和、带偶宇称参考的diamond下界，以及原一粒子码全部矩阵单位误差为0。
3. 实际有限码含原模式和漏出模式时的输出见证，核有界源期望与通道对偶一致。
4. 支撑正交隐藏模式的代数作用：原端口全部矩阵单位不变，而隐藏变化范数按两个声明尺度成比例。

代码中的有限生成元矩阵仅独立核交换子和完整CAR恒等式；不是连续Dirac谱截断。实际表示、无界域、同球光滑缩放和任意参考量词由上述解析证明承担。数值复用1121径向诊断的$\sigma_j^2=11/(3R^2)$及选$\beta\chi=\chi$时$\mu_0=m,\sigma_0^2=11/R^2$；该径向函数属$H_0^1\cap H^2$而非$C_c^\infty$，可在相关$H^1$范数内光滑逼近，不增加科学校准。

默认复算只读比较[results.json](results.json)；初次保存才用显式写入选项。本轮没有证明局部应力$T_{\mu\nu}(x)$、任意来源微分、噪声／Hessian、四动量完整谱分布或实际总能源代价受(13)控制。它只关闭本操作合同下**四个完整物理平移生成元的一阶变化误差**。
