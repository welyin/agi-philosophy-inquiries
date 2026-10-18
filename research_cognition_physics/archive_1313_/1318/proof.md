# 1318：同一331作用的Dirac谱、跨味流与完整中性交换和式

2026-10-19。[准入](../_admission/round1318_admission.md)。作者证明，状态以独审与验收为准。新增公理0；本轮只签结构层、树级及明示的一圈有限运行合同。

## 1. 净问题与继承范围

本轮采用用户允许的双向研究法：不同家族承担不同规范表示，须共同消去反常，作为“不同职责共同组成一致整体”的物理候选。此动机不推出331群、表示、真空或运行正号。

[Pisano—Pleitez](https://arxiv.org/pdf/hep-ph/9206242)、[Frampton](https://doi.org/10.1103/PhysRevLett.69.2889)、[Long §II](https://arxiv.org/pdf/hep-ph/9504274)的代际消反常以及[331LHN §II](https://arxiv.org/pdf/1206.5779)的重中性伙伴、三标量直接复用。原文的中性菜单及附加对称与本轮不同，不能整包继承。本轮不声称首次发现331、三代计数、味非普适中性流或Schur补。

项目新增的共同判定是：在保留r允许的全部中性质量交叉的同一作用中，给完整Dirac极点及有限界；从同一标量背景和规范动能推出跨夸克部门的流张量与**包含两个有质量中性模的交换权和式**；所有依赖参数的流、质量和混合一起取来源。并给一个有界势的稳定局部真空及非零交叉／相位实例，避免拼接不同父模型。

1060的原SM联合全局反常不自动运输；1303的原三味实际仪器不自动运输；1314的旧六味群不适用于不同规范表示的家族。这里保多Higgs，不证明单Higgs退耦、全UV、完整仪器或现实参数。

## 2. 字段、全部中性质量项与同一作用

取$Q=T_3-T_8/\sqrt3+X$，规范动能标准正号，$D_\mu=\partial_\mu-ig_LW^a_\mu T_a-ig_X X B_\mu$。x=3X只用于整数荷计数，实际耦合仍为$g_XX$，没有额外$\sqrt6$约定。字段为

|字段|$SU(3)_C,SU(3)_L,X$|份数|全左r|
|---|---|---:|---|
|L|$(1,3,-1/3)$|n|i|
|$Q_+$|$(3,3,1/3)$|a|i|
|$Q_-$|$(3,\bar3,0)$|b|i|
|$u^c,d^c$|$(\bar3,1,-2/3),(\bar3,1,1/3)$|各n|i|
|$e^c,N^c$|$(1,1,1),(1,1,0)$|各n|i|
|$F^c$|$(1,1,0)$|n|−i|
|$U^c,D^c$|$(\bar3,1,-2/3),(\bar3,1,1/3)$|a、b|−i|

$a+b=n>0$。$L=(\nu,e,F)_L$；$Q_+=(u,d,U)_L$，$Q_-=(d,-u,D)_L$。三标量$\eta,\chi\sim(1,3,-1/3),\rho\sim(1,3,2/3)$，r分别为−1、+1、−1。规范场r平凡。普通Dirac门户用$\eta,\rho$或其共轭，重门户用$\chi$或其共轭；例如四分量中有$\bar L\eta Y_\nu\nu_R,\bar L\rho Y_e e_R,\bar L\chi Y_FF_R$。夸克三／反三按Long的相应标量共轭收缩，每个允许行的系数一般，故两个普通质量矩阵可任意满秩。

还须保留

$$
\epsilon_{ijk}h_{ab}(L_a^iL_b^j)\rho^k+\mathrm{h.c.},
\quad h^{\mathsf T}=-h,
\qquad N^{c\mathsf T}\mu F^c+\mathrm{h.c.}
\tag{1}
$$

左旋量收缩后两费米双线性在综合指标下对称，结合SU(3)的反对称张量给h反对称。两项的X和r均中性。$N^cN^c,F^cF^c$则被r禁止。允许项分类限这些字段的维数不超过4；更高阶项和其它字段须另登记。本轮不要求每个允许系数非零，但(1)在共同例中两者均非零。

同一父作用含这些规范、标量动能、下节势和全部显示Yukawa／混合；无另设的有效质量或交换系数。取实际标量幅（不是$v/\sqrt2$符号）

$$
\langle\eta\rangle=w_\eta e_1,\quad
\langle\rho\rangle=w_\rho e_2,\quad
\langle\chi\rangle=w_\chi e_3,\qquad w_i>0.
\tag{2}
$$

真空保持$r_{\rm eff}=r\operatorname{diag}(-1,-1,+1)_{SU(3)_L}$。中性全左字段$\Psi_-=(\nu_L,F^c)$和$\Psi_+=(N^c,F_L)$分别取−i、+i。因此二次Weyl系数恰为$\Psi_-^{\mathsf T}\mathcal D\Psi_++\mathrm{h.c.}$，无同相位Majorana块。本文全部谱和物理左流使用

$$
M=\mathcal D^*=\begin{pmatrix}A&B\\C&D\end{pmatrix},
\quad A=w_\eta Y_\nu^{\rm phys},\quad
B=w_\rho H,\ H^{\mathsf T}=-H,\quad C=\mu^{*\mathsf T},\quad D=w_\chi Y_F^{\rm cc},\qquad
Y_F^{\rm cc}=(Y_F^{\rm orig})^{\mathsf T}.
\tag{3}
$$

H、μ的转置及(1)的常数按(3)定义吸收，Weyl系数统一取相应共轭；物理质量为M奇异值，物理左平方核为$MM^\dagger$。这里物理左字段为$(\nu_L,F^c)$，右字段为$(N^c,F_L)^c$，包含对原重F粒子／反粒子标签的重排。不得把$\mathcal D\mathcal D^\dagger$的相位号直接冒充物理左号。带电普通与重夸克在该真空的$r_{\rm eff}$不同，二次轻重混合被禁，普通上／下质量块各为n维。保r不等于全部中性轻重混合为零，也不保证完整连续轻子数对称；下节允许的$(\eta^\dagger\chi)^2$已说明后者不能免费领取。

这里$Y_F^{\rm orig}$专指前述四分量门户$\bar F_LY_F^{\rm orig}F_R$，$Y_F^{\rm cc}$为重排后的电荷共轭Dirac图册；两者不能同名误认。代码以独立物理块C和H为坐标，原始μ及h通过这项可逆线性字典一起变动，未冻结原始来源。实对角D例中转置不改变数值。

## 3. 有界势和局部真空确实非空

令$\Phi=(\eta,\rho,\chi)$，$x_i=\phi_i^\dagger\phi_i$，$z_{ij}=\phi_i^\dagger\phi_j$，$W=w_\eta w_\rho w_\chi$。给一个允许的可重整化势

$$
V=\sum_i\lambda_i(x_i-w_i^2)^2+\sum_{i<j}k_{ij}|z_{ij}|^2
+2\epsilon\operatorname{Re}(z_{\eta\chi}^2)
+\sum_i\frac{fW}{w_i^2}(x_i-w_i^2)-2f(\operatorname{Re}\det\Phi-W).
\tag{4}
$$

所有项规范及r不变；ε实，f>0，$\lambda_i>0,k_{ij}>0,k_{\eta\chi}>2|\epsilon|$。正四次支配三次，V有下界。(2)为驻点，数值系数是选定的输入，不是认知最小化原则。稳定性只要求这个真空的有限邻域，未证明其为唯一全局最低点。

在该满秩真空附近，对八个已破缺规范方向取局部横截面，可用Gram矩阵$G=\Phi^\dagger\Phi$的九实分量及$\theta=\arg\det\Phi$参数化十个物理涨落。部分$z_{ij}$仍受未破缺电磁U(1)作用，不能称它们全部为完整规范商的不变量；该U(1)在真空处不另产生切向零模。恒有$\operatorname{Re}\det\Phi=\sqrt{\det G}\cos\theta$。x方向Hessian为

$$
(H_x)_{ii}=2\lambda_i+\frac{fW}{2w_i^4},\qquad
(H_x)_{ij}=-\frac{fW}{2w_i^2w_j^2}\ (i\ne j).
\tag{5}
$$

$z_{ij}$的二次系数为$k_{ij}+fW/(w_i^2w_j^2)$，ηχ两实方向再加减$2\epsilon$；θ二次系数为fW。只要(5)正定，全部十物理方向严格正，余八方向为规范轨道。共同例取$w=(1,1,100)$、$\lambda_i=.2,k_{ij}=.3,f=.01,\epsilon=.02$；(5)严格对角占优，故该真空严格局部稳定。这里并未删去其余物理标量以假称只有一个Higgs。

式(4)把二次势系数写成所选真空的参数化。取真正的作用来源时，应先恢复独立二次、四次、三次系数再变分，不能把驻点条件代回后删源。下文w导数是已声明背景响应，不单独称为额外能量或引力源。

## 4. 代际与保护的条件计数

全左计数给$\mathcal A_{L^3}=2(3a-n),\mathcal A_{L^2x}=3a-n,\mathcal A_{x^3}=6(n-3a)$；$\mathcal A_{C^3},\mathcal A_{C^2x},\mathcal A_{\mathrm{grav}^2x}$恒为零。故n=3a、b=2a。全新彩色阈值之上有3n个Dirac夸克，色一圈系数$b_{0C}=11-2n$。额外要求该有限全活跃弱耦合窗中$b_{0C}>0$，即得n=3、a=1、b=2。该正号不等于一般有限窗口可控，不声称无截断UV完成。

r的两个SU(3)单位瞬子指数均为4n，Jacobian为1；$N_+-N_-=16n,\sum sx^2=24n,\sum s^2x=0$（s=±1）。纯Spin–Z4指标为零，但这些检查未证明任意联合全局背景异常平凡。共同中心核为高层Z3；选忠实商可对应低层G6，群核不替代全局异常定理。本轮只领取上述明示范围的资格。

本类别内的n=3结论不排除原SM字段类别的四／五代共同模型；类别选择、331表示和运行正号仍是独立输入。

## 5. 完整2n Dirac极点与有限近似

取同一M，写

$$
MM^\dagger=\begin{pmatrix}E&F\\F^\dagger&G\end{pmatrix},\quad
E=AA^\dagger+BB^\dagger,\ F=AC^\dagger+BD^\dagger,\ G=CC^\dagger+DD^\dagger.
\tag{6}
$$

设$\gamma^2=\lambda_{\min}(G)>0,z_*=\|E\|<\gamma^2$，M满秩。主子矩阵交错给恰有n个低极点$0<z_j\le z_*$，其余不低于$\gamma^2$。低极点**精确**满足

$$
\det S(z)=0,\quad S(z)=E-zI-F(G-zI)^{-1}F^\dagger.
\tag{7}
$$

令$K_0=E-FG^{-1}F^\dagger>0,C_0=FG^{-2}F^\dagger,\delta=z_*/\gamma^2<1$，则

$$
S(z)=K_0-z\{I+FG^{-1}(G-zI)^{-1}F^\dagger\},
$$

$$
\boxed{\lambda_j(K_0,I+C_0/(1-\delta))\le z_j\le\lambda_j(K_0,I+C_0).}
\tag{8}
$$

证明：在$0\le z\le z_*$，花括号中新增项介于$C_0$与$C_0/(1-\delta)$；$S'(z)=-I-F(G-zI)^{-2}F^\dagger\prec0$。比较零点惯性给逐序界。相应任意$|z|<\gamma^2$的重核相对零点余量至多$|z|/(\gamma^2-|z|)$。这不是将$A-BD^{-1}C$奇异值直接当完整物理质量，也不把轻极点附近的传播误差无界地遗漏。

完整源由$\partial_x(MM^\dagger)=M_xM^\dagger+MM_x^\dagger$给出。固定z时(7)的导数为

$$
S_x=E_x-F_xRF^\dagger-FRF_x^\dagger+FRG_xRF^\dagger,\quad R=(G-zI)^{-1}.
\tag{9}
$$

不得在变化h、μ、Y或背景时冻结这些依赖。代码同时核全谱Hellmann–Feynman和消元来源。

## 6. 对角中性规范块的完整混合与跨味流

331另外含由$W^4,W^5$组成的电中性非厄米场$X^0$。本节3×3块仅指Cartan–U(1)对角中性部门，不声称整个331只有两个有质量中性粒子。各标量真空仅沿一个正交基向量，故完整8+1规范质量核中，$W^{4,5}$与本节块无交叉；二者的质量平方都是$g_L^2(w_\eta^2+w_\chi^2)/2$。生成元只连接第一与第三分量，对夸克即普通—重伙伴，其普通字段子空间内的双线性投影为零。残余r同时禁止普通—重夸克质量混合。因此在下节所声明的普通u/d横向非普适树级二源响应中，$X^0$不额外贡献。其它外腿、圈阶及一般非守恒源不在这个结论内。

按基$(W^3,W^8,B_X)$定义三个真实标量荷向量

$$
v_\eta=(g_L/2,g_L/(2\sqrt3),-g_X/3),\quad
v_\rho=(-g_L/2,g_L/(2\sqrt3),2g_X/3),\quad
v_\chi=(0,-g_L/\sqrt3,-g_X/3).
\tag{10}
$$

它们之和为0。标量动能的完整中性质量平方为

$$
M_V^2=2\sum_i w_i^2v_iv_i^{\mathsf T}.
\tag{11}
$$

系数2来自$\frac12 A_\mu^{\mathsf T}M_V^2A^\mu$且使用实际幅w。对$g_L,g_X>0$秩为2；光子方向为$e(1/g_L,-1/(\sqrt3g_L),1/g_X)$，其中$e^{-2}=4/(3g_L^2)+1/g_X^2$。其余两个正本征值$m_1^2,m_2^2$及正规本征向量$o_1,o_2$来自完整(11)，没有先删去Z混合。

设P为夸克三重态家族的秩a投影。普通u、d的左耦合向量分别为

$$
v_u^-=(g_L/2,-g_L/(2\sqrt3),0),\quad
v_d^-=(-g_L/2,-g_L/(2\sqrt3),0),\quad
d=(0,g_L/\sqrt3,g_X/3)=-v_\chi.
$$

因此每个真实中性模s在物理质量基的耦合为

$$
\Gamma_u^{(s)}=c_u^{(s)}I+\kappa_sU_u^\dagger PU_u,\quad
\Gamma_d^{(s)}=c_d^{(s)}I+\kappa_sU_d^\dagger PU_d,
\quad c_f^{(s)}=v_f^-\cdot o_s,\ \kappa_s=d\cdot o_s.
\tag{12}
$$

右手耦合按X普适。光子有$\kappa_\gamma=0,c_u=2e/3,c_d=-e/3$。设$V=U_u^\dagger U_d$、$P_f=U_f^\dagger PU_f$，则对**两个有质量模同时**有

$$
P_d=V^\dagger P_uV,\quad P_f^2=P_f,\ \operatorname{rank}P_f=a.
\tag{13}
$$

n=3条件将其变为rank-one。这不是把SU(3)活动像改称原1312的单参考SU(2)：高层采用了不同字段、多个标量参考和不同家族职责。

## 7. 同一真空的完整交换和式与有限窗口

将两部门的**非普适部分**组成同一流$j_P=j_{P_u}+j_{P_d}$。两个有质量中性模都耦合为$\kappa_s j_P$。本节核的是树级横向源响应／横向交换系数；一般不同质量的手征味流不守恒，实际完整散射还含纵向、Goldstone与普适流，未把这一系数冒充整个FCNC散射率。

定义类时解析变量z和横向权

$$
\mathcal F(z)=\sum_{s=1}^2\frac{\kappa_s^2}{m_s^2-z},\qquad
\mathcal F_0=d^{\mathsf T}(M_V^2)^+d.
\tag{14}
$$

在所声明横向源的领先响应中，零点非普适局部项按约定为$-\frac12\mathcal F_0 j_{P\mu}j_P^\mu$；特别是uu、ud、dd三块的味系数不能独立输入，它们分别是同一$\mathcal F_0$乘$P_u\otimes P_u,P_u\otimes P_d,P_d\otimes P_d$。这不是删去一般非守恒流的纵向与Goldstone项后仍声称完整散射。

**完整和式：**

$$
\boxed{\mathcal F_0=\frac1{2\left(w_\chi^2+\dfrac{w_\eta^2w_\rho^2}{w_\eta^2+w_\rho^2}\right)}.}
\tag{15}
$$

证明：令$B_v=(v_\eta,v_\rho)$满列秩，$\alpha=(1,1)^{\mathsf T}$，$a_0=w_\eta^2,b_0=w_\rho^2,c_0=w_\chi^2$。则$d=B_v\alpha$，

$$
M_V^2=2B_v\begin{pmatrix}a_0+c_0&c_0\\c_0&b_0+c_0\end{pmatrix}B_v^{\mathsf T}.
$$

限制到其秩2像后求逆，(14)为$\frac12\alpha^{\mathsf T}D_v^{-1}\alpha=(a_0+b_0)/[2(a_0b_0+c_0(a_0+b_0))]$。因此$g_L,g_X$在**两模总和**中消去，不能冻结混合或只保重Z后仍领取此式。

同一完整8+1核还给$m_W^2=g_L^2(a_0+b_0)/2$、$m_X^2=g_L^2(a_0+c_0)/2$、$m_Y^2=g_L^2(b_0+c_0)/2$。因此(15)也可以写成同作用的质量—耦合字典

$$
\mathcal F_0=\frac{g_L^2}{2m_X^2+2m_Y^2-m_W^2-(m_X^2-m_Y^2)^2/m_W^2}.
$$

这只是原和式的等价表达，未新增独立预测或圈级物理质量认证。

对任意$|z|<m_{\min}^2=\min(m_1^2,m_2^2)$，正残数直接给

$$
|\mathcal F(z)-\mathcal F_0|\le\frac{|z|}{m_{\min}^2-|z|}\mathcal F_0.
\tag{16}
$$

这个界按较轻中性质量，不能换成重Z质量。在只积分重Z时应单独保留轻Z传播。$\mathcal F_0$也不是固定正动量下的成功率或资源证书；当耦合趋零时质量窗口随之收缩，零点和式不提供统一有限分辨任务。

令$K=(M_V^2)^+$，d始终在其像内。全部参数／背景来源满足

$$
\partial_x\mathcal F_0=2d_x^{\mathsf T}Kd-d^{\mathsf T}K(M_V^2)_xKd.
\tag{17}
$$

伪逆微分的零模项被两侧d消去；这不是假定零模及混合与参数无关。特别地$\partial_{g_L}\mathcal F_0=\partial_{g_X}\mathcal F_0=0$，而

$$
\partial_{w_\chi}\mathcal F_0=-4w_\chi\mathcal F_0^2,\quad
\partial_{w_\eta}\mathcal F_0=-4w_\eta\left(\frac{w_\rho^2}{w_\eta^2+w_\rho^2}\right)^2\mathcal F_0^2
\tag{18}
$$

ρ式交换η、ρ。若连Yukawa一起变，(13)的$P_f,V$也随之变化，完整味系数导数还须含它们。程序核有限非简并矩阵上的这种变化，不只核一个冻结P的标量源。

## 8. 同父非空点、相位与剩余输入

使用第3节同一w与势，$g_L=.65,g_X=.4$。令U为三阶单位DFT，取物理块

$$
A=U\operatorname{diag}(.01,.02,.03)U^\dagger,\quad
B=10^{-4}\begin{pmatrix}0&1&0\\-1&0&1\\0&-1&0\end{pmatrix},\quad C=10^{-4}I,\quad D=\operatorname{diag}(10,12,15).
\tag{19}
$$

带电轻子取正、非简并对角质量；普通夸克取$M_u=\operatorname{diag}(.10,.20,.30),M_d=U\operatorname{diag}(.15,.25,.35)$，重夸克质量10、12、15。所有Yukawa由同一w登记，皆有限。高阈值色运行窗与低交换窗口是该作用的不同声明范围，不把两处截断混用。

(19)相对$\operatorname{diag}(A,D)$的扰动范数为$\sqrt2\times10^{-4}$，满秩与重隙有严格裕量。完整轻谱和有限界见结果。物理左低模投影给$K_{\rm low}=W_{\nu,\ell}\operatorname{diag}(m_\ell^2)W_{\nu,\ell}^\dagger$，其三元乘积虚部约$1.15470051\times10^{-11}\ne0$。这是非零相位结构证据；$W_{\nu,\ell}$一般非酉，未因此自动认证原三味仪器。

零交叉点有显式非零DFT相位；正质量隙、势的物理Hessian正性、非零流／相位均在小邻域保持。故共同实现不是孤立的精细零交叉点。相同w下改变允许H、μ、Y可以改变完整质量和味取向；相同表示下改变wχ改变(15)。三代没有决定现实质量、混合、尺度或规范耦合。

四删弱直接保留：n=6、a=2、b=4仍消局部反常但全活跃$b_{0C}=-1$；若六代新增彩色态不活跃，只计12普通Dirac夸克则$b_{0C}=3$；删r允许单态同类Majorana；保r但改H、μ改变完整谱。它们只判相应额外条件，不宣称六代在任何有限物理类都无效。

## 9. 核验及停止

程序核精确荷计数、残余r、势的局部正性、完整中性质量及(8)、物理CP、完整中性规范核、两部门流和(15)—(18)，以及父参数来源与删弱。独审另用独立实现。成熟计数／一般矩阵工具计入复用，本轮只领取所列同对象连接。

保持普通局部／瞬子范围；完整331联合全局背景未证明。多Higgs与新重物质保留，未称原单HiggsSM精确恢复；横向核没有代替全部实际散射；也未把用户的全部认知公理由此实现。完成后停止该候选参数扫描，回到表示／运行／真空等剩余角色的独立来源；阶段目标仍未结项。
