# 第571轮研究：连续共同源到精确Gauss初值的受控连接

日期：2026-10-01。接续570。代码joint_gauss_continuum_sampling.py，结果joint_gauss_continuum_sampling_results.json；正式完成以research_round_571_checks.json为准。八组复算，十四式。

## 0. 合并了哪一处条件

568的光滑场能量取样与569的完整规范约束尚未共同成立；570又证明，小Gauss残差不能无条件保证靠近精确物理初值。本轮给出一条限定清楚的正面连接：

**对固定周期盒内、同一非零Higgs的光滑连续Gauss初值，可以构造满足图上完整电弱Gauss的精确初值；在声明的背景条件下，动量修正的真实动能范数为O(epsilon²)，并保留Higgs径向动量以及圆电场的横向／谐和部分。**

这里不增加物质种类或相互作用，使用569同一Hamiltonian及其568尺度族。几何、空间维数、光滑性、场及耦合依然输入；本轮完成的是经典初值与约束的连接，不是完整动力学、量子连续极限或时空生成。数值同时检查d=2、3，也不选择三维。

## 1. 成熟方法与去重

核读[Berchenko-Kogan与Stern，Charge-conserving hybrid methods for the Yang–Mills equations，2003.10054v2](https://arxiv.org/pdf/2003.10054)导言、§3.4—3.5及命题3.6。该文证明其混合有限元变量的逐单元荷守恒，并明确普通全局守恒不能代替局部条件、投影某一个守恒量也可能破坏其它结构。本文沿用“先保持约束再比较来源”的成熟思路；本项目的紧群链路、Higgs荷和圆群不同，不把该文的混合变量或时间结论直接搬过来。

568只给能量一致性，570只给固定图反例和双环构造；549的无规范场引力约束也不能代替本轮非Abelian采样。中点规则、周期Poisson和正交分解均是成熟工具，不作为新基础理论。增量在于它们能在当前共同模型中联立，同时控制全局补荷、物质来源及原电磁自由数据。

|层次|本轮地位|
|---|---|
|认知动机|接入共同约束时保留原有状态信息，补偿资源在整体内部|
|额外输入|固定平直周期盒、d≥2、全局平凡化和单位规范、光滑初值、h≥h_min>0|
|继承物理|569群／荷／正势和同尺度参数；568的v、边势、电磁系数族|
|解析结果|精确Gauss补全、固定背景下二阶动量误差、原自由数据保持|
|数值用途|独立源导数、群矩阵、二维／三维格图、原全能量复核|
|未签收|Higgs零点或非平凡丛、所有背景的一致常数、时间收敛、量子相、引力约束|

## 2. 同一连续来源和符号

令Omega=(R/LZ)^d。所有源在固定盒上光滑、周期，独立于格距。取X=h n，n=(0,1)^T，h≥h_min>0，另有原singlet s。t_a=σ_a/2；连接a_mu已吸收g_w，共轭电密度E_mu相应为568规范化动量Pi_mu/g_w。圆连接a_0和E_0使用Q=6Y；其耦合g_0=g_Y/6，E_0=Pi_Y/g_0。下标0在此标记圆群，**不是时间分量**。颜色初值连接与动量为零；不宣称这一静默选择在未来保持。

完整连续电弱Gauss写成

$$
r=\sum_{\mu=1}^d(\partial_\mu E_\mu-a_\mu\times E_\mu),\quad
q_a(P)=\operatorname{Re}\langle P,it_aX\rangle,\quad
q_0(P)=\operatorname{Re}\langle P,3iX\rangle,\qquad
q(\Pi_X)+r=0,\quad q_0(\Pi_X)+\operatorname{div}E_0=0 .
\tag{1}
$$

负号由W=exp(+i epsilon a·t)和靶动量−Ad(W)^T pi共同固定。本文不与568反Hermitian生成元记法混用。单位规范下3iX=−6it_3X，因此q_0=−6q_3，周期性推出积分r_3为零。径向Pi_X、s及其动量没有被Gauss固定，必须原样保留。

## 3. 对称取样与初始误差

规则周期格N^d、epsilon=L/N，边从x到x+epsilon e_mu；D取源正靶负。边中点m=x+epsilon e_mu/2。采用明确的中点连接及半边框架：

$$
W_e=e^{i\epsilon a_\mu(m)\cdot t},\quad
V_e=e^{i\epsilon a_\mu(m)\cdot t/2},\quad z_e=e^{i\epsilon a_{0,\mu}(m)},\qquad
\pi_e=\epsilon^{d-1}\operatorname{Ad}(V_e)E_\mu(m),\quad
\pi_{0,e}=\epsilon^{d-1}E_{0,\mu}(m),\quad
P_i=\epsilon^d\Pi_X(x_i).
\tag{2}
$$

X_i、s_i及P_s也按同一格取样。W是所选中点近似，不冒充任意变连接的精确路径有序输运；其给定平凡化和单位规范明确输入。生成的精确Gauss数据可以再作任意节点规范变换，但本采样算法本身不宣称对任意连续规范变换精确交织。

定义B_W pi为源pi、靶−Ad(W)^T pi的节点和。由Ad(exp(i theta a·t))=exp(−theta[a×])，每方向节点项等于epsilon^(d−1)乘
exp(−epsilon[a_+×]/2)E_+−exp(+epsilon[a_-×]/2)E_-。这是epsilon的奇对称展开，故

$$
B_W\pi=\epsilon^d r(x_i)+O(\epsilon^{d+2}),\qquad
D^T\pi_0=\epsilon^d\operatorname{div}E_0(x_i)+O(\epsilon^{d+2}),\qquad
G_i^{\rm trial}=O(\epsilon^{d+2}).
\tag{3}
$$

余项在固定光滑源上逐节点一致。下面仍须构造精确修正；式(3)单独不足以签收，570反例继续有效。

## 4. 先关闭全局总荷

令e_3=(0,0,1)。沿所有边求和，输运使弱第三分量不再普通望远镜相消：

$$
c_e=e_3-\operatorname{Ad}(W_e)e_3,\quad
S_\epsilon=\sum_e c_e\cdot\pi_e=\sum_i(B_W\pi)_{i3},\quad
M_\epsilon=\sum_e|c_e|^2,\qquad
\Gamma=\int_\Omega\sum_\mu|a_\mu\times e_3|^2\,dx .
\tag{4}
$$

连续Gauss和周期积分给S_epsilon=O(epsilon²)；周期中点求积及光滑性足够得到此阶。若固定背景Gamma>0，则epsilon^(d−2)M_epsilon=Gamma+O(epsilon²)，细格下正。取

$$
\delta\pi_e=-\frac{S_\epsilon}{M_\epsilon}c_e,\qquad
\widehat\pi=\pi+\delta\pi,\qquad
\sum_i(B_W\widehat\pi)_{i3}=0 .
\tag{5}
$$

这是消去一个明确全局线性泛函的最小普通边动量修正，不是整个Gauss的能量最小投影，也不是实际状态制备。

不能只由delta pi小推其散度小。直接恒等式给

$$
(B_Wc)_i=\sum_\mu\left[
2e_3-\operatorname{Ad}(W_{i,\mu})e_3
-\operatorname{Ad}(W_{i-\mu,\mu})^Te_3\right]=O(\epsilon^2),
\quad
B_W\delta\pi=O(\epsilon^{d+2}).
\tag{6}
$$

其中两个一阶项相消；S/M=O(epsilon^d)，所以这里没有损失一阶导数。

## 5. 用原Higgs补局部荷，保留径向数据

记r_hat=B_W pi_hat，b_i=q(P_i)+r_hat_i。物质生成元的实Gram矩阵为h_i² delta_ab/4，故

$$
\delta P_i=-\frac4{h_i^2}\sum_{a=1}^3 b_i^a\,it_aX_i,\qquad
\widehat P_i=P_i+\delta P_i,\qquad
q(\widehat P_i)+\widehat r_i=0,\quad
\operatorname{Re}\langle X_i,\delta P_i\rangle=0 .
\tag{7}
$$

这一步只修正角动量，保留任意给定的径向动量；同一h_min使逆Gram有界。式(3)、(6)给delta P_i=O(epsilon^(d+2))。新的圆物质荷严格是6r_hat_i3。

## 6. 修正原圆电场，保留横向与谐和数据

由式(5)，rho_i=−6r_hat_i3−(D^T pi_0)_i严格零均值。用周期图Laplace的零均值逆：

$$
\rho=-6\widehat r_3-D^T\pi_0,\qquad
\delta\pi_0=D(D^TD)^+\rho,\qquad
\widehat\pi_0=\pi_0+\delta\pi_0,\qquad
q_0(\widehat P)+D^T\widehat\pi_0=0 .
\tag{8}
$$

delta pi_0属于im D，因而与ker D^T正交：原电场的离散无散部分（含周期谐和）精确保留。不是从零重建一个满足同荷的圆电场。后者会无端删除原物理自由度。

周期格最小非零特征值是4 sin²(pi/N)，故

$$
\|\delta\pi_0\|_2\le
\frac{\|\rho\|_2}{2\sin(\pi/N)}
\le\frac{L}{4\epsilon}\|\rho\|_2,\qquad
\rho_i=O(\epsilon^{d+2}),\quad
\|\delta\pi_0\|_2=O(\epsilon^{d/2+1}) .
\tag{9}
$$

这里的逆谱隙常数基于固定物理盒；不是固定未缩放图的统一常数。

## 7. 共同能量范数及适用边界

沿用568／569的v=epsilon^d、k_h=k_s=epsilon^(d−2)、b_w(epsilon)=b_w(1)epsilon^(2−d)、b_0同尺度，磁势整体系数乘epsilon^(d−4)。定义

$$
\|\delta p\|_{\rm kin,\epsilon}^2=
\sum_i\frac{|\delta P_i|^2}{\epsilon^d}
+2\epsilon^{2-d}\sum_e
\left[b_w(1)|\delta\pi_e|^2+b_0(1)|\delta\pi_{0,e}|^2\right].
\tag{10}
$$

三种修正分别满足Euclidean大小O(epsilon^(d/2+2))、O(epsilon^(d/2+1))、O(epsilon^(d/2+1))。代入正确尺度即得

$$
G_{\rm discrete}(\widehat q,\widehat p)=0,\quad \widehat q=q,\qquad
\|\widehat p-p\|_{\rm kin,\epsilon}\le C\epsilon^2,\qquad
|H_\epsilon(q,\widehat p)-H_\epsilon(q,p)|
\le\|p\|_{\rm kin,\epsilon}\|\delta p\|_{\rm kin,\epsilon}
+\tfrac12\|\delta p\|_{\rm kin,\epsilon}^2=O(\epsilon^2).
\tag{11}
$$

全H所有配置势能不变，动量来源有界；s动量亦不变。由此568给定光滑源的经典能量一致性可以在**精确Gauss初值**上继承，569的磁函数及同尺度Hessian照旧保留。此处第二个O(epsilon²)是补全过程造成的能量变化，不能冒充原整套离散化误差全部已升为二阶。

C依赖L、固定源的导数、h_min、耦合及正Gamma分支的Gamma下界；不承诺所有状态或任意扩大盒的统一常数。若Gamma=0，则每个a_mu严格沿e_3，c_e和S_epsilon恒零；令delta pi=0，式(7)—(11)仍成立。严格中性分支可以连接，不是零隙失败。跨越Gamma趋零的一般背景族则须另外估计，不能用这两个分支宣称统一界。

该动能L²估计足以控制相应有界权重的积分二次动量读数。本文没有证明点态完整Hilbert应力收敛，也没有把本初值送入Einstein初始约束；图—引力连接仍未关闭。

## 8. 非平凡来源与不能丢的自由数据

数值用L=2pi，A=.41、B=.29、u=.38、E*=.17，f=1+u cos x。二维主来源为

$$
a_x=Af\,e_1,\quad a_y=B(1+.2\sin y)e_2,\qquad
E_x=E_*\!\left[f^2-(1+\tfrac32u^2)\right]e_2,\quad
E_y=.11\cos x\,e_2,\qquad
r_3=-Af(E_x)_2 .
\tag{12}
$$

该连接有非零非Abelian曲率，非零规范电场和物质角荷。取圆电场

$$
E_{0x}=6AE_*\left[(2u-\tfrac34u^3)\sin x
+\tfrac34u^2\sin2x+\tfrac1{12}u^3\sin3x\right]+.13,\qquad
E_{0y}=.07\cos x .
\tag{13}
$$

其散度为−6r_3。h=h_star[1+.05 sin(x+y)]，s=s_star[1+.06 cos y]，径向Pi_X=.04 cos(x+y)n；角动量由连续式(1)确定，Pi_s=.025 cos x。另取圆连接a_0x=.03 cos y、a_0y=.08 sin x。三维例增加a_z=.12 cos z e_3、E_z=.09 sin y e_3。原共同势及参数不变；这些光滑状态是指定初值，不是宇宙初态推导。

连续总r_3为零，但对称取样仍留下真实的二阶总荷：

$$
S_\epsilon=
\frac{\epsilon^2 A^3E_*|\Omega|}{24}
\left(2u^2-\frac38u^4\right)+O(\epsilon^4).
\tag{14}
$$

可直接由每条x边的配对−2epsilon^(d−1)sin(epsilon Af/2)(E_x)_2展开；充分细的周期三角求积无低阶混叠。说明补全并非在一个自动满足Gauss的特例上重算零。

反例控制：若从零解圆Poisson代替修正原pi_0，上例会丢掉harmonic .13和横向.07 cos x，二维动能范数损失恒为0.052262895，细化不消失。若把物质动量整个重置为角向补荷，还会丢掉径向部分，损失恒为0.177715318。两者都可能满足Gauss，却不是原来源的受控逼近。

## 9. 八组复算

|诊断|结果|
|---|---|
|连续源的独立中心差分与真实荷|最大残差4.82e-11；不用定义r的表达式单独自证|
|真实SU(2)矩阵与Rodrigues输运|最大差1.34e-15；弱散度密度误差/epsilon²约0.0100—0.0104|
|精确采样，d=2、3共七格|完整Gauss密度残差最大2.57e-16|
|径向、圆横向与谐和保持|径向变化0，横向测试配对3.47e-18，谐和均值6.70e-19|
|严格中性背景，三格|Gamma=0仍精确补全，残差不超过1.78e-16|
|全局总荷、实际gap、Poisson界|二维Gamma_epsilon趋向10.5020021；S/epsilon²趋向0.005414|
|完整569势及真实动能|四格全能量变化均满足式(11)，未把变化称为最小制备功|
|删来源的对照|两种非零范数损失在三格保持不变|

二维N=8、12、18、26时，补全动能范数依次为0.1165515、0.0528808、0.0237215、0.0114137，除epsilon²后为0.18895、0.19289、0.19468、0.19544。三维N=8、12、16的对应比为0.47362、0.48349、0.48703。二阶一般性由前述估计承担，有限样本只交叉复核。

## 10. 对统一条件的实际贡献与下一步

同一既定物质模型的连续约束源、离散Gauss、内部补荷和真实能源现在可在初值层同时成立。没有为补全再增加新物质或新相互作用；所用h非零和背景正Gamma是分支条件，必须保留。该成果排除了“Gauss成立就可以任意重置自由数据”的错误合并方式。

下一步优先核**同一映射能否承接实际演化及共同应力**。先查成熟Hamilton格点规范方法和已有549约束工作，列明能直接复用的定理与剩余模型映射，不继续精修本轮常数。不能从能量一致与精确初始约束直接签收有限时间解收敛；若完整569磁势在消去高能自由度或共同几何中增加独立输入，先登记其具体作用。维数、手征费米子、实际装置、量子连续与引力来源继续开放。


**完整证据：** [代码](joint_gauss_continuum_sampling.py)、[结果](joint_gauss_continuum_sampling_results.json)、[核验](research_round_571_checks.json)、[条件账](unified_physics_condition_ledger_571.md)、[预审](round571_drafts/sampling_pre_review.txt)、[终审](round571_drafts/final_review.txt)。入口与历史稿件保留。
