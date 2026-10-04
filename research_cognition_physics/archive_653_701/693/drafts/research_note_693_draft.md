# 第693轮：完整规范平均的近零谱控制与原物理来源误差

日期：2026-10-02。接[692](../../research_note_692.md)和[已执行有限表示筛查](finite_character_entry.md)。[代码](../joint_gauge_sublevel_control.py)、[结果](../joint_gauge_sublevel_control_results.json)、[核验](../research_round_693_checks.json)、[条件账](../unified_physics_condition_ledger_693.md)。两组新核验、十六式；主代理审查，无独立代理审查。

## 1. 接续问题与新旧边界

上一目标轮完成692、核验发布并执行693入口，属于进展。本轮先读导航、最新结果与入口核验，未发现运行中的Python。回查614、669、673、676—677、686及相关谱隙／规范截断历史。

677／686已经证明原全部平均的支配收敛；673已经证明Wilson零集Haar零测度。**本轮新增的是对所有原空间规范背景统一的近零谱测度幂上界，并把它接成全部指定物理来源的显式全平均误差式。** 不增加硬可容许域，不删近零谱配置，不把“几乎处处收敛”重复当新结果。

该界极松，尚不能提供实际正负号证书。它关闭“是否还须额外假定小谱概率衰减”的逻辑缺口，没有关闭原完整动态RP、非零归一、H_F、共同连续或量子GR。

|层次|地位|
|---|---|
|认知动机|统一模型的表示近似须保原来源、全配置与可审计误差|
|继承输入|原2空间点×2AP时间片、两传播方向、16内部通道、m0=1、原H_b与全平均|
|成熟工具|Cayley下的Haar密度、Brudnyi–Ganzburg多项式次水平不等式|
|解析增量|单位时间边界的统一见证隙；原表示次数与全Haar近零谱界；同一来源的热平均误差|
|验证|原完整群空间链路、原表示清分母、多项式次数账、512维B来源及奇异Pfaffian控制|
|保留限制|固定有限盒；热矩常数未数值求出；界不实用，不代表RP或实际宇宙尺度|

空间按[679旧接口表](../../research_note_679.md)复用382—386、425、522—523；384消去的输入不恢复，386与425保持替代路线，523共同实现不重开。以下谱参数不是物理时间或实际空间位移，E仍非s记录。

## 2. 673的非零点可给统一下界

先固定任意原空间链路，但把两个时间边界设为单位。沿673，反周期两片移位T0酉且反厄米；空间Wilson部分可写X_s=I-U_s，其中U_s由gamma1的两个正交投影分别乘T_s、T_s†，故U_s酉。对任意单位psi：

$$
X=X_s+\gamma_0T_0,\quad W_s=\operatorname{Re}X_s\succeq0,
\qquad \|X_s\psi\|^2=2\langle\psi,W_s\psi\rangle.
\tag{1}
$$

令epsilon=||X psi||。由于时间项反厄米，Re〈psi,Xpsi〉=〈W_s〉≤epsilon；三角不等式给

$$
1=\|\gamma_0T_0\psi\|\le\epsilon+\sqrt{2\epsilon},
\qquad \sigma_{\min}(H(U_s,I,I))\ge g_*:=2-\sqrt3>0.
\tag{2}
$$

这是原单位边界上的见证下界，不声称其它边界没有零谱。原H维数n=256、全域范数h≤3，故

$$
|\det H(U_s,I,I)|\ge g_*^{256},\qquad
\|H(U_s,a,b)\|\le3
\quad\text{对全部原 }U_s,a,b.
\tag{3}
$$

因此下一节的常数不再依赖任选空间配置上的未知最小行列式。

## 3. 原整个商群的有理坐标，不以U1代替全群

对Hermitian A_j，先构造Haar U(j)的Cayley变量，再去掉行列式得到SU(j)：

$$
Q_j=(I+iA_j)(I-iA_j)^{-1},\qquad
C_j=Q_j\operatorname{diag}((\det Q_j)^{-1},1,\ldots,1),
\qquad z=(1+it)/(1-it).
\tag{4}
$$

去行列式映射对左SU(j)作用等变，所以Haar U(j)被送到Haar SU(j)。独立取j=3、2和Cauchy t，再投到原Z6商群，恰得原归一Haar；不是修改规范群。每个群元使用9+4+1=14个实坐标，允许表示冗余。

采用Hermitian矩阵对角实分量、上三角实虚分量的Lebesgue测度，Cayley密度与尾界为

$$
\rho_j(A)=c_j\det(I+A^2)^{-j},\quad
c_j=\frac{2^{j(j-1)}\prod_{k=1}^{j-1}k!}{\pi^{j(j+1)/2}}<1\ (j=1,2,3),
\qquad \Pr(\|A_j\|>R)\le\frac{2j}{\pi R}.
\tag{5}
$$

密度及归一来自[Assiotis等，式(1)—(5)的s=0情形](https://www.pure.ed.ac.uk/ws/portalfiles/portal/213863459/2009.04752.pdf)。尾界另由Haar特征角单点密度均匀和联合界推出：|tan(theta/2)|>R的单角概率为2 arctan(1/R)/pi。没有假定特征角独立。

记delta_j=det(I+A_j²)、delta_1=1+t²。每个SU(j)矩阵元及其共轭可用正分母delta_j，分子次数≤2j。614的原模块依次为z C⊗W、z^-4 C*、z² C*、z^-3 W、z^6、1，所以统一清分母只需

$$
d_g=\delta_3\delta_2\delta_1^6>0,\qquad
\deg d_g=22,\qquad \deg(d_gR(g)_{uv})\le22.
\tag{6}
$$

R†也用同一实正分母，原固定内部字典J不改变次数。本项用原具体荷表，不采用更大的虚拟表示。

## 4. 全双Haar近零谱概率的统一幂界

当前共有K=4个时间边界群元，故Cayley实坐标数M=56。固定任意空间链路，令d为这四个d_g之积，则

$$
p(\xi)=d(\xi)^n\det H(U_s,a(\xi),b(\xi))
\in\mathbb R[\xi_1,\ldots,\xi_M],\qquad
\deg p\le D:=22Kn=22528,\quad |p(0)|\ge g_*^n.
\tag{7}
$$

p实值由H厄米及实正分母保证；空间链路只改变其系数，不改变次数或原点下界。在坐标立方体[-R,R]^M，||A_j||²≤2j²R²，因而d^n≤(1+18R²)^(D/2)。若原谱隙g<delta，则

$$
|p|\le3^{n-1}\delta\,(1+18R^2)^{D/2}.
\tag{8}
$$

多变量Remez不等式的一个直接推论是：对次数≤D的实多项式及立方体Q，

$$
\frac{|\{\xi\in Q:|p(\xi)|\le t\}|}{|Q|}
\le4M\left(\frac{t}{\sup_Q|p|}\right)^{1/D}.
\tag{9}
$$

使用[Ganzburg，定理1.2(a)](https://www.impan.pl/shop/publication/transaction/download/product/92006)，以Chebyshev上界和1-(1-lambda)^(1/M)≥lambda/M得到式(9)。仅借这个多项式工具，不把文献当本模型的谱定理。

Cayley联合密度≤1，立方体体积(2R)^M；落到盒外的概率由四组U3、U2、U1尾界控制。式(3)、(8)、(9)给对所有原空间背景一致的归一Haar界

$$
\Pr_{a,b}\{g<\delta\}\le\mathcal P(\delta;R):=
\min\!\left\{1,\frac{12K}{\pi R}
+4M(2R)^M\sqrt{1+18R^2}
\left(\frac{3^{n-1}\delta}{g_*^n}\right)^{1/D}\right\}.
\tag{10}
$$

对0<delta≤1选R=delta^-alpha≥1，得到简明的显式模量

$$
\Pr_{a,b}\{g<\delta\}\le P_*(\delta):=\min(1,C_*\delta^\alpha),\qquad
\alpha=\frac1{D(M+2)}=\frac1{1306624},\quad
C_*:=\frac{12K}{\pi}+4M2^M\sqrt{19}
\left(\frac{3^{n-1}}{g_*^n}\right)^{1/D}.
\tag{11}
$$

这一幂次极弱，但常数和量词都明确。673的零测度是delta趋零的直接后果；当前新结果是有限delta的统一上界。delta只是证明中拆分积分域的参数，两部分都保留，没有施加物理截止。

## 5. 把这个界接到原来源，而非只控制谱矩阵

对677的2k个源列Z，增广反对称矩阵尺寸2m=2(n+k)。复反对称矩阵的酉合同标准形给Pfaffian微分界；沿两矩阵的线段积分，得

$$
|\operatorname{Pf}A-\operatorname{Pf}B|
\le m\|A-B\|\max(\|A\|,\|B\|)^{m-1},\qquad
|C_{\varepsilon,Z}-C_{\varepsilon',Z}|\le B_Z\|\varepsilon-\varepsilon'\|.
\tag{12}
$$

具体取c=(1+sqrt5)/2、a_Phi=sqrt5/4，令R_Z=2+|lambda|c²||P_phi||+c||Z||，L_Z=1/2+2a_Phi |lambda|c||P_phi||+a_Phi||Z||，则B_Z=m L_Z R_Z^(m-1)。标准形中每个微分项为一个扰动二阶元乘其它m-1个奇异值，故这个界在奇异端点也成立，不借逆传播子。

来源为有限多项式或固定有限来源族时，逐项取上述界并求有限和。B_Z由原w(q_i)+w(q_j)的有限次多项式支配。记原热核对角为k_tau(q,q)，则669已有热矩保证

$$
\mathcal B_Z:=\int B_Z(q_i,q_j)k_\tau(q_i,q_i)k_\tau(q_j,q_j)
d\mu(q_i)d\mu(q_j)<\infty.
\tag{13}
$$

这里可用B_Z≤A(w_i+w_j)^p及669的M_p、M_0写成2^p A M_p M_0的明确上界。没有把这些热矩的“已知有限”冒称为“已经数值求出”。

复用677的实际调节epsilon_(a,L)，在g≥delta且delta_a<delta时，有

$$
\|\varepsilon_{a,L}-\operatorname{sign}H\|
\le e_{a,L}(\delta):=2e^{-aL(\delta-\delta_a)}+
\frac{\delta_a}{\delta-\delta_a},\qquad
\delta_a=\frac{9a}{2(1-3a/2)},\quad a<2/3.
\tag{14}
$$

坏谱域上两Hermitian收缩的差范数≤2；原热核乘积受两对角核乘积控制。式(11)—(14)因此给完整原控制测度dvarpi_tau上的误差

$$
\boxed{\displaystyle
\int|C_{a,L,Z}-C_{\infty,Z}|d\varpi_\tau
\le\mathcal B_Z\,[e_{a,L}(\delta)+2P_*(\delta)] .}
\tag{15}
$$

保留S9、全空间字段及双Haar平均，包含零权重来源，未改成只有电Casimir热核。固定盒、固定正tau及固定来源的677收敛现在有明确上界；仍不提供体积／格距／几何一致性或指定配分函数非零性。

例如沿原a趋零、aL趋无穷的任意序列，可取delta=sqrt(delta_a)+(aL)^(-1/2)。足够后delta_a<delta且式(15)右侧趋零。不需要把La²趋零重新加进原精确调节合同。

## 6. 与692的九点规则怎样共同使用

692只对精确投影的原泛函证明九点身份。软epsilon一般不幂等，不能对其自动套用675解耦或九点次数界。

但可以明确另外定义F_(a,L)^(9)：在692的同一九点正权辅助和中，逐项用677软来源近似。因为式(12)—(15)对单位E统一，而九点权归一，故对Gauss不变、E独立物理来源F，

$$
|\mathfrak I_{a,L}^{(9)}(F)-\mathfrak I_\infty(F)|
\le\mathcal B_F[e_{a,L}(\delta)+2P_*(\delta)].
\tag{16}
$$

这是九点调节表达的误差，不宣称有限a,L时它等于原连续S9软泛函。两种表达都趋向同一精确原目标；完整规范／玻色积分仍要保留。正权及收敛不证明目标正。

## 7. 复算及为什么还不能判定符号

第一组把三个幅度的原完整SU3／SU2／U1空间链路接到单位时间边界，核式(1)—(3)；数值不承担全背景证明。另按原左手各荷模块直接清分母，核同一16维表示和固定内部字典，保存22次上界及Cayley密度常数。

第二组在677原非平坦背景上计算512维B的九项六场来源，检验实际软／精确来源差服从式(12)。用另两份复反对称控制矩阵（含奇异端点）交叉检验Pfaffian界。保存模量的对数值，防止超小delta浮点下溢被误报为精确零。

本例式(11)要保证P_*≤1/2，竟要求log10(delta)≤约-26341846.4。它是严格但非常保守的控制；以此直接判定675原积分的小正负值不可行，且热矩系数尚未求值。**不继续靠缩小a、增加L或优化一两个常数凑研究轮次。**

本轮实际关闭的是近零谱衰减需不需要另加假设，以及该控制能否与同一物理来源／全热平均相接；答案分别为“不需要另加”和“可以，固定有限盒”。原动态RP、H_F、连续物质与量子几何仍未统一。

接[694](../../694/drafts/STATUS.md)：优先检验原来源自身在临界谱附近的消去／抵消是否比裸行列式次水平更强，或寻找能直接控制物理二次型的结构证书。须明确这些性质是否由当前原表示产生，不能偷加横截性、谱密度有界或固定拓扑假设。
