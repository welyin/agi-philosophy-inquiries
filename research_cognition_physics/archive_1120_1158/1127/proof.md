# 1127证明：外给平滑脉冲下的逐点材料装载

2026-10-11。接[1126](../1126/proof.md)的区域CAR门权限，本轮给出一个**单端、有限窗、逐点局域材料交换**实例：相同自由Dirac谱的两种实际材料场，经紧时空支撑的外给时间分量耦合，精确交接完整双轨码与任意参考；正常Fock过程保留全部余部，并有有限末端平均场能量。

这不是空间分离通信整链，也没有实现正常输入的制备仪器、外脉冲来源或自主开关。新增认知公理0、科学校准0，阶段仍active。成熟外场与CAR工具直接复用；净连接及未闭合项见[selection_audit.md](selection_audit.md)、[NEXT.md](NEXT.md)。

## 1. 实际材料、有限窗与点作用

取自然单位，$m>0$，

$$
\mathcal K=L^2(\mathbb R^3,\mathbb C^4)\otimes\mathbb C^2_{AB}\otimes\mathbb C^2_{rail},\quad
h=(\boldsymbol\alpha\cdot\boldsymbol p+\beta m)\otimes I_{AB,rail},\quad
\omega=|h|=\sqrt{|p|^2+m^2},\quad\epsilon=h/\omega,\quad S=(I-\epsilon)/2.
\tag{1}
$$

A、B是两个可区分的实际材料场，各有两条rail；它们均携带(1)的正谱粒子／负谱空穴能源，B不是1126的零能静止辅助。本模型明示这些材料标签及其可访问偶观测量；没有另外施加Gauss商，也没有声称已构造电磁场或某个内部规范群。

选$\tau,r>0$、归一$f\in C_c^\infty(B_r,\mathbb C^4)$，再取实函数

$$
\chi\in C_c^\infty((0,\tau)),\quad \lambda\int_0^\tau\chi(t)dt=\pi/2,
\qquad\eta\in C_c^\infty(\mathbb R^3),\quad \eta=1\text{ 于 }B_{r+\tau}.
\tag{2}
$$

可先取非负非零光滑bump $b$，令$\chi=b/\int b$、$\lambda=\pi/2$；积分归一是解析定义，无需把数值求积充作面积证明。空间平台可在$B_{r+\tau+\delta}$内平滑降至0，$\delta>0$固定有限。

声明一个共同惯性钟与位置字典，给出单粒子局部方程

$$
i\partial_t\psi=(h+V(t))\psi,\qquad
V(t,x)=\lambda\chi(t)\eta(x)I_{spin}\otimes\sigma_x\otimes I_{rail}.
\tag{3}
$$

此处是时间分量矩阵外势：Hamiltonian中的旋量因子为$I_{spin}$，不是Lorentz标量质量混合$\beta\otimes\sigma_x$。有界光滑自伴势给出唯一酉传播$U(t,s)$，保持Sobolev正则性；Dirac主部未变，传播速度界为1。场乘积密度只是(3)的局域耦合记法，**不据此把裸瞬时正常序密度的积分自动定义为真空域上的Fock生成元**。完整Fock实施在§3另证。

## 2. 因果平台内的精确交换

令$u_t=e^{-iht}f$（此处$h$只作用旋量空间）、$\theta(t)=\lambda\int_0^t\chi(s)ds$。有限传播给$\operatorname{supp}u_t\subset B_{r+t}$。于是

$$
\psi_i(t)=u_t\otimes\bigl(\cos\theta(t)e_A-i\sin\theta(t)e_B\bigr)\otimes e_i
\tag{4}
$$

在其完整支集上满足$\eta\psi_i=\psi_i$。直接代入(3)，(4)满足同一初值，唯一性遂给

$$
\boxed{U(\tau,0)F_i=-iG_i,\quad
F_i=f\otimes e_A\otimes e_i,\quad
G_i=u_\tau\otimes e_B\otimes e_i,\quad
\operatorname{supp}G_i\subset B_{r+\tau}.}
\tag{5}
$$

这只是指定两模式上的精确等式。空间平台外仍存在真实演化，$U$在正交余部不被设为恒等。整个作用区有限，窗口长度$\tau>0$；B的目标偶矩阵单位由该有限球内的实际B场给出。A、B在同一端点区域共存，(5)没有声称把载荷传到空间分离的接收端。后续接1126仍须匹配模式、等待、提取和资源。

## 3. 两谱、标准正常Fock与实施相位

对$\sigma_x$作一次恒定酉对角化，(3)成为两个独立Dirac分支，外势为$A_0^{\pm}=\pm\lambda\chi\eta$，空间分量$\boldsymbol A^{\pm}=0$；两rail复制该过程。这是调用外场定理的数学分解，不另赋物理规范权限。

[Deckert—Merkl，Theorem 1.2、1.5及Corollary 1.4](https://arxiv.org/html/1505.06039)给出：实施条件是相应两谱非对角块为Hilbert–Schmidt；等时超平面的极化类只取决于空间外势。故本例每一时间的标准自由极化类可用，存在标准正能Fock上的酉lift $\widehat U$，满足

$$
\widehat U\,c(v)\widehat U^\dagger=c(Uv),\qquad
\widehat U\,c^\dagger(v)\widehat U^\dagger=c^\dagger(Uv),\qquad
\{c(v),c^\dagger(w)\}=\langle v,w\rangle.
\tag{6}
$$

$c$反线性、$c^\dagger$线性；$c(v)$包括正谱粒子湮灭和负谱空穴产生。真空协方差为$S$，不把紧支模式投到正频部分，不把真空局部边缘当空模式。完整lift保留所有Bogoliubov块与场余部。一个固定外势的lift只差整体相位，该共同相位不改变本轮任何$Q/R$通道；但外势变分、电流／应力或量子控制支路比较会涉及其相位导数，不能沿用此处的相位无关性删除来源账。[同文§1.2](https://arxiv.org/html/1505.06039#S1.SS2)

## 4. 正常源码准备与全部未知输入

设$s=\langle f,P_-f\rangle<1$，$n_i=c^\dagger(F_i)c(F_i)$，

$$
e_{empty}=(I-n_1)(I-n_2),\qquad
p_{empty}=\|e_{empty}\Omega\|^2=(1-s)^2>0,\qquad
\Omega_0=e_{empty}\Omega/\sqrt{p_{empty}},\qquad
J_i=c^\dagger(F_i)\Omega_0.
\tag{7}
$$

两rail的原真空独立，给出所示概率。非空例可取实归一$\phi\in C_c^\infty(B_r)$和$\alpha_z\xi=\xi$的归一旋量，$f=\phi\xi$。反对易关系使$\xi^\dagger\beta\xi=\xi^\dagger\alpha_{x,y}\xi=0$，实$\phi$的动量密度在$p\mapsto-p$下偶，故$s=1/2$、$p_{empty}=1/4$；两谱确实都保留。

$c(F_i)\Omega_0=0$，CAR给$\langle J_i,J_j\rangle=\delta_{ij}$。$J:\mathbb C^2\to\mathcal F$是正常等距编码；两列同为奇宇称，逻辑矩阵代数由偶元实现，任意逻辑相干及物理允许的参考纠缠均可保留。(7)是一份**新增正常准备输入**，不是已经实现的概率制备仪器或已证明局部可用的源；$p_{empty}$只证明该数学输入非零。

令$\Xi=\widehat U\Omega_0$。由(5)—(7)，

$$
c(G_i)\Xi=0,\qquad
\boxed{\widehat UJ_i=-i c^\dagger(G_i)\Xi.}
\tag{8}
$$

在B的两目标模式上，$M_{ij}=c^\dagger(G_i)c(G_j)$是单占据码的完整偶矩阵单位，严格有

$$
J^\dagger\widehat U^\dagger M_{ij}\widehat UJ=|i\rangle\langle j|,
\qquad
J^\dagger\widehat U^\dagger(n_{G_1}+n_{G_2})\widehat UJ=I_2.
\tag{9}
$$

因此B的单占据双轨解码对所有$\rho_{QR}$给$\rho_{BR}=\rho_{QR}$，不是只验两个基态或有限参考样本。有限CAR因子化把$\Xi$写成目标空模与某个固定余部向量；(8)使整个编码输出等于目标码乘以这个与$Q$无关的余部（保共同费米相位）。余部包含原真空关联与实际对产生等变化，$\Xi$一般不等于$\Omega_0$，也不等于$\Omega$。互补输出的$Q$独立性不等于余部恒等。

## 5. 完整真空输出的有限末端均能

实施性本身不足以推出有限能量。以下给出额外加权证书，所有算符、迹与HS范数均在**完整旋量×A/B×两rail**空间计算。

### 5.1 空间快降符号与有限阶消除

在$\mathbb R^3$使用$\Psi^q_{dec}$：其矩阵符号$a(x,p)$对任意多重指标$\alpha,\beta$及任意$L\ge0$满足

$$
\|\partial_x^\alpha\partial_p^\beta a(x,p)\|
\le C_{\alpha\beta L}\langle x\rangle^{-L}\langle p\rangle^{q-|\beta|},
\tag{10}
$$

并包括具有同样空间快降的平滑余项。时间族在紧窗上各符号半范数一致。$V(t)\in\Psi^0_{dec}$；$\omega,\omega^{-1},\epsilon$是平移不变的阶$1,-1,0$符号，$m>0$避免低频奇点。

标准全局伪微分复合给：空间快降类是与这些乘子及有界零阶符号复合的理想；乘法相加阶数且保留(10)；$[\omega,a]\in\Psi^q_{dec}$、$\omega^{\pm1}a\omega^{\mp1}\in\Psi^q_{dec}$，$\omega^{-1}a,a\omega^{-1}\in\Psi^{q-1}_{dec}$。此处$\omega$的主符号是标量，交换子降低一阶，恰抵其阶1；不对一般矩阵主符号无条件使用降阶。

还需明确指数闭合：若$\alpha^*=-\alpha\in\Psi^{-k}_{dec}$、$k>0$，则$e^\alpha-I\in\Psi^{-k}_{dec}$，并保时间光滑性。可用标准带阶数Beals刻画／伪微分指数演算：$e^{t\alpha}$在每个带多项式空间权重的Sobolev空间有界；$e^\alpha-I=\int_0^1e^{t\alpha}\alpha\,dt$；对$x_j,D_j$的各阶交换子使用Duhamel，把交换落在$\alpha$上，保住全部空间权重和所需动量阶。这里没有以算符范数Taylor收敛冒充全部符号半范数收敛。

在以上理想内使用[Langmann—Mickelsson §2，Theorem 2.10及Lemma 2.12](https://arxiv.org/html/hep-th/9509034)：逐阶酉消除使$[\epsilon,A_N]\in\Psi^{-N}_{dec}$；变换只依赖当时势及有限时间导数，故势在端点邻域消失时$T_N=I$。这些成熟递推用到的复合、交换及指数已在上面验明空间条件。

为避免原文变换方向与负号生成元记号混用，本稿固定自洽约定：在每步令$a_j=\frac18\{\omega^{-1},[\epsilon,A_j]\}$、$W_j=e^{a_j}$，取$U_{j+1}=W_jU_jW_j^\dagger$，$A_{j+1}=i\dot W_jW_j^\dagger+[W_j,h]W_j^\dagger+W_jA_jW_j^\dagger$。首阶为$A_j-[h,a_j]$，恰消去相应非对角主阶；这是同一消除方法的方向一致写法。总$T_3=W_2W_1W_0$。

具体地，$[\epsilon,A_j-[h,a_j]]=\frac14[\omega^{-1},[\omega,[\epsilon,A_j]]]$，右端降阶；其余项由以上理想与指数闭合控制。取三步即可得到

$$
U_3(t,s)=T_3(t)U(t,s)T_3(s)^\dagger,\quad
i\partial_tU_3=(h+A_3(t))U_3,\quad
B_3(t)=[\epsilon,A_3(t)]\in\Psi^{-3}_{dec},\quad T_3(0)=T_3(\tau)=I.
\tag{11}
$$

### 5.2 加权Duhamel界

$\omega^{1/2}B_3\in\Psi^{-5/2}_{dec}$。符号对$x$快降、对$p$为阶$-5/2<-3/2$，在$(x,p)$上平方可积，故为HS；符号积分或核的Plancherel公式给出该结论。**仅有负阶或仅有所有Schatten性质不够，空间可积及能量权重均已另验。**

$A_3$为光滑有界零阶族，其$H^{1/2}$共轭有界；自由群在该范数酉。Duhamel／Gronwall在有限窗上给

$$
M(t):=\|\omega^{1/2}U_3(\tau,t)\omega^{-1/2}\|<\infty,
\qquad\sup_{0\le t\le\tau}M(t)<\infty.
\tag{12}
$$

因$[\epsilon,h]=0$，交换子Duhamel式及右酉因子给

$$
[\epsilon,U_3(\tau,0)]=-i\int_0^\tau
U_3(\tau,t)B_3(t)U_3(t,0)\,dt,
$$
$$
\boxed{\|\omega^{1/2}[\epsilon,U(\tau,0)]\|_{HS}
\le\int_0^\tau M(t)\|\omega^{1/2}B_3(t)\|_{HS}dt=:L<\infty.}
\tag{13}
$$

端点$T_3=I$使此式回到原实际过程；不把变换后过程偷换成新物理控制。$L$是本固定脉冲的解析有限资源表达，不声称给出最优数值。

正能物理场Hamiltonian $H_F$以自由真空正常序，真空末态协方差为$USU^\dagger$。按正粒子及负谱空穴分别计能，

$$
E_{vac}:=\langle\widehat U\Omega,H_F\widehat U\Omega\rangle
=\|\omega^{1/2}P_+UP_-\|_{HS}^2+
\|\omega^{1/2}P_-UP_+\|_{HS}^2
=\frac14\|\omega^{1/2}[\epsilon,U]\|_{HS}^2\le L^2/4.
\tag{14}
$$

$1/4$来自两个非对角块中的因子2；A/B及rail的重数已经包括在整空间范数中，不能再漏乘或重复乘。此处只签末端平均能量，不追加全部能源矩或窗内裸生成元的域声明。

## 6. 准备与未知码的有限能量

令$F:\mathbb C^2\to\mathcal K$以$F_1,F_2$为列，$E=I-FF^\dagger$，$R=I-F^\dagger SF=(1-s)I_2$。条件空态(7)仍为纯规范不变准自由态，其协方差为

$$
C_0=E\bigl[S+SF R^{-1}F^\dagger S\bigr]E.
\tag{15}
$$

单rail在$\mathbb Cf\oplus f^\perp$下，若$S=\begin{pmatrix}s&v^\dagger\\v&A\end{pmatrix}$，则
$C_0=\begin{pmatrix}0&0\\0&A+vv^\dagger/(1-s)\end{pmatrix}$；正号由空条件的CAR收缩确定。故$C_0-S$有限秩，所有腿在$\operatorname{span}\{F_i,SF_i\}$内。对任意逻辑密度$\rho$，编码态的协方差为

$$
C_\rho=C_0+F\rho F^\dagger,\qquad \Delta_\rho=C_\rho-S.
\tag{16}
$$

对纠缠$QR$取$Q$边缘即可计算场能量，参考不增加该场能预算。$f$光滑紧支，$S$为阶0乘子，$F_i,SF_i\in H^{1/2}$；实际$U$由光滑有界势保持$H^{1/2}$，所以$U\Delta_\rho U^\dagger$的有限秩腿仍在能量形式域。相对(14)的修正严格有限：

$$
\boxed{E_\rho-E_{vac}=\operatorname{Tr}_{form}
\bigl(h\,U\Delta_\rho U^\dagger\bigr).}
\tag{17}
$$

这里形式迹以$\langle\omega^{1/2}v,\epsilon\omega^{1/2}w\rangle$解释每个秩一腿，未先假设无界算符乘积为迹类。更显式地，把$C_0-S=\sum_{a=1}^d d_a|w_a\rangle\langle w_a|$作有限维谱分解，并令$K_\tau=\|\omega^{1/2}U\omega^{-1/2}\|$，则

$$
0\le E_\rho\le L^2/4+K_\tau^2\left[
\sum_{a=1}^d|d_a|\|\omega^{1/2}w_a\|^2+
\|F^\dagger\omega F\|\right]<\infty,
\tag{18}
$$

统一覆盖全部$\rho_{QR}$。同一论证在输入端取$U=I$也给有限准备均能。它是全部材料场的实际正能账，不是1126零能辅助的四P缺陷公式，也不是外源所供能量、守恒总能源或应力反冲的证明。

## 7. 验证范围与停止线

[check.py](check.py)仅作三组有区分力的有限代数诊断：相关有限CAR环境的条件准备、全部码矩阵单位与纠缠参考；保两谱的$h$与时间分量混合及错误$\beta$质量混合对照；准备协方差有限秩修正与加权HS能量恒等式。默认只读比较[results.json](results.json)，仅显式`--write-results`登记本轮结果。

有限矩阵不是连续真空、局域因果传播或无限维能量界的证明，也不把数值样本说成覆盖任意参考。相应量词分别由§2、§3—4、§5—6承担。通过后只领取本例外给脉冲下的一端实际材料码交接。正常准备来源、自主控制与定位、作用来源和相位／应力、另一端及完整D-P仍开放。
