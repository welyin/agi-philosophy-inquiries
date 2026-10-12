# 1128证明：固定三维Wilson材料的自治局部装载

2026-10-11。接[1127](../1127/proof.md)的一端外给脉冲接口，本轮另声明一份**有限物理格点有效模型**：有非零三维最近邻传播的两材料Dirac CAR、正常负谱真空、有限交换平台，以及每站独立的正常钟。一个固定有限实例的完整未知双轨码与任意参考交接误差小于1%；同一自治Hamiltonian保留钟反作用、端点交互能和局部能源流。

这是新增的格点材料实例，不是1127连续Fock过程的已证离散近似。格点、固定支撑、正常准备与读取权限是输入；共同Lorentz、定位架的空间动量及连续匹配未闭合。新增认知公理0、科学校准0，阶段active。一般钟Duhamel、有限CAR及格点尾界直接复用，不另记成熟工具为新定理。

## 1. 固定三维材料及全部两谱

取自然单位，固定

$$
a=m=r_W=T=1,\qquad\Lambda=\{-22,\ldots,22\}^3,\qquad |\Lambda|=91125.
\tag{1}
$$

$T_j$是沿第$j$轴的单位格移，移出盒子的分量置0，即开放边界。选$\alpha_j=\tau_x\otimes\sigma_j$、$\beta=\tau_z\otimes I_2$，单材料单rail空间为$\mathcal K_0=\ell^2(\Lambda)\otimes\mathbb C^4$，

$$
h=\beta+\sum_{j=1}^3\left[
\alpha_j\frac{T_j-T_j^\dagger}{2i}
+\beta\left(I-\frac{T_j+T_j^\dagger}{2}\right)\right].
\tag{2}
$$

这是最近邻Wilson选择，Wilson项与空间传播一起保留；不把naive格点的额外低能物种静默解释为连续原材料。两可区分材料A/B、两条rail均用同一$h$，完整Cauchy空间$\mathcal K=\mathcal K_0\otimes\mathbb C^2_{AB}\otimes\mathbb C^2_{rail}$有$1458000$个CAR模式。

写$h=4\beta+K$。每条有向邻键的旋量系数为$-(\beta+i\alpha_j)/2$或其伴随，其范数为1，故$\|K\|\le6$、$\|h\|\le10$。三轴邻键都非零，不是独立站点qubit。

另令$L_W=\sum_j[I-(T_j+T_j^\dagger)/2]\ge0$。$\beta h=I+L_W+\sum_j\beta\alpha_j(T_j-T_j^\dagger)/(2i)$最后一项反Hermitian，因此

$$
\operatorname{Re}\langle\psi,\beta h\psi\rangle\ge\|\psi\|^2,
\qquad\|h\psi\|\ge\|\psi\|,\qquad |h|\ge I.
\tag{3}
$$

这在开放边界成立，无需假设各同轴差分与Wilson项对易。取$S=1_{h<0}\otimes I_{AB,rail}$，$\Omega$为填满负谱的正常真空。完整物理场能为

$$
H_D=d\Gamma(h\otimes I_{AB,rail})-\operatorname{Tr}(Sh)I,
\quad 0\le H_D,\quad\|H_D\|\le16|\Lambda|\,10=14580000.
\tag{4}
$$

迹中及以下$h$的重数均按完整空间理解。常数减除是本固定材料的同一能源基准，未随钟态或未知Q变化。全CAR正谱粒子、负谱空穴及余部均保留；没有以一粒子Hilbert空间冒充整个物理态空间。

## 2. 正常紧支源码及实际接收区域

令$\Gamma=\tau_y\otimes I_2$，它与(2)的$h$反对易。取归一$\xi=(1,0,i,0)^T/\sqrt2$，$\Gamma\xi=\xi$，并令$f=\delta_0\xi$。由(3)，$\Gamma\operatorname{sgn}(h)\Gamma=-\operatorname{sgn}(h)$，故

$$
s=\langle f,1_{h<0}f\rangle=1/2.
\tag{5}
$$

$F_i=f\otimes e_A\otimes e_i$，$n_i=c^\dagger(F_i)c(F_i)$；完整CAR满足$\{c(v),c^\dagger(w)\}=\langle v,w\rangle$。采用1127同一正常准备结构：

$$
e_0=(I-n_1)(I-n_2),\quad p_0=\|e_0\Omega\|^2=(1-s)^2=1/4,
\quad\Omega_0=2e_0\Omega,\quad J_i=c^\dagger(F_i)\Omega_0.
\tag{6}
$$

两rail真空独立，$c(F_i)\Omega_0=0$，$J^\dagger J=I_2$。任意$\rho_{QR}$经该固定奇宇称双轨码输入，逻辑代数为偶元。原局部真空不是空模；条件投影保留相应全场关联。式(6)是新增正常准备输入，非已造出的制备仪器。有限模型中其场能由(4)一致控制。

交换及读取平台固定为

$$
C=\{x\in\Lambda:|x|_1\le20\},\quad
\eta_x=1_C(x),\quad |C|=(4\cdot20^3+6\cdot20^2+8\cdot20+3)/3=11521.
\tag{7}
$$

最终B读模将严格支撑在$C$。格点上的阶跃$\eta$是逐站硬件部署，不声称它是连续空间光滑场。

## 3. 跳跃尾与局部脉冲的模式误差

取内窗光滑函数

$$
b(s)=e^{-1/[s(1-s)]}1_{(0,1)}(s),\quad I_b=\int_0^1b(s)ds,
\quad\chi(t)=2b(2t-1/2)/I_b,\quad\lambda=\pi/2.
\tag{8}
$$

$\chi\in C_c^\infty((0,1))$，支撑$[1/4,3/4]$、非负且积分为1。$I_b\ge e^{-16/3}/2$，$\|b\|_\infty=e^{-4}$；令$y=1/[s(1-s)]\ge4$，有$|b'|\le y^2e^{-y}\le16e^{-4}$，于是

$$
\|\chi\|_\infty\le16,\qquad \|\chi'\|_\infty\le512.
\tag{9}
$$

这里采用$e^{4/3}<4$的保守界；面积归一由解析定义承担，不靠数值求积。

在现场项$4\beta$的相互作用绘景，$K(t)$仍是最近邻且$\|K(t)\|\le6$。从中心到$C$外至少21跳，故Dyson级数前20项在该外部为0；对全部$0\le t\le1$，

$$
\|(I-\Pi_C)e^{-iht}f\|\le
\sum_{n\ge21}\frac{(6t)^n}{n!}\le
p_*:=\frac{6^{21}}{21!}\frac1{1-6/22}
<0.000590384633.
\tag{10}
$$

这是固定盒上解析尾界；有限格点没有严格连续光锥支撑。质量／Wilson现场项及所有重复往返跳跃都已包含。

同一源码确实传播：令$P_N$投到中心的六个近邻，由各轴两方向块平方和为$I$，$\|P_Nhf\|=\sqrt3$。在$t_0=1/100$，谱Taylor余项给
$\|P_Ne^{-iht_0}f\|\ge\sqrt3/100-1/200>0.01232$。这时$\chi(t)=0$，外给比较过程为自由传播；取近邻归一模式$P_Ne^{-iht_0}f$，其A材料rail1占据数在两个源码基输入间的概率差大于$0.01232^2$。§5的短时Duhamel积分不超过式(21)右侧明确上界的$1/100$，实际自治过程的同一概率差仍大于$1.36\times10^{-4}$；这里没有从实际终时误差反推短时误差。因此这里不仅写了非零邻键，而且同一正常源码在实际过程中有可区分的空间传递。

此下界只验本实例确有非平凡空间传播，不签严格光锥因果性、连续速度1或空间四动量；相关范围仍由固定格点合同限定。

先比较外给格点脉冲$h+\lambda\chi(t)\eta\otimes\sigma_x$的单粒子演化$U_\eta$。仅作为数学比较，令$\eta$换成全盒1，得到$U_1(t)=e^{-iht}\otimes e^{-i\theta(t)\sigma_x}$，$\theta(t)=\lambda\int_0^t\chi$。$U_1(1)F_i=-i(u\otimes e_B\otimes e_i)$，$u=e^{-ih}f$。Duhamel及(10)给

$$
\|U_\eta(1)F_i-U_1(1)F_i\|\le(\pi/2)p_*.
\tag{11}
$$

定义实际接收模式

$$
g=\Pi_Cu/\|\Pi_Cu\|,\quad G_i=g\otimes e_B\otimes e_i.
\tag{12}
$$

分母正；若外尾范数为$r\le p_*$，则$\|u-g\|^2=2(1-\sqrt{1-r^2})\le2r^2$。因此

$$
\|U_\eta(1)F_i-(-iG_i)\|\le\delta:=(\pi/2+\sqrt2)p_*.
\tag{13}
$$

全盒旋转没有被授予实际控制器。实际控制和B读模只在有限平台$C$内。

## 4. 完整CAR校正：不随Fock维数放大

需要把(13)升为全量子码的界，不能仅比较两条一粒子振幅。以下有限CAR引理给维数无关常数。

对任意单位向量$v,w$，记$a=\langle v,w\rangle$、$b=\sqrt{1-|a|^2}$。$b>0$时以$v,e=(w-av)/b$为基定义

$$
R\big|_{\operatorname{span}\{v,e\}}=
\begin{pmatrix}a&-b\\b&\overline a\end{pmatrix},
\qquad R=I\text{ 于正交补}.
\tag{14}
$$

这是SU(2)，$Rv=w$；特征值为$e^{\pm i\vartheta}$，$\cos\vartheta=\operatorname{Re}a$，故$\|R-I\|=\sqrt{2-2\operatorname{Re}a}=\|v-w\|$。其**完整两模式外代数**上的$\Gamma(R)$具有特征值$1,e^{i\vartheta},e^{-i\vartheta},1$，因而

$$
\boxed{\|\Gamma(R)-I\|=\|v-w\|.}
\tag{15}
$$

全部余模只张量恒等，范数不增加。若$b=0$，在任一额外正交模式上配以$\overline a$，仍有(14)—(15)；$a=1$可直接取恒等。不能把复overlap默认为实，也不丢配对相位。

本例两个rail正交，分别令$v_i=U_\eta(1)F_i$、$w_i=-iG_i$，取得偶CAR校正$\mathcal R$，三角不等式给$\|\mathcal R-I\|\le2\delta$。记完整外给脉冲Fock演化为$\widehat U_\eta$；其所有负谱块和余部按同一有限CAR Hamiltonian保留。令$\Xi=\mathcal R\widehat U_\eta\Omega_0$，则

$$
c(G_i)\Xi=0,\qquad
\mathcal R\widehat U_\eta J_i=-i c^\dagger(G_i)\Xi.
\tag{16}
$$

这给B码的全部偶矩阵单位及任意参考精确恒等。其余输出是固定的$Q$独立态，却可以与原环境不同。$\mathcal R$只作证明中的比较，不是实际设备或额外被领取的区域门。

实际解码按B的单占据两模式读逻辑码，其余占据输出送到正交失败标记，是确定性CPTP映射。式(15)—(16)控制完整输出等距的范数差，故对所有$QR$，

$$
\varepsilon_{\rm platform}:=\tfrac12\|\mathcal N_\eta-\mathrm{id}\|_\diamond
\le2\delta=(\pi+2\sqrt2)p_*<0.003525440.
\tag{17}
$$

有限失败支路、场余部和相位没有被后选择删去。

## 5. 每站正常自治钟及完整比较界

每个$x\in C$配置独立$L^2(\mathbb R)$钟，自由参数及初态固定为

$$
M=10^{24},\quad\sigma=10^{-11},\qquad
\varphi(X)=(2\pi\sigma^2)^{-1/4}
\exp[-X^2/(4\sigma^2)+iMX].
\tag{18}
$$

$\langle X\rangle=0$、$\langle P\rangle=M$、$\operatorname{Var}X=\sigma^2$、$\operatorname{Var}P=1/(4\sigma^2)$，自由钟均轨道为$X=t$。钟坐标是每站的内部指针坐标，$P$是其共轭动量，**不是三维空间动量**；$M$是该钟的惯性参数，不被替代成定位架已认证的物理质量。各钟在未知输入之前独立准备／部署，共同定时字典是明确输入，不用一个中心钟算符远控全部站点。

定义本地偶交换流

$$
J_x=\sum_{s=1}^4\sum_{i=1}^2
(c^\dagger_{A,i,s,x}c_{B,i,s,x}+c^\dagger_{B,i,s,x}c_{A,i,s,x}),
\qquad\|J_x\|=8.
\tag{19}
$$

八个独立hopping对的完整CAR谱相加给该范数，不能用一粒子范数1取代。实际自治模型为

$$
\boxed{H=H_D+H_C+V,\quad H_C=\sum_{x\in C}P_x^2/(2M),
\quad V=\lambda\sum_{x\in C}\chi(X_x)J_x.}
\tag{20}
$$

$V$有界，$H_D$有界且正，故$H$在$D(H_C)$自伴并半有界。总输入是$J\rho_{QR}J^\dagger\otimes|\varphi\rangle\langle\varphi|^{\otimes|C|}$，对未知Q不重选H或钟。

参考过程为$\widehat U_\eta(t)\otimes e^{-itH_C}$。自由钟的均方宽度为$\sigma_t=\sqrt{\sigma^2+t^2/(4M^2\sigma^2)}$。标量Lipschitz界通过谱演算给$\|[\chi(X)-\chi(t)]\varphi_t\|\le512\sigma_t$。Duhamel作用于该完整参考态，利用真实传播酉性，得到对全部场输入及任意参考一致的

$$
\varepsilon_C\le\lambda\sum_{x\in C}\|J_x\|\,512\int_0^1\sigma_tdt
\le\frac\pi2\,8\cdot11521\cdot512\,
10^{-11}\frac{80001}{80000}
<0.000741567.
\tag{21}
$$

用了$\sigma_t/\sigma\le\sqrt{1+1/40000}\le80001/80000$。实际钟后态及其反作用留在左端完整自治传播中；自由钟只是误差比较，不是实际冻结源。式(21)复用有界作用的钟工具，不将该一般引理单算本轮增量。

合并(17)、(21)，同一平台解码满足事前固定1%目标：

$$
\boxed{\varepsilon_{\rm handoff}
\le\varepsilon_{\rm platform}+\varepsilon_C
<0.004267007<0.01.}
\tag{22}
$$

所有显示的向上小数预算采用$\pi<22/7$、$\sqrt2<99/70$及(21)的有理平方根上界；故比直接浮点代入略宽。组合数、尾级数及这些有理预算在代码中以精确分数复算。

精确正常准备和声明解码是本实例输入；若另有准备／读出误差，分别按完整通道合同另加。本式不包含任何$a=1$到1127连续场的匹配误差。

## 6. 同一能源、真实反作用和端点交互

每钟初始动能及其总和明确为

$$
E_{C,x}=M/2+1/(8M\sigma^2)=5\cdot10^{23}+1/800,
\qquad E_C=11521 E_{C,x}\simeq5.7605\cdot10^{27}.
\tag{23}
$$

该资源巨大但有限，不作经济性或最优性声明。由(9)、(19)，$\|V\|\le(\pi/2)\,16\cdot8\cdot11521<2317367$；$H\ge-\|V\|$，输入具有有限能源及$D(H_C)$正则性。自治传播保同一完整H，严格有

$$
\Delta\langle H_D\rangle+\Delta\langle H_C\rangle+
\Delta\langle V\rangle=0.
\tag{24}
$$

能源交换不是由接收误差推得。尤其控制器实际受力

$$
\dot P_x(t)=-\lambda\chi'(X_x(t))J_x(t),\qquad
\|\dot P_x(t)\|\le B:=\lambda\,512\cdot8<6436.572.
\tag{25}
$$

场与钟的瞬时算符都由同一H运输；$J_x$一般不守恒。Heisenberg积分式在共同Schwartz芯成立，因力有界延拓到所用正常输入，给

$$
X_x(t)=X_x(0)+tP_x(0)/M
-\frac\lambda M\int_0^t(t-s)\chi'(X_x(s))J_x(s)ds,
$$
$$
\|(X_x(t)-t)\Psi_{in}\|\le\sigma_t+Bt^2/(2M).
\tag{26}
$$

这包括未知Q对钟的真实反作用，界对全部$QR$一致。在$t=0,1$，$\chi(t)=0$。通过**同一标量点$t$**的谱演算，$\|\chi(X_x(t))\Psi_{in}\|\le512\|(X_x(t)-t)\Psi_{in}\|$；无需对两个不对易位置算符滥用一般算符Lipschitz。于是实际端点有

$$
\boxed{|\langle V(0)\rangle|,|\langle V(1)\rangle|
\le\lambda\,8\cdot11521\cdot512
\left[10^{-11}\frac{80001}{80000}+\frac B{2M}\right]
<0.000741567.}
\tag{27}
$$

Gaussian尾部使实际交互没有严格时间截断；(27)给有限端点预算，未将尾部删零。它来自钟力与位移估计，不是用普通通道误差控制无界钟能源。式(4)、(23)、(24)、(27)共同给本固定实例的有限输入／输出能源及交换账。

## 7. 同H的局部能源平衡

将$H_D$拆为现场偶项$d_x$与无向邻键项$k_{xy}$，固定正常序常数按站分配。定义

$$
E_x=d_x+1_C(x)\left[P_x^2/(2M)+\lambda\chi(X_x)J_x\right]
+\tfrac12\sum_{y:\,y\sim x}k_{xy},\qquad H=\sum_{x\in\Lambda}E_x.
\tag{28}
$$

不同站点的偶CAR项及独立钟满足局部对易关系。按此半键分配，$[E_y,E_x]$只可能在格距至多2时非零；并非只准最近邻唯一电流。取

$$
j_{y\to x}=i[E_y,E_x]=-j_{x\to y},\qquad
\dot E_x=\sum_yj_{y\to x}.
\tag{29}
$$

对互异$x,y$，钟动能与对方各项对易，故这些站间流是有界局部项；站内钟力作功及交互变化保在同一$E_x$中。区域求和后内部流严格抵消，余下跨边界流与区域能源变化一致。这个账不是另拟合的压力模型，也没有用形式Ward身份取代具体装置。

固定格点和支撑属于本模型输入。(25)是内部钟反作用，不是完整机械反冲；本轮没有连续空间平移生成元、动态定位架或共同Lorentz证书，不能把(24)、(29)升格为完整连续四动量／应力闭合。

## 8. 复算、净连接与停止线

[check.py](check.py)默认只读复算并比较[results.json](results.json)。代码核三类：固定实例的精确组合数与有理预算；缩小但真正三维的开放Wilson正负谱／传播；完整有限CAR旋转与有限矩阵钟的实际反作用和共同能源流。后两类仅校准代数，有限矩阵钟不具精确CCR，不替代(18)—(27)的解析Gaussian证明。显式`--write-results`才登记本轮结果。

本轮净增量是上述三维有传播材料、正常源码、实际有限区域、全部Q/R、分站自治控制与同H能源的共同有限实例；不把422／524的钟界、525—527的反馈、944—947的自治能源工具重报一遍。连续场匹配、固定支撑与正常准备的来源、实际提取器及完整D-P仍开放。详见[selection_audit.md](selection_audit.md)、[NEXT.md](NEXT.md)。
