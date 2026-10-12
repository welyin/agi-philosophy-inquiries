# 1155：同一有限记录源的分站四交换与共同运动描述

2026-10-12。依[事前准入](admission.md)和[1154后方向资格](../_shared/1120_1156/after1154_covariant_parent_screen.md)。新增认知公理0，科学校准0。证明使用1144—1145的正常线性源工具，但选择新的QND有限材料合同；不把本轮有限CAR源改名为1153的全海局部复合材料，也不继承1154的固定实例。

## 1. 命题、适用范围与净新增

在以下明列模型中，同一正常KG过程同时具有：

1. 对两个独立未知材料码及任意被动参考统一的有限角门和完整场输出界；准备独立指针后成为强记录。
2. 从每份真实源的变分及退迟解取得的分站四交换，满足同一场的弱能动量收支。
3. 固定有限期限内的空间冲量、一般非零的时间交换及同一标量时空核在第二惯性描述中的有限预算。
4. 可声明的非简并有限CAR钟，以及一次有因果缓冲的非对易继续；场不重置。

1145已经证明总KG末四动量的共同运输。本轮新增的是**各站分别承担的区域交换**：总末场能和总空间动量可以很小或相消，而每站的空间交换及运动描述中的能源交换仍为有限量。该分摊由源的实际变分固定，不能任意从总四P倒定义。

本轮不构造正常自治支撑、开关或控制器，不声称全装置四P闭合，不认证冲量的全部测量分布，也不提出新的空间或Lorentz公理。

## 2. 同一实际模型与初态

取已给定的 $3+1$ Minkowski背景、度规 $(+---)$、$\hbar=c=1$ 和质量1正常KG Fock场。空间Fourier采用酉归一，时间Fourier采用直接积分；令

$$
\omega=\sqrt{-\Delta+1},\qquad H_\phi=d\Gamma(\omega),\qquad
P_\phi^\mu=d\Gamma(k^\mu),\quad k^\mu=(\omega_{\mathbf k},\mathbf k).
$$

材料是完整四模CAR，模式为 $a_0,a_1,b_0,b_1$。两站各自的一个粒子双轨空间构成声明未知码

$$
\mathcal Q_A=\operatorname{span}\{a_0^\dagger|0\rangle,a_1^\dagger|0\rangle\},
\qquad
\mathcal Q_B=\operatorname{span}\{b_0^\dagger|0\rangle,b_1^\dagger|0\rangle\}.
$$

保留整个16维材料空间；不删真空或双占据态。取

$$
N_A=a_1^\dagger a_1,\qquad N_B=b_1^\dagger b_1,\qquad
N_i^2=N_i,\quad [N_A,N_B]=0.
$$

这两个偶算符属于不同材料，且其谱可独立变化。初态为任意未知材料码/被动参考R联合态，场为同一正常真空 $\Omega$。记 $W_0\psi=\psi\otimes\Omega$；以下等距界实际延伸到完整16维输入，物理准备权限只需上述双轨码。

沿z部署 $z_A=-3/2,z_B=3/2$，每站实际源核为1144原非负径向单位积分bump：

$$
q_i(\mathbf x)=Z^{-1}\exp\!\left[-\frac1{1-|\mathbf x-\mathbf y_i|^2}\right]
\mathbf1_{|\mathbf x-\mathbf y_i|<1}.
$$

固定支持与调用同一个 $N_i$ 的整站权限为本轮明确输入，不认证管内任意点态材料微因果。源相互作用绘景Hamiltonian为

$$
H_I(t)=\sum_{i=A,B}N_i s_T(t)\phi_t(q_i),\qquad
s_T(t)=\alpha T^{-1/2}f(t/T),
$$

$$
f(u)=30u^2(1-u)^2\mathbf1_{0<u<1},\quad
\|f\|_1=1,\quad F_2=\|f\|_2^2=\frac{10}{7},\quad
N_2:=\|f''\|_1=\frac{40\sqrt3}{3}<24.
$$

该开关是紧支C1且二阶弱导数可积，不误称C∞。空间核C∞。记

$$
C(d)=\langle q_A,\omega^{-2}q_B\rangle,\qquad
\alpha^2=\frac{\pi}{F_2C(3)},\qquad T=10^{10}.
$$

$\alpha$由已知部署一次确定，与未知材料/R无关。以后对中心位置求导时，固定 $\alpha,T$，不微分重校准实验。

1144冻结界直接给

$$
\|q_i\|_2^2<3,\qquad C(3)\ge\frac{e^{-5}}{20\pi}>10^{-4},
\qquad \alpha^2<24000.
\tag{1}
$$

核范数界来自半径1/2球上的归一积分下界及最大值；$C$的下界来自全部点间距离不超过5的正Yukawa核。没有为本轮重拟合源。

## 3. 正常完整过程及强记录界

各材料荷相互对易，场交换子为数。因此Magnus展开三阶及以上精确消失；这是线性源的完整结果，不是有限耦合下截断微扰级数。令

$$
R_T(\omega)=\int dt\,s_T(t)\int_{u<t}du\,
\frac{\sin[\omega(t-u)]}{\omega}s_T(u),
$$

$$
b_i=\frac12\langle q_i,R_Tq_i\rangle,\qquad
\Theta_T=\langle q_B,R_Tq_A\rangle.
$$

实际散射酉为

$$
S_T=\sum_{a,b\in\{0,1\}}\Pi_{ab}\,
e^{i(b_Aa+b_Bb+\Theta_Tab)}W(ah_A+bh_B),
\tag{2}
$$

其中 $h_i(\mathbf k)=-i(2\omega_{\mathbf k})^{-1/2}
\widehat q_i(\mathbf k)\int s_T(t)e^{i\omega_{\mathbf k}t}dt$，$\Pi_{ab}$是整个材料空间的共同谱投影。相位中的两条交叉项抵消二次型的1/2，所以 $\Theta_T$没有额外半因子。

源自相位 $b_iN_i$始终保留。它对未知味叠加不是全局相位，不能从记录比较中删掉。

### 3.1 精确受迫解和有限余项

对任一空间向量q，记

$$
v_q(t)=\int_0^t\frac{\sin[\omega(t-u)]}{\omega}s_T(u)q\,du.
$$

利用 $s_T(0)=s_T'(0)=0$ 两次分部积分，得到精确式

$$
v_q(t)=\omega^{-2}s_T(t)q
-\omega^{-3}\int_0^t\sin[\omega(t-u)]s_T''(u)q\,du.
\tag{3}
$$

又有 $\|s_T\|_1=\alpha\sqrt T\|f\|_1$、
$\|s_T''\|_1=\alpha T^{-3/2}N_2$，且 $\omega\ge1$，故

$$
\|R_T-\alpha^2F_2\omega^{-2}\|
\le\frac{\alpha^2\|f\|_1N_2}{T}.
\tag{4}
$$

于是

$$
|\Theta_T-\pi|
\le\frac{\alpha^2\|q_A\|_2\|q_B\|_2\|f\|_1N_2}{T}
<\frac{1728000}{T}
=:\delta=0.0001728.
\tag{5}
$$

时间Fourier同样两次分部积分给

$$
\|h_i\|\le\frac{\alpha N_2\|q_i\|_2}{\sqrt2T^{3/2}},
\qquad
\|\omega^{1/2}h_i\|\le\frac{\alpha N_2\|q_i\|_2}{\sqrt2T^{3/2}}.
\tag{6}
$$

后一个是独立的能源加权界，不由前一个态距推出。由(1)及 $2\sqrt{24000}\sqrt3\,24/\sqrt2<10000$，

$$
\|h_A\|+\|h_B\|<e:=\frac{10000}{T^{3/2}}=10^{-11}.
\tag{7}
$$

### 3.2 全输入、全部场后态与实际指针

定义保全部实际自相位的参考门

$$
V_T=e^{i(b_AN_A+b_BN_B)}e^{i\pi N_AN_B}.
$$

对每个正交材料谱块，
$\|W(h)\Omega-\Omega\|\le\|h\|$。结合
$|e^{i\Theta_T}-e^{i\pi}|\le|\Theta_T-\pi|$，得到

$$
\|S_TW_0-W_0V_T\|<\epsilon:=\delta+e
=0.00017280001.
\tag{8}
$$

张量任意R不变；未知混态可纯化。全部场、材料及参考关联保留，不以末读概率替代完整输出。声明的同一后续仪器对完整cq后态的半迹距误差至多 $\epsilon$。

具体地，B独立准备 $|+\rangle_B=(|0\rangle_B+|1\rangle_B)/\sqrt2$ 后，A的两条数分支分别生成

$$
|+_{b_B}\rangle_B,\qquad|-_{b_B}\rangle_B,
$$

二者正交；A的实际相对相位 $b_A$及任意R仍保留。这是相干记录；不是复制两个任意未知态。已知自相位确定末读基，不要求实际撤销自相位。在Schrödinger绘景中，还须同时保留第8节的实际自由钟相位。

## 4. 分站交换由真实源取得

设 $J_i(x)=s_T(t)q_i(\mathbf x)$，$j(x)=\sum_iN_iJ_i(x)$。
采用相互作用拉氏项 $-j\phi$，实际Heisenberg场满足

$$
(\Box+1)\phi_H=-j,\qquad
\phi_H(x)=\phi_{\rm in}(x)-\sum_jN_j(G_{\rm ret}J_j)(x).
\tag{9}
$$

定义各站的**源侧交换，收到为正**：

$$
I_i^\mu=\int d^4x\,N_iJ_i(x)\partial^\mu\phi_H(x).
\tag{10}
$$

它的来源有两个相符的身份：

- 对第i份实际源作平移变分 $J_i(x)\mapsto J_i(x-a)$，相互作用作用量的一阶响应为 $-\delta a_\mu I_i^\mu$。
- 在静止参考中，先对相互作用绘景Hamiltonian的源部署参数求导，再运输到Heisenberg绘景。以 $S(t)$ 表示第一段截至t的实际源酉，

$$
F_{i,H}^z(t)
=S(t)^\dagger[-\partial_{y_i^z}H_I(t)]S(t)
=N_i s_T(t)\int d^3x\,(\partial_zq_i)\phi_H
=N_i\int d^3x\,J_i\partial^z\phi_H.
$$

其时间积分正是 $I_i^z$。这里源位置导数只作用 $H_I(t)$ 的既定核，不再对已随来源改变的Heisenberg场作全导数；不能混用两个绘景的势与力。

紧支开关和核允许分部积分：

$$
I_i^\mu=-\int d^4x\,N_i(\partial^\mu J_i)\phi_H.
\tag{11}
$$

所以时间分量的相反数是开关对场作的净功；不能把时间交换事先设为零。

### 4.1 正常定义及同一场的弱收支

对每个材料谱块，(9)是自由正常场加一份实退迟解。(10)因此是实线性场算符加一个有限数值二次型。其测试向量因q的空间C∞紧支有任意所需大动量加权矩；各分量在有限粒子Schwartz核心上定义，并延伸到这里实际产生的正常相干态和有限材料叠加。场自由四P及以下弱应力积分也在该域上。

以入真空Wick序定义实际场应力，或等价地将自由应力、线性移位项、源平方项分开。由(9)按分布微分得

$$
\partial_\mu:T_\phi^{\mu\nu}[\phi_H]:
=-j\,\partial^\nu\phi_H.
\tag{12}
$$

先用光滑时空测试及大空间截止积分，再取包住源与实际传播输出的极限；源前后面处耦合为0。自由场入/出四动量及真实(2)给

$$
\sum_i I_i^\nu
=-\big(S_T^\dagger P_\phi^\nu S_T-P_\phi^\nu\big).
\tag{13}
$$

(12)是地区域的弱来源身份，(13)是其全空间归一；分站项在此之前已由各源变分固定。没有把总差额任意一分为二。也没有把(13)称为已构造自治源的守恒总四P。

以下 $\bar I_i^\mu:=W_0^\dagger I_i^\mu W_0$ 是初始场真空压缩后的有限材料期望算符。其范数不等于无界全场 $I_i^\mu$ 的算符范数，不承诺该量的全部分布或一个已制造的测量仪器。

## 5. 有限T空间冲量：不需要新的核导数常数

由(9)，$\phi_H=\phi_{\rm in}-\sum_jN_jv_{q_j}$；入真空的一阶场均值为0，$\partial^z=-\partial_z$，故

$$
\bar I_B^z
=\sum_jN_BN_j\int dt\,s_T(t)
\langle q_B,\partial_zv_{q_j}(t)\rangle.
\tag{14}
$$

$R_T$是实的自伴谱乘子，与空间平移和导数对易。实q的自项
$\langle q_B,\partial_zR_Tq_B\rangle=0$。将A固定、B中心沿正z移动，则
$\partial_dq_B=-\partial_zq_B$，因而

$$
\bar I_B^z=N_AN_B\,\partial_d\Theta_T(d),\qquad
\bar I_A^z=-\bar I_B^z.
\tag{15}
$$

径向部署还给 $\bar I_i^x=\bar I_i^y=0$。这些是此真实有限时间部署的精确压缩身份。

### 5.1 导数余项

由(3)取出 $\omega^{-3}$ 后，空间导数可以直接作用在谱乘子上，而不估计q的导数：

$$
\|\partial_z\omega^{-3}\|
=\sup_{\mathbf k}\frac{|k_z|}{(1+|\mathbf k|^2)^{3/2}}\le1.
$$

正弦谱乘子与其对易且范数至多1。所以

$$
\|\partial_z(R_T-\alpha^2F_2\omega^{-2})\|
\le\frac{\alpha^2\|f\|_1N_2}{T},
$$

$$
|\partial_d\Theta_T-\alpha^2F_2C'(d)|
\le\frac{\alpha^2\|q_A\|_2\|q_B\|_2\|f\|_1N_2}{T}<\delta.
\tag{16}
$$

此法只用已冻结的核L2范数。所有位置导数均固定本次 $\alpha,T$。

### 5.2 静态比较值与有限数值证书

非重叠相同径向核的Yukawa外部解可由径向方程或球面平均得到

$$
C(d)=\frac{M_q^2}{4\pi}\frac{e^{-d}}d,\qquad
M_q=\int d^3x\,q(\mathbf x)\frac{\sinh|\mathbf x|}{|\mathbf x|},
\qquad d>2.
$$

这只是计算(16)中的比较系数；实际传播仍是(9)的退迟过程。因此

$$
\frac{C'(3)}{C(3)}=-\frac43,\qquad
\alpha^2F_2C'(3)=-\frac{4\pi}{3}.
$$

结合(15)—(16)，

$$
\boxed{\left\|\bar I_B^z+\frac{4\pi}{3}N_AN_B\right\|<0.0001728.}
\tag{17}
$$

两站11分支的交换分别约为 $\pm4\pi/3$；其总和为0。这不是高能极限、无限精度外推或仅本征态Hellmann–Feynman值，而是实际裸初态的有限期限积分压缩，作用于全部声明未知码/R。

## 6. 时间分量：严格保留实际开关交换

令

$$
F_T(\omega)=\int dt\,s_T(t)e^{i\omega t},
$$

$$
E_{\rm self}=\frac12\|F_T(\omega)q_A\|_2^2,\qquad
E_{\rm cross}=\frac12\operatorname{Re}
\langle F_T(\omega)q_B,F_T(\omega)q_A\rangle.
$$

同形平移核有相同 $E_{\rm self}$，Cauchy–Schwarz给
$|E_{\rm cross}|\le E_{\rm self}$。由
$\partial_tv_{q_j}(t)=\int_{u<t}\cos[\omega(t-u)]s_T(u)q_jdu$，
将对称余弦核的半三角积分补成全时间平方，直接得到

$$
\boxed{\bar I_B^0=-N_B^2E_{\rm self}-N_AN_BE_{\rm cross},}
$$

$$
\bar I_A^0=-N_A^2E_{\rm self}-N_AN_BE_{\rm cross}.
\tag{18}
$$

特别是10、01单源分支，开关向场输入的正能为 $E_{\rm self}$，源侧收到 $-E_{\rm self}$。由于f和q非零，紧时开关的Fourier整函数不可能在q所见的整个质量壳频率区间为零，因此 $E_{\rm self}>0$。时间分量一般非零，不可从静止中心不动推成零。

两次分部积分给

$$
E_{\rm self}
\le\frac{\alpha^2N_2^2\|q_A\|_2^2}{2T^3}
<\frac{20736000}{T^3},
$$

$$
\boxed{\|\bar I_i^0\|<\beta:=\frac{41472000}{T^3}
=4.1472\times10^{-23}.}
\tag{19}
$$

完整真空输出的KG能源算符压缩恰为

$$
W_0^\dagger S_T^\dagger H_\phi S_TW_0
=(N_A^2+N_B^2)E_{\rm self}+2N_AN_BE_{\rm cross};
$$

它正好等于 $-(\bar I_A^0+\bar I_B^0)$。总空间动量压缩为0。最大末态场能源小于 $2\beta$，但(17)的分站冲量保持有限。该加权计算独立于(8)，没有从小态距推无界能源，也没有把中途场能、相互作用储能或外部开关责任删去。

## 7. 第二惯性描述：先给真实核，再核分站工作

取主动boost

$$
\binom{t'}{z'}=
\begin{pmatrix}\gamma&\kappa\\\kappa&\gamma\end{pmatrix}
\binom{t}{z},\qquad
\gamma=\frac54,\quad\kappa=\frac34,\quad v=\frac35.
$$

本合同将 $J_i$和内部荷作为标量来源；这是已声明物理类型，不是从认知原则排除其它四流或张量来源。实际新描述使用同一时空源的图像

$$
J_i^\Lambda(x')=J_i(\Lambda^{-1}x'),\qquad
\zeta=\gamma(z'-vt'),\quad\tau=\gamma(t'-vz'),
$$

$$
J_i^\Lambda(t',\mathbf x')
=\alpha T^{-1/2}f(\tau/T)\,q_i(x',y',\zeta).
\tag{20}
$$

保有限核内的完整时间纹理，不能将f移出空间积分后重新按实验室体积归一。沿中心 $\zeta=z_i$，

$$
\tau=t'/\gamma-vz_i,\qquad \Delta\tau=\Delta t'/\gamma.
$$

原来已同步的两站时钟在新参考有相应原点偏移；不能将它们同时重置成 $\tau=t'/\gamma$。初末准备面采用原面之Lorentz图像，或加入实际自由传播走到包住全部源的常时面，不在新 $t'=0$ 重准备。

### 7.1 同一过程和分站矢量身份

自由正常场有实际Lorentz酉 $\mathscr U_\Lambda$，满足
$\mathscr U_\Lambda\phi(x)\mathscr U_\Lambda^\dagger=\phi(\Lambda x)$且保真空。这里密度的材料荷两两对易，类空点处自由场交换子为0；正时向boost改变的仅是这些可交换密度的排序。等价地，直接运输(9)的Lorentz不变退迟核，可得真实源酉及Heisenberg场共同身份。

把(20)和已运输的退迟解代入(10)，以不变四体积换元，导数按四矢量变换，得到

$$
I_i^\mu[J^\Lambda]=
\Lambda^\mu{}_\nu\,
(I_{\rm mat}\otimes\mathscr U_\Lambda)
I_i^\nu[J]
(I_{\rm mat}\otimes\mathscr U_\Lambda^\dagger).
\tag{21}
$$

不是以右端倒定义左端。每站固定源自身的变分标签随之运输，不能把总交换运输正确而重新分摊站点。

正常性及所需加权域保持：在本有限boost下 $\omega_{\mathbf k'}\le2\omega_{\mathbf k}$，真空不变，源的正频位移由真实一粒子表示运输。完整等距界(8)及实际自相位保持。管内任意点态材料操作仍不由此获得；本轮仅使用已声明的标量QND密度。

### 7.2 有限区域能源—动量预算

真空压缩(21)给

$$
\bar I_B'^0=\frac54\bar I_B^0+\frac34\bar I_B^z,\qquad
\bar I_B'^z=\frac34\bar I_B^0+\frac54\bar I_B^z.
$$

所以由(17)、(19)，

$$
\boxed{\|\bar I_B'^0+\pi N_AN_B\|
<\frac54\beta+\frac34\delta
<0.000129600001,}
\tag{22}
$$

$$
\boxed{\left\|\bar I_B'^z+\frac{5\pi}{3}N_AN_B\right\|
<\frac34\beta+\frac54\delta
<0.000216000001.}
\tag{23}
$$

11分支的B站源侧时间交换约为 $-\pi$，而两站合计仍只承担小的末态场四动量。A站相反机械交换在运动参考产生相反的领先能源项；这与(13)相容。末态总能源小不能替代分站责任，也不能令支持反作用免费消失。若改为被动boost，应一起改变符号约定，而非改源。

## 8. 有限CAR实际钟：给定SR下的共同权限

可在同一有限材料上加入实际非简并自由Hamiltonian

$$
H_{D,i}=n_{i0}+1.001\,n_{i1},\qquad
\nu=0.001>0.
$$

它与 $N_i$对易，故所有源式在共同自由绘景中原样成立。味相干算符有实际Heisenberg演化

$$
e^{iH_{D,i}\tau}c_{i0}^\dagger c_{i1}e^{-iH_{D,i}\tau}
=e^{-i\nu\tau}c_{i0}^\dagger c_{i1}.
$$

读口可取完整CAR效应 $c_{i,\varphi}^\dagger c_{i,\varphi}$，
其中 $c_{i,\varphi}=(c_{i0}+e^{i\varphi}c_{i1})/\sqrt2$。
在单粒子码上即为相干相位读取；在全部有限材料空间仍满足 $0\le E\le1$。这不是无Hamiltonian的标签钟。

在相互作用结束后，实际味相位还含源自相位 $b_i$；记录指针的Schrödinger绘景末读基必须包含自由相位 $-\nu\Delta\tau_i$及 $b_B$，二者均已知并保留。有限周期相位钟不单独提供无歧义任意长历时计数，源开启的期限T仍是外部部署。

本合同直接输入：各外部世界管中的该材料按SR proper时间演化；boost按(20)运输同一读事件及钟原点。没有从这个有限模型重新导出SR，没有领取原CO位置、原全Dirac钟、正常运动支撑或自主开关触发。场与钟共享一份已声明proper部署，而非各取一份任意时间字典。

## 9. 同一场上的一次非对易继续，不重准备环境

为给一个明确完整历史，可采用第二份相同源窗

$$
J_i^{(2)}(t,\mathbf x)=s_T(t-T-20)q_i(\mathbf x),
$$

其散射酉为 $S_2$；第一段为 $S_1=S_T$。在A站源管内声明一次偶CAR旋转权限

$$
Y_A=-i(a_0^\dagger a_1-a_1^\dagger a_0),\qquad
U_A=e^{-i\pi Y_A/4},\qquad [U_A,N_A]\ne0.
$$

可以在proper时段 $[T+9,T+11]$用一个面积 $\pi/4$ 的平滑控制系数实现该有限材料旋转。若按共同自由绘景定义上述U，Schrödinger控制Hamiltonian须带实际自由钟的已知相位补偿；该控制及其有限味能源变化为明示外部权限，不冒充自治控制源。

两站全部空间支持的直径不超过5。第一段止于T，中插始于T+9；中插止于T+11，第二段始于T+20。故每对相关支持点之间时间间隔至少9而空间间隔至多5：前段整体位于控制的因果过去，后段整体位于其因果未来。这是实因果缓冲，不是只在某坐标中不重叠。整站扩展控制的内部点态实现并未因此证明。

完整过程为

$$
\mathcal S=S_2(U_A\otimes I_\phi)S_1.
$$

第二段与第一段具有相同自相位、交叉相位及残余范数，源时间平移只给一粒子位移自由相位。设相应理想有限材料门为 $V_1,V_2$。插入并减去共同中间比较项，利用实际后续酉范数为1，

$$
\|\mathcal S W_0-W_0V_2U_AV_1\|
\le
\|S_1W_0-W_0V_1\|
+\|S_2W_0-W_0V_2\|
<2\epsilon=0.00034560002.
\tag{24}
$$

这一步从未把实际第二段输入重置成真空：$W_0$只出现在三角比较的参考过程，第一段真实场的全部残余与相关性在 $\mathcal S$ 中保留。未知材料/R、正常全场以及非对易继续共同满足(24)。若只需第一段后末端U，则误差仅为 $\epsilon$。

第一段的(10)—(23)是这同一历史截至T的实际Heisenberg积分，不由之后的U和第二段改写。**本轮只认证第一段交换的数值。** 第二段实际入场并非真空，不能直接套(18)—(19)；若未来增加第二段交换菜单，必须以精确分支和加权交叉项另算。中插时源关闭，材料能源变化可单独列入外部控制功，未获其正常自治空间动量。

每一QND段内部的类空重排由第7节保证，段间与中插的顺序由上述因果缓冲保证；运输同一控制区域、proper时程、材料相位与场后态，得到同一完整历史。没有宣称任意扩展不对易探针的任意切片都可免费排序。

## 10. 因果范围、验证责任与结果边界

所有实际场解与相位均使用正常退迟/Pauli–Jordan核。若两个源支持整体类空，交叉源相位为0；本次 $T=10^{10}$部署允许因果传播。Yukawa核仅是带显式有限T余项的比较系数，不被当成瞬时信号传播。

本轮解析取得：

- 原核、有限T及预选强度下的完整记录误差 $\epsilon<0.00017280001$；
- 同一源的B站空间交换误差 $\delta<0.0001728$；
- 真正非零的时间交换及独立加权上界 $\beta<4.1472\times10^{-23}$；
- 实际boost源的区域能源与动量误差分别小于0.000129600001、0.000216000001；
- 声明有限CAR钟和一次非对易继续的完整历史界(24)。

工具复算应分别检查CAR/全输入代数、时间卷积与源交换符号、有限常数与boost预算；有限模诊断不能替代上述正常连续场证明。本文不以尚未产生的数值结果预先宣布检查通过。

仍未领取：1153局部复合材料的共同运动、全海身份、任意区域仪器、开关和支持的正常自治实现、总装置守恒四P、所有冲量分布或实际经验真实性。这些范围限制对应所选候选，不是对整个认知研究纲领的否定，也不要求先完成无截断UV问题。

相对于1145，本轮的净成果是：**同一有限记录过程不仅运输末态场资源，还给出由实际源固定的每站四交换；非零区域反作用可在两个惯性描述中共同计算，且与实际有限钟及不重置场的继续相容。**
