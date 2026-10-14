# 第1189轮解析证明：正TT有限带与同一Dirac来源的受驱共同过程

2026-10-14。作者证明；代码、独立审阅及最终验收另行交付。本轮只完成一个有限有效域的共同实现桥，不产生新的引力作用或认知公理。

## 1. 净命题及继承对象

[1180](../research_note_1180.md)已给正的高阶局域关系成本，在固定非零频带逼近线性TT频率，并证明**自由**振子的完整Fock后态误差。[940](../../archive_935_955/research_note_940.md)已给同一Dirac来源、正物理传播模式与非传播约束势的共同量子Hamiltonian。不能直接拼贴两项结论：源与物质Hamiltonian、不同源之间一般不对易；持续交互会把有限占据初码带到任意高占据，1180的自由初码界不能原样套用。

本轮在两模型的同一物质因子、同一$Q,P$、同一源及原约束势上，只改变TT势的频率系数。证明完整受驱后态与原无界来源插入具有共同的显式误差界。没有将输出裁回初码，也不将量子源替换成经典规定函数。

实际Dirac顶点、共同工作背景、TT物理部门及$V_g$继续按940输入。1180的正成本限制属于其声明的张量/TT类别；本轮没有生成其物理约束，没有证明该有限调节器的严格微观空间局域性，也没有证明实际内部钟资格。把关系自由度独立放入内部整体是允许的模型选择，不要求几何只能由物质复合。

## 2. 同一Hilbert空间、生成元和误差参数

取940声明的有限原物质空间$\mathcal H_F$，及有限个物理玻色模式的**完整**Fock空间：

$$
\mathcal H=\mathcal H_F\otimes\bigotimes_{a=1}^ML^2(\mathbb R),\qquad
[Q_a,P_b]=i\delta_{ab},\qquad \omega_a>0.
\tag{1}
$$

取$\hbar=1$。全部$J_a=J_a^\dagger$只作用于$\mathcal H_F$，故有界；不要求$[J_a,J_b]=0$或$[H_F,J_a]=0$。定义

$$
\begin{aligned}
H_0&=H_F+V+\sum_a\frac{P_a^2+\omega_a^2Q_a^2}{2}
       -\sum_a g_aQ_aJ_a,\\
H_+&=H_F+V+\sum_a\frac{P_a^2+\omega_{+,a}^2Q_a^2}{2}
       -\sum_a g_aQ_aJ_a,\\
\delta_a&=\omega_{+,a}^2-\omega_a^2,\qquad
H_+-H_0=\frac12\sum_a\delta_aQ_a^2.
\end{aligned}
\tag{2}
$$

只在TT模式允许$\delta_a\ne0$；保留电磁模式时，其频率和全部源不变。所有$\omega_{+,a}>0$。$V$包含940同一$V_e,V_g$、原正常序/对称乘积及原能源基线。$H_0$是声明的弱场目标，$H_+$是接入正TT成本的比较过程；两者均非完整非线性GR。

1180的频率误差可以给

$$
\left|\frac{\omega_{+,a}}{\omega_a}-1\right|\le\epsilon_a
\quad\Longrightarrow\quad
|\delta_a|\le(2\epsilon_a+\epsilon_a^2)\omega_a^2.
\tag{3}
$$

若将1180表达式缩放到940实际动量带，缩放、模式和时间单位必须共同登记；不能把1180的64模数值直接当作另一调节器的证书。以下证明只需(2)中实际系数，不要求无限UV或包含零波数的统一极限。

以$\omega_a$定义同一参考产生/湮灭算符：

$$
Q_a=\frac{a_a+a_a^\dagger}{\sqrt{2\omega_a}},\quad
P_a=-i\sqrt{\frac{\omega_a}{2}}(a_a-a_a^\dagger),\quad
N=\sum_a a_a^\dagger a_a,\quad A=N+1.
\tag{4}
$$

初码投影为$\mathsf P=I_F\otimes\mathbf1_{N\le N_{\max}}$，$m=N_{\max}+1$。若逐模式上限为$n_a\le L_a$，只有取$N_{\max}\ge\sum_aL_a$时，该逐模式码才是此总占据码的子空间并可使用相同上界。不能将“每模$\le1$”误写成“总占据$\le1$”。$\mathsf P$不是演化后的测量或后选择。取

$$
\begin{aligned}
j_a&=\|J_a\|,\qquad
b=\sqrt2\sum_a\frac{|g_a|j_a}{\sqrt{\omega_a}},\\
d&=\frac54\sum_a\frac{|\delta_a|}{\omega_a},\qquad
r=\frac32\sum_a\frac{|\delta_a|}{\omega_a},\\
z_0(t)&=\sqrt m+\frac b2t,\qquad
z_+(t)=e^{rt/2}z_0(t),\qquad t\ge0.
\end{aligned}
\tag{5}
$$

$b$包括全部受驱模式，$d,r$只计频率实际改变的模式。三者均为可从原矩阵和实际频率计算的保守常数，不是新增动力参数。

## 3. 自伴性、下界与共同域

有限模式、全部正频率给标准振子自伴算符。有限物质矩阵有界；$Q_aJ_a$相对任一正频振子Hamiltonian有任意小相对界。因此两Hamiltonian在共同振子图域$D(A)$上自伴。Schwartz函数乘有限物质空间是共同核；不同固定正频振子的图范数等价。这一结论可由展开$\|(P_a^2+\omega_a^2Q_a^2)\psi\|^2$及分部积分控制各$P_a^2,Q_a^2$取得，不由形式上每项有下界替代域论证。

完成平方给

$$
H_\nu\ge\lambda_{\min}(H_F+V)
-\sum_a\frac{g_a^2j_a^2}{2\omega_{\nu,a}^2},
\qquad \nu\in\{0,+\}.
\tag{6}
$$

不同$J_a$不对易不影响逐项正平方下界。此处“正”指Hilbert内积及关系模式的正成本；总$H$有下界不自动等于$H\ge0$。若数值对象需要总正谱，应核原参数下界，不能免费更换原能源基线或分别减去真空能。零点能始终保留。

## 4. 驱动后数矩：非对易源无需冻结

对归一$\psi\in D(A)$，产生/湮灭算符给

$$
\|a_a\psi\|\le\|A^{1/2}\psi\|,\qquad
\|a_a^\dagger\psi\|\le\|A^{1/2}\psi\|.
\tag{7}
$$

$H_F,V$与$A$对易。原受驱过程$U_0(t)=e^{-itH_0}$满足

$$
\|[A,H_0]\psi\|\le b\|A^{1/2}\psi\|.
\tag{8}
$$

这里用的是当时联合态上的算符不等式，$J_a$无需是常数或守恒源。令$y_0(t)=\|AU_0(t)\psi\|$。由Schrödinger方程及Cauchy–Schwarz，

$$
\frac{d^+y_0}{dt}\le b\sqrt{y_0},\qquad
\sqrt{y_0(t)}\le\sqrt{y_0(0)}+\frac b2t.
\tag{9}
$$

因此对全部初码向量及任意关联参考一致有

$$
\boxed{\|AU_0(t)\mathsf P\|\le z_0(t)^2.}
\tag{10}
$$

为比较模型的来源也取图范数，注意

$$
[A,Q_a^2]=\frac{(a_a^\dagger)^2-a_a^2}{\omega_a},\qquad
\|a_a^2\psi\|\le\|A\psi\|,\quad
\|(a_a^\dagger)^2\psi\|\le2\|A\psi\|.
\tag{11}
$$

故$\|[A,H_+-H_0]\psi\|\le r\|A\psi\|$。同样令$y_+(t)=\|AU_+(t)\psi\|$，有

$$
\frac{d^+y_+}{dt}\le b\sqrt{y_+}+ry_+,
\qquad
\boxed{\|AU_+(t)\mathsf P\|\le z_+(t)^2.}
\tag{12}
$$

最后使用$\sqrt{y_+}\le e^{rt/2}(\sqrt{y_+(0)}+bt/2)$，无需除以$r$，因此$\delta=0$时无奇异定义。

这些微分不等式可先在有限数截断上严格计算：截断投影与$A$对易，交换子上界不变。将交互及差分二次项压到截断，保自由振子在补空间上；截断算符自伴并在共同核上收敛到(2)。共同核的强预解/酉极限、$A$的闭性和下半连续性给(10)、(12)，随后推广到全部$D(A)$初态。该截断仅为域与极限证明，不是最终过程的占据裁剪。

## 5. 完整Fock后态的统一界

由

$$
Q_a^2=\frac{a_a^2+(a_a^\dagger)^2+2a_a^\dagger a_a+1}{2\omega_a}
$$

及(11)，保守地有

$$
\|Q_a^2\psi\|\le\frac5{2\omega_a}\|A\psi\|,
\qquad
\|(H_+-H_0)\psi\|\le d\|A\psi\|.
\tag{13}
$$

Duhamel公式保留同一完整后态：

$$
(U_+(t)-U_0(t))\mathsf P
=-i\int_0^tU_+(t-s)(H_+-H_0)U_0(s)\mathsf P\,ds.
\tag{14}
$$

利用(10)、(13)，得到本轮主界

$$
\boxed{
\|(U_+(t)-U_0(t))\mathsf P\|
\le\eta(t):=d\left[mt+\frac{b\sqrt m}{2}t^2+\frac{b^2}{12}t^3\right].
}
\tag{15}
$$

亦可取与2的较小者作为向量范数上界。$\eta$中的源项反映真实受驱出带，不能把$b$置零后继续声称已接入物质。有限$\mathcal H_F$并不意味着源可相互对角化；本证明没有使用这种加强。

算符范数在张量任意$I_R$后不变。因此对初码内任意未知态、任意外部关联参考，两完整输出的半迹距离不超过$\min(1,\eta)$。接同一完整CP仪器并保全部分支/记录，距离仍不增。这个结论不授予新仪器的实际权限，不认证未知初态的免费准备，也不保证任意后续无界源在没有能源条件时仍有同样误差。

## 6. 原引力耦合来源：必须加权，不能仅用迹距

只改变TT频率、保持其与$\kappa$无关；原质量、动量、参考、正常序处方也固定。继承940来源公式、在本轮声明的参数点取值，两模型中有同一个算符：

$$
V=V_e+\kappa^2V_g^{(0)},\qquad
S_\kappa=\partial_\kappa H_0=\partial_\kappa H_+
=2\kappa V_g^{(0)}-\frac12\sum_{a\in\mathrm{TT}}Q_aJ_a.
\tag{16}
$$

这里$J_a$取940对应的原TT应力模式，Fourier及实模归一仍沿同一字典。$V_g^{(0)}$必须包含原完整密度、流和应力约束势；不能为比较而只留TT辐射项。式(16)没有重新生成该势。$S_\kappa$是本轮有限Hamiltonian的完整$\kappa$参数导数，包含全部已声明$V_g$项；它不是已认证的完整场论Hilbert应力张量$T_{\mu\nu}$，也不认证相互作用Ward恒等式、引力自能或非线性GR来源。

更一般地，在当前固定参数点记$S=V'-\sum_ag'_aQ_aJ_a$，$J_a$对该参数无显式依赖。定义

$$
u=\|V'\|,\qquad
w=\sqrt2\sum_a\frac{|g'_a|j_a}{\sqrt{\omega_a}}.
\tag{17}
$$

$Q$线性项由(7)给$\|Q_a\psi\|\le\sqrt{2/\omega_a}\|A^{1/2}\psi\|$。所以

$$
\|S U_\nu(t)\mathsf P\|\le u+wz_\nu(t),\qquad \nu\in\{0,+\}.
\tag{18}
$$

将两份Heisenberg来源作同一码上的比较，使用$SU_\nu\mathsf P$及其伴随的有界延拓，得

$$
\boxed{
\left\|\mathsf P\left(U_+^\dagger S U_+-U_0^\dagger S U_0\right)\mathsf P\right\|
\le\eta(t)\left[2u+w\big(z_+(t)+z_0(t)\big)\right].
}
\tag{19}
$$

证明在两项差中分别放置$U_+-U_0$；另一因子分别是$SU_+\mathsf P$或$(SU_0\mathsf P)^\dagger$，不能直接估计无界$S$的全空间范数。张量参考后同界成立，故它控制任意初码输入的来源期望及其有界码上二次型。对区间上同一加权来源积分，可积分右侧再乘实际权函数绝对值。

这里的$\mathsf P$是将输入范围写成算符证书，不代表在来源之前投影后态。高占据激发及原约束势一直保留。若实际来源参数还改变$\omega$、$J$、测度或准备，必须共同微分；那不属于(16)—(19)已证明的对象。

**证明强度止于同一来源的Heisenberg插入/期望界。** 没有据此签收完整$\partial_\kappa U$通道的范数界、所有高阶来源插入或任意新仪器后来源；这些通常需要更高图范数与独立误差责任。

## 7. 固定实例及可复算有理包络

本实例继承940的Dirac顶点、实模提升、二次量子化及对称正常序公式，**不是重放940原诊断的数值参数**。新参数一次固定为

$$
M=4,\quad k=\omega_a=1,\quad m_D=1,\quad
\mathbf p=(1,0,0),\quad\mathbf p'=(1,0,1),\quad
\kappa=\frac1{10},\quad g_a=\frac\kappa2,\quad 0\le t\le1.
\tag{20}
$$

物质空间是四个CAR动量—自旋模式的二粒子部门，维数6；四个玻色实模是两种TT极化各自的cos/sin，均保完整Fock。电磁耦合在本实例固定为零；没有通过扔掉引力约束项换取预算。取$N_{\max}=1$，故初码是总占据不超过1、维数$6(1+4)=30$的未知态空间，$m=2$。参考$R$仍任意。

### 7.1 源和完整约束势的解析范数

940使用归一Dirac旋量$u^\dagger u=I$。令$j^\mu=u(\mathbf p')^\dagger\alpha^\mu u(\mathbf p)$，则各$\|j^\mu\|\le1$。同一顶点满足

$$
T^{\mu\nu}=\frac14\big[(p+p')^\mu j^\nu+(p+p')^\nu j^\mu\big],
\qquad p^0+p'^0=\sqrt2+\sqrt3<\frac{16}5.
\tag{21}
$$

因此本参数点的raw TT矩阵为$j^x/\sqrt2,j^y/\sqrt2$。对raw矩阵$B$定义940提升

$$
\ell=\begin{pmatrix}0&0\\B&0\end{pmatrix},\qquad
L_c(B)=\frac{\ell+\ell^\dagger}{\sqrt2},\quad
L_s(B)=\frac{\ell-\ell^\dagger}{i\sqrt2},\qquad
\|L_{c,s}(B)\|=\frac{\|B\|}{\sqrt2}.
\tag{22}
$$

二粒子空间上$\|d\Gamma(A)\|\le2\|A\|$，所以全部四个$J_a=d\Gamma(L_{c,s}(T^{TT}))$满足$\|J_a\|\le1$。

令正常序双线性为

$$
\mathcal N(A,B)=\frac12\{d\Gamma(A),d\Gamma(B)\}
-d\Gamma\!\left(\frac12\{A,B\}\right).
\tag{23}
$$

在二粒子部门，它是$A\otimes B+B\otimes A$限制于反对称空间，因此$\|\mathcal N(A,B)\|\le2\|A\|\|B\|$。保留940的全部引力约束势：

$$
\begin{aligned}
V_0&=-\frac18C_g,\qquad V_g=\kappa^2V_0,\\
C_g&=\sum_{r_0=c,s}\left[
\frac12\mathcal N(\rho_{r_0},\rho_{r_0})
+\mathcal N(\rho_{r_0},\sigma_{r_0})
-2\mathcal N(v^x_{r_0},v^x_{r_0})
-2\mathcal N(v^y_{r_0},v^y_{r_0})
-\frac32\mathcal N(v^z_{r_0},v^z_{r_0})\right],
\end{aligned}
\tag{24}
$$

其中$\rho=L(T^{00})$，$\sigma=L(\sum_iT^{ii})$，$v^i=L(T^{0i})$。由(21)，对应raw矩阵的范数上界可分别取

$$
R=\frac85,\qquad S=\frac32,\qquad
v_x=\frac{13}{10},\quad v_y=\frac45,\quad v_z=\frac{21}{20}.
\tag{25}
$$

每个提升贡献$1/\sqrt2$，两个实模全部计入，故

$$
\|V_0\|\le\frac14\left[\frac{R^2}{2}+RS+2v_x^2+2v_y^2+\frac32v_z^2\right]
=\frac{1599}{640}<3.
\tag{26}
$$

这不是只计算辐射TT部分。$H_{\rm mat}=H_m+\kappa^2V_0$及全部(24)在$H_0,H_+$中完全相同。原940来源家族在此参数点定义了本轮实际矩阵；并未要求新模式、质量、耦合等于旧诊断值。数值可检查非对易性，但(26)不依赖浮点谱估计。

### 7.2 固定1180成本与完整Fock误差

1180已选定的32阶正成本及$10^{-12}$正底不再调高精度。复用其全带解析相对频率误差

$$
\epsilon=\frac{364325745847155802203}{2384185791015625000000000000},\qquad
\Delta=2\epsilon+\epsilon^2.
\tag{27}
$$

四个模式均为$k=1$，属于1180已认证频带。用全带包络而非对该谱点单独拟合，有$|\delta_a|\le\Delta$。记大写$D$为$d$的统一上界，则

$$
b\le\frac25,\quad D=5\Delta,\quad r\le6\Delta,
\qquad u=\|2\kappa V_0\|\le\frac35,\quad w\le4.
\tag{28}
$$

这里更直接的$b\le\sqrt2/5,w\le2\sqrt2$已放宽为有理界。取$\sqrt m=\sqrt2\le3/2$，$z_0(1)\le17/10$，得到

$$
\eta(t)\le\eta(1)\le
D\left(2+\frac3{10}+\frac1{75}\right)
=D\frac{347}{150}
\approx3.534988588683327\times10^{-6}<4\times10^{-6}.
\tag{29}
$$

(27)—(29)由有理数直接比较完成；小数只是显示。共同过程在原同一能源基线上也有正谱下界：由(6)、$H_m\ge2\sqrt2$，

$$
H_\nu\ge2\sqrt2-\frac3{100}-\frac{1}{200(1-\epsilon)^2}>2,
\tag{30}
$$

其中$\omega_{\nu,a}\ge1-\epsilon$，没有分别减去各自真空能以制造谱匹配。共同表达式(2)中的振子零点能保持不动。

完整Fock输出不要求无限初始资源。对固定30维初码内任意归一输入及任意参考，$\|H_m\|\le2\sqrt3<7/2$，原自由振子能源至多3，$\kappa^2\|V_0\|<3/100$；每模$\|Q_a\mathsf P\|\le\sqrt{2m}=2$，线性耦合的初态期望绝对值总和不超过$2/5$。频差项的初态期望绝对值不超过$dm=2d<4\times10^{-6}$。因此有共同的保守有限资源界

$$
0<\langle H_\nu\rangle_{t=0}<8,\qquad
\langle H_\nu\rangle_t=\langle H_\nu\rangle_{t=0},\quad \nu\in\{0,+\}.
\tag{30a}
$$

守恒针对各自同一时不变Hamiltonian，未声称两模型初始能源期望逐态精确相等，也未删去实际激发的高占据尾。任意参考指允许已有相关性；此处有限能源界只计声明的物质—几何系统，不给任意外部参考的总能源免费上界。


### 7.3 同一参数导数来源误差

本实例对$\kappa$求导时，固定$k,m_D,\mathbf p,\mathbf p'$及两套频率，故

$$
S_\kappa=2\kappa V_0-\frac12\sum_{a=1}^4Q_aJ_a
\tag{31}
$$

精确相同。定义$\bar r=6\Delta$，利用$e^{\bar r/2}\le1/(1-\bar r/2)$及(19)，在$0\le t\le1$上

$$
\begin{aligned}
\|\mathsf P(U_+^\dagger S_\kappa U_+-U_0^\dagger S_\kappa U_0)\mathsf P\|
&\le\left(D\frac{347}{150}\right)
\left[\frac65+4\frac{17}{10}\left(1+\frac1{1-\bar r/2}\right)\right]\\
&\approx5.231785315184274\times10^{-5}<6\times10^{-5}.
\end{aligned}
\tag{32}
$$

这界仍对全部未知初码及任意关联参考成立，控制本有限过程的完整$\kappa$来源插入，不扩称全场应力来源。

### 7.4 数值诊断的责任

[check.py](check.py)与[results.json](results.json)计算(27)—(32)的精确有理包络，并用声明的有限占据诊断检查非零驱动、$[H_m,J_a]\ne0$、源间非对易、完整$V_0$保留及固定频率的来源导数。它们继承940公式并公开上述新参数。有限矩阵只显示具体交互确实发生；全Fock与任意输入/参考量词由第3—6节解析证明承担，不能由单输入演化、Ward数值残差或截断误差小替代。

本轮无需增加多项式阶数、优化占据、重调耦合或模拟无限体积。自由$b=0$时(15)回到自由数矩界；真实本例$b\ne0$，且允许超出初码的占据被保留。

## 8. 结果定位与停止线

本轮若数值对象通过验收，取得的是：**正TT有限带成本与原物质—约束势能够在同一受驱量子过程内比较，完整未知输入/参考及原耦合来源均有有限域误差。** 它消除的是“自由频率逼近一接实际非对易物质源就没有任何共同证书”的接口空缺。

未减少的输入是：原TT物理约束、两导数目标、Dirac顶点、$V_g$的Einstein系数、工作背景及材料参数。本轮既不从一般认知原则选择它们，也不领取358/1029条件自耦的前提已经由新母模型产生。尤其静态Newton交换来自完整$V_g$，不是正TT频率拟合的输出。

591的正诱导惯量不适配未约束六分量ADM签名，并不阻碍这里已经约化的正模式；但本轮也没有从591外参数族生成几何变量。不得用纯材料复合分支的困难排除允许内部独立关系量的总目标。

关闭这一次受驱共同实现接口后停止常数、占据、谱点与高阶多项式精调。后继若要增加科研轮次，须真正选择原本自由的约束/自耦输入或建立另一项有物理选择力的同源机制；不能继续扩大本比较器规模代替引力生成。
