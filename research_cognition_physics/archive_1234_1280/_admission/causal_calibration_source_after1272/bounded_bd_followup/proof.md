# 同一最小标量BD参考：有界局域读口的有限正响应窗口

2026-10-17。非编号同材料接续；原目标、CT待验证身份及正式计数不变。本稿只证明原场的谱侧非空性及有限比较过程的正概率、正容差存在；不做场谱数值，不认证实际设备或固定非零引力耦合下的误差。

## 1. 原对象与直接复用

沿[原正常脉冲](../../finite_horizon_transfer_after1272/proof.md)保四维de Sitter静态背景、BD参考及实际最小耦合物质作用，单位$c=\hbar=k_B=1$：

$$
ds^2=-F(r)dt^2+F(r)^{-1}dr^2+r^2d\Omega_2^2,\quad
F(r)=1-H^2r^2,\quad
x=H^{-1}\operatorname{atanh}(Hr),
$$

$$
I_\Phi=-\frac12\int\sqrt{-g}\big[(\nabla\Phi)^2+m^2\Phi^2\big]d^4x,
\qquad H>0,\quad m^2=2H^2.
\tag{1}
$$

$m$是固定材料参数；度规变分时不把它更新为$R/6$。本背景上方程和BD二点函数恰与无质量共形场相同，**不意味着采用共形Hilbert应力或物质Wald接触项**。

原脉冲在$t\ge0$为$\varphi=f(x-t)/r$，$f\in C_c^\infty([a,b])$、$0<a<b$，故中心管$x<a$内位移恒零。$\partial_t$在$r=0$单位归一，中心固有时为$\tau=t$。BD在静态patch代数上的KMS逆温度为

$$
\beta_0=2\pi/H.
\tag{2}
$$

扩展读口其它位置的局部固有时有不同lapse；本稿使用统一中心时间，不把整个空间包称为同一随身钟。

|旧成果|本稿直接复用|没有重复领取|
|---|---|---|
|[1232 §2—5](../../../../archive_1219_1233/1232/proof.md)|正谱测度、有限开关及完整Dyson概率定理|原Dirac数值谱下界或实际设备合同|
|[dS参考候选 §4.2—5](../../thermal_mechanism_admission_after1270/actual_thermal_reference_candidate_after1270.md)|BD温标及同源扰动的候选范围|实际测温认证|
|[因果校准来源包 §2—4、6](../README.md)|合格局域过程的非选择来源差与条件几何运输|由抽象有界门自动取得完整来源|
|[原脉冲 §2—3](../../finite_horizon_transfer_after1272/proof.md)|正常Weyl输入、最小来源和中心空管|另一材料或新几何误差|

本次净增：**原BD标量的非零紧支Cauchy读口具有连续非零正频带，因而这同一材料满足1232有限正响应窗口的谱前件。** KMS、Bochner、Gaussian与Dyson工具均为成熟工具。

## 2. Cauchy涂抹与全部有限时差

固定$0<x_c<a$，置

$$
R_c=H^{-1}\tanh(Hx_c),\qquad F_c=1-H^2R_c^2>0.
\tag{3}
$$

在静态Cauchy片$\Sigma_0$取实非零$h\in C_c^\infty$，支撑严格位于$r<R_c$。采用$d\Sigma=F^{-1/2}r^2dr\,d\Omega_2$，定义

$$
A=\Phi_{\Sigma_0}(h)=\int h\Phi\,d\Sigma,\quad
v=\omega_{\rm BD}(A^2),\quad
B=\sin A=\frac{e^{iA}-e^{-iA}}{2i}.
\tag{4}
$$

这是正规自由KG场的Cauchy Weyl构造，不是点场的正弦。$A$采用标准自伴场实现，$\|B\|\le1$。Hadamard二点函数可拉回类空Cauchy片；光滑紧支空间涂抹给$v<\infty$。取共轭动量测试$k=h$，则

$$
[\Phi_\Sigma(h),\Pi_\Sigma(h)]
=i\int h^2d\Sigma\ne0.
$$

共轭动量同样光滑涂抹且方差有限，不确定性不等式因而给$v>0$。

记静态流为$\alpha_t$，$A_t=\alpha_t(A)$、$B_t=\alpha_t(B)$，令$w_h(t)=\omega_{\rm BD}(A_tA)$。BD表示中的静态流强连续且保参考，向量$A\Omega_{\rm BD}$存在，所以

$$
w_h\text{连续},\qquad |w_h(t)|\le v\quad(t\in\mathbb R).
\tag{5}
$$

这是Hilbert空间矩阵系数的结论，不要求Weyl算符范数连续，也不把有限时差的点核奇性作未经许可的绝对值积分。

零均值准自由态的Weyl乘积公式为

$$
\omega(e^{isA_t}e^{iuA})
=\exp\!\left[-\frac{s^2+u^2}{2}v-su\,w_h(t)\right],
\qquad s,u\in\mathbb R.
\tag{6}
$$

$w_h$包含对易子的虚部，不能只用对称二点。展开两正弦给精确恒等式

$$
\boxed{W_B(t)=\omega(B_tB)=e^{-v}\sinh w_h(t),\qquad
\omega(B)=0,\qquad W_B(0)=\frac{1-e^{-2v}}2>0.}
\tag{7}
$$

## 3. 原BD核与显式大时间界

[Fröb—Much—Papadopoulos，附录A式(A20)—(A21)](https://arxiv.org/pdf/2310.12185)给BD二点函数。本稿取$n=4,m^2/H^2=2$，对应$\mu=1/2$；由$\,{}_2F_1(2,1;2;z)=(1-z)^{-1}$得到

$$
G^+(X,X')=\frac{H^2}{8\pi^2(1-Z_\epsilon(X,X'))}.
\tag{8}
$$

$\epsilon$指同一BD边界值处方。原文§III给楔流KMS归一，换为中心时间即(2)。以下阈值自行由该核推导。

静态时间差$t$及空间夹角$\gamma$满足

$$
Z=\sqrt{F(r)F(r')}\cosh(Ht)+H^2rr'\cos\gamma.
\tag{9}
$$

两个点均在$r,r'\le R_c$时，

$$
Z-1\ge F_c\cosh(H|t|)-(1+H^2R_c^2).
$$

事前定义

$$
t_*=\frac1H\operatorname{arcosh}
\!\left(\frac{2(1+H^2R_c^2)}{F_c}\right),\qquad
C_h=\frac{H^2\|h\|_{L^1(d\Sigma)}^2}{2\pi^2F_c}.
\tag{10}
$$

当$|t|\ge t_*$，有

$$
Z-1\ge\frac{F_c}{2}\cosh(H|t|)
\ge\frac{F_c}{4}e^{H|t|}>0,\qquad
|w_h(t)|\le C_he^{-H|t|}.
\tag{11}
$$

此时紧支区域之间已无类光奇点，可以对普通函数核作绝对值估计。有限时差则始终使用(5)，不穿越点核奇性作此估计。

由$|\sinh z|\le\sinh|z|$、$|w_h|\le v$及$\sinh y\le y\cosh v$（$0\le y\le v$），

$$
|W_B(t)|\le e^{-v}\cosh(v)|w_h(t)|\le|w_h(t)|.
$$

故有明确有限界

$$
\boxed{\|W_B\|_{L^1(\mathbb R)}
\le 2t_*v+\frac{2C_h}{H}e^{-Ht_*}<\infty.}
\tag{12}
$$

$v,C_h,t_*$是给定$h,H,x_c$的解析常数；本稿没有声称已数值计算该包的方差或谱。

## 4. 连续非零正频带

按1232的Fourier符号，Bochner正测度满足

$$
W_B(t)=\int e^{i\nu t}\mu(d\nu),\qquad
\mu\ge0,\quad \mu(\mathbb R)=W_B(0)\le1.
\tag{13}
$$

(12)使这份有限正测度具有连续、有界、趋于零的密度

$$
\mu(d\nu)=p(\nu)d\nu,\qquad
p(\nu)=\frac1{2\pi}\int e^{-i\nu t}W_B(t)\,dt\ge0.
\tag{14}
$$

$W_B(0)>0$保证$p$不恒零，并排除仅零频原子承担全部质量的例外。KMS给

$$
p(-\nu)=e^{\beta_0\nu}p(\nu).
\tag{15}
$$

该式先在测度意义成立，再由连续性延为每点成立。因此存在$E>0$使$p(E)>0$：非零点若为负频可用(15)反射，若为零频则用连续性取相邻正频点。于是存在$c>0$和$0<r_0<E$，使

$$
p(\nu)\ge c\quad\text{当 }|\nu-E|\le r_0.
\tag{16}
$$

这是同一非零$h$的谱侧存在结论，不是另输入一份替代谱；没有指定数值$E,c,r_0$。

## 5. 复用1232：有限概率及正容差

取$0<r\le r_0$、有限$T\ge1/r$，由(16)有

$$
m_E=\mu([E-1/T,E+1/T])\ge\frac{2c}{T}.
\tag{17}
$$

沿原tent开关$\chi_T(\tau)=(1-|\tau|/T)_+$，令$c_0=(23/24)^4$。1232 §3给可选的领先响应下界与外带比

$$
A_T=2c_0cT,\qquad
u_T=\frac{16}{r^4T^2A_T}=\frac8{c_0cr^4T^3}.
\tag{18}
$$

先减小$r$使$\beta_0r$够小，再取充分大但有限的$T$使$u_T$够小，即可令谱误差$\beta_0r+\log(1+u_T)$任意小。这没有把无限时跃迁率当有限记录。

作为明确的有界比较过程，复用纯翻转交互

$$
V_I(\tau)=\lambda\chi_T(\tau)
(\sigma_+e^{iE\tau}+\sigma_-e^{-iE\tau})\otimes B_\tau.
\tag{19}
$$

有限两能级自由能隙为$E$，场后态完整保留。$\|B_\tau\|\le1$保证Dyson酉存在；tent在此有界比较中合法。若实际控制采用其它平滑开关，其差另计操作误差。

令$x_\lambda=|\lambda|T$、$R_3=\sinh(x_\lambda)-x_\lambda$。原1232 §4奇数翻转阶数证明原样给出

$$
q_\pm\ge\lambda^2A_T,\qquad
|P_\pm^\chi-q_\pm|\le d_D=2x_\lambda R_3+R_3^2.
\tag{20}
$$

任取$0<\rho<1$，可选非零有限$\lambda$满足

$$
0<\lambda^2\le\min\{T^{-2},\,\rho A_T/T^4\}.
\tag{21}
$$

原$x_\lambda\le1$护栏给$d_D\le x_\lambda^4/2\le\rho\lambda^2A_T/2$，所以两份完整比较概率严格正，且留有正的操作误差容差：

$$
0\le\epsilon_{\rm op}\le\rho\lambda^2A_T/2.
\tag{22}
$$

**若**同对象实际过程另外证明$|P_\pm^{\rm act}-P_\pm^\chi|\le\epsilon_{\rm op}$，则

$$
\left|\log\frac{P_-^{\rm act}}{P_+^{\rm act}}-\beta_0E\right|
\le\beta_0r+\log(1+u_T)+\log\frac{1+\rho}{1-\rho}.
\tag{23}
$$

$\rho$亦可取有限正小量。因此有限窗口、非零比较概率和正容差能够共存。式(22)不证明已有设备达到该容差；小耦合时容许绝对误差也缩小。没有沿用1232原Dirac数值表或其具体核心谱质量。

## 6. 中心因果域、来源及实际实现边界

整个比较窗可平移到$[t_c-T,t_c+T]\subset(0,\infty)$，例如$t_c=T+H^{-1}$，以$\tau=t-t_c$写相位与开关。这只是共同时间原点改变，时间间隔仍为中心固有时；BD平稳性使谱和比较概率不变。

$h_t$由静态流运输，始终位于$x<x_c<a$。Cauchy Weyl元属于相应球的局域代数，time-slice性质允许在任意足够薄的时间邻域表示。$h$支撑与$x_c,a$之间有正间距，窗口起点也严格大于0，因此可以统一选择薄邻域，使整个有限调用仍落在原中心空管。任意上述有限$T$均可如此；双侧$t$的相关函数只是平稳谱工具，不要求在脉冲前的负时刻执行校准。

原位移在该管内为零。若实际实现另满足[因果校准来源包§2](../README.md)的共同作用、时间前缀、控制、来源变分及定义域前件，则其非选择装置记录在脉冲／BD两输入中相同，且可继承该包§4的完整平均来源增量。这里只运输已经资格化的过程。

必须区分：

1. $B=\sin\Phi_\Sigma(h)$的局域代数归属不自动提供许可的局部材料耦合。**空间整体涂抹的正弦怎样实际执行**仍是待证的有限有效合同，不由(19)造出装置。
2. 抽象比较酉不提供内部开关、完整Hilbert来源或实际末读。全部绝对参考、探针、控制、记录及交互来源保持；条件选择分支不自动继承非选择来源身份。
3. 来源始终取(1)的最小耦合定义，不能利用固定背景传播相同而换成共形源。
4. 加入几何反作用后，传播支撑、参考、实际钟及谱可能改变。本稿未取得实际温标、Wald热身份或固定非零$G$的共同误差。

## 7. 交付与停止线

本次已证同一原BD材料、中心时间和中心因果域中存在非零连续正频带，并直接接成有限比较概率与正容差。有限时差由正规Cauchy涂抹和强连续矩阵系数控制，大时间由明确BD核控制。

没有谱数值、最短等待时间、新材料、实际设备、制备工厂或新认知公理。本项到此停止；不展开谱扫描、探针优化或几何回授工程。

