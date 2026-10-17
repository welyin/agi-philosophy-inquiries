# 1264来源证明：非恒定共形KG响应、正频率运输与原记录的有限共同界

2026-10-16。正式预约1264后的来源交付，source_author_groups=1；新增认知公理0，不改变现目标、1263或其它冻结结果。最终联合验收以本轮acceptance记录为准。本次只使用1262已公开的原864维最小KG材料，实现一份**给定静态非均匀几何中的共同来源—传播—记录连接**。几何不是由材料自生的；不沿用1263 AW散射为同一实现，也不要求它们精确相等。

[作者复算](source_check.py)默认只读比对[source_results.json](source_results.json)，首次写入须显式使用 --write-results，以排他模式创建。代码核明确profile导数上界与共同预算；不模拟场传播，不作Fourier数值积分。

## 1. 查重与净增量

- [1186](../../archive_1161_1200/1186/proof.md)已证明同母作用的背景接触与有限记录相容，但没有本次原质量材料的非恒定静态KG来源运输。
- [1262](../1262/proof.md)提供同一最小KG一粒子应力核及有限涂抹；[1263](../1263/trace_source.md)把leading物质迹插入运到已给散射态。二者未计算本次同母作用的平方根生成元对非恒定度规的响应。
- 942原生场嵌入、1194的能源定位与[Moretti §3.2、§5.1](https://arxiv.org/pdf/2304.02133)用于身份/范围复用，不重新计作结果。紧支的一粒子位置振幅不等于紧支KG场数据或实际局域仪器。
- 1263后已冻结的常lapse/conformal方向筛查只给全局身份。本次不是任意扩展旧H[g]，而从1262同一个最小作用固定非恒定背景族及其来源。

新增连接是：正确区分固定Cauchy场图与固定正频率L²图以后，两种来源响应核有明确运输差；这份差与有限背景幅度余项共同控制**相同准备和原记录效果**的预测，且可以小于一个已指定的记录概率差。

## 2. 同一作用给出的静态正生成元

保原全部M及简并、跨质量相干、未知资料和任意参考：

$$
M=\sum_\alpha m_\alpha P_\alpha,\qquad 26<M<32.
$$

对静态度规族和同一最小作用，

$$
g_\lambda=e^{2\lambda f(x)}\eta,\qquad
S=-\int\sqrt{-g}\left(g^{\mu\nu}\partial_\mu\phi^\dagger
\partial_\nu\phi+\phi^\dagger M^2\phi\right)d^4x,
\quad \xi=0,
$$

令场变量 $\psi=e^{\lambda f}\phi$。空间分部积分给

$$
S=\int\left(\|\dot\psi\|^2-\langle\psi,K_\lambda\psi\rangle\right)dt,
$$

$$
K_\lambda=-\Delta+M^2e^{2\lambda f}
+\lambda\Delta f+\lambda^2|\nabla f|^2,\qquad
H_\lambda=\sqrt{K_\lambda}.
\tag{1}
$$

这里 $\psi$ 是场变量；一粒子L²振幅另记u，初始正频率场数据为
$\psi=(2H_\lambda)^{-1/2}u$。不能混用两者的支撑。

对 $0\le f\le1$、$\lambda\ge0$，

$$
K_\lambda=(\nabla-\lambda\nabla f)^\dagger
(\nabla-\lambda\nabla f)+M^2e^{2\lambda f}\ge M^2.
\tag{2}
$$

各质量部门的势是有界实函数，K的共同域为H²，正平方根H的域为H¹；相对 $H_0$ 的差有界。式(2)不依赖删掉正常态动量尾。

## 3. 两个不同的正确来源核

以下在一个质量m部门计算。记

$$
d=p-q,\quad \omega=E_p-E_q,\quad
s=\sqrt{E_pE_q},\quad b=(E_p+E_q)/2,\quad
V=2m^2-|d|^2,\quad E_p=\sqrt{m^2+|p|^2}.
$$

Fourier约定沿用1262：$\widehat f(d)=(2\pi)^{-3}\int e^{-id\cdot x}f(x)d^3x$。
固定L²图的 $A_f=H'_0$ 由Sylvester方程给出：

$$
H_0A_f+A_fH_0=2m^2f+\Delta f,\qquad
A_f(p,q)=\widehat f(d)\frac{V}{2b}.
\tag{3}
$$

固定原Cauchy正则场数据 $(\phi,\Pi)$ 的Hamiltonian是

$$
\mathcal H_\lambda=\int
\left[e^{-2\lambda f}\Pi^\dagger\Pi+
e^{2\lambda f}|\nabla\phi|^2+
m^2e^{4\lambda f}|\phi|^2\right]d^3x.
$$

其导数的一粒子核正是1262物理负迹 $Q_f=-T^\mu{}_\mu(f)$：

$$
Q_f(p,q)=\widehat f(d)
\frac{\boldsymbol p\cdot\boldsymbol q+2m^2-E_pE_q}{s}
=\widehat f(d)\frac{V+\omega^2}{2s}.
\tag{4}
$$

因此不能写成 $A_f=Q_f$。例如 $p=0,q=me_1$ 时，去掉共同Fourier权，

$$
A=m(\sqrt2-1),\qquad Q=m(2-\sqrt2)/2^{1/4},
$$

两者不等。等能源矩阵元 $\omega=0$ 才相等。

## 4. 差项来自场重标定与正频率运输

固定重标定场图 $(\psi,\pi)$ 的场Hamiltonian导数，一粒子核先给
$\widehat f V/(2s)$。于是两步差分别为

$$
Q_f-\frac{\widehat f V}{2s}
=\widehat f\frac{\omega^2}{2s},\qquad
\frac{\widehat f V}{2s}-A_f
=\widehat f\frac{V\omega^2}{8bs(b+s)}.
$$

第一项来自原场的正则重标定；第二项来自正频率归一化图。合写为

$$
Q_f-A_f=[H_0,[H_0,C_f]],\qquad
C_f(p,q)=\widehat f(d)
\left[\frac1{2s}+\frac{V}{8bs(b+s)}\right].
\tag{5}
$$

这不是任意添加的修正。等时将新正频率数据投影到旧正频率分量的映射为

$$
T_\lambda=\frac12\left[
H_0^{1/2}e^{-\lambda f}H_\lambda^{-1/2}
+H_0^{-1/2}e^{\lambda f}H_\lambda^{1/2}
\right],
$$

在共同光滑核上直接微分得到

$$
T'_0=-[H_0,C_f]=:-S_f,\qquad
A_f+[T'_0,H_0]=Q_f.
\tag{6}
$$

式(6)最后一式只陈述基点微分的图运输身份，不要求另证T的全幅度可逆或酉性质；**T本身不是已证明精确酉的比较**，完整正负频率分解还含负频部分。$S_f^\dagger=-S_f$；若使用 $e^{-\lambda S_f}$，它只是一个明确的一阶酉运输选择，不能免费宣布等于T的全部有限幅度作用。

本次下节实际概率界不依赖把T当酉，也不换准备：直接比较同L²中的两个有界差生成元。若另使用图运输，准备、原效果和后态都必须一起运输，不能暗中重准备。这里是不同**静态背景**之间的比较，不是 $\lambda(t)$ 开关过程；后者的额外时间导数和正负频率混合未获授权。

## 5. 紧支测试的统一来源算符界

定义

$$
C_j(f)=\int |d|^j|\widehat f(d)|d^3d.
$$

由 $b,s\ge m$、$|\omega|\le|d|$，式(5)给

$$
|Q-A|/|\widehat f|
\le \frac{5|d|^2}{8m}+\frac{|d|^4}{16m^3}.
$$

Schur检验因此给

$$
\boxed{\ \|Q_f-A_f\|\le
D_f:=\frac{5C_2(f)}{8\cdot26}
+\frac{C_4(f)}{16\cdot26^3}.\ }
\tag{7}
$$

这是对全部原质量、全部机械L²空间与任意参考统一的界，不是低动量压缩。各m的直和保持相干。另有
$\|C_f\|\le5C_0(f)/(4\cdot26)$、
$\|S_f\|\le5C_1(f)/(4\cdot26)$：
用 $|d|^2\le4b^2,m^2\le b^2$ 即得式(5)中括号的绝对值不超过 $5/(4s)$。

$A_f$ 有界，且 $\|A_f\|\le32+\|\Delta f\|_\infty/52$；故本次负迹Q也有界。这是这一特定质量正、有光滑空间权的迹插入结论，不推广为每个应力分量均有界。

## 6. 同一个明确平顶profile及有限常数

与记录候选共用唯一函数

$$
\eta(t)=\begin{cases}e^{-1/t},&t>0,\\0,&t\le0,\end{cases}
\qquad
\rho(s)=\frac{\eta(4-s)}{\eta(4-s)+\eta(s-1)},\qquad
F(z)=\rho(|z|^2),\quad f_R(x)=F(x/R).
\tag{8}
$$

$F=1$ 于单位球，支撑包含于半径2的闭球；$F\in C_c^\infty$，$0\le F\le1$。过渡内分母大于 $1/2$。原生作用只要求此明确光滑profile。下面直接由同一Fourier预算控制梯度和Hessian；不添加另一份profile正则性门槛。

不使用Gaussian的 $C_0=1$。缩放只给 $C_j(f_R)=R^{-j}C_j(F)$。
这些C是明确profile的有限可计算常数。对 $j=2,4$，Cauchy–Schwarz、Plancherel给

$$
C_j(F)^2\le\frac1{8\pi}
\sum_{r=0}^{j+2}\binom{j+2}{r}
\sum_{|\alpha|=r}\frac{r!}{\alpha!}
\|\partial^\alpha F\|_2^2.
\tag{10}
$$

其使用 $\int_{\mathbb R^3}(1+|d|^2)^{-2}d^3d=\pi^2$；
将 $|d|^{2j}(1+|d|^2)^2$ 上界为 $(1+|d|^2)^{j+2}$。
只需要到六阶的有限导数，不要求计算最优Fourier积分。

[来源程序](source_check.py)还给完全有限的有理上界 $\overline C_j\ge C_j(F)$：
写 $\eta^{(n)}=e^{-v}P_n(v),v=1/t$，递推
$P_{n+1}=v^2(P_n-P'_n)$；
逐项用 $\sup_{v>0}v^k e^{-v}=(k/e)^k<(k/2)^k$ 控制。
从 $(\eta(4-s)+\eta(s-1))\rho=\eta(4-s)$ 逐次求导，
再用 $\partial_i\rho^{(k)}(|z|^2)=2z_i\rho^{(k+1)}(|z|^2)$ 生成所有多指标上界。
支撑在 $[-2,2]^3$，故L²平方可由64倍相应sup平方控制。
将式(10)右侧记为 $S_j/(8\pi)$，取 $\overline C_j=S_j/24+1$ 即足。
这个选择很保守，只为证明一份明确有限任务存在，不做常数优化。递推实际给出

$$
\overline C_2=\frac{91062363423308015531}{3},\qquad
\overline C_4=\frac{610242449441917720898626279056227}{3}.
$$

这些是解析导数递推的上界，不是Fourier积分的数值估计。式(13)因此可明确取
$L=19100606975$、$R=38201213950$。保守常数造成的大尺度只属于这份充分证书；没有证明更小尺度不可行。

## 7. 有限幅度余项，不把一阶式冒充完整传播

采用下一节同一L条件后，来源差的第一项满足

$$
\frac{5C_2(f_R)}{208}\le\frac1{2000},
\qquad C_2(f_R)\le\frac{13}{625}=0.0208=:h.
\tag{9}
$$

Fourier逆变换因此同时给
$\|\Delta f_R\|_\infty\le h$ 和 $\|\operatorname{Hess}f_R\|_{\rm op,\infty}\le h$。
沿梯度方向作一维Taylor，使用 $0\le f_R\le1$ 得
$\|\nabla f_R\|_\infty^2\le2h=26/625=:g$。
这是同一来源预算的推论，不是额外强假设。

固定一个 $26<m<32$ 部门，$0\le\lambda\le\log(1.01)<1/100$。令

$$
B_1(m)=2m^2(1.01)^2+h+2(1/100)g,\qquad
B_2(m)=4m^2(1.01)^2+2g.
$$

有 $\|K'_\lambda\|\le B_1$、$0\le K''_\lambda\le B_2$。
Lyapunov逆
$\mathcal L_{H}^{-1}(X)=\int_0^\infty e^{-tH}Xe^{-tH}dt$
是正映射且范数不超过 $1/(2m)$。对
$HH'+H'H=K'$、
$HH''+H''H=K''-2(H')^2$，
先得到

$$
\|H'_\lambda\|
\le m(1.01)^2+\frac{h+2(1/100)g}{2m}<33.
$$

由于 $K''\ge0$，上下算符序分别给

$$
-\frac{33^2}{26}I < H''_\lambda
\le \left(2m(1.01)^2+\frac gm\right)I
<66I.
$$

下界绝对值小于42，上界小于66；所以 $\|H''_\lambda\|<66$。固定同一正频率L²图的Taylor余项满足

$$
\boxed{\ \|H_\lambda-H_0-\lambda A_f\|\le33\lambda^2.\ }
\tag{11}
$$

这是有界势族正平方根的算符范数微分：共同正下界和Sylvester逆控制H的一、二阶导数；不把形式链式法则用于非对易算符。

因此在同一准备和同一终末效果下，

$$
\|e^{-itH_\lambda}-e^{-it(H_0+\lambda Q_f)}\|
\le |t|(\lambda D_f+33\lambda^2).
\tag{12}
$$

两个生成元自伴，后者只用于有限比较；没有把它宣称为全幅度精确KG生成元。
式(12)覆盖原Mr、机械效果及任意参考的联合有限效应：纯化后任何 $0\le E\le I$ 的概率差不超过右侧。它不制造新的实际当地仪器。

## 8. 一份统一的非空有限参数合同

只作这一组准入，不再扫描profile或设备：

$$
r=10^{-4},\qquad \lambda=\log(1+r),\qquad T=\pi/2,\qquad R=2L,
$$

$$
L=\left\lceil\max\left\{
10^4,\
\sqrt{\frac{625\,\overline C_2}{52}},\
\left(\frac{125\,\overline C_4}{281216}\right)^{1/4}
\right\}\right\rceil .
\tag{13}
$$

$\overline C_j$是上一节同一个F的有限有理上界；所以式(13)定义实际有限L。式(7)的两项各不超过 $1/2000$，总和 $D_{f_R}\le1/1000$。

使用本轮[记录证明](record_proof.md)中的两份**已知准备任务**：三维正常cos⁴位置包半宽L，中心分别0与 $10Le_1$，同一原内部编码 $\delta=\pi/4$。前者全在 $f=1$ 平顶内，后者在 $f=0$ 区；两份任务各自保全部原未知载荷与参考。它们不是一个没有资源代价的免费路径标签效果。

同轮记录证明给出的同母传播对照为：每支包误差不超过 $18/(13L^2)$；
原Mr两份理想读数之差大于 $9r/16$。
对每支再使用式(12)，两任务的线性来源比较概率差下界为

$$
\Delta p_Q>
\frac{9r}{16}-\frac{36}{13L^2}
-4\left(rD_{f_R}+33r^2\right)
>\frac r2=\frac1{20000}.
\tag{14}
$$

最后一式只用 $L\ge10^4,D\le1/1000,\lambda<r,T<2$。
同一物理负迹Q、同一原材料、相同准备/效果的两种演化在这里形成了一条实际数学连接，而非两张独立的来源表与钟表。概率差、来源能源差和时间系数在式(12)—(14)中明确转换，没有直接相加不同单位。

这里引用本轮针对同一H和同一包的[记录证明](record_proof.md)；原生终末操作权限与全轮独审共同登记在验收中。来源程序只核本稿有限常数与预算，不重复积分或重做记录证明。

## 9. 交付范围与停止线

本轮交付的**限定中间结果**是：同一非恒定静态母作用固定来源与传播，补上正频率图差项，并在一份有限两任务菜单上留下严格正的共同误差余量。它避免1263后只证明global lapse或任意H[g]扩展的止线。

结论限于事前给定的静态共形几何及其一粒子有效实现。未证明动态源生成度规、完整自引力/场存能、普适物质或GR；未证明原Mr在场论意义上的实际当地测量与测后来源，也未授权时间脉冲。不能把这些未签项目抬成整个研究纲领的普遍前置门槛。

后继以本轮[NEXT](NEXT.md)和最终验收为准。不要据此继续制造测量仪器、优化巨大L或重启全QFT工程。
