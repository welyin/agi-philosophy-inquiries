# 固定温度重KG参考：平均源与完整热方差的有限界

2026-10-18。1291作者证明；[准入](../_shared/round1291_admission.md)、[原候选与查重](../_admission/heavy_thermal_screen/selection.md)。正式独审与验收另行登记。下述 $\beta>0$、测试、长度 $\ell$ 均事前固定；只改变候选物种质量 $M>0$。

## 1. 同一场、参考、来源与定义域

自然单位 $\hbar=c=k_B=1$，$\eta=\operatorname{diag}(-1,1,1,1)$。最小实KG作用与Hilbert应力为

$$
S_M=-\frac12\int\sqrt{-g}\,[g^{\mu\nu}\partial_\mu\phi\partial_\nu\phi+M^2\phi^2],
\qquad
T_{\mu\nu}=\partial_\mu\phi\partial_\nu\phi
-\frac12g_{\mu\nu}[(\partial\phi)^2+M^2\phi^2].
\tag{1}
$$

始终用各质量的Minkowski真空正常序 $:T:_M$。这是公开一致的真空均值校准规则，不是跨质量固定裸反项；更不是随温度另扣热均值。

记 $d\mu(p)=d^3p/(2\pi)^3$、$E_p=\sqrt{|p|^2+M^2}$、$n_p=(e^{\beta E_p}-1)^{-1}$。选相对惯性时间平移的标准零均值准自由KMS态，其二点为

$$
W_{M,\beta}(x-y)=\int\frac{d\mu(p)}{2E_p}
\left[(1+n_p)e^{-iE_p(t_x-t_y)+ip\cdot(x-y)}
+n_pe^{iE_p(t_x-t_y)-ip\cdot(x-y)}\right].
\tag{2}
$$

空间指数中的 $x-y$ 指三维部分。式(2)的正性、CCR及KMS边界可逐频率检查。$n_p$的指数高能衰减使它与真空核之差光滑，因此保Hadamard奇性。所选态的局部正常性也可直接复用 [Sanders，1209.6068v2，定理6.2(ii)、(iv)及其证明](https://arxiv.org/pdf/1209.6068)。不假设无限空间上存在真空Fock正常的全局Gibbs密度矩阵。

对实对称 $f^{\mu\nu}\in C_c^\infty(\mathbb R^4)$，定义Wick域上的对称二次形式及其GNS二阶矩

$$
A_M(f)=\sum_{\mu,\nu}\int f^{\mu\nu}:T_{\mu\nu}:_M\,d^4x,\quad
\mu_\beta(f)=\omega_{M,\beta}(A_M(f)),\quad
V_\beta(f)=\omega_{M,\beta}((A_M(f)-\mu_\beta(f))^2).
\tag{3}
$$

以下积分给出有限二阶矩的直接证明，不从一个未经证明的全局矩阵迹推断。双指标求和包含两次非对角项；单个 $T_{12}$ 测试 $h$ 应取 $f^{12}=f^{21}=h/2$。数学源测试未被授予实际测量仪器。

## 2. 热均值：不因真空减除而为零

式(2)和完整(1)给平移不变的静止系均值

$$
\omega(:T_{00}:_M)=\rho_{M,\beta}
=\int d\mu(p)\,E_pn_p,\qquad
\omega(:T_{ij}:_M)=p_{M,\beta}\delta_{ij},\qquad
p_{M,\beta}=\frac13\int d\mu(p)\frac{|p|^2}{E_p}n_p.
\tag{4}
$$

混合分量为零，$0\le p_{M,\beta}\le\rho_{M,\beta}/3$。设 $d_M=1-e^{-\beta M}$。对 $E\ge M$，

$$
n(E)\le\frac{e^{-\beta E}}{d_M},\qquad
p^2dp=pE\,dE\le E^2dE.
\tag{5}
$$

故

$$
\rho_{M,\beta}\le
\frac{e^{-\beta M}}{2\pi^2d_M}
\left(\frac{M^3}{\beta}+\frac{3M^2}{\beta^2}
+\frac{6M}{\beta^3}+\frac6{\beta^4}\right)
=\frac{e^{-\beta M}}{2\pi^2d_M}P_3(M,\beta).
\tag{6}
$$

这是分部积分的上界，未把精确密度积分误写为该多项式。令 $B_0(f)=\sum_{\mu,\nu}\|f^{\mu\nu}\|_1$，则 $|\mu_\beta(f)|\le B_0(f)\rho_{M,\beta}$。

## 3. Wick方差的两类不同通道

令未来壳协矢 $p_\mu=(-E_p,p)$，$p\cdot p=-M^2$。定义

$$
t_{\mu\nu}(p,q)=-\frac12(p_\mu q_\nu+p_\nu q_\mu)
+\frac12\eta_{\mu\nu}(p\cdot q-M^2),
$$

$$
s_{\mu\nu}(p,q)=t_{\mu\nu}(p,-q)
=\frac12(p_\mu q_\nu+p_\nu q_\mu)
-\frac12\eta_{\mu\nu}(p\cdot q+M^2).
\tag{7}
$$

这里 $t$ 是成对创建/湮灭顶点，$s$ 是数保持顶点。两者都含完整质量项，且直接满足

$$
(p+q)^\mu t_{\mu\nu}=0,\qquad
(p-q)^\mu s_{\mu\nu}=0,\qquad
|t_{\mu\nu}|,\ |s_{\mu\nu}|\le3E_pE_q.
\tag{8}
$$

因此来源守恒是同一场方程的恒等式，不是另加噪声守恒假设。定义
$\widehat f(E,P)=\int f(t,x)e^{iEt-iP\cdot x}\,dtd^3x$，并令

$$
F_+(p,q)=\sum_{\mu,\nu}t_{\mu\nu}(p,q)
\widehat f^{\mu\nu}(E_p+E_q,p+q),\qquad
F_-(p,q)=\sum_{\mu,\nu}s_{\mu\nu}(p,q)
\widehat f^{\mu\nu}(E_p-E_q,p-q).
\tag{9}
$$

按 $[a(p),a^\dagger(q)]=(2\pi)^3\delta^3(p-q)$ 和模系数 $1/\sqrt{2E_p}$，对称成对项的核是 $F_+/(2\sqrt{E_pE_q})$；数保持项有两份交叉收缩，其核是 $F_-/\sqrt{E_pE_q}$。Wick展开并减去均值平方，得到

$$
\begin{aligned}
V_\beta(f)&=V_{\mathrm{pair},\beta}(f)+V_{\mathrm{sc},\beta}(f),\\
V_{\mathrm{pair},\beta}(f)
&=\frac12\int d\mu(p)d\mu(q)\,
\frac{|F_+(p,q)|^2}{E_pE_q}
\big[(1+n_p)(1+n_q)+n_pn_q\big],\\
V_{\mathrm{sc},\beta}(f)
&=\int d\mu(p)d\mu(q)\,
\frac{|F_-(p,q)|^2}{E_pE_q}\,n_p(1+n_q).
\end{aligned}
\tag{10}
$$

两类交叉项为零；$a^\dagger a$项的断开收缩已经进入被扣除的均值平方。成对收缩的因子2没有再次乘进散射项。单振子检查为
$\operatorname{Var}(A a^2+\bar A a^{\dagger2}+B a^\dagger a)
=2|A|^2[(n+1)^2+n^2]+B^2n(n+1)$。

**散射项不具 $2M$ 频率阈值。** 当 $p=q$ 时 $E_p-E_q=0$，而 $s_{00}(p,p)=E_p^2$。其邻域可对固定低频测试贡献非零。不能把整个(10)直接套入真空的总频率估计。

## 4. 成对项保留真空界

令 $n_M=(e^{\beta M}-1)^{-1}$，则

$$
V_0(f)\le V_{\mathrm{pair},\beta}(f)
\le(1+2n_M+2n_M^2)V_0(f).
\tag{11}
$$

这里 $V_0$ 就是(10)令 $n_p=0$的真空项。对任意整数 $r>4$，令 $B_r(f)=\sum_{\mu,\nu}\|\partial_t^r f^{\mu\nu}\|_1$。由(8)及时间分部积分，

$$
V_0(f)\le \frac{9B_r(f)^2}{8\pi^4}
\int_M^\infty dE\int_M^\infty dF\,
\frac{E^3F^3}{(E+F)^{2r}}
\le\frac{9B_r(f)^2}{1120\pi^4(2r-8)}(2M)^{8-2r}.
\tag{12}
$$

末步使用 $\int_0^1x^3(1-x)^3dx=1/140$。这复用1290的真空估计而保其完整归一。温度固定不使真空项乘上 $e^{-\beta M}$；小于等号的热增量另由(11)及下一节控制。

## 5. 无能隙散射项的空间平滑与热占据界

选择任一事前固定长度 $\ell>0$，定义

$$
C_\ell(f)=\sum_{\mu,\nu}
\|(1-\ell^2\Delta_{\mathbf x})^2f^{\mu\nu}\|_1.
$$

在空间分部积分且不依赖时间频率，得

$$
|F_-(p,q)|\le
\frac{3E_pE_qC_\ell(f)}
{(1+\ell^2|p-q|^2)^2}.
\tag{13}
$$

另一方面

$$
n_p(1+n_q)\le\frac{e^{-\beta E_p}}{d_M^2},\qquad
E_q\le E_p+|q-p|.
\tag{14}
$$

后一式来自 $|\nabla_pE_p|\le1$，因此对任意两个动量成立。置 $k=q-p$，径向Beta积分给

$$
J_0(\ell)=\int\frac{d^3k}{(2\pi)^3}
\frac1{(1+\ell^2k^2)^4}=\frac1{64\pi\ell^3},
\qquad
J_1(\ell)=\int\frac{d^3k}{(2\pi)^3}
\frac{|k|}{(1+\ell^2k^2)^4}=\frac1{24\pi^2\ell^4}.
\tag{15}
$$

所以

$$
V_{\mathrm{sc},\beta}(f)\le
\frac{9C_\ell(f)^2}{d_M^2}
\left[J_0\int d\mu(p)E_p^2e^{-\beta E_p}
+J_1\int d\mu(p)E_pe^{-\beta E_p}\right].
\tag{16}
$$

定义
$P_j(M,\beta)=\sum_{a=0}^j j!M^{j-a}/[(j-a)!\beta^{a+1}]$。
沿(5)的径向上界及有限次分部积分，

$$
\int d\mu(p)E_p^ke^{-\beta E_p}
\le\frac{e^{-\beta M}}{2\pi^2}P_{k+2}(M,\beta).
$$

合并即得完整显式散射预算

$$
\boxed{
V_{\mathrm{sc},\beta}(f)\le
\mathcal S_{M,\beta,\ell}(f):=
\frac{9C_\ell(f)^2e^{-\beta M}}{2\pi^2d_M^2}
\left[
\frac{P_4(M,\beta)}{64\pi\ell^3}
+\frac{P_3(M,\beta)}{24\pi^2\ell^4}
\right].}
\tag{17}
$$

其中
$P_4=M^4/\beta+4M^3/\beta^2+12M^2/\beta^3+24M/\beta^4+24/\beta^5$。
本证明保留高动量、近零能量转移和刺激因子 $1+n_q$；没有假设热散射的全部动量都被某个固定频带截去。

## 6. 固定有限菜单的量词和共同热参考

由(6)、(11)、(12)、(17)，对每一固定 $\beta>0$、固定有限测试族 $\{f_j\}$ 和固定 $\ell>0$，

$$
\mu_\beta(f_j)\longrightarrow0,\qquad
V_\beta(f_j)\longrightarrow0
\qquad(M\longrightarrow\infty).
\tag{18}
$$

对任意指定正均值/方差容差，右侧显式函数给有限质量门槛，使所有有限测试同时满足要求。不是先改变测试再选结果。更具体地，

$$
0\le V_\beta(f)-V_0(f)
\le(2n_M+2n_M^2)V_0(f)+\mathcal S_{M,\beta,\ell}(f).
\tag{19}
$$

连接交叉协方差由GNS的Cauchy–Schwarz不等式控制：
$|\omega((A_i-\mu_i)(A_j-\mu_j))|\le\sqrt{V_iV_j}$。
未经居中的二阶矩也满足 $\omega(A_i^2)=V_i+\mu_i^2\to0$。没有给不对易源测试赋予一份免费联合经典概率分布。

若母参考为同一物理时间的 $\beta$-KMS态，零门户下总自同构为两因子的乘积，则在解析元素上KMS边界函数相乘，给乘积参考仍为同一 $\beta$ 的KMS态。实际母输入及其未知 $R$ 可保持原合同；新增参考独立附加时，固定背景、零门户的母记录过程不变，新增源为(1)，均值和居中协方差按乘积态相加。没有把未知非平衡输入称作KMS，也没有据此授予所有热操作。

## 7. 来源完整性与不作的外推

- 真空与热态用同一 $:T:_M$。热均值为(4)而非零；$\partial^\mu T_{\mu\nu}=0$来自KG方程，(8)同时核两类顶点。确定的c-number匹配项不改连接方差，但可改变均值，须保所声明处方。
- $\beta$及惯性时钟归一固定。不把本结论搬为曲域、加速观察者、HG或瞬时热参考结论；也不允许用 $\beta(M)\to0$ 保持热占据后仍套(18)。
- 无限空间热参考一般不是全局有限能源准备；此处只给局部正常参考和固定紧时空源矩。若任务要求有限体积制备/边界，须另核该任务，不能直接用本式给其设备资格。
- 没有相互作用、开关或装置来源被隐去：候选本来是自由场零门户。若加入门户/控制/仪器，应把其全部应力和相互作用矩重新纳入；本式没有验证这些项。
- 新增正热源不使Minkowski自动成为自引力解。此处是固定背景源资格，不是全Einstein共同反作用、所有高阶噪声或任意细分辨率的预算。
- 若公开有效类别另固定质量上限 $\Lambda$ 或总资源界，有限门槛仍须与它相容；不得声称一个固定 $\Lambda$ 的EFT类别允许任意小误差。数学质量极限不删除该约束。

这为同温参考的新增物种基线提供具体可检验后继，不决定实际物种数、质量、三代或门户允许性。本轮按已登记的1291准入提交作者证据；该结果不替代主线的结构选择判断或长程验收。

## 8. 一份真正紧时空、固定同温的有限预算

沿1290已认证的盒核方法，令 $u_h=h^{-1}1_{[-h/2,h/2]}$、$h=1/8$，事前固定非负归一 $b\in C_c^\infty(\mathbb R)$，支撑在 $[-1/8,1/8]$。例如取归一后的 $\exp[-1/(1-64t^2)]$ 于开区间内，区间外取零。置

$$
g=u_h^{*6}*b,\qquad v=u_h^{*4}*b,\qquad
s(x)=v(x_1)v(x_2)v(x_3),\qquad f^{00}(t,x)=g(t)s(x),
\tag{20}
$$

其余张量分量为零。这里 $g,s$ 非负、归一且光滑；时间支撑为 $[-1/2,1/2]$，各空间坐标支撑为 $[-3/8,3/8]$。这不是Gaussian诊断或截断后的近似紧窗。

把各导数落在不同盒核上，$u_h'$ 的总变差为 $2/h=16$，因此

$$
B_0(f)=1,\qquad B_6(f)\le16^6,\qquad
\|\partial_{x_i}^k v\|_1\le16^k\quad(0\le k\le4).
$$

末项的 $v$ 是单轴函数，$\partial_{x_i}$指对应轴导数；对 $s$ 的混合导数取三因子L¹范数的积。令固定 $\ell=1$，展开空间算符给

$$
C_1(f)\le1+6\cdot16^2+9\cdot16^4
=(1+3\cdot16^2)^2=769^2=591361.
\tag{21}
$$

其中 $\Delta^2$ 有9个按顺序计的四阶项，交叉项没有漏掉。固定 $\beta=1$、$M=2^{20}$，令 $q_*=2^{-256}$。只用 $M\ge256$ 与 $e>2$，就有 $e^{-M}<q_*$、$d_M>1/2$、$n_M<2q_*$，且 $1+2n_M+2n_M^2<2$。无需计算巨大指数或依赖浮点underflow。

用 $\pi>3$，(6)、(12)、(17)分别给严格上界

$$
\rho_{M,1}<
\frac{q_*P_3(M,1)}9<10^{-59},
\qquad
V_0(f)<\frac1{40320\,2^{36}},
$$

$$
\mathcal S_{M,1,1}(f)<
2(769^2)^2q_*
\left[\frac{P_4(M,1)}{192}+\frac{P_3(M,1)}{216}\right]
<10^{-43},
$$

$$
\boxed{\quad|\mu_1(f)|=\rho_{M,1}<10^{-59},\qquad
V_1(f)<
\frac2{40320\,2^{36}}+
2(769^2)^2q_*
\left[\frac{P_4(M,1)}{192}+\frac{P_3(M,1)}{216}\right]
<10^{-15}.\quad}
\tag{22}
$$

全部比较可由有理数作严格核验。最后方差上界约为 $7.2183\times10^{-16}$；小数仅展示量级，严格结论用(22)。热均值没有设零，散射正项没有丢失。

测试、温度及长度由(20)—(21)固定后才核此质量。若另有事前固定有效上限 $\Lambda$ 且要求 $M\le\Lambda$，本具体实例只在 $2^{20}\le\Lambda$ 时可用；更一般门槛同样须落域，不能把固定 $\Lambda$ 类别说成允许任意小容差。此处单位是已声明的场/时间单位，不额外指定Planck质量或弱引力参数。

## 9. 复算和范围

[check.py](check.py)默认只读重算并比较[results.json](results.json)，首次显式使用--write排他创建。作者一组检查包括：

1. 从完整张量顶点核两种Ward身份、反向动量校准，以及删除质量项的失败见证。
2. 从单振子数态跃迁与几何分布矩核Wick因子；删除真空项、相同粒子因子或热散射项均有非零误差。零频数保持例明确不能套 $2M$ 阈值。
3. 用Beta/Gamma有理系数核(15)，以分部积分递推核 $P_3,P_4$。
4. 盒核导数的有限有符号测度、混合空间导数计数和(22)严格有理预算。

这些有限检查辅助核对实现，连续量词由§1—7解析证明承担；不作质量扫描或另造仪器。Sanders原文定理6.2及其局部正常性证明已实际核读。先前候选的内存常数核对不另计正式组；本轮独立代码/审阅及最终验收另行登记。
