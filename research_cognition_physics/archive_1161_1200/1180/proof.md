# 1180证明：点态正成本的阶数限制、反例和同一量子后态

2026-10-13。三组有限证书；内部独立审阅另列。本文检验新增候选条件，不改变现行公理，也不把传播比较器称为生成的GR。

## 1. 对象和准确量词

在三维平直空间上取实对称场 $h_{ij}$ 及canonical动量 $p^{ij}$，只要求线性空间规范变换

$$
\delta_\xi h_{ij}=\partial_i\xi_j+\partial_j\xi_i,\qquad
\delta_\xi p^{ij}=0.
\tag{1}
$$

低阶类 $\mathcal L_2$ 的Hamiltonian为

$$
H=\frac12\int d^3x\left[p:A:p+\sum_\alpha |L_\alpha h|^2\right].
\tag{2}
$$

$A$ 是严格正、超局域、常系数动能；$L_\alpha$ 为至多二阶空间导数的线性算子；每份平方严格逐点规范不变。总二次型平移、SO(3)不变，无 $p$—$h$ 混合；TT势非退化。没有加入微分动量、辅助约束消元、非局域核或额外变量。TT指沿非零波矢的横向无迹部分。它是本轮必要传播检验的扇区，不自动等于完整健康引力的物理态空间。

这里甚至不必额外假设每个线性 $L_\alpha h$ 已不变：若有限平方和在每点对全部 $\xi$ 不变，将 $h$ 换成 $h+t\delta_\xi h$，$t^2$ 系数是 $\sum_\alpha|L_\alpha\delta_\xi h|^2$；正性使每项为零。正平方之间不能相消规范变差。但**同一个不变量内不同导数阶可以交叉**；第4节正利用这一自由。

## 2. jet的完整线性审计

在一点，$\partial_i\xi_j$ 的对称部分任意，故可任意改变6个 $h_{ij}$。记

$$
A_{kij}=\partial_kh_{ij},\qquad
B_{i,jk}=\partial_j\partial_k\xi_i=B_{i,kj}.
$$

18维一阶jet的变差和逆映射为

$$
\delta A_{kij}=B_{j,ki}+B_{i,kj},\qquad
B_{i,jk}=\frac12\bigl(\delta A_{jik}+\delta A_{kij}-\delta A_{ijk}\bigr).
\tag{3}
$$

右边是线性Christoffel组合。零阶和一阶规范jet可独立指定，所以没有非零的零／一阶线性严格点态不变量。无边界条件限制的光滑 $\xi$ 可用局部Taylor多项式乘紧支函数实现这些jet；这不是把全空间非正规多项式当作物理参考。

二阶jet共 $6\times6=36$ 维，规范三阶jet共 $3\times10=30$ 维。定义

$$
C_{ij}(h)=2R^{(1)}_{ij}(h)
=\partial_l\partial_i h_{lj}+\partial_l\partial_j h_{li}
-\Delta h_{ij}-\partial_i\partial_jh_{ll}.
\tag{4}
$$

交换偏导直接给 $C(\delta_\xi h)=0$。三维Riemann张量由Ricci及其迹恢复，故 $C=0$ 时二阶曲率为零。令 $\Gamma_{i,jk}$ 为式(3)对应的场连接；零曲率使 $\partial_l\Gamma_{i,jk}$ 对 $j,k,l$ 全对称，它正可作为 $\partial_j\partial_k\partial_l\xi_i$，将该二阶jet消去。因此二阶不变量空间维数为6，由 $C_{ij}$ 承载。

代码独立列出全部基底，精确有理消元给规范映射秩 $6,18,30$、$C$ 的秩6，并检查 $C$ 乘规范映射精确为零及式(3)的矩阵逆。数值秩阈值未参与证明。于是式(2)的每个 $L_\alpha$ 均因子化通过 $C$；零／一阶项不能与二阶规范jet交叉抵消。

## 3. 低阶类的有限两点障碍

在TT扇区，式(4)为 $C=-\Delta h$。SO(3)不变常系数二次型在对称张量的迹与无迹部分分别为标量；TT动能及势遂为

$$
H_{\mathrm{TT}}(\mathbf k)=\frac12\left(A_T|p|^2+B_Tk^4|h|^2\right),
\quad A_T,B_T>0,\qquad
\omega(k)=b k^2,\quad b=\sqrt{A_TB_T}>0.
\tag{5}
$$

若允许TT零支，$\omega=0$ 也不能近似一个给定正 $ck$；它不属于这里的非退化比值公式。对任何预先共同标定的 $c,k>0$，令 $u=bk/c$，两点最大相对误差为

$$
E(u)=\max\{|u-1|,\ |2u-1|\},\qquad
\inf_{u>0}E(u)=\frac13,\qquad u_{\mathrm{opt}}=\frac23.
\tag{6}
$$

证明：$E<1/3$ 将同时要求 $u>2/3$ 和 $u<2/3$；等号在 $u=2/3$ 达成。这是两个固定非零波数的类内限制，不需UV或 $k\to0$。它是必要传播标准的数学拒绝，不是原物质端口已经实施的频率仪器证书。

仅要求积分规范不变会改变限制。线性Einstein张量 $G^{(1)}_{ij}=R^{(1)}_{ij}-\delta_{ij}R^{(1)}/2$ 满足 $\partial_iG^{(1)}_{ij}=0$，故

$$
\int h_{ij}G^{(1)}_{ij}(h)\,d^3x
\tag{7}
$$

在周期边界或妥当边界收支下规范不变，TT中等于 $\frac12\int|\nabla h|^2$，允许 $k^2$ 势。然而它**不是全未约束相空间的正成本**：$h_{ij}=\phi\delta_{ij}$ 给式(7)为 $-\int|\nabla\phi|^2$。处理标量、约束及边界是额外工作，不能从TT正性领取完整健康GR。

## 4. 跨阶混合的严格正两点逃逸

仍用同一 $h,p,C$ 和正超局域动能，取

$$
H_{\mathrm{esc}}=\frac12\int d^3x\left[
|p|^2+\frac97\left(\left|\left(1+\frac{\Delta}{6}\right)C\right|^2
+\frac1{12}|C|^2\right)\right].
\tag{8}
$$

每个成本仍是局域、严格点态规范不变的正平方。TT频率为

$$
\omega^2(k)=\frac97 k^4\left[\left(1-\frac{k^2}{6}\right)^2+\frac1{12}\right].
\tag{9}
$$

$k=1,2$ 处精确有 $\omega^2=1,4$，即与同一 $c=1$ 的线性目标完全一致。正基底 $1/12$ 保证全部 $k\ne0$ 的TT频率严格正，不靠额外零点。该反例不在 $\mathcal L_2$ 内：一个平方含二阶和四阶 $h$ 导数的混合。它证明“正平方”本身不足以给所有有限两点的公共障碍；不能把展开后的负交叉系数误判为负成本。

## 5. 整个有限频带的显式局域近似

不是只配两点。固定 $s=k^2\in[1,4]$，令

$$
c_j=\frac{\binom{2j}{j}}{4^j},\qquad
S_n(s)=\sum_{j=0}^{n}c_j\left(1-\frac{2s}{5}\right)^j,\qquad
P_n(s)=\sqrt{\frac25}\,S_n(s).
\tag{10}
$$

二项式展开 $(1-x)^{-1/2}=\sum c_jx^j$，且 $0<c_j\le1$。因 $|x|=|1-2s/5|\le3/5$，整个带上有

$$
|kP_n(k^2)-1|
\le\sqrt{\frac85}\frac{(3/5)^{n+1}}{1-3/5}
<\frac{16}{5}\left(\frac35\right)^{n+1}
 =:\varepsilon_{0,n}.
\tag{11}
$$

固定 $n=32$，再取 $\delta=10^{-12}$，势为

$$
V_{32}=\frac12\int d^3x\left[
\frac25\left|S_{32}(-\Delta)C\right|^2+\delta|C|^2\right],
\qquad
\omega_{32}^2(k)=k^4\left[\frac25S_{32}(k^2)^2+\delta\right].
\tag{12}
$$

它在全部非零TT波数严格稳定。平方内最高为66阶空间导数，Euler算子最高132阶；时间仍为普通二阶canonical振子，不能把高阶空间导数混称高阶时间不稳定。也不能由此声称严格因果锥、LR界或局域格点实现。

令 $r=kP_{32}(k^2)\ge1-\varepsilon_{0,32}$。加入正基底后相对频率为 $\sqrt{r^2+s\delta}$，增量至多 $2\delta/(1-\varepsilon_{0,32})<3\delta$。于是全带

$$
\left|\frac{\omega_{32}(k)}k-1\right|
\le\varepsilon:=\frac{16}{5}\left(\frac35\right)^{33}+3\cdot10^{-12}
<1.529\cdot10^{-7}.
\tag{13}
$$

对任意固定 $0<k_{\min}<k_{\max}<\infty$，把 $5/2$ 换成 $\mu=(k_{\min}^2+k_{\max}^2)/2$，收敛比为 $(k_{\max}^2-k_{\min}^2)/(k_{\max}^2+k_{\min}^2)<1$；同一尾界给任意所需精度的有限阶多项式，正基底可相应选小。因此没有全体有限阶正成本在所有固定非零带上的统一拒绝精度。

这不提供一个固定有限阶模型在含 $k=0$ 的整个红外域中的线性传播。对固定系数和正基底，$k\to0$ 仍有 $\omega\sim b_0 k^2$。所选系数依赖声明的有效域，是建模输入；有限域近似不是选择这些系数的认知机制。

## 6. 同一canonical对象的完整量子后态

为使全带结论超出谱点拟合，采用一份公开有限调节器：边长 $2\pi$ 的三维周期盒，所有整数波矢 $0<|\mathbf k|\le2$，每个实场的 $\pm\mathbf k$ 合成余弦／正弦坐标，每方向有两份TT极化。各壳整数向量数为

$$
\begin{array}{c|rrrr}
|\mathbf k|^2&1&2&3&4\\ \hline
\#\mathbf k&6&12&8&6
\end{array}
\qquad M=64\text{ 个实canonical振子}.
\tag{14}
$$

两个Hamiltonian都作用于同一 $\mathcal H=\bigotimes_{a=1}^{64}L^2(\mathbb R,dQ_a)$：

$$
H_0=\frac12\sum_a(P_a^2+\omega_{0a}^2Q_a^2),\qquad
H_{32}=\frac12\sum_a(P_a^2+\omega_{32,a}^2Q_a^2),
\quad \omega_{0a}=|\mathbf k_a|.
\tag{15}
$$

Schwartz空间给共同核，两个正二次振子Hamiltonian均有自伴实现。所有 $Q_a,P_a$ 身份固定。只用 $H_0$ 定义一次参考升降算符 $a_a$ 及输入投影 $\Pi_N=\bigotimes_a1_{0\le n_a\le N}$，固定 $N=1$；不把 $H_{32}$ 的另一个粒子数识别成 $n_a$。任意未知输入和其任意参考R可占据这个 $2^{64}$ 维初码。

在单模参考码上，

$$
Q^2=\frac{a^2+a^{\dagger2}+2a^\dagger a+1}{2\omega_0},
\qquad
\|Q^2\Pi_N\|
\le\frac{\sqrt{N(N-1)}+\sqrt{(N+1)(N+2)}+2N+1}{2\omega_0}
\le\frac{2N+3}{\omega_0}.
\tag{16}
$$

$N=0$ 时第一个根式按零理解。本式允许输出 $N+2$ 占据，未将 $Q^2$ 压回初码。真正共同差是

$$
H_{32}-H_0=\frac12\sum_a(\omega_{32,a}^2-\omega_{0a}^2)Q_a^2.
\tag{17}
$$

它通常会产生挤压／成对激发；不能用 $\sum_a(\omega_{32,a}-\omega_{0a})n_a$ 替换。

由于 $e^{-itH_0}$ 保 $\Pi_N$，prepared-Duhamel恒等式在该有限码上给

$$
\begin{aligned}
\|(e^{-itH_{32}}-e^{-itH_0})\Pi_N\|
&\le |t|\left(\varepsilon+\frac{\varepsilon^2}{2}\right)
\sum_a\omega_{0a}(2N+3)\\
&\le640\left(\varepsilon+\frac{\varepsilon^2}{2}\right)
=:\eta<10^{-4},\qquad |t|\le1 .
\end{aligned}
\tag{18}
$$

最后仅用64模、$\omega_{0a}\le2$、$N=1$，未优化壳权重。向上展示的实际界为 $9.779795519123902\times10^{-5}$；精确有理值见结果文件。初码参考能源至多 $\frac32\sum_a\omega_{0a}\le192$，候选初始能源至多 $192(1+\varepsilon)^2$，并由其自主H守恒。零点能保留，没有通过各自真空减能来混同身份。

张量 $I_R$ 不改变式(18)范数；对码内任意纯化，输出向量距离不超过 $\eta$，故两完整输出的迹距离 $D=\|\rho-\sigma\|_1/2\le\eta$。混态由纯化／收缩继承。再接同一任意完整CP仪器、保其所有分支和细记录，距离仍不增；这是限制初码的半diamond界。该结论并未将 $H_{32}$ 输出投影回 $\Pi_N$，保留64模的**全部Fock占据**和R。

调节器内所有非零模都在声明带内，没有带外真空被免费删除。若附加零模或旁观参考，双方使用同一自由演化因子，结论张量继承。无限空间／无限UV、其它 $h$ 扇区、物质相互作用及其后续不是本次量词；“完整Fock”不指已解决这些未纳入对象。Fourier截止也不等于严格局域的实际空间准备；有限微分表达式的局域性与有限调节器的性质分开登记。连续表达式可SO(3)不变，立方周期盒及有限整数波矢集只保其相应离散旋转，不能把这份调节器写成全SO(3)实际模型。

## 7. 三组检查、成熟来源和决定

[代码](check.py)仅用标准库精确分数、整数平方根和有向展示值；没有采样拟合承担全带证明，没有截断演化Fock矩阵。默认只读重算，[结果](results.json)只由首次显式写入创建。三组为：jet及低阶两点限制；严格正两点逃逸；全带余项、64模计数及完整未知输入误差。源码、文本绑定和精确有理字段为权威；浮点仅有向展示。它是数学算术证书，不是硬件浮点形式验证或真实装置精度。

[Gu—Wen，0907.1203，§IV.2—IV.4、VI—VII](https://arxiv.org/html/0907.1203v3)已讨论曲率平方的二次色散、附加约束后的三次色散以及仅积分规范不变路线；其可靠L-type为三次色散，N-type线性色散的晶格推断被原文标为未受控。本轮二次色散属于式(2)自定的超局域动能类，不是对其L-type的直接复算。我们不据此领取健康量子GR，也不沿用“只有引力模才许可”的额外排除。jet、二项式逼近及Duhamel均是成熟方法；本轮增量是此候选条件的精确删弱范围及同canonical后态的有限证书。

因此：**拒绝把资源完整性自动加强为低阶逐点正不变成本；同时拒绝把该低阶障碍扩大为全部有限有效域的不可能性。** 高阶局域拟合没有选择作用，也未接原共同钟／来源，不能作为生成完成。区域边界、规范表示、动量导数、辅助场、非局域消元、额外物质／关系模式都是另有条件的逃逸，未被本轮排除；暗部门更没有被排除或确认为暗物质、暗能量。下一步不再为增加多项式阶数另立轮次。
