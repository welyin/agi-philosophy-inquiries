# 1046 解析证明：有限波包与同一偶极接收器

2026-10-08。所有精确等式针对下述有效作用；不声称已界定它与完整 QED、真实某种原子的误差。

## 1. 物理采用及正规共同对象

采用三维 Euclidean 空间、$\hbar=c=1$，并取能量单位 $\mu$，下述频率与时间分别以 $\mu$、$\mu^{-1}$ 计。原子钉扎于原点；基态和激发态是相反宇称、分别变换为 $j=1/2$ 的二重态。能隙 $\Omega=3/2$，偶极的约化矩阵元并入 $g=1/10$。从矢量协变关系

$$
U_R D_j U_R^\dagger=\sum_iR_{ij}D_i
$$

在两份 $j=1/2$ 间展开一般 $2\times2$ 矩阵，标量项不能作矢量，三个 Pauli 分量的交织唯一至一个共同系数。因此选择相位后 $D_i=\sigma_i$ 是该已给多重态的偶极形式；原子的物种、对称性与二重态选择并未由认知推出。

横向单光子空间与测度为

$$
\mathcal F=\{\psi\in L^2(\mathcal B,d^3k;\mathbb C^3): k\cdot\psi(k)=0\},\quad
\mathcal B=\{1<|k|<2\},\quad \omega(k)=|k|,
$$

$$
P_k=I-\hat k\hat k^T,\qquad
v(k)=\sqrt{\frac2{3\pi}}\frac{\sin^2[\pi(|k|-1)]}{|k|}\quad(k\in\mathcal B). \tag{P1}
$$

域外延为零。$\int|v|^2d^3k=1$，径向概率为 $r(\omega)d\omega=(8/3)\sin^4[\pi(\omega-1)]d\omega$，其均值为 $3/2$。有限频带和径向形状是所选有效耦合的输入，不从该函数推出微观最小长度或严格局域场。

单激发扇区包括“基态原子加一个光子”及“激发原子加真空”：

$$
\mathcal H_1=(\mathbb C_g^2\otimes\mathcal F)\oplus\mathbb C_e^2,
\quad H_0=(I_2\otimes\omega)\oplus\Omega I_2,
$$

$$
A\Big(\sum_a|a\rangle\otimes\psi_a\Big)
=\sum_{a,i}\sigma_i|a\rangle\int v(k)\psi_{a,i}(k)d^3k,
\qquad H=H_0+g\begin{pmatrix}0&A^\dagger\\A&0\end{pmatrix}. \tag{P2}
$$

这是固定的 Schrödinger 图景旋波偶极作用，不随入射方向改变。$A^\dagger$ 包含 $P_k$ 横向投影。因 $P_k$ 为实对称、$\operatorname{tr}P_k=2$，Pauli 反对易给

$$
AA^\dagger=\int |v|^2\sum_{ij}(P_k)_{ij}\sigma_i\sigma_jd^3k=2I_2,
\quad\|A\|=\sqrt2. \tag{P3}
$$

所以 $H$ 是整个 $\mathcal H_1$ 上的有界自伴算符，且

$$
\frac67I<H<\frac{15}7I. \tag{P4}
$$

没有把有限离散模式的数值矩阵当作连续场自伴性证明。有限频带是在本父模型中已采用的正有效扇区，不认证对任意域外脉冲仍成立。RWA、频带形状、钉扎及单激发准备是独立输入；没有模拟机械反冲或无截断 QED。

## 2. 正常方向包与同一准备资源

给 $n\in S^2$，选实正交 $u,v$ 满足 $u\times v=n$，令 $\epsilon_+(n)=(u+iv)/\sqrt2$。改变这份局部框架只乘整体相位；以下光子态和效果与之无关。定义

$$
w(c)=\begin{cases}(2c-1)^2,&c>1/2,\\0,&c\le1/2,\end{cases}
\qquad f_n(k)=\frac{v(k)w(\hat k\cdot n)P_k\epsilon_+(n)}{\sqrt{a_2}}, \tag{P5}
$$

其中

$$
a_j=\frac14\int_{1/2}^1(1+c^2)w(c)^jdc,
\quad a_1=\frac{71}{960},\quad a_2=\frac{31}{672}. \tag{P6}
$$

故 $\|f_n\|=1$。这是有限角锥、有限频带的普通 $L^2$ 波包，非方向或动量本征态。角窗及径向窗在支撑边界的值和一阶导数均为零；延拓向量场具有平方可积的一阶导数，因而其普通 Fourier 位置分布有有限二阶矩。没有紧位置支撑或精确事件局域化的声明。

轴对称与角支撑给

$$
\langle k\rangle_{f_n}=\frac32\zeta n,\qquad \frac12<\zeta<1. \tag{P7}
$$

因此 $n$ 是波包的真实平均传播动量方向，不是给内部 Bloch 球重新命名。所有方向具有同一径向能量分布、同一角宽、同一准备类型和一致的全 $H$ 能量上界。源与未知原子／被动参考初始独立；各方向准备的可用性与接收器支撑是输入，不宣称可免费从任意旧光子制造这些包。

由于 $U_R$ 作用在两原子二重态，场作用为 $(\mathcal U_R\psi)(k)=R\psi(R^{-1}k)$，式(P2)严格旋转不变。$f_{Rn}=\mathcal U_Rf_n$ 至多差整体相位。各方向使用相同 $H$、相同等待及同一个激发投影 $Q_e$。

## 3. 完整演化的秩一性，不依赖最低阶截断

记 $R(t)=\int |v(k)|^2e^{-i\omega(k)t}d^3k$。式(P3)逐频率仍成立，给

$$
A(I_2\otimes e^{-i\omega t})A^\dagger=2R(t)I_2. \tag{P8}
$$

同时角积分给

$$
A(\xi\otimes e^{-i\omega t}f_n)
=\frac{a_1}{\sqrt{a_2}}R(t)\,\sigma\cdot\epsilon_+(n)\xi,
\qquad \sigma\cdot\epsilon_+(n)=\sqrt2|+n\rangle\langle-n|. \tag{P9}
$$

将完整 Schrödinger 方程的基态—光子分量消去，激发分量 $e(t)$ 满足

$$
\dot e(t)=-i\Omega e(t)-ig\frac{a_1}{\sqrt{a_2}}R(t)\sigma\cdot\epsilon_+\xi
-2g^2\int_0^t R(t-s)e(s)ds,\quad e(0)=0. \tag{P10}
$$

这是有连续有界核的精确 Volterra 方程。其解的唯一性及标量自能核使

$$
Q_e e^{-iHt}J_n=c(t)|+n\rangle\langle-n|,
\quad J_n\xi=\xi\otimes f_n,\quad
E_n(t):=J_n^\dagger e^{iHt}Q_e e^{-iHt}J_n=p(t)P_{-n},
\quad p(t)=|c(t)|^2. \tag{P11}
$$

所有方向的 $c(t)$ 相同。这也直接说明暗自旋的激发振幅对所有时间严格为零，而非取小 $g$ 后把高阶暗吸收丢掉。自旋翻转来自原偶极矩阵，与光子偏振本身是否有自主 qubit 因子不同。

## 4. 一个有限时间的有理正对比证书

取共同 $t=1$。相互作用图景中 $\|gV_I(s)\|=g\sqrt2<1/7$；基态到激发态只含奇数阶。第一阶亮态振幅满足

$$
|c_1(1)|\ge g\alpha\left(1-\frac18\right),
\qquad \alpha=\frac{\sqrt2a_1}{\sqrt{a_2}}>\frac{12}{25}. \tag{P12}
$$

证明：在转去 $e^{-i\Omega t}$ 后，第一阶径向积分的实部为 $\int_0^1\int r(\omega)\cos[(\omega-\Omega)s]d\omega ds$；$|\omega-\Omega|\le1/2$，故逐点 $\cos\ge1-1/8$。$\alpha$ 下界由式(P6)的有理数平方直接核实。

设 $x=g\sqrt2<1/7$，全部余下奇数阶满足

$$
|c-c_1|\le\sinh x-x
\le\frac{x^3}{6(1-x^2/20)}
<\frac{980}{2014782}<\frac1{2000}. \tag{P13}
$$

后续奇数阶与上一项的比值至多 $x^2/20$。因此

$$
|c(1)|>\frac{21}{500}-\frac1{2000}=\frac{83}{2000},
\quad p(1)>\frac{6889}{4000000}=.00172225. \tag{P14}
$$

这是真正整个连续模型的保守下界，不是数值积分差或有限模式收敛拟合。对任何 $n$，$\|E_n-E_{-n}\|=p(1)$，同一输入 $|-n\rangle$ 见证该反向概率差。正 $p$ 还使方向差空间张成全部 $\operatorname{Herm}_0(2)$，且 $E_{Rn}=U_RE_nU_R^\dagger$；变化率为 $\|E_n-E_m\|=p|n-m|/2$。这项上界只用于有限误差，不把384已消去的下界 Lipschitz 当回必要条件。

若实际来源、传播与末读共同给每份方向效果算子误差至多 $\varepsilon$，则反向 gap 至少 $6889/4000000-2\varepsilon$。误差预算必须另证；本轮不假造仪器实测精度，也不对无限多重复任务许诺无资源成本。

## 5. 未知输入、参考和实际后态

对任意未知 $\rho_{QR}$，$R$ 被动，真实两个未归一化输出为

$$
\mathcal I_{b,n}(\rho_{QR})=(K_{b,n}\otimes I_R)\rho_{QR}(K_{b,n}^\dagger\otimes I_R),
\quad K_{e,n}=Q_e e^{-iHt}J_n,\quad
K_{g,n}=(I-Q_e)e^{-iHt}J_n. \tag{P15}
$$

两者的输出空间保留完整原子／场扇区，$\sum_bK_b^\dagger K_b=I_2$。若引入记录寄存器，则取 $\sum_b|b\rangle\langle b|\otimes\mathcal I_b$。没有将完整仪器换成只作用于旧 qubit 的 $\sqrt{E_n}$：吸收支路原子位于激发二重态且自旋从 $-n$ 变为 $+n$，未吸收支路仍可有原子—光子关联。

最终 $Q_e$ 是采用的固定 PVM 权限。$[H,Q_e]\ne0$，其算子范数为 $g\sqrt2$。非选择测量 $\mathcal D(\rho)=Q_e\rho Q_e+(I-Q_e)\rho(I-Q_e)$ 满足

$$
\mathcal D^*(H)=H_0,
\quad\Delta\langle H\rangle=-g\langle V\rangle,
\quad|\Delta\langle H\rangle|\le g\sqrt2. \tag{P16}
$$

故自主的相互作用阶段不等于末端测量装置能源闭合。上界是整个单激发扇区的预算；本轮特定对称准备在某些时刻可恰好零平均能变，仍不改变非对易事实。测后正能量及平方矩有共同有界范围，由(P4)保持，不声称读出器免费执行。

## 6. 真正无偏振的正常反例

不能将有限角宽的 $P_k\epsilon_\pm(n)$ 错称每个 $k$ 的纯 helicity 态，也不能将两包等混合自动称为逐动量无偏振。下面给一份满足后者的有限秩态。

令实验室实基 $e_j$，$h_{j,n}(k)=v(k)w(\hat k\cdot n)P_ke_j$（未归一化），并定义

$$
j_2=\frac12\int_{1/2}^1w(c)^2dc=\frac1{20},
\qquad \rho_n^{\rm unpol}=\frac1{2j_2}\sum_{j=1}^3|h_{j,n}\rangle\langle h_{j,n}|. \tag{P17}
$$

它是正常、迹1、秩3的光子态，点态偏振核为

$$
\rho_n^{\rm unpol}(k,k)=\frac{|v(k)|^2w(\hat k\cdot n)^2}{2j_2}P_k. \tag{P18}
$$

因此每个所含动量上两个横向偏振等权；不是不可归一化的精确动量经典混合。仍只有正向角锥支持、$\langle k\rangle$ 为正倍数 $n$，径向能量分布与前例相同。

设

$$
B_n=a_1(I-nn^T)+b_1nn^T,\qquad b_1=3/160.
$$

角积分使 $Ae^{-i\omega t}(\xi\otimes h_{j,n})=R(t)\sigma\cdot(B_ne_j)\xi$。式(P10)的同一标量自能表明精确激发映射与 $\sigma\cdot(B_ne_j)$ 成共同系数。各向量 $B_ne_j$ 为实，故

$$
\sum_j[\sigma\cdot(B_ne_j)]^\dagger[\sigma\cdot(B_ne_j)]
=\operatorname{tr}(B_n^2)I_2. \tag{P19}
$$

无偏振包在同一固定仪器下的效果为 $q(t)I_2$，与 $n$ 无关。该光仍有明确平均传播方向，接收器却不再给反向对比。这是本仪器的失败判据，不是所有探测器不能测无偏振光的方向。

完整后态亦可保留三个准备分量的 Kraus 映射，且不把隐藏混合标签赠给读取者。控制反例还可用两相反轴向圆包的等混合，所得效果为 $p(t)I/2$；它与(P17)的真正逐动量无偏振态需区别。

## 7. 数值有限表示与不承担的证明

径向星形算符是连续场的精确角分解：令 $T$ 为角积分的偶极映射，$TT^\dagger=2I$，则 $B=T^\dagger/\sqrt2$ 给两份正交亮角模式。每份原准备 $J_n$ 在它们上的投影是

$$
C_n=B^\dagger J_n=\frac{a_1}{\sqrt{a_2}}|+n\rangle\langle-n|.
$$

其余部分为暗模式，Gram 矩阵 $I-C_n^\dagger C_n$。有限径向求积后，数值可在“激发2维＋每径向两亮两暗”中精确保留该已声明准备的全部原子／场后态。不同方向所用暗基可以不同；不把各自压缩擅自当成保任意方向**相干叠加**的一份共同小编码。共同原父对象始终是(P2)。

有限求积、完整矩阵指数及独立角积分只校准实现。原连续模型的自伴性、全时间 rank1、全部方向和有限误差均由上述证明承担，网格加密差不是完整 QED 的物理误差证书。

## 8. 历史接口及停止线

所得 $E_n$ 给旧382—384的固定载体、固定迹、连续反向可分、全族旋转协变和对比完整的**方向仪器部分**。数学旋转协变不自动实现位置端口的真实可逆重定向；全族入射准备及其参考资源仍是采用权限。它不提供真实位置端口全部邻域，也未交付386强缩放或425一致半幅／成本合同。我们已采用三维 Maxwell 空间，因此不得反过来引用旧维数上界宣传成认知生成三维。

本轮也没有给1041的事件同一性、任意可替换中继或闭环接线权限；带限耦合不自动是严格时空局域信号。准备、支撑、时间标定和末读仍需其任务合同。到这份实际方向桥与偏振删除反例为止，停止器件优化。

成熟背景：[Wang等，1010.4661，§II、式(1)—(5)、§III.1](https://arxiv.org/html/1010.4661)给固定位置原子、横向传播波包及偶极RWA接口；[Stobińska等，0808.1666](https://arxiv.org/html/0808.1666)给自由空间单光子激发背景。本轮没有采用其最优吸收结论、Markov／Weisskopf–Wigner近似或实验数值作为自己的证书。此处双二重态、有限频带作用及全参考仪器由(P1)—(P19)自行定义和证明。
