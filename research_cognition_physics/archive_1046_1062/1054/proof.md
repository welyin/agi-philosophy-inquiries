# 1054证明：同一偶极作用的位置概率坐标

2026-10-08。科学基线1053，父物理作用为[1046](../1046/proof.md)。下文是已声明有限频带、单激发、旋波偶极模型内的精确结论；不是该模型相对完整QED的误差定理。特别是 $t_0=.01\mu^{-1}$ 短于载频周期，不能据模型内部的严格余项声称忽略反旋转项已获现实物理认证。

## 1. 采用对象、位置身份和权限

沿1046采用3D Euclidean空间、$\hbar=c=1$，频率单位 $\mu$、长度单位 $\mu^{-1}$。接收原子的两份二重态为基态 $\mathbb C_g^2$ 和激发态 $\mathbb C_e^2$，$\Omega=3/2$、$g=1/10$。正常单光子空间、横向投影和耦合函数为

$$
\mathcal F=\{\psi\in L^2(\mathcal B,d^3k;\mathbb C^3):k\cdot\psi(k)=0\},\quad
\mathcal B=\{1<|k|<2\},\quad P_k=I-\hat k\hat k^T,
$$


$$
v(k)=\sqrt{\frac2{3\pi}}\frac{\sin^2[\pi(|k|-1)]}{|k|},\qquad
\int |v|^2d^3k=1,\qquad \omega=|k|. \tag{1}
$$

在 $\mathcal H_1=(\mathbb C_g^2\otimes\mathcal F)\oplus\mathbb C_e^2$ 上，令接收器的钉扎位置为经典参数 $x$，同一偶极有效作用在该点的表达是

$$
A_x(\xi\otimes\psi)=\int v(k)e^{ik\cdot x}\sigma\cdot\psi(k)\xi\,d^3k,
\qquad H_x=(I_2\otimes\omega)\oplus\Omega I_2+
g\begin{pmatrix}0&A_x^\dagger\\A_x&0\end{pmatrix}. \tag{2}
$$

若 $T_a\psi(k)=e^{-ik\cdot a}\psi(k)$，$V_x=(I_2\otimes T_x)\oplus I_2$，则 $H_x=V_xH_{x=0}V_x^\dagger$。1046的 $A_xA_x^\dagger=2I$ 与谱界直接运输：

$$
\frac67I<H_x<\frac{15}7I. \tag{3}
$$

这里 $x$ 不是接收器的量子位置算符，不包含机械反冲；该酉等价不是可主动执行的位移操作。不把多种位置下的同一作用表达叫作“一个自主装置能自由移动原子”。位置参数的3D身份及其准备参照已采用，以下新增的是其实际读数字典。

定义三个正常横向向量与混合光源

$$
h_{j,a}(k)=e^{-ik\cdot a}v(k)P_ke_j,\qquad
\rho_a=\frac12\sum_{j=1}^3|h_{j,a}\rangle\langle h_{j,a}|.
\tag{4}
$$

角平均 $\int P_kd\Omega/(4\pi)=2I/3$ 给 $\langle h_{i,a},h_{j,a}\rangle=2\delta_{ij}/3$。故 $\rho_a$ 为正常秩3态，非零谱均为 $1/3$。其自由光子能源分布为 $(8/3)\sin^4[\pi(k-1)]dk$，均值 $3/2$，与 $a$ 无关；初始全 $H_x$ 平均亦为 $3/2$，因为初始位于基态—光子部门，作用交叉项期望为零。没有声称所有全 $H_x$ 高阶矩均与 $x-a$ 无关；共同有界谱已给统一资源上界。

径向窗在频带端点一阶消失，角投影在 $|k|>1$ 上光滑，所以这些包具有有限普通Fourier位置二阶矩。它们无紧位置支撑，不是方向或位置本征态。三份准备参照取

$$
a_i=r_0e_i,\qquad r_0=\frac14,\quad i=1,2,3;\qquad t_0=\frac1{100}. \tag{5}
$$

源与任意未知基态自旋 $Q$、其被动参考 $R$ 初始独立。不同 $a_i$ 是事先声明的正常光包准备权限，未构造三台发射器。锚点框架、钉扎、计时、基态部门准备、新鲜光包和原激发末读 $Q_e$ 均采用。三个读数是三个备选实验的概率，不是把同一未知输入复制三份，也不是免费将吸收后的原子复位。重复估计还须相应稳定性与统计资源。

## 2. 完整演化的标量概率

允许保留光源的纯化标签 $E\simeq\mathbb C^3$，共同等距为

$$
J_a\xi=\frac1{\sqrt2}\sum_j(\xi\otimes h_{j,a})\otimes|j\rangle_E. \tag{6}
$$

对 $r=x-a$，消去基态光场时自能和注入分别为

$$
A_xe^{-i\omega s}A_x^\dagger=2R(s)I_2,\qquad
R(s)=\int |v|^2e^{-i\omega s}d^3k,
$$


$$
B_s(r)=\int |v(k)|^2e^{-i\omega s}e^{ik\cdot r}P_k\,d^3k. \tag{7}
$$

奇角部分由 $P_{-k}=P_k$ 消去；余下为实对称矩阵角核乘复的径向系数。绕 $r$ 的旋转遂给

$$
B_s(r)=b_T(s,|r|)P_T+b_L(s,|r|)P_L,\quad
P_L=\hat r\hat r^T,\quad P_T=I-P_L. \tag{8}
$$

这是具体横向矩阵的结论。一般旋转协变的qubit效果仍可含 $\hat r\cdot\sigma$，不能仅从旋转宣布下文的标量性。

完整激发解由共同标量Volterra方程唯一决定，故第 $j$ 个未除 $\sqrt2$ 的注入给振幅

$$
Q_ee^{-itH_x}(\xi\otimes h_{j,a})=\sigma\cdot[C_t(r)e_j]\xi,
\qquad C_t(r)=c_T(t,|r|)P_T+c_L(t,|r|)P_L. \tag{9}
$$

所有高阶吸收和重发射均包含在 $C_t$ 中。将基转至沿 $r$ 的实正交基，再用 $\sigma_i^2=I$，得到

$$
E_{a,x}(t)=J_a^\dagger(e^{itH_x}Q_ee^{-itH_x}\otimes I_E)J_a
=p_t(|x-a|)I_2,
\qquad p_t(r)=\frac12(2|c_T|^2+|c_L|^2). \tag{10}
$$

所以该概率对任意未知自旋及其旧参考相同；并非假定它们已经重准备成某个校准态。所有证明在连续光子空间成立。

## 3. 时间偶性与最低非零空间斜率

令 $B(r)=B_0(r)$。沿 $r=e_3r$，其本征值为 $B_T(r),B_T(r),B_L(r)$，均实，定义

$$
F(r)=\frac12\operatorname{tr}B(r)^2=\frac12(2B_T^2+B_L^2). \tag{11}
$$

需要证明完整 $p_t$ 的三阶系数确实为零。对于每个整数 $m\ge0$，$A_x\omega^mA_x^\dagger=2\langle\omega^m\rangle I$ 是实标量；$A_x\omega^m(\,\cdot\,\otimes h_{j,a})$ 为Pauli乘实对称角矩。展开 $Q_eH_x^n$ 的块乘积，所有完整往返贡献这些实标量，剩下一个这样的注入。因此 $C_t=\sum_{n\ge1}(-it)^n C_n/n!$，各 $C_n$ 为实对称、具有式(8)的轴结构。$H_x$ 有界使该级数在任意有限 $t$ 收敛，故

$$
C_{-t}=C_t^*,\qquad p_{-t}(r)=p_t(r),\qquad
p_t(r)=g^2t^2F(r)+O(t^4). \tag{12}
$$

不能只根据最低阶吸收概率删去奇数阶；式(12)提供全阶依据。

取 $r_0=1/4$，$|k\hat k_3r_0|\le1/2$，用 $\cos z\ge1-z^2/2\ge7/8$ 和 $\sin z/z\ge1-z^2/6\ge23/24$。径向二阶矩 $\mu_2=\langle|k|^2\rangle\ge1$；角矩为

$$
\langle\hat k_3^2(1-\hat k_1^2)\rangle=4/15,\qquad
\langle\hat k_3^2(1-\hat k_3^2)\rangle=2/15.
$$

遂有

$$
B_T,B_L\ge7/12,\quad -B_T'\ge23r_0/90,\quad -B_L'\ge23r_0/180,
$$


$$
F(r_0)\ge49/96,\qquad
-F'(r_0)=2B_T(-B_T')+B_L(-B_L')\ge\frac{161}{1728}. \tag{13}
$$

## 4. 全部阶数的有限有理证书

有限频带给方向导数的统一算子界（任意单位方向均成立）

$$
\|H_x\|<M=15/7,\qquad\|\partial_iH_x\|<2/7,\qquad
\|\partial_i\partial_jH_x\|<4/7. \tag{14}
$$

例如 $\partial_iA_x$ 在 $A_x$ 的核中插入 $ik_i$，其范数至多 $2\sqrt2$；乘 $g$ 小于 $2/7$。高阶导数同理存在，有界的整个表达对 $x$ 光滑。

在 $e^{itH_x}Q_ee^{-itH_x}$ 的嵌套对易子展开中，第 $n$ 项位置导数有 $n$ 个插入位置，故其范数至多 $n2^nM^{n-1}(2/7)$。结合式(12)，$z=2Mt<1$ 时

$$
|\partial_rp_t-g^2t^2F'|\le\frac2{15}\sum_{n\ge4}\frac{n z^n}{n!}
\le\frac2{15}\frac{z^4e^z}6\le\frac2{15}\frac{z^4}{6(1-z)},
$$


$$
|p_t-g^2t^2F|\le\frac{z^4}{24(1-z)}. \tag{15}
$$

只保偶数项会更紧，但式(15)宁可累加全部 $n\ge4$。没有把算符阶数截断当成完整动力学。

在固定 $t_0=1/100$，$z=3/70$，得到

$$
-p'_{t_0}(r_0)>\kappa:=\frac{161}{1728000000}-\frac9{114905000}
=\frac{589541}{39711168000000}>10^{-8}, \tag{16}
$$


$$
p_{t_0}(r_0)>\frac{49}{96000000}-\frac{27}{183848000}
=\frac{802069}{2206176000000}>\frac1{3000000}. \tag{17}
$$

严格性可沿式(3)、(14)的严格界保留；不影响以下保守非严格常数。固定时间有真实非零吸收，未让任务在 $t\to0$ 极限退成恒等。数值显示约 $6.36\times10^{-7}$，概率较小；不能隐去达到某个统计精度所需的重复资源。

## 5. 同一实际读口给局部概率坐标

定义三备选实验的概率图

$$
\mathcal P(x)=\big(p_{t_0}(|x-a_1|),p_{t_0}(|x-a_2|),p_{t_0}(|x-a_3|)\big).
\tag{18}
$$

在 $x=0$，$D\mathcal P(0)=-p'_{t_0}(r_0)I_3$。因此逆函数定理给局部开图；不是把预先存好的位置数字作为输出。

为明确一个共同窗口，不只援引存在性。对 $K_x=(Q_e\otimes I_E)(e^{-itH_x}\otimes I_E)J_a$，有

$$
\|K_x\|\le t/7,\quad \|\partial_iK_x\|\le2t/7,\quad
\|\partial_i\partial_jK_x\|\le4t/7+4t^2/49. \tag{19}
$$

第一界由初始在 $Q_g$ 和 $\|[H_x,Q_e]\|=g\sqrt2<1/7$ 得到；其余由有界Duhamel公式及二阶有序积分得到。令

$$
C=16t_0^2/49+8t_0^3/343=1401/42875000.
$$

$\partial_i\partial_j(K_x^\dagger K_x)$ 的四项乘积给每份概率Hessian各分量至多 $C$，其算子范数至多 $3C$；三个读口的Jacobian变化模小于 $6C|x-y|$。取

$$
\mathcal U=\{x:|x|<10^{-5}\},\qquad
6C10^{-5}=4203/2143750000000<2\times10^{-9}. \tag{20}
$$

对 $x,y\in\mathcal U$ 沿凸线段积分Jacobian，并将 $D\mathcal P(0)$ 分离，可得

$$
\|\mathcal P(x)-\mathcal P(y)\|_2\ge(\kappa-6C10^{-5})|x-y|
>\frac{|x-y|}{2\times10^8}. \tag{21}
$$

故 $\mathcal P$ 在该球上单射、局部微分同胚且逆图具有显式误差界。若两个位置都与同一估计概率向量逐项相容至 $\eta$，则

$$
|x-y|\le4\sqrt3\times10^8\eta. \tag{22}
$$

这不是断言任意有噪概率向量都落在像中，也不是数值网格认证邻域。解析Hessian、凸域与式(16)承担所有点的量词。若另给独立重复Bernoulli试验，每个设置 $N$ 次，则通常Hoeffding并合界给同时误差 $\eta$ 的失败概率至多 $6e^{-2N\eta^2}$；这种稳定重复、基态准备与资源并未由当前 $H_x$ 自行产生。

## 6. 完整后态、参考与读出资源

保源纯化时实际支路算子必须写成

$$
K_b=(Q_b\otimes I_E)(e^{-itH_x}\otimes I_E)J_a
=\sum_jK_{b,j}\otimes|j\rangle_E,
\quad K_{b,j}=\frac1{\sqrt2}Q_be^{-itH_x}(\,\cdot\,\otimes h_{j,a}). \tag{23}
$$

对未知 $\rho_{QR}$，完整cq输出是

$$
\sum_b|b\rangle\langle b|\otimes(K_b\otimes I_R)\rho_{QR}(K_b^\dagger\otimes I_R). \tag{24}
$$

这是 $j$ 的**相干和**。只有将源纯化标签偏迹，才得到 $\sum_jK_{b,j}\rho K_{b,j}^\dagger$；两者不可混称。式(10)仅给效果，不能用 $\sqrt{p_t}I_Q$ 替代真实后态。校准中输入 $|+z\rangle$ 的激发支路自旋Z期望约 $-.320524$，与这种伪替换的 $+1$ 明显不同。

所有输入及任意被动参考的保迹由 $\sum_bK_b^\dagger K_b=I_Q$ 精确保证。比较同一个源准备下的位置误差，Duhamel给

$$
\|e^{-itH_x}-e^{-itH_y}\|\le(2t/7)|x-y|.
$$

接固定等距、末读和偏迹不增加含参考误差，故实际仪器的半diamond距离亦至多 $(2t/7)|x-y|$。这是真实全部后态的连续性，未只保概率。

末读沿1046，$\|[H_x,Q_e]\|=g\sqrt2$。非选择末读删除交叉作用能，最大平均能源变化至多 $g\sqrt2$；检测器、支撑、源准备与计时未闭合。短时有限带模型仍有非紧局域尾，不提供严格光锥、无截断QED或完整实验误差。尤其本轮不认证 $t_0$ 下RWA反旋转项可忽略。

## 7. 旧方向任务在同一作用族中的位置

1046的正常圆偏振源族 $f_n$、$x=0,t=1$ 原效果直接复用：

$$
E_n^{(0)}=p(1)(I-n\cdot\sigma)/2,\quad
p(1)>6889/4000000.
$$

它与本轮使用同一 $H_x$ 作用族、同一个激发末读、不同的预声明源和等待菜单。对未经重定中心、但 $x\in\mathcal U$ 的方向任务，保守算子比较给

$$
\|E_n^{(x)}-E_n^{(0)}\|\le(4/7)|x|,\qquad
\|E_n^{(x)}-E_{-n}^{(x)}\|>6889/4000000-8\cdot10^{-5}/7
=47903/28000000>0. \tag{25}
$$

这只证明同一位置窗仍有统一反向对比；不把有限误差协变冒领384的精确完整群作用。要在任意其它位置精确得到原方向效果，需要另有共同平移准备／标定；从位置概率图不推出免费控制。更不将 $n$ 光包参数直接当成真实位置壳，二者的同一性、主动位移和完整邻域生成仍需对应合同。

## 8. 前提删除与结论强度

一个中心的概率恒为径向函数，等半径的不同位置精确给相同读数；它不能构成三维局部坐标。三个共线中心在共线基准位置的梯度均沿同一轴，Jacobian秩至多1。二者是当前读口类别的明确失败控制，不是否定所有其它光学定位方法。

本轮可签收的是：采用3D光场后，同一有界RWA作用提供实际方向效果和一个严格可反演的位置概率字典，且每份任务的真实后态与旧参考未删除。它没有从认知生成3D，没有提供425的位移群、半幅和成本，也没有提供386的真实缩放。有限Lie商1051是另一条件工具，本轮不声称当前无束缚场的自然能源具有紧预解。整体M3A与ROADMAP保持开放。
