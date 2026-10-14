# 固定Hartree候选的有限时间残差证书

2026-10-13。辅助准入证明与只读接口；不另编号，不修改已冻结的[复制误差证明](mean_field_bound.md)或[小K检查](mean_field_bound_check.py)。对象为[固定候选](mean_field_candidate/candidate.md)及其[preview](mean_field_candidate/preview.py)。本页处理数值轨道误差，不替代父模型、格点、实际仪器的其他预算。

## 1. 不以步长差充当误差

将每步保存的数值端点及导数视为精确的有限二进有理数，定义连续三次Hermite曲线。对它计算**整段**方程残差，再使用稳定性不等式。生成端点的RK4算法只负责给出候选曲线；它的截断及舍入误差均可反映在残差里，不要求先证明RK4的理论阶数。

在一步 $t=t_j+\Delta u,\ 0\le u\le1$，定义

$$
\widehat y(u)=c_0+c_1u+c_2u^2+c_3u^3,
\quad
\begin{aligned}
c_0&=y_j,&c_1&=\Delta f_j,\\
c_2&=3(y_{j+1}-y_j)-\Delta(2f_j+f_{j+1}),\\
c_3&=2(y_j-y_{j+1})+\Delta(f_j+f_{j+1}).
\end{aligned}
\tag{1}
$$

这些是保存数据所定义的**精确实数公式**；不能先把系数舍入后又假称所得幂基曲线严格连接原端点。辅助接口以区间基本运算包住式(1)的系数，保留这个区别。相邻段因此在其数学定义中连续，没有隐去网格接缝误差。

## 2. 消去快速旋转后的稳定性

写归一权重为 $\omega_s$，并令

$$
f_s=(\nu s,0,\epsilon m),\quad
d_s(g)=(s,0,m)+g_zf_s,\quad
b(r)=2\sum_s\omega_s f_s\cdot r_s,\quad
a(r)=(-\Omega,0,b(r)).
\tag{2}
$$

精确系统为
$\dot r_s=2d_s(g)\times r_s$、
$\dot g=2a(r)\times g$。
从指定单位初态出发，$|r_s(t)|=|g(t)|=1$。近似曲线残差定义为

$$
\rho_s=\dot{\widehat r}_s-2d_s(\widehat g)\times\widehat r_s,\qquad
\rho_g=\dot{\widehat g}-2a(\widehat r)\times\widehat g,\qquad
R^2=|\rho_g|^2+2\sum_s\omega_s|\rho_s|^2.
\tag{3}
$$

对误差 $e_g=g-\widehat g,\ e_s=r_s-\widehat r_s$，取

$$
E^2=|e_g|^2+2\sum_s\omega_s|e_s|^2,\qquad
F^2=\sum_s\omega_s|f_s|^2 .
\tag{4}
$$

用下列分解，两个快速旋转项与自身误差的内积都精确为零：

$$
\begin{aligned}
\dot e_g
&=2a(\widehat r)\times e_g
+4\left(\sum_s\omega_s f_s\cdot e_s\right)e_z\times g-\rho_g,\\
\dot e_s
&=2d_s(g)\times e_s+2(e_g)_z f_s\times\widehat r_s-\rho_s .
\end{aligned}
\tag{5}
$$

因此，若整段上 $|\widehat r_s|\le r_{\max}$，

$$
\dot E\le L E+R,\qquad
L=\sqrt2F(1+r_{\max}).
\tag{6}
$$

证明只用加权Cauchy及 $2xy\le x^2+y^2$；$E=0$ 时取右上Dini导数。这里没有把最大动量或最大单源旋转频率粗暴塞进指数。还可将第二行改写为近似轴旋转加真实 $r_s$ 的forcing，得到更紧的 $L=2\sqrt2F$；为不追逐小改进，当前接口保留合法的式(6)。

每段取常数上界 $L_j,R_j$，记
$\varphi_1(z)=(e^z-1)/z$、
$\varphi_2(z)=(e^z-1-z)/z^2$，在零点分别续为1和1/2，则

$$
\begin{aligned}
E_{j+1}&\le e^{L_j\Delta}E_j+
\Delta\varphi_1(L_j\Delta)R_j,\\
\int_{t_j}^{t_{j+1}}E(t)\,dt
&\le\Delta\varphi_1(L_j\Delta)E_j+
\Delta^2\varphi_2(L_j\Delta)R_j .
\end{aligned}
\tag{7}
$$

初始误差为零；若使用其他准备，须另填实际 $E_0$，不能自动置零。

## 3. Bernstein全段残差，不漏采样点之间

Bloch右端为二次多项式，三次Hermite代入后的残差为至多6次。设幂基残差为 $d(u)=\sum_{k=0}^6d_ku^k$，其Bernstein控制向量为

$$
b_j=\sum_{k=0}^j
\frac{\binom jk}{\binom6k}d_k,\qquad
d(u)=\sum_{j=0}^6 b_j\binom6j u^j(1-u)^{6-j}.
\tag{8}
$$

所有Bernstein权重非负且和为1。将 $(\rho_g,\sqrt{2\omega_s}\rho_s)_s$ 视作一个向量，整段 $R$ 因而不超过七个控制向量范数的最大值。同理，三次曲线自身的Bernstein控制点给整段 $r_{\max}$。这是全段上界，不是“残差采样很小”的经验判断。

## 4. 传到同一移动探针、钟标定及平直对照

在声明的有限动量准备域中，同一探针生成元为

$$
H_p(t)=(1+\epsilon g_z(t))\beta M+
(1+\nu g_z(t))\boldsymbol\alpha\cdot\ell_p .
\tag{9}
$$

若 $m_{\max}$、$p_{\max}$ 分别控制质量和动量，
由Duhamel及 $|g_z-\widehat g_z|\le E$，

$$
\eta_{\rm BG\to P}
\le (|\epsilon|m_{\max}+|\nu|p_{\max})
\int_0^T E(t)\,dt .
\tag{10}
$$

这对整个声明准备空间及任意参考同时成立，不只对某个末端概率。新增质量零直和部门仍被相同 $m_{\max}=0.5$ 界覆盖。

探针的数值 $U$ 另用其保存端点、导数构成矩阵Hermite曲线，计算
$\dot{\widehat U}+iH_p(\widehat g)\widehat U$。残差同为至多6次；接口用各Bernstein矩阵的Frobenius范数控制算符范数，再积分。记所得界为 $\eta_{\rm P,res}$，则

$$
\|U_p(T)-\widehat U(T)\|
\le\eta_{\rm BG\to P}+\eta_{\rm P,res}.
\tag{11}
$$

背景误差与探针自身积分误差在这里分开，没有重复将背景数值方程当作精确物理。

静止钟标定 $\tau=\int(1+\epsilon g_z)\,dt$ 也不能直接取RK4累计值为真值。Hermite曲线的 $g_z$ 积分可以精确写成每段 $\Delta\sum_{k=0}^3c_{k,z}/(k+1)$；接口将它与原累计值的差另计，得到

$$
|\tau-\tau_{\rm stored}|
\le|\epsilon|\int E+
\left|T+\epsilon\int\widehat g_z-\tau_{\rm stored}\right|.
\tag{12}
$$

同一 $\tau$ 的平直对照因此再付出
$(m_{\max}+p_{\max})|\delta\tau|$。$\xi=\int(1+\nu g_z)$ 同理；质量零传播可直接用
$p_{\max}(|\delta\xi|+|\delta\tau|)$。

**从向量误差到概率要保规范化。** 精确比较器是酉的；Hermite终点矩阵未必精确酉。可将其极分解后取酉因子，费用不超过 $\eta_{\rm P,res}$，因为它已距离同一Hermite背景下的精确酉矩阵至多该值。固定见证也可直接规范化末向量并计规范化差。若沿用原未规范化数值二次式，则另外计其范数误差对概率的影响；不能无声把非酉浮点矩阵当成精确通道。这里 $\eta_{\rm P,res}$ 约为 $7\times10^{-11}$，额外保守支付其数倍不影响下列预算。

## 5. 算术假设及接口

只读接口：[mean_field_numerical_bound.py](mean_field_numerical_bound.py)，调用：

~~~python
from mean_field_numerical_bound import run
certificate = run(steps=2048)
~~~

它重新生成固定preview端点，并输出数值包络；不写任何文件，不修改preview或其结果。当前网格要求二进时间步，确保 $T/\mathrm{steps}$ 精确可表示。

接口的包络层采用：

- 基本加减乘除后向外取相邻float；假设基础float64基本算术符合通常的IEEE舍入模型，无溢出或破坏包络的下溢。
- 向量求和用 $\gamma_{n-1}$ 绝对误差界，故不依赖NumPy具体采用串行还是树形加法；不使用点积的经验误差代替此界。
- $\pi$ 用明确的两个二进端点夹住；格点正弦用22项Taylor及小于 $4\times10^{-34}$ 的截断余项，权重用整数重数除以 $65^3$ 包络。源输入不依赖未经说明的libm正弦精度。
- 平方根采用正确舍入的通常sqrt假设并向外扩一格。小正指数及 $\varphi_1,\varphi_2$ 用正项Taylor和几何尾界，不靠大数相减。
- 残差端点数据取保存的二进有理数；物理参数则另作区间。对 $\epsilon=3/20,\nu=7/20$ 及探针质量的十进有理值，在相应float两侧各扩一格，覆盖它们与储存常量的差。积分轨道仍原样使用preview的float数据，不因包络而改参数或重拟合。

**这是公开算术模型下的区间算法，不是对CPU、NumPy构建及每条实际机器指令的形式化验证。** 同样，任意给普通浮点结果加一份padding，也不因此成为逐指令证明。正式结论应明确对所声明算术包络成立；独立复算检查实现错误，但不替代这些假设。需要其他算术保证时，可以保留本页数学证书，只替换残差上界的计算器。

## 6. 固定2048步的结果与有限预算

固定物理参数不变，$T=2$、源质量2、$\epsilon=0.15,\nu=0.35$。当前向外包络给出约：

|量|上界或包络值|
|---|---:|
|$F$|1.767086|
|整窗最大 $L_j$|4.998101|
|整窗最大Hermite Bloch残差|$3.090\times10^{-7}$|
|$\int_0^T E$|$2.525\times10^{-4}$|
|主动移动探针：背景数值误差|$5.364\times10^{-5}$|
|探针自身Hermite残差|$6.6\times10^{-11}$|
|同一钟标定后的平直对照误差|$3.381\times10^{-5}$|
|主动加平直数值算符误差|$8.746\times10^{-5}$|

因此可以给这两份预测合计预留 $9\times10^{-5}$ 的数值预算，并把剩余约 $2.5\times10^{-6}$ 中所采用的末端效果计算、规范化及额外算术包络逐项列清。它不是自动授予这些末端计算任意误差的许可。若正式check未能满足该预算，应报告失败或只加积分分辨率；不改物理参数、见证、端口来追逐通过。

本页结果只支付数值比较误差。有限 $K$、格点到连续带、实际准备和末读合同仍分别计账。尤其不把区间积分通过当作完整GR、源自身满足Einstein或现实仪器已构造的证明。
