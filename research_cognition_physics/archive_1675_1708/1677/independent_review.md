# 第1677轮独审　无阈值输运的方向与有限阈值的反例

在主证明认证的完整初值盒和整个一圈标准模型运行历史内，高能到低能的有向 $\gamma$ 漂移与 Jarlskog 不变量 $J$ 同号。$J>0$ 的 $67.070235^\circ$ 见证在这种演化下继续增大。另构造了一个带两个重矢量型下夸克的有限匹配，将同一见证的角度移到 $63.2^\circ$，同时在维四轻理论中保持六个夸克质量和三个混合角。该构造说明实际阈值和无阈值输运需要分别约束。[主算全史界](results.json) · [独审数值](independent_results.json)

## 有向四元组

令 $a_i=y_{u_i}^2$、$b_j=y_{d_j}^2$ 按质量升序排，$\tau=\log\mu/(16\pi^2)$。在瞬时质量基中取对角生成元为零，非简并夸克本征态的一圈非对角生成元为

$$
(T_u)_{ik}=-\frac32\frac{a_i+a_k}{a_k-a_i}(V b V^\dagger)_{ik},\quad
(T_d)_{lj}=-\frac32\frac{b_l+b_j}{b_j-b_l}(V^\dagger a V)_{lj},\quad
\frac{dV}{d\tau}=-T_uV+VT_d .
$$

以 $k=i+1,l=i+2$ 和 $m=j+2,n=j+1$ 循环编号，定义 $a'_i=a_k-a_l$，$b'_j=b_m-b_n$，$A_i=(a_k-a_i)(a_l-a_i)$，$B_j=(b_m-b_j)(b_n-b_j)$。利用三阶幺正矩阵的四元组虚部均为有向的 $\pm J$，其中 $J=\operatorname{Im}(V_{ud}V_{cs}V_{us}^*V_{cd}^*)$，得到

$$
\operatorname{Im}\frac{d V_{ij}/d\tau}{V_{ij}}
=-3J\frac{a'_i b'_j}{|V_{ij}|^2}\left(\frac{a_i}{A_i}+\frac{b_j}{B_j}\right).
$$

可直接复算其中的系数差，例如对固定 $i$ 的两个非对角贡献相减，
$(a_i+a_k)/(a_k-a_i)-(a_i+a_l)/(a_l-a_i)=-2a_i a'_i/A_i$。$d\gamma/d\tau$ 从 $ud,cb$ 取正号，从 $ub,cd$ 取负号，任意外部夸克重定相的导数随之消去。代码用完整反厄米基底生成元逐项对照全部九个读数，在七个正、负、零初始相位上最大绝对差小于 $1.9\times10^{-16}$。这些点检支持公式实现，以下不等式负责覆盖全域。[独审程序](independent_check.py)

## 全历史符号证明

把括号两项的有向系数分别求和，记作 $C_u,C_d$，即 $d\gamma/d\tau=-3J(C_u+C_d)$。主算给出上型相邻平方比上界 $\epsilon_u,\epsilon_c$，下型为 $\epsilon_d,\epsilon_s$，并给出 $|V_{ij}|$ 的全程管界。另需一个下界。由 Yukawa 本征值方程中普适项在对数差中严格相消，

$$
\left|\frac{d\log(b_s/b_b)}{d\tau}\right|\le3(a_t+b_b).
$$

高尺度两宽最小比乘以 $\exp[-3\int(a_t+b_b)|d\tau|]$ 给出 $b_s/b_b\ge0.000315393$。这个数字与主计算来自相同的一圈历史约束。

在 $C_u$ 中选择来自 $cd$ 的正项为单位。$ub$ 项同样有利，两个不利项 $ud,cb$ 合计相对于该单位至多为

$$
\frac{\epsilon_u}{|V_{ud}|_{\min}^2}
+\frac{\epsilon_s}{(1-\epsilon_s)|V_{cb}|_{\min}^2}
\le0.3977175 .
$$

第一项使用 $|a'_u a_u/A_u|/|a'_c a_c/A_c|=(a_u/a_c)[(a_t-a_c)/(a_t-a_u)]^2$，第二项使用 $b'_b/b'_d=(b_s-b_d)/(b_b-b_s)$，两项均用 $|V_{cd}|\le1$。

在 $C_d$ 中选择来自 $ub$ 的正项为单位。$cb$ 项有利，两个不利项 $ud,cd$ 合计相对于单位至多为

$$
\frac{\epsilon_d |V_{ub}|_{\max}^2}{(b_s/b_b)_{\min}(1-\epsilon_d)^2}
\left[\frac1{|V_{ud}|_{\min}^2}
+\frac1{(1-\epsilon_c)|V_{cd}|_{\min}^2}\right]
\le0.0396859 .
$$

这里使用 $|b'_d b_d/B_d|/|b'_b b_b/B_b|=(b_d/b_b)[(b_b-b_s)/(b_s-b_d)]^2$。因此 $C_u,C_d$ 在整条运行管上均为正。实 Yukawa 的 CP 守恒子空间在标准模型一圈方程下封闭，解的唯一性使非零 $J$ 不能穿越该子空间。上述符号结论适用于连续选取的有向角；最终与实验相比较时再按 $180^\circ$ 取模。

## 规范允许的有限匹配

给定 $V_0=V(\theta_{12},\theta_{13},\theta_{23},\delta_0)$ 及 $Y_0=V_0\operatorname{diag}(y_d,y_s,y_b)$。保持三个 PDG 混合角，单独换取能使树级 $\gamma=63.2^\circ$ 的 $\delta_1$，得 $Y_1$。PDG 分解中相位只影响 $13$ 旋转的两个非对角元，所以 $\Delta=Y_1-Y_0$ 的秩至多为二，且

$$
\frac{\|\Delta\|_2}{y_b}\le2\sin\theta_{13}\left|\sin\frac{\delta_1-\delta_0}{2}\right|.
$$

对 1675 两宽见证，实际左边为 $0.00029518$。为检验此改变可由实际规范允许的端口承担，引入两个 $(\mathbf3,\mathbf1,-1/3)$ 矢量型 Dirac 下夸克。令 $R$ 的两行为 $\Delta$ 的两个单位右奇异行向量，各乘 $\kappa=0.01$，$C=(I+R^\dagger R)^{1/2}$。取

$$
L=(Y_0-Y_1 C)R^\dagger/\kappa^2,\qquad
\mathcal M=\begin{pmatrix}vY_0&vL\\M R&M I_2\end{pmatrix}.
$$

这里 $v$ 取轻质量等于 $vy$ 的约定。上排右块可由 $\bar Q_L H D_R$ 产生，下排均为规范允许的单态质量。$P=R^\dagger R/\kappa^2$ 投影到改变的右行空间，满足 $(Y_0-Y_1C)(I-P)=0$。在电弱破缺前，下排作用于正则归一的轻右子空间 $\binom{I}{-R}C^{-1}$ 恰为零。因而维四匹配关系精确给出

$$
Y_{\rm light}=(Y_0-L R)C^{-1}=Y_1.
$$

$\|L\|_2=1.779\times10^{-4}$。完整破缺质量矩阵在 $M/v=100,1000$ 的轻态 SVD 上分别给出 $63.19999999^\circ$ 和 $63.20000000^\circ$，各轻质量的最大相对差小于 $2\times10^{-10}$。这些有限值是对 EFT 近似的数值检查；精确相等的是维四匹配。该重块的耦合按目标相位选取，展示了阈值自由可达到的量级，没有独立选定这些耦合，也未将它与高尺度带电轻子 Clebsch 关系的运行作统一拟合。若在高尺度加这种匹配，再在低能采用无阈值一圈 SM，尚需计入主算给出的不足 $0.060^\circ$ 后续漂移。第 1675 轮已计算三条具体无阈值轨迹，本轮新量是全域漂移的符号与可显式实现的有限阈值竞争例。[此前的点检](../1675/independent_review.md)
