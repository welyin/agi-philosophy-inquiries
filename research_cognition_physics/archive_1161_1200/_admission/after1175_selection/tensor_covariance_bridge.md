# 双关系自旋候选的十分量协方差桥

2026-10-13。计算前辅助，不编号，不新增认知公理或科学组。只给代数省算与小格点校准，不签大格点概率、数值残差包络、实际仪器或GR。

配套：[独立小矩阵代码](tensor_covariance_check.py)、[只写一次结果](tensor_covariance_checks.json)。本页不改变作者候选参数；校准使用另列的小格点窗口。

## 1. 符号固定及精确闭合

依1175取四个Hermitian Clifford矩阵，

$$
\Gamma=(\beta,\alpha_x,\alpha_y,\alpha_z),\qquad
\{\Gamma_a,\Gamma_b\}=2\delta_{ab}I,\qquad
L_{ab}=i\Gamma_a\Gamma_b\quad(a<b).
$$

$L_{ab}$ 的符号不能换成 $-i\Gamma_a\Gamma_b$ 后继续沿用下面的方程。取每个动量处的占据协方差

$$
P=\frac12\left(I+\sum_{a=0}^3r_a\Gamma_a+
\sum_{a<b}s_{ab}L_{ab}\right),\qquad s_{ba}=-s_{ab},\quad s_{aa}=0.
$$

十个系数均实，提取式为

$$
r_a=\frac12\operatorname{tr}(\Gamma_aP),\qquad
s_{ab}=\frac12\operatorname{tr}(L_{ab}P).
$$

这不是任意 $4\times4$ 协方差的全参数化；它覆盖本候选的初态，并在以下生成元下精确不变。不需要三重或四重Clifford分量。

对 $h=d_a\Gamma_a$，由 $\dot P=-i[h,P]$ 得

$$
\boxed{\quad
\dot r_a=-2\sum_b s_{ab}d_b,\qquad
\dot s_{ab}=2(d_b r_a-d_a r_b)\quad(a<b).
\quad}
$$

证明仅用

$$
[\Gamma_a,\Gamma_b]=-2iL_{ab},\qquad
[\Gamma_c,L_{ab}]
=2i(\delta_{ca}\Gamma_b-\delta_{cb}\Gamma_a).
$$

例如 $d_0>0,r_1=1$ 时 $\dot s_{01}=-2d_0$；这对应通常 $\beta$ 旋转产生正的相应Pauli分量，因为本页 $L_{01}=i\beta\alpha_x$ 已含负号。代码直接比较完整复矩阵对易子，不靠记忆判断此符号。

## 2. 两份准备及完整来源

记 $k_i=\sin(\ell p_i)/\ell$，

$$
\begin{aligned}
N&=1+\epsilon z_s,\\
v_x&=1+\nu z_s+g z_t,\qquad
v_y=1+\nu z_s-g z_t,\qquad
v_z=1+\nu z_s,\\
d&=(Nm_S,v_xk_x,v_yk_y,v_zk_z).
\end{aligned}
$$

源A在 $k_x\ne0$ 时为

$$
P_A(p)=\frac12[I+\operatorname{sgn}(k_x)\alpha_x],
\qquad r_x=\operatorname{sgn}(k_x),\quad r_{a\ne x}=s_{ab}=0.
$$

在整个零面 $k_x=0$（不只零动量点）取 $P_A=(I+\beta)/2$。源B把 $x$ 换成 $y$，并在 $k_y=0$ 面采用同一 $\beta$ 投影。这是明确的准备输入；离散零面占有有限权重，不能删掉。原65的奇数格点中可以用整数动量标签判断正负和零，避免把舍入的 $\sin0$ 当物理分支。

每个 $P$ 均为秩2投影。完整源Slater占据数是 $2V$，不是取一个正能支或两能级近似。由Clifford迹身份，

$$
E_m=2m_S\langle r_0\rangle,\qquad
E_{t_i}=2\langle k_i r_i\rangle .
$$

角括号是全部格点动量的归一和，不是只取低带。若另外合并等价动量类，必须保对应权重及双旋转后的初态，不能直接复用1175各向同性 $|\mathbf k|$ 分类。

关系Bloch向量 $\mathbf g_s,\mathbf g_t$ 满足

$$
\begin{aligned}
b_s&=\epsilon E_m+\nu(E_{t_x}+E_{t_y}+E_{t_z}),&
\dot{\mathbf g}_s&=2(-\Omega_s,0,b_s)\times\mathbf g_s,\\
b_t&=g(E_{t_x}-E_{t_y}),&
\dot{\mathbf g}_t&=2(-\Omega_t,0,b_t)\times\mathbf g_t .
\end{aligned}
$$

同一Hartree轨道的每站能量是

$$
\mathcal E=NE_m+v_xE_{t_x}+v_yE_{t_y}+v_zE_{t_z}
-\Omega_s x_s-\Omega_t x_t .
$$

固定系数时每个 $\operatorname{tr}(hP)$ 的交换子贡献为零；系数变化贡献 $b_s\dot z_s+b_t\dot z_t$，恰被两份关系自能变化抵消。因此 $\dot{\mathcal E}=0$。投影性、迹2、关系Bloch长度、完整CAR的粒子数／奇偶与全局粒子相位对称也保持。十分量的二次不变量

$$
\sum_a r_a^2+\sum_{a<b}s_{ab}^2=1
$$

由投影性及演化推出，但它单独不是任意十分量对应投影的充分条件。也不把这些守恒误称为已经建立局部规范约束或一般GR约束。

## 3. A/B对称及物理比较方式

取绕 $z$ 的正 $\pi/2$ 旋转

$$
R(p_x,p_y,p_z)=(-p_y,p_x,p_z),\qquad
U_R=\exp[-\pi\alpha_x\alpha_y/4].
$$

则 $U_R\alpha_xU_R^\dagger=\alpha_y$、$U_R\alpha_yU_R^\dagger=-\alpha_x$，并且

$$
P_B(p,t)=U_RP_A(R^{-1}p,t)U_R^\dagger,\quad
\mathbf g_{s,B}=\mathbf g_{s,A},\quad
\mathbf g_{t,B}=(x_{t,A},-y_{t,A},-z_{t,A}).
$$

最后一个变换来自实际关系自旋的 $X_t$ 共轭。初始两份 $+X$ 均保持。标量来源不变，张量来源变号，$E_{t_x},E_{t_y}$ 交换；初始与后续背景能源相同。立方格点及各向同性标量规则是这条对称的前提。

这是不含未知probe反馈的Hartree背景对称。完整父规则在连probe、端口及记录也共同旋转时有相应对称；若固定probe和末读，仅旋转源准备，比较的是不同相对方向，不能由共同换表把它们宣布等价，也不能保证任意读口都有非零差。

$g=0$ 确实恢复1175的标量耦合规则，并使初始 $+X_t$ 的额外关系自旋无效；**新的A/B源Slater不是1175原来逐站 $\beta$ 正占据的准备**，故不能继承旧具体轨道或概率。

## 4. 有限K证明直接复用；不得只算单粒子范数

令每份cluster含原完整源CAR和两份关系自旋。取

$$
\begin{aligned}
h_r&=H_{S,0,r}-\sum_i(\Omega_sX_{s,ir}+\Omega_tX_{t,ir}),\\
W_{ru}&=\sum_{i,\sigma=s,t}
(Z_{\sigma,ir}F_{\sigma,iu}+F_{\sigma,ir}Z_{\sigma,iu}),\\
F_{s,i}^a&=\epsilon m_i^a+\nu\sum_b t_{i,b}^a,\qquad
F_{t,i}^a=g(t_{i,x}^a-t_{i,y}^a).
\end{aligned}
$$

这里 $t_{i,b}$ 是该方向两相邻边的各半份，$r=u$ 仍排除。所有跨复制因子在CAR上是偶算符；两种关系的标签只扩大旧证明的求和指标。因此[原有限K证明](../after1174_direction/mean_field_bound.md)逐项适用，无需新计数定理。完整Fock的充分界为

$$
\begin{aligned}
b_S&=2|\epsilon|m_S+\frac{6|\nu|+4|g|}{\ell},\\
b_P&=2|\epsilon|(0.3+0.5)+\frac{18|\nu|+12|g|}{\ell},\\
w&=2Vb_S\ge\|W\|,\qquad
\mathcal L=Vb_P\ge\sum_{i,\sigma}\|F_{\sigma,i}^P\|\,\|Z_{\sigma,i}\|.
\end{aligned}
$$

源的一个方向站点动能界为 $2/\ell$；同父探针有三份Dirac味（含质量零），所以为 $6/\ell$。不能用质量零使其跳跃Fock范数消失。两方向差用三角不等式得 $4/\ell$ 和 $12/\ell$，已足够，无需优化。

对于固定已知跨复制乘积准备、任意未知probe／参考，

$$
\|(U_{\rm act}-U_{\rm cmp})E_K\|_{\rm op}
\le\frac{4\mathcal L T e^{6wT}}{\sqrt K}.
$$

比较器保留**精确BG联合演化**，实际反作用没有删除；有限K的少一份自耦已在旧证明中计入。新的 $w,\mathcal L$ 须重算，不能沿用1175已选的K充分值。$\Omega_s,\Omega_t$ 属单体 $h$，不进入该粗指数，但全部资源与正谱平移仍要计账。双关系来源也不因同一总H而自动成为全部协变应力或普适几何自耦。

## 5. 小格点核验与停止线

独立测试只取 $5^3$ 模、$\ell=.8,T=.4$，$\epsilon=.15,\nu=.35,g=.12,\Omega_s=\Omega_t=.5,m_S=2$，160步RK4。它不是作者候选的大格点参数选择或严格误差包络。

- 随机37份十分量对完整对易子：最大差 $2.51\times10^{-16}$；十维闭合差 $1.12\times10^{-16}$。
- 自洽十分量与原 $4\times4$ 复矩阵协方差共同推进：最大差 $1.45\times10^{-15}$。
- A/B正旋转差 $2.23\times10^{-16}$；同源能源相同，RK4能量漂移约 $6.33\times10^{-12}$。
- 投影误差约 $3.76\times10^{-12}$。小局部单粒子谱经完整Fock范数公式给 $2.5,7.5$，均小于来源／三味probe的充分界 $5,15$。

数值只是符号、归一、A/B和完整Fock责任交叉检查。正式信号仍需作者的实际读口、有限K预算、格点／连续字典、数值包络及仪器误差。到此停止辅助计算，不用小格点给大格点签收。

## 6. 双关系系统的误差稳定常数不必相加

令 $q=(r_0,\ldots,r_3,s_{01},\ldots,s_{23})\in\mathbb R^{10}$，上述源方程写成 $\dot q=A(d)q$。$A(d)$ 为实斜对称矩阵，并满足

$$
\|A(\delta d)\|_{\rm op}=2|\delta d|
\quad(\delta d\ne0).
$$

证明可用正交Clifford旋转把 $\delta d$ 化成 $(|\delta d|,0,0,0)$。此时三个 $(r_i,s_{0i})$ 平面各以频率 $2|\delta d|$ 转动，$r_0$ 与三个 $s_{ij}$ 不动。旋转在四向量加六二形式的系数范数上正交，所以该算符范数在原坐标相同。

记

$$
f_s=(\epsilon m_S,\nu k_x,\nu k_y,\nu k_z),\qquad
f_t=(0,gk_x,-gk_y,0),\qquad
C_{\lambda\mu}=\langle f_\lambda\cdot f_\mu\rangle,\qquad
F_{\rm op}=\sqrt{\lambda_{\max}(C)}.
$$

平均与第2节一致，是完整动量网格的非负归一权；若按对称类合并，必须精确运输权重，不能以任意截样本替代。设精确解与一份可微近似轨道之差为 $e_q,e_{g_s},e_{g_t}$，定义

$$
G^2=|e_{g_s}|^2+|e_{g_t}|^2,\quad
Q^2=\langle|e_q|^2\rangle,\quad
E^2=G^2+2Q^2 .
$$

取近似轴对应的斜对称旋转放在齐次误差项，精确态放在强迫项。源侧强迫为 $A(\sum_\lambda\delta z_\lambda f_\lambda)q_{\rm exact}$；两份关系侧强迫为 $4\langle f_\lambda\cdot e_r\rangle e_z\times g_{\lambda,\rm exact}$。精确 $|q|=|g_\lambda|=1$，而Gram矩阵同时给

$$
\left\langle\left|\sum_\lambda u_\lambda f_\lambda\right|^2\right\rangle
\le F_{\rm op}^2|u|^2,\qquad
\sum_\lambda|\langle f_\lambda\cdot e_r\rangle|^2
\le F_{\rm op}^2 Q^2 .
$$

因此斜对称部分消去后，

$$
\frac12\frac{dE^2}{dt}
\le 8F_{\rm op}GQ+E R
\le 2\sqrt2F_{\rm op}E^2+ER,\qquad
D^+E\le LE+R,\quad L=2\sqrt2F_{\rm op}.
$$

这里 $D^+$ 可处理 $E=0$，且近似方程的残差按同一权重合为

$$
R^2=|\rho_{g_s}|^2+|\rho_{g_t}|^2+2\langle|\rho_q|^2\rangle.
$$

于是

$$
E(t)\le e^{Lt}E(0)+\int_0^t e^{L(t-u)}R(u)\,du .
$$

若改用近似 $q$ 而非精确 $q$ 控制源侧强迫，也可得较粗的 $L=\sqrt2F_{\rm op}(1+\sup|\widehat q|)$；没有必要为了两份关系自旋直接把两个独立常数相加。

完整均匀立方周期网格（每方向不少于3点）上，

$$
C_{st}=0,\qquad
C_{ss}=\epsilon^2m_S^2+\frac{3\nu^2}{2\ell^2},\qquad
C_{tt}=\frac{g^2}{\ell^2}.
$$

采用1175的 $\ell=16/65,\epsilon=.15,\nu=.35,m_S=2$，并对拟议的 $g=.2$，有 $F_t=.8125<F_s\approx1.7670859553$；因而 $F_{\rm op}=F_s$，稳定指数无需超过原单关系值。此数值是解析常数的代入，不是新轨道或概率认证。

实际近似轨道的初值误差 $E(0)$、十实分量残差、完整权重、参数舍入与算术包络仍须另核。方程仍为二次，逐步三次Hermite轨道的残差仍最高六次，因此旧Bernstein方法可复用；**旧已算出的残差数值不能原样复用。** 这只提供新的省算稳定界，不承诺候选留出必能超过总误差。
