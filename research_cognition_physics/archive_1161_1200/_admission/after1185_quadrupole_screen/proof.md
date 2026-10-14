# 1185后未编号筛查：完整约束核的双紧支四极配对

2026-10-13。只签数学源与核；不新增1186，不认定实际仪器、局部制备权限或动态系数已被选择。固定约定与计算见[check.py](check.py)、[results.json](results.json)。内部代理曾交叉手算，不是外部盲审。

## 1. 完整核与守恒源

采用[940 §3](../../../archive_935_955/research_note_940.md)的符号，$D=k^2-\omega^2$，传播分母统一取退迟处方。比较族为

$$
K_\xi(T_A,T_B)=\frac{T_{A,\mathrm{TT}}:T_{B,\mathrm{TT}}}{D_{\mathrm{ret}}}
+\xi\frac{C_g(T_A,T_B)}{k^2},
$$

$$
C_g=\frac12\rho_A\rho_B+\frac12(\rho_A\sigma_B+\sigma_A\rho_B)
-2\mathbf v_A\cdot\mathbf v_B+\frac12v_{A,L}v_{B,L}.
\tag{1}
$$

这里写的是实源的双线性积分；Fourier表示中第一份源按同一约定取复共轭。式(1)与横动量系数 $-2$、纵动量系数 $-3/2$ 的原式等价，未丢应力迹或动量项。

对任一对称紧支 $S_{ij}$，取

$$
\rho=\partial_i\partial_jS_{ij},\qquad
v_i=-\partial_t\partial_jS_{ij},\qquad
T_{ij}=\partial_t^2S_{ij}.
\tag{2}
$$

交换偏导即得 $\partial_t\rho+\partial_iv_i=0$、$\partial_tv_i+\partial_jT_{ji}=0$。这是数学守恒身份；不说明这些函数已由同一真实材料、控制器和参考制备。

## 2. 双四极与精确空间系数

固定 $e=\operatorname{diag}(1,-1,0)$。两形状 $f_A,f_B$ 非负、归一、径向紧支，半径分别为 $r_A,r_B$，中心沿z轴相距 $R>r_A+r_B$。令

$$
S^A_{ij}=a(t)e_{ij}f_A,\quad S^B_{ij}=b(t)e_{ij}f_B,\qquad
m_2=\int |x|^2f_A(x)\,d^3x+\int |x|^2f_B(x)\,d^3x.
$$

本例 $\sigma_A=\sigma_B=0$ 来自 $\operatorname{tr}e=0$，不是普遍省略应力迹。设 $\mathscr D=\partial_x^2-\partial_y^2$，

$$
G_1=(-\Delta)^{-1}=\frac1{4\pi R},\qquad
G_2=(-\Delta)^{-2}=-\frac R{8\pi}.
$$

所需 $G_2$ 总伴四阶空间导数；标准基本解的低阶多项式零模歧义不影响这些配对。所有空间积分使用相同的 $d^3k/(2\pi)^3$ 约定。

对式(1)分部积分，完整约束配对为

$$
\mathcal I_C=\frac12\mathscr D^2\overline G_1\int ab\,dt
+\left[2(e^2):\nabla\nabla\overline G_1
+\frac12\mathscr D^2\overline G_2\right]\int a'b'\,dt.
\tag{3}
$$

上横线表示双形状平均。动量直接项的正导数系数来自两个空间梯度在分部积分后的负号，不能仅计算 $\rho_A\rho_B$。在z轴上，未经形状平均的导数为

$$
\mathscr D^2G_1=\frac3{\pi R^5},\quad
2(e^2):\nabla\nabla G_1=-\frac1{\pi R^3},\quad
\frac12\mathscr D^2G_2=\frac1{4\pi R^3}.
\tag{4}
$$

### 径向平均为何精确

若函数在整个平均球内调和，其球面平均等于中心值；若双调和，则球面平均恰为中心值加 $r^2\Delta u/6$。后一式可由球面平均的径向双调和方程和原点正则性得到，不是截断Taylor近似。用径向形状对球面半径再积分，由 $\Delta G_2=-G_1$、$\Delta G_1=0$，并连续作两次平均，得到

$$
\overline G_1=G_1,\qquad
\overline G_2=G_2-\frac{m_2}{6}G_1.
\tag{5}
$$

分离支撑保证平均区域避开原点；各导数可与平均交换。两径向形状的一阶矩为零，所以其差分布的二阶矩就是两二阶矩之和。结合(3)—(5)，

$$
\boxed{\displaystyle
\mathcal I_C=\frac1{4\pi R^5}
\left[6\int ab\,dt-(3R^2+m_2)\int a'b'\,dt\right].}
\tag{6}
$$

## 3. 有限窗的严格非零与有理见证

令 $a=b\ne0$，支撑在同一长度T的时间窗，且 $T<R-r_A-r_B$。两整个时空支撑因此处处类空分离。对 $a\in H^1_0(0,T)$，

$$
\int a'^2\,dt\ge\frac{\pi^2}{T^2}\int a^2\,dt,
\qquad
6-(3R^2+m_2)\frac{\pi^2}{T^2}<0.
\tag{7}
$$

故 $\mathcal I_C<0$。这覆盖任意非零 $C^\infty$ 紧支同包络，不依赖一次偶然非零采样。这里只使用径向非负形状的 $m_2\ge0$；$T<R$ 与 $\pi>3$ 已足以严格判号。

一份可精确积分的有限见证是

$$
a(t)=b(t)=q_0p(t/T),\qquad p(u)=u^3(1-u)^3\quad(0<u<1),
$$

外延为零。$p,p',p''$ 在两端均为0，所以这是 $C^2$ 紧支源；不把它称为 $C^\infty$。直接有理多项式积分给

$$
\int a^2dt=\frac{q_0^2T}{12012},\qquad
\int a'^2dt=\frac{q_0^2}{770T},
$$

$$
\boxed{\displaystyle
\mathcal I_C=\frac{q_0^2}{4\pi R^5}
\left[\frac T{2002}-\frac{3R^2+m_2}{770T}\right]<0.}
\tag{8}
$$

代码另固定一份数学实例：$R=3,r_A=r_B=1/4,T=1,q_0=1$；空间形状为归一的 $(1-|x|^2/r_A^2)^4_+$，是 $C^3$ 紧支函数。每形状二阶矩为 $3r_A^2/13$，故 $m_2=3/104$。代码精确给出 $\pi\mathcal I_C$ 的负有理数，不引入点源、无限时长或数值积分误差。这里没有宣称源满足正能条件或可由量子材料制备。

## 4. 真正接触退化与读口角节点

取另一紧支源

$$
S_{ij}=(\partial_i\partial_j-\delta_{ij}\Delta)\chi.
$$

则 $\rho=\mathbf v=0$、$\sigma=-2\Delta\ddot\chi$。空间张量正比于横向投影的纯迹，因此 $T_{\mathrm{TT}}=0$。与任意该核允许的守恒B源配对时，只剩

$$
\mathcal I_C=\int\ddot\chi(t,x)\rho_B(t,x)\,dt\,d^3x.
\tag{9}
$$

这是局部接触项，空间支撑分离便严格为零，对任意ξ仍如此。不能把所有双散度源都归为这种退化类；第2—3节已给反例。

另有较弱的仪器盲点：对第2节的e，在z轴处 $\mathscr DG_1=0$。一个静止、仅含密度的点测试源因此读不到此约束配对；但双四极读口由 $\mathscr D^2G_1$ 等项给非零。这里的静止点测试不是一份可任意启停且仍守恒的完整钟仪器；更不能由这个角节点推断源本身只有接触效应。

## 5. 完整退迟相消与下一接口

在同一守恒源和退迟处方下，$\xi=1$ 的完整核等于协变GR领先交换核，仅含局域导数作用于退迟波核。其对处处类空分离的两份紧支源配对为0。成熟退迟交换的约定见[Hertzberg—Litterer—Sandora §II](https://arxiv.org/html/2005.01744#S2)；本页没有从该文领取两传播模或实际仪器的认知推导。

因此在这里的数学源类中

$$
\mathcal I_{\mathrm{TT,ret}}=-\mathcal I_C,\qquad
\mathcal I_\xi=(\xi-1)\mathcal I_C.
\tag{10}
$$

TT空间投影本身非局部；单独TT项或单独约束项都不能作为GR超光速通信证据。式(10)也不是一份测得的概率差：从交换核到实际局部动态四极、控制/触发、相同初始约束资料、完整量子响应及原记录，还须在同一个过程中建立误差受控的映射。一般末端门编译或运输算符不会自动保持完整来源插入。

本页仅排除“该双紧支四极方向全部是接触项”的停止理由。没有选择ξ，没有构造健康自旋2或Einstein动力学，没有完成未知输入/R的实际通道，也不预占下一轮。代码默认只读比较，只有显式 `--write` 才能独占创建结果；本三件供后续独立数学审阅。
