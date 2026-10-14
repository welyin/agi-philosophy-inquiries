# 1193证明：完整分布应力到有限时间质心轨迹的条件连接

2026-10-14。本文不输入Einstein方程，也不把几何运动方程作为物质的独立轨迹假设。它从指定共同度规下的完整体应力平衡及公开的有限矩合同，给出应力质心与该度规测地线之间的有限时间界。该质心尚不自动等于原CO读口；比较使用的坐标时也尚不自动等于原内部钟的固有时。

## 1. 准入、旧结果和成熟来源

[1192](../research_note_1192.md)把有限钟束端点接到弱Ricci判据，但其“原自由腿按同一度规作仿射测地运动”仍是额外有效合同。本轮只提供其中一段可独立检验的充分连接：

$$
\text{同一度规下的完整体应力平衡和有限矩控制}
\ \Longrightarrow\
\text{应力质心在有限坐标时内接近该度规的测地线}.
\tag{1}
$$

[937](../../archive_935_955/research_note_937.md)的尘模型、[953](../../archive_935_955/research_note_953.md)的给定世界线来源、[941](../../archive_935_955/research_note_941.md)—[942](../../archive_935_955/research_note_942.md)的完整内部质量，以及[1119](../../archive_1106_1119/1119/proof.md)的完整实际材料运输均直接复用。没有重新证明质量比例，也没有把给定世界线再次变分称为分布物质运动的来源。

成熟入口是Geroch与Weatherall的[《The Motion of Small Bodies in Space-time》§1—3](https://arxiv.org/html/1707.04222)。该文在固定时空中，以能量条件、守恒与tracking/分布极限联系小体运动和测地线。本轮采用其“先保完整应力，再判断整体运动”的思路，但下述有限界另有显式矩假设，不是该文极限定理的直接有限误差版本，也不声称原创普遍小体运动定理。

## 2. 固定对象及公开合同

### 2.1 同一图、同一度规、同一材料体

在固定操作坐标图 $(t,x^1,x^2,x^3)$ 中取一个光滑有效Lorentz度规 $g$，号差为 $+---$。所用 $t=\mathrm{const}$ 切片为类空，比较测地线沿未来具有 $dt/ds>0$。令

$$
\tau^{\mu\nu}(t,x)=\sqrt{|\det g(t,x)|}\,T_B^{\mu\nu}(t,x)
\tag{2}
$$

为指定体 $B$ 的对称密度化应力。$B$ 包括事前分配给它的内部钟、程序、记录、结合能、控制及相互作用承担者。体与其余环境的交换以 $f^\nu$ 公开保留：

$$
\partial_\mu\tau^{\mu\nu}
=-\Gamma^\nu_{\alpha\beta}\tau^{\alpha\beta}+f^\nu.
\tag{3}
$$

式(3)是 $\nabla_\mu T_B^{\mu\nu}=F^\nu$ 的坐标密度形式，$f^\nu=\sqrt{|g|}F^\nu$。这里尚未证明原实际量子父过程具有这份共同度规及应力字典；式(3)是一项明示的有效平衡合同。若从给定协变物质作用导入，还须满足相应物质方程并保留非最小耦合的完整来源。

以下空间积分采用 $d^3x$。允许非紧支尾部，但要求全部出现的绝对矩可积，分部积分的无穷远边界通量消失，积分与时间导数可合法交换。有限窗口的边界交换必须保留到 $f$ 或显式补入相同预算，不能丢弃。可用光滑场证明后按受控弱极限推广；不免费假定点粒子或严格紧支。

设在 $0\le t\le T$ 上

$$
\begin{aligned}
E&=\int\tau^{00}\,d^3x>0,\\
P^\nu&=\int\tau^{0\nu}\,d^3x,\qquad P^0=E,\\
X^i&=\frac{\int x^i\tau^{00}\,d^3x}{E},\\
w^0&=1,\qquad w^i=\frac{P^i}{E}.
\end{aligned}
\tag{4}
$$

$X$ 是这份固定坐标与应力分配下的能量质心。$E>0$ 本身不要求 $\tau^{00}$ 逐点非负，也不保证 $X$ 在材料支撑内；后述有效域必须显式包含 $X$。不得根据候选结果重新分配体/环境或重选坐标。

### 2.2 三项矩合同与交换合同

对所有时刻声明

$$
\Sigma^{\alpha\beta}
=\frac1E\int\tau^{\alpha\beta}\,d^3x-w^\alpha w^\beta,
\qquad
\sum_{\alpha,\beta=0}^{3}|\Sigma^{\alpha\beta}|\le\sigma.
\tag{5}
$$

由对称性与式(4)，$\Sigma^{00}=\Sigma^{0i}=\Sigma^{i0}=0$，但以下统一对16个有序分量求和。式(5)约束剩余集成应力、内部速度离散及非尘项；不能从“整体正能”或“所有资源均已记账”自动领取。

另声明

$$
\frac1E\int |x-X|_\infty
\sum_{\alpha,\beta=0}^{3}|\tau^{\alpha\beta}|\,d^3x\le r,
\tag{6}
$$

$$
\frac1E\int |x-X|_\infty|f^0|\,d^3x\le r_f,
\qquad
\max_{\nu=0,1,2,3}
\frac1E\int|f^\nu|\,d^3x\le\eta.
\tag{7}
$$

$r$ 是完整应力的加权绝对一阶矩，不只是 $\tau^{00}$ 的几何半径；$r_f$ 控制交换能源的作用位置。均匀小体尺寸或小总力可作为这些合同的某些充分实现，但不是其定义，也不在此证明其可实现。

全体系守恒不等于子体 $f=0$。环境、介质、规范场和暗部门的动量交换均需纳入 $f$ 及全体系另一侧的反向收支；不能把它们删去后宣称自由落体。本体内所有应力也须保留在式(5)—(6)。

### 2.3 几何和轨迹域

采用分量最大范数，假定

$$
|\Gamma^\nu_{\alpha\beta}(t,x)|\le B,
\qquad
|\Gamma^\nu_{\alpha\beta}(t,x)-\Gamma^\nu_{\alpha\beta}(t,y)|
\le L|x-y|_\infty.
\tag{8}
$$

后一式是相同坐标时下的空间Lipschitz界。界须覆盖材料积分所在区域、$X$、比较测地线及它们之间的连接线段；无限尾的使用必须有相应全尾域控制或另付尾误差。不能各点另换正规坐标使 $\Gamma$ 消失。

令 $v\ge0$ 控制三个空间分量：

$$
\max_i|w^i(t)|\le v,\qquad \max_i|V^i(t)|\le v,\qquad S=1+3v.
\tag{9}
$$

$w^0=V^0=1$ 不受 $v<1$ 的空间界限制。$|w|_1,|V|_1\le S$。这里 $v<1$ 也不独立保证曲度规中的类时性；比较初值的类时资格单列。所有界若只在局部域成立，以下结论仅到首次退出，除非另外给出缓冲区和不退出的bootstrap证明。

## 3. 精确矩方程：不对误差再求导

写

$$
s^\nu=-\Gamma^\nu_{\alpha\beta}\tau^{\alpha\beta}+f^\nu.
\tag{10}
$$

由式(3)及边界条件，

$$
\dot P^\nu=\int s^\nu\,d^3x,\qquad
\dot E=\int s^0\,d^3x.
\tag{11}
$$

对一阶空间矩分部积分，并使用 $\tau^{i0}=\tau^{0i}$：

$$
\frac d{dt}\int x^i\tau^{00}\,d^3x
=P^i+\int x^is^0\,d^3x.
\tag{12}
$$

于是精确得到

$$
\dot X^i=w^i+d^i,\qquad
d^i=\frac1E\int(x^i-X^i)s^0\,d^3x,
\tag{13}
$$

且

$$
|d|_\infty\le Br+r_f=:D.
\tag{14}
$$

**$\dot X$ 一般不等于 $P/E$。** 能源在不同位置获得或流失会移动能量质心。本文不要求 $\dot d$ 存在或很小，因而不把式(13)偷偷求导成一个已有控制的二阶 $X$ 方程。

在 $(t,X(t))$ 处定义

$$
C^\nu=\Gamma^\nu_{\alpha\beta}(t,X)w^\alpha w^\beta.
\tag{15}
$$

利用式(5)—(8)分解

$$
\begin{aligned}
\frac1E\int\Gamma^\nu_{\alpha\beta}(t,x)\tau^{\alpha\beta}\,d^3x-C^\nu
={}&\Gamma^\nu_{\alpha\beta}(t,X)\Sigma^{\alpha\beta}\\
&+\frac1E\int[
\Gamma^\nu_{\alpha\beta}(t,x)-\Gamma^\nu_{\alpha\beta}(t,X)]
\tau^{\alpha\beta}\,d^3x.
\end{aligned}
\tag{16}
$$

右端绝对值至多 $B\sigma+Lr$。故有

$$
\frac{\dot P^\nu}{E}=-C^\nu+\epsilon^\nu,
\qquad
|\epsilon^\nu|\le B\sigma+Lr+\eta=:A.
\tag{17}
$$

对比值 $w^i=P^i/E$ 求导得到

$$
\dot w^i
=-C^i+w^iC^0+e^i,
\qquad e^i=\epsilon^i-w^i\epsilon^0,
\qquad |e|_\infty\le(1+v)A=:W.
\tag{18}
$$

因此整段轨迹满足一阶受迫系统

$$
\boxed{
\dot X=w+d,\quad
\dot w=F_\Gamma(t,X,w)+e,\quad
|d|_\infty\le D,\quad |e|_\infty\le W,
}
\tag{19}
$$

其中

$$
F_\Gamma^i(t,x,z)
=-\Gamma^i_{\alpha\beta}(t,x)z^\alpha z^\beta
+z^i\Gamma^0_{\alpha\beta}(t,x)z^\alpha z^\beta,
\qquad z^0=1.
\tag{20}
$$

这一推导只使用体应力平衡与矩合同，没有使用Einstein方程、点质量世界线、物质最小耦合作用或量子能量条件。

## 4. 同一度规的坐标时测地线与有限比较界

### 4.1 坐标时的非仿射项必须保留

设 $Y(t)$ 是同一 $g$ 的类时测地线，$V=dY/dt$。仿射参数写作 $s$，消去 $d^2t/ds^2$ 后得到

$$
\dot Y=V,\qquad \dot V=F_\Gamma(t,Y,V).
\tag{21}
$$

式(20)的第二项正是改用坐标时后的非仿射项。不能略去它，也不能据式(21)认领原内部钟已经是 $s$。

默认比较初值为 $Y(0)=X(0)$、$V(0)=w(0)$。若原实验提供的是 $\dot X(0)$，则除仪器误差外，还须把式(14)给出的至多 $D$ 差异计入初速匹配。一般初值误差记为

$$
z_x(0)\le e_x,\quad z_v(0)\le e_v,\qquad
z_x=|X-Y|_\infty,\quad z_v=|w-V|_\infty.
\tag{22}
$$

### 4.2 向量场的显式Lipschitz常数

在式(8)—(9)的共同域内，对空间位置变化有

$$
|F_\Gamma(t,x,z)-F_\Gamma(t,y,z)|_\infty
\le A_x|x-y|_\infty,\qquad
A_x=(1+v)LS^2.
\tag{23}
$$

对两个空间速度，利用 $\delta z^0=0$ 和
$\sum_{\alpha=0}^{3}|\delta z^\alpha|\le3|\delta z|_\infty$，

$$
|\Gamma_{\alpha\beta}(z^\alpha z^\beta-\widetilde z^\alpha\widetilde z^\beta)|
\le6BS|z-\widetilde z|_\infty.
\tag{24}
$$

再保留式(20)外乘的 $z^i$ 之变化，得到

$$
|F_\Gamma(t,x,z)-F_\Gamma(t,x,\widetilde z)|_\infty
\le A_v|z-\widetilde z|_\infty,
\qquad
A_v=B[S^2+6(1+v)S].
\tag{25}
$$

将位置与速度之差分开，式(19)、(21)给出积分不等式

$$
\begin{aligned}
z_x(t)&\le e_x+\int_0^t[z_v(s)+D]\,ds,\\
z_v(t)&\le e_v+\int_0^t[A_xz_x(s)+A_vz_v(s)+W]\,ds.
\end{aligned}
\tag{26}
$$

这一步不需要可微的最大范数，也不需要 $\dot d$。

### 4.3 有限时间定理

令非负非对角元的矩阵及驱动为

$$
M=\begin{pmatrix}0&1\\A_x&A_v\end{pmatrix},
\qquad b=\begin{pmatrix}D\\W\end{pmatrix},
\qquad e_0=\begin{pmatrix}e_x\\e_v\end{pmatrix}.
\tag{27}
$$

由逐次代入式(26)或标准正比较系统，

$$
\boxed{
\begin{pmatrix}|X(t)-Y(t)|_\infty\\|w(t)-V(t)|_\infty\end{pmatrix}
\le e^{Mt}e_0+\int_0^t e^{M(t-s)}b\,ds.
}
\tag{28}
$$

不等号按两个分量理解。矩阵指数与积分不要求 $M$ 可逆；$B$、$L$ 为零等退化情形亦适用。非恒定可积误差预算可保留在积分中；本轮使用整窗统一上界。

式(28)是真正有限尺寸、有限时间的条件误差界。它没有先取“小体极限”再把结论当有限器件。其对象是坐标能量质心和积分动量比；并非所有量子位置输出、单次实际轨迹或任意内部材料都已证明服从同一测地线。

## 5. 固定有限预算及解析可复算形式

事前固定

$$
\begin{gathered}
B=\frac1{100},\quad L=\frac1{50},\quad v=\frac12,\quad
r=\frac1{10000},\quad r_f=\frac1{10^7},\\
\sigma=\frac1{10000},\quad\eta=\frac1{100000},\quad T=1,
\qquad e_x=e_v=0.
\end{gathered}
\tag{29}
$$

直接代入得

$$
S=\frac52,\quad A_x=\frac3{16},\quad A_v=\frac{23}{80},
\quad D=\frac{11}{10^7},\quad W=\frac{39}{2000000}.
\tag{30}
$$

预定判据为第一分量 $<2\times10^{-5}$、第二分量 $<3\times10^{-5}$，不是按计算结果选择阈值。

有理复算可用

$$
Z_{12}=\sum_{n=0}^{12}\frac{M^nb}{(n+1)!}.
\tag{31}
$$

因为 $\|M\|_\infty=1$、$M,b$ 逐项非负、$e<3$，余项各分量满足

$$
0\le \int_0^1e^{M(1-s)}b\,ds-Z_{12}
\le \frac{3\|b\|_\infty}{14!}
\begin{pmatrix}1\\1\end{pmatrix}.
\tag{32}
$$

其中余项上界约为 $6.711\times10^{-16}$。精确有理运算给出

$$
\begin{aligned}
|X(1)-Y(1)|_\infty
&\le
\frac{30117024078568791167039466180833}
{2496186064167763968000000000000000000}
<1.207\times10^{-5}<2\times10^{-5},\\
|w(1)-V(1)|_\infty
&\le
\frac{359750129290938154233532896012643}
{15361145010263162880000000000000000000}
<2.342\times10^{-5}<3\times10^{-5}.
\end{aligned}
\tag{33}
$$

每个分数均为对应的 $Z_{12}$ 分量加式(32)尾界。代码与结果见[check.py](check.py)、[results.json](results.json)；独立审阅应分别核矩方程、常数和有理不等式。

这些数值只是有效合同预算，不证明所有前提已由某个实际共同量子过程同时实现；既不是实测数据，也没有认证原材料尺寸、正能源或测量精度。若材料/关系域不能支持此预算，应报告该实现失败，不将其升级为认知纲领失败。

## 6. 删除条件后的诊断意义

- 不控制 $\eta$，子体可以通过环境交换动量而加速；全体系总应力守恒对此没有禁止作用。
- 不控制 $\sigma$，积分压力与内部速度离散可改变式(16)，不能把任意完整应力直接当作 $Ew^\alpha w^\beta$。
- 不控制 $r$，相隔位置上的应力可以采样明显不同的连接；点上几何不足以替代整个体的加权响应。
- 不控制 $r_f$，即使净能源交换很小，远离质心的能源转移仍可能改变 $X'$。
- 不保同一图、同一时间和共同域，式(28)不能用于跨模型或各自正规坐标的拼接。
- 不保实际钟/CO读口匹配，式(28)仍只是一份有效应力模型内的结论。

这些是证明各误差项的角色定位。某个代数或坐标变换诊断不自动构成完整量子父模型或全公理反例。特别是非惯性坐标中的非零 $\Gamma$ 不等于真实曲率，也不是等效原理违例。

## 7. 接1192的范围和下一项停止线

本轮可以替换1192自由运动合同中的**一部分数学输入**：不再把应力质心的近测地运动独立假设，而以式(3)、(5)—(9)和有限时间误差推出。尚有四项独立匹配：

1. **原CO读口与应力质心。** [1117 §2.2及§7](../../archive_1106_1119/1117/proof.md)已区分实际CO位置、载体裸质心和完整仪器。即便另证读口均值与 $X$ 相差至多 $\epsilon_{\rm CO}$，也只得到读口均值距 $Y$ 至多 $\epsilon_{\rm CO}+Z_x$。它不直接给单次端点误差、完整位置分布、读后态或同一量子参考上的仪器距离；分布集中与端点统计需另证。
2. **实际内部钟与仿射固有时。** 本轮使用事前固定坐标时 $t$。式(21)给坐标时测地线，不能直接替代1192的原钟仿射参数合同；钟率/间隔误差和统一单位仍须独立匹配，不能通过重新参数化消去失败。
3. **共同实际历史与完整来源。** 多束材料需同处于一份历史并共享同一 $g$。所有本体、探针及其环境的来源须进入总账；额外暗部门的交换若未被 $\eta,r_f$ 控制，不能称这些体为本轮意义的近自由体。其对总弱来源下界的影响仍按1192另列，不能用本轮替代。
4. **量子量词。** 可把式(2)解释为指定态族中的有效期望应力，但须对该态族统一证明正 $E$、所有绝对矩及几何/交换合同。量子Hamiltonian正谱不自动给逐点经典能量条件或式(5)。本轮没有任意未知输入/参考 $R$ 的通道界，没有认证全Fock尾、仪器或完整量子几何。对多个态逐一成立的平均轨迹界，不等于完整过程的共同量子字典。

因此，本轮新增采用公理数为0；Einstein方程生成、实际钟束认证、暗物质/暗能量身份及整体目标交付均未取得。下一步须真正验证上述匹配在原实际过程是否成立，或由有限可观测任务拒绝某项合同；不以继续压小本组误差、调整任意候选参数或重做小体极限作为进展。
