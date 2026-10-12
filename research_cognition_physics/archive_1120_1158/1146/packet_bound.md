# 1146：同一正常紧支运动包的钟速与标量操作——完整两谱有限误差界

2026-10-12。本文是1146的包级证明，配合本轮主证明的完整CAR、正常量子来源与实际读取。旧941—943/953的固定动量或世界线公式、1126的紧支Dirac传播、1129的有限核与完整来源比较直接复用；新增的是**同一正常紧支运动材料、同一未知内部态和同一实际移动读口**同时给出钟系数与标量操作系数的有限证书。没有把一粒子概率当成完整场论证书。

## 1. 模型、真实模式和两项对照

自然单位。旋量与两味内部空间的次序为 $\mathbb C^4_{\rm spin}\otimes\mathbb C^2_{\rm flavor}$。标准Dirac矩阵满足Clifford关系。取

$$
P_1=|1\rangle\langle1|,\quad X=\sigma_x,\quad
m=10^{10},\quad \delta=10^{-3},\quad
p_0=\frac{3m}{4}e_z,\quad E_0=\frac{5m}{4},\quad
r_0=\frac45,\quad v_0=\frac35,\quad T=1 .
\tag{1}
$$

同一自由材料的两个味质量是 $m,m+\delta$。质量分裂、Dirac材料、标量型或时间分量矢量型交互、部署切片与有限分辨率都是本候选的物理输入，并非新增认知公理。

严格采用1129的实际归一包：

$$
q(s)=e^{-1/(1-s^2)}1_{|s|<1},\qquad
a(x)=\frac{q(x_1)q(x_2)q(x_3)}{\|q\|_2^3},\qquad
f(x)=e^{ip_0\cdot x}a(x)u_0 ,
\tag{2}
$$

其中 $\|u_0\|=1$，$(\alpha\cdot p_0+\beta m)u_0=E_0u_0$。因此 $\|f\|=1$，$f\in C_c^\infty(B_2,\mathbb C^4)$。这里只在载波中心选正能旋量；**真实Cauchy模式 $f$ 未作正能投影，保有完整负频尾部**。

令 $\kappa_1=\|\nabla a\|_2$、$\kappa_2=\|\Delta a\|_2$。第5节给解析界

$$
\kappa_1<20,\qquad \kappa_2<600<1000 .
\tag{3}
$$

门函数逐字复用1129：

$$
b(s)=e^{-1/[s(1-s)]}1_{0<s<1},\quad I_b=\int_0^1 b(s)\,ds,\quad
F(t)=\frac{2b(2t-1/2)}{I_b},\quad \lambda=\frac{\pi}{2}.
\tag{4}
$$

$F\ge0$、支撑于 $[1/4,3/4]$、$\int_0^1F=1$、$\|F\|_\infty\le16$、$\|F'\|_\infty\le512$。它先单调增、后单调减，所以

$$
\int_0^1|F'|\le32,\qquad \int_0^1F^2\le16 .
\tag{5}
$$

先定义两个仅供比较的空间均匀Hamiltonian，完整作用在 $L^2(\mathbb R^3,\mathbb C^4)\otimes\mathbb C^2$：

$$
H_s(t)=\alpha\cdot p+\beta\otimes\big(mI+\delta P_1+\lambda F(t)X\big),
$$
$$
H_v(t)=(\alpha\cdot p+\beta m)\otimes I+
       \beta\otimes\delta P_1+I_{\rm spin}\otimes\lambda F(t)X .
\tag{6}
$$

下标 $s$ 是本轮欲连接的标量质量型；$v$ 是改变真实交互物理类型的对照，不能据此排除真正的矢量理论。二者质量、初态、读口、来源门函数和有限空间核均相同；不把后者称为同一个标量材料的“另一种坐标”。

共同实际目标模式为

$$
g_t(x)=e^{ip_0\cdot x-iE_0t}a(x-v_0t e_z)u_0 .
\tag{7}
$$

$\|g_t\|=1$，支撑于 $B_{2+v_0t}$。两个内部比较酉分别满足

$$
i\dot R_s(t)=r_0\big(\delta P_1+\lambda F(t)X\big)R_s(t),\qquad
i\dot R_v(t)=\big(r_0\delta P_1+\lambda F(t)X\big)R_v(t),\qquad R_j(0)=I .
\tag{8}
$$

**两支的真实钟质量分裂共同乘 $r_0$，只有源型不同；目标读口始终是相同的 $g_t\otimes e_0,g_t\otimes e_1$。** 不随味质量或比较对象重选模式。这两列严格正交，任意 $g_t\otimes R_j(t)e_i$ 也严格正交。

## 2. 完整两谱有隙比较：不删负能块、不假定内操作可交换

先证一个可复用的算符界。固定动量 $p$，写

$$
h_m=\alpha\cdot p+\beta m,\quad E_p=\sqrt{|p|^2+m^2},\quad
P_\pm=\frac12(I\pm h_m/E_p),\quad Q=P_+-P_- .
$$

Clifford恒等式给

$$
P_\pm\beta P_\pm=\pm r(p)P_\pm,\qquad r(p)=m/E_p .
\tag{9}
$$

对一般绝对连续Hermitian内部矩阵 $V_s(t),V_v(t)$，令

$$
H=h_m+\beta\otimes V_s+I\otimes V_v=D+O,
$$
$$
D=h_m+D_1,\quad D_1=rQ\otimes V_s+I\otimes V_v,\quad
O=(\beta-rQ)\otimes V_s,\quad
\mu(t)=\|V_s(t)\|+\|V_v(t)\| .
$$

$O$ 是纯两谱交叉块。$\|D_1\|\le\mu$、$\|O\|\le\|V_s\|\le\mu$。定义反自伴算符

$$
S(t)=\frac{P_-O(t)P_+-P_+O(t)P_-}{2E_p}.
\tag{10}
$$

有 $[h_m,S]=-O$，$\|S\|\le\|V_s\|/(2E_p)$、$\|\dot S\|\le\|\dot V_s\|/(2E_p)$。令 $U,U_D$ 是 $H,D$ 的精确传播；直接微分得到

$$
\frac{d}{dt}\{U(T,t)S(t)U_D(t,0)\}
 =
U(T,t)\{\dot S-iO+i[D_1,S]+iOS\}U_D(t,0).
\tag{11}
$$

将Duhamel式 $U(T,0)-U_D(T,0)=-i\int U(T,t)OU_D(t,0)dt$ 代入并分部积分；用所有传播均为酉，得到完整算符范数界

$$
\|U-U_D\|\le
\frac{\|V_s(0)\|+\|V_s(T)\|+
\int_0^T\|\dot V_s\|\,dt+3\int_0^T\mu(t)^2\,dt}{2E_p}.
\tag{12}
$$

常数3来自 $2\|D_1\|\|S\|+\|O\|\|S\|$。这是有边界项的精确恒等式所给上界，不是仅取有效Hamiltonian一阶项；不要求 $[V_s(t),V_s(t')]=0$，也不要求 $\delta P_1$ 与 $X$ 对易。若 $V_v\ne0$ 而 $V_s$ 固定，式(12)可以收紧；本轮为保持两支统一只用上述安全粗界。

式(6)两种情况分别为 $(V_s,V_v)=(\delta P_1+\lambda FX,0)$ 和 $(\delta P_1,\lambda FX)$。二者同取

$$
\mu(t)=\delta+\lambda F(t),\quad
L(t)=\int_0^t\mu(u)du,\quad
\mathcal B(t)=\frac{2\delta+\lambda\int_0^t|F'(u)|du+
3\int_0^t\mu(u)^2du}{2}
\tag{13}
$$

即可上界相应的分子一半。这里式(13)仅在 $F(0)=F(t)=0$ 的终点使用（本轮 $t=1/4,1$），避免遗漏一般终点的源边界项。因 $E_p\ge m$，两种完整传播与其两谱对角比较的差均不超过 $\mathcal B(t)/m$。

## 3. 从每个动量到同一真实紧支移动包

四维Clifford向量 $(p,m)/E_p$ 给投影Lipschitz界；沿动量线段积分：

$$
\|P_+(p)-P_+(p_0)\|\le\frac{|p-p_0|}{2m},\qquad
\|P_- f\|_2\le\frac{\kappa_1}{2m}.
\tag{14}
$$

把完整两谱对角传播与只写作比较的正带内酉输出相减，负带项以其**真实保留的范数**计入，代价至多 $2\|P_-f\|\le\kappa_1/m$；没有把它设为0。

正带内演化的生成元分别为 $r(p)(\delta P_1+\lambda FX)$ 和 $r(p)\delta P_1+\lambda FX$。$\|\nabla r\|\le1/m$；再一次精确Duhamel给正带内部酉到式(8)的包范数误差不超过 $L(t)\kappa_1/m$。矢量支只需 $\delta t\kappa_1/m$，但同一粗界仍有效。

自由色散满足 $\|\nabla^2 E_p\|\le1/m$，因而

$$
|E_p-E_0-v_0(p_z-p_{0z})|\le\frac{|p-p_0|^2}{2m}.
\tag{15}
$$

Plancherel将它变成 $t\kappa_2/(2m)$。式(7)恰是线性色散对应的实空间平移；载波的空间相位已保留，未把不同位置上的相位比较成一个错误钟。

综上，对 $j=s,v$ 以及任意内部态和任意被动参考 $R$，有**同一输入码空间上的完整算符／等距界**

$$
\boxed{
\left\|U_j(t)(f\otimes I)-g_t\otimes R_j(t)\right\|
\le d_{\rm packet}(t):=
\frac{\mathcal B(t)+(1+L(t))\kappa_1+t\kappa_2/2}{m}.
}
\tag{16}
$$

张量 $I_R$ 不改变该范数。式(16)包括全部负带、离模式与相位误差，而非只比较一次测量均值。它仍是单Cauchy空间定理，必须经本轮主证明的完整CAR实施与正常源比较，才成为整个物理输出的Q/R证书。

## 4. 相同有限空间核：明确计算载波成本

复用1129的实、偶、非负径向归一平滑核 $\rho_\epsilon$，$K_\epsilon=\rho_\epsilon*$；$\|K_\epsilon\|\le1$、$\|(K_\epsilon-I)v\|\le\epsilon\|\nabla v\|$。固定

$$
\epsilon=10^{-16},\qquad
0\le\eta\le1,\quad \eta|_{B_4}=1,\quad
\operatorname{supp}\eta\subset B_5 .
$$

实际外源比较一粒子Hamiltonian为

$$
H_{s,\epsilon}=h_m\otimes I+\beta\otimes\delta P_1+
\lambda F(t)(K_\epsilon\eta K_\epsilon)\otimes\beta X,
$$
$$
H_{v,\epsilon}=h_m\otimes I+\beta\otimes\delta P_1+
\lambda F(t)(K_\epsilon\eta K_\epsilon)\otimes I_{\rm spin}X.
\tag{17}
$$

张量写法中 $K_\epsilon\eta K_\epsilon$ 只作用空间、$\beta$ 只作用旋量、$X$ 只作用味；首项同样按此扩张。

均匀完整Dirac传播 $U_j(t)f$ 的支撑在 $B_{2+t}$。因此 $K_\epsilon U_j(t)f$ 落在 $B_{3+\epsilon}\subset B_4$，平台函数在这条**比较轨道**上恰为1。$K_\epsilon$ 和梯度与均匀 $U_j(t)$ 对易。故

$$
\|(K_\epsilon\eta K_\epsilon-I)U_j(t)(f\otimes I)\|
\le2\epsilon\|\nabla f\|
\le2\epsilon(|p_0|+\kappa_1).
$$

沿均匀比较轨道使用Duhamel即可得到

$$
\boxed{
d_{\rm kernel}(t)\le
2\lambda\epsilon(|p_0|+\kappa_1)\int_0^tF(u)du .
}
\tag{18}
$$

这里不能只用旧零载波的 $\kappa_1$；$|p_0|$ 的巨大但有限资源已计入。对实际有核传播不假设其严格光锥支撑；有限核是声明的有效分辨率，不冒称点态局域场密度或原点态Ward壳流。

## 5. 解析包常数与固定有限预算

旧1129已证 $\|q\|_2^2\ge e^{-8/3}>1/16$、$|q'|<2$，给 $\kappa_1^2<384<400$。令 $y=1-s^2$，则

$$
q''(s)=\left(-\frac2{y^2}-\frac{8s^2}{y^3}+\frac{4s^2}{y^4}\right)e^{-1/y}.
\tag{19}
$$

取 $z=1/y\ge1$，用 $\sup_{z\ge0}z^ne^{-z}=n^ne^{-n}$ 及 $e>8/3$，

$$
|q''|\le2(2/e)^2+8(3/e)^3+4(4/e)^4
<\frac{2097}{64}<33 .
$$

所以 $\|q''\|_2<33\sqrt2<50$、$\|q''\|_2/\|q\|_2<200$。对归一乘积包以三角不等式估计三个二阶导数，得 $\kappa_2<600$。以下仍使用较松的1000，数值积分不承担此证书。

取 $\lambda<11/7$、$\int F^2\le16$、$\int|F'|\le32$。终时 $T=1$：

$$
\mathcal B(1)\le
\delta+16\lambda+
\frac32(\delta^2+2\delta\lambda+16\lambda^2),
\quad L(1)\le\delta+\lambda.
$$

完整精确有理预算为

$$
d_{\rm packet}(1)<6.358625\times10^{-8},\qquad
d_{\rm kernel}(1)<2.357142864\times10^{-6},
$$
$$
\boxed{d(1):=d_{\rm packet}(1)+d_{\rm kernel}(1)
<2.420729109\times10^{-6}<3\times10^{-6}.}
\tag{20}
$$

不存在对 $\epsilon\to0$、$m\to\infty$ 的物理假设；这是固定有限 $m,\epsilon$ 下的误差结论。正常KG来源所用 $A=10^{100}$ 及其全部有限能源、后态、反作用由主证明登记。本包定理不将它们视为免费背景。

实际自由海为各味 $S_i=1_{\alpha p+\beta(m+\delta_i)<0}$，$\delta_i=0,\delta$。沿质量方向同样用投影Lipschitz估计，

$$
\|S_i f\|\le\frac{\kappa_1+\delta_i}{2m}
<1.00005\times10^{-9}.
\tag{21}
$$

因此任何有限个正交味／rail实际模式的真空占据都很小但非零。若共有 $N$ 个，联合空模投影 $e_0=\prod(1-n_i)$ 满足

$$
\|e_0\Omega\|^2\ge1-\sum_i\langle\Omega,n_i\Omega\rangle
\ge1-N\left(\frac{20+\delta}{2m}\right)^2>0
$$

（本轮只需两个味模式）。可以按真实范数精确归一；不能沿用旧1129另一旋量的 $1/4$，也不能称局部模式恰为空真空。

## 6. 同一材料的两个实际有限读窗

### 6.1 源脉冲开始前的实际钟

取 $t_c=1/4$，同一初始味态 $|+\rangle$，终止时在共同移动模式 $g_{t_c}$ 的味空间读取 $E_{-Y}=(I-\sigma_y)/2$。这是可选终止实验，未声称先读钟而后续操作不受扰。

因为 $F|_{[0,t_c]}=0$，外源比较中两支都满足

$$
R(t_c)=e^{-ir_0\delta t_c P_1},\qquad
p_c=\frac{1+\sin(r_0\delta t_c)}2.
\tag{22}
$$

此窗 $L=\delta t_c$、$\mathcal B=\delta+\frac32t_c\delta^2$、$d_{\rm kernel}=0$，给

$$
d_{\rm packet}(t_c)
\le\frac{\delta+\frac32t_c\delta^2+(1+\delta t_c)20+500t_c}{10^{10}}
=1.45006000375\times10^{-8}.
\tag{23}
$$

错误实验室钟字典在同一读口预测 $(1+\sin(\delta t_c))/2$。因 $\sin x\ge x-x^3/6$、$\sin y\le y$，

$$
\Delta p_{\rm clock}^{\rm ideal}
\ge\frac{\delta t_c-r_0\delta t_c-(\delta t_c)^3/6}{2}
>2.4999998697\times10^{-5}.
\tag{24}
$$

这个差须减主证明给的完整CAR与量子源误差，不能直接领取整场结论。真实量子来源在早窗并非精确关闭；式(23)针对外源比较，来源涨落尾由同一全输出Duhamel界另计。该早窗比用终时总界有意义：实际钟信号没有被后来源型操作的误差吞没。

### 6.2 同一钟下，错误源类型仍在终末实际读口可辨

终时 $T=1$ 初始味 $|0\rangle$，读取同一 $g_1$ 上的 $P_1$。在 $\delta=0$ 时，式(8)给

$$
p_s^{(0)}=\sin^2(2\pi/5),\qquad p_v^{(0)}=1,\qquad
p_v^{(0)}-p_s^{(0)}=\cos^2(2\pi/5)=\frac{3-\sqrt5}{8}.
\tag{25}
$$

真实 $\delta>0$ 与 $X$ 不对易，未把该项删去。每支由Duhamel离 $\delta=0$ 的酉差至多 $r_0\delta T\le\delta$，故各支任一效果概率差不超过 $\delta$。于是同一实际内钟下的有效源型分离至少

$$
\Delta p_{\rm source}^{\rm ideal}\ge\frac{3-\sqrt5}{8}-2\delta
>0.0934915.
\tag{26}
$$

再减两支实际完整输出误差，才是本轮主证明的最终合同；单粒子包层可先各减式(20)。该分离只拒收“把 $I_{\rm spin}X$ 改叫 $\beta X$ 而不变预测”的配方。它不排除矢量物质，不宣称所有源必须是标量。

## 7. 复算内容与明确边界

同目录 packet_check.py 给三组独立检查：

1. 精确有理预算、解析包常数链及归一bump的独立积分诊断；
2. 可精确求解的连续旋转非交换内部源，在完整 $8\times8$ Dirac矩阵中核两谱分部积分结构与传播上界，既测标量也测矢量；
3. 同一实际味效果的早期钟与终末源型读取，独立核矩阵读口符号、严格分离预算及保留质量分裂的平滑脉冲诊断。第2组同时包含固定偏移动量的两谱投影与色散余项检查。

矩阵诊断不模拟三维连续场；积分诊断也不替代解析支撑或严格常数。完整矩阵诊断采用适中质量，实际 $m=10^{10}$ 的微小质量分裂与相位由有理／解析界负责，不用两个巨大浮点能源相减认证精度。最终预算使用有理不等式而不是浮点观察。全部资源是固定参数，未作优化扫描。

本证书领取同一实际移动模式上的共同系数 $r_0\delta$ 与 $r_0\lambda$，以及它们同矢量源的可辨差。它不领取：原1144独立双站KG闭环、实际定位架与准备读取的自治完备、点态壳通量、精确Lorentz协变、共同总四动量闭合、或全部认知公理已经唯一选出Dirac/Yukawa材料。正常源与完整CAR接口须与本轮主证明共同读取。
