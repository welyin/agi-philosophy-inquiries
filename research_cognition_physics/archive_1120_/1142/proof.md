# 1142证明：实际BT动力学的固定区域源能弱收支

2026-10-12。正式有限连接；验收状态见本轮独审及acceptance.json。继承1137原正常材料、b=0、固定实验室w、T=5/4及原参数，不新增认知公理、材料权限或科学校准。

本证明将1141的同一分辨率区域源能、独立波流及裸交互功放到实际H的同一连续账上。来源为[直接资格](../_shared/after1141_defect_direct_qualification.md)、[裸功资格](../_shared/after1141_BT_power_qualification.md)及[弱方法资格](../_shared/after1141_defect_weak_qualification.md)。正式稿放宽配对核中间常数，并补写外部弱Sylvester及混合图域论证；原资格保留其研究时点。

实际源子账的弱缺陷绝对期望积分小于 $2\times10^{11}$，同一裸交互功率期望绝对值积分小于 $10^{23}$。这些是原固定正常准备、有限分辨率和短窗的解析上界，不是实验数据。未迁移1141相对区域能源或相对流的单项旧精度，未给整个装置的局部能动量或能源仪器。

## 1. 固定字典与主恒等式

对象及约定逐字取前次缺陷筛查：

$$
S=\sqrt{m_0^2+P_f^2},\quad S_Y=\sqrt{m_0^2+(K-P_f)^2},\quad m_0=M+H_D+1,
$$
$$
\mathsf M=H_F+H_\phi+S+2M+H_Z+H_O+V,\quad
H=\sqrt{\mathsf M^2+K^2},\quad
H_o=\mathsf M-S+S_Y.
$$

E 始终是 1141 原**固定实验室** $E_{\epsilon,w}$ 连同输入共同作 LLP 变换的像；不是在 LLP 中重新固定一个相对区域。$[S_Y,E]=0$，而

$$
C_j=i[K_j,E]=i[P_{f,j},E]=E_{\epsilon,\partial_j w},
\quad C_o=i[H_o,E]=-\Phi_w+P_w,
\quad C_R=i[S,E],\quad g=H-\mathsf M.
$$

前次已独核的 Sylvester 式可以进一步精确重排为

$$
\boxed{
\mathcal D_w=C_R+
\int_0^\infty e^{-sH}
\left(\sum_j\{K_j,C_j\}-\{g,C_o\}-\{g,C_R\}\right)e^{-sH}\,ds .}
\tag{1}
$$

这里 $\mathcal D_w=i[H-H_o,E]$；积分先按弱形式解释。理由是将 $\mathsf M=H-g$ 代回[缺陷筛查](../_shared/after1141_balance_defect_screen.md)式 (6)，把 $\{H,C_R\}$ 的正半群积分精确积分掉。这一步消去原先需要先估 $\mathsf M e^{-sH}\Psi$ 的主交叉项，不是近重质量展开。下面证明右端在原实际向量上有可积支配，并解释弱极限；不能只凭式 (1) 的形式写法签收。

## 2. 一个统一的平滑源图界

以下极粗的常数只用于带 g 的小余项：对本轮的 $C_j$、$E_{\epsilon,\partial_l\partial_j w}$ 与 $\Phi_w$，均可取

$$
\|O\psi\|\le c\|H_\phi\psi\|+d\|\psi\|,
\qquad c=10^{27},\quad d=10^{30}.
\tag{2}
$$

这是独立核界，不从小态差推出。简要核证明如下。把四个平滑场写成 $a(f_{x,a})+a^\dagger(f_{x,a})$。旧核给

$$
\sum_a\|f_{x,a}\|^2=\nu_\epsilon<100\epsilon^{-4},\quad
\sum_a\|f_{x,a}/\sqrt\omega\|^2\le8\epsilon^{-3},\quad
\sum_a\|f_{x,a}/\omega\|^2<18\epsilon^{-2},
$$

其中 $\omega=|k|$；最后一个数是 1129 的原测试 φ 真空界。对权重 h，令 $V_h=\int|h|$。正常序二次型分成数项 $d\Gamma(A_h)$ 及配对项。三角与 Cauchy 给数核的单粒子图界 $\|A_h\omega^{-1}\|\le50V_h\epsilon^{-3}$、形式界 $\|\omega^{-1/2}A_h\omega^{-1/2}\|\le8V_h\epsilon^{-3}$。同粒子项用前者、异粒子项用后者，故它们提升为全部 Fock 上的同类 H_φ 图界，而不是把单粒子范数直接当全 Fock 有界性。

配对核 B 满足

$$
\|B/(\sqrt\omega\sqrt{\omega'})\|_2\le8V_h\epsilon^{-3},
\quad \|B\omega^{-1/2}\|_2\le30V_h\epsilon^{-7/2},
\quad\|B\|_2\le100V_h\epsilon^{-4}.
$$

双 annihilation 的加权 Cauchy 由 $H_\phi^2-d\Gamma(\omega^2)\le H_\phi^2$ 控制；双 creation 的 CCR 恒等式另带 $4d\Gamma(B^*B)+2\|B\|_2^2$。由 $\|H_\phi^{1/2}\psi\|\le(\|H_\phi\psi\|+\|\psi\|)/2$，得到总图系数可取

$$
c_h\le V_h(66\epsilon^{-3}+30\epsilon^{-7/2}),\quad
 d_h\le V_h(30\epsilon^{-7/2}+142\epsilon^{-4}).
$$

原 w 的一阶导不超过 16；由 $\|R''\|\le128$、径向 η 的 Hessian 界 128，二阶导不超过 $128+128+128<400$。因此 $V_h<400\cdot524<210000$。对 Φ 的 π/梯度混合项作实场极化，再安全付 12 倍，ε=10^-5 下仍严格小于式 (2)。这些中间系数从显示的松界直接推出：配对贡献使用 $2\kappa+\eta$ 与 $\eta+\sqrt2\beta$，不将更紧但未显示的半因子带入。这是刻意保守的界，无须优化常数。源动量与源能强对易，故这些图界在 $e^{-rS}$ 上保持指数权。

## 3. 实际 H 的二阶共同动量域可以补上

令 $\mathsf M_0=\mathsf M-V$、$H_0=\sqrt{\mathsf M_0^2+K^2}$。沿既有有界质量差引理，

$$H=H_0+W,\qquad \|W\|\le2v,\qquad v=16(11/7)C_*+558/7<10^{21}.$$

$H_0$ 与 $P_f$、K、H_φ 强对易。P_f 平移只改变 V；1138 已给 V 的共同一阶界 $C_Y$ 和 Laplacian 界 $C_{2,x}=9\times10^{21}$，并证明各混合二阶导存在。

固定 K 纤维及一个平移方向，记 A 为质量、F=√(A²+k²)。二次微分精确平方关系得到

$$
FF''+F''F=AA''+A''A+2(A'^2-F'^2).
\tag{3}
$$

原一阶 Sylvester 映射作用于任意有界自伴增量的范数至多为 2。正下界 $m_*=3M-v>2M$ 使最后一项的半群积分范数至多 $(\|A'\|^2+\|F'\|^2)/m_*$。逐固定 k 的有界差商先给 C² 范数可微，统一界再进入直积分；没有对无界 A 直接积分发散范数。三轴安全地取

$$
C_1=4C_Y<3\times10^{21},\qquad
C_2=2C_{2,x}+20C_Y^2/m_*<2\times10^{22}.
\tag{4}
$$

4 和 20 已付三轴列范数，未把单轴的 5 直接当向量常数。现在 W 是一个具有这些 C²(P_f) 界的**有界**扰动。因 H_0 群与 P_f 强对易，1138 的 D(P_f²) 图空间 Dyson/闭性证明可用于实际 H 本身。因此在原完整准备上

$$
\|P_f\Psi_t\|\le p_1+C_1t<10^{36},\quad
\|P_f^2\Psi_t\|\le p_2+2C_1p_1t+C_1^2t^2+C_2t<1.201\times10^{72}.
\tag{5}
$$

同样，对 $\Phi_s=e^{-sH}\Psi_t$ 有

$$
\|P_f^2\Phi_s\|\le e^{-m_*s}
[P_2+2C_1P_1s+C_1^2s^2+C_2s],
\quad P_1=10^{36},\ P_2=1.201\times10^{72}.
\tag{6}
$$

这里使用实际 H 的有界一、二阶共同力，不再误用只在 L 轨道证明过的结论。K² 与 H 强对易，原 $\|K^2\Psi_0\|=A_2<1.201\times10^{72}$ 原样保留。

## 4. 平台相消是一个对任意合法双向量成立的身份

令 W_t 为原完整自由正常相干源的 Weyl 位移，定义全局涨落能源

$$N_t=W_tH_\phi W_t^\dagger=H_\phi-L_t+E_{cl},\qquad E_{cl}\le1073A^2.$$

L_t 是原相干均值对应的线性 Segal 项。因为 $5+T+\epsilon<7$，w 支持中的自由均值精确为 π=A、∇φ=0。因此

$$
W_t^\dagger C_jW_t=E_{\epsilon,\partial_jw}+A\Pi_K(\partial_jw).
\tag{7}
$$

经典常数 $A^2\int\partial_jw/2$ 精确为零。式 (7) 是二次形式的 Weyl 共轭恒等式，**不要求被代入的实际/半群双向量是相干态**。前次筛查的真正缺口是涨落列的估计，不是平台常数在交叉型中不能消去。

对任意 u，记 $e(u)=\langle u,N_tu\rangle$。在 B_5（安全体积 600）内，四分量涨落场列满足

$$Z(u):=\|\mathbf1_{B_5}(\pi_K-\pi_{cl,K},\nabla\phi_K-\nabla\phi_{cl,K})u\|
\le2\sqrt{e(u)}+4\times10^{12}\|u\|.
\tag{8}
$$

它来自全空间 annihilation 列与有限区域 creation 真空项；保留正常序真空责任。由式 (7)，$\|\partial_jw\|_\infty\le16$、$\|\partial_jw\|_2<400$，

$$
|C_j(u,v)|\le8Z(u)Z(v)+400A\|u\|Z(v).
\tag{9}
$$

线性项是一个实际单场算符作用于 v；不必对它再使用一个会引入大混合图量的对称 Cauchy 界。

## 5. 不作重质量展开的实际低涨落证书

### 5.1 先处理同一自由 BT 加原 V

取 $\widetilde H=H_0+V$，仅作数学比较。既有质量差资格给

$$
\|(e^{-itH}-e^{-it\widetilde H})\Psi_0\|
\le\frac{tvA_2}{2m_*^2}<1.581\times10^{-68}.
\tag{10}
$$

K 与两生成元均强对易，所以此界保原 K/场纠缠，不换正常准备。

H_0 是自由强对易电荷 $Q=(E_f,P_f,K)$ 的函数

$$h(Q)=\sqrt{(E_f+\sqrt{m_0^2+P_f^2})^2+K^2},\quad E_f=H_F+H_\phi+2M+H_Z+H_O.$$

a(k) 始终指原物理源 annihilator 的共同 LLP 像；它在 LLP 中包含所需的 Y 相位，不是另选的相对场裸 annihilator。因此它同时改变 P_f 与 K、保持 K−P_f 不变，原 Y/源关联没有重置。其电荷增量是 $q=(|k|,k,k)$。自由物理 $H_{o,0}=H_o-V$ 在该增量下恰增加 |k|，因为 K−P_f 不变。令 $R_0=H_0-H_{o,0}$；在 Fock 光滑芯有精确有限差分身份

$$[R_0,a(k)]=[R_0(Q)-R_0(Q+q)]a(k).\tag{11}$$

方括号是自伴乘子。沿这个真实正自由谱段作标量均值定理，质量至少 3M，得到

$$
|\Delta R_0|\le\frac{|k|}{M}
\left(|P_f|+\frac23|K|+\frac53|k|\right).
\tag{12}
$$

这是平方根的全余项不等式，不是截断 1/M 展开。

令 $e(t)=\langle\Psi_{\widetilde H}(t),N_t\Psi_{\widetilde H}(t)\rangle$，$b_t(k)=a(k)-a_t(k)$。涨落能源导数中的自由差项只取
$\int\omega\overline{a_t}\langle\Delta R_0 a\rangle$ 的虚部。写 a=b_t+a_t 后，$|a_t|^2\langle\Delta R_0\rangle$ 精确为实数，因而整个经典项消失。加权 Cauchy 得

$$|e'_{R_0}(t)|\le2\sqrt{e(t)}\,r(t),\qquad
r(t)\le\frac AM[300(p(t)+k_0)+9000]<10^{-7}.\tag{13}$$

常数由原 β_0 的 $\|\nabla\beta_0\|_2/\sqrt2<300$、$\|\Delta\beta_0\|_2/\sqrt2<5000$ 给出；$p(t)\le p_1+C_Yt$、$k_0=A_1<10^{36}$。未冻结支撑或切掉源的任何正常尾。

V 的同一源导数给

$$|e'_V(t)|\le\frac DA\sqrt{2e(t)+600\nu_\epsilon},\qquad D=1.2\times10^{21}.$$

这里使用未卷积力括号的原 L² 列界，以及 Y 纤维中实际 π_ε 涨落列；全部保持、接收和继续尾均包含。若 e≤1，右端连同式 (13) 小于 0.295。初态 e(0)=0，因此第一次越界法给

$$\sup_{0\le t\le5/4}e_{\widetilde H}(t)<0.4.\tag{14}$$

上述导数按时间弱形式理解。H_0 与 H_φ 强对易，V 是有界且源实谱导数有统一列界的扰动；先对线性 Segal 测试和有限自由电荷芯证明式 (11)，再以式 (12) 的加权支配、原 P_f 一阶域和 N_t 形式域延拓。全过程 H_φ 图界独立成立，故有限差分积分、时间弱导数和第一次越界极限均有支配；没有假定实际 H² 或更高自由能源域。

### 5.2 移回实际 H：只在极小差处使用大图界

实际 H 已有 $\|H_\phi\Psi_H(t)\|<B_*=1.074\times10^{71}$。比较 $\widetilde H$ 的相应图界由 $H_0-3M$ 的正自由项和 $\widetilde H$ 自身守恒、加 2v 另证，同样小于 B_*，不能从式 (10) 推出。

原源 $\|\omega a_t\|\le94A$，因此对两条过程分别有

$$
\|N_t\Psi_\bullet(t)\|\le B_*+2(33A)(33A)+1073A^2+94A
<B_c:=4.4\times10^{71}.
\tag{15}
$$

现在式 (10) 可以合法地乘**两侧已证明的** N_t 图界：

$$e_H(t)<0.4+2(1.581\times10^{-68})B_c<1.5\times10^4.\tag{16}$$

这一步没有把原 $10^{-8}$ 态差乘区域能源。使用的是交互专属 $10^{-68}$ 比较及独立涨落能源图界；其放粗误差约 $1.4\times10^4$，足以服务下述 1/M 缺陷，不声称直接保持 1141 的相对区域能源单项精度。

## 6. 正半群所需的列界

自由 S 的源电荷差不超过 |k|。由与式 (11) 同样的精确 annihilator 身份，

$$
\sqrt{e(e^{-rS}\Psi_t)}\le e^{-Mr}
[\sqrt{e_H(t)}+300Ar].
\tag{17}
$$

对 H_0，电荷增量的绝对值不超过 3|k|，且 H_0(Q) 与 H_0(Q+q) 均至少为 3M；不需要另假定该增量非负。故自由正半群相应界为 $e^{-3Ms}[\sqrt e+900As]$。实际 $H-H_0=W$ 有界≤2v；Duhamel 给

$$\|(e^{-sH}-e^{-sH_0})\Psi_t\|\le2vs e^{-m_*s}.$$

两侧的 N_t 图界分别由正自由图域与实际 $G=H-3M+v$ 守恒给出 $e^{-m_*s}B_c$。对差向量作正 N_t 的 Cauchy，得到

$$
\sqrt{e(e^{-sH}\Psi_t)}\le e^{-m_*s}
[\sqrt{e_H(t)}+900As+2\sqrt{vB_cs}].
\tag{18}
$$

因此相应局部场列安全满足

$$
Z(e^{-rS}\Psi_t)\le e^{-Mr}[5\times10^{12}+600Ar],
$$
$$
Z(\Phi_s)\le e^{-m_*s}[5\times10^{12}+1800As+10^{47}\sqrt s].
\tag{19}
$$

后者用了 $4\sqrt{vB_c}<10^{47}$。这些界明确支付了半群对相干状态的改变，没有假定半群仍保持平台态。

另一方面，原实际 $P_f^2$ 及 H_φ 图界、二者强对易给

$$\|H_\phi^{1/2}P_f\Psi_t\|^2\le\|P_f^2\Psi_t\|\,\|H_\phi\Psi_t\|.$$

连同 $\sqrt{E_{cl}}<33A$ 和真空列，得到

$$Z(P_{f,j}e^{-rS}\Psi_t)<1.5\times10^{72}e^{-Mr}.$$

对外 K 方向，用 $H_\phi\le G$ 的**形式序**（不将其误作算符范数序），K/G 强对易与守恒给

$$\|\sqrt G K\Psi_t\|^2=\langle K^2\Psi_0,G\Psi_0\rangle
\le A_2\|G\Psi_0\|.$$

故同样有

$$Z(K_j\Phi_s)<1.5\times10^{72}e^{-m_*s}.\tag{20}$$

这些是当前实际及半群交叉向量的独立混合形式界，正是此前缺少的对象。

## 7. 四项缺陷的固定预算

置 $L_c=1.6\times10^{73}$。式 (9)、原 $\|P_f\Psi_t\|,\|K\Psi_t\|<10^{36}$ 及式 (20)使每个 C_j 交叉项的系数至多为 L_c。

### 7.1 内部反冲

用原精确 S-Sylvester 积分与式 (19)，三轴及反对易付 6 倍，

$$
|C_R(\Psi_t,\Psi_t)|
\le6L_c\left[\frac{5\times10^{12}}{2M}+\frac{600A}{4M^2}\right]
<3\times10^6.\tag{21}
$$

### 7.2 外部 K 项

$2m_*>4M$；$\Gamma(3/2)<1$。因此

$$
\left|\int_0^\infty\sum_j\{K_j,C_j\}(\Phi_s,\Phi_s)ds\right|
\le6L_c\left[\frac{5\times10^{12}}{4M}
+\frac{1800A}{16M^2}+\frac{10^{47}}{8M^{3/2}}\right]
<2\times10^6.\tag{22}
$$

### 7.3 gC_o 项

$g\le K^2/(2m_*)$，故 $\|g\Phi_s\|\le g_0e^{-m_*s}$，可取 $g_0=10^{-8}$。式 (2)、原源图界及已核 P 图界给

$$\|C_o\Phi_s\|,\ \|C_j\Phi_s\|,
\ \|E_{\epsilon,\partial_l\partial_jw}\Phi_s\|
\le N_c e^{-m_*s},\qquad N_c=2\times10^{98}.$$

所以

$$\int_0^\infty|\{g,C_o\}(\Phi_s,\Phi_s)|ds
\le\frac{2g_0N_c}{4M}<2\times10^{10}.\tag{23}$$

### 7.4 gC_R 项：把动量移到左侧

不能在这里用 $H_\phi g\Phi_s$ 的 1/s 粗界，那会产生假的对数发散。在内部 S 积分中使用

$$\{P_{f,j},C_j\}=2P_{f,j}C_j-[P_{f,j},C_j],\qquad
[P_{f,j},C_j]=-iE_{\epsilon,\partial_j^2w}.$$

这使全部源图算符作用于右侧 Φ_s，而动量作用于左侧 gΦ_s。又因 g≤|K|，且 $[P_f,g]$ 的共同范数≤5C_Y，

$$\|P_fg\Phi_s\|
\le\|KP_f\Phi_s\|+5C_Y\|\Phi_s\|.$$

K 与 P_f 强对易，故 $\|KP_f\Phi_s\|^2\le\|K^2\Phi_s\|\,\|P_f^2\Phi_s\|$。式 (6)于是给

$$\|P_fg\Phi_s\|
\le e^{-m_*s}[2\times10^{72}+10^{65}\sqrt s+4\times10^{57}s].\tag{24}$$

只用了已证明的实际二阶动量图域，没有新增混合高阶能源输入。由 S≥M 的积分及式 (2)，

$$
2\int_0^\infty|C_R(g\Phi_s,\Phi_s)|ds
\le\frac{6N_c}{M}\left[
\frac{2\times10^{72}}{4M}+\frac{10^{65}}{8M^{3/2}}+\frac{4\times10^{57}}{16M^2}\right]
+\frac{3N_cg_0}{4M^2}<7\times10^{10}.\tag{25}
$$

因此四项之和严格小于 $10^{11}$，原时间窗给有限证书

$$
\boxed{\int_0^{5/4}|\langle\mathcal D_w\rangle_{H,t}|dt<2\times10^{11}.}
\tag{26}
$$

## 8. 弱形式身份、正则化及不能省略的范围

式 (2)给 E 及源图量的合法域；实际 D(H) 与原自由/比较共同域等价，原轨道保该域。标量 $\langle E\rangle_{H,t}$ 的导数可以按 $i\langle H\Psi,E\Psi\rangle-i\langle E\Psi,H\Psi\rangle$ 的弱形式定义，不预设 HEΨ 或 EHΨ。

内部 S 的 Sylvester 身份先对 S 的谱截断双向量证明。S 的谱截断与 H_φ、P_f 强对易，保平方关系和所需源图域；它不是对 E 的朴素谱截断。式 (21)/(25) 的可积图支配随后识别其极限为同一弱交换子 $i[S,E]$。

外部弱身份可直接证明，避免先写未经定义的强乘积。令

$$
C_g=i[g,E],\qquad C_{\mathsf M}=i[\mathsf M,E]=C_o+C_R.
$$

E、g、S都为H图有界：H=H_0+W且W有界，H_0是共同正自由电荷的函数，E的H_φ图界由式(2)给出。因此上列交换子及C_R均为D(H)上的连续双形式，定义为 $i(\langle Au,Ev\rangle-\langle Eu,Av\rangle)$，无须假定HEu存在。

H、g、$\mathsf M$和K强对易，且 $2\mathsf M g+g^2=K^2$。在D(H²)双向量上按上述弱积展开，E的对称性取消中间项，得到

$$
\{H,C_g\}=\sum_j\{K_j,C_j\}-\{g,C_{\mathsf M}\}.
$$

其中 $\mathsf M u$、gu、Ku属于D(H)，所以E只作用于已经合法的向量。对 $u_s=e^{-sH}u$ 微分连续形式 $C_g(u_s,u_s)$，其负导数就是右侧；正质量下界及H图界使无穷端边界为0。积分后用 $\mathcal D_w=C_R+C_g$ 得到式(1)。全过程不用HEu或EH²u强积。

为使内部S积分与此弱展开在同一密集域上相遇，取图Banach空间 $\mathcal X=D(P_f^2)$，范数为 $\|u\|+\|P_f^2u\|$。自由H_0与P_f强对易，其正半群在X上解析；W由第3节C²(P_f)界在X上有界。有界扰动后的解析半群在X上等于Hilbert空间中e^{-sH}的限制，故原u∈X在s>0时的各阶H半群导数仍在X。这补足Pf与H混合图域；不是仅从Hilbert空间H平滑性推断这些混合域。内部S谱截断仍与H_φ、P_f强对易，源核图支配及第7节右作用重排分别负责移除截断。

先用 $e^{-\delta H}\Psi_t$ 检验，再移除δ；不对E作会产生未交代[V,截断]项的朴素谱截断。C_o始终是源方程已独立给出并延拓的 $-\Phi+P$。

δ→0 时，H 图收敛给源图/普通场形式收敛；式 (5)/(6)控制动量图。对 $0<\delta\le1/M$，式 (18)给额外涨落列仍远小于 $5\times10^{12}$，式 (19)—(25)可用相同安全常数支配。因此可移除数学正则化，识别式 (1)右端为原实际轨道的同一弱缺陷。没有遗漏对 E 截断的 [V,截断] 项，因为没有作那种截断。

由此实际来源子账为

$$
\langle E\rangle_H(t)-\langle E\rangle_H(0)
=-\int_0^t\langle\Phi_w\rangle_Hds
+\int_0^t\langle P_w\rangle_Hds
+\int_0^t\langle\mathcal D_w\rangle_Hds.
$$

Φ 由原源场独立定义，P 由原 V 的源导数独立定义，D 是经过估计的实际差；没有用总 H 的交换子反向造流。原 BT 裸功率资格可以接入；其迁移误差 $4\times10^{15}$ 大于本稿缺陷预算，但仍远小于原交互功的安全上界。

**不随式 (26)自动成立的事项：** 实际 H 的相对区域能/相对流单项分别具有 1141 的同样上界；源未分辨能源和其余材料的完整局部账；实际能源仪器；任意 boost、任意区域或其他正常输入的同一误差；全物质的普遍局域性。本资格只处理原固定区域、原短窗、原有限材料合同。

## 9. 与同一裸交互量资格合并

[已核裸功资格](../_shared/after1141_BT_power_qualification.md)对同一P_w建立完整源图界，两过程各自的能源资源独立受控。未分辨源及全部量子接收尾仍在实际态中。原固定短窗结论为

$$
\int_0^{5/4}|\langle P_w\rangle_H-\langle P_w\rangle_o|\,dt
<3.526\times10^{15}<4\times10^{15}.
$$

1141的未放粗预算小于 $6.958\times10^{22}$；相加仍有余量，给

$$
\int_0^{5/4}|\langle P_w\rangle_H|\,dt<10^{23}.
$$

这里是功率期望绝对值的积分，非 $\langle|P_w|\rangle$，不是正功下界。联合式(26)，实际同一对象满足

$$
\left|\langle E\rangle_H(t)-\langle E\rangle_H(0)
+\int_0^t\langle\Phi_w\rangle_Hds-\int_0^t\langle P_w\rangle_Hds\right|
<2\times10^{11}.
$$

有界时间测试也可由缺陷L¹界控制；本轮不以时间平滑或另换仪器掩盖原误差。

## 10. 净增量、解释与下一方向

1141只有H_o来源子账，裸P资格只有一个无界量的同名迁移。本轮新增实际H的低涨落资源、实际二阶动量域、交叉半群列及同一连续弱收支身份；标准Fock、Sylvester和Weyl工具本身不另计物理发现。

该连接说明在已声明的短窗和材料上，BT改变动力学没有破坏原源能—流—交换功的有限精度核对。它没有证明现行认知公理唯一选择本模型，亦未消除既有质量、核、正常准备和有效读取等物理输入。

下一项先核源子账能否与剩余材料、支撑及交互储能的同一有限切面账相接，并明确哪些现有总四P及旧规范／引力条件定理前提真正得到同一对象。不改A/M/ε、不施工全部BT应力或1140汇集网络。有限域失败只限相应对象/合同，不能放大为整个认知纲领的反证。

[代码](check.py) · [结果](results.json) · [选择审计](selection_audit.md) · [下一步](NEXT.md)。解析证明承担连续结论，复算只核有理预算、冻结来源及有限代数诊断，不是连续物理模拟。
