# 1299：受保护 Dirac 质量、低维多体交换与整数超荷

2026-10-18。作者证明；状态由后续验收记录确定。本轮按[事前准入](../_shared/round1299_admission.md)给结构定理与声明阶四点矩阵元，不把其升级为实际制备／接收合同。本目录仅一组有限代数及归一核验；全整数、全低维结论由下述解析证明承担。

## 1. 同一比较类与旧结果

直接继承[1035解析稿 §§1–5](../../archive_1009_1043/1035/drafts/joint_selection_derivation.md)：三代左手 Weyl 字段 $Q,u^c,d^c,L,e^c,N$，单复 Higgs $H$，固定

$$
G_6=(SU(3)\times SU(2)\times U(1))/\mathbb Z_6,
\qquad Y_k=Y_0+3k(B-L),\quad k\in\mathbb Z.
\tag{1}
$$

全部左手 Weyl 在 $r$ 下乘 $i$，$H$ 乘 $-1$，$r^2=(-1)^F$。[1060证明 §§2–5](../../archive_1046_1062/1060/proof.md)已给每个整数 $k$ 的精确全局群嵌入及每代 16 的反常拉回；本轮沿用其声明的 $\mathrm{Spin}^{\mathbb Z_4}\times G_6$ 背景范围。它不等于全部相互作用量子场论或全部边界仪器的构造。

本轮局域顶点是 Lorentz、$G_6$ 和 $r$ 不变的多项式，正则维数 $d\le6$；允许共轭字段、任意分配的协变导数和规范场强，固定系数不携带规范／离散荷。没有新增带荷字段、反 Higgs 幂、外部会变换的系数。$d$ 是字段和导数的正则维数，不能用系数中的逆质量把高维顶点重新称为低维。

记 $\mathfrak b=\Delta(B-L)$，与1035表中的整数 $b=3(B-L)$ 区分。各基本字段资料为

|字段|$Y_0$|$B-L$|正则维数|
|---|---:|---:|---:|
|$Q$|$1/6$|$1/3$|$3/2$|
|$u^c$|$-2/3$|$-1/3$|$3/2$|
|$d^c$|$1/3$|$-1/3$|$3/2$|
|$L$|$-1/2$|$-1$|$3/2$|
|$e^c$|$1$|$1$|$3/2$|
|$N$|$0$|$1$|$3/2$|
|$H$|$1/2$|$0$|$1$|

共轭取相反荷；$D$ 与规范场强的两种荷为零，维数为1和2。一个多项式可按 $\mathfrak b$ 分解，规范及 $r$ 作用保此分解。因此只需逐个非零荷单项证明；不同荷分量不能靠相互相消获得不变性。

**定理。** 在上述类内，允许非零 $B-L$ 破坏顶点当且仅当 $k=0$。其非零荷部分必为 $d=6$ 的四 $N$ Lorentz 标量及其共轭的线性组合，$\mathfrak b=\pm4$。至少两个味时该空间非零，原三代菜单满足。定理只给许可；不迫使许可系数实际非零。

## 2. 完备必要性，不靠有限扫描

### 2.1 不采用 $r$ 已能证明的部分

颜色中心不变要求净夸克数为3的倍数，故 $\mathfrak b\in\mathbb Z$。Lorentz 标量含偶数个费米字段 $n_F$：旋量中心对每个 Weyl 或其共轭给负号，对导数、场强和标量给正号。每个费米的 $3(B-L)$ 都是奇整数，所以

$$
3\mathfrak b\equiv n_F=0\pmod2,
\qquad \mathfrak b\in2\mathbb Z.
\tag{2}
$$

逐字段有 $|Y_0|\le(2/3)d_{\rm field}$。三角不等式给所有允许的导数、场强及共轭单项

$$
|Y_0(\mathcal O)|\le\frac23d\le4.
\tag{3}
$$

若 $k\ne0$、$\mathfrak b\ne0$，规范中性要求

$$
|Y_0(\mathcal O)|=3|k\mathfrak b|\ge6,
\tag{4}
$$

与(3)矛盾。因此**低维非零 $B-L$ 破坏要求 $k=0$，本身不需要 $r$**。不能在删弱比较中声称“删去 $r$ 就恢复非零 $k$ 的低维通道”。

### 2.2 精确 $r$ 将首个通道限定为四 $N$

1035的同一荷字典有

$$
X_k=5(B-L)-4Y_k,
\qquad r=\exp(i\pi X_k/2).
\tag{5}
$$

它对 Weyl 为 $i$，对 $H$ 为 $-1$；对共轭、协变导数及场强也与上述指定作用一致。规范中性单项有 $X_k=5\mathfrak b$。结合 $\mathfrak b\in\mathbb Z$，$r$ 中性等价于

$$
\mathfrak b\in4\mathbb Z.
\tag{6}
$$

另一方面 $d\le6$ 迫使 $n_F\le4$，且每个费米 $|B-L|\le1$。非零情形只能为 $\mathfrak b=\pm4$，并同时饱和：四个费米、各取同号单位 $B-L$、维数恰6。没有剩余维数放置 $H,D$ 或场强。这也排除了“加导数／共轭强度可能漏过”的余项。

考虑 $\mathfrak b=4$。仅有 $N,e^c,L^\dagger$ 能提供单位正荷，$Y_0$ 分别为 $0,1,1/2$，所以 $0\le Y_0\le4$。$Y_k=Y_0+12k=0$ 对整数 $k$ 只能为 $k=0,Y_0=0$，即四个 $N$。负荷情形取共轭。四个左旋量的 Lorentz 单态由两两反对称收缩张成，故其空间可写为

$$
(N_iN_j)(N_kN_l),\qquad
N_iN_j:=\epsilon_{\alpha\beta}N_i^\alpha N_j^\beta,
\quad \epsilon_{12}=1,
\tag{7}
$$

模去 Schouten／Grassmann 关系后的线性张成空间。本结论不把所有配对系数当作线性独立。空间非空由下一节及其在壳矩阵元证明，而不是只数荷。

## 3. 两味非零、质量保护与真空前件

取独立 Grassmann 生成元 $a=N_1^1,b=N_1^2,c=N_2^1,d=N_2^2$。则

$$
N_1N_1=2ab,\quad N_2N_2=2cd,\quad
\mathcal O_4:=\frac14(N_1N_1)(N_2N_2)=abcd\ne0.
\tag{8}
$$

相反 $(N_iN_i)^2=0$。两味已经足够；第三味不被删除，也不凭此选定代际数。$k=0$ 时 $N$ 是规范单态，故(8)为真正规范及 Lorentz 标量，不依靠额外参考字段。其 $r$ 相位为 $i^4=1$。

原 Dirac Yukawa $LHN$ 中性；$NN$ 和 $(LH)^2$ 的 $r$ 相位为 $-1$。在非零 Higgs 真空上应使用规范补偿后的

$$
\widetilde r=r\,z_w,\qquad z_w=-I_2\in SU(2),
\quad \widetilde r H=H,\quad
\widetilde r N=iN,\quad\widetilde r L=-iL,
\quad\widetilde r^2=(-1)^F.
\tag{9}
$$

故 Higgs 本身不破坏这个保护。这里仍须选择未另行自发破坏 $\widetilde r$ 的相，并用保对称的有效重整化条件；这不是证明所有参数／背景都具有同一有隙真空。在此分支中精确离散对称排除 Majorana 二次项，而不排除(8)。满秩 Dirac Yukawa 和不同正质量仍可自由输入。

成熟的“Dirac 中微子伴随四单位轻子数破坏”已见 [Heeck–Rodejohann §II、式(6)及图1](https://arxiv.org/html/1306.0580v2#S2)。其模型借新增标量产生有效四费米项；本轮只借此事实作来源，未采用其额外字段或声称重新发现该机制。[1034](../../archive_1009_1043/research_note_1034.md)及[1062共同分支账](../../archive_1046_1062/_shared/common_branch_scope_after1062.md)的离散保护范围继续成立。

## 4. 同一作用的来源与连续 Ward 式

在 $k=0$ 的原局部规范—Higgs—Weyl作用上选定

$$
S=S_0+\int d^4x\sqrt{-g}\,\mathcal L_4,
\qquad \mathcal L_4=\frac{C\mathcal O_4+C^*\mathcal O_4^\dagger}{\Lambda^2}.
\tag{10}
$$

$C$ 是固定中性无量纲参数，$\Lambda>0$ 是声明的有效尺度。$C=0$ 和 $C\ne0$ 都允许；本轮不从认知公理选择 $C$。该项改变父作用，不能无误差地继承1293—1297精确二次模型的全部 CP 准备／传播结论。

定义局部 $B-L$ 相位变化为 $\delta\psi=i q_\psi\alpha(x)\psi$，并选电流约定使未破坏部分 $\delta\mathcal L_0=-j^\mu\partial_\mu\alpha$。于是

$$
\delta\mathcal L=-j^\mu\partial_\mu\alpha
 +\frac{4i\alpha}{\Lambda^2}(C\mathcal O_4-C^*\mathcal O_4^\dagger),
\qquad
\nabla_\mu j^\mu=-\frac{4i}{\Lambda^2}
 (C\mathcal O_4-C^*\mathcal O_4^\dagger)
\tag{11}
$$

为在壳 Ward 式。积分分部及场方程给第二式；同样，$H_4=-\int\mathcal L_4$ 与 $[Q_{B-L},\mathcal O_4]=-4\mathcal O_4$ 给相同的 $\dot Q$ 符号。保留标准含 $N$ 的无反常连续 $B-L$ 电流归一，式(11)显示显式四单位破坏插入；不是把 Dirac 质量错误地当作连续 $B-L$ 守恒。

完整来源统一取 $S$ 的变分，例如签名 $(+---)$、$T_{\mu\nu}=2(-g)^{-1/2}\delta S/\delta g^{\mu\nu}$，以及规范、Higgs和参数来源。固定局域旋量分量时，(10)没有显式 vierbein／导数，故其 Hilbert 项为

$$
(T_4)_{\mu\nu}=-g_{\mu\nu}\mathcal L_4,
\qquad \frac{\partial\mathcal L}{\partial C}=\frac{\mathcal O_4}{\Lambda^2},
\quad \frac{\partial\mathcal L}{\partial C^*}=\frac{\mathcal O_4^\dagger}{\Lambda^2}.
\tag{12}
$$

同一变分还必须保 $S_0$ 的全部动能、规范、Higgs及 Yukawa 来源；在物理规范中重写变量不等于删去它们。(11)只是一个 Ward 身份，不是已经数值认证全 Hilbert 来源、曲域参考或实际装置。一般背景上的总守恒用完整场方程；本轮树级在壳例不提供装置来源预算。

## 5. 一个非零质量、有限动量的同作用四点例

选择原三代 Dirac 质量 $(m_1,m_2,m_3)=(3,4,5)$，其余 Dirac Yukawa 也满秩，$C=1,\Lambda=100$。这些是公开示例输入，未重新选择群／质量规律。两种参与味不简并。取质量本征场 $\Psi_i=(N_i,\ell_i^\dagger)^T$，其中 $\ell_i$ 是配对的左中微子。

使用手性 Hamiltonian 和 Pauli 矩阵

$$
h_m(\mathbf p)=
\begin{pmatrix}-\boldsymbol\sigma\cdot\mathbf p&mI_2\\mI_2&\boldsymbol\sigma\cdot\mathbf p\end{pmatrix}.
\tag{13}
$$

协变单粒子归一为 $\langle p,s|q,s'\rangle=(2\pi)^3 2E_p\delta_{ss'}\delta^3(p-q)$。取单位二旋量 $\xi,\eta$，外腿为

$$
u_m(p,\xi)=
\binom{\sqrt{EI_2-\boldsymbol\sigma\cdot\mathbf p}\,\xi}
 {\sqrt{EI_2+\boldsymbol\sigma\cdot\mathbf p}\,\xi},
\quad
v_m(p,\eta)=
\binom{\sqrt{EI_2-\boldsymbol\sigma\cdot\mathbf p}\,\eta}
 {-\sqrt{EI_2+\boldsymbol\sigma\cdot\mathbf p}\,\eta}.
\tag{14}
$$

直接有 $h(p)u=Eu$、$h(-p)v=-Ev$ 以及 $u^\dagger u=v^\dagger v=2E$。负频外腿的空间 Fourier 指数与正频相反，所以第二式用 $h(-p)$，不能混用 $h(p)$。

考虑由 $\mathcal O_4$ 给出的两味交叉通道 $\bar\nu_1\bar\nu_1\to\nu_2\nu_2$。等价地，这是 $\Psi_1$ 的两粒子入射与 $\Psi_2$ 的两反粒子出射，连续 $B-L$ 总荷改变 $-4$。在质心系固定每粒能量 $E=5$：

$$
\mathbf p_1=4\hat z,\quad\mathbf p_2=-4\hat z,
\qquad\mathbf p_3=3\hat z,\quad\mathbf p_4=-3\hat z.
\tag{15}
$$

各对第一腿取 $\downarrow$，第二腿取 $\uparrow$。$N$ 对应(14)上半外腿，因此

$$
x_1=3\downarrow,\quad x_2=3\uparrow,
\qquad y_3=\sqrt8\downarrow,\quad y_4=\sqrt8\uparrow;
\qquad x_1^T\epsilon x_2=-9,\quad y_3^T\epsilon y_4=-8.
\tag{16}
$$

式(8)中 $1/4$ 被两对相同字段的 $2!\,2!$ 消去。采用 $S=\mathcal T\exp(i\int\mathcal L_4)$，有序 amputated 顶点为 $iC\epsilon_{\alpha\beta}\epsilon_{\gamma\delta}/\Lambda^2$。定义 $i\mathcal M$ 为该 S 矩阵去掉整体动量守恒 delta 后的系数，则

$$
\mathcal M=\frac{C}{\Lambda^2}
 (x_1^T\epsilon x_2)(y_3^T\epsilon y_4)
 =\frac{72C}{\Lambda^2}=\frac9{1250},
\qquad |\mathcal M|^2=\frac{81}{1562500}>0.
\tag{17}
$$

整体相位随 S 矩阵约定改变不影响非零性。费米外腿交换使振幅反对称，故这不是把两相同外态当成可区别粒子；选定的反向动量／自旋态本就正交。

不只一个孤立点：保持同一配对及方向、$5\le E\le21/4$ 时

$$
|\mathcal M(E)|=\frac{|C|}{\Lambda^2}
 \big(E+\sqrt{E^2-9}\big)\big(E+\sqrt{E^2-16}\big)
 \ge\frac{72|C|}{\Lambda^2}.
\tag{18}
$$

两个因子在该窗递增；两体相空间的入／出速度分别至少 $4/5$、$3/5$，都非零。自旋随末对方向共同旋转时反对称收缩保不变，也可由连续性取得非零角邻域。适当支持于这些动量邻域的正常波包具有非零领先阶矩阵元。

这是(10)在 $C/\Lambda^2$ 首阶的树级过程，不给循环／高阶余项上界。总能量至多 $21/2<\Lambda$ 只检查所选有效输入的分离，不自动证明误差小；平面波与有限动量窗也不等于已实现紧时空准备、有限探测概率或完整记录装置。

## 6. 删除比较与不能推广的量词

1. **删除 $r$：** 第2.1节仍排除全部非零 $k$ 的 $d\le6$ 非零 $B-L$ 顶点；但 $k=0$ 的 $NN$（维3）及 Weinberg 项（维5）不再被离散关系禁掉。于是 Dirac 保护及“第一类只能是四 $N$”失效，低维 $k$ 必要性不失效。
2. **加入额外带荷来源：** 例如新复字段 $S$ 取 $Y(S)=-6k,r(S)=-1$，可使 $SNN$ 中性；需把其动能、来源、相与可能凝聚纳入另一个父作用。若同时赋 $B-L(S)=-2$，它保持的是扩大后总 $B-L$，不能误称扩大后总荷被该项破坏。若赋其它 $B-L$，则必须重做全荷账。该例表明原字段菜单必要，不免费移植其量子／全局资格。
3. **放开维数：** 在原字段菜单内定义

$$
J_\mu=H^T\epsilon D_\mu H,
\quad Y_k(J)=1,\quad B-L(J)=0,\quad r(J)=1,\quad d(J)=3.
\tag{19}
$$

$H$ 与 $D_\mu H$ 是不同 jet，故反对称收缩不恒零。局部取 $H=(1,0)^T,D_0H=(0,1)^T,D_iH=0$ 即 $J_0=1,J_i=0$。对任意 $k\ne0$，有非零高维局域许可

$$
\mathcal O_{4,k}=
\begin{cases}
\mathcal O_4(J_\mu^\dagger J^{\dagger\mu})^{6k},&k>0,\\
\mathcal O_4(J_\mu J^\mu)^{6|k|},&k<0,
\end{cases}
\quad d=6+36|k|,
\quad Y_k=0,\quad \mathfrak b=4.
\tag{20}
$$

这只是高维许可非空，不认证其有效率、系数、实际制备或在任何背景上可用。因此不能把低维条件结论推广为永久排除全部非零 $k$。

4. **令 $C=0$：** 保留全部低维、反常及 Dirac 条件而关闭所选通道。认知机制并未证明该系数必须非零；如果未来明确要求这类非零低维集体交换，它才给 $k=0$ 的条件性结构筛选。

## 7. 净增量、复算及停止线

成熟离散选择规则、一般四单位轻子数破坏、旧全局反常与一般树级 Wick 规则均直接复用。本轮新增联合证书限于：**原无限整数超荷族的完整 $d\le6$ 菜单 → 允许的四 $N$ 结构 → 非零 Grassmann 顶点 → 同一作用、非简并正质量的非零四点过程**；并精确拆开低维 $k$ 筛选与 Dirac 离散保护的前件。

[作者代码](check.py)从有理字段荷枚举必要条件（不充当完备性证明）、独立外代数乘法及显式四分量外腿核验；[结果](results.json)保存一组，默认和 `--check` 只读复算，`--write` 仅首次排他创建。源哈希绑定原件与本证明，数值只校验归一，式(2)–(6)、(18)、(20)的全称量词由解析承担。

交付是结构 S 和声明阶矩阵元；反常 Q 仅继承1060原背景范围。没有借此认证一般非微扰 EFT、实际任意伙伴仪器 O、完整热几何 M4 或长程目标完成；也不把1293的二次仪器结果与此新四费米作用并列为同一个已认证装置。后续若只换味数、优化小概率或再次做荷相加，不再新开轮。
