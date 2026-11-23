# 已命名的破缺标量能否保护 Nelson–Barr 质量结构？

## 1. 耦合的对称性与正向、反向作用

取三个普通下型味和一个矢量下型弱单态 $D$。在 CP 基中普通下型 $m$ 可逆实、$\det m>0$，上型实正，裸强角取零支。单态重质量 $M>0$ 是对称允许的常数；一个或多个标量 $S$ 以实 CP 系数产生复行 $B(S)$。上右列记 $c$。以下结论属于共同背景和树级匹配，沿用第1364轮已得的精确 Schur 补

$$
\bar\theta=\arg(1-z),\qquad z=B m^{-1}c/M .
$$

取 $-\bar\psi_L\mathcal M\psi_R+\mathrm{h.c.}$ 的手征约定，已计入完整四个下型夸克；积分掉重态还要匹配其相位。第1364轮已研究一个 $Z_2$ 门、一列任意小扰动和保持左 Gram 的单独右相位；第1365轮在同一个径向 Higgs 下给出物理 Higgs 耦合；第1388轮说明真实来源需要在另一个过程获得标定。这里进一步问普通的内部对称能否保护零列，或至少保护强相位。[1364](../1364/proof.md) · [1365](../1365/proof.md) · [1388](../1388/NEXT.md) · [Nelson](https://doi.org/10.1016/0370-2693%2884%2992025-2) · [Barr](https://doi.org/10.1103/PhysRevLett.53.329)

令紧群 $G$ 在正则动能的有限维场空间幺正作用。已有的普通 Yukawa 给出 $Y(H):d_R\to Q_L$，标量混合给 $B(S):d_R\to D_L$，重质量使 $J:D_R\to D_L$ 成为常数的等距交织，背景值 $MJ$。三个网络均按左右端口的表示变换，例如 $Y(gH)=U_QY(H)U_d^\dagger$。紧致性保证可取不变内积，$B(gS)^\dagger=U_d B(S)^\dagger U_{D_L}^\dagger$。于是

$$
\mathcal C(H,S)=Y(H)\,B(S)^\dagger J:\ D_R\longrightarrow Q_L
$$

自动协变。其标量含 $H$ 以及与底行相应的共轭标量；若 $B$ 本来同时含 $S,S^\dagger$，反向式相应含两种。两费米场和两个标量给规范维数五，系数 $\kappa/\Lambda$，模型数例取 $\langle H\rangle=1$ 个质量单位。这一证明包括非阿贝尔内部作用及非平凡 Higgs 表示。一般在具体背景的非零条件是 $\operatorname{im}B^\dagger J\not\subset\ker Y(H)$，这里 $m$ 满列秩、$B\ne0$ 保证非零。每个指定味分量仍可因零元或对称扇区消失。若 $M$ 也必须由新带荷场产生，就没有假设的常数交织；手性规范拓展可以把危险项推至更高维。幺正对称只给许可，不能确定 Wilson 系数及其是否为零。[手性 Nelson–Barr 保护及维五问题](https://doi.org/10.1103/PhysRevD.107.115012)

## 2. 对齐与全系数判据

对最直接的 $c=(\kappa/\Lambda)mB^\dagger$，实 $\kappa$ 给 $z=\kappa\|B\|^2/(\Lambda M)$，完全为实。$1-z>0$ 时强角精确为零，等于零时行列式为零，不能沿用有质量相位公式，小于零时在 $\pi$ 支。因而逆向门存在本身不破坏树级强 CP。

更一般取 $c=(\kappa/\Lambda)m X B^\dagger$，$X$ 为已允许的无量纲右味算符。对固定 $B$，
$\operatorname{Im}(BXB^\dagger)=B(X-X^\dagger)B^\dagger/(2i)$。对所有允许的小实 Wilson 系数独立保护，相当于这一实线性泛函在全部许可插入的实张成空间上恒零；大幅扰动另须保持 $1-\operatorname{Re}z>0$。对全部复行 $B$ 都保护的充要条件是每个 $X$ Hermitian，由 Hermitian 二次型的极化恒等式得到。如果物理上只允许 $B$ 在某真子空间，须在该空间压缩后检验，不能擅自推广。

**一个已有的左插入已经可能破保护。** 令 $L$ 为现成的 Hermitian 左味张量，$c=(\kappa/\Lambda)Lm B^\dagger$。则 $X=m^{-1}Lm$，并有

$$
X=X^\dagger\quad\Longleftrightarrow\quad [L,mm^\dagger]=0 .
$$

即使 $L$ 本身 Hermitian，把它通过非等距 $m$ 转到右边也未必 Hermitian。普通上型 Yukawa 可给 $L=\sum_\alpha Y_{u,\alpha}Y_{u,\alpha}^\dagger$，希格斯分量按不变内积收缩。这是维数零的味张量，不额外插入标量场，整个反门仍是维五。对于单个通常的 Higgs 双态，去掉共同的弱恒等因子即为所用 $Y_uY_u^\dagger$。一般内部表示须用完整分量收缩，不能把任一真空分量平方当作原群的不变量。

**固定破缺背景也可判所有既有味插入。** 在普通 Yukawa 全实的 CP 基中，令 $\mathcal A_R$ 为已有实右味张量的实结合代数。例如由 $X_d=\sum Y_{d,\alpha}^\dagger Y_{d,\alpha}$ 与 $Z=\sum Y_{d,\alpha}^\dagger L Y_{d,\alpha}$ 生成。每个词 $W$ 均给独立许可的复合映射 $Y(H)W B(S)^\dagger J$，其 CP 合法性容许实 Wilson 系数。若 $\mathcal A_R=M_3(\mathbb R)$，矩阵单位 $E_{ij}$ 都在实张成中。所有这些门在正行列式邻域内保护强相位，当且仅当 $\operatorname{Im}(B_i B_j^*)=0$ 对所有 $i,j$ 成立，即 $B=e^{i\alpha}b$ 且 $b$ 为实行（零分量不妨碍）。可用重费米子左右的相同相位消去 $e^{i\alpha}$ 而保留正 $M$，于是当 $B$ 是此质量部门唯一 CP 来源时，四列父质量与普通弱流都可取实。**在满实味代数并实际具有弱 CP 时，总有一种利用已有味张量的反向维五作用会对强角贡献虚部。** 存在可选作用没有给出其最低非零系数，也没有推出辐射修正的大小。

一个可查的满代数充分条件是 $X_d$ 实对称非简并，在它的谱基中另一个实对称 $Z$ 的非零边图连通。谱投影是 $X_d$ 的实多项式，以投影夹住 $Z$ 选出边矩阵单位，沿路径连乘得到全部单位。对称分块或味简并的模型可另有保护；一般应回到对 $\mathcal A_R$ 的线性泛函判据。

## 3. 同一物质、弱 CP 与可留出的强相位

用第1364轮的 $m$ 和 $B=(1,i,1+i)$、$M=4$、$\langle S\rangle=1+i$，另登记 $\Lambda=10$（同一质量单位）、$L=\operatorname{diag}(1,4,9)/9$。此 $L$ 来自一个实非简并上型作用；保持旧例的弱质量基。$\kappa=1$ 时，对齐门给 $z=0.1$、$\bar\theta=0$，对角化后的轻弱四元积为 $0.005341$。加 $L$ 的门给

$$
z=\kappa(161/2160+i/180),\qquad
\bar\theta=-0.00600293\quad(\kappa=1),
$$

弱四元积为 $0.005550$。同一 $m,B,L$ 即可判别两种方向，不需发明新的复杂味张量。数例取 $X_d=m^Tm,\ Z=m^T(9L)m$。$X_d^iZX_d^j$（$i,j=0,1,2$）的九个展平行有非零**整数**行列式 $300321699195100225536$，直接证得满实代数。符号和尺寸没有假定现实夸克质量等于数例。[精确代数与数值](check.py) · [结果](results.json)

另外还有一份**质量、全部左弱流及同一常数正 $M$ 完全相同**的树级对照。在零反门的有限母矩阵右端乘酉算符，只在第一个普通右手和重右手之间用下列块，其余单位：

$$
R=\begin{pmatrix}t&s e^{i\beta}\\-s e^{i(\phi-\beta)}&t e^{i\phi}\end{pmatrix},
\quad t=\cos\alpha,\ s=\sin\alpha,\quad
\cos\phi=\frac{1+t^2-|B_1/M|^2s^2}{2t},\quad
e^{i\beta}=\frac{M(1-te^{i\phi})}{B_1s}.
$$

例如 $\alpha=0.1$ 给 $\phi=0.02451517$。新上左块是 $m\operatorname{diag}(t,1,1)$，实且正行列式；重块仍为 $M$。另外两个方向可用同一 $Z_2$ 下 $S,S^*$ 的不同实系数产生：任意复向量 $w$ 取 $g=(\operatorname{Re}w+\operatorname{Im}w)/2$、$h=(\operatorname{Re}w-\operatorname{Im}w)/2$，则 $gS+hS^*=w$ 在该背景成立。上右系数另带 $H/\Lambda$。这次同读竞争容许**独立**反向张量，不声称它就是唯一的 $LmB^\dagger$。两理论选同一裸角零支，$\det R=e^{i\phi}$ 使强角相差 $\phi$；因 $\mathcal M\mathcal M^\dagger$ 不变，质量和左流相同，第1365轮的径向配接条件下相应 Higgs 顶点也相同。此为比较两个 UV CP 对称的不同耦合实现；若只是被动手征换基，必须连同反常移动裸角，物理强角不变。具名标量的涨落顶点或其它高能读数能区别两者。

## 4. 中子端口的条件定量

PSI 的 $90\%$ 实验限为 $|d_n|<1.8\times10^{-26}\ e\,\mathrm{cm}$。采用一份历史格点计算 $d_n=-0.00148(14)(31)\,\bar\theta\ e\,\mathrm{fm}$；若把两项误差线性展开一倍作敏感性框，斜率模的条件下界为 $1.03\times10^{-16}\ e\,\mathrm{cm}$。假定近零支、小角线性有效，阈值、环修正、独立拓扑源和其它 EDM 贡献已控制，本例把剩余 EDM 预算取零，得 $|\bar\theta|\le1.75\times10^{-10}$。若独立预算为 $b$，分子须换为 $d_{\rm lim}+b$。实验限和格点误差的组合没有新的联合置信水平。格点的连续、手征外推等系统误差属于该文的响应输入。[PSI 实验](https://arxiv.org/abs/2001.11966) · [格点 QCD 斜率](https://arxiv.org/abs/2301.04331)

对单独的 $L$ 插入，写上式 $z=\kappa(r+ij)$，$j>0$，令 $t_\theta=\tan\theta_{\max}$。在正支恰有 $-t_\theta/(j-r t_\theta)\le\kappa\le t_\theta/(j+r t_\theta)$，即约 $\pm3.15\times10^{-8}$。取 $\kappa=10^{-7}$ 给 $\bar\theta=-5.56\times10^{-10}$，用中央斜率时 $|d_n|=8.22\times10^{-26}\ e\,\mathrm{cm}$。这是对指定无量纲物质比、唯一系数及斜率范围的条件留出；若多项抵消，约束实际只落在合计相位。没有从弱 CP 或这一数例算出现实强角。
