# 1181证明：固定EPRL边界的配对、形状记录与完整输出

2026-10-13。作者草稿，三组。所有结果均在本节明确的有限静态接口内成立。

## 1．对象和成熟来源

令 $\hbar=1$。取Euclidean Spin(4)的标准张量积内积，

$$
\gamma=\frac13,\quad k_i=\frac32,\quad
j_i^+=\frac{1+\gamma}{2}k_i=1,\quad
j_i^-=\frac{1-\gamma}{2}k_i=\frac12 .
\tag{1}
$$

$\mathcal H_0=\operatorname{Inv}_{SU(2)}(V_{3/2}^{\otimes4})$ 维4；$\mathcal H_1=\operatorname{Inv}(V_1^{\otimes4})\otimes\operatorname{Inv}(V_{1/2}^{\otimes4})$ 维6。每腿最高自旋CG嵌入 $Y:V_{3/2}\to V_1\otimes V_{1/2}$ 等距，定义

$$
F=P_{\mathrm{Spin}(4)}Y^{\otimes4}|_{\mathcal H_0},\qquad
G=F^\dagger F,\qquad \mathcal S=\operatorname{im}F.
\tag{2}
$$

采用Condon–Shortley实CG相位及耦合 $(12)\to a,(34)\to a\to0$ 的正交归一基，原码 $a=0,1,2,3$；目标行次序 $(a_+,a_-)=(0,0),(0,1),(1,0),(1,1),(2,0),(2,1)$。不变张量逐项为

$$
I^{(j,a)}_{m_1m_2m_3m_4}
=\frac{(-1)^{a-M}}{\sqrt{2a+1}}\,
C^{aM}_{jm_1,jm_2}C^{a,-M}_{jm_3,jm_4},
\quad M=m_1+m_2,\quad \sum_i m_i=0.
\tag{3}
$$

代码直接对四腿磁指标收缩(3)与四个单腿CG；先独立检验三套基正交归一和 $Y^\dagger Y=I$。

EPRL非等距及以逆Gram表达映像正交投影，归属[KKL，0912.0540v2，§II.0.4–6、式20–21](https://arxiv.org/html/0912.0540v2)。注入性的一般适用条件见[1109.5023，定理2](https://arxiv.org/html/1109.5023)；本例 $\sum j_i^+=4$ 为整数，两个目标不变空间均非零，另有本轮直接正谱证书。因此不把早期仅SO(3)版本无条件套到半整数输入。[Barbieri量子四面体](https://arxiv.org/html/gr-qc/9707010)提供闭合不变空间与非对易形状记录的成熟几何解释。

## 2．直接CG结果与完整形状记录

精确收缩给

$$
F=
\begin{pmatrix}
2/3&0&0&0\\
0&10/27&0&0\\
0&5/9&0&0\\
1/9&0&2\sqrt5/9&0\\
0&0&1/3&0\\
0&2\sqrt5/135&0&\sqrt{105}/15
\end{pmatrix},
\quad
G=
\begin{pmatrix}
37/81&0&2\sqrt5/81&0\\
0&181/405&0&2\sqrt{21}/405\\
2\sqrt5/81&0&29/81&0\\
0&2\sqrt{21}/405&0&7/15
\end{pmatrix}.
\tag{4}
$$

$$
\operatorname{spec}G=\left\{\frac13,\frac{35}{81},\frac{13}{27},\frac{13}{27}\right\},
\qquad \det G=\frac{5915}{177147}>0.
\tag{5}
$$

所以 $F$ 注入但非等距，也不是标量倍等距。原SU(2)码中定义

$$
q_{12}=\operatorname{diag}(0,2,6,12),\qquad
q_{13}=
\begin{pmatrix}
15/2&-5\sqrt3/2&0&0\\
-5\sqrt3/2&13/2&-4\sqrt{15}/5&0\\
0&-4\sqrt{15}/5&9/2&-3\sqrt{35}/10\\
0&0&-3\sqrt{35}/10&3/2
\end{pmatrix}.
\tag{6}
$$

$q_{13}$ 由实际交换腿2、3在不变基中的矩阵 $R_{23}$ 给出：$q_{13}=R_{23}q_{12}R_{23}^\dagger$。这是对同一磁基张量的独立重排，不靠拟合Gram。记录合同包含各自四个谱结果及读后态；两份读口在不同任务选择，未声称非对易形状同时锐测。

两者与 $G$ 均不对易。若矩阵 $M$ 与 $q_{12}$ 对易，简单谱使它在该基中对角；再与 $q_{13}$ 对易，由非零链 $0-1-2-3$ 迫四个对角元相同。因此

$$
\{q_{12},q_{13}\}'=\mathbb C I_4 .
\tag{7}
$$

这也适用于由它们的谱投影生成的全部代数。

## 3．同边界无新事件身份唯一选择配对

固定上述原内积、映像与全部输入支。此次候选合同把一份插入明确识别为**同一边界、不引入新事件的无损交接**。它要求振幅算符在 $\mathcal S$ 上恒等，而不是任意物理更新必须恒等。对于待选配对矩阵 $W$，

$$
P_W=FWF^\dagger,\qquad
P_WF=F
\ \Longleftrightarrow\
F(WG-I)=0
\ \Longleftrightarrow\
W=G^{-1}.
\tag{8}
$$

最后一步使用注入性，未丢弃任何Gram支，也不需要先假定 $W$ 正定。解自动正定，给正交投影

$$
P=FG^{-1}F^\dagger=P^\dagger=P^2,\quad PF=F,\quad
\operatorname{rank}P=4.
\tag{9}
$$

整体共同相位按无新事件的约定固定为1；若只比较无相位参考的通道，则同一结论至多留下共同相位，它不恢复标量权。式(8)是KKL逆Gram配对在本次完整边界合同下的唯一性，不是新的EPRL非等距定理，更不选择顶点动力学。

对实标量 $W=cI$，$\mathcal S$ 中 $FF^\dagger$ 的谱就是(5)。故

$$
\min_{c\in\mathbb R}\|(cFF^\dagger-I)|_{\mathcal S}\|
=\frac{\lambda_{\max}-\lambda_{\min}}
{\lambda_{\max}+\lambda_{\min}}
=\frac2{11},\qquad c_*=\frac{27}{11}.
\tag{10}
$$

端点处误差等幅反号，中间特征值处于区间内，证明全谱的极小极大结论。该最优算符有最大特征值 $13/11>1$；(10)是范数证书，不能称某实际成功仪器的概率误差。

若欲让原 $F$ 通过更换原码内积而变等距，需用 $\langle x,y\rangle_G=x^\dagger Gy$。但保留同一 $q_{12},q_{13}$ 的Hermitian身份要求

$$
q_\alpha^\dagger M=Mq_\alpha\quad(\alpha=12,13).
\tag{11}
$$

原矩阵Hermitian及(7)迫 $M=cI$；原归一化再固定 $c=1$。非标量 $G$ 不满足这个必要条件。改内积是一条新物理／读口合同，不能免费解释成原两份完整记录不变。

## 4．原未知码／参考与极分解运输

逆Gram修正边界配对，不会同时把原 $F$ 变为等距。原码若需要确定性无损编码，取正平方根规范的polar映射

$$
V=FG^{-1/2},\qquad V^\dagger V=I_4,\quad VV^\dagger=P.
\tag{12}
$$

对任意未知 $\rho_{0R}$，$\rho_{0R}\mapsto(V\otimes I_R)\rho_{0R}(V^\dagger\otimes I_R)$ 保全部相干和与任意 $R$ 的关联；在支撑 $\mathcal S$ 上由 $V^\dagger$ 恢复。此处是精确有限维恒等式，不是只在几个已知态上验真。

原形状谱投影 $E_{\alpha,a}$ 对应

$$
\widetilde E_{\alpha,a}=V E_{\alpha,a}V^\dagger,\quad
\sum_a\widetilde E_{\alpha,a}=P,\quad
\widetilde E_{\alpha,a}V=VE_{\alpha,a}.
\tag{13}
$$

它们在 $\mathcal S$ 上是完整PVM；若需全 $\mathcal H_1$ 仪器，另加 $I-P$ 的域外结果，编码输入上其概率严格为0。带所有结果的Lüders读后态也由(13)完全运输，任意 $R$ 同样成立。更一般原有限仪器可在映像上用 $V A_rV^\dagger$ 运输并在补空间补全。

这保的是**已声明的原几何记录及其运输**，不宣称它们自动等于未经运输的Spin(4)原生面积／角算符。$V$ 改变原 $F$ 的非均匀输入权；$FqF^{-1}|_{\mathcal S}$ 一般又不Hermitian。故不能把polar归一化说成原嵌入、原生读口和全部振幅全都不变。

## 5．保留失败与成功条件化的固定实例

为把范数与概率严格分开，另声明一份数学滤波仪器，作用在 $\mathcal H_1$：

$$
K=\frac{27}{13}FF^\dagger,\quad 0\le K\le I,\qquad
L=(I-K^\dagger K)^{1/2},\qquad
\mathcal J\phi=|s\rangle K\phi+|f\rangle L\phi.
\tag{14}
$$

$\mathcal J^\dagger\mathcal J=I_6$，所以保完整相干旗标及后态时，对所有输入／$R$ 是等距。$L$ 保含 $\mathcal S^\perp$ 的两维，声明输入码仍为 $\mathcal S$。此处滤波对象是裸配对 $FF^\dagger$，不是原 $F$ 单次滤波；后者在polar识别下给平方根特征值，不能套用下述概率。旗标实际读出后的instrument有两个结果，满足 $K^\dagger K+L^\dagger L=I$；它是合法完整仪器，但已发生一次有信息的读出，不能把测量后的经典旗标自动认作相干可逆交接。尤其不得丢掉失败后再称确定性。

取 $\mathcal S$ 中正交归一的最小／最大谱向量 $u_-,u_+$，令

$$
\psi=\frac{u_-+u_+}{\sqrt2},\qquad
K u_-=\frac9{13}u_-,\quad Ku_+=u_+.
\tag{15}
$$

随后以原已标定 $\psi$ 投影回读，则精确有

$$
p(s)=\frac{125}{169},\quad
p(s,\psi)=\frac{121}{169},\quad
p(\psi\mid s)=\frac{121}{125},\quad
1-p(\psi\mid s)=\frac4{125}.
\tag{16}
$$

理想无事件恒等输出的 $\psi$ 回读为1。这里(16)只诊断成功条件化造成的明确差异，不能把条件概率与(10)的范数、或与不带成功设置的另一个完整联合实验混为同一误差。归一成功更新一般随输入成功率而非线性；其未归一成功CP支与完整两支instrument仍线性。整套仪器由已知固定矩阵给出，不复制未知输入，也未假定实际硬件已实现。

## 6．计算与来源冻结范围

三组检查依次覆盖直接CG／谱、完整配对／记录、polar／保旗标输出。精确算术以有理系数平方根和保存所有矩阵；谱投影由 $G$ 的三个不同有理特征值构造，精确核平方根、等距和完整仪器恒等式。NumPy谱及残差仅作展示，不承担证明。

首次结果显式 `--write` 创建，已存在则拒绝；默认从头重算并只读比对。作者代码不导入独审或预览结果；预览源码的通用CG实现复制后发展，预览本身不改。来源绑定与正式值见 `results.json`，最终六核SHA由冻结交接另报。

本轮只选择固定边界图、标签、简单性嵌入与原内积下的配对身份。没有证成实际动态承担者、钟、物质、有限能量实施、来源反作用、边界 $T_{21}$、Lorentz接续、$\gamma$ 选择、顶点或Einstein方程。EPRL约束映射本来不是CPTP时间演化；本轮未将二者混同，也不从这一接口结果宣布EPRL或整个GR纲领失败。
