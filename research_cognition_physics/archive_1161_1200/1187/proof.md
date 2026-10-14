# 1187 解析证明：复合材料的真实钟与裸传播主符号

本文件是第1187轮的证明附件。研究对象限于 Batz–Erlich–Mrini 行列式候选在 $D=4$ 的**正实材料、正势、正则经典分支**。结论不针对其重整化集体真空，也不把该分支的失败提升为认知研究纲领或全部涌现引力的失败。

## 1. 复用、输入与新命题

原文入口：[Batz–Erlich–Mrini，Background-Independent Composite Gravity，式(1)、(5)、(9)及§IV](https://arxiv.org/html/2011.06541)。本轮从其行列式作用在所声明分支的表达式出发；以下正则分析、真实材料规约与线性诊断由本轮推导，不借用文献的量子真空结论。

历史接口：

- [366轮](../../archive_342_369/research_note_366.md)已证明完整重切片约束仍可对应固定 Minkowski 背景。本轮不重证形变代数，也不以约束数目生成引力。
- [858轮](../../archive_854_872/research_note_858.md)已要求材料参考由约化后的实际动力学计算。本轮遵守这条规则，直接用原作用的实际材料字段作钟，计算同一约化作用的传播。
- [359轮](../../archive_342_369/research_note_359.md)已限定 Weinberg–Witten 前提。本轮不使用该定理，也不预设或排除复合引力子的量子态。
- 旧 Gram 符号审计只判断有限正实材料矩阵的符号。本轮的新增量是：保留正则约束并落实真实材料钟后，裸分支的**有限波数线性响应**仍不等于 Lorentz 波传播。
- [1186证明](../1186/proof.md)区分条件协变背景身份与独立正量子见证。本轮不将二者合并，也不重复接触项推导。

认知动机是同一材料能否同时承担实际参考与共同传播。额外建模输入是所选行列式作用、正实字段、正则局部材料坐标和所声明的诊断支。解析结论不需要无限高能、无限精细测量或量子引力的完整紫外定义。

**待证命题。** 在这些输入下，选四个实际材料字段为局部坐标后，任意有限数目剩余字段的裸图像作用具有椭圆主符号；平坦图像线性化给出明确的传播诊断。常正势分支在固定有限波数和有限材料时间上给出双曲余弦响应，而不是 Lorentz 波方程的余弦响应。仅改变时钟方向或作用的总体符号不能消除此差别。

## 2. 原作用与适用域

设实际材料为 $\phi^a$，$a=1,\ldots,N$，参数坐标为 $x^\mu=(\tau,x^i)$，$i=1,2,3$。取

$$
S[\phi]=\int d^4x\,F(\phi)\sqrt{\det A},\qquad
F=\frac1{V(\phi)}>0,\qquad
A_{\mu\nu}=\sum_a\partial_\mu\phi^a\partial_\nu\phi^a.
\tag{1}
$$

原文行列式带绝对值。本轮限制 $A$ 为正定实 Gram 矩阵的正则分支，故这里 $\det A>0$。式(1)中的材料指标以正的 Euclidean 内积收缩；这不是已经推得的物理时空度规。参数 $\tau$ 起初只是分解变量，不能直接当成实际时钟。$N\ge4$ 才能满足非退化四维材料条件；要有剩余图像自由度，需 $N\ge5$。

定义 $N$ 维列向量 $e_i=\partial_i\phi$、$v=\partial_\tau\phi$，以及

$$
E=(e_1,e_2,e_3),\qquad
h=E^TE>0,\qquad
P_\perp=I-Eh^{-1}E^T,\qquad
u=P_\perp v\ne0,\qquad
W=\frac{\sqrt{\det h}}{V}>0.
\tag{2}
$$

$P_\perp$ 是材料内积下投向三个空间切向量之正交补的投影。以下所有正则式均限于式(2)成立的开集；不把 $h$ 退化、$u=0$ 或 $V=0$ 的边界混入同一坐标图。

## 3. Legendre 秩、四个主约束与零正则 Hamiltonian

Gram 矩阵的 Schur 补给出

$$
\det A
=\det h\left(v^Tv-v^TEh^{-1}E^Tv\right)
=\det h\,v^TP_\perp v,\qquad
\mathcal L=W\|u\|.
\tag{3}
$$

固定字段及其空间导数求速度导数，有

$$
\pi=\frac{\partial\mathcal L}{\partial v}
=W\frac{u}{\|u\|},\qquad
C_i=\pi\cdot e_i=0,\qquad
C_\perp=\frac12(\pi^2-W^2)=0.
\tag{4}
$$

令 $n=u/\|u\|$，速度 Hessian 为

$$
\frac{\partial^2\mathcal L}{\partial v\,\partial v^T}
=\frac W{\|u\|}\left(P_\perp-nn^T\right),\qquad
\operatorname{rank}\frac{\partial^2\mathcal L}{\partial v\,\partial v^T}=N-4.
\tag{5}
$$

其核由 $e_1,e_2,e_3,u$ 张成。约束关于动量的梯度分别为 $e_i$ 与 $\pi$，在本正则分支上线性独立，故得到四个独立主约束。速度的一次齐次性同时给出

$$
\mathcal H_{\rm can}=\pi\cdot v-\mathcal L=0.
\tag{6}
$$

本轮仅使用这些直接计算出的主关系，不另行声称已完成完整量子约束代数或物理态空间构造。式(6)也不意味着没有关系时间或物理变化：约束系统选择实际内部钟后可以有非零约化 Hamiltonian，第5节明确给出它。

## 4. 原文辅助坐标不能直接充当实际材料钟

原文式(5)将辅助 $X^M$ 项加上再减去，式(9)也说明总作用与 $X$ 无关。这种引入可以组织背景坐标下的展开；它本身没有提供参与材料动力学的时钟。

严格地说，若把该原作用扩展到独立的 $X$ 变量，则有独立局部冗余

$$
X^M(x)\longmapsto X^M(x)+\epsilon^M(x),\qquad
\Pi_M^X=0.
\tag{7}
$$

固定 $X^M=x^M$ 可固定这份新增的坐标冗余，但不会使式(4)的原材料约束凭空消失。也可用“参数微分同胚加补偿 $X$ 平移”的组合保持 $X=x$，同时仍改变原材料的参数描述。因此，不能仅根据辅助 $X$ 的固定，就把它解释成实际物质钟或已经完成物理 Hilbert 空间构造。

这里不判定文献的背景坐标组织错误，也不据此否定其散射推导。需区分的只是：**辅助坐标方便性与本项目所需的实际材料参考，是不同的接口。**

## 5. 四个真实材料字段的局部规约

现在选实际材料中的四个字段作坐标。在其 Jacobian 非零的局部区域，经材料指标重排，写为

$$
T=\phi^N=\tau,\qquad
\phi^i=x^i\quad(i=1,2,3),\qquad
\text{其余字段为 }\psi^A,\quad A=1,\ldots,N-4.
\tag{8}
$$

这是实际字段的局部非退化假设，不声称全域存在这样的材料坐标。定义

$$
J_{Ai}=\partial_i\psi^A,\qquad
h=I_3+J^TJ,\qquad M=I_{N-4}+JJ^T.
\tag{9}
$$

令 $p_A$ 为 $\psi^A$ 的动量，$p_T$ 为真实钟的动量。空间主约束直接解出

$$
\pi_i=-\sum_Ap_A\partial_i\psi^A.
\tag{10}
$$

代入法向约束得到

$$
p_T^2=\mathcal R,\qquad
\mathcal R=\frac{\det h}{V(T,x,\psi)^2}-p^TMp.
\tag{11}
$$

在 $\mathcal R>0$ 的固定代数支上，$p_T=\sigma\sqrt{\mathcal R}$，$\sigma=\pm1$。由辛项 $p_T\,dT+p_A\,d\psi^A$ 的拉回可知

$$
H_{\rm red}=-p_T=-\sigma\sqrt{\mathcal R}.
\tag{12}
$$

**方向限制。** 对原正作用、式(8)的正向参数化 $\dot T=1$，空间切向量的 $T$ 分量为零，故 $(P_\perp v)_T=1$。式(4)于是要求 $p_T>0$；这一具体定向分支给出

$$
H_{\rm red}=-\sqrt{\mathcal R}.
\tag{13}
$$

不能未经说明就把相反根当成同一正向原作用的另一个已实现物理 Hamiltonian。改变钟方向或作用总体符号会相应改变动量及 Hamiltonian 的约定，第7节说明它们不改变传播方程的类型。

未规约前的正则括号给出

$$
\{T(x),C_\perp(y)\}=p_T(x)\delta(x-y).
\tag{14}
$$

所以在有限任务区域内若 $\mathcal R\ge\kappa^2>0$，选定根和该局部钟是正则的。到 $\mathcal R=0$，这份时钟规约失去可逆性，需换图或另作分析；这不是“所有内部钟都不存在”的结论。

式(11)也不是现成的正自伴量子平方根构造。若进一步量子化，算符域、谱限制、时间依赖投影与相应接续尚须给出。本轮不把这些尚未研究的问题设为整个目标必须先完成的紫外门槛。

## 6. 图像作用及其主符号

把式(8)直接代入同一原作用，得到

$$
A_{\mu\nu}=\delta_{\mu\nu}
+\sum_A\partial_\mu\psi^A\partial_\nu\psi^A,\qquad
\mathcal L_{\rm graph}
=F(T,x,\psi)\sqrt{\det(I_4+ZZ^T)},\quad
Z_{\mu A}=\partial_\mu\psi^A.
\tag{15}
$$

此处 $\delta_{\mu\nu}$ 来自四个正实材料坐标在原正内积中的 Gram 矩阵，不是外加 Lorentz 度规。真实钟已经代入，不能再以辅助参考替换它而宣称同一规约。

### 6.1 单剩余字段：非线性主符号严格椭圆

只留一个剩余字段，矩阵行列式引理给出精确表达式

$$
\mathcal L_{\rm graph}=F\,s,\qquad
s=\sqrt{1+\partial_\mu\psi\,\partial_\mu\psi}
=\sqrt{1+\psi_T^2+|\nabla\psi|^2}.
\tag{16}
$$

Euler–Lagrange 方程为

$$
\partial_\mu\left(F\frac{\partial_\mu\psi}{s}\right)
-F_\psi s=0.
\tag{17}
$$

$F$ 不含字段导数。所有最高二阶导数的系数为

$$
\mathsf A^{\mu\nu}
=\frac F s\left(\delta^{\mu\nu}
-\frac{\partial^\mu\psi\,\partial^\nu\psi}{s^2}\right).
\tag{18}
$$

这里升指标只是使用式(15)中的 $\delta$。对任一非零实协向量 $\xi$，由 Cauchy–Schwarz 不等式，

$$
\xi_\mu\mathsf A^{\mu\nu}\xi_\nu
=\frac F s\left(|\xi|^2-\frac{(\partial\psi\cdot\xi)^2}{1+|\partial\psi|^2}\right)
\ge\frac F{s^3}|\xi|^2>0.
\tag{19}
$$

故在声明的正则分支，单字段非线性方程的主符号正定；不存在实的非零特征协向量。这是椭圆方程，不能仅把材料坐标 $T$ 政名为物理时间就变成双曲的波方程。

### 6.2 多剩余字段：平坦图像的线性化

先给出平坦图像附近的明确线性结论，再于第6.3节核对一般图像的 Legendre–Hadamard 主符号；不把后者混同于面积泛函对全部梯度矩阵的普通凸性。

令 $\psi^A=0$，并要求 $F_A(T,x,0)=0$，使该图像满足剩余字段的背景方程。写 $F_0=F(T,x,0)>0$、$F_{AB}=\partial_A\partial_BF|_0$。二次作用与线性方程分别为

$$
S^{(2)}
=\frac12\int d^4x
\left[F_0\sum_A\partial_\mu\eta^A\partial_\mu\eta^A
+\sum_{A,B}F_{AB}\eta^A\eta^B\right],
\tag{20}
$$

$$
\partial_\mu(F_0\,\partial_\mu\eta^A)
-\sum_BF_{AB}\eta^B=0.
\tag{21}
$$

主符号为 $F_0\delta^{\mu\nu}\delta_{AB}$，仍为正定。$F_0$ 的变化和 $F_{AB}$ 改变一阶及零阶项，不会把该裸主符号改成 Lorentz 符号。

例如原文势在此局部规约下为
$V=V_0'+m^2(T^2+|x|^2+\sum_A(\psi^A)^2)/2$，若 $V>0$，则 $F_A|_0=0$、$F_{AB}|_0=-m^2\delta_{AB}/V(T,x,0)^2$。此例仍满足上述主符号结论。但零阶项能够改变某些有限波数的具体行为；不能由主符号独自断言原固定有质量支的所有低频模式都指数增长。

### 6.3 一般多字段图像的主符号

仍限于实图像和正 $F$。记 $Z_{\mu A}=\partial_\mu\psi^A$，沿用式(15)的梯度矩阵，并置

$$
g=I_4+ZZ^T,\qquad Q=g^{-1}Z,\qquad
\mathsf N=I-Z^Tg^{-1}Z=(I+Z^TZ)^{-1}>0.
\tag{21a}
$$

对图像 Lagrangian 的梯度变量求两次导数，得到

$$
\frac{\partial^2\mathcal L_{\rm graph}}
{\partial Z_{\mu A}\partial Z_{\nu B}}
=F\sqrt{\det g}\left[
g^{\mu\nu}\mathsf N_{AB}
+Q_{\mu A}Q_{\nu B}-Q_{\mu B}Q_{\nu A}\right].
\tag{21b}
$$

$F$ 对梯度变量无依赖。最后两项与 $\xi_\mu\xi_\nu$ 收缩后相消，故系统主符号为

$$
\sigma_{AB}(\xi)
=F\sqrt{\det g}\,(\xi^Tg^{-1}\xi)\,\mathsf N_{AB}.
\tag{21c}
$$

对所有非零实协向量 $\xi$ 和非零剩余字段极化向量 $z$，

$$
z^T\sigma(\xi)z
=F\sqrt{\det g}\,(\xi^Tg^{-1}\xi)(z^T\mathsf Nz)>0.
\tag{21d}
$$

因此不仅平坦背景，一般正则实图像也满足该系统的强椭圆主符号条件。紧致正则任务域上的系数界可以给一致椭圆界；本轮不据此宣称把 $T$ 当时间的完整非线性初值问题已适定。式(21b)对一般梯度矩阵变化不必处处正定，故上述结论没有偷换为面积泛函的普通凸性。

## 7. 同一实际钟下的 Hamiltonian 交叉检查

单剩余字段时，式(16)给出

$$
p=\frac{F\psi_T}{s},\qquad
H_T=-\sqrt{(1+|\nabla\psi|^2)(F^2-p^2)}
\tag{22}
$$

于原正向分支。式(22)既来自直接 Legendre 变换，也等于式(13)，因为此时 $\det h=1+|\nabla\psi|^2$、$M=1+|\nabla\psi|^2$。

在常 $F$ 的平坦图像附近，二次展开是

$$
H_T=-F+\frac{p^2}{2F}
-\frac F2|\nabla\psi|^2+O\!\left((p,\nabla\psi)^4\right).
\tag{23}
$$

因此线性 Hamilton 方程给出
$\partial_T\psi=p/F$、$\partial_Tp=-F\Delta\psi$，即 $\psi_{TT}=-\Delta\psi$，与式(21)一致。Hamiltonian 并非因式(6)就没有演化；它产生的是所选裸分支本身的椭圆传播方程。

将完整作用乘以负号不改变 Euler–Lagrange 方程的零集合；把 $T$ 翻向也不改变二阶主符号的惯性指数。一般实可逆坐标变换对主符号作合同变换，保持正定或负定的类型。因此这类符号或坐标选择不能把式(19)改成 Lorentz 光锥。

这不等价于声称“出现增长就不能量子化”：正常量子系统也可有不稳定经典解。本轮只检验同一裸材料钟接口能否直接提供目标所需的稳定 Lorentz 波传播。

## 8. 固定有限波数与有限时窗的线性诊断

为得到不依赖无限频率外推的明确见证，选原作用族中的 $m=0$、$V=V_c>0$ 常势诊断支，并取平坦图像。该选择不是文献用于诱导 Planck 系数的固定有质量量子真空，也不把此诊断代替该真空。此时式(20)精确给出线性化作用

$$
S^{(2)}=\frac{F_c}{2}\int dT\,d^3x
\left[(\partial_T\eta)^2+|\nabla\eta|^2\right],
\qquad F_c=1/V_c,
\qquad \eta_{TT}+\Delta\eta=0.
\tag{24}
$$

选一个固定有限非零波向量 $\mathbf k$，$k=|\mathbf k|>0$。可在局部平移不变图像作模态诊断；若使用有限盒的 Fourier 规范化，则周期性只用于扰动的谱诊断，不声称四个全局材料坐标也周期。对线性初值

$$
\eta(0,x)=\epsilon\cos(\mathbf k\cdot x),\qquad
\eta_T(0,x)=0,
\tag{25}
$$

裸线性解与相应模态的推迟逆分别为

$$
\eta_{\rm raw}(T,x)
=\epsilon\cosh(kT)\cos(\mathbf k\cdot x),\qquad
G_{{\rm raw},R}(T,k)=\theta(T)\frac{\sinh(kT)}k.
\tag{26}
$$

$G_R$ 的含义仅为固定模态算子 $\partial_T^2-k^2$ 的带推迟初值逆，并非已经构造出具有有限因果域的完整椭圆场仪器。

对照采用同一 $T$、$x$、初值归一和单位传播系数的 Lorentz 波方程
$\eta_{TT}-\Delta\eta=0$，有

$$
\eta_{\rm L}(T,x)=\epsilon\cos(kT)\cos(\mathbf k\cdot x),\qquad
G_{{\rm L},R}(T,k)=\theta(T)\frac{\sin(kT)}k.
\tag{27}
$$

令有限时窗端点为 $T_*=1/k$。两个模型在原 Fourier 模式振幅上的线性响应导数之差为

$$
\left.\frac{d}{d\epsilon}
\left(A_{\rm raw}(T_*)-A_{\rm L}(T_*)\right)\right|_{\epsilon=0}
=\cosh1-\cos1
=1.002778328947104\ldots.
\tag{28}
$$

对于模态推迟逆，同一点亦有
$k[G_{{\rm raw},R}(T_*,k)-G_{{\rm L},R}(T_*,k)]=\sinh1-\sin1>0$。
这些均取固定有限 $k$ 和固定有限 $T_*$，没有使用 $k\to\infty$ 或无限等待。因此“差异只存在于无限高能外推”不能解释这个已声明常势支的线性见证。

需要严格区分：

1. 式(26)是**线性化方程的精确解**，不是完整非线性图像作用的精确有限振幅解。
2. 可选小 $\epsilon$ 控制形式上的扰动斜率，但本轮没有因此取得完整非线性解的存在时间或严格有限振幅误差界。
3. Fourier 振幅是数学响应诊断，尚不是已经内部实现的量子探针、准备和记录。
4. 常势支的见证没有自动转化为原有质量重整化真空的定量预测。对后者，本轮提供的是裸实际材料分支的主符号限制。

## 9. 可签结论与停止线

可签的条件结论是：

$$
\begin{gathered}
\text{正实行列式裸作用}
+\text{正则真实四材料坐标}
+\text{以该材料钟读取同一约化过程}
\\
\Longrightarrow\
\text{正则实图像的椭圆主符号及其平坦线性化};
\\
\text{再取常正势诊断支}
\Longrightarrow\
\text{固定有限波数、有限时窗的线性响应不等于 Lorentz 波响应}.
\end{gathered}
\tag{29}
$$

这补上了旧 Gram 审计没有单独完成的“实际钟—正则约化—同作用传播”连接。它排除了把这个裸正实分支直接接成既有 Lorentz 钟与信号接口的捷径。

它没有证明全部复合几何候选失败，没有反驳原文的重整化集体真空，也没有证明现行认知公理不能产生有效 GR。若候选改用集体钟、重整化复合量、不同的内部收缩符号或受控解析延拓，必须说明新的共同材料过程怎样改变有效主符号并保留实际参考、资源与记录；那是另一个有明确输入的机制问题，不能由辅助 $X$ 的坐标选择代替。

下一步不继续修饰辅助钟、不重复正 Gram 数值验证，也不把完整紫外理论当成准入要求。对本候选的后续投入，应由一个可检验的**同材料集体参考与有效双曲响应机制**决定。涉及真空来源时仍可与暗部门共同生成猜想衔接，但本轮未得到暗能量、暗物质或 Einstein 动力学。
