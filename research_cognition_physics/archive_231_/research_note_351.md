# 第351轮：超曲面形变能否选择 Einstein 动力学——纯度规约束括号的条件性定理

## 0. 本轮问题与结论

本轮继承截至348轮的冻结基线，主要接续[第344轮的协变引力唯一性审计](research_note_344.md)，不依赖本批其他未完成轮次。先读项目 README、本研究目录的 README、research_direction、RESEARCH_STATE 及既有共同度规／约束结果；未发现已经完成本轮纯度规超曲面形变括号的推导。

**有明确的正向结果：** 在下列已给定的局域纯度规正则模型族里，要求任意局域 lapse 的约束括号在整个相空间上精确表示 Lorentz 型超曲面形变，确实唯一选出 Einstein 动能的迹系数。在三维空间中，本轮 Hamiltonian 的系数是 λ＝1/2，文献常用的 Lagrangian 系数则是 λ_L＝1；两者不能混用。

**这不是“操作可以组合，因此完整 GR 已推出”。** 认知操作的组合律尚未给出平滑空间度规、纯度规正则相空间、指定导数阶的 Hamiltonian 或超曲面形变代数。本轮将这条既有研究路线变成可复算的条件接口，并给出每一项额外输入的账目。它也没有证明所有 λ≠1/2 的理论不一致。

| 层次 | 本轮处理 |
|---|---|
| 认知动机 | 同一整体的描述不应依赖无可观测意义的处理顺序 |
| 额外建模输入 | 平滑正定空间度规、局域正则变量、特定动能／曲率势、任意局域 lapse、Lorentz 型形变括号 |
| 解析证明 | 完整计算异常项；在指定模型族内证明迹系数的必要性与充分性 |
| 数值复算 | Fourier 数据；非线性曲率作用量的独立方向差分；13项检查 |
| 物理解释 | 条件性恢复 ADM／Einstein 纯引力；不能将条件解释为已由认知生成 |

## 1. 对接原始文献：两种“路径独立”不是同一个条件

[Teitelboim（1973）](https://doi.org/10.1016/0003-4916(73)90096-1)研究的是：给定初始与最终空间超曲面，经过不同中间超曲面的形变应给出一致的演化。这已经涉及超曲面、法向与切向、诱导度规以及可嵌入时空的几何要求。两个法向形变的交换一般产生切向形变，不是要求它们彼此对易。

[Hojman、Kuchař、Teitelboim（1976），《Geometrodynamics regained》](https://doi.org/10.1016/0003-4916(76)90112-3)是从此接口重建纯度规动力学的原始入口。本轮不冒充复证完整 HKT 唯一性定理：其出版社全文本轮访问受限，下面的结论由明确受限的 ansatz 独立证明，并用可读取的 λ-R 原始论文核对符号与例外。

“先做操作 A 再做 B 的复合有良好定义”“不同算法实现同一个信道”均不足以自动给出下面的度规依赖括号。尤其是：局部量子操作的组合本身，没有指定一个代表空间切片的 q，也没有指定哪些生成元应解释为法向变形。

## 2. 额外输入、记号与待满足的代数

取维数已给定的紧致无边界空间 Σ，d≥2；或采用足以令全部分部积分边界项消失的支撑条件。q_ij 是光滑正定空间度规，q 是其行列式，π^ij 是唯一共轭动量、权重为1的对称张量密度，π 是其度规迹。采用

$$
\begin{aligned}
\{F,G\}
&=\int_\Sigma d^dx\left(
\frac{\delta F}{\delta q_{ij}}\frac{\delta G}{\delta\pi^{ij}}
-\frac{\delta F}{\delta\pi^{ij}}\frac{\delta G}{\delta q_{ij}}
\right),\\
\pi&=q_{ij}\pi^{ij},\qquad
\pi_{ij}=q_{ik}q_{jl}\pi^{kl},\\
D[v]&=\int_\Sigma d^dx\,\pi^{ij}{\cal L}_v q_{ij}
=-2\int_\Sigma d^dx\,v^i q_{ij}\nabla_k\pi^{jk}.
\end{aligned}
\tag{1}
$$

整个模型族预先限制为：局域、动量二次、动能中无空间导数、势仅含标量曲率与空间常数项。暂设两个系数的单位归一化为1：

$$
H_\lambda[N]=
\int_\Sigma d^dx\,N\left[
\frac{\pi^{ij}\pi_{ij}-\lambda\pi^2}{\sqrt q}
-\sqrt q\,R+2\Lambda\sqrt q
\right].
\tag{2}
$$

这里 λ 乘在 Hamiltonian 的迹平方上。N 是给定的任意光滑 smear，即 lapse；Λ 是空间常数。两者不依赖相空间变量。本轮未从认知原则证明式(2)的函数形式，也未加入物质或额外引力场。

待要求的 Lorentz 型形变代数是

$$
\begin{aligned}
\{D[v],D[w]\}&=D[[v,w]],\\
\{D[v],H_\lambda[N]\}&=H_\lambda[{\cal L}_vN],\\
\{H_\lambda[N],H_\lambda[M]\}
&=D[v_{N,M}],\qquad
v_{N,M}^{\,i}=q^{ij}(N\partial_jM-M\partial_jN).
\end{aligned}
\tag{3}
$$

右侧含 q^ij，因此其“结构系数”实际上依赖度规。前两式由空间微分同胚的正则生成与 Hamiltonian 密度的变换性质满足，不挑选 λ。第三式是关键。

**本轮选择条件是强恒等式：** 式(3)最后一式须对整个未约束的相空间、任意 N、M 成立，不能只在某组解上成立，也不能在右侧另加尚未申明的迹约束。它比“某个扩展的 Dirac 约束算法最终一致”更强。

## 3. 解析推导：多出来的正是迹梯度

动量变分与曲率势中含 lapse 导数的度规变分分别为

$$
\begin{aligned}
\frac{\delta H_\lambda[N]}{\delta\pi^{ij}}
&=\frac{2N}{\sqrt q}
\left(\pi_{ij}-\lambda\pi q_{ij}\right),\\
\left.
\frac{\delta\left(-\int_\Sigma N\sqrt q\,R\,d^dx\right)}
{\delta q_{ij}}\right|_{\nabla N}
&=\sqrt q\left(q^{ij}\Delta N-\nabla^i\nabla^jN\right),\\
\Delta&=q^{ij}\nabla_i\nabla_j.
\end{aligned}
\tag{4}
$$

动能、Λ项及曲率变分的其余部分都局域正比于 N；与另一动量变分相乘后含 NM，反对称化即消去。势与势的括号本来为零；无导数动能之间的括号也反对称消去。因此式(4)已包含全部非消去项，而不是只保留线性近似。

将其代入式(1)，收缩指标可得

$$
\begin{aligned}
C_d(\lambda)&=1-(d-1)\lambda,\\
\{H_\lambda[N],H_\lambda[M]\}
&=2\int_\Sigma d^dx\,
\Big[
N\pi^{ij}\nabla_i\nabla_jM
-M\pi^{ij}\nabla_i\nabla_jN\\
&\hspace{38mm}
+C_d(\lambda)\pi(M\Delta N-N\Delta M)
\Big].
\end{aligned}
\tag{5}
$$

对称 π^ij 使含两个一阶 lapse 梯度的反对称项消去。再分部积分，得到本轮的完整括号：

$$
\boxed{
\{H_\lambda[N],H_\lambda[M]\}
=D[v_{N,M}]
+2C_d(\lambda)\int_\Sigma d^dx\,
v_{N,M}^{\,i}\nabla_i\pi
}.
\tag{6}
$$

这里不可把 π 当普通标量。它是标量密度，所以

$$
\nabla_i\pi
=\partial_i\pi-\Gamma^k{}_{ki}\pi
=\sqrt q\,\partial_i\!\left(\frac{\pi}{\sqrt q}\right).
\tag{7}
$$

式(6)直接证明充分性：若 C_d＝0，第三个形变括号对所有相空间数据成立。Λ 没有出现在异常项中，也没有出现在普通动量项的系数中。

## 4. 必要性不是扫描猜测：一个散度为零的精确见证

取平坦 d 维环面，所有坐标周期2π。记归一化体积平均为尖括号，并对所有泛函与正则配对一致采用此归一化。构造

$$
\begin{aligned}
q_{ij}&=\delta_{ij},&
\pi^{22}&=\sin x^1,&
\pi^{ij}&=0\quad (ij\ne22),\\
N&=1,& M&=2+\sin x^1,&
\nabla_j\pi^{ij}&=0,\\
\pi&=\sin x^1,&
v_{N,M}^{\,1}&=\cos x^1,&
D[v_{N,M}]&=0.
\end{aligned}
\tag{8}
$$

两个 lapse 都严格为正；见证不依赖允许负 lapse。动量张量的散度为零，但其迹沿第一个坐标变化，所以

$$
\{H_\lambda[N],H_\lambda[M]\}
=2C_d(\lambda)\langle\cos^2x^1\rangle
=C_d(\lambda).
\tag{9}
$$

于是，在已经限定的模型族中，

$$
\boxed{
\text{式(3)对所有相空间数据及局域 lapse 强成立}
\quad\Longleftrightarrow\quad
\lambda=\frac1{d-1}.
}
\tag{10}
$$

这是真正的必要性反例，而不是有限参数搜索给出的推测。**见证是离壳数据**：它一般不满足 Hamiltonian 约束。它正确地反驳错误的强恒等式，但不能单凭此断言某个增加了约束、固定了 lapse 的理论没有一致解。

在 d＝2、3、4 中分别得到 λ＝1、1/2、1/3。这个定理对每个已给定的 d 成立，因而**没有选出 d＝3**。

## 5. 它如何有条件地回到 Einstein 动力学

对式(10)选出的系数，使用总 Hamiltonian H[N]＋D[β]，定义切片的外挠曲率，并作逆 Legendre 变换：

$$
\begin{aligned}
K_{ij}
&=\frac{\dot q_{ij}-{\cal L}_\beta q_{ij}}{2N}
=\frac1{\sqrt q}\left(\pi_{ij}-\frac{\pi}{d-1}q_{ij}\right),\\
\pi^{ij}&=\sqrt q\left(K^{ij}-Kq^{ij}\right),\\
S_{\rm ADM}
&=\int dt\,d^dx\,N\sqrt q
\left(K_{ij}K^{ij}-K^2+R-2\Lambda\right).
\end{aligned}
\tag{11}
$$

在给定 Lorentz 时空重建、边界约定和本轮单位后，这就是 Einstein–Hilbert 作用量的 ADM 形式，相差相应边界项；变分得到含 Λ 的真空 Einstein 方程。这一步表明本轮不只是得到一个抽象数值 λ，而是接上已知纯度规 GR。

不过，把式(11)认作某个时空的 ADM 分解，本来就在式(3)的几何解释之内。此处没有由认知生成连续时空，也没有从量子操作证明经典度规的唯一性，更未选定 matter action。

### 5.1 文献 λ 约定的核对

[Loll–Pires（2014）](https://arxiv.org/html/1407.1259)的式(7)、(9)分别使用外挠曲率动能系数及其逆超度规系数。用不同字母可避免把该文“λ＝1”错套进本轮 Hamiltonian：

$$
\begin{aligned}
{\cal L}_{\rm kin}
&=N\sqrt q\left(K_{ij}K^{ij}-\lambda_LK^2\right),\\
\pi^{ij}&=\sqrt q\left(K^{ij}-\lambda_LKq^{ij}\right),\\
\lambda&=\frac{\lambda_L}{d\lambda_L-1},
\qquad
\lambda_L=\frac{\lambda}{d\lambda-1},\\
\lambda_L=1&\quad\Longleftrightarrow\quad
\lambda=\frac1{d-1}.
\end{aligned}
\tag{12}
$$

转换只在对应 Legendre 映射可逆时成立。d＝3 时，式(5)的异常部分转换成该文式(30)的系数和 lapse 次序，符号一致。

## 6. 几个必须保留的分支

### 6.1 常平均曲率及新增约束

若额外限制

$$
\nabla_i\pi=0
\quad\Longleftrightarrow\quad
\pi=a(t)\sqrt q
\quad\text{（每个连通空间切片）},
\tag{13}
$$

则式(6)的异常在这个子空间上消失，并不需要式(10)。这不是对强恒等式的反例，而是换成了附加约束／规范条件。保持该条件可进一步限制 lapse，并改变约束的第一类／第二类性质；不能据此认为任意 lapse 下的原代数已在整个相空间成立。

[Bellorín–Restuccia（2010）](https://arxiv.org/html/1004.0055)在其渐近平坦 λ-R 分析中给出了迹约束与相应 GR 规范描述的联系。[Loll–Pires](https://arxiv.org/html/1407.1259)对闭空间明确区分 π＝0 与一般 π＝a(t)√q，以及 projectable 部门。因此本轮不能把 λ≠GR 值一概称为“不一致”，也不能把全部边界条件下的 λ-R 都宣称已等价于 GR。

### 6.2 Projectable lapse

若只准 N＝N(t)、M＝M(t)，则空间梯度为零，v＝0，两侧括号都为零，对 λ 不产生选择。它放弃的是任意局部法向形变的要求；是否构成可用理论要另作约束与物理分析。

### 6.3 曲率／动能系数、Euclidean 与超局域分支

进一步允许空间常数 a、b、u，并将式(2)推广为

$$
\begin{aligned}
H_{\lambda,a,b,u}[N]
&=\int_\Sigma d^dx\,N
\left[
\frac{a}{\sqrt q}(\pi^{ij}\pi_{ij}-\lambda\pi^2)
-b\sqrt q\,R+2u\sqrt q
\right],\\
\{H[N],H[M]\}
&=ab\left[D[v_{N,M}]
+2C_d(\lambda)\int_\Sigma d^dx\,v_{N,M}^{\,i}\nabla_i\pi\right].
\end{aligned}
\tag{14}
$$

在相同归一化下，指定式(3)的正号要求 ab＝1，并由式(8)继续选出式(10)。a＝2、b＝1/2 与 a＝b＝1 具有同一个括号；单凭此不能确定 Newton 常数或物质相对归一化。

若指定相反号的形变代数，则 ab＝−1 对应相反的签名分支；其异常消除条件相同。本轮已输入 Lorentz 符号，不能再将其当成输出。若 b＝0、a≠0，则势无空间导数，两个 H 的括号对任意 λ 都为零，是超局域分支；a＝0 的退化势模型也使括号为零。零积不实现本轮指定的非平凡式(3)，但不能仅凭这一点说它们在自己的动力学意义下“不一致”。两种非零系数同时反号也保持 ab，额外的能量／时间取向要求不由这个括号自动决定。

### 6.4 两种不同的退化点与 d＝1

Lagrangian 在 λ_L＝1/d 时迹方向的 Legendre 映射奇异，通常必须从新的初级约束开始分析；不能用式(12)除以零。另一方面，直接给出的 Hamiltonian 动能在 λ＝1/d 时，对动量的迹 Hessian 有零特征值，是另一退化点，也不是前一个点经可逆转换得到的有限像。本轮的括号公式仍可评价后者，见证给出异常1/d；但这不替代其完整约束分类。

d＝1 的内禀标量曲率恒为零，且式(10)分母为零，不在本轮定理适用范围。

### 6.5 额外引力自由度

第344轮的健康 R＋αR² 模型可在适当分支转成包含额外标量及其共轭动量的正则理论。对完整变量和全部物质一致处理时，它仍可具有共同的时空形变代数；标量提供的约束贡献不能删掉。

因此它不满足本轮“唯一正则变量就是 q 与 π”的选择输入，不与式(10)矛盾。本轮也不能借“纯度规”把这个输入偷偷证明掉。排除额外引力自由度仍是额外结构选择；若它们在低能未被激发、可被有效积分掉，GR 作为有效低能近似是否涌现仍须独立研究。

## 7. 可复算核验与实际结果

代码采用 NumPy，不依赖符号软件。不是把式(6)编码两遍互相印证，而有三个层次：

1. 用完整的未积分梯度收缩计算 Poisson 括号，包含动能、曲率 lapse 二阶导数和 Λ项；与分部积分后的动量项＋迹异常比较。
2. 独立从 q、Christoffel、Ricci 计算非线性 Hamiltonian 泛函，在正定的 q±εh 上取真实作用量差分。该路径不调用异常公式；非平坦共形度规的 Ricci 还与解析共形恒等式交叉检查。
3. 用式(8)的散度零见证、一般对角／非对角 Fourier 动量、多个维数、两种网格及不同分支复算。

差分针对的正则方向为

$$
\begin{aligned}
A_N&=\left.\frac{\delta H[N]}{\delta\pi}\right|_{q=I},\\
B_\varepsilon
&=\frac{H[N](I+\varepsilon A_M,\pi)-H[N](I-\varepsilon A_M,\pi)}
{2\varepsilon}\\
&\quad-\frac{H[M](I+\varepsilon A_N,\pi)-H[M](I-\varepsilon A_N,\pi)}
{2\varepsilon},\\
B_{\rm extrap}
&=\frac{4B_{\varepsilon/2}-B_\varepsilon}{3}.
\end{aligned}
\tag{15}
$$

本轮一般 Fourier 数据在 d＝3、λ＝0.17 时，解析括号为−0.1840489。ε＝0.002、0.001、0.0005 的差分绝对误差依次约为1.129×10⁻⁷、2.822×10⁻⁸、7.055×10⁻⁹，按二阶差分预期缩小；外推误差约1.2×10⁻¹⁴至4.0×10⁻¹⁴。式(8)的异常同样由实际曲率作用量差分独立恢复。

在 GR 迹系数下，一般数据的动量项为−0.066955。Λ＝−3、0、0.7、12得到同一个括号，误差在舍入范围；a＝2、b＝1/2保持此值，a＝1、b＝−1反号，零积给零。后者验证的是符号／归一化边界，不是增加物理模型的可接受性证明。

13项检查全部通过：谱导数／分部积分；非平坦 Ricci；完整作用量变分；梯度与积分括号；必要性见证；差分收敛；非零动量项下的充分性；Λ自由；ab签名；超局域／projectable；CMC子空间；两种 λ 约定及退化；独立作用量与网格。数值只覆盖一个坐标依赖的张量配置，任意维平滑场的结论来自第3—4节解析恒等式，不能由这些样本取代。

## 8. 本轮完成与下一步

关闭的问题是：**一个明确的几何可组合性要求，能否在已限制的纯度规模型族中实际选出 Einstein 的动能？能。** 代价也清楚：这里的可组合性是带度规结构函数、任意 lapse、离壳成立的超曲面形变代数。

仍未关闭的桥接包括：为何认知操作应有这种形变解释；为何有效引力仅含度规正则变量；为何该低阶局域 ansatz 足够；为何三维、Lorentz 符号、特定 Λ 和物质参数成立。一般协变并不单独排除第344轮的额外自由度；本轮的纯度规限制正是新增选择条件之一。

因此下一步应检验具体涌现模型能否提供这组生成元和约束，或把额外场／高阶项的低能消退给出可控制论证。不能在操作层只有组合律时，就把式(3)写成认知定理。物质部门的共同形变要求应单独完整推导，本轮不提前假定它已经完成。G、Λ及初边值可以作为GR方程族的自由数据；其具体数值尚未计算，不单独构成本项目“推导GR方程族”的未完成理由。

[代码](hypersurface_constraint_selection_audit.py)、[结果](hypersurface_constraint_selection_audit_results.json)、[本轮核验](research_round_351_checks.json)。

    python -B -X utf8 research_cognition_physics/archive_231_/hypersurface_constraint_selection_audit.py

保存结果使用同一命令加 --write-results；先通过检查，再保存，拒绝覆盖不同的既有结果。无图像检查。
