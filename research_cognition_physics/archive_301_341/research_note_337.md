# 第337轮：耦合扰动的真实光锥与稳定范围

日期：2026-09-23。接续[336轮](research_note_336.md)。本轮12项，第三阶段累计971项科学检查。

## 1. 问题、输入及已有工作

认知动机仍是：共同测量几何需要能实际传播信息的模式共同趋同。336轮只证明了两份指定度规的差异趋零，没有把非零物质背景中的动能混合算进去。本轮检验该缺口；不增加新的认知公理。

输入沿用336轮：Einstein引力、正则φ、通过g_B＝g＋D(φ)dφdφ耦合的无势χ、3＋1维平直FLRW均匀背景、正D和φ̇＞0。D＝0.6exp(√6φ)及初态区域仍是指定模型。当前计算是经典线性扰动的最高二阶导数部分，即决定局部传播特征和高频梯度稳定性的“主部”。

采用[van de Bruck、Koivisto、Longden（2016）§2.2.2与§3.1，式39—45](https://arxiv.org/html/1510.01650)的现成扰动接口，独立变分复算。该文的C是共形因子，须取为1；本笔记的C＝1−q是另一符号。文献已给出一般声速，本轮贡献是核对336吸引区域、统一误差界及边界反例，不把已有声速公式称为新发现。

## 2. 从完整标量作用量核对主部

在引力度规的局部惯性系，记X＝(∂φ)²、Y＝(∂χ)²、Z＝∂φ·∂χ。由逆度规及行列式直接得到标量拉氏密度

$$
\mathcal L=-\frac X2-\frac{\sqrt{1+DX}}2Y+
\frac{DZ^2}{2\sqrt{1+DX}},\qquad
v=\dot\phi,\quad w=\dot\chi,\quad q=Dv^2,\quad C=1-q,\quad
r=D\rho_\chi=\frac{Dw^2}{2C^{3/2}}=\frac{fq}{2(1-f)}.
\tag{1}
$$

对u＝(δφ,δχ)ᵀ展开二阶主部，令A＝Dvw/C：

$$
\mathcal L_{\rm prin}^{(2)}
=\frac12\dot u^{\mathsf T}K\dot u-\frac12(\nabla u)^{\mathsf T}G(\nabla u),
\quad
K=\begin{pmatrix}1+r+3rq/C&A/\sqrt C\\A/\sqrt C&1/\sqrt C\end{pmatrix},
\quad
G=\begin{pmatrix}1+r(2q-1)&A\sqrt C\\A\sqrt C&\sqrt C\end{pmatrix}.
\tag{2}
$$

FLRW坐标中空间导数需另除a，整体作用量乘a³；上述速度按引力局部钟尺计。代码同时对完整式(1)的速度、空间梯度作独立有限差分Hessian，核对式(2)，不是只比较两份同公式。

本模型物质作用量不含度规导数。规范固定后的Einstein方程主部与标量主部分块；标量引力约束消去不改最高导数系数，亦与所引文献式39—40一致。D的导数、H、背景加速度和Ȧ会影响较低导数项，本轮没有将这些项宣称为零。

## 3. 同时配方及真实特征速度

冻结背景系数，作可逆变量δθ＝δχ＋Aδφ，便有

$$
\dot u^{\mathsf T}K\dot u=(1+r/C)\delta\dot\phi^2+C^{-1/2}\delta\dot\theta^2,
\qquad
(\nabla u)^{\mathsf T}G(\nabla u)
=(1-r)|\nabla\delta\phi|^2+\sqrt C|\nabla\delta\theta|^2.
\tag{3}
$$

因此广义本征问题G e＝c²K e给出

$$
c_\theta^2=C,\qquad c_\phi^2=\frac{1-r}{1+r/C},\qquad
K>0\quad(C>0,\ D\ge0),\qquad
G>0\Longleftrightarrow r<1.
\tag{4}
$$

这说明物质占据也改变φ模式的传播。K正定表示没有负时间动能模；G正定排除该均匀背景上的高频梯度失稳。r＝1是零梯度边界，不能称为严格正定。时间变化的变量变换只对主部同时对角化，不意味着完整扰动方程解耦。

## 4. 与336的吸引证明拼接

在336的前向不变区域0≤f≤1/2、0≤q≤0.4中，r≤q/2≤0.2，C≥0.6。直接相减可得

$$
c_\phi^2-C=\frac{q-2r}{1+r/C}\ge0,\qquad
0.6\le c_\theta^2\le c_\phi^2\le1,\qquad
0\le1-c_\phi^2=\frac{r(1+C)}{C+r}\le q.
\tag{5}
$$

Einstein张量模式速度为1。故336的解析q(N)→0同时控制两条真实标量特征锥相对引力锥的偏差，并保持全程严格正的主部。不是从“度规看起来相同”直接跳到结论。

原变量的混合也消失：

$$
A^2=\frac{2rq}{\sqrt C}\le\frac{q^2}{\sqrt C},
\qquad |A|\le qC^{-1/4},\qquad
q\longrightarrow0\Longrightarrow K,G\longrightarrow I
\quad\text{在上述受控区域内}.
\tag{6}
$$

## 5. 数值复算与边界反例

| 初始f₀ | N | 当时q | c_φ² | c_θ² |
|---|---:|---:|---:|---:|
| 0.1 | 0 | 0.4 | 0.94285714 | 0.6 |
| 0.1 | 100 | 2.79476×10⁻⁷ | 0.99999998649 | 0.99999972052 |
| 0.01 | 100 | 0.08392940 | 0.99953272 | 0.91607060 |

N＝log a，不是固有时间秒数。复用336轨迹，没有重复运行其全套背景研究。f₀＝0.5初始两个标量均有c²＝0.6，仍不同于引力张量的1；两个标量彼此等速不等于整体共同几何。

完整作用量Hessian核验中，动能矩阵误差约2.50×10⁻⁸、空间矩阵误差约2.53×10⁻⁹，步长减半改善误差；独立Cholesky广义本征计算与式(4)一致。

两个反例限制推广范围：

- q＝0.4、f＝0.9时，两份度规仍为Lorentz型、K仍正，但r＝1.8，c_φ²＝−0.2，存在高频梯度失稳。因此度规正则不足以保证传播稳定。
- 取局部状态族r＝1/2、f＝1/(1＋q)，则

$$
q\to0,\quad f\to1,\qquad c_\theta^2\to1,\qquad c_\phi^2\to\frac13.
\tag{7}
$$

后一例仅是局部状态的代数反例族，未声称它是本模型的一条完整未来轨道。它证明q趋零单独不够；336中对能量比例的控制确实承担了证明责任。

## 6. 本轮结论与下一步

**在指定作用量及336初态区域内，度规吸引提升为真实标量主部的稳定共同传播极限。** 已补上动能混合这一明确缺口。

尚未证明低频无增长、非线性稳定、任意背景、全物种普适性或量子有效理论的适用范围。“高频主部”是该经典方程的特征分析；任意高频实际可用性及UV截止尺度仍需独立审计，不能由此宣称完整操作等价。作用量和耦合函数也尚未从认知原则选出。

下一轮检查对耦合函数参数的稳定性：β必须精确等于√6，还是存在参数开区间也会趋同；同时保留反作用和本轮主部稳定条件。其余输入按[GR依赖清单](333/gr_assumption_dependency_review.md)追踪。

[代码](337/disformal_characteristics_audit.py)、[结果](337/disformal_characteristics_audit_results.json)、[核验](337/research_round_337_checks.json)。

    python -B -X utf8 research_cognition_physics/archive_231_/disformal_characteristics_audit.py
