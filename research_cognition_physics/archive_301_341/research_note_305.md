# 第305轮：量子熵第一定律、有限余项与物理能量的识别

日期：2026-09-22。承接[304轮](research_note_304.md)，检查视界热力学哪些部分已经包含在前两阶段的量子状态结构中。本轮是有限维量子计算，不预装几何、视界、场论真空或Unruh温度。

## 1. 问题、已有工具与输入

**认知动机：** 子系统对整体的可访问范围受限时，其状态变化与熵变化之间能否提供物理生成所需的能量接口？

**待检验假设：** 量子相对熵给出普适的一阶熵关系及精确有限余项；要将其中的模生成元解释为指定物理能量，仍需可检验的相容条件。

复用[Casini，2008，§3—4](https://arxiv.org/abs/0804.2182)中的相对熵与模能量接口。模生成元和几何流的关系另参见[Witten的量子场论纠缠讲义](https://arxiv.org/abs/1803.04993)。这些是成熟结果，本轮的工作是有限矩阵复算、误差控制及与物理生成元的识别条件。

**输入：** 有限维复矩阵状态空间、归一化密度矩阵、von Neumann熵、满秩参考态σ，以及指定可实现的态扰动。对数为自然对数、k_B＝1。采用这个熵作为待比较信息量是明示选择；参考态、物理Hamiltonian和可访问扰动集合尚未由认知原则唯一决定。

## 2. 精确恒等式与一阶关系

定义无量纲模Hamiltonian及相对熵：

$$
K_\sigma=-\log\sigma,\qquad
S(\rho)=-\operatorname{Tr}\rho\log\rho,\qquad
D(\rho\Vert\sigma)=\operatorname{Tr}\rho(\log\rho-\log\sigma).
\tag{1}
$$

令Δ表示与同一参考σ比较，直接整理迹得到：

$$
D(\rho\Vert\sigma)=\Delta\langle K_\sigma\rangle-\Delta S\ge0.
\tag{2}
$$

其中K始终固定为参考态的生成元。不能同时将它换成新态的−logρ再套用同一差分公式。熵和相对熵的定义及非负性不要求有时空几何。

对归一化可微路径\(\rho_\epsilon=\sigma+\epsilon X\)，\(X=X^\dagger\)、TrX＝0，且在σ的正定邻域内，D在ε＝0取最小值，因此：

$$
\left.\frac{dS(\rho_\epsilon)}{d\epsilon}\right|_0
=\operatorname{Tr}(XK_\sigma),\qquad
\Delta S=\epsilon\operatorname{Tr}(XK_\sigma)-D(\rho_\epsilon\Vert\sigma).
\tag{3}
$$

第一式就是纠缠第一定律的有限维形式；有限变化的等号缺口恰为相对熵。它不要求X与σ交换。

## 3. 二阶项与有限误差界

在σ的本征基中记本征值为pᵢ，矩阵对数的Fréchet导数给出：

$$
D(\sigma+\epsilon X\Vert\sigma)
=\frac{\epsilon^2}{2}\sum_{ij}\lvert X_{ij}\rvert^2 L_{ij}+O(\epsilon^3),\qquad
L_{ij}=\frac{\log p_i-\log p_j}{p_i-p_j},\quad
L_{ij}=\frac1{p_i}\ \text{当}\ p_i=p_j.
\tag{4}
$$

这是对矩阵函数的导数，不把非交换问题退化成对角概率变化。若整条线段的最小本征值下界为\(\mu=p_{\min}(\sigma)-\lvert\epsilon\rvert\lVert X\rVert_{\rm op}>0\)，则：

$$
\frac12\lVert\epsilon X\rVert_1^2
\le D(\sigma+\epsilon X\Vert\sigma)
\le\frac{\epsilon^2\lVert X\rVert_F^2}{2\mu}.
\tag{5}
$$

左边用自然对数的量子Pinsker界；右边用对数导数的正算子积分，将沿路径的二阶导数控制在\(\lVert X\rVert_F^2/\mu\)内，再积分两次。二阶界依赖参考谱；当最小本征值趋于0时，它不提供统一控制。

数值例取H＝diag(0,1,2,4)、σ∝exp(−0.7H)，再选含非对角项的迹零X，保证所有测试态正定。σ与X的交换子范数约0.00411747。

| ε | 相对熵D | D／ε² |
|---:|---:|---:|
| 0.8 | 0.000434378598934 | 0.000678716561 |
| 0.4 | 0.000108193493336 | 0.000676209333 |
| 0.2 | 0.000027000588527 | 0.000675014713 |
| 0.1 | 0.000006744317869 | 0.000674431787 |

二阶解析极限为0.000673858319；逐次减半的D比值为4.0148、4.0071、4.0035。精确恒等式误差约1.52×10⁻¹⁶以内，有限上下界均满足。纯相干扰动还给出ΔK＝0而ΔS＝−D＜0的例子：一阶能量变化为零，不代表有限熵变化也为零。

## 4. 模能量何时等于指定物理能量

给定一个非恒等的物理Hamiltonian H和候选逆温度β。若要求**所有迹零Hermitian方向**的一阶熵变化都满足物理热关系，则：

$$
\operatorname{Tr}[X(K_\sigma-\beta H)]=0\quad\forall X=X^\dagger,\ \operatorname{Tr}X=0
\quad\Longleftrightarrow\quad
K_\sigma=\beta H+cI
\quad\Longleftrightarrow\quad
\sigma=\frac{e^{-\beta H}}{\operatorname{Tr}e^{-\beta H}}.
\tag{6}
$$

证明：满秩保证每个这样的X在足够小邻域内对应合法态变化；迹零Hermitian空间的正交补恰为I的线性张成。指数化并归一化完成最后一步。

这给出可操作的识别条件，而非将“任何σ能写成exp(−K)”直接解释成真实温度。对K、H去掉单位阵部分后，最小二乘拟合β，完整残差应为零。正例恢复β＝0.7，残差3.19×10⁻¹⁶以内。

反例采用同一H和\(\sigma_b=\operatorname{diag}(0.45,0.30,0.20,0.05)\)。两者交换，状态在H下平稳，却不满足单一温度的Gibbs关系：最佳β约0.5534159，残差约0.2382234。平稳性不能代替温度识别。

在ℏ＝1的数值约定下，正例还满足模流与物理流重参数化一致：

$$
e^{-iuK_\sigma}Oe^{iuK_\sigma}
=e^{-i\beta uH}Oe^{i\beta uH}.
\tag{7}
$$

实际物理时间和单位须由物理H及其标定给出；式(7)本身不推出微观时间性质。

## 5. 可访问方向不足的反例

若只允许对角准备变化，式(6)只约束对角投影。构造：

$$
\sigma_r=Z^{-1}e^{-[0.7H+0.3X_{01}]},\qquad
X_{01}=\lvert0\rangle\langle1\rvert+\lvert1\rangle\langle0\rvert.
\tag{8}
$$

K−0.7H的非平凡部分纯非对角，故所有对角迹零方向都通过β＝0.7的测试，完整残差却为√0.18≈0.424264。用相干方向X₀₁即可得到0.6的非零见证。这里讨论的是指定准备集合的不足；若还能观察非平稳演化，也可能另行排除该候选。

零本征值参考态的−logσ在全空间发散；代码拒绝用本轮满秩定理处理它，而不静默添加小数正则化。恒等Hamiltonian也无法从此接口唯一标定β。

## 6. 接向视界，还差哪一步

在适当的局域、正能量、Poincaré协变量子场论真空和Rindler楔代数条件下，Bisognano–Wichmann结果给出几何模流。其通常写法是：

$$
K_R=\frac{2\pi}{\hbar}\int_{x>0}x\,T_{00}\,d^3x+cI.
\tag{9}
$$

该几何识别使用额外QFT、真空及区域代数条件；不是有限矩阵恒等式的直接推论。连续场的局域代数通常不能直接照搬有限密度矩阵的全空间迹，需用适当代数定义或有控制的截断。本轮未证明有限矩阵序列已有该连续极限。

**已经补上的接口：** 指定熵与满秩参考以后，一阶熵关系及精确相对熵余项无需新增“热力学公理”。**仍需补上：** 参考态来源、物理能量／几何模流的对应，以及面积熵与几何变化的联系。

后者由[306轮](research_note_306.md)先对接Jacobson的小球方法，核验固定体积的面积响应。9项检查，代码不渲染图像。

[代码](305/modular_energy_interface_audit.py)、[结果](305/modular_energy_interface_audit_results.json)、[核验](305/research_round_305_checks.json)。

    python -B -X utf8 research_cognition_physics/archive_231_/modular_energy_interface_audit.py
