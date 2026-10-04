# 第354轮：集体量子自由度怎样给出经典 Poisson 流——受控极限、时间窗口与几何边界

日期：2026-09-23。独立完整研究；先读取项目 README、研究目录三个导航及已完成[351轮](research_note_351.md)，不依赖未冻结353轮。代码14项检查，无图像核验。

## 1. 问题、原始文献与本轮增量

351轮给出了“已有度规正则相空间＋指定超曲面形变代数”如何条件性选择 Einstein 动能。本轮往前检验一个不同接口：**有限量子系统内部的集体变量，能否在明确误差范围内给出经典 Poisson 相空间和动力学？若能，得到的是不是空间度规？**

采用成熟的集体自旋路线，而不重新创造经典极限理论：

- [Kitagawa–Ueda（1993），《Squeezed spin states》原文](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevA.47.5138/fulltext)第II节给出自旋相干态的同向自旋乘积表示，第III节及附录计算单轴扭转的矩。其赤道初态平均自旋公式与本轮结果在参数换算后相符。本轮独立用单个自旋的环境因子推出任意极角公式。
- [van de Ven（2020）的原始论文](https://arxiv.org/pdf/2007.03390)将对称量子自旋部门与球面的严格变形量子化连接起来。其一般量子化、谱及特定态极限结果，不能直接冒充本轮的有限时间误差定理；本轮另行证明所用动力学子类。
- [Lieb（1973）的作者机构原始摘要](https://collaborate.princeton.edu/en/publications/the-classical-limit-of-quantum-spin-systems/)说明该文主要控制配分函数、自由能及基态能量的经典极限。本轮不把这个静态结果误引成所有动态／所有输入态的收敛保证。

**结论：** 对指定的同向自旋乘积初态和单轴扭转 Hamiltonian，归一化集体均值在固定时间区间以 O(1/N) 逼近球面 Poisson 流，集体方差也为 O(1/N)。但在时间随 √N 增长时，均值与集中性可以发生有限偏离。更强的边界是：即使固定时间，完整量子态也未必接近沿经典轨道移动的相干乘积态。

这给出真实的经典动力学接口；**球面是内部自旋相空间，尚不是物理空间，更不是 q_ij 的相空间。**

| 层次 | 本轮内容 |
|---|---|
| 认知动机 | 大量内部自由度能否支持稳定、可预测的集体描述 |
| 额外输入 | N个已标记二能级系统、集体 Pauli 代数、相干乘积初态、特定耦合及时间单位 |
| 解析证明 | 精确矩、误差／集中界、√N窗口失效、局部正则坐标与全态边界 |
| 数值验证 | 对称部门矩阵、完整张量积矩阵两条实现，固定时间与增长时间对照 |
| 物理解释 | 指定态族／集体观测的经典近似，不是完整 GR 或所有量子通道的经典化 |

## 2. 量子模型与每项新增选择

设 ℏ＝1，N个二能级寄存器的 Pauli 算子为 σ_i^(a)。定义集体自旋及无量纲平均：

$$
J_i=\frac12\sum_{a=1}^{N}\sigma_i^{(a)},\qquad
\widehat s_i=\frac{2J_i}{N},\qquad
[\widehat s_i,\widehat s_j]
=\frac{2i}{N}\epsilon_{ijk}\widehat s_k.
\tag{1}
$$

N趋大时这些平均量的交换子范数趋零，但保留缩放后的括号，仍有非平凡的动力学结构。选取的初态是

$$
\begin{aligned}
|\boldsymbol n;N\rangle
&=\left[
\sqrt{\frac{1+z}{2}}\,|\uparrow\rangle
+e^{i\phi}\sqrt{\frac{1-z}{2}}\,|\downarrow\rangle
\right]^{\otimes N},\\
\boldsymbol n
&=\left(\sqrt{1-z^2}\cos\phi,\sqrt{1-z^2}\sin\phi,z\right),
\qquad -1\le z\le1,\\
\langle\widehat s_i\rangle_0&=n_i,\qquad
\frac12\langle\{\Delta\widehat s_i,\Delta\widehat s_j\}\rangle_0
=\frac{\delta_{ij}-n_in_j}{N}.
\end{aligned}
\tag{2}
$$

这是最大总自旋 j＝N/2 的对称部门，其维数为 N＋1；完整 N 比特空间维数仍为2^N。初态选择已保证集体涨落小，而非对所有量子态都如此。选取内部自主 Hamiltonian

$$
\begin{aligned}
H_N&=\frac{\chi}{N}J_z^2
=\frac{\chi}{4}I+\frac{\chi}{2N}
\sum_{a<b}\sigma_z^{(a)}\sigma_z^{(b)},\\
\hbar_{\rm eff}&=\frac2N,\qquad
\widehat h_N=\hbar_{\rm eff}H_N
=\frac{\chi}{2}\widehat s_z^2,\\
\|H_N\|&=\frac{|\chi|N}{4},\qquad
\widehat{\boldsymbol s}^{\,2}
=\left(1+\frac2N\right)I
\quad\text{（最大自旋部门）}.
\end{aligned}
\tag{3}
$$

1/N是明示的耦合缩放，χ及时间尺度也是输入。Hamiltonian 不随演化时刻改变，并保持总自旋部门、J_z及能量；没有依赖外部随时脉冲。但制备、识别同一自旋轴、选择耦合和读取仍是资源与建模输入。

式(3)在寄存器标签上是全连接两体耦合。它没有给出空间邻接、局域距离或传播光锥，不能将 N 当成已推导的物理体积。

## 3. 相空间是球面，而非平面或三维物理空间

对式(1)以 iℏ_eff 归一化，坐标括号与球面 Lie–Poisson 括号一致：

$$
\{s_i,s_j\}_{S^2}=\epsilon_{ijk}s_k,\qquad
\{f,g\}_{S^2}
=\boldsymbol s\cdot(\nabla f\times\nabla g),
\qquad
|\boldsymbol s|^2=1,\qquad
h(\boldsymbol s)=\frac{\chi}{2}s_z^2.
\tag{4}
$$

在球面上延拓 f、g 时，沿径向延拓的自由不改变这个切向括号。单位长度是 Casimir 叶条件，来自所选最大自旋及相干态分支，不是全部量子态的均值恒有单位长度。

在避开两极并选定方位角分支的局部区域，取 q＝φ、p＝z，则

$$
\{q,p\}=1,\qquad
\omega=d\phi\wedge dz,\qquad
\int_{S^2}\omega=4\pi.
\tag{5}
$$

这给出局部 Darboux 正则坐标，但没有一个覆盖全球面的普通实平面坐标对：φ有角度周期与极点奇性；若存在全局光滑实函数 q、p 使 ω＝dq∧dp，则 ω＝d(q\,dp) 在闭球面积分为零，与式(5)矛盾。

**三个 s 分量不是三个空间维度。** 它们受一个 Casimir 约束，所选经典相空间只有两维，即一个正则自由度。它不同于每个空间点上具有多个度规分量及共轭动量的 ADM 场相空间。

由式(4)得到经典流：

$$
\dot z=0,\qquad
\dot\phi=\chi z,\qquad
s_+^{\rm cl}(t)
=\sqrt{1-z^2}\,e^{i\phi}e^{i\chi zt},
\qquad s_+=s_x+is_y.
\tag{6}
$$

## 4. 精确量子均值：无需数值拟合

记单自旋升算子 σ_+＝(σ_x＋iσ_y)/2。对第a个自旋，翻转所产生的能量差只取决于其余自旋，所以 Heisenberg 演化精确为

$$
\begin{aligned}
e^{itH_N}\sigma_+^{(a)}e^{-itH_N}
&=\sigma_+^{(a)}
\exp\!\left(\frac{i\chi t}{N}\sum_{b\ne a}\sigma_z^{(b)}\right),\\
\left\langle e^{iu\sigma_z}\right\rangle_0
&=\cos u+iz\sin u,\qquad u=\frac{\chi t}{N},\\
\boxed{\;
\langle\widehat s_+(t)\rangle
=\sqrt{1-z^2}\,e^{i\phi}
[\cos u+iz\sin u]^{N-1}\;},\qquad
\langle\widehat s_z(t)\rangle=z.
\end{aligned}
\tag{7}
$$

第一式可直接在 σ_z 的共同本征态中按翻转能量差验证；第二行因初态为乘积而因子化。N−1来自其余自旋的个数，不能误换成N。N＝1时 Hamiltonian 仅为常数，式(7)确实不发生旋转。

这也澄清经典平均场方程在哪一步出现：有限 N 的精确方程含反对易乘积

$$
\frac{d\widehat s_+}{dt}
=\frac{i\chi}{2}
\left(\widehat s_z\widehat s_+
+\widehat s_+\widehat s_z\right).
\tag{8}
$$

把乘积期望近似成期望之积不是无条件代数恒等式；下面的界才说明本态族何时允许经典描述。

## 5. 固定时间的显式误差与集中性

设 X取±1且均值为z，则 cos u＋iz sin u＝E e^(iuX)。利用中心化变量、Taylor积分余项和 |E e^(iuX)|≤1，有

$$
\begin{aligned}
\left|\cos u+iz\sin u-e^{izu}\right|
&\le\frac{1-z^2}{2}u^2,\\
\left|a^{N-1}-b^{N-1}\right|
&\le(N-1)|a-b|\qquad (|a|,|b|\le1).
\end{aligned}
\tag{9}
$$

再计入 N−1 与经典相位 Nzu 的差，得对全部 |t|≤T 的界：

$$
\begin{aligned}
E_N(t)
&:=\left|\langle\widehat s_+(t)\rangle-s_+^{\rm cl}(t)\right|\\
&\le B_N(T,z):=
\sqrt{1-z^2}\left[
\frac{|z\chi|T}{N}
+\frac{(N-1)(1-z^2)\chi^2T^2}{2N^2}
\right].
\end{aligned}
\tag{10}
$$

这个证明不要求先作小角截断；对任意有限 N、T 都是上界。只有在说误差趋零时才取固定T、N趋无穷。三分量均值的欧氏误差与式(10)相同，因为z分量没有误差。

仅平均值接近还不够。由 Casimir 和式(7)，总集体方差精确为

$$
\begin{aligned}
V_N(t)
&:=\sum_i\operatorname{Var}(\widehat s_i)
=\frac2N+(1-z^2)
\left[1-\left(1-(1-z^2)\sin^2u\right)^{N-1}\right]\\
&\le\frac{2+(1-z^2)^2\chi^2T^2}{N},\\
\sum_i\left\langle
\left(\widehat s_i-s_i^{\rm cl}(t)I\right)^2
\right\rangle
&=V_N(t)+
\left|\langle\widehat{\boldsymbol s}(t)\rangle
-\boldsymbol s^{\rm cl}(t)\right|^2
\le\frac{2+(1-z^2)^2\chi^2T^2}{N}+B_N(T,z)^2.
\end{aligned}
\tag{11}
$$

对固定集体分量的单次测量，可以进一步用其方差与 Chebyshev 不等式控制偏离概率。这没有假设三个不对易量可无扰动联合精确测量，也没有给所有微观／全局可观察量一个统一误差保证。

因此经典轨道同时得到均值与集体集中性的支持。这里不需要先让单个比特彻底退相干；完整总态始终纯态酉演化，内部相关性会增长。本轮不依赖任何单比特退相干模型。

## 6. 时间随资源增长后，哪个结论失效？

假设χ≠0、|z|＜1，并取 t_N＝τ√N/|χ|。由 log(cos u＋iz sin u)的局部展开，得到去除经典旋转后的极限：

$$
\begin{aligned}
\frac{\langle\widehat s_+(t_N)\rangle}
{s_+^{\rm cl}(t_N)}
&\longrightarrow
\exp\!\left[-\frac{1-z^2}{2}\tau^2\right],\\
E_N(t_N)
&\longrightarrow
\sqrt{1-z^2}
\left(1-e^{-(1-z^2)\tau^2/2}\right),\\
V_N(t_N)
&\longrightarrow
(1-z^2)\left(1-e^{-(1-z^2)\tau^2}\right).
\end{aligned}
\tag{12}
$$

这是经典点轨道在增长时间上的有限偏离，不是数值不稳定。赤道 z＝0、τ＝1 时，平均横向对比度趋向0.6065307，均值误差趋向0.3934693，总方差趋向0.6321206。两极没有横向自旋，不能套用除以 s_+ 的式子；它们也不是这个失效见证。

式(10)实际上保证 T_N＝o(√N) 时该均值误差趋零，式(11)也给集中性；本轮并不宣称 √N 是所有初态、Hamiltonian 和全部观测量的普适经典时间界。

有限 N 模型还是可逆的。t＝2πN/|χ|使全部 J_z本征分量恢复同一整体相位，故有完整态回归；“均值塌缩”不能解释成已证明不可逆熵产生。

## 7. 更强的范围检查：固定时间也未必有全态收敛

### 7.1 同一相干初态，宏观经典并不等于完整态仍是乘积

取赤道初态 |＋⟩^⊗N。其经典轨道在φ不动；定义**平方保真度**为精确量子态对这个经典相干乘积态的重叠平方。初态中的 J_z是N个独立±1/2变量之和，因此

$$
\begin{aligned}
F_N(t)
&=\left|
\langle+|^{\otimes N}e^{-i\chi tJ_z^2/N}|+\rangle^{\otimes N}
\right|^2,\\
\frac{J_z}{\sqrt N}
&\ \Longrightarrow\ {\cal N}(0,1/4),\\
F_N(t)&\longrightarrow
\left|\mathbb E_{Y\sim{\cal N}(0,1/4)}
e^{-i\chi tY^2}\right|^2
=\frac1{\sqrt{1+(\chi t/2)^2}}.
\end{aligned}
\tag{13}
$$

此处中央极限定理只用于初态的经典二项分布；指数是有界连续函数，故可直接传递期望极限，再计算 Gaussian 积分。不是把整个量子态先假设成经典概率分布。

χt＝2时，F极限为1/√2，两个纯态的迹距趋向约0.541196；与此同时，式(10)—(11)中的宏观误差仍趋零。内部微观相关和放大的涨落可保留量子差别。

这反驳的是“经典点轨道可直接提升为完整相干乘积态近似”，并非反驳经典有效理论。一般量子到经典通道极限还须定义观察代数、嵌入与比较范数；本轮没有完成或宣称这种更强极限。

### 7.2 任意输入态也不自动集中

同一对称部门、同一 Hamiltonian 还允许

$$
|{\rm GHZ}_N\rangle
=\frac{|\uparrow\rangle^{\otimes N}
+|\downarrow\rangle^{\otimes N}}{\sqrt2},
\qquad N\ge2:\quad
\langle\widehat{\boldsymbol s}\rangle=0,\quad
\operatorname{Var}(\widehat s_z)=1.
\tag{14}
$$

两端分量具有相同的 J_z²，故这个态仅变化整体相位；方差不随N消失。它不接近单一球面点。若只看固定阶集体矩在N趋于无穷时的极限，它可以对应两极的经典混合；本例不否定经典概率分布型极限。

但完整量子 GHZ 与去掉交叉项的混合态仍有迹距1/2，而且全局 X^⊗N 的期望分别为1与0。不能用“集体量可以有经典描述”删除全部量子相干，或说已经对所有未知输入完成经典化。

## 8. 可复算实现与结果

代码用两条独立有限矩阵实现核对式(7)：

1. N＋1维自旋j＝N/2矩阵，按 J_z本征基作精确相位演化；另算均值、方差和保真度。
2. 小N时构造完整2^N维 Pauli 张量积 Hamiltonian 与同向乘积初态，直接核验同一公式。前者的对称性压缩没有隐藏物理系统所需的N个寄存器。

固定χ＝0.8、z＝0.4、φ＝0.3、T＝1.2：

| N | 实际均值误差 | 式(10)上界 | 总集体方差 | 式(11)方差上界 |
|---:|---:|---:|---:|---:|
| 16 | 0.0296888 | 0.0427829 | 0.162263 | 0.165643 |
| 64 | 0.00770843 | 0.0109556 | 0.0411928 | 0.0414106 |
| 256 | 0.00194574 | 0.00275513 | 0.0103389 | 0.0103527 |
| 1024 | 0.000487611 | 0.000689798 | 0.00258731 | 0.00258817 |

赤道增长时间 t＝√N、χ＝1，N＝64到4096时均值误差由0.389499变为0.393408，趋向式(12)的非零常数，不能被固定时间表格掩盖。

另一组固定t＝2、χ＝1：N＝2048时宏观均值误差约0.000975610，但全态平方保真度约0.70710678118，迹距约0.54119610015。两种结论可在同一个模型、同一个初态上同时成立。

14项检查覆盖：精确交换关系／有效ℏ；Casimir及相干协方差；对称矩阵精确矩；完整张量实现；固定时间界；方差与集中性；能量及z分布守恒；√N失效；态回归；球面Poisson／Darboux局部坐标；全态保真度边界；GHZ及全局测量；耦合资源恒等式；有限N的精确Heisenberg方程。解析结论不依赖有限N采样的外推。

## 9. 对 GR 主线究竟增加了什么？

本轮给出一个有限量子操作结构内部的**条件性正例**：指定集体变量、相干态族及耦合后，可以定量得到经典 Poisson 动力学，并控制其时间窗口。经典变量不必作为与量子世界无关的另一种基本材料直接加入。

但351轮需要的对象并没有由此出现：

- 球面上的一个正则坐标对不是空间点上的 q_ij 与 π^ij；本轮没有生成空间指标、图、平滑空间体积或正定度规场。
- Poisson 括号不等于超曲面形变代数。这里的结构系数是自旋 s；ADM 约束括号中的 q^ij、lapse、shift 与局域导数尚未构造。
- 阶段二允许实现这些有限量子状态和酉变换，不等于认知原则自然选择了同向乘积态、最大自旋部门、J_z²或1/N缩放。
- N增加同时增加实际寄存器和耦合资源；对称矩阵较小仅是可复算优势，不是无成本生成连续自由度。

下一实质接口应把候选集体变量与关系／几何可观察量联系起来，证明其有效正则代数和约束，并检查额外内部模式在所需时空／能标下如何被控制。仅把自旋球面重新命名为“空间”或“度规”，不能补上这个缺口。

本轮没有因G、Λ或初边值的数值尚未计算而宣告 GR 方程族不完整；开放项是量子集体结构与有效几何／约束之间的生成关系。

[代码](354/collective_poisson_limit_audit.py)、[结果](354/collective_poisson_limit_audit_results.json)、[核验](354/research_round_354_checks.json)。

    python -B -X utf8 research_cognition_physics/archive_231_/collective_poisson_limit_audit.py

加 --write-results 可在检查通过后保存；已有不同结果拒绝覆写。本轮未修改旧轮、导航或其他并行文件。
