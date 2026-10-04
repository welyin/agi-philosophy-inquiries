# 第358轮：已给定无质量自旋2后，非线性自耦合何时成为 Einstein 理论

## 0. 本轮问题、基线与结论

本轮以截至356轮的完成结果为基线，独立于同批尚未冻结的357、359轮。已读取根 README、本研究三个导航文件，检索既有自旋2／自耦合研究，并重读[326轮软自旋2接口](../archive_301_341/research_note_326.md)、[344轮引力唯一性审计](research_note_344.md)、[351轮超曲面形变约束](research_note_351.md)。此前326轮处理软极限下的共同耦合，没有完成本轮的一阶非线性自耦合；351轮则从另一组已给定的正则几何条件出发。

**本轮补齐一条有条件的正向接口：给定健康的单一无质量 Pauli–Fierz 自旋2自由场，并限定局域、Poincaré 不变、最多两导数及逐阶规范一致的非平凡相互作用，已有分类定理将其限定到 Einstein–Hilbert 类；Deser 一阶构造可以明确给出该类中的完整非线性作用量。**

本轮独立复算构造中的关键环节：自由场的辅助背景变分、完整应力源、非线性联络消元、Palatini 作用量恒等式和变分。没有仅把已知 Einstein 作用量展开，然后宣称已从认知原则导出它。

**从认知到上述自由自旋2、局域 Lorentz 场论及导数限制的入口仍未证明。** 本轮没有使这一总缺口消失，也没有否定 GR 作为低能有效理论的可能性。

| 层次 | 本轮内容 |
|---|---|
| 认知动机 | 相互影响的系统应把自身携带的能量也纳入共同动力学 |
| 额外输入 | 四维 Minkowski 背景、单一 Pauli–Fierz 场、局域作用量、导数阶、非平凡耦合及规范一致性 |
| 既有定理 | 受限自旋2一致形变分类；不是本项目原创 |
| 本轮解析复算 | 一阶三次作用量、完整源的产生、联络消元和边界等价 |
| 本轮数值验证 | 15项独立或交叉检查；包含实际作用量方向差分 |
| 物理解释 | 已有自旋2后的非线性完成；不是认知原则到 GR 的最终证明 |

## 1. 对接原始研究，并区分“构造”与“唯一性”

[Deser（1970），Self-Interaction and Gauge Invariance](https://arxiv.org/abs/gr-qc/0411023)的核心技术是使用独立联络的一阶变量，使相互作用在作用量中只需增加一个三次项。自由场的完整应力源还含联络消元产生的导数项，不能把它简化成“往作用量加上 h 乘应力”。

[Deser（2009预印本、2010发表），Gravity from self-interaction redux，v3全文](https://arxiv.org/pdf/0910.2975)明确承认导数阶限制是输入，并解释应力、边界项和耦合常数重标度。其 PDF 排版日期不是论文提出该结果的年份。本轮只对下面明确构造的局域场重定义作结论，不将文中“任意守恒改进均可移除”的泛称当作已证明的局域定理。

更精确的唯一性范围来自[Boulanger、Damour、Gualtieri、Henneaux（2001），原文定理1.1及第2—4页](https://arxiv.org/pdf/hep-th/0007220)：在局域、Poincaré 不变、Pauli–Fierz 自由极限、Lagrangian 最多两导数及可逐阶展开的规范一致形变条件下，允许的理论是独立 Einstein–Hilbert 部门或仍自由的 Pauli–Fierz 部门，模去相应场重定义；宇宙项亦在分类内。该文也给出放宽导数限制后的曲率三次相互作用。

因此，要说“得到相互作用的 Einstein 理论”，还必须排除完全自由的分支。**允许自旋2自由场本身，不会迫使非零引力相互作用。** 本轮显式构造检验的是指定自耦合接口的充分性；受限类别中的唯一性使用上述既有分类结果，不能由一次数值拟合或一次成功构造替代。

## 2. 自由模型、记号与假设账目

采用四维背景及紧支撑变分，或足以使本轮分部积分成立的边界条件。所有式中的时空维数、微分结构和 Lorentz 符号均已给定：

$$
\eta_{\mu\nu}=\operatorname{diag}(-1,1,1,1),\qquad
h^{\mu\nu}=h^{\nu\mu},\qquad
C^\alpha{}_{\mu\nu}=C^\alpha{}_{\nu\mu}.
\tag{1}
$$

h 是相对于背景的逆变二阶张量密度，权重为1；C 在自由背景变分中作为普通张量。它们是独立变量。定义对称 Ricci 的线性和二次部分：

$$
\begin{aligned}
P_{\mu\nu}(C)
&=\partial_\alpha C^\alpha{}_{\mu\nu}
-\frac12\bigl(\partial_\mu C^\alpha{}_{\alpha\nu}
+\partial_\nu C^\alpha{}_{\alpha\mu}\bigr),\\
Q_{\mu\nu}(C)
&=C^\alpha{}_{\mu\nu}C^\beta{}_{\alpha\beta}
-C^\alpha{}_{\mu\beta}C^\beta{}_{\alpha\nu},\\
\mathcal L_2&=h^{\mu\nu}P_{\mu\nu}(C)
+\eta^{\mu\nu}Q_{\mu\nu}(C).
\end{aligned}
\tag{2}
$$

引入度规密度的联络作用算子：

$$
\begin{aligned}
(\mathcal M_G C)_a{}^{mn}
&=C^m{}_{ar}G^{rn}
+C^n{}_{ar}G^{mr}
-C^r{}_{ar}G^{mn},\\
C_L(h)&=-\mathcal M_\eta^{-1}(\partial h),\\
B(h,C)&=-\mathcal M_\eta^{-1}(\mathcal M_h C).
\end{aligned}
\tag{3}
$$

这里逆的是每一点有限维的代数算子，不是微分方程的非局域 Green 函数。四维、对称下指标的 C 有40个分量。为明确其局域性，令 F＝M_η C，将 F 后两个指标及 C 首指标用 η 降下，可直接反解：

$$
\begin{aligned}
T_a&=-\frac12\eta^{mn}F_{a mn},\\
C_{a mn}
&=\frac12\bigl[
F_{m an}+F_{n am}-F_{a mn}
+T_m\eta_{an}+T_n\eta_{am}-T_a\eta_{mn}
\bigr].
\end{aligned}
\tag{4}
$$

对式(2)作联络变分即得 C＝C_L(h)。把密度扰动转为协变度规扰动 k 后，得到通常的线性联络：

$$
\begin{aligned}
k_{\mu\nu}(h)
&=-\eta_{\mu\alpha}\eta_{\nu\beta}h^{\alpha\beta}
+\frac12\eta_{\mu\nu}\eta_{\alpha\beta}h^{\alpha\beta},\\
C_L^\alpha{}_{\mu\nu}
&=\frac12\eta^{\alpha\rho}
\bigl(\partial_\mu k_{\rho\nu}
+\partial_\nu k_{\rho\mu}-\partial_\rho k_{\mu\nu}\bigr),\\
P_{\mu\nu}(C_L(h))&=R^L_{\mu\nu}(k).
\end{aligned}
\tag{5}
$$

因而自由 h 方程就是线性真空 Ricci 方程，等价于四维线性 Einstein 方程。通常线性规范变换给出零线性曲率：

$$
\delta_\xi k_{\mu\nu}
=\partial_\mu\xi_\nu+\partial_\nu\xi_\mu,
\qquad
R^L_{\mu\nu}(\delta_\xi k)=0.
\tag{6}
$$

这里是对给定 Pauli–Fierz 自由理论的等价写法，不是从有限维量子操作权限推出空间张量、自旋或上述 Hamiltonian。

## 3. 独立计算完整应力源：为什么只取 Q 不够

用辅助逆变度规密度 γ 替代背景 η，K 为其 Levi-Civita 联络。把式(2)的导数协变化，h 仍为权重1密度、C 为张量。辅助作用量取

$$
\begin{aligned}
I_2[\gamma;h,C]
&=\int d^4x\,
\left[h^{mn}\bigl(P_{mn}(C)+A_{mn}(K,C)\bigr)
+\gamma^{mn}Q_{mn}(C)\right],\\
A_{mn}(K,C)
&=K^a{}_{ar}C^r{}_{mn}
-K^r{}_{am}C^a{}_{rn}
-K^r{}_{an}C^a{}_{mr}
+K^r{}_{mn}C^a{}_{ar}.
\end{aligned}
\tag{7}
$$

γ 在此只是计算源的辅助变量，不新增物理背景自由度。K 对 γ 在 η 处的一阶变化为 C_L(δγ)。把式(4)—(5)代入式(7)，逐项收缩并对 δγ 的导数分部积分，得到

$$
\begin{aligned}
\int d^4x\,
h^{mn}A_{mn}\bigl(C_L(\delta\gamma),C\bigr)
&=\int d^4x\,
\delta\gamma^{mn}P_{mn}\bigl(B(h,C)\bigr),\\
\left.\delta I_2\right|_{\gamma=\eta}
&=\int d^4x\,\delta\gamma^{mn}\tau_{mn},\\
\tau_{mn}
&=Q_{mn}(C)+P_{mn}\bigl(B(h,C)\bigr).
\end{aligned}
\tag{8}
$$

这些恒等式对任意光滑 h、C 成立，不要求先解真空方程。式(7)列出全部四个协变化项，式(4)给出代数逆，因而式(8)可以直接逐指标复算。τ 是本作用量归一化下、以逆变密度变分定义的迹调整源；转换成通常度规应力时还需相应迹反转及归一化。不能把 τ 自身不加区别地当作通常 Einstein 方程右边的 T。

数值实现的两侧独立：左侧从辅助度规、Christoffel 符号和实际 I_2 方向差分计算；右侧才使用 M_η 的代数逆和 B。两者不是调用同一个“应力公式”。

现在加入三次项，而非先要求完整作用量等于 Einstein：

$$
\mathcal L=\mathcal L_2+\kappa h^{mn}Q_{mn}(C),
\qquad \kappa\ne0.
\tag{9}
$$

这里 κ 是场重标度采用的耦合参数，与通常 Einstein 方程中乘 T 的常数不能直接同名等同。新项中没有 η、辅助 γ 或导数；在这套协变化选择中，它没有新的辅助背景应力项。由式(9)直接变分得到

$$
\begin{aligned}
P(C)+\kappa Q(C)&=0,\\
\partial h+\mathcal M_\eta C+\kappa\mathcal M_h C&=0,\\
C&=C_L(h)+\kappa B(h,C).
\end{aligned}
\tag{10}
$$

第二行是经取迹整理后的联络方程。对第三行施加线性算子 P，再用第一行消去 P(C)，即得所需完整源：

$$
R^L_{mn}(k)
=-\kappa\left[Q_{mn}(C)+P_{mn}\bigl(B(h,C)\bigr)\right]
=-\kappa\tau_{mn}.
\tag{11}
$$

**关键在这里：三次作用量只显含 hQ，方程中的完整源却是 Q＋P(B)。** 后者的导数部分来自联络方程的非线性修正。若直接把 Q 当作全部应力，就会漏掉可直接计算的作用量变分。

“一轮完成”指独立变量下的作用量仅需该三次项，不意味着消去 C 后度规作用量只有三次，也不意味着解的非线性演化只需一次迭代。

## 4. 识别完整 Palatini 作用量及边界范围

定义新变量：

$$
G^{mn}=\eta^{mn}+\kappa h^{mn},
\qquad
\Gamma^\alpha{}_{\mu\nu}=\kappa C^\alpha{}_{\mu\nu},
\qquad
R_{(mn)}(\Gamma)=P_{mn}(\Gamma)+Q_{mn}(\Gamma).
\tag{12}
$$

直接展开，而无需场方程，得到逐点恒等式：

$$
\begin{aligned}
\kappa^{-2}G^{mn}R_{(mn)}(\Gamma)
&=\mathcal L_2+\kappa h^{mn}Q_{mn}(C)
+\kappa^{-1}\partial_\alpha V^\alpha,\\
V^\alpha
&=\eta^{mn}C^\alpha{}_{mn}
-\eta^{\alpha n}C^\beta{}_{\beta n}.
\end{aligned}
\tag{13}
$$

因此式(9)的体内 Euler–Lagrange 方程就是 Palatini 形式。表面的 κ 的负幂并不要求自由理论以负幂展开：式(13)中该部分是总散度，剩余作用量在 κ 中是显式有限多项式。

必须固定**非退化、一个时间方向、与 η 连通的分支**。四维中可定义

$$
g_{\mu\nu}
=\sqrt{-\det G}\,(G^{-1})_{\mu\nu},
\qquad
G^{\mu\nu}=\sqrt{-g}\,g^{\mu\nu}.
\tag{14}
$$

仅检验 det G＜0 不足以排除三个时间方向，代码另检验负本征值的个数。对称联络变分给出

$$
\begin{aligned}
-\nabla_\lambda G^{\mu\nu}
+\delta_\lambda^{(\mu}\nabla_\rho G^{\nu)\rho}&=0
\quad\Longrightarrow\quad \nabla_\lambda G^{\mu\nu}=0,\\
\Gamma&=\Gamma_{\mathrm{LC}}(g),\\
R_{\mu\nu}(g)=0
&\quad\Longleftrightarrow\quad
R_{\mu\nu}(g)-\frac12 g_{\mu\nu}R(g)=0.
\end{aligned}
\tag{15}
$$

于是局域真空体内作用量就是 Einstein–Hilbert，外加已明确的边界项。这里采用无挠独立联络；若允许挠率、物质自旋直接耦合联络或度规退化，需要重新处理方程和分支。

**边界等价不是“边界物理无关”。** 式(13)在紧支撑变分／适当无边界设定下保证体内方程一致；有限边界的良定变分、边界应力、守恒荷与作用量值，须连同边界条件及所需边界项一起匹配。本轮不由丢掉总散度推出这些对象自动相同。

## 5. 受限唯一性、局域改进与不唯一项

### 5.1 可明确证明的局域改进

若改进项能写成线性 Einstein 算子作用于某个局域张量，它可以在相应微扰阶通过局域场移位处理。一个显式例子为

$$
\begin{aligned}
\Delta_{\mu\nu}(F)
&=(\partial_\mu\partial_\nu-\eta_{\mu\nu}\Box)F,\\
\partial^\mu\Delta_{\mu\nu}(F)&=0,\\
G^L_{\mu\nu}(-\eta F)&=\Delta_{\mu\nu}(F).
\end{aligned}
\tag{16}
$$

若 F 是字段的局域二阶小量，这给出保留自由极限的近恒等微扰场重定义。代码验证任意测试 F 的算子恒等式；完整场重定义还需满足相应微扰可逆性。

**本轮没有证明“散度恒为零”足以保证任意改进都有局域可逆的原像。** 求一个微分算子的逆通常可能引入非局域性。不能用形式上解出一个场移位，就把任意高导数项、额外传播自由度或不同边界问题全部消去。

### 5.2 两导数限制确实有内容

考虑非线性度规作用量中的高导数候选：

$$
\begin{aligned}
\Delta S
&=\beta\int d^4x\,\sqrt{-g}\,
R_{ab}{}^{cd}R_{cd}{}^{ef}R_{ef}{}^{ab},\\
g_{\mu\nu}&=\eta_{\mu\nu}+\varepsilon k_{\mu\nu}
\quad\Longrightarrow\quad
\Delta S=O(\varepsilon^3).
\end{aligned}
\tag{17}
$$

它从平直背景的三阶开始，因而不改二次 Pauli–Fierz 自由作用量；它含六导数，已经超出本轮引用的两导数分类。也可直接用线性曲率构造保持线性规范不变的六导数顶点，原始分类文献对此给出明确例子。

为排除“这只是代数恒等为零或全是 Ricci 项”的误解，构造一个代数曲率张量：

$$
\begin{aligned}
E&=\operatorname{diag}(-2,1,1),\\
\mathcal R_{0i0j}&=E_{ij},\qquad
\mathcal R_{ijkl}=-\epsilon_{ijm}\epsilon_{kln}E_{mn},
\qquad \mathcal R_{0ijk}=0,\\
\eta^{ac}\mathcal R_{abcd}&=0,\qquad
\mathcal R_{ab}{}^{cd}\mathcal R_{cd}{}^{ef}
\mathcal R_{ef}{}^{ab}=96.
\end{aligned}
\tag{18}
$$

它满足 Riemann 的指标对称性和代数 Bianchi 恒等式。此检查证明曲率三次收缩并非代数零式，也不是纯 Ricci 多项式；**单点代数检查不等于完整时空解、边界非平凡性分类或健康的高能理论证明。** 本轮仅用它展示导数限制不能从“已有同一自由自旋2”中省略。作为低能有效项时，其影响还依赖系数、截止尺度与态，并不排除 GR 主导红外行为。

344轮 R＋αR² 在平直附近带额外标量，所以不属于“只有单一 Pauli–Fierz 自由场”的本轮类别；两轮没有矛盾。本轮曲率三次例子则进一步说明：即使固定同一二次自由作用量，若放开导数限制仍不能仅凭自由谱宣称精确唯一性。

### 5.3 物质、宇宙项和数值常数

- 本轮完整解析构造是纯引力真空部门。加入物质需为整个物质作用量施加与最终规范对称性相容的耦合条件；由 Noether 恒等式得到的协变守恒，并不唯一决定物质谱、势、非最小曲率耦合或标准模型。
- 326轮软极限给出的是其假设下的普适软耦合约束，不能自动替代完整硬动量相互作用及局域物质作用量分类。
- 本轮显式构造采用无宇宙项、平直真空的分支；这不证明宇宙学常数必须为零。已有分类允许宇宙项，Deser 路线也可从常曲率背景展开。
- G 与 Λ 的数值不在本轮输出中；它们未被计算并不单独构成本项目接受有效 GR 的否决条件。真正未证的桥梁是自由自旋2及其所需结构、相互作用、尺度分离与导数截断的形成。

## 6. 可复算结果与独立性

代码：[spin2_self_coupling_audit.py](358/spin2_self_coupling_audit.py)；结果：[spin2_self_coupling_audit_results.json](358/spin2_self_coupling_audit_results.json)；轮次交付检查：[research_round_358_checks.json](358/research_round_358_checks.json)。

采用确定随机种子358，96个周期相位采样点。场依赖相位 u＝ℓ_μ x^μ，ℓ＝(0.4, 1, −0.35, 0.2)，因此四个微分方向均参与完整四维指标运算；这不是宇宙周期时间模型。它是一般恒等式的可复算样本，不能承担全部场配置的量词。一般性由上面的代数和变分推导承担。

| 检验对象 | 结果 |
|---|---:|
| M_η 的秩／条件数 | 40／2.61803 |
| 三次作用量与 Palatini 逐点恒等式残差 | 3.89×10⁻¹⁶ |
| 总散度密度最大值／周期积分 | 0.21020／−5.28×10⁻¹⁸ |
| 完整辅助应力方向变分 | 0.001341810417426 |
| 只用 Q 的方向变分 | 0.001106564063181 |
| 漏掉的 P(B) 贡献 | 0.000235246354245 |
| 辅助作用量差分 Richardson 误差 | 8.25×10⁻¹⁶ |
| 非线性度规相容性残差 | 1.67×10⁻¹⁴ |
| C＝C_L＋κB 残差 | 1.13×10⁻¹⁴ |
| 联络非线性修正最大值 | 0.01552 |
| 真 Levi-Civita 联络处作用量方向变分 | −1.30×10⁻¹⁴ |
| 人为偏离联络处同一方向变分 | −0.06659 |
| 消元后的完整 Ricci／源恒等式残差 | 7.28×10⁻¹³ |
| 代数 Ricci 残差／曲率三次收缩 | 0／96 |

辅助背景中央差分的步长依次为0.01、0.005、0.0025；误差依次约7.07×10⁻¹¹、1.77×10⁻¹¹、4.41×10⁻¹²，符合二阶收敛。实际用于联络消元的度规具有非零 Ricci，最大分量约0.21370：因此源恒等式验证没有靠选择平直或真空解令两侧同时变零。

15项检查还覆盖：密度—度规往返及错误符号分支拒绝、代数逆与直接线性 Christoffel 比较、纯规范零曲率、h 的作用量方向变分、明确局域改进、曲率代数对称性与三次缩放。数值部分没有重新求解完整量子引力或重做旧轮测试。

复算方式，在本目录运行：

~~~powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 -m unittest spin2_self_coupling_audit.Checks -v
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 spin2_self_coupling_audit.py --write-results
~~~

先独立执行 Checks 全部通过，再经既有 growing_stream_audit.main 再次检查并写入结果。运行时为 Python 3.12.14、NumPy 2.3.5，无新增依赖。

## 7. 对总目标的实际推进与下一问题

本轮把326轮之后的一个问题从“自旋2大概会给出引力”推进为**有明确输入、已知分类依据及完整一阶构造的非线性接口**。一旦某个认知—操作模型真正产生符合条件的低能单一无质量自旋2部门，就无需重新猜测每一阶 Einstein 相互作用。

但不能倒过来把这条条件定理解释为自由自旋2已经由认知原则得到。与351轮对照，两条路线的后半段都已有可靠工具：一条给定超曲面形变与纯度规变量，一条给定自由自旋2及规范一致性；两者共同指向 Einstein 类，同时共同暴露入口的额外物理结构。

**下一项有实质意义的检验**应定位在入口：给定可复算的集体模型，是否出现受保护的无质量横向无迹张量自由度、相应局域冗余及足够的尺度分离；为什么其他低能模不会破坏单一自旋2有效描述。若还需要人工指定 Pauli–Fierz 动力学，就应明确登记该输入，而不是继续复述已知自耦合定理来增加轮次。
