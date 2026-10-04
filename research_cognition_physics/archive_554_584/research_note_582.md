# 第582轮：同一规范质量、量子连接与共同应力

日期：2026-10-01。接[581条件压缩](research_note_581.md)、[580标量有效项](research_note_580.md)及[572共同规范来源](research_note_572.md)。链接按archive_231_解析；完成以独立终审和冻结为准。

## 1. 为什么必须把规范背景接回

581在零规范曲率的标量部门压缩了局部有效项，但572的同一物质来源含非零规范场。本轮保留原Higgs＋singlet目标、原SU(2)×U(1)荷及同尺度归一，把非零规范连接和应力接回同一标量环系数。

**有两条可合并的新联系：** 规范质量的目标Gram矩阵同时固定一部分场依赖规范热核系数；使用完整Einstein—标量—规范领先作用压缩时，又必须保留确定的物质—规范应力交叉项。两者都不能从“只保留规范能量大小”得到。

这是一项指定模型、指定一圈部门的共同实现与条件限制。未包括规范场、ghost、费米子及引力内部环，也未构造量子修正后的572完整在壳来源；物理维数、群、作用和尺度依然有输入。

## 2. 成熟映射与归一

[Alonso–Jenkins–Manohar，1511.00724，式(73)、(79)](https://arxiv.org/html/1511.00724)给带规范背景的sigma模型协变Hessian和连接曲率。580已经映射其零规范部门，本轮只补同一模型的非零规范曲率，不重述一般热核作为新理论。

保留四维Euclidean背景、无边界局部密度及580的极点约定。五实分量按φ=(ReX0,ReX1,ImX0,ImX1,s)排列，F=M−φ²/6、M=2；目标K、U来自574。弱连接已吸收物理耦合，圆连接a0对应Q=6Y：

$$
D_\mu X=\partial_\mu X+i a_\mu^a\frac{\sigma_a}{2}X+3i a_{0\mu}X,
\qquad f^a_{\mu\nu}=\partial_\mu a^a_\nu-\partial_\nu a^a_\mu-\epsilon^{abc}a^b_\mu a^c_\nu,
\qquad f^0_{\mu\nu}=\partial_\mu a_{0\nu}-\partial_\nu a_{0\mu}.
\tag{1}
$$

圆荷3由实生成元承担，不能再给f0多乘一次g0。规范动能权w_A=(1/g_w²,1/g_w²,1/g_w²,36/g_Y²)直接取569保存的同尺度参数；规范动能是Σ_A w_A f_A²/4。Higgs不带颜色，故本轮标量环的色生成元为零；色场若作为额外经典来源，仍进入总应力。

|层次|本轮地位|
|---|---|
|认知动机|共同物质、规范作用与几何来源不能分别更换|
|继承输入|原目标、势、弱与整数圆表示、同尺度规范系数、四维背景|
|量子范围|只积分五标量，冻结规范及度规背景；本轮不作完整量子规范固定／ghost计算|
|压缩范围|以动态Einstein—标量—规范领先作用定义一阶EFT重定义，未加入引力环|
|解析增量|质量Gram与连接平方关系、真实局部规范资料、规范应力修正后的完整离壳分解|
|未完成|完整环总和、独立可观测算符基、连续图匹配、来源制备及在壳量子几何|

## 3. 同一Killing场连接质量与量子规范系数

令T_A是iσ_a/2及3iI在上述四实Higgs坐标上的矩阵，在singlet方向为零。它们反对称，k_A=T_Aφ。Klein坐标的Christoffel及Killing导数为

$$
\Gamma^a{}_{bc}=\frac{\delta^a_b\phi_c+\delta^a_c\phi_b}{6F},\qquad
J_A{}^a{}_b:=\nabla_b k_A^a=T_A{}^a{}_b+\frac{k_A^a\phi_b}{6F},\qquad
KJ_A+J_A^\top K=0.
\tag{2}
$$

由φ·k_A=0，角向K作用为Kk_A=k_A/F。定义原连接变量下的质量Gram和常数表示指标：

$$
\mathcal M_{AB}=K(k_A,k_B)=\frac{k_A\cdot k_B}{F},\qquad
\mathcal I_{AB}=-\operatorname{tr}(T_AT_B)=\operatorname{diag}(1,1,1,36)_{AB},\qquad
\boxed{-\operatorname{tr}(J_AJ_B)=\mathcal I_{AB}+\frac13\mathcal M_{AB}.}
\tag{3}
$$

证明：展开J_AJ_B的两个交叉迹，各为−k_A·k_B/(6F)；二次外积迹因φ·k=0而消失。因此同一目标的规范质量Gram固定场依赖的连接平方。它不是另设一份质量或量子耦合。

物理规范场须先标准化。记W_g=diag(w_A)，则

$$
\mathcal M_{\rm can}=W_g^{-1/2}\mathcal M W_g^{-1/2},\qquad
X=(0,h_*)\ \Longrightarrow\ 
\operatorname{spec}\mathcal M_{\rm can}
=\left\{0,\frac{g_w^2h_*^2}{4F_*},\frac{g_w^2h_*^2}{4F_*},\frac{(g_w^2+g_Y^2)h_*^2}{4F_*}\right\}.
\tag{4}
$$

这是当前Einstein目标的常场二次质量矩阵；与569平直目标相比有明确1/F因子，不能把两个不同目标的数值直接等同。它尚非完整量子极点质量、低能实验拟合或新粒子预测。

## 4. 规范曲率补入同一标量行列式

在目标正交标架记X_mu=D_muφ，M_X=ΣX_mu X_muᵀ，A=∇²_K U。规范曲率只在连接交换子中新增项；在所列规范动能系数不依赖φ且冻结规范背景的处方下，P没有额外直接f项：

$$
P=A+\frac16(SI-M_X),\qquad
\Omega_{\mu\nu}=\Omega^{(0)}_{\mu\nu}+\mathcal G_{\mu\nu},\quad
\Omega^{(0)}_{\mu\nu}=-\frac16(X_\mu X_\nu^\top-X_\nu X_\mu^\top),\quad
\mathcal G_{\mu\nu}=\sum_A f^A_{\mu\nu}J_A.
\tag{5}
$$

J已在同一正交标架表示。由580的tr Ω²/12，额外局部系数是

$$
\Delta_g b_4=
\frac16\sum_{\mu,\nu}\operatorname{tr}(\Omega^{(0)}_{\mu\nu}\mathcal G_{\mu\nu})
+\frac1{12}\sum_{\mu,\nu}\operatorname{tr}\mathcal G_{\mu\nu}^2
=\Delta_{\rm mix}
-\frac1{12}\sum_{\mu,\nu,A,B}f^A_{\mu\nu}f^B_{\mu\nu}
\left(\mathcal I_{AB}+\frac{\mathcal M_{AB}}3\right).
\tag{6}
$$

反对称矩阵平方迹一般非正，不可替换成正Frobenius平方。μν双求和已计两个方向。式(6)是局部热核系数而非有限物理耦合；完整量子群和全部物质环可以改变净系数。

原群变换同时作用φ、Dφ及曲率，Ω与G在切空间共轭，式(6)的两项分别规范不变。本轮数值使用原SU(2)变换验证，没有用不属于规范群的一般目标旋转代替。

## 5. 同规范应力不能决定这个局部系数

在Abelian子部门固定φ、X，令f→−f。纯规范平方与经典应力不变，混合项变号：

$$
\tau_{\mu\nu}^{(0)}=w_0\left(f_{\mu\rho}f_\nu{}^\rho-\frac14g_{\mu\nu}f_{\rho\sigma}f^{\rho\sigma}\right),\quad
\tau[f]=\tau[-f],\qquad
b_4[f]-b_4[-f]=2\Delta_{\rm mix}.
\tag{7}
$$

保存的局部例有Δmix=.0181983387223，差=.0363966774446。它们不是同一规范轨道，因为规范不变系数不同。结论只限制“同φ、同X及同规范应力足以确定这一局部系数”的合同；不声称局部层析失效、完整作用积分不同或该混合项已被证明是独立可观测耦合。

这些局部资料不是不能实现的任意矩阵。对给定Abelian f及目标正交标架e，取

$$
a_{0\mu}(x)=-\frac12f_{\mu\nu}x^\nu,\qquad
\phi(x)=\phi_0+V_\mu x^\mu,\quad V_\mu=eX_\mu.
\tag{8}
$$

原点a0=0、Dφ=V、da0=f；足够小邻域保持F>0。以−f替换亦可。代码直接从这个局部连接差分计算交换子，独立核式(5)。它没有强迫这两份局部资料满足完整Gauss、物质或Einstein方程，也没有完成572的实际来源转换。

## 6. 规范应力必须进入同一压缩

以κ=1的共同Euclidean领先作用及总规范应力定义EFT冗余：

$$
I_0^{g}=\int\sqrt g\left[-\frac R2+\frac S2+U+\frac14\sum_A w_A(f^A)^2\right],\qquad
\tau_{\mu\nu}=\sum_Aw_A\left(f^A_{\mu\rho}f^A_\nu{}^\rho-\frac14g_{\mu\nu}(f^A)^2\right),\quad
\operatorname{tr}_g\tau=0.
\tag{9}
$$

tracefree是四维领先Yang–Mills性质，不是宣称完整量子无迹反常。由度规Euler导数有

$$
E^g_{\mu\nu}=R_{\mu\nu}-B_{\mu\nu}-Ug_{\mu\nu}-\tau_{\mu\nu},\quad
e=R-S-4U,\qquad
H^g_{\mu\nu}=\frac1{\sqrt g}\frac{\delta I_0^g}{\delta g^{\mu\nu}}
=-\frac12(E^g_{\mu\nu}-e g_{\mu\nu}/2).
\tag{10}
$$

因此Ric=B+Ug不是当前共同源的领先关系；遗漏τ会在同一阶漏交叉项。记b4,red为581式(8)，则完整约化候选为

$$
\boxed{b_{4,\mathrm{red}}^{g}=b_{4,\mathrm{red}}+\Delta_g b_4
+\frac16 B_{\mu\nu}\tau^{\mu\nu}+\frac1{12}\tau_{\mu\nu}\tau^{\mu\nu}.}
\tag{11}
$$

证明不需要预设E^g=0。令

$$
D^g_{\mu\nu}=\frac{R_{\mu\nu}+B_{\mu\nu}+Ug_{\mu\nu}+\tau_{\mu\nu}}{12}
+g_{\mu\nu}\left[\frac{R+S+4U}{24}-\frac{\operatorname{tr}A}{6}-\frac S9\right],\quad
T^g=-2D^g+(\operatorname{tr}_gD^g)g.
\tag{12}
$$

展开Ric²的平方差；trτ=0使新增部分恰为B:τ/6+τ²/12，R关系不变。原U规范不变且w_A为常数，标量链式法则仍成立，所以

$$
b_4^{g}=b_{4,\mathrm{red}}^{g}+H^g:T^g
-\frac16\langle\operatorname{grad}U,E_\phi^g\rangle_K-\frac16\Box U,
\quad E_\phi^g=\operatorname{grad}U-\mathscr D^\mu D_\mu\phi.
\tag{13}
$$

度规与标量共同替换仍为

$$
g_{\rm old}^{\mu\nu}=g_{\rm new}^{\mu\nu}-\lambda(T^g)^{\mu\nu},\qquad
\phi_{\rm old}=\phi_{\rm new}+\lambda\operatorname{grad}U/6,\qquad
a_{\rm old}=a_{\rm new}.
\tag{14}
$$

gradU按原群等变，因而不破坏表示结构。无需使用Yang–Mills EOM；仍须保留全作用的规范背景，因Dφ及度规变分已经含其作用。与581相同，仅是一阶微扰、有限系数和有界光滑资料下的局部形式逆；源、观察量及边界必须共同拉回，未证明完整量子处方或非扰动等价。

## 7. 六组实际复算

[代码](582/joint_gauged_effective_geometry.py)、[结果](582/joint_gauged_effective_geometry_results.json)复用原Python与NumPy，默认完整复算比较，不覆盖。两个入口检查首次纳入正式计数，另加四组独立验证。

|检查|结果|
|---|---|
|原Abelian目标连接|荷3、Killing反对称误差2.22×10⁻¹⁶，trJ²=−36.7952605201|
|同应力局部见证|f与−f连接系数−.68182693069、−.71822360814，应力相同；Δmix=.01819833872|
|全电弱Gram与同一质量|三种场点的式(3)误差≤2.36×10⁻¹⁵；共同真空标准化质量平方约0、.024799208、.024799208、.032392908，均含1/F|
|原群规范协变|完整弱变换后新增系数误差≤5.33×10⁻¹⁵，混合项单独保持|
|共同规范应力压缩|原权w_w=2.3788490133、w_0=279.675086406；三份局部资料式(13)的度规部分误差≤5.56×10⁻¹⁷；领先度规例若漏新增应力项，局部系数差.001729987945|
|连接交换子独立差分|步长2×10⁻⁴、10⁻⁴、5×10⁻⁵的误差1.87×10⁻¹⁰、4.67×10⁻¹¹、1.17×10⁻¹¹，符合直接局部连接实现|

一般局部张量核验不等于构造完整在壳解。规范源中的更多物种、全部环、场重定义后的约束相空间和实际量子记录尚未完成。

## 8. 对共同条件账的改变

[582条件账](582/unified_physics_condition_ledger_582.md)把规范质量、连接量子修正和共同应力放进同一明确归一。式(3)减少了把场依赖规范系数任意单独选择的余地；式(11)排除了“保留非零规范来源却原样沿用纯标量压缩”的接法。

没有证明所列局部项均为独立实验参数，也没有把由额外物理作用推得的关系称为新的认知公理。580—582仍只含标量内部环；553固定Jordan、574固定Einstein的量子变量匹配并未因经典或EFT变量替换自动解决。

下一步优先核这个量子处方接口：哪些变量必须联合量子化，哪些可以在同一受控近似中作为背景，来源与有效作用怎样共同匹配。接入成熟的框架等价及联合Hessian方法，先核前提，避免继续以孤立规范算符扫描代替整体推进。统一目标、连续尺度和物理预测仍开放。
