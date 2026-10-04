# 第643轮：完整Gauss物质、共同热参考与来源的有序历史表示

日期：2026-10-01。接[642](../../research_note_642.md)。[代码](../joint_gauss_fermion_influence.py)、[结果](../joint_gauss_fermion_influence_results.json)、[核验](../research_round_643_checks.json)、[条件总账](../unified_physics_condition_ledger_643.md)。三组检查、十四式；主代理审查，无新增独立代理审查。

## 1. 从限定反例转向共同构造

642说明，原完整模型的内部CAR不能迹掉后只靠剩余等时态重建来源。本轮给一条可用的替代：**保留被消去物质在整个历史上的有序作用，并将Gauss投影放在同一历史的闭合处。** 原完整非线性有限图满足这种热核表示的适用条件；正常热参考、费米消元、原标量记录和全部几何来源可共用该表示。

没有把“可能存在影响泛函”当作结果。下文给出原曲目标与非对角电动能上的实际积分对象、矩阵势可积界、Gauss闭合因子以及原32模式的精确条件权重。原整体相互作用未作Gaussian化；只有**给定玻色历史后**的原有限CAR生成元为二次型。

这提供共同对象及精确表示，不提供廉价算法、空间连续极限、手征重建或GR。数值只核原单节点完整32模式的条件历史因子；没有执行全图玻色路径积分、Gauss平均或计算完整物理配分函数。

## 2. 成熟工具、对象映射与去重

[Güneysu，定理2.11及2.13](https://arxiv.org/html/1109.0151v3)给向量值Schrödinger半群的Feynman–Kac表示和标量支配；所需矩阵正部局部可积、负部满足Kato条件。本项目采用平凡有限Fock纤维，原动能作为配置流形Laplace型算子；下节核实完整矩阵势的负部有界，因此无需另假定质量项有界。

[Klich，§II、§IV式(37)](https://arxiv.org/html/1403.7824)给含配对的有限CAR指数乘积迹公式，并明确平方根符号须由原Fock提升与解析性确定。本轮使用这一成熟迹工具，保留原质量、原商群和闭合插入；不将它称为新数学定理。

历史上，[529](../../../archive_467_530/research_note_529.md)是另一受限关系模型的标量路径界；[603](../../../archive_585_628/research_note_603.md)只用标量比较证明全热迹有限，未把原矩阵势保留为精确路径权；[615](../../../archive_585_628/research_note_615.md)计算的是另一手征候选的辅助Pfaffian测度。[622](../../../archive_585_628/research_note_622.md)—[625](../../../archive_585_628/research_note_625.md)处理完整过程的二阶来源及有限近似。本轮新增的是**原完整非线性Gauss热态的精确有序费米消元与同一来源表示**，不是复述路径积分或再做Gaussian记录扫描。

|层次|本轮地位|
|---|---|
|认知动机|物质被消去后，参考、记录及反作用仍应来自同一过程|
|继承输入|原固定有限三维图、H⁵目标、G商群、32CAR／节点、Y及跳跃、严格正外部γ、β>0|
|既有结果|598点态质量界、603全热迹、623共同域；全部常数可依赖图|
|解析增量|原全模型的矩阵桥表示、Gauss扭闭合、绝对可积及来源／记录共用|
|数值验证|完整原32模式条件因子；规范协变、排序、配对及条件来源插入|
|未证明|宇宙参考选择、连续手征、统一物理极限、动态引力与自主认知设计|

## 3. 原完整矩阵势满足适用条件

记q=(φ,U)，配置空间Q=(H⁵)^V×G^E，测度沿642式(1)。固定正γ时，589的Tγ包含全部曲目标动能及原非对角电动能；其配置度量完整。紧群各出边簇具有正右不变度量，体积与原Haar乘积只差常数，核以下均相对于**原固定dμ**定义。取ℏ为既有正值，令扩散生成元为−Tγ，不擅自增减1/2。

$$
H_\gamma=T_\gamma\otimes I+V_\gamma(q)I+B_\gamma(q),\qquad
V_\gamma=W_\gamma+V_{{\rm grad},\gamma}+V_{{\rm mag},\gamma}\ge W_\gamma,
\qquad W_\gamma=\sum_vw_vU(\phi_v).
\tag{1}
$$

B含598原Dirac、Majorana及规范协变跳跃，Fock维数d_F=2ⁿ，n=32|V|（一代）。它是关于q的光滑有限矩阵，但靠近F=0可以无界。598／603的点态界给

$$
\|B_\gamma(q)\|\le\tfrac12 W_\gamma(q)+C_\gamma,\qquad
\mathcal V_\gamma(q):=V_\gamma(q)I+B_\gamma(q)
\ge(\tfrac12W_\gamma-C_\gamma)I.
\tag{2}
$$

于是完整矩阵势的负部全局有界，正部在配置流形内部光滑且局部可积。常数属于标量Kato类，满足引用定理的前提。只检验B而忘记原束缚势会误判适用性。F=0在完整H⁵目标的无穷远，有限时间连续扩散路径的像几乎处处紧；原正动能对应的扩散不爆炸。所有几何只在固定正紧参数域内比较。

## 4. 精确矩阵桥与绝对可积

用Dγ,β^{y→x}表示−Tγ从y到x的**未归一**桥测度，其总质量是动能核kγ,β(x,y)。路径中的矩阵演化按晚时在左的约定定义：

$$
\frac{d\mathcal U_F(t)}{dt}=-B_\gamma(q_t)\mathcal U_F(t),\qquad
\mathcal U_F(0)=I,\qquad
\mathcal U_F[q]=\mathcal T\exp[-\int_0^\beta B_\gamma(q_t)dt].
\tag{3}
$$

每条非爆炸路径上系数有界，故该有限矩阵ODE正常定义。真正的全H热核为

$$
K_{\gamma,\beta}(x,y)=\int D_{\gamma,\beta}^{y\to x}(dq)\,
e^{-\int_0^\beta V_\gamma(q_t)dt}\mathcal U_F[q].
\tag{4}
$$

原非线性标量、所有链路和跳跃均仍在积分中；没有固定Higgs方向、删除磁势或把玻色态换成经典背景。规范协变来自原T、V的等变性与B(gq)=Γ(g)B(q)Γ(g)†。

Gronwall及式(2)给逐路径支配

$$
e^{-\int V_\gamma}\|\mathcal U_F[q]\|
\le e^{\beta C_\gamma}e^{-\frac12\int W_\gamma(q_t)dt},\qquad
\|K_{\gamma,\beta}(x,y)\|
\le e^{\beta C_\gamma}k^{(1/2)}_{\gamma,\beta}(x,y),
\tag{5}
$$

其中k^(1/2)是Tγ+Wγ/2的标量核。603的原势尾和热迹比较保证它对角可积。完整Fock有限这一事实控制维数，但不是忽略无界质量的理由。

令局部规范群𝒢=G^V，Γ(g)为原全Fock作用。群作用约定为(U(g)Ψ)(q)=Γ(g)Ψ(g⁻¹q)。原Gauss热迹遂化为

$$
Z_{\rm phys}(\beta;\gamma)
=\int_{\mathcal G}dg\int_Qd\mu(q)
\int D_{\gamma,\beta}^{q\to g^{-1}q}(d\omega)\,
e^{-\int V_\gamma(\omega_t)dt}\,
\underbrace{\operatorname{Tr}_{\mathcal F}[\Gamma(g)\mathcal U_F[\omega]]}_{\mathcal W_F[\omega,g]}.
\tag{6}
$$

g既决定玻色桥的端点，也出现在同一费米迹中。不能只保闭环q→q，再把费米迹预先投成各点独立singlet；也不能只删除Γ(g)而保原端点。

标量核的半群Cauchy–Schwarz和规范不变性给k^(1/2)(g⁻¹q,q)≤k^(1/2)(q,q)，所以整个式(6)的绝对值积分满足

$$
\int dg\,d\mu(q)\int D^{q\to g^{-1}q}
e^{-\int V_\gamma}|\mathcal W_F|
\le d_F e^{\beta C_\gamma}
\operatorname{Tr}e^{-\beta(T_\gamma+W_\gamma/2)}<\infty.
\tag{7}
$$

因此矩阵迹、规范平均和桥积分可以按此支配交换。原P_G与H对易，603又给物理空间非零，最终0<Z_phys<∞。证明不要求每条条件路径权重都是实的非负概率；这不是把整个量子系统改成经典随机模型。

## 5. 保留原配对、闭合及历史排序的有限矩阵字典

沿固定玻色历史，B始终为CAR二次型，尽管整体耦合模型不是自由场。以分段常量近似路径系数，记h_j、Δ_j为原单粒子块，Δᵀ=−Δ，允许z_j为实热时或复轮廓参数。Nambu矩阵与闭合矩阵为

$$
\mathsf H_j=\begin{pmatrix}h_j&\Delta_j\\\Delta_j^\dagger&-h_j^T\end{pmatrix},\qquad
\mathsf R(g)=\operatorname{diag}(R_F(g),R_F(g)^*),\qquad
\mathsf M=\mathsf R(g)e^{-z_m\mathsf H_m}\cdots e^{-z_1\mathsf H_1}.
\tag{8}
$$

当前原h仅含跨模块质量与不同节点跳跃，Tr h_j=0；完整32模式R_F的行列式为1，见642。按原Γ(g)的连续Fock提升选择迹，Klich公式给

$$
\left(\operatorname{Tr}_{\mathcal F}
[\Gamma(g)e^{-z_mB_m}\cdots e^{-z_1B_1}]\right)^2
=\det(I_{2n}+\mathsf M),\qquad
\left(2^{-n}\mathcal W_F\right)^2=\det\frac{I_{2n}+\mathsf M}{2}.
\tag{9}
$$

此式先在有限分段成立，再按矩阵ODE的收敛取连续路径极限。它不是允许每条路径任取主平方根：符号／相位由原Fock算符乘积固定，在零点附近应延拓实际迹。若另加有迹的h或改变真空常数，需补相应标量因子，不能直接沿用式(9)。本轮没有改变原排序或真空能源。

给定源a，在𝒲_F≠0且矩阵可逆的局部可写

$$
\partial_a\log\mathcal W_F
=\tfrac12\operatorname{tr}[(I+\mathsf M)^{-1}\partial_a\mathsf M],\qquad
\partial_a\mathcal U_F(\beta,0)
=-\int_0^\beta\mathcal U_F(\beta,t)(\partial_a B_t)
\mathcal U_F(t,0)dt.
\tag{10}
$$

零权重点用右侧无除法插入式定义导数。源若影响多段，所有段均须求导；非对易乘积不可替换成质量平均后的单个指数。

## 6. 全部来源和实际记录的共同范围

式(10)只给条件费米贡献。空间γ还改变原Tγ、Vγ、可能的跳跃系数及桥核，不能把条件费米得分当完整应力。尤其，不把不同扩散系数的Brownian测度未经证明地视为有可微概率密度比。

全部来源按623共同域、603热矩控制下的Duhamel式处理，P_G与γ无关：

$$
\partial_a Z_{\rm phys}
=-\int_0^\beta d\tau\operatorname{Tr}\big[
P_Ge^{-(\beta-\tau)H}(\partial_aH)e^{-\tau H}\big]
=-\beta Z_{\rm phys}\langle\partial_aH\rangle_{\beta,\rm phys},\qquad
\partial_aH=\partial_aT+\partial_aV+\partial_aB.
\tag{11}
$$

微分型来源必须保相应核导数或形式插入；式(6)若仅微分𝒲_F会漏项。该表示接通参考与来源，但没有给动态几何或Einstein约束。

原标量仪器L_r(s)=√(.5+r sin s/4)为有界规范不变乘法，作用在Fock上为恒等。在Euclidean算符字中，它直接作为桥上的标量插入，与同一𝒲_F相乘。这类Euclidean关联本身不能叫实际测量概率；真实历史仍用

$$
C_{\boldsymbol r}=L_{r_m}e^{-it_{m-1}H}\cdots L_{r_2}e^{-it_1H}L_{r_1},\qquad
p_{\boldsymbol r}=Z_{\rm phys}^{-1}
\operatorname{Tr}[P_G C_{\boldsymbol r}e^{-\beta H}C_{\boldsymbol r}^\dagger].
\tag{12}
$$

它由同一H、P_G和β确定。对每条实时间段先加正实部ε，则热半群的算符解析延拓及迹类连续性给

$$
e^{-(\varepsilon+it)H}\xrightarrow[\varepsilon\downarrow0]{\rm strong}e^{-itH},\qquad
\operatorname{Tr}[P_G C_{\boldsymbol r,\varepsilon}e^{-\beta H}
C_{\boldsymbol r,\varepsilon}^\dagger]
\longrightarrow Z_{\rm phys}p_{\boldsymbol r}.
\tag{13}
$$

下界移位使规整因子一致有界，L有界且e⁻ᵝᴴ迹类，故极限成立。实时端不主张普通Wiener概率测度；有限CAR复轮廓因子则仍可由式(8)—(9)求取。本轮没有数值完成整个式(12)。若改为几何依赖的区域仪器，还需639的显式仪器导数。

## 7. 原32模式的可复算检验

数值取原单节点的四组内部标量值，均满足F>0，完整复Y及原Majorana保留。给定一份非平凡原C、W、z闭合元素。因为这一**条件单节点因子**的夸克部门数守恒且与轻子部门分离，可独立计算

$$
\mathcal W_F=
\det[I_{24}+R_q(g)e^{-z_mh_{q,m}}\cdots e^{-z_1h_{q,1}}]
\operatorname{Tr}_{\Lambda\mathbb C^8}
[\Gamma(R_l(g))e^{-z_mB_{l,m}}\cdots e^{-z_1B_{l,1}}].
\tag{14}
$$

轻子迹使用实际256维Fock矩阵；Γ(R_l)通过外代数子式直接构造，独立于64维Nambu计算。夸克24模式全保留，未用一个中性16维因子替代整代。原全图跳跃在解析式(6)—(10)保留；数值未计算有空间跳跃的多节点因子。

所有权重除以2³²只为显示和数值稳定，不是物理配分函数归一。四段Euclidean时长为(.4,.3,.5,.25)，另一复轮廓为(.7,.4i,−.6i,.25i)：

|条件因子或改法|归一权重|
|---|---:|
|完整原Euclidean有序因子|.754683209868|
|只省掉闭合Γ(g)|1.377037153935|
|改成积分质量的单个指数|.731154764743|
|删除原Majorana|.752250309003|
|完整原复轮廓因子|.377920862696＋.112617834303i|

三种改法分别改变条件权重约.622354、.0235284和.00243290。不是完整Gauss积分的误差下界，也不是物理实验概率。原Nambu不同片段的对易子范数最大约.683092，排序影响真实存在。

共同变换全部φ和g→kgk⁻¹，权重误差低于5.5×10⁻¹⁵；Z₆ lift残差低于1.1×10⁻¹⁵。完整32模式迹平方与64维行列式残差低于1.3×10⁻¹⁴。另以原四模式实际Fock确认配对惯例，残差低于1.8×10⁻¹³。

将每段B_j乘eᵃ，只检查**条件费米lapse贡献**：a=0时权重导数约.527304505130，对数导数约.698709734410。矩阵插入与行列式导数差低于5.6×10⁻¹⁶；直接差分步长从2×10⁻⁴减到10⁻⁴，误差由3.11×10⁻⁸降至7.78×10⁻⁹。完整lapse还会改变T、V，见式(11)，本表不是其绝对来源。

## 8. 对联合条件的实际贡献

本轮将C01量子对象、C14／C15的原规范与CAR、C19参考、C03记录和C22来源放进同一份有序历史表示。**不用再为被消去的原物质独立猜测一份噪声或记忆作用；它已由原B、Gauss闭合与同一玻色桥确定。** 代价是保留其历史依赖和原完整积分，不能据此宣称获得局部、无记忆的粗Hamiltonian。

原规范群、物种、Y、三维图、正γ与β仍是输入；式(6)不是连续手征路径测度，也没有补足615的一般局域测度问题。固定图存在性不能替代统一的尺度界。全部27项已回填[条件总账](../unified_physics_condition_ledger_643.md)，整理本身不另计轮次。

接[644入口](../../644/drafts/STATUS.md)：利用同一表示检查区域拼接／区域态与几何来源的共同信息，优先对接617、636—639及635的几何熵条件。若仅重述一般复制技巧或再改一组质量参数，则不计新轮次；不继续条件因子的精度扫描，不进入认知系统设计。
