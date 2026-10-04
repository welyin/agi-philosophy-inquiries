# 第618轮：原非最小耦合、几何边界与标量来源的共同匹配

日期：2026-10-01。接[617](research_note_617.md)及[618入口](618/drafts/STATUS.md)。[代码](618/joint_geometric_boundary_matching.py)、[结果](618/joint_geometric_boundary_matching_results.json)、[核验](618/research_round_618_checks.json)、[条件账](618/unified_physics_condition_ledger_618.md)。三组复算、十四式；主代理审查，无新增独立代理审查。

## 1. 本轮关掉的具体缺口

617保持的是同一有限图的规范区域表示；它没有给动态几何的区域拼接。[358](../archive_342_369/research_note_358.md)明确保留有限边界，[583](../archive_554_584/research_note_583.md)使用周期或已处理边界，[600](research_note_600.md)—[601](research_note_601.md)的共同作用论证也没有补出具体边界变分。

本轮将**原五标量目标、非最小引力耦合、区域边界动量和标量来源**接通：

1. 原F同时决定几何边界动量和标量法向通量，二者不能独立选。
2. 原F给共同边界混合块的Schur量恒为M=2；在F>0分支，若没有表面来源，连续边界值还必须配合连续法向一阶数据。不能仅凭两边的度规、标量值相同就签收无源拼接。
3. Jordan—Einstein转换同时变换边界作用及标量来源。明确的表面张力诊断显示：Jordan中不显含标量的表面作用，在Einstein中必有标量来源；遗漏它会破坏同一变分。
4. 使用原势的非周期实际作用验证上述边界项与变分，未靠周期边界令所有差异消失。

**结论限给定二导数标量—引力作用的经典类时界面。** Einstein项、四维、原F及势都是继承输入；不是从认知原则生成引力，也未完成规范／费米边界或量子拼接。

## 2. 成熟结果、对象映射与保留输入

[Avilés、Maeda、Martínez，2019](https://arxiv.org/pdf/1910.07534)在§3及附录B给标量—张量非类光界面的接合与变分；其单场记号fR对应本轮FR/2。[Bhattacharya、Bamba，2023](https://arxiv.org/pdf/2307.06674)讨论两种框架的良定边界项和Brown—York数据。成熟junction理论和GHY边界项不算本项目新发现。

本轮保留原五场的非平直Einstein目标，直接推导五场边界混合，不套单场的全局canonical坐标。真正新增是将上述条件落实到574／583的原F与K、验证其特殊非退化恒等式，并补出有限边界的共同来源。

|层次|本轮采用|
|---|---|
|认知动机|同一整体切成区域时，不应凭切法增添未记账来源|
|继承输入|四维、领先二导数作用、M=2、原势V、F>0及无挠联络|
|边界约定|签名−+++，光滑类时Σ，向外单位法向n²=+1，K_ab=h_a^μh_b^ν∇_μn_ν|
|场范围|原五标量与度规的变分；规范／费米场在本轮诊断中取零，不签收其边界问题|
|额外诊断|一项明示测试用表面张力，检验表示变换；不把它设成宇宙新物质|
|未完成|类光边界、角点、所有物种、完整EFT边界反项、量子态与尺度映射|

假设界面两侧诱导度规和φ连续且分片足够光滑，变分沿Σ切向紧支撑以避开角点。下文使用两侧各自的向外法向之和，不使用未说明取向的跳跃符号。无源条件排除的是界面δ来源；体内方程、约束和全局解存在仍需另核。

## 3. 原目标与带边界的共同作用

令f_A=∂_A F，场指标用δ_AB升降。原模型为

$$
F=M-\frac{\phi\cdot\phi}{6}>0,\qquad f_A=-\frac{\phi_A}{3},\qquad
\mathcal K_{AB}=\frac{\delta_{AB}}F+\frac{3f_Af_B}{2F^2},\qquad
U=\frac V{F^2},\quad M=2.
\tag{1}
$$

K正是574完整双曲目标；并未换成新的平直场空间。在当前单位g_E=F g_J，Lorentz领先作用连同类时边界写为

$$
S_J=\int_{\mathcal M}\sqrt{-g}\,
 \left[\frac12 FR-\frac12\delta_{AB}\partial\phi^A\cdot\partial\phi^B-V\right]
 +\int_{\partial\mathcal M}\sqrt{|h|}\,F K ,
\qquad
S_E=\int_{\mathcal M}\sqrt{-g_E}
 \left[\frac12 R_E-\frac12\mathcal K_{AB}\partial\phi^A\cdot\partial\phi^B-U\right]
 +\int_{\partial\mathcal M}\sqrt{|h_E|}\,K_E .
\tag{2}
$$

V直接继承原二标量势及全部常数。F K边界项来自给定FR作用的Dirichlet变分，不是额外独立耦合；暂不加参考背景减项或真实表面物质。两侧分别加入此项，再共同变化界面数据，才得到拼接条件。

## 4. 两种边界动量必须一起匹配

取协变h_ab为边界变量，设v_A=n^μ∂_μφ_A。共同变分的边界一形式为

$$
\delta S_J\big|_{\partial\mathcal M}
 =\int \left(\Pi_J^{ab}\delta h_{ab}+\pi_{J,A}\delta\phi^A\right)
 =\int\sqrt{|h|}\left(\frac12P^{ab}\delta h_{ab}+j_A\delta\phi^A\right).
\tag{3}
$$

对曲率项分部积分，F K消去法向度规变分的导数；剩余系数为

$$
P^{ab}=F(Kh^{ab}-K^{ab})+h^{ab}f_Av^A,\qquad
j_A=f_AK-v_A,\qquad
\Pi_J^{ab}=\frac{\sqrt{|h|}}2P^{ab},\quad \pi_{J,A}=\sqrt{|h|}\,j_A.
\tag{4}
$$

j中的f_AK来自F K的标量变分；P中的nF来自非恒定F的度规变分。丢掉其中一项会破坏同一作用的变分。符号与改用逆度规变量时相反，不能混用两种记号。

两侧共同h、φ下，定义κ_ab=K_ab^++K_ab^-、κ=h^{ab}κ_ab、ν_A=v_A^++v_A^-，并令p^{ab}=P_+^{ab}+P_-^{ab}、j_A=j_{+,A}+j_{-,A}、p=h_ab p^{ab}。则

$$
p^{ab}=F(\kappa h^{ab}-\kappa^{ab})+h^{ab}f\cdot\nu,\qquad
p=2F\kappa+3f\cdot\nu,\qquad
j=f\kappa-\nu.
\tag{5}
$$

有表面作用时，该p、j须抵消表面作用的对应变分；无表面作用且界面数据共同自由变化时，p^{ab}=j=0。若将界面当作两个永不共同变化的固定外边界，则没有得到无源整体拼接的变分问题。

## 5. 原F使混合块可逆，不允许无来源的一阶折角

消去ν后，唯一混合分母为2F+3|f|²。原模型恰满足

$$
F+\frac32|f|^2=M,\qquad
\kappa=\frac{p+3f\cdot j}{2M},\qquad
\nu=f\,\frac{p+3f\cdot j}{2M}-j .
\tag{6}
$$

无迹部分再给出

$$
\kappa_{\rm TF}^{ab}=-\frac{p_{\rm TF}^{ab}}F,\qquad
p^{ab}=0,\ j=0\ \Longrightarrow\ \kappa_{ab}=0,\ \nu_A=0 .
\tag{7}
$$

于是，在共同适配的Gaussian法向坐标中，h及φ连续再加无来源界面方程，就要求法向一阶数据匹配。反之匹配这些数据不会产生这类界面δ来源；这并不保证体内方程或全局解存在。

式(6)没有新参数：与原目标K的正性及逆矩阵使用同一个F、M。若以法向共形速率w及五场速率v表示，法向动能Hessian为块矩阵[[6F,3fᵀ],[3f,−I₅]]，行列式恒−6M；这说明非退化，**不是**正性或量子引力稳定性。完整无迹块还依赖F，故不能将F=0也纳入结论。

成熟标量—张量理论允许某些退化耦合上的无源薄壳；本轮原正F分支排除相应混合退化，并没有反证文献中的其他分支，也不外推到类光界面或高阶引力。

## 6. 共同框架转换包含边界标量来源

保留五场坐标不变。沿同一向外法向，Weyl转换给

$$
h_{E,ab}=Fh_{ab},\qquad
K_{E,ab}=\sqrt F\left(K_{ab}+\frac{h_{ab}}{2F}f\cdot v\right),
\qquad v_E=F^{-1/2}v,\qquad
\sqrt{|h_E|}=F^{3/2}\sqrt{|h|}.
\tag{8}
$$

将式(8)代入Einstein动量，并使用式(1)的K，可直接核共同边界一形式：

$$
\Pi_J^{ab}=F\Pi_E^{ab},\qquad
\pi_{J,A}=\pi_{E,A}+f_A h_{ab}\Pi_E^{ab},\qquad
\Pi_J^{ab}\delta h_{ab}+\pi_{J,A}\delta\phi^A
=\Pi_E^{ab}\delta h_{E,ab}+\pi_{E,A}\delta\phi^A .
\tag{9}
$$

第二式的混合项不能删除。原因是“固定Jordan诱导度规时改变φ”会改变Einstein诱导度规。583已核体内固定变量的差异；本轮补的是其有限边界的实际共轭数据。

### 有来源诊断：原框架中的常张力并非另一框架中的常张力

只作检验，另给测试表面作用

$$
S_\Sigma=-\tau\int_\Sigma\sqrt{|h_J|}
=-\tau\int_\Sigma\sqrt{|h_E|}\,F^{-3/2}.
\tag{10}
$$

τ是此诊断的新输入，不用于声称底层必有膜或新增认知资源。Jordan中表面变分为P_Σ,J^{ab}=−τh_J^{ab}、j_Σ,J=0。Einstein中则为

$$
P_{\Sigma,E}^{ab}=-\tau F^{-3/2}h_E^{ab},\qquad
j_{\Sigma,E,A}=\frac{3\tau}{2}F^{-5/2}f_A .
\tag{11}
$$

需要由体侧抵消的Jordan数据是p^{ab}=τh^{ab}、j=0。式(6)—(7)唯一给

$$
\kappa=\frac{3\tau}{2M},\qquad
\kappa_{ab}=\frac{\tau}{2M}h_{ab},\qquad
\nu_A=\frac{3\tau}{2M}f_A,\qquad
j_{E,A}=-\frac{3\tau}{2}F^{-5/2}f_A .
\tag{12}
$$

所以Jordan表面不显含φ，仍可能因非最小耦合要求φ的法向导数跳变；在Einstein变量中，同一效应由表面标量来源显式表示。设j_Σ,E=0会漏来源，不能凭相同“张力”名称视为同一物理边界。这只是局部接合数据及表示变换，未构造满足全部体方程的完整薄壳时空。

## 7. 非周期实际作用验证

为避免只在同一公式的两种写法间比对，代码从两框架的原曲率和GHY分别计算作用。取法向x∈[0,1]，每单位切向坐标体积，g_J=dx²+e^{2a(x)}η_ab dy^a dy^b，η=diag(−1,1,1)，a及五场均为非恒定二次多项式；端点不周期。R_J=−6a″−12a′²，向外K分别为−3a′(0)、+3a′(1)。直接分部积分应给

$$
S_J^{\rm bulk+GHY}
=\int_0^1 e^{3a}\left(3F(a')^2+3F'a'
-\frac12|\phi'|^2-V\right)dx
=S_E^{\rm bulk+GHY}.
\tag{13}
$$

g_E的法向lapse为√F，代码用其独立曲率表达式积分，保留原V，不通过目标答案生成S_E。一般这些剖面是离壳的；没有将任意五场剖面称为完整规范物质方程的解。

对端点也非零的线性变化δa、δφ，独立的变分检验是

$$
\delta S
=\int_0^1\left(\mathcal E_a\delta a+\mathcal E_A\delta\phi^A\right)dx
+\left[e^{3a}\left((6Fa'+3F')\delta a
 +(3f_Aa'-\phi'_A)\delta\phi^A\right)\right]_0^1 .
\tag{14}
$$

体Euler项由该原作用直接微分；变分左侧对原曲率加GHY作中心差分。64／96点Gauss积分核积分精度，三个变分步长核二阶收敛。

|三组复算|保存结果|
|---|---|
|五场混合、反解及框架边界一形式|四个原半径位置，含F=0.195；Schur量恒2，六维Hessian行列式−12；最大残差7.11×10⁻¹⁵|
|非周期完整作用与变分|S_J=S_E=−0.000517877225632，差小于1.53×10⁻¹⁶；漏GHY时两体作用差0.06536904006；变分误差6.13×10⁻⁸、1.53×10⁻⁸、3.83×10⁻⁹，按步长平方收敛|
|原真空的来源诊断|F₀=1.8765364152，τ=.07时κ=.0525；转换后的标量表面来源范数0.006244829989，不能置零；共同匹配残差6.94×10⁻¹⁸|

只保持h、φ连续而取κ_ab=.04h_ab、ν=0，会留下度规动量差0.2600205131及标量动量差0.0344274660。它是错误拼接数据的局部变分见证，不是两份已解完整体方程的全局反模型。

## 8. 条件压缩与边界

本轮把C02区域拼接、C08目标几何、C10共同框架、C12给定引力作用、C22来源的一部分条件归为**同一完整带边界作用的共同变分**。原F确定混合系数与可逆性；不必再独立规定一套与之无关的标量拼接法则。

C11约束得到的只是必要界面一致性，不是完整约束代数、初边值问题或量子物理态。规范与费米表面项、高阶EFT边界项、类光面、角点和实际区域操作仍未覆盖。617固定图酉拼接与本轮经典时空接合属于不同层面的对象，尚无共同连续尺度映射将二者等同。

[619入口](619/drafts/STATUS.md)继续核原规范／费米物质的边界变分怎样加入同一结构，优先明确场、法向电流、作用与来源的共享条件；回查既有手征与spin范围，不把形式边界数据当作认知控制装置。总目标与“先整合条件、再研究设计”的顺序不变。

