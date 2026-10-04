# 第587轮：同一量子物质的能源流、切向来源与记录涨落

日期：2026-10-01。接[586](research_note_586.md)、[路径积分方法对接](586/closed_time_path_bridge_review_586.md)，完成原[能源流入口](587/drafts/STATUS.md)。[代码](587/joint_matter_energy_current.py)、[结果](587/joint_matter_energy_current_results.json)、[核验](587/research_round_587_checks.json)。主代理完成推导与核验，未取得新的独立代理审查。

## 1. 问题、输入及实际增量

同一物质既承担记录又作为几何来源时，只核总能源不够。本轮问：原完整H能否给出与既有切向来源相容的能源流，并同时保留真实记录过程的反作用？

**结论：** 在585已选局部能源分配下，原曲目标、测地边及全部商群磁势给一份保Gauss、有限支撑的一阶量子流。在固定光滑源上，其经典弱流共同趋向572原物质的切向动量来源，无需另加标量与规范两套流模型。574物理波包将固定图量子期望接到同一经典符号，但不提供固定ℏ的量子连续极限。586两种同报告仪器都保持瞬时平均流，流的二阶量却一般不同。

|层次|本轮地位|
|---|---|
|认知动机|记录、资源转移与整体来源应属于同一过程，平均报告和平均流不能替代全部关联。|
|继承输入|574完整H、569群及系数、585节点／端点／面能源分配、586实际仪器。|
|连续比较输入|给定三维规则周期图、固定光滑正ψ、正F紧范围、光滑经典源及同一规范平凡化；单位lapse、零shift。|
|解析增量|完整曲目标／规范流，原来源的弱连续匹配，固定图半经典连接，实际仪器的流二阶差。|
|未完成|量子应力全部分量、局部几何耦合唯一性、量子形变代数、动态几何及自主仪器。|

### 去重及成熟成果

[318](../archive_301_341/research_note_318.md)已有连续改进应力问题；[352](../archive_342_369/research_note_352.md)已有给定形变代数下的条件性共同传播；[572](../archive_554_584/research_note_572.md)已定义本模型连续动量密度；[585](research_note_585.md)保留lapse接法不唯一性。通用对易子连续性本身不是新定理。

[Boguslavski、Lappi、Peuron、Singh，EPJC 84, 368（2024）](https://pmc.ncbi.nlm.nih.gov/articles/PMC11001687/)从实际格点Hamiltonian构造经典Yang–Mills能源流，强调直接套连续应力不保证离散守恒，并区分能源与动量连续性的精度。本轮采用这一检查思路；其纯规范Wilson作用结果不能替代当前曲目标物质与忠实商群磁势的推导，更不能当作量子应力重整化定理。

## 2. 从完整H和既定分配得到流

写H=Σ_A T_A+Σ_Z W_Z。A遍历节点及全部规范链路因子，Z遍历现场、测地边和磁面。各因子度量沿574、几何权固定；不同T因子对易，W均为配置乘法。在紧支撑光滑规范不变共同核上，定义

$$
h_i=\sum_A a_{iA}T_A+\sum_Z b_{iZ}W_Z,\qquad
\sum_i a_{iA}=\sum_i b_{iZ}=1,\qquad
D_AW_Z=\frac{i}{\hbar}[T_A,W_Z]
=-i\hbar\left(\nabla_AW_Z\cdot\nabla_A+\frac12\Delta_AW_Z\right).
\tag{1}
$$

梯度／Laplacian用原加权因子度量。节点T及现场W全归本节点，边电能及测地W端点均分，磁W四顶点均分。于是

$$
J_{ij}:=J(i\leftarrow j)=\frac{i}{\hbar}[h_j,h_i]
=\sum_{A,Z}(a_{jA}b_{iZ}-a_{iA}b_{jZ})D_AW_Z,\qquad
\dot h_i=\sum_jJ_{ij},\qquad J_{ij}=-J_{ji}.
\tag{2}
$$

这是完整算符恒等式，不删原势、电场或颜色部门。所有项规范不变，流也保持物理共同核。这里只建立对称微分算符及其核上矩，不未经证明宣称每份流在全空间本质自伴或已有实际读取装置。

任意实节点测试函数f，记bar f_A=Σ_i a_iA f_i、bar f_Z=Σ_i b_iZ f_i。弱形式为

$$
H[f]=\sum_i f_i h_i,\qquad
\dot H[f]=\sum_{A,Z}(\bar f_Z-\bar f_A)D_AW_Z
=\frac12\sum_{i,j}(f_i-f_j)J_{ij}.
\tag{3}
$$

节点动能—本地势、同边电能—测地势分别因分配相同而不给成对流。后一项有真实内部交换，只是在本节点分配中抵消。对e=(i,j)及磁面f，剩余为

$$
J^{\rm scalar}_{ij}=\frac12(D_j-D_i)W_e,\qquad
J^{\rm mag}_{ij}=\frac18\sum_f\sum_{e\subset\partial f}
\bigl(1_{j\in e}1_{i\in f}-1_{i\in e}1_{j\in f}\bigr)D_eW_f.
\tag{4}
$$

磁项可连面内对角顶点，最大图距离2。仅凭连续性可加散度零的面环流；选定h_i和式(2)才确定本轮流，也不自动给全部空间应力。

## 3. 曲目标及完整磁势的实际导数

沿574，F_x=M−|x|²/6、M=2，c=(M−x·y/6)/√(F_xF_y)，a(c)=arcosh(c)/√(c²−1)，且a(1)=1。半距离平方的端点导数为

$$
\partial_x\frac{d_{\mathcal K}(x,y)^2}{2}
=a(c)\left(\frac{cx}{F_x}-\frac{y}{\sqrt{F_xF_y}}\right)
=-\mathcal K_x\log_x y.
\tag{5}
$$

设y=R(g_e)φ_j、k_e=εψ_e²、v_i=w_i⁻¹K_i⁻¹P_i，则式(4)的经典符号是

$$
j^{\rm scalar}_{ij}
=\frac{k_e}{2}\left\langle v_i+\mathsf P_{y\to\phi_i}R(g_e)v_j,\log_{\phi_i}y\right\rangle_{\mathcal K_i}.
\tag{6}
$$

mathsf P沿唯一目标测地线平行移动；正号指从j流入i。原完整磁势沿569字符和正补项，不换成基本迹。正向链路上第a规范因子的电动能使

$$
\dot g_{e,a}=i\Omega_{e,a}^At_Ag_{e,a},\qquad
\Omega_{e,a}^A=2b_a\epsilon^{-1}\psi_e^{-2}P_{e,a}^A,\qquad
(D_{e,a}W_f)_{\rm cl}
=\epsilon^{-1}\psi_f^{-2}\,dV_{\rm mag}(g_f)[\delta_{e,a}g_f].
\tag{7}
$$

圆群取t=1。δg_f对四个有序因子逐项求导；逆边用δ(g⁻¹)=−g⁻¹δg g⁻¹。代码对原569函数作独立有限差分，包含非零颜色、弱及圆群holonomy，另核独立lift变化和节点规范变化。下降性与共轭不变性由原函数继承，不另选磁模型。

## 4. 同一流接到原切向来源

以下为**固定光滑经典源上的弱采样命题**。周期域上F、ψ有固定正下界，字段及测试f有足够有界导数；沿571—574取节点／中点链路。其精确Gauss动量修正保留O(ε²)密度误差。

测地边给log_φ R(g)φ_next=εD_aφ+O(ε²)，v=ψ⁻⁶K⁻¹p+O(ε²)，所以

$$
(\dot H_\epsilon[f])_{\rm scalar,cl}
=-\int d^3x\,\psi^{-4}(\partial_a f)p_A D_a\phi^A+O(\epsilon).
\tag{8}
$$

证明：每条边j=ε²ψ⁻⁴p·D_aφ+O(ε³)，弱权f_i−f_j=−ε∂_a f+O(ε²)，按ε⁻³条边求和。K与K逆在原动能及边势中抵消，未改用平直目标。

取向μν的面满足g_f=exp(iε²F_μν+O(ε³))，磁函数Hessian为K_a I；电动速度2b_a εψ⁻²E_μ+O(ε²)。面平均与两条平行边平均之差为±ε∂_ν f/2+O(ε²)。四边合起来给

$$
(\dot H_\epsilon[f])_{\rm gauge,cl}
=-\int d^3x\,\psi^{-4}\sum_a(2b_aK_a)
(\partial_\mu f)E_a^{\nu A}F^A_{a,\mu\nu}+O(\epsilon).
\tag{9}
$$

单位元处混合Hessian为零，高阶字符交叉项进入余项，非Abelian输运已含在F中。每面局部余项至多O(ε⁴)，总和O(ε)。规范约束不是此符号恒等式的必要条件；物理态连接时仍必须满足Gauss。

569既定匹配给2b_aK_a=1。令572原规范协变动量密度为mathcal D_a=p_A D_aφ^A+Σ_c E_c^{bA}F^A_c,ab，则

$$
(\dot H_\epsilon[f])_{\rm cl}
=-D_m[\psi^{-4}\nabla f]+O(\epsilon),\qquad
D_m[\xi]=\int d^3x\,\xi^a\mathcal D_a,\qquad
S^a=-\psi^{-4}\delta^{ab}\mathcal D_b.
\tag{10}
$$

S是单位lapse、零shift的向外能源流密度，dot H[f]=∫S·∇f。规范协变D在Gauss面上与普通空间变换生成元相差规范生成项；不能直接将有限图流称为空间微分同胚生成元。若单独改某规范磁Hessian而保电能及几何校准，式(9)保留2bK，不能沿用同一D来源。这是既有匹配在流层的检验，非新推出耦合值。

574波包对固定图、固定紧源管内的一阶不变符号适用。在同一来源及稳定子前提下，

$$
\langle\dot H_\epsilon[f]\rangle_{\Psi_{\hbar,\epsilon}}
=(\dot H_\epsilon[f])_{\rm cl}+O_\epsilon(\sqrt\hbar),\qquad
\lim_{\epsilon\to0}\lim_{\hbar\to0}\langle\dot H_\epsilon[f]\rangle
=-D_m[\psi^{-4}\nabla f].
\tag{11}
$$

常数依赖图、源和管，未交换极限，也未给固定ℏ的连续量子应力。578—579的能源細化障碍不被解决。本轮只核单位lapse与一个测试函数的关系，未完成任意两lapse括号或全部形变代数。

## 5. 实际记录保平均流，却不保流涨落

流可写J=−iℏ(V·∂+div_μ V/2)，V是完整实配置向量场。对586正乘法Kraus族，ΣL_r²=1给ΣL_r VL_r=0。在共同核上，

$$
\Phi^\dagger(J)=\sum_r L_rJL_r=J,\qquad
\sum_r\|JL_r\Psi\|^2-\|J\Psi\|^2
=\hbar^2\left\langle\sum_r(VL_r)^2\right\rangle_\Psi.
\tag{12}
$$

证明：J(LΨ)=LJΨ−iℏ(VL)Ψ，交叉项按完备性消失。这是核上的二阶矩／二次型恒等式，不额外假定J的谱测量已实现。对586两次细读与直接奇偶读，

$$
\langle J^2\rangle_{\rm fine}-\langle J^2\rangle_{\rm direct}
=\hbar^2\left\langle(Vs_B)^2\Delta A(s_B)\right\rangle,\qquad
\Delta A=\frac{\cos^2s_B}{8(1+\sin^2s_B/4)}.
\tag{13}
$$

对B、C的测地边，Vs_B=k_e(log_B Rg φ_C)^s/(2w_B)。磁流只对链路微分，对本次singlet读取的直接贡献为零；磁动力学仍在H中。系数在允许配置上不恒零；586同类独立径向紧包及Haar链路在非零开集有正权，所以存在严格Gauss、有限矩来源使差严格为正。两种读后均值相同，此差也为方差之差。

路径积分对接所保留的非对角信息由此进一步接到**流涨落**。相同平均报告和平均流不足以签收同一来源；本式也不是完整协变应力噪声核。

## 6. 有源连续性与共同过程

原H等待、中途真实非选择仪器下，在具有所需有限矩和可微性的状态范围内，读取时刻t_k的源注能s_i,k给

$$
\frac{d}{dt}\langle h_i\rangle
=\sum_j\langle J_{ij}\rangle+\sum_k s_{i,k}\delta(t-t_k),\qquad
s_{i,k}=\langle\Phi_k^\dagger(h_i)-h_i\rangle_{t_k^-},\qquad
s_{i,k}=\delta_{iB}\Delta E_k
\quad\text{对本轮B配置读取}.
\tag{14}
$$

这是非选择平均账，选择分支另有条件化变化。一个形式对易子不自动解决全域算符域传播；共同核上使用瞬时导数及二次型，有限区间推广须维持上述矩条件。理想读取的δ项不是自主装置功率，装置、控制及记录总成本仍开放。

CTP实际分支核可保留这些插入及涨落；仅替成同POVM会删改注能与二阶量。CTP不自动决定h_i、局部lapse耦合或空间维数。

## 7. 可复算结果及范围

五组检查：完整分配入口复算、原曲目标导数、全部商群磁势微分、同一严格Gauss来源弱流、原仪器流二阶差。

- 4³图1408项：反对称误差0，连续性残差2.66×10⁻¹⁵，最大支撑距离2。
- 测地导数最细中心差分误差1.21×10⁻⁸，按二阶缩小；规范协变差3.47×10⁻¹⁸。
- 非零全规范面势9.1530731，四边实际导数最细误差6.25×10⁻⁶；独立lift／节点规范变化误差至多3.55×10⁻¹⁴。
- 原574修正Gauss源，ψ=1.1+.05cos x+.02sin y，测试sin x、cos(x+y)、sin x cos y。N=8、12、20、32的总弱流误差0.10321、0.04777、0.01759、0.006924。观察到约二阶，解析仅用保守O(ε)界。Gauss密度残差最大2.30×10⁻¹⁶。源的弱／圆群非零、颜色为零；完整颜色导数由前组另测。
- 实际曲边流s主系数的局部微分诊断，二阶差系数0.00018472076231，最细展开误差3.40×10⁻⁹。诊断固定邻居，不冒充全图量子传播；Gauss二阶矩结论由式(12)—(13)承担。

没有重新求引力初值，指定光滑ψ不是新Einstein解。一般群微分、共同核与Taylor展开承担证明；有限格收敛只核所列来源，不能替代统一量子误差界。

## 8. 条件合并与下一步

[587条件账](587/unified_physics_condition_ledger_587.md)接通原物质、能源流、连续切向来源及记录涨落，减少另立流模型的需要；保留585局部耦合选择、578—579量子细化障碍及自主仪器缺口。

下一项[588共同局部几何源族](588/drafts/STATUS.md)：从同一h_i检验任意两lapse的量子括号及经典来源，核能否保留真实记录插入和涨落并接入共同几何约束。只定义一个对易子并改名“shift”不算完成，须与独立切向作用及其组合比较。不继续优化单读口，不改统一目标。
