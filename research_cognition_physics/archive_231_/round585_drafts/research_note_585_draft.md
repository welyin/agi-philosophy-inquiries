# 第585轮：同一量子演化为何仍未确定局部几何来源

日期：2026-10-01。接[584](research_note_584.md)和[合并清单](joint_condition_compression_update_584.md)。本轮五组复算；主代理核验，未新增独立代理审查。正式核验见[检查记录](research_round_585_checks.json)。

## 1. 结果、输入及去重

本轮定位并检验了一个具体连接缺口：**574的单位lapse有限图Hamiltonian，即使连同584的全部既有记录与空间权响应一起保留，也没有唯一确定局部Hamiltonian约束所需的能源分布。**

构造只改变“如何把局部时间推进率耦合到原Hamiltonian”的扩展。它不新增物质、不删除原势、不改Gauss和单位时间演化。对于任意正lapse，候选仍有正闭二次型。原严格Gauss态与原相位准备给出显式的不同局部能源来源。

同时得到限定的正向筛选：若沿用已有连续理论的lapse几何意义、切向生成元和法向形变关系，这种只平滑物质lapse的改动不能任意加入。这个筛选是原有匹配条件的落实，不是已经从认知原则新增推出了广义相对论。

|层次|本轮地位|
|---|---|
|认知动机|同一来源、记录与几何必须共同定义，不能只比较总能源|
|继承|574有限图、正F曲目标、原规范物质与LB排序；577相位准备；584表示字典|
|候选输入|有限图局部正能源分配、正lapse扩展，以及比较时固定的lapse校准|
|解析增量|同一单位lapse理论的正局部扩展非唯一；原Gauss源见证；光滑采样与量子资源的不同极限；受限法向括号筛选|
|未完成|量子引力约束、shift完成、协变应力、完整量子连续极限及物理结构选择|

[351](research_note_351.md)—[352](research_note_352.md)已经研究任意lapse的条件性形变关系，[366](research_note_366.md)已证明这种关系本身不等于动态引力；本轮不重复宣称这些是新结果。[350](research_note_350.md)的状态依赖均值反馈反例也直接复用。本轮新增对象是574原完整有限图及其正规Gauss来源。

## 2. 为什么需要局部时间推进率

574明确取单位lapse、零shift。已有同一经典连续作用可以给出完整经典来源，但它没有自动选定有限截止下的量子局部耦合。空间共形权的偏导也不同于时间方向的变分。

沿用标准正则约定，连续物质生成元和来源分别为

$$
H_m[N,\beta;\gamma]
=\int d^3x\,[N\,\mathcal H_m+\beta^i\mathcal D_i],\qquad
\mathcal H_m=\sqrt\gamma\,\rho,\qquad
\mathcal H_m(x)=\frac{\delta H_m}{\delta N(x)}.
\tag{1}
$$

这里N是lapse，控制相邻空间切片之间的法向时间推进；beta是切向推进。rho不是一个单独的总能量数。若只知道N=1、beta=0上的函数值及其空间权偏导，不能直接作未定义方向的泛函微分。

已核对[Gourgoulhon，§4.1.2及§10.2](https://arxiv.org/pdf/gr-qc/0703035)：法向能源、动量与空间应力分别进入初始约束和约束传播。[Teitelboim，式(6)、(28)及§4](https://www.fis.uc.cl/~mbanados/Cursos/TopicosRelatividadAvanzada/Teitelboim.pdf)提供不同切片推进的几何相容条件。两者都以几何对象为前提；本轮仅采用其接口，不借此声称生成时空。

## 3. 原完整图上的两种正局部扩展

取已有周期三维图，每方向至少五个节点以避免邻居重合。将574原节点动能与现场势记T_i、V_i；原边梯度及电能记L_e；原面磁势记B_f。这些是非负二次型。按端点／面顶点均分，定义

$$
h_i=T_i+V_i+\frac12\sum_{e\ni i}L_e
          +\frac14\sum_{f\ni i}B_f,\qquad
\sum_i h_i=H_{\epsilon,\psi},\qquad h_i\ge0.
\tag{2}
$$

没有删除颜色或弱规范动能，没有改成平直目标。每项都严格内部规范不变；原全图正形式域记为D。式(2)是一种明确的局部分配，并非唯一性定理。

对0≤alpha≤1，定义对称的邻居平均

$$
(P_\alpha f)_i=(1-\alpha)f_i+\frac\alpha6\sum_{j\sim i}f_j,
\qquad P_\alpha\mathbf1=\mathbf1,
\qquad P_\alpha^T=P_\alpha,
\qquad h_i^{(\alpha)}=(P_\alpha h)_i.
\tag{3}
$$

本轮数值取alpha=0.4，不扫描最优参数。平均增加了有限的支撑范围，仍与总图大小无关；没有声称保持原来每项的精确最小支撑半径。候选lapse扩展为

$$
H_\alpha[N;\psi]=\sum_iN_i h_i^{(\alpha)}[\psi]
                 =\sum_i(P_\alpha N)_i h_i[\psi].
\tag{4}
$$

对任意固定严格正N，令N_-、N_+为有限图上的最小与最大值。非负形式给

$$
N_- H_{\epsilon,\psi}\le H_\alpha[N;\psi]
\le N_+H_{\epsilon,\psi}\quad\text{（二次型意义）}.
\tag{5}
$$

因此加入范数后的形式范数与原闭形式范数等价，D上得到正闭形式及其自伴算符，规范群仍约化该算符。这里只主张共同**形式域**，不擅自断言所有lapse下闭算符域相同。固定N的演化是与输入态无关的酉演化，不存在350那种将输入均值代入Hamiltonian的步骤。

## 4. 保留了什么，以及哪里不同

对每个固定图、全部正空间权psi，严格有

$$
H_\alpha[\mathbf1;\psi]=H_{\epsilon,\psi},\qquad
\partial_{\psi_k}H_\alpha[\mathbf1;\psi]
=\partial_{\psi_k}H_{\epsilon,\psi},\qquad
\frac{\partial H_\alpha[N;\psi]}{\partial N_i}=h_i^{(\alpha)}.
\tag{6}
$$

前两式意味着原单位lapse下的谱、全部准备和记录统计、能源以及已定义的空间权响应均相同，适用于任意原允许态。最后一式说明，同一个态的局部法向能源一般不同。各自除以原节点体积才是候选密度；不是改了体积定义。

在固定几何校准N下，局部差为

$$
\Delta\langle H[N]\rangle
=\sum_i[(P_\alpha N)_i-N_i]\langle h_i\rangle.
\tag{7}
$$

这不是对“原完整协变作用能决定其应力”的反驳。它证明有限图单位lapse资料没有保留全部协变耦合信息。若同时重新定义N、引力部门和钟的校准，应当作为另一套完整字典检验；式(7)比较的是同一N，而不是排除所有重定义。

## 5. 原Gauss来源与原准备的严格见证

每节点采用584的实紧支撑径向态，链路取Haar常函数。曲测度波函数由该轮的半密度字典取得。全图态正规、严格Gauss，处于所需有限次H幂域。对A节点施加577原准备

$$
U_\theta=\exp(i\theta h_A^2/(2\hbar)),\qquad h_A=|X_A|.
\tag{8}
$$

该准备不改变乘法势期望，也不改变链路动能。实幅度使节点动能的线性相位交叉项为零；唯一能源增量来自A节点：

$$
\delta E=\frac{\theta^2}{2w_A}
\left\langle F_A\left(h_A^2-\frac{h_A^4}{6M}\right)\right\rangle>0,
\qquad w_A=\epsilon^3\psi_A^6.
\tag{9}
$$

支持h∈[0.42,0.62]、s∈[0.25,0.55]处F>0，且h²<6M，故正性由解析界保证。式(2)将增量全部记在A；式(3)则给

$$
\delta\langle h_A^{(\alpha)}\rangle=(1-\alpha)\delta E,
\qquad
\delta\langle h_j^{(\alpha)}\rangle=\frac\alpha6\delta E\quad(j\sim A),
\qquad\sum_i\delta\langle h_i^{(\alpha)}\rangle=\delta E.
\tag{10}
$$

数值theta=0.3、hbar=0.7、w=0.512，得delta E≈0.0460104368949481；alpha=0.4后，A为0.0276062621369689，六个邻居各为0.00306736245966321。正增量可由直接积分带相位动能独立复算，不是任意指定的经典正数组。

在5³图、N_i=1+0.2 cos(2πi_x/5)下，两种处方对此次准备给出的lapse能源增量差约−0.000847798132688。总能源和原单位时间记录完全相同。这里比较的是局部来源响应，不宣称已造出测量局部lapse的自主装置。

## 6. 相同领先经典极限不能替代量子尺度控制

对同一光滑N的网格采样，中心差分的精确余项估计给

$$
P_\alpha N=N+\frac{\alpha\epsilon^2}{6}\Delta N+O(\epsilon^4),\qquad
|\Delta\langle H[N]\rangle|
\le\frac{\alpha\epsilon^2}{6}
\left(\sum_{a=1}^3\|\partial_a^2N\|_\infty\right)\langle H\rangle.
\tag{11}
$$

不等式在正有限能量态成立，使用所有h_i非负。固定光滑经典能源密度的有界总能采样给O(epsilon²)差；所以574的领先经典极限不选择alpha。该界也明确指出，必须同时控制能源。

例如固定psi=1、theta和同一节点包形状，记C=delta E·epsilon³>0。在周期长度2π、N=1+nu cos x、准备节点x=0的细化族上，式(10)精确给

$$
\delta E=\frac C{\epsilon^3},\qquad
\Delta\delta\langle H[N]\rangle
=-\frac{\alpha\nu C}{3\epsilon^3}(1-\cos\epsilon)
\sim-\frac{\alpha\nu C}{6\epsilon}.
\tag{12}
$$

这是一族每张图上都正规、Gauss且有限能的态，但准备成本没有网格一致上界。它不是有界能源量子连续极限的反例。它说明“平滑误差是二阶”单独不够；必须把578—579的能源限制带入同一比较。代码只积分单节点包并利用全图解析恒等式，不模拟指数维全图波函数。

## 7. 既有法向形变条件如何筛选该扩展

此节是独立的**经典连续接口诊断**。沿原singlet测地子部门取X=0，令s=√(6M)tanh(chi/√6)，目标动能成为规范化的一标量；保留原U(s)，不是改成自由粒子。以平直周期切片、单位横截面积写

$$
H[N]=\int dx\,N\left[\frac{\pi^2}2+\frac{(\chi')^2}2+U(s(\chi))\right],
\qquad D[v]=\int dx\,v\pi\chi',\qquad
\{H[N],H[M]\}=D[NM'-MN'].
\tag{13}
$$

最后的成熟关系由352继承；直接变分中势项及无lapse导数项抵消。设P是保持常量、与导数交换的固定对称平移平均，定义H_P[N]=H[PN]。固定原N的几何意义与原D，则

$$
\{H_P[N],H_P[1]\}-D[-N']
=-D[(PN-N)'].
\tag{14}
$$

若对所有离壳chi、pi及N都要求原关系，右侧只能恒零，即(PN−N)'=0。对本轮邻居平均的非零Fourier模，乘子为1−alpha(1−cos a)/3；当平均半径a不等于周期整数倍时，alpha>0给违反见证。代码取chi=0.3+0.15 sin x、pi=0.2 sin 2x、N=1+0.2 cos x，两lapse均正；未改动原势。a=0.8时缺陷约−0.000381129589519。

这只排除“保持原切向生成元与lapse校准，却单独平均物质法向生成元，并要求原精确连续关系”的接法。它不是完整引力—物质交叉括号计算，也不证明alpha=0的有限图具有精确量子形变代数。有限截止若只要求有效近似，仍需误差匹配；共同重定义几何或钟属于另一个分支。366关于HDA不足以生成引力的结论完整保留。

## 8. 五组可复算证据

[代码](joint_local_lapse_source.py)默认重算并比对完整[结果](joint_local_lapse_source_results.json)，复用既有Python与NumPy。

|检查|结果与范围|
|---|---|
|原完整图能源分配|4³、6³原572—574经典来源；逐部门和总和误差≤1.43×10⁻¹⁴，空间权方向偏导差≤7.11×10⁻¹¹；局部lapse来源不同。全量子结论由正形式恒等式承担|
|原Gauss包及相位准备|48、80、128点求积；直接动能增量与式(9)最细差2.09×10⁻¹⁷；六邻居分配和非恒定lapse差按式(10)独立汇总|
|光滑经典采样|8³至32³；lapse差从−0.387479降到−0.0254198，差／epsilon²趋−0.661467；仅诊断声明的固定光滑采样|
|量子来源的资源限制|同包相位成本随细化从0.0486247增至3.11198；lapse增量差乘epsilon趋−0.000314098。无有界资源连续声明|
|法向括号|保留原singlet势的泛函导数直接求括号；未平均值0.00942477796077等于原D；非零平均缺陷与式(14)差小于10⁻¹³。非量子图约束检查|

初版四组结果和代码保存在round585_drafts，新增资源检查后才冻结正式五组，旧四组不另计轮次。全部检查由主代理运行与核验；没有声称取得新的独立代理审查。

## 9. 合并与下一步

[585条件账](unified_physics_condition_ledger_585.md)把缺口收束为同一个函数族：需要共同的H[N,beta;gamma]及对应状态和观察量，而不只H[1,0;gamma]。能源、流、应力与约束传播应从它共同取得。这有望减少分开指定来源的自由度；本轮尚未构造完整函数族。

下一项在原物质和有限截止内推进局部能源流：从已经选定的局部正形式与同一Hamiltonian推导离散连续性关系，核所得流能否与切向几何生成元匹配。复用351—352的成熟条件，明确有限图误差及量子范围，不继续扫描alpha，也不重做“形变代数是否足以生成GR”的旧问题。
