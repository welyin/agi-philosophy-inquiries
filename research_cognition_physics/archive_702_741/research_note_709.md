# 第709轮：保原量子过程的荷粗化、跨区相位与反常接口

日期：2026-10-03。接[708](research_note_708.md)与[709入口](709/drafts/conditional_reference_entry.md)。[代码](709/joint_charge_quantum_coarse.py)、[结果](709/joint_charge_quantum_coarse_results.json)、[核验](709/research_round_709_checks.json)、[条件账](709/unified_physics_condition_ledger_709.md)。

## 1. 从资料不足转向可行的量子描述

708及709入口说明，经典配置统计不能承接所有后续量子操作。此次不再追加矩，而在**原完整有限图H及原Gauss空间**中构造一份真正较小、仍非对易的可观测代数，检验同一参考、演化、读口、质量和资源能否一起保留。

得到三个共同结论：

1. 原有限图H严格保持总夸克占据数。只要求保留对此荷中性的任务时，荷扇区之间的相干可以从描述中去掉；扇区内部的量子态、原热参考、全部允许多时历史和原能量均精确保留。
2. **这一全局构造不能分别在各区域独立执行。** 原跨区夸克跳跃需要相同总荷下、不同局部荷之间的相干。明确的正常Gauss态给非零输运；逐区去相干将其删掉。
3. **不能把该连续U(1)原样提升为整个统一模型的基本认知对称。** 629已经给出同一重定相方向在连续手征测度上的非零指数Jacobian。参考／粗化选择与反常部门因此必须共同验收。

|层次|准确地位|
|---|---|
|认知动机|某份信息只有对所允许的未来过程确实不可用时才可从整体描述中删除|
|继承输入|598完整有限图H、全部CAR、给定Y、规范协变跳跃及原Gauss；629另列的连续手征分支|
|新增任务约定|先要求保所有总夸克荷中性的操作；这不是所有可想象的Gauss操作|
|成熟工具|守恒荷分解、群平均、条件期望、夹断相对熵身份|
|项目增量|原全质量／Gauss上实际严格粗化，正常物理态的真信息损失，原跨边流的区域障碍及与629的共同范围条件|
|未完成|空间自由度的粗化、最小充分态、原手征过程正性／连续映射、实际装置及GR生成|

此构造保留了很多微观变量，不能当成宏观时空涌现的完成。它提供的是可行的非对易内部描述及其明确组合限制。

## 2. 去重、文献与原对象映射

598已经有原质量和CAR；629有重子方向对质量及连续测度的不同作用；642有正常Gauss粒子—空穴态及删CAR失败；705有区域边界支持；707—708有实际标量分块和资源连接。它们均直接复用。

[Bartlett—Rudolph—Spekkens，§II.3、III.3及IV.1](https://arxiv.org/html/quant-ph/0610030)给群平均、表示重数空间和共同／各自相位参考的成熟框架。本轮把其中群作用具体替换成**原有限CAR的总夸克占据数**，逐项核原H、Gauss和跨边项。群平均公式及全局／局部区别本身不当新理论发现，也不从参考限制推导该群是宇宙基本对称。

对接连续部门时直接使用629冻结的指数矩阵。标准模型中重子数破坏过程具有具体动力学研究，例如[D’Onofrio等的sphaleron速率计算](https://arxiv.org/abs/1404.3565)；这是理论模拟而非本项目的实验验证。本轮不重算该速率。

## 3. 原有限模型的荷与严格分解

沿598物理右手约定，每代每节点32个模式。Q_L、u_R、d_R共24个，记其一粒子投影为q，定义

$$
q=\operatorname{diag}(I_{24},0_8),\qquad
N_q=\sum_{v,a,b}c_{v,a}^\dagger q_{ab}c_{v,b},\qquad
V_\alpha=e^{i\alpha N_q},\quad 0\le N_q\le24n_g|V|.
\tag{1}
$$

这里是原有限CAR占据数，不直接等同于连续重整化的重子流。经614粒子—空穴字典后，uᶜ、dᶜ带相反权并出现无关共轭作用的常数；629的左手重定相方向由此一致对应。

原质量与链路表示满足

$$
[q,h(\phi)]=0,\qquad q\Delta(\phi)+\Delta(\phi)q^{\mathsf T}=0,
\qquad [q,R(g)]=0.
\tag{2}
$$

Yukawa只在夸克内部或轻子内部变换；Majorana仅作用于ν_R，因此满足第二式。原跳跃按同一表示及spin／代指标输送，不把夸克变为轻子，所以总N_q守恒。一般**总费米数**却不守恒，不能用它替代式(1)。所有纯玻色项也与N_q对易。

N_q是有限谱有界算符，与全部节点Gauss作用对易，保持原闭能量形式及核；各谱投影于是约化自伴H。设P_k为它在物理Hilbert空间中的非零谱投影，定义

$$
\mathfrak A_q=\{X\in\mathfrak M_{\rm phys}:[X,N_q]=0\},\qquad
\mathcal E(X)=\sum_kP_kXP_k
=\frac1{2\pi}\int_0^{2\pi}V_\alpha X V_\alpha^\dagger\,d\alpha.
\tag{3}
$$

可取物理偶可观测代数为M_phys；以下见证都在偶部门。E是正常、幺正保持、完全正、幂等的条件期望，双模性质直接由P_k给出。其迹类对偶E_*有同一夹断公式，正常CPTP且保Gauss支持。全部空间／规范／CAR重数资料仍在各荷块内。

## 4. 原参考、演化及真实历史共同保留

对原U_t=e⁻ⁱᵗᴴ/ℏ、σ_β=e⁻ᵝᴴ/Z，有

$$
[H,P_k]=0,\qquad
\mathcal E_*\big(U_t\rho U_t^\dagger\big)=U_t\mathcal E_*(\rho)U_t^\dagger,
\qquad \mathcal E_*(\sigma_\beta)=\sigma_\beta.
\tag{4}
$$

等式针对所有实t及正常态，未换Hamiltonian、丢掉相互作用或在每次操作后重设热态。粗描述上的动力学仍由原块内H产生；这一路不需要为所列任务另配记忆核，因为那些内部资料被完整保留了。

取任意有限历史，每个Kraus算符与N_q对易，中间等待仍为原U_t。整个有序历史算符记K_h，则

$$
[K_h,N_q]=0,\qquad
\mathcal E_*(K_h\rho K_h^\dagger)=K_h\mathcal E_*(\rho)K_h^\dagger,
\qquad
\operatorname{Tr}(K_h\rho K_h^\dagger)
=\operatorname{Tr}(K_h\mathcal E_*(\rho)K_h^\dagger).
\tag{5}
$$

因此保真实概率、条件态对A_q的限制和整个经典记录／粗后态，而非只保单次期望。认知任务约定没有宣称必须限制为中性；加入其它操作须重新验收。附加被动中性参考时，E_*⊗id也给同一身份；若参考可交换该荷，须纳入总荷及共同操作重新分析，不能假装免费外部相位标准。

全部原节点质量、跳跃、708的P₂／预算和709入口相位V=e^(iαsinZ)都在这份保留代数内。尤其入口的相位准备不会再被配置统计压掉：

$$
\mathcal E_*(V\sigma_\beta V^\dagger)=V\sigma_\beta V^\dagger.
\tag{6}
$$

对所有存在的原能量矩，以及固定原参数族内的荷中性来源，有

$$
\operatorname{Tr}\big(\mathcal E_*(\rho)H^m\big)=\operatorname{Tr}(\rho H^m),
\qquad
\operatorname{Tr}\big(\mathcal E_*(\rho)\partial_\gamma H_\gamma\big)
=\operatorname{Tr}(\rho\partial_\gamma H_\gamma).
\tag{7}
$$

先在有限谱／原核上验证，再以已有形式与热矩界延拓。E与γ无关，所以对**已经存在**的迹类响应导数也可交换；不借此补签707／708新读口尚未证明的所有动力二阶来源域。

## 5. 内部资源账与一份真的非对易真子代数

有限能量态在原Gibbs约束下有有限熵。因σ_β与P_k对易，夹断的相对熵分解为

$$
D(\rho\Vert\sigma_\beta)
=D(\rho\Vert\mathcal E_*\rho)+D(\mathcal E_*\rho\Vert\sigma_\beta),\qquad
D(\rho\Vert\mathcal E_*\rho)=S(\mathcal E_*\rho)-S(\rho)\le\log r,
\tag{8}
$$

r是非零荷扇区数，≤24n_g|V|+1。证明用logσ及log(E_*ρ)的荷块对角性，并用ρ≤rE_*ρ控制最后一项。该界随系统规模变动，不声称资源统一界。

E_*在此是描述的粗化，不是已实现的物理遗忘装置。原能量保存不能被解释成记录擦除免费；若实施随机相位，控制器／相位记录仍须在整体内部处理。保留代数所列真实仪器的原能量代价则由式(7)共同保留。

证明它确实丢了一份原物理信息：在某节点取两个不同spin标签，原右手模式构成

$$
B_\sigma^\dagger=\sum_{a,b,c=1}^3\epsilon_{abc}
u_{a\sigma}^\dagger d_{b\sigma}^\dagger d_{c\sigma}^\dagger,
\qquad |b\rangle=\mathcal N B_\uparrow^\dagger B_\downarrow^\dagger|0\rangle.
\tag{9}
$$

每个udd的超荷为4−2−2=0，颜色用ε收缩，弱群作用平凡，故|b>为严格Gauss态；两份spin模式给非零六粒子偶态。乘以574的共同规范不变紧支撑光滑玻色波包，得到原完整模型中的正常有限能量态。令|v>为同波包乘CAR真空。

$$
|\psi\rangle=\frac{|v\rangle+|b\rangle}{\sqrt2},\qquad
\mathcal E_*(|\psi\rangle\langle\psi|)
=\frac{|v\rangle\langle v|+|b\rangle\langle b|}{2},\qquad
\frac12\big\||\psi\rangle\langle\psi|-\mathcal E_*(|\psi\rangle\langle\psi|)\big\|_1=\frac12.
\tag{10}
$$

两态能被荷不守恒但仍Gauss／偶的探针区分，故没有把荷守恒偷偷升为全部物理操作的超选择公理。另一方面，ν_↑†ν_↓†及n_ν↑均保N_q、保Gauss、为偶，且交换子在真空上非零；真空与中性ν对之间的量子相干保留。故A_q确是真非对易子代数，不是只剩经典荷标签。

## 6. 同一构造在区域拼接处有额外要求

对区域A定义N_A为其中夸克占据数。跨边e取原非零夸克跳跃通道，记由A端v向另一端w的规范不变偶算符为T_e：

$$
T_e=\sum_{a,b}c_{w,a}^\dagger R(g_e)_{ab}c_{v,b},\qquad
[N_A,T_e]=-T_e,\quad [N_{\bar A},T_e]=T_e,\quad [N_q,T_e]=0.
\tag{11}
$$

因此局部荷夹断E_A、E_Ā会删T_e，而全局E保T_e。原H含t_eT_e+t_e* T_e†；选实t_e的既有通道作校准时，切口流由i(T_e−T_e†)给出。若spin矩阵有其它分量，逐通道使用同样身份，不假设所有图边相同。

这不是只在非物理单夸克态中的现象。取v全部32模式填满、w真空的CAR态Ω，其它节点选合法真空，乘原共同正常玻色波包。原完整表示det R(g)=1，故全满海也是局部Gauss singlet。对一个固定spin的u_R颜色通道，d=3，令

$$
|\chi\rangle=\frac{T_e|\Omega\rangle}{\sqrt d},\qquad
\|T_e\Omega\|^2=d,\quad
(N_A,N_{\bar A})\Omega=(24,0)\Omega,\quad
(N_A,N_{\bar A})\chi=(23,1)\chi.
\tag{12}
$$

T_e保Gauss，波包乘光滑链路矩阵仍属正常有限能量核；这复用642粒子／空穴的实际方法。Ω、χ都具有总N_q=24及同一偶宇称，且正交。其明确相位态满足

$$
|\psi_\pm\rangle=\frac{\Omega\pm i\chi}{\sqrt2},\qquad
I_e=i(T_e-T_e^\dagger),\qquad
\langle I_e\rangle_{\psi_\pm}=\pm\sqrt d.
\tag{13}
$$

逐区夹断令两态都变为(ΩΩ†+χχ†)/2，电流变成0；全局夹断则保持两态本身：

$$
\mathcal E_*\rho_\pm=\rho_\pm,\qquad
(\mathcal E_A)_*(\mathcal E_{\bar A})_*\rho_\pm
=\frac{\Omega\Omega^\dagger+\chi\chi^\dagger}{2}.
\tag{14}
$$

这也可从实际原演化验收：E_A保持初态ΩΩ†，但原H的该跳跃在一阶产生Ω—χ非对角项，E_A会删掉它，所以E_A不与原演化交换。数值所用两向量不是H的不变子空间，未把2×2矩阵当作完整动力学。

在已有617／637边界荷匹配下，区域总荷k的载体还须保留

$$
\mathcal H_k\subseteq\bigoplus_{a+b=k}\mathcal H_{A,a}\mathbin{\widehat\otimes}\mathcal H_{\bar A,b}
\quad\text{中的不同 }(a,b)\text{ 项之间的相干}.
\tag{15}
$$

帽号提醒仍有原规范边界配对和CAR分次，未假定物理区域空间普通独立因子化。全局不可用的相位与区域之间可用的相对相位必须分别处理。

## 7. 与旧连续反常条件共同验收

629的全左手模块按(Q,uᶜ,dᶜ,L,eᶜ,νᶜ)排列，原质量相位矩阵C和指数矩阵A已被冻结。对重子方向b=(1,−1,−1,0,0,0)，有

$$
Cb=0,\qquad -A^{\mathsf T}b=(0,0,3,0).
\tag{16}
$$

它说明质量保持不等于指定连续手征测度保持。在n_g代、允许629全部整数弱拓扑数的同一参数描述中，共同相位α的Jacobian含exp(i3n_gαc₂(E₂))；符号沿629约定。此处c₂(E₂)为拓扑整数，不是本轮有限图上的自由新噪声。

若要求固定其它参数时，该重定相对**每个**所列拓扑部门的测度都不变，至少必须满足

$$
e^{i3n_g\alpha}=1.
\tag{17}
$$

所以有限图的整个连续U(1)不能仅凭质量对称就在这条连续分支中被签为相同固定参数的精确对称。这里未重新证明指数定理，也未声称所有有限图模型不可能有手征极限；空间/时间物理映射尚未建立，699的特定反射正性问题也没有被修复。

一个可检验的失败判据是：若某共同Hilbert表示下的有界酉族U_a(t)、V_a(α)强收敛，且每个尺度保持可比的同一个荷动作，那么

$$
[U_a(t),V_a(\alpha)]=0,\quad U_a(t)\xrightarrow{\rm s}U(t),\quad V_a(\alpha)\xrightarrow{\rm s}V(\alpha)
\quad\Longrightarrow\quad[U(t),V(\alpha)]=0.
\tag{18}
$$

乘积的强收敛因酉族一致有界而成立。要与具有不同重子Ward身份的目标连接，至少须明确哪个条件改变：物理流的字典、UV／拓扑资料、对称实现、允许任务或极限的强度。不能一面坚持上述全U(1)合同原样成立，一面直接借629的反常破坏而不交代连接。此为条件性接口限制，不是全统一目标的反证。

## 8. 复算、联合限制及下一项

四组复算通过：

- 原32模式质量、Majorana与规范链路的荷身份及完整海det R=1；总费米U(1)由非零Majorana排除。
- 原六夸克正常态的9个非零CAR分量，在所有81个允许规范变换输出上核不变；荷夹断迹距离1/2，中性ν对仍有非对易相干。
- 原64模式位标上的满海—空穴—Wilson链路态，切口电流为±√3，局部夹断后为0；仅检查实际矩阵元，不模拟二向量封闭H。
- 直接复用629的整数矩阵核质量方向与拓扑相位，未重新登记反常计算或给sphaleron率。

本轮把C01／C03／C04／C17／C19／C21接到同一全局量子描述，并给C02／C20的区域相干要求；C16／C18同时限制这份描述能否沿用到连续手征目标。读口、参考和动力学由同一条件期望承接，省去为它们各自拟定恢复规则的自由度，但未选出原群、Y、代数目或宇宙几何。

下一项优先联合核查：在629允许的拓扑部门与实际区域合并下，哪些较小对称或内部参考描述还能共同保留？先检查原全局群中心、CAR宇称及残余离散动作，区分非平凡物理作用与本来就是规范的动作；同时确认它是否保跨区荷相干。不能把找出一个有限图可行粗化直接当作全尺度认知结构。

旧空间382—386、425、522—524按原条件复用；384删除的Lipschitz不恢复，386／425是替代桥，实际端点映射仍须接到当前模型。统一目标未完成。
