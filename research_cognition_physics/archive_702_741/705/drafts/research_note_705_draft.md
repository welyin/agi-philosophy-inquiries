# 第705轮：区域压缩、边界电荷与Gauss支持的共同限制

日期：2026-10-02。接[704](../../research_note_704.md)与[705入口](regional_compression_entry.md)。回查[360](../../../archive_342_369/research_note_360.md)—[365](../../../archive_342_369/research_note_365.md)、[617](../../../archive_585_628/research_note_617.md)、[625](../../../archive_585_628/research_note_625.md)、[636](../../../archive_629_652/research_note_636.md)—[638](../../../archive_629_652/research_note_638.md)。[代码](../joint_regional_charge_compression.py)、[结果](../joint_regional_charge_compression_results.json)、[核验](../research_round_705_checks.json)、[条件账](../unified_physics_condition_ledger_705.md)。三组检查、二十式；主代理审查，无新增独立代理审查。

## 1. 新结论与去重

704的共同热准备、真实历史及来源极限使用全局谱压缩。若改为“各区域独自压缩成有限维，远处不变，原Gauss匹配精确保持”，这些要求是否仍相容？

**在本原固定图的非平凡切口与完整边界支持上，不相容。** 原商群有无限多个实际出现的边界表示。有限输出丢掉其中一些，而严格区域通道不能改变另一侧所携带的匹配表示。本轮给出原模型上的误差权衡：**Gauss支持损失＋远端状态改变＋有限输出泄漏，至少等于原态落在遗漏边界扇区的概率。**

同时给出正向替代：保留所有边界表示及其载体，只压缩每个表示内的重数空间，可构造严格局部、保持原Gauss、保持任意远端／参考边缘、并在原参考能源控制下趋向恒等的CPTP族。它每个扇区有限，但整体仍无限维；不把这一代价隐藏起来。

|层次|地位|
|---|---|
|认知动机|区域、压缩、共同记录与约束必须使用同一状态和访问范围|
|继承模型|原固定图、Z₆商群、全部CAR与标量、617切边、603热态、637比较能源|
|受检验的额外要求|严格单侧操作、另一侧完全不变、固定有限维输出、原Gauss精确支持|
|解析增量|原无限边界支持下的容量／支持／远端误差权衡；保边界的局部重数压缩族及误差|
|数值范围|原字符环态、原非Abel表示与明确的有限重数诊断；不求完整H热谱|
|保留缺口|本局部族的实际动态／来源匹配、空间连续、手征统一及量子引力|

**归属补正：** 入口式(1)的压缩交换子恒等式已在363式(4)证明。[补正记录](attribution_correction_363.md)明确其旧归属，冻结入口不改。入口热尾是旧身份的推论，未计轮次；本轮不以它作为新结果。360—365的共享中心和有限Z₂桥、617精确切边、636表示熵、637—638区域能源及共同粗化均直接复用。

成熟输入为[Donnelly，§II、式(18)—(20)](https://arxiv.org/html/1109.0036)的边界表示分解，以及[Van Acoleyen等，式(1)—(5)](https://arxiv.org/html/1511.04369)的局部物理操作／扇区结构。它们支持下节的表示语言，不把本轮的原热支持、截断误差或完整物理统一归为文献已证明。

## 2. 固定原切口，而非任意改变区域定义

先在区域内部施加原顶点Gauss，保留617边界扩展。边界群为原G的有限乘积，λ是其允许表示元组；方向由另一侧取对偶。忽略没有配对对象的空块，有

$$
\mathcal H_A=\bigoplus_\lambda M_{A,\lambda}\otimes V_\lambda,\qquad
\mathcal H_B=\bigoplus_\lambda M_{B,\lambda}\otimes V_\lambda^*,\qquad
\mathcal H_{\rm phys}^{\rm cut}
=\bigoplus_\lambda M_{A,\lambda}\otimes M_{B,\lambda}\otimes\mathbb C\Omega_\lambda ,
\quad
\Omega_\lambda=\frac1{\sqrt{D_\lambda}}\sum_{i=1}^{D_\lambda}|i,i\rangle .
\tag{1}
$$

这是同一原物理空间J的像。重数M包含原区域内的标量、链路、CAR等自由度，V是边界匹配载体；不是把边界当作新增独立物种。保持原CAR分级，本轮所有通道取偶操作。

设正常CPTP通道只作用A，即Φ_A⊗id_B，要求对全部匹配物理输入仍有物理支持。对每个纯匹配输入，输出是正Kraus项之和；总输出没有码外支持，迫使每个Kraus项都没有码外分量。将任意重数向量与Ωλ代入，再用不可约表示的Schur引理，得到配对部门上的形式

$$
K_\alpha=\bigoplus_\lambda k_{\alpha,\lambda}\otimes I_{V_\lambda},
\qquad
\sum_\alpha k_{\alpha,\lambda}^\dagger k_{\alpha,\lambda}=I_{M_{A,\lambda}}.
\tag{2}
$$

反之，式(2)直接使每个分支保持匹配。因此局部物理操作可以改变区内重数，却不能把边界λ改成别的表示，也不能随意毁掉与B的载体配对。本轮在固定扩展和同一边界约束下讨论这一结论；跨边操作、改变边界、带补偿荷的更大系统是不同合同。

若把带荷装置或环境留在A内，其边界表示也必须计入整个A。把荷移到被排除的环境中，只说明改了系统边界，并未完成所要求的有限输出。

## 3. 原模型确实具有无限边界支持

636已有完整原图中的字符环态。取跨越A、B的一块原面，标量为原规范不变光滑紧支撑径向包，费米子为原空态。对

$$
\lambda_n=(0,0,0,6n),\quad n=0,1,2,\ldots,\qquad
\Psi_n=\chi_{\lambda_n}(U_\square)\prod_v f(\phi_v)\otimes|0_F\rangle,
\qquad
\langle\Psi_n,\Psi_m\rangle=\delta_{nm}.
\tag{3}
$$

Q=6n满足原Z₆商条件，dλ=1、Cλ=36n²。每个Ψ_n都是完整Gauss态，属于原共同光滑核；不是完整H的本征态或不变部门。它们是原链路电表示，不是额外粒子。切开两条跨边后，边界元组为λ_n与相应对偶，互不相同。

原正热算符在物理Hilbert空间无核。令S为任意有限个边界元组，Π_B,S为对应远端投影。式(3)提供其补空间内非零物理向量，因此603的忠实Gibbs态给

$$
p_{\bar S}(\rho_\beta)
=\operatorname{Tr}\!\left[J\rho_\beta J^\dagger(I_A\otimes\Pi_{B,\bar S})\right]>0,
\qquad \rho_\beta=e^{-\beta H_F}/Z,\quad \beta>0 .
\tag{4}
$$

不需要假定边界标签与H_F对易，也不需要把H换成纯电Hamiltonian。标签可被原跨区域相互作用改变；“局部受限操作不能改标签”不是全动力学的新守恒律。本轮没有计算式(4)在完整热态中的数值。

## 4. 有限容量、远端不变和Gauss支持不能同时精确

令F_A为候选有限输出的边界不变有限秩投影，其出现的边界元组集合为S。有限维连续群表示只含有限个不可约分量。给定原物理输入ρ与任意归一输出Σ，定义

$$
\epsilon_A=\operatorname{Tr}[(I-F_A)\otimes I\,\Sigma],\qquad
\eta_G=\operatorname{Tr}[(I-P_{\rm phys}^{\rm cut})\Sigma],\qquad
\delta_B=\frac12\|\Sigma_B-\rho_B\|_1 .
\tag{5}
$$

F_A内没有与Π_B,barS配对的载体；这些投影与切口匹配投影相容，故在扩展空间上

$$
I_A\otimes\Pi_{B,\bar S}
\ \le\ I-P_{\rm phys}^{\rm cut}+(I-F_A)\otimes I,\qquad
\boxed{\eta_G+\delta_B+\epsilon_A\ge p_{\bar S}(\rho).}
\tag{6}
$$

证明第二式：第一式给输出远端补扇区概率≤η_G+ε_A；同一个投影在输入、输出边缘上的概率差≤δ_B。这个界不要求Σ来自通道，是对任何声称满足三份误差合同的输出的必要条件。

严格单侧保迹操作给δ_B=0，严格有限输出给ε_A=0。结合式(4)，它对原忠实热输入必有η_G>0。即使不要求所有未知态，只要求同一完整热态也不能三项同时精确。全局625／704压缩不受此反证，因为它允许改变两侧边缘；617等距也不受此反证，因为它没有有限维化。

这不是禁止近似。若S_L保留边界总Casimir不超过L的全部元组，637原区域比较能源含每个半边的c_T/2 Casimir，因此

$$
p_{\bar S_L}(\rho)
\le \frac{2\,\operatorname{Tr}(\rho_AK_A)}{c_TL}
\le \frac{2\,\operatorname{Tr}\rho(H_F+C_*)}{c_TL}.
\tag{7}
$$

固定图c_T>0，任意有限能源态的尾可趋零；但有限L的热尾严格正。L是近似标签，不是认知公理或物理硬阈值。有限标签集合也不等于有限重数维数；式(7)不独自解决全部容量。细化时c_T和能源预算必须另核。

## 5. 规范协变或只记标签都不够

规范协变仅要求Φ(U X U†)=UΦ(X)U†，并不等于每个Kraus满足式(2)。把A重置成其平凡表示真空是协变通道，但会把B原有非平凡边界荷留在一侧。

非Abel情形即使λ标签不变，载体配对也不能随意去掉。对原636两切口环态，载体维数Dλ=dλ²，在A载体施加完全去极化通道后，

$$
(\mathcal D_A\otimes{\rm id})(|\Omega_\lambda\rangle\langle\Omega_\lambda|)
=\frac{I_{D_\lambda}}{D_\lambda}\otimes\frac{I_{D_\lambda}}{D_\lambda},
\qquad
\operatorname{Tr}(P_{\rm phys}\,\Sigma)=D_\lambda^{-2}=d_\lambda^{-4}.
\tag{8}
$$

输出密度算符与联合群作用对易，但它大部分不在不变态向量的子空间。原L、u、Q表示的d为2、3、6，存活概率分别1/16、1/81、1/1296。保留经典标签不能冒充保留整个匹配载体。

## 6. 正向替代：保边界，只压缩重数

直接使用637构造的区域比较能源K_A。它不是原区域独立物理Hamiltonian；其热迹有限、预解紧，且与边界群对易。在式(1)中

$$
K_A=\bigoplus_\lambda K_{A,\lambda}\otimes I_{V_\lambda},\qquad
K_{A,\lambda}e_{\lambda j}=a_{\lambda j}e_{\lambda j},\qquad
\mu_\lambda=\min_j a_{\lambda j},\quad g_\lambda=e_{\lambda0}.
\tag{9}
$$

每个非空重数块有最低本征向量，允许简并时固定一个；K_A≥0。对N>0，pλ,N保留aλj≤N，P_A,N=1_{K_A≤N}。定义一份保持不同已保留扇区相干性的共同主Kraus，再加各扇区尾重置：

$$
K_{0,N}=P_{A,N}=\bigoplus_\lambda p_{\lambda,N}\otimes I_{V_\lambda},
\qquad
K_{\lambda j,N}
=|g_\lambda\rangle\langle e_{\lambda j}|\otimes I_{V_\lambda}
\quad(a_{\lambda j}>N),
\qquad
\Phi_{A,N}(X)=\sum_\alpha K_{\alpha,N}XK_{\alpha,N}^\dagger .
\tag{10}
$$

尾Kraus只在所标λ块作用。所有块可数，Kraus平方和强收敛给

$$
\sum_\alpha K_{\alpha,N}^\dagger K_{\alpha,N}=I_A,\qquad
[K_{\alpha,N},U_A(g)]=0,\qquad
(\Phi_{A,N}\otimes{\rm id}_B)(P_{\rm phys}\rho P_{\rm phys})
\ \hbox{仍有完整物理支持}.
\tag{11}
$$

这是真正作用A的正常CPTP通道，对任意B及被动参考保持其联合边缘。保留扇区之间的相干项都由单个K₀传递；没有先测出λ再整体去相干。尾的变化受下面误差控制，而非宣布它物理不可见。

输出支撑包含

$$
\mathcal F_{A,N}
=\bigoplus_\lambda
\left(p_{\lambda,N}M_{A,\lambda}+\mathbb C g_\lambda\right)\otimes V_\lambda .
\tag{12}
$$

每块有限，整体一般无限。特别是μλ>N时仍保留gλ⊗Vλ，不能把它删掉后继续声称通道保迹且保原边界。选择gλ和N属于分析性近似处方，未证明由原自主装置自动执行。

对p=1,2，直接在K_A本征基比较可得

$$
\Phi_{A,N}^*(K_A^p)
=P_{A,N}K_A^pP_{A,N}
+\bigoplus_\lambda\mu_\lambda^p(1-p_{\lambda,N})\otimes I_{V_\lambda}
\le K_A^p.
\tag{13}
$$

控制的是同一原比较能源，不能宣称完整相互作用H_F必然下降。原电剪切、磁项、跨边及全部CAR没有从H_F中删除。

## 7. 区域误差、参考与组合共用一个界

对任意包含B及未知参考R的正常联合态，设ε_N=Tr[(I−P_A,N)ρ_A]。gentle投影界加正尾的迹给

$$
\|(\Phi_{A,N}\otimes{\rm id}_{BR})(\rho)-\rho\|_1
\le2\sqrt{\epsilon_N}+\epsilon_N,\qquad
\epsilon_N\le\frac{\operatorname{Tr}\rho_AK_A}{N}.
\tag{14}
$$

可把右边再截到2。这里的参考一致性仅在给定区域能源预算的态类上；不声称全输入diamond范数趋零。所有正常固定态因P_A,N强趋I都逐态收敛，能源界给其中明确的定量范围。

另一侧按同样规则取Φ_B,N。两图表示上的实际单侧通道交换，每个Kraus均保匹配，所以

$$
\Phi_{A,N}\Phi_{B,N}=\Phi_{B,N}\Phi_{A,N},\qquad
\|\Phi_{A,N}\Phi_{B,N}(\rho)-\rho\|_1
\le2\sqrt{\frac{2E}{N}}+\frac EN,\quad
E=\operatorname{Tr}\rho_AK_A+\operatorname{Tr}\rho_BK_B.
\tag{15}
$$

采用恒等扩展省略下标。三角界中A通道不改变B的边缘，故使用同一输入预算；E由637同一H_F控制。式(13)还使K_A+K_B及其平方的期望不增加：平方中的交叉项为正的K_AK_B，逐侧Heisenberg不等式适用。没有据此证明新有效Hamiltonian或空间传播速度。

同一套gλ使N≤M时Φ_A,N Φ_A,M=Φ_A,N。故对真实扰动态与同一参考，利用已有638的数据处理和联合下半连续性，

$$
D(\Phi_{A,N}\rho\Vert\Phi_{A,N}\sigma)
\uparrow D(\rho\Vert\sigma)\quad(N\uparrow\infty),
\tag{16}
$$

有限相对熵时为有限值极限，无穷时为扩展值结论。此处是局部、保完整边界的版本；不是638全谱粗化的新证明，更不提供有限矩阵的相对熵算法，因为输出保留无限边界。真实态与参考必须用同一通道。

## 8. 三组复算及其边界

### 8.1 原商群的精确尾见证

对式(3)取原非Gibbs混态

$$
\rho_{\rm loop}=\sum_{n=0}^\infty 2^{-(n+1)}
|\Psi_n\rangle\langle\Psi_n|,\qquad
p_{n>N}=2^{-(N+1)}.
\tag{17}
$$

各固定次H作用在Ψ_n上的范数至多随n多项式增长：紧规范方向微分带来n次幂，标量包保持紧支撑，原系数和有限Fock矩阵在其支撑上有界。因此此态有全部固定阶能源矩，不借用非正规或无限能源输入。

保留A的n≤N并将其余重置到n=0，而不动B；该协变通道的Gauss泄漏恰为式(17)的尾。数值取0—12标签，将剩余几何尾收在12以使有限夹具归一；N=0、1、3、7分别精确给1/2、1/4、1/16、1/256。原四链路局部群变换核对字符不变。此夹具不是完整H的平衡态；Gibbs结论由式(4)。

### 8.2 原非Abel载体

对原L、u、Q构造两个切口群表示R(g₁)⊗R(g₂)*，实际匹配向量在联合动作下不变。L、u还用完整密度矩阵独立复核式(8)，Q用同一归一迹公式；不能把仅与群对易的密度当成Gauss支持。此组保持原Z₆允许表示，没有新增规范群。

### 8.3 重数通道及未知参考

取原ν、L、u边界载体，D分别1、4、9，每块给三维重数，另有二维被动参考。总区域维数42。重数能源明确声明为诊断输入：

$$
a_{\lambda j}=1+2C_\lambda+\Delta_j,\qquad
(\Delta_0,\Delta_1,\Delta_2)=(0,1,7).
\tag{18}
$$

原Casimir按636计算；Δ不是原完整K_A的已求本征值。以任意复重数振幅和匹配载体组成含参考的纯态，核所有Kraus、远端／参考边缘、匹配支持、两侧顺序及重数能源。

|N|低能主投影秩|真实迹范数误差|式(14)上界|比较能源〈K_A+K_B〉|
|---|---:|---:|---:|---:|
|10|3|1.70452818|2.15146503|38.09894060|
|23|11|1.37925274|1.57996107|38.22504450|
|36|24|1.12266529|1.21561546|38.64757798|
|44|42|约9.14×10⁻¹⁸|0|39.44453443|

这些中间误差并不小，不把它们当高精度物理近似。初始比较二阶矩约2399.11486281，每次通道后均不增。Gauss泄漏≤8.28×10⁻³³，远端／参考边缘误差≤1.73×10⁻¹⁷；Kraus交织和两侧顺序在该表示中为零残差。

N≥23的一个被保留跨扇区相干系数约.0037403743+.0041267068i，证明夹具没有先把全部λ去相干。有限夹具中N=44恢复身份；无限原族的收敛与每扇区有限、整体无限的区别由解析证明承担。

## 9. 对整体条件的真实削减

本轮把C02区域、C14原Gauss、C19同一热参考、C20压缩和C21内部载体接为共同合同：

$$
\{\text{严格单侧、有限总输出、原完整热支持、精确Gauss}\}
\quad\text{不能同时成立；}\qquad
\text{式(6)给允许近似时的必要误差账}.
\tag{19}
$$

替代族同时保区域组合、Gauss、未知参考及原比较能源，但必须留下边界表示与载体。因此实际可用的正向接口是

$$
\{\text{原区域参考能源、保边界重数压缩}\}
\Longrightarrow
\{\text{局部CPTP、精确匹配、参考保持、受控逐态近似}\},
\quad
\text{整体有限维与实际动态来源极限仍不在此结论内}.
\tag{20}
$$

这没有增加“宇宙必须有限维压缩”的认知原则，也没有反证规范场论、连续极限或统一计划。跨边协调、保边界无限载体、有限精度支持以及扩大内部装置是不同竞争方案，须分别核资源与来源。

旧空间382上界、383自由反向、384已消去的额外Lipschitz、386真实邻域桥或425半幅／成本桥、522—523热参考与实际方向仪器、524局域探针及UV有限噪声均直接复用。本轮没有改它们的前提，也没有把边界群的内部表示维数认作空间维数。辅助E仍不当记录s。

接[706](../../706/drafts/STATUS.md)：检验该保边界局部族如何接入原真实记录与来源，先核比较能源图范数和误差合同，不能直接照搬704中全H谱投影的对易证明。保留四分支、699反例范围和当前统一目标；不再重复有限容量障碍或调参优化式(14)。

