# 第665轮：原Gauss量子过程、接触作用与辅助历史的共同条件

日期：2026-10-02。接[664](research_note_664.md)、[实际入口](round665_drafts/auxiliary_current_entry.md)。[代码](joint_contact_gauss_history.py)、[结果](joint_contact_gauss_history_results.json)、[核验](research_round_665_checks.json)、[全条件账](unified_physics_condition_ledger_665.md)。三组正式检查、十八式；主代理推导与审查，无独立代理审查。

## 1. 本轮实际接通的范围

664指出最小独立联络产生总流接触作用，单节点原质量不再保持二次费米计算形式。本轮将其接回**原固定有限图的完整玻色、规范、CAR量子系统**：在保留原标量目标的声明分支中，新的有界接触项与原共同域、Gauss热态、实际标量记录及几何来源相容。并给出接回原有序Gaussian条件历史的辅助表示，明确原行列式公式需要补的提升因子。

保留664的分支区分：本轮讨论C_E，或已经补齐标量匹配项而保留原K的接触扩展；不把C_J的新标量目标直接塞进旧完整性证明。无挠基线H仍保留，扩展Hᶜ通常是不同模型；反向接触匹配项才使有限调节下两者相同。

上一目标轮完成并发布664、执行665入口，属于进展。本次先核README、方向、状态、664报告与结果及保存入口；无在运行的Python进程。回查[603](research_note_603.md)、[623](research_note_623.md)—[624](research_note_624.md)、[643](research_note_643.md)、[655](research_note_655.md)和[662](research_note_662.md)。不重做旧热核、原记录定义或一般Trotter定理；新增的是接触作用加入后这些条件能否共同成立，以及实际辅助字典的来源。

|层次|本轮内容|
|---|---|
|认知动机|状态、相互作用、读取和反作用须由同一过程承载|
|继承输入|固定有限图、原H⁵目标和束缚势、原商群／32CAR、完整质量与跳跃、正给定几何|
|新增分支输入|664指定正常排序的接触作用；辅助Gaussian分解仅为表示|
|解析连接|原完整域／Gauss热态／记录与新增来源共同存在；归一辅助表示与原有序权重相接|
|数值范围|实际32模式的条件提升；八轻子的辅助平均；原两节点中性不变扇区|
|仍未证明|原联络路径测度唯一性、全图连续手征极限、动态量子几何与GR|

## 2. 原完整域与Gauss热态继续适用

固定一代32模式／节点，令Vᶜ为664的正常排序总流接触项，w_x为原正单元体积：

$$
H_\gamma^c=H_\gamma+V_\gamma^c,\qquad
V_\gamma^c=\sum_x k_x:\widehat Q_x:,
\quad k_x=\frac3{16w_x},\qquad
\|V_\gamma^c\|\le204\sum_xw_x^{-1}=:v_\gamma.
\tag{1}
$$

粗界由∥Ĵ₀∥≤16、∥Ĵ_i∥≤16、∥N∥≤32及:Q̂:=−Ĵ₀²+ΣĴ_i²−2N给出1088，再乘3/16。全部物种均在内，不是ν扇区的界。正几何紧参数域使Vᶜ及其指定一、二阶参数导数一致有界。

原内部群逐点与四个流系数对易，故Vᶜ保Gauss及宇称。有界自伴扰动保623的共同算符域和原紧预解；在同一物理子空间上用min–max比较本征值，得到

$$
\mathcal D(H_\gamma^c)=\mathcal D(H_\gamma)=\mathcal D,
\qquad [P_G,H_\gamma^c]=0,\qquad
0<e^{-\beta v_\gamma}Z_{\rm phys}
\le Z^c_{\rm phys}\le e^{\beta v_\gamma}Z_{\rm phys}<\infty.
\tag{2}
$$

不使用一般不成立的算符指数单调性。Hᶜ全部热能量矩由同一正β谱估计有限；原核在有界扰动后仍为算符核。改变原标量目标的未匹配C_J、w可趋0或图规模无界，均不在此命题内。

## 3. 记录不必另配一套能源，但概率会改变

原L_r(s)=√(1/2+r sin s/4)在Fock上为身份，与新增接触项对易。因此

$$
[V_\gamma^c,L_r]=0,\qquad
\sum_rL_rH_\gamma^cL_r-H_\gamma^c
=\sum_rL_rH_\gamma L_r-H_\gamma,
\qquad L_r\mathcal D\subset\mathcal D.
\tag{3}
$$

这是同一瞬时非选择注能算符相同，绝不意味着两模型的平衡态、时间历史或读数概率相同。以Hᶜ代入624的实际实时间历史、并采用同一初态或明确的新Gibbs态，仍得到归一CP仪器及有限历史域控制。存在仪器的数学表示不等于自主装置和资源机制已完成。

同一几何来源需补新增项。对保原canonical CAR识别的体积参数，

$$
\partial_\lambda H^c=\partial_\lambda H+\sum_x(\partial_\lambda k_x):\widehat Q_x:,
\qquad
\partial_\lambda\log Z^c_{\rm phys}
=-\beta\langle\partial_\lambda H^c\rangle_{\beta,c}.
\tag{4}
$$

若还改变旋量识别，应另保663的字典导数。式(4)的完整原∂H包含玻色动能、势和跳跃，不能用条件费米得分替代全部几何应力。有界的新增导数保623的共同域时间演化及已有标量二阶响应论证。

## 4. 接触作用的共同辅助分解

令Qχ=c†χc、S_i=c†σ_i c/2，均对原全部物种求和。N、Qχ与各S_i对易，而不同S_i一般不对易。664正常排序给

$$
:\widehat Q:=4\mathbf S^2-Q_\chi^2-2N,\qquad
e^{-\delta V_x^c}=e^{2aN}e^{aQ_\chi^2}e^{-4a\mathbf S^2},
\quad a=\delta k_x>0.
\tag{5}
$$

[Assaad的辅助场研究](https://arxiv.org/pdf/cond-mat/9806307)示范了相互作用的辅助二次表示及保对称性的差别；这里不套其Hubbard物种或具体离散恒等式，而对原共同流直接作谱证明。标准正态Dz给

$$
e^{aQ_\chi^2}=\int Dz\,e^{\sqrt{2a}\,zQ_\chi},\qquad
e^{-caS_i^2}=\int Dz\,e^{i\sqrt{2ca}\,zS_i},\qquad
Dz=(2\pi)^{-1/2}e^{-z^2/2}dz.
\tag{6}
$$

使用对称spin分割

$$
P(a)=e^{-2aS_1^2}e^{-2aS_2^2}e^{-4aS_3^2}
e^{-2aS_2^2}e^{-2aS_1^2},\qquad
P(a)>0,\quad P(a)=e^{-4a\mathbf S^2}+O(a^3).
\tag{7}
$$

有限CAR上余项由有界矩阵乘积展开控制；这是成熟对称乘积公式的使用。有限a的P不具精确全Spin旋转不变性，连续细分恢复指定总Casimir。内部标准模型规范群则逐辅助配置都保持，因为它与N、Qχ、S_i分别对易。

一个共同charge字段和五个spin字段给实际二次Fock算符

$$
\mathcal A_a(z)=e^{2aN}e^{\sqrt{2a}z_0Q_\chi}
e^{i\sqrt{4a}z_1S_1}e^{i\sqrt{4a}z_2S_2}
e^{i\sqrt{8a}z_3S_3}e^{i\sqrt{4a}z_4S_2}e^{i\sqrt{4a}z_5S_1},
\quad \mathbb E\mathcal A_a=e^{2aN}e^{aQ_\chi^2}P(a).
\tag{8}
$$

这不是六种新物质或六个外部控制器。辅助平均与原质量、跳跃按真实时间顺序交错；入口已核质量与Qχ、N不对易，两节点又核局部spin与原跳跃不对易，不能将所有辅助量预先合为一次静态修正。

## 5. 原有序桥、Gauss闭合与可积性

各spin因子酉，且与N、Qχ对易，因此条件于spin字段，正态charge平均给

$$
\mathbb E(\mathcal A_a^\dagger\mathcal A_a)
=e^{4aN+4aQ_\chi^2}\le e^{1152a}I,
\qquad
\mathbb E\|\mathcal A_a B\|_{\rm HS}^2
\le e^{1152a}\|B\|_{\rm HS}^2.
\tag{9}
$$

不同节点同片仍独立积分，但流不能在同一节点按物种分别积分。逐片迭代式(9)的指数与Σa=βΣk有关，不随片数发散。这比直接逐片E∥A∥的粗界更强，后者可能产生虚假的√δ累积。

固定一条原非爆炸玻色桥q_t，质量B(q_t)有界；保原时间顺序的辅助乘积均值收敛到生成元B(q_t)+Vᶜ的矩阵演化。外面的原标量权及598质量界共同提供

$$
e^{-\int_0^\beta V_\gamma(q_t)dt}
\mathbb E_z\left|\operatorname{Tr}\{\Gamma(g)\mathcal U_{F,z}[q]\}\right|
\le d_F e^{C\beta+576\beta\sum_x k_x}
e^{-\frac12\int_0^\beta W_\gamma(q_t)dt}.
\tag{10}
$$

这里用HS二阶矩、Cauchy–Schwarz及∥B∥≤C+cW^{1/4}，再以Young界吸收进W/2。对实际有序ODE可先沿路径作子步极限；不将端点Riemann和未经支配直接当连续泛函。603／643的标量桥支配保证接回原规范扭闭合积分后可积。因而同一完整热迹可写为

$$
Z^c_{\rm phys}
=\int_{\mathcal G}dg\int_Qd\mu(q)
\int D_{\gamma,\beta}^{q\to g^{-1}q}(d\omega)\,
e^{-\int V_\gamma(\omega_t)dt}
\lim_{\delta\to0}\mathbb E_z
\operatorname{Tr}_{\mathcal F}[\Gamma(g)\mathcal U_{F,z}^{(\delta)}[\omega]].
\tag{11}
$$

式(11)保留643同一个群元对应的玻色端点及Fock闭合。等价性按辅助平均后的矩阵极限及支配定义；不声称存在处处正的单路径权，也不将Euclidean桥当实际记录概率。实时间记录仍用式(3)的原仪器和Hᶜ。

## 6. 旧行列式公式必须补上真实提升因子

原643质量／跳跃的正常h无迹，公式可不带额外标量。本辅助分解新增e^{2aN}。令r为原单粒子辅助乘积，非酉情况下canonical Nambu作用必须写成

$$
\Gamma(r)=\exp(c^\dagger\log r\,c)\ \text{（连续提升）},\qquad
\mathsf R(r)=\operatorname{diag}(r,r^{-T}),\qquad
\det r=e^{2an},\quad n=32\ \text{每节点}.
\tag{12}
$$

一般乘积以实际外代数表示定义，不能任取矩阵对数的分支。spin与手性部分迹为0，e^{2aN}贡献全部det r。原配对热段的正常h仍无迹，所以共同历史迹满足

$$
\mathcal W_{F,z}^{,2}
=\exp\!\left(2n\sum_{x,j}a_{xj}\right)
\det(I+\mathsf M_z),\qquad
\mathcal W_{F,z}=\exp\!\left(n\sum_{x,j}a_{xj}\right)
\operatorname{Lift}\sqrt{\det(I+\mathsf M_z)}.
\tag{13}
$$

Lift沿原Fock乘积固定符号和相位，不是主平方根，也不对权重取绝对值。若再加其它有迹正常块，其提升因子也要算入。此处不是修改原物理真空；它是式(8)这个辅助表示不可省的标量。

当所有相关w=e^{6θ}倍增时，同一非零条件权的来源为

$$
\partial_\theta\log\mathcal W_{F,z}
=\frac12\operatorname{tr}[(I+\mathsf M_z)^{-1}\partial_\theta\mathsf M_z]
-6n\sum_{x,j}a_{xj}.
\tag{14}
$$

原全部32模式、完整复质量及非平凡商群闭合的直接“24夸克行列式×八轻子Fock迹”，与式(13)平方的相对误差1.93×10⁻¹⁴。a=0.035时漏提升会漏掉trace因子3.06485420及来源−6.72；共同来源与独立差分误差9.98×10⁻⁸。原643已说明有迹块需要补项，本轮给出实际新增块及其原物质结果，不宣称一般提升公式是新定理。

## 7. 辅助测度归一与几何导数必须一起输送

采用标准z时Dz与θ无关，只微分式(8)幅度。若改用物理辅助变量y=√(b a)z，测度宽度随几何改变，a′=−6a给

$$
d\nu_a(y)=\frac{dy}{\sqrt{2\pi b a}}e^{-y^2/(2ba)},\qquad
\left.\partial_\theta\log d\nu_a(y)\right|_y
=3-\frac{3y^2}{ba}=3-3z^2.
\tag{15}
$$

Gaussian积分分部证明该得分与标准z下振幅求导相同。每个字段漏掉归一导数中的+3，都会产生−3倍原平均；本时间片共有六个字段，故

$$
(\partial_\theta\log Z_{\rm step})_{\rm omitted}
-(\partial_\theta\log Z_{\rm step})_{\rm full}=-18.
\tag{16}
$$

这是固定本分解的漏项误差，不是新的宇宙学常数。八轻子24点Gaussian积分验证两种正确求导误差≤7.54×10⁻¹⁵，平均算符误差≤5.98×10⁻¹⁴，漏项确为−18。平均单片最小特征值0.71377；spin分割相对精确接触指数还有0.00141222误差，保留而不冒称有限片精确。

逐片得分身份不保证其连续随机得分有统一绝对变差界。完整来源仍由式(4)的算符Duhamel导数定义，辅助平均后再取已控制的极限；不以未经证明的路径测度微分代替它。

## 8. 原空间跳跃要求补接触能量流

在655实际两节点中性扇区，取零Higgs、s₁=0.62、s₂=−0.41和原允许的几何无关跳跃0.23I。每节点保留Majorana质量；令J为整条跳跃，H_x⁰=M_x+J/2，则

$$
\frac{[H_1^0+V_1^c,H_2^0+V_2^c]}i
=\frac{[H_1^0,H_2^0]}i
+\frac{[V_1^c,J/2]+[J/2,V_2^c]}i.
\tag{17}
$$

不同节点的接触及质量偶算符对易。最后一项不是旧二次Nambu能源流的重命名，实际Fock范数为0.92091749；原矩阵恒等式误差1.12×10⁻¹⁶。局部S₁与跳跃对易范数0.46，而此指定spin标量跳跃和质量的全局S₁+S₂仍守恒。

θ₁=0.13、θ₂=−0.08、β=0.83时，同一四模式热态给

$$
\log Z=3.93131499206,\qquad
\partial_{\theta_1}\log Z=-1.26651608984,\qquad
\langle j_{12}\rangle\simeq0,\quad
\langle j_{12}^2\rangle=0.07542254653.
\tag{18}
$$

该条件源的独立差分误差6.46×10⁻⁹；在本允许的诊断跳跃下，费米几何变化来自同一V₁。没有把这四模式计算当完整玻色或Gauss热积分；它核的是原实际空间块的新流与来源。

## 9. 整合成果与下一步

三组正式检查分别核全32模式辅助提升、同测度的平均／来源、实际空间接触流；665入口的Casimir和时间片检查只复算、不重复计组。固定图下，无需为声明接触分支另选热态规则、记录域、空间来源或Gauss闭合：它们由同一个扩展H和实际辅助字典固定。

这一整合没有强制引入接触项，更没有将辅助变量认作宇宙底层实体。有限调节下采用同处方的相反接触项，会在算符层恢复原H、态和所有原记录；连续独立联络测度的行列式、计数、边界及匹配仍另待证明。正几何下的常数随w⁻¹增长，不能用于声明跨图连续或强引力控制。

四个主分支仍未完全统一：原固定图量子过程、声明连续物质、给定作用的经典动态几何、原手征辅助测度。下一步回到它们的实际字典，接[666原CAR—全左手表示的排序与共同测度](round666_drafts/STATUS.md)：回查614粒子—空穴变换及615—661原手征权重，确认本接触与真空排序在同一表示中如何输送。先处理原对象映射，不再扫描辅助精度或设计认知装置。
