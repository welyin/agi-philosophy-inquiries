# 第717轮：实际记录的质量反馈与有界物理读数

日期：2026-10-03。接[716](research_note_716.md)、[717入口](717/drafts/postrecord_reference_entry.md)。[代码](717/joint_record_mass_feedback.py)、[结果](717/joint_record_mass_feedback_results.json)、[核验](717/research_round_717_checks.json)、[条件账](717/unified_physics_condition_ledger_717.md)。四组检查、二十二式；主代理推导与审查，无独立代理审查。

## 1. 实际推进与去重

上一目标轮完成716及717入口，属于有效进展。本轮先核导航、716报告／结果和入口，无在运行的Python研究进程。回查577、598、600、623—624、633—637、669、712及716：600已有质量力；623已有原域与目标多项式；669已有梯度增长；712已有二阶交换子工具。以下不把这些工具重新计为定理。

新增连接是：**原CAR记录在不改变即时玻色边缘的情况下，通过同一质量力改变后续玻色二阶响应。读取资源界不能统一控制所有无界响应，但项目原有的有界sin s读口具有明确、与输入态无关的二阶反馈界。** 后者提供可继续接有限时间物理记录的方向，无须先追加一串全部矩闭合公理。

|层次|本轮地位|
|---|---|
|认知动机|保留一份记录后，继续认知须用真实更新的整体态|
|继承输入|原固定图、全部H_b／CAR／Gauss、原K与U、原Y及跳跃、给定正几何|
|操作输入|716的配置无关sterile占据数仪器；577／598已有singlet读口|
|解析增量|原记录后质量反馈、正常Gauss反模型族、旧有界效果的全配置二阶上界|
|数值验证|原96模式、原H⁵正常波包求积、物理权重及原有界效果|
|未完成|统一有限时间余项、跨图真实参考／传播极限、因果装置、动态度规|

原动力学不被冻结，也不在读取后重新设为Gibbs。下文的加速度是某个量子可观测期望的二阶时间导数，不是物理空间中粒子的加速度，亦不是Einstein方程。

## 2. 真实原过程的二阶响应

沿716令ρ'=Φ_f(ρ)，R_f=1−2n_f，D_f=Φ_f*(H)−H。对于纯玻色配置函数A，

$$
A(t)=e^{itH/\hbar}Ae^{-itH/\hbar},\qquad
\Delta_A(t)=\operatorname{Tr}[(\rho'-\rho)A(t)],\qquad
\Delta_A(0)=\dot\Delta_A(0)=0 .
\tag{1}
$$

初始ρ先取原光滑紧支撑Gauss核的有限秩混合。其像仍在该核。A取下文P_q或有界效果，623保证A对H图范数有界；核向量属所有H幂域，因此原演化在D(H)中二次可微，给实际矩阵元的Taylor展开。没有宣称所有有限能量态都有同样的二阶展开。

这与成熟[Ehrenfest定理的严格域要求](https://arxiv.org/pdf/1003.3372)一致：Friesecke—Schmidt定理1要求相关算符域在演化下保持，结论为期望的一阶可微。本轮另用上述更强核正则性支持二阶步骤，不把其定理直接当成所有二阶结论。

写H=H_b+B，B为原完整矩阵乘法势。原节点动能为−ℏ²Δ_K/(2w_v)。因为[B,A]=0，Laplace乘积法则给

$$
[B,[H_b,A]]=\hbar^2\mathscr F_A,\qquad
\mathscr F_A=\sum_v\frac1{w_v}
K_v^{ab}(\partial_{v,a}A)(\partial_{v,b}B).
\tag{2}
$$

原规范电项对本轮A的导数为零；其余玻色势不参与这个内外交换子。原跳跃不依赖φ，所以此处只剩原质量导数，并不是删去跳跃演化。

由于[H_b,[H_b,A]]是纯玻色算符，其在ρ和ρ'中的初始期望相同，得到

$$
\ddot\Delta_A(0)
=-\operatorname{Tr}\rho\,\delta_f\mathscr F_A,\qquad
\delta_f X:=\Phi_f^*(X)-X,\qquad
\Delta_A(t)=-\frac{t^2}{2}\operatorname{Tr}\rho\,\delta_f\mathscr F_A+o_\rho(t^2).
\tag{3}
$$

ℏ²由二阶Heisenberg公式与(2)相消；这不代表没有量子效应，ρ、CAR相干与Φ_f仍在。余项依赖实际状态及固定图，尚无统一时长结论。

## 3. 原质量力的同一系数

复用623的内部双曲坐标x=φ/√F，t=|x|²、a=|x_H|²、b=x₅²：

$$
U=\frac14\delta^TL\delta,\quad
\delta=\left(I-\frac{u{\bf1}^T}{6M}\right)
\begin{pmatrix}a\\b\end{pmatrix}-\frac uM,\qquad
K_x^{-1}=I+\frac{xx^T}{6},\quad B_{\rm mass}=\sum_{j=1}^5x_j\mathcal B_j .
\tag{4}
$$

每个𝓑_j是原完整32模式的Dirac／Majorana常系数二次CAR算符。此处x是内部目标坐标，不能当物理五维空间。

取A=P_q=∑q_vU_v，q_v≥0、∑q_v=1。设J=(I+xxᵀ/6)∇_x U：

$$
J_H=2x_H\!\left[U_a+\frac{aU_a+bU_b}{6}\right],\qquad
J_5=2x_5\!\left[U_b+\frac{aU_a+bU_b}{6}\right],\qquad
\mathscr F_{P_q}=\sum_{v,j}\frac{q_v}{w_v}J_j(x_v)\mathcal B_{v,j}.
\tag{5}
$$

这是600中grad U收缩质量力的原Hamiltonian实现；新增的是实际记录前后的δ_f及其后续响应。其系数、体积和原Y不能再独立指定。

对716模式f，p_v=|f_v|²，δ_f𝓕的活跃单粒子支撑仍至多三模。其正常／配对交叉范数满足

$$
\|\delta_f\mathscr F_{P_q}(q)\|
\le |Y_\nu|\left[\sum_vp_v\left(\frac{q_v}{w_v}\right)^2|J_H(x_v)|^2\right]^{1/2}
+|Y_s|\left[\sum_vp_v\left(\frac{q_v}{w_v}\right)^2J_5(x_v)^2\right]^{1/2}.
\tag{6}
$$

这里式左侧括号q表示配置，权重q_v另列。可改称配置Q以避免混淆；数值代码始终分开存储。

623的1+t²≤C_c(1+U)、|U_a|,|U_b|≤ℓ+ht给全域充分界

$$
|J|\le C_J(1+U)^{5/4},\qquad
C_J=2(\ell+h)(2C_c)^{5/4},\qquad
\|\delta_f\mathscr F_{P_q}\|
\le (|Y_\nu|+|Y_s|)C_J r_q
\left[\sum_vp_v(1+U_v)^{5/2}\right]^{1/2},
\quad r_q=\max_v(q_v/w_v).
\tag{7}
$$

证明是|J|≤2√t(1+t/6)(ℓ+ht)以及(1+t)²≤2(1+t²)。这是宽松充分界，未称5/2是最佳矩指数。固定图还可用节点和给∥δ_f𝓕∥≤C_graph(1+W)^(5/4)，W=∑w_vU_v。

在原完整Gibbs中，623的W图范数界与热H²矩保证TrρβW²有限，故这个瞬时质量力期望有限；读后玻色边缘相同，相关配置矩仍有限。单凭这一点不领取ρβ'的全部高阶实时间Taylor性质或任何跨图一致性。

## 4. 有限读取资源不能控制所有势响应

这节给**同一固定原图内的状态族反例**，不是对任意参考、连续极限或统一计划的否定。取f在节点v有p_v>0。其余节点用固定光滑Gauss包、规范链路取常数Haar波函数。v节点取规范不变的正常包

$$
\chi_R(x)=N_R\,b_H(|x_H|)\,b_s(x_5-R),\qquad
b_H(r)=e^{-1/(1-r^2/d^2)}{\bf1}_{r<d},\quad
b_s(y)=e^{-1/(1-y^2)}{\bf1}_{|y|<1},\quad d=0.7 .
\tag{8}
$$

二者零延拓均光滑。归一使用原H⁵测度d⁵x／√(1+|x|²/6)，不使用平直测度。x_H径向性保证内部规范不变。

全部其余CAR为空，在v处使用原偶sterile相干

$$
|\eta_+\rangle=\frac{|0\rangle+e^{i\arg Y_s}
c_{v,30}^\dagger c_{v,31}^\dagger|0\rangle}{\sqrt2},\qquad
\rho_R=(1-R^{-4})\rho_0+R^{-4}|\chi_R\eta_+\rangle\langle\chi_R\eta_+|,
\quad R\ge4 .
\tag{9}
$$

ρ₀用R=0的同型包和相同费米向量。每个ρ_R都是原正常、偶、Gauss且核光滑的有限秩态；它们不是热态。

令P_mat=I−u1ᵀ/(6M)、u₄=(P_mat[:,2]ᵀLP_mat[:,2])/4>0，其中第2列对应b=x₅²。χ_R支持上x_H有界、x₅=R+O(1)，故

$$
\langle U_v\rangle_{\chi_R}=u_4R^4(1+o(1)),\qquad
\langle J_5\rangle_{\chi_R}=\frac{2u_4}{3}R^5(1+o(1)),
\qquad
\langle T_v\rangle_{\chi_R}=O(R^2).
\tag{10}
$$

动能阶来自原逆度量和固定宽度包的导数；归一后的坐标密度趋于固定紧支持密度。边势仅随目标距离平方增长，为O((log R)²)；规范磁势有界、电能为零，原质量为O(R)。因此在同一固定图、某个固定分析移位c下，

$$
\sup_R\operatorname{Tr}\rho_R(H+c)<\infty,\qquad
\sup_R\operatorname{Tr}\rho_R P_f<\infty,\qquad
\sup_R\operatorname{Tr}\rho_R|D_f|^4<\infty .
\tag{11}
$$

末式直接由716，非新的独立资源假设。非选择ρ_R'也有相同势矩和有界平均总能量。

在|η₊〉中，Dirac项期望为零，Majorana力期望为|Ys|J₅。原全局模式反射使本节点配对系数的δ_f为−p_v倍；这是二维反对称块与rank-one反射的行列式身份，不要求f仅在该节点。对q=p，

$$
\ddot\Delta_{P_f,R}(0)
=\frac{p_v^2}{w_v}|Y_s|
\left[(1-R^{-4})\langle J_5\rangle_{\chi_0}
+R^{-4}\langle J_5\rangle_{\chi_R}\right]
\sim \frac{2u_4p_v^2|Y_s|}{3w_v}R\longrightarrow+\infty .
\tag{12}
$$

χ₀关于x₅对称，首项为零。相同初始玻色边缘与一阶变化，连同(11)，仍不足以给无界势期望的统一二阶响应。

这不证明任意固定非零时刻有能量或势爆炸。各态的Taylor余项不统一，不能交换R→∞与t→0；也没有证明原Gibbs序列违反(7)。不继续制造更高阶矩反例。

## 5. 项目原有的有界读数有更强的正面结果

回到577／598已经使用的效果e_+(s)=1/2+sin s/4。对固定正报告权q，定义A_q=∑q_ve_+(s_v)。它可对应先选节点v、做原instrument再只报告正负的方案，Kraus为√q_v L_{r,v}；不是悄悄换成sqrt(A_q)的另一仪器。节点选择与实体装置的资源仍需内部实现，此处只算同一已声明效果的概率。

原逆度量与原质量坐标满足712身份的完整版本

$$
\sum_bK_\phi^{ab}\partial_b(\phi_j/\sqrt F)=\sqrt F\,\delta^a_j,\qquad
\nabla_K e_+(s)\cdot\nabla_K x_j
=\frac{\sqrt F\cos s}{4}\,\delta_{j5}.
\tag{13}
$$

因此Dirac质量力在这份二阶差量中精确消去，只剩原Majorana：

$$
\mathscr F_{A_q}
=\sum_v\frac{q_v\sqrt{F_v}\cos s_v}{4w_v}
\left(Y_s c_{v,30}^\dagger c_{v,31}^\dagger+
\bar Y_s c_{v,31}c_{v,30}\right).
\tag{14}
$$

这与712的singlet二阶消去不冲突：712追踪的是原Q³L创建张量的第一Yukawa通道，而这里是sterile占据数读取对原Majorana相干的改变。对象不同，原逆度量一致。

反对称配对块Δ_A给b=Δ_A bar f⊥f。δ_f𝓕只有f∧b的配对项，范数恰为∥b∥。由F≤M、|cos s|≤1，得到完整配置空间及物理Gauss子空间上的精确范数

$$
\|\delta_f\mathscr F_{A_q}\|_{\rm phys}
=\frac{|Y_s|\sqrt M}{4}
\left[\sum_vp_v\left(\frac{q_v}{w_v}\right)^2\right]^{1/2}
\le\frac{|Y_s|\sqrt M}{4}\,r_q .
\tag{15}
$$

上确界在各φ_v趋0时达到；此时f与b都属规范平凡sterile模块，取相应偶配对本征态及原点附近正常Gauss包可逼近该值。无需用带荷非物理向量充当物理范数见证。

对上述核上的全部输入态，真实概率差满足

$$
|\ddot\Delta_{A_q}(0)|\le
\frac{|Y_s|\sqrt M}{4}
\left[\sum_vp_v(q_v/w_v)^2\right]^{1/2}.
\tag{16}
$$

右侧不含势矩或态。相关二阶系数作为有界算符可延伸到所有正常态，但这仍不意味着所有正常态的实际时间概率都有统一二阶余项。定理的导数范围和系数延拓范围须分开。

## 6. 物理体积与几何来源继续共用

若q_v按同一物理体积取q_v=w_v g(x_v)/∑w g，g有界非负且分母有正极限，则q_v/w_v一致有界。由(15)可得到该实际效果二阶反馈系数的跨图统一界，条件为

$$
\sup_a \frac{\|g_a\|_\infty}{\sum_vw_vg_a(x_v)}<\infty
\ \Longrightarrow\
\sup_a\|\delta_{f_a}\mathscr F_{A_{q_a}}\|_{\rm phys}<\infty .
\tag{17}
$$

这不需要先求原全图Gibbs，也不证明该图序列真的逼近一个连续物理场。716的体积归一直接复用；维数没有由此选定。

沿均匀w_v=ε³exp(6σ)族，固定物理测试权的归一形状、φ、Y以及单位lapse，p、q不变。于是无界势力和有界效果力均有

$$
\partial_\sigma\delta_f\mathscr F_A=-6\,\delta_f\mathscr F_A .
\tag{18}
$$

非均匀几何时还需716的dot f、dot q项；不将(18)当全部应力，也不把背景体积变分当动态GR。共同点是：原记录、质量、曲目标、体积与来源使用同一组系数，没有新增“认知反馈强度”。

## 7. 可复算结果

四组检查使用现有Python／NumPy：

1. 原三节点96模式，在四组内部幅度中保留全部原质量。CAR夹断前后质量力差与(12)精确一致，原φ坐标方向差分相对误差随步长缩小。计算同时保留整个原H的解析消去项，不用点采样替代完整演化。
2. 原H⁵规范不变正常波包。以48、80阶求积比较，所报告矩的相对变化小于7.14×10⁻⁹；不是严格积分区间证书。R=4到128时，

$$
\begin{array}{c|ccc}
R&\langle U_v\rangle_{\rho_R}&
\langle w_vU_v+T_v+B_{M,v}\rangle_{\rho_R}&
\ddot\Delta_{P_f,R}(0)\\ \hline
4&0.06245205&26.44255368&0.02298402\\
128&0.06256744&26.41789516&0.52170533
\end{array}
\tag{19}
$$

表中第二能源列只计该节点三项；整个图的额外项依(10)后给出的解析界控制，未声称数值积分全图。实际渐近斜率为0.0040742202876，R=128的比值为0.0040758228566。

3. 同一物理紧支撑权在12³、20³、32³单元上，原势力交叉范数上界约0.00390167、0.00393981、0.00394219；这只是给定光滑配置的物理求积。均匀几何导数误差小于1.43×10⁻¹¹，未模拟量子参考。
4. 原sin s效果在全96模式中的配对差量只需两活跃模。精确全配置范数

$$
C_{\rm effect}=0.05589302508195645,
\qquad
\|\delta_f\mathscr F_{A_q}\|_{\phi=0}=C_{\rm effect}.
\tag{20}
$$

五组配置的原逆度量收缩误差小于4.1×10⁻¹⁷；非零配置的实际范数均在该界内。由全域证明及φ→0见证得到精确范数，非采样推测。

早期三组势力结果保存在[草稿结果](717/drafts/force_results_before_bounded_record.json)。范围元数据中一处非热标识在正式保存前修正；当前四组结果与代码一致，未覆盖任何既有冻结轮次。

## 8. 结论与下一项

本轮同时给出一条限制和一条可用连接：

$$
\text{有界总能量＋716读取资源}
\not\Longrightarrow\text{任意无界势响应一致受控};
\qquad
\text{原曲目标＋Majorana＋物理报告权}
\Longrightarrow\text{原有界读数的二阶反馈界}.
\tag{21}
$$

新结果约束C03／C15／C19／C20／C22，不消除原Y、物种或给定几何的输入。633瞬时分布式实现、634来源拼接及699固定候选反例保持。

下一步不继续高矩反例或波包优化。接[718入口](718/drafts/STATUS.md)，将原实际有界读数接到有限时间共同历史。因R_f与A_q对易，可以先复用精确身份

$$
\Delta_{A_q}(t)=\frac12\operatorname{Tr}\rho
\left[e^{it(H+2D_f)/\hbar}A_qe^{-it(H+2D_f)/\hbar}
-e^{itH/\hbar}A_qe^{-itH/\hbar}\right],\qquad
H+2D_f=R_fHR_f .
\tag{22}
$$

右侧只是相同物理读取的比较表示，不增设一套实际自然规律。须回查623—625、634和704已有真实响应，核有限变化、真实状态与几何来源，不以二阶系数代替全历史。旧空间382—386、425、522—524逐项复用：384删除的Lipschitz不恢复，386／425保持替代，523／524不重复证明；共同物理端点身份仍待落实。统一目标继续开放。
