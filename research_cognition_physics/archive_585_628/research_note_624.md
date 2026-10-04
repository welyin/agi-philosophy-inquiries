# 第624轮：实际记录、条件响应与共同来源过程

日期：2026-10-01。接[623](research_note_623.md)，回查[586](research_note_586.md)、[592](research_note_592.md)、[593](research_note_593.md)、[598](research_note_598.md)、[617](research_note_617.md)及[586双历史方法桥](586/closed_time_path_bridge_review_586.md)。[代码](624/joint_record_source_response.py)、[结果](624/joint_record_source_response_results.json)、[核验](624/research_round_624_checks.json)、[条件账](624/unified_physics_condition_ledger_624.md)。三组检查、十四式；主代理审查，无新增独立代理审查。

## 1. 这次合并什么，哪些不重做

623已将原完整固定图的无记录来源接到实际过程。本轮把原基本仪器插入同一过程，得到**记录概率、条件量子过程和来源响应可以共同定义的有限历史二阶展开**。筛选某份记录时，记录概率本身也随参数改变，不能把无记录公式中的态直接替换为“读后态”并忽略这个变化。

592—593已经给原仪器的算符域与H²预算，586已经给配置效果的Fisher条件分解；这些均继承，不计新发现。新增接口是：完整CAR模型的共同来源过程带上实际记录后，怎样保留归一、条件选择项和来源信息的尺度验收。

|层次|本轮地位|
|---|---|
|认知动机|记录与物理响应不能来自互不相容的过程|
|继承输入|623原完整H及指定参数族；原L_±；固定有限记录时刻；同一初态|
|状态要求|Tr(A²ρ)<∞；包括603 Gibbs态及原光滑紧支撑Gauss包|
|解析新增连接|有限记录双历史的二阶展开、条件归一与原来源选择项|
|数值核对|原H生成的瞬时报告／条件响应，以及报告合并后的来源分辨性|
|未完成|原仪器的自主装置实现、长期统一界、连续尺度及动态引力|

成熟方法上，[Gammelmark—Julsgaard—Mølmer的式(2)](https://arxiv.org/pdf/1305.0681)用前向态和未来记录的反向效果计算过去测量的条件概率。这里复用其前向／反向区分，以原有限历史算符直接定义对象，不引入独立Markov环境。反向效果是推断工具，未使物理时间倒退。

## 2. 原记录确实能进入623的域

固定节点的原仪器及共同移位A=H+c≥1为

$$
L_r(s)=\sqrt{\tfrac12+\tfrac r4\sin s}\,I_{\mathcal F},\quad r=\pm1,
\qquad \sum_rL_r^2=I,\qquad \tfrac14I\le L_r^2\le\tfrac34I.
\tag{1}
$$

继承592的具体交换子：只有该节点动能参与，[H,L_r]=−ℏ²(2L_r′grad s·grad+Δ_KL_r)/(2w)。全部原矩阵质量、跳跃、标量势与L_r对易，因此新增CAR没有增加交换子项。598的半有界势控制又使节点梯度形式≤C〈A〉；623已经给共同核与算符域，于是

$$
\sum_r\|[H,L_r]\Psi\|^2\le c_1\langle\Psi,A\Psi\rangle+c_0\|\Psi\|^2,
\qquad L_r\mathcal D\subset\mathcal D,
\qquad \|AL_r\Psi\|\le C_L\|A\Psi\|.
\tag{2}
$$

这是将旧交换子预算接到完整CAR族，未重新证明592全部结论。来源参数不改变L_r定义；若将来改为参数依赖仪器，必须另加其导数。

设固定k次读、U_j为623在相邻记录间的实际传播子，包含最后一次读后的可选等待。把最后等待并入K即可，记

$$
K_{\boldsymbol r}[\lambda]=L_{r_k}U_k[\lambda]\cdots L_{r_1}U_1[\lambda],
\qquad p_{\boldsymbol r}[\lambda]=\operatorname{Tr}(K_{\boldsymbol r}\rho K_{\boldsymbol r}^\dagger),
\qquad p_{\boldsymbol r}\ge4^{-k},\quad \sum_{\boldsymbol r}p_{\boldsymbol r}=1.
\tag{3}
$$

下界来自每步L_r的最小奇异值≥1/2，适用于任意归一输入及任意等待。各K及其伴随保持D，有固定有限历史的图范数界。故每个条件后态仍有有限A²矩；没有与k无关的长期下界，不能借此省略原记录的能源成本。

## 3. 带实际记录的同一双历史对象

对同一结果串的两条源历史定义

$$
Z_{\boldsymbol r}[\lambda_-,\lambda_+]
=\operatorname{Tr}(K_{\boldsymbol r}[\lambda_-]^\dagger K_{\boldsymbol r}[\lambda_+]\rho),
\qquad Z_{\boldsymbol r}[\lambda,\lambda]=p_{\boldsymbol r}[\lambda].
\tag{4}
$$

共同历史的单分支对象等于该记录概率，**不是1**。只有对记录求和才有Σ_r Z_r[λ,λ]=1。即使对记录求和，两个不同源历史下的Σ_r Z_r也一般不等于完全未执行测量的623对象，因为非选择记录仍有反作用。

将有限乘积逐项作Duhamel差：一个来源插入时，右侧基点历史保持D，左侧都是有界算子；两个来源插入时，将较晚的一项移到bra侧，之间的U与L保持有界，两端基点历史及其伴随保持D。623式(16)的同一证明因此给

$$
\left|\langle\Xi,(K_\epsilon-K_0-\epsilon K_1)\Psi\rangle/\epsilon^2\right|
\le C_k\|A\Xi\|\|A\Psi\|,
\qquad (K_\epsilon-K_0)\Psi/\epsilon\longrightarrow K_1\Psi.
\tag{5}
$$

弱二阶极限包括同段两插入、跨段插入和实际Hessian接触项。记录算子按原时序保留，不交换到来源之外。对具有有限A²矩的ρ作谱分解，非零权重的向量属于D且Σ_jp_j‖Aψ_j‖²有限；以式(5)支配可将二阶展开移入Z_r的迹。这一步不要求ρ与H对易，也不要求ρ为Gaussian。

因此623的实际来源、原记录与有限历史在同一模型内相容。此处不另宣称完整路径积分测度、强二阶传播子或任意无界读数的二阶响应；已经足够定义记录概率和双历史来源系数。

## 4. 条件归一需要保留记录概率

设某个有界末端读数O的未归一分支均值为M_r=Tr(O K_rρK_r†)，μ_r=M_r/p_r。沿同一源参数的实际一、二阶系数满足

$$
\mu_r'=\frac{M_r'}{p_r}-\mu_r\frac{p_r'}{p_r},\qquad
\mu_r''=\frac{M_r''}{p_r}-\mu_r\frac{p_r''}{p_r}
-2\frac{p_r'}{p_r}\mu_r'.
\tag{6}
$$

这里对M的二阶结论要求O及O†保持D，或其插入满足同样弱展开条件；不由“有界”一词自动保证二阶性质。原光滑s读数在此固定目标域上有界，且有受控交换子，满足要求。对单独p不需这个附加读数条件。

对623源求导，未来记录由效果E(t)向过去输送，过去记录形成未归一σ(t)。时刻t不位于瞬时读取点时，一阶Duhamel密度是

$$
p_r=\operatorname{Tr}(E_r(t)\sigma_r(t)),\qquad
\frac{\delta p_r}{\delta\lambda_a(t)}
=-\frac i\hbar\operatorname{Tr}\big(E_r(t)[G_a(t),\sigma_r(t)]\big),
\qquad \partial_a\log p_r=\frac{\partial_ap_r}{p_r}.
\tag{7}
$$

G_a(t)在这里指该截面的Schrödinger来源，时间标签不是重复作Heisenberg输送。所有过去、未来操作都已放进σ及E。公式来自允许Hamiltonian族的光滑源变化，不擅自假定每个对称G_a都有可执行的独立瞬时酉exp(−iλG_a)。域／矩控制使相应迹配对有意义。

对全部未来结果求和，E之和是I，交换子迹消失，恢复过去边缘不受以后未筛选操作影响。固定某个未来结果时，其条件推断可变化；这不是逆因果传播，也不能套用无记录平衡态的纯迟致核来替代全部条件项。

## 5. 条件过程的正性与可分辨变化

在Hilbert–Schmidt空间令v_r(λ)=K_r[λ]ρ^{1/2}。则Z_r为v的Gram核。若要一个对角恒为1的条件比较核，必须明确除去两端各自的记录概率：

$$
\widetilde Z_r[\lambda_-,\lambda_+]
=\frac{Z_r[\lambda_-,\lambda_+]}{\sqrt{p_r[\lambda_-]p_r[\lambda_+]}},
\quad |\widetilde Z_r|\le1,\quad \widetilde Z_r[\lambda,\lambda]=1.
\tag{8}
$$

不能在使用它后丢掉p_r；两者共同描述记录及分支过程。单源方向记dot v为一阶导数，内积用HS迹。归一向量v/√p的二阶横向衰减给

$$
g_r=\frac{\|\dot v_r\|^2}{p_r}-\frac{|\langle v_r,\dot v_r\rangle|^2}{p_r^2}\ge0,
\qquad \operatorname{Re}\log\widetilde Z_r[-\epsilon/2,+\epsilon/2]
=-\tfrac12\epsilon^2g_r+o(\epsilon^2).
\tag{9}
$$

这是当前原Kraus表示和固定初态纯化的Gram几何，不直接认作可访问的混态量子Fisher信息，更不是物理时空度规。正性是同一完整过程的结果，不能分别拟合记录概率和噪声而不检验Gram兼容。

令s_r=p_r′/p_r为记录得分，a_r=Im〈v_r,dot v_r〉/p_r。直接分解复内积的实虚部：

$$
I_{\rm record}:=\sum_rp_rs_r^2\le4\sum_r\|\dot v_r\|^2<\infty,
\qquad \sum_r\|\dot v_r\|^2
=\tfrac14I_{\rm record}+\sum_rp_r(g_r+a_r^2).
\tag{10}
$$

有限性由式(5)的一阶图范数界及初态A²矩保证。该身份不是普遍量子Fisher链式定理；它说明同一源的记录概率、条件幅度与相位变化共用同一有限变化预算。

## 6. 原H和原Gauss包上的非零选择项

数值使用590的正规紧支撑Gauss包Ψ₀(h,s)，给singlet相位Ψ=e^{iκs}Ψ₀，κ=3。其它节点取原光滑Gauss因子，链路常数及有限费米真空是合法旁观输入。原全部势和CAR质量与f(s)对易，其它动能也不参与，因此完整H在该读前截面的真实瞬时生成元给

$$
\left.\frac d{dt}\langle f(s)\rangle_{e^{-itH/\hbar}\rho e^{itH/\hbar}}\right|_{t=0}
=\int p(h,s)\,v(h,s)f'(s)\,dhds,
\qquad v=\frac{\hbar\kappa}{w_B}\mathcal K^{ss}.
\tag{11}
$$

这里v是目标singlet概率流的系数，不是宇宙空间传播速度。计算的是原实际时间演化在指定正常态处的局部导数，没有数值求全图传播，也没有将其它部门从模型删除。读取在这个小时间变化之后进行；“先读后演化”是另一个过程。

对e_r=L_r²，原s条件均值的响应为

$$
p_r'=\int pve_r',\qquad
\mu_r'=\frac{\int pve_r}{p_r}
+\frac{\int pv(s-\mu_r)e_r'}{p_r},
\qquad
\frac d{dt}\sum_rp_r\mu_r=\sum_rp_r\mu_r'+\sum_rp_r'\mu_r.
\tag{12}
$$

第一项是若只用固定读后态计算s速度所得结果；第二项是时间变化发生在筛选**之前**带来的实际选择修正。κ=3、旧ℏ=.7、w_B=.8，64／96／128点求积给：

- p_+=.5972287741、p_−=.4027712259，p_±′=±1.1473920758。
- 两分支的选择修正为−.004691355606和−.0000884041323，均非零。
- 不筛选的s速度4.9886128134；Σp_rμ_r′=4.9857753942，缺少的Σp_r′μ_r=.002837419198恰好补齐。
- 96与128点对报告的关键组合差≤1.78×10⁻¹⁵。独立使用原曲Laplace对[H,e_+]作有限差分，步长.0005到.000125的当前系数误差从1.38×10⁻⁵降至8.65×10⁻⁷，包含原测度漂移。

这排除“保留各条件态的响应但固定记录权重即可恢复共同响应”这一具体接法；没有排除正确保留概率变化的条件模型。

## 7. 对尺度匹配增加一个实际验收条件

将若干实际细历史r合并为a，P_a=Σ_{r∈a}p_r，s_a=P_a′/P_a。586已证明经典Fisher分解；此处仅对**原H产生的参数变化**使用它：

$$
I_{\rm fine}-I_{\rm coarse}
=\sum_a\sum_{r\in a}p_r(s_r-s_a)^2\ge0.
\tag{13}
$$

只要同一粗类的得分不相同，就不能仅从该粗报告恢复原细记录对这一源的全部一阶分辨性。若在同一无等待记录串后只报告奇偶，原包诊断给

$$
\begin{array}{c|rr}
\text{读取次数}&I_{\rm fine}&I_{\rm parity}\\\hline
1&5.47298845&5.47298845\\
2&10.93975712&.79224342\\
3&16.40031659&.06880932\\
4&21.85467743&.00486078
\end{array}
\tag{14}
$$

各串都在同一次原H微小演化之后读，串内等待为零；不是把每次结果当独立样本。这里只合并实际报告，未用直接奇偶平方根仪器替换原过程，也未从Fisher差额推断新的能源下界。若尺度目标只需粗报告，则无需保留细报告的全部分辨性；如果声称保留原细来源信息，则必须满足对应得分／误差合同。

## 8. 核验与下一步

三组检查覆盖原条件响应与归一、独立曲生成元诊断、实际源下的历史合并。首次执行有一处不必要的数值量级预期失败：要求每个选择修正都超过10⁻⁴，而第二支是8.84×10⁻⁵。模型参数不改，改核非零幅度超过计算量级100个浮点间距后全部通过，详见[诊断](624/drafts/development_diagnostic.md)。没有把成熟Fisher身份或旧H²预算另计为新增定理。

本轮合并C03实际记录、C19状态矩、C22来源及C20报告压缩的部分接口。留下的主要问题不再是怎样额外设计一个记录器，而是：**同一尺度映射怎样同时保留记录、来源均值、噪声和迟致相位，并明确允许丢弃什么。** 接[625](625/drafts/STATUS.md)先回查592谱近似、612全时障碍与622虚跃迁结果，寻找共同的有限尺度匹配条件。原仪器的能源、控制、存储及引力约束问题仍按593—597保留，统一目标未完成。
