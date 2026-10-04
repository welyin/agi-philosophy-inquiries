# 第608轮：投影补全的独立组合、局域性与共同过程

日期：2026-10-01。接[607](research_note_607.md)及[608入口](608/drafts/STATUS.md)。[代码](608/joint_projection_local_composition.py)、[结果](608/joint_projection_local_composition_results.json)、[核验](608/research_round_608_checks.json)、[条件账](608/unified_physics_condition_ledger_608.md)。三组复算、十二式；主代理审查，无新增独立代理审查。

## 1. 新问题和结论

607在给定投影下接通了约束演化、完整历史、热态与几何来源，但使用“整个系统的允许部门P／其余全部Q”两块补全。本轮检验**同一规则能否与独立组合及局域传播共同成立**，不另造记忆装置。

结果：全局两块规则一般不保独立组合，甚至对严格在位的约束与原相互独立的动力学，也能引入距离无关的影响。这个反例针对607规则在整个P＋Q空间上的普遍局域性，不是针对只允许P内态和操作的物理理论。

若实际约束是独立在位投影，可改为逐个约束分块。这样不扩大有界相互作用的支持，保原相互作用范数上界，并与独立张量组合相容；在共同允许部门，607的演化、完整记录、任意参考和来源全部保留。代价是放弃“保留整块QHQ”的额外要求。真实GW投影一般不是在位投影乘积，因此本轮没有签收其局域性。

|层次|内容|
|---|---|
|认知继承|早期独立组合与230完整有限过程；比较同一对象的组合和演化|
|建模输入|有限在位因子、独立在位正交投影、有界局部相互作用；指定保部门仪器|
|解析结果|全局补全的组合缺陷及定时远程影响；逐约束补全的支持、范数及共同过程恒等式|
|数值验证|非平凡二维允许部门的qutrit链、带相互作用的逻辑部门、反馈记录与同一来源|
|不能外推|当前Gauss区域自动因子化、GW投影在位分解、无界电动能的LR界、维数或GR|

去重复：[360](../archive_342_369/research_note_360.md)已给Gauss区域非因子化，[363](../archive_342_369/research_note_363.md)已给交织泄漏，[466](../archive_429_466/research_note_466.md)已强调逐原子不变部门；本轮没有重做这些结论。成熟[Nachtergaele—Ogata—Sims传播界](https://arxiv.org/abs/math-ph/0603064)用于有界局部相互作用的后续推论；新增工作是审计607补全是否满足这些前提。[Hernández—Jansen—Lüscher](https://arxiv.org/abs/hep-lat/9808010)的光滑规范背景下指数局部性针对一粒子核，不能替代整Fock空间补全的组合审计。

## 2. 全局两块补全一般不保独立组合

记607有限矩阵版本为

$$
\mathcal C_P(H)=PHP+(I-P)H(I-P).
\tag{1}
$$

两个独立系统取H=H_A⊗I＋I⊗H_B、P=P_A⊗P_B，定义O_A=H_A−C_{P_A}(H_A)，B同理。直接按两个约束的四个联合部门展开：

$$
\mathcal C_{P_A\otimes P_B}(H)
-\left[\mathcal C_{P_A}(H_A)\otimes I+
I\otimes\mathcal C_{P_B}(H_B)\right]
=O_A\otimes Q_B+Q_A\otimes O_B .
\tag{2}
$$

当B违反自己的约束时，A在其允许／不允许部门之间的跃迁被保留；B满足约束时，同一跃迁被删除。虽然原H_A完全不依赖B，新规则却依赖B的状态。这是整个Q块内部仍保原跃迁所导致，与约束是否严格在位、原H是否有长程作用无关。

## 3. 非平凡允许部门中的明确诊断

每点取三维空间，两个允许态|0〉、|1〉和一个补态|2〉。按一条长度可增大的链赋予距离，定义

$$
p=|0\rangle\langle0|+|1\rangle\langle1|,\quad
X=|0\rangle\langle2|+|2\rangle\langle0|,\quad
Z=|0\rangle\langle0|-|2\rangle\langle2|,\quad
P=\bigotimes_{i=0}^{N-1}p_i,\quad H=NI-X_0 .
\tag{3}
$$

原H正且无节点间作用；允许部门维数2^N，不是把整个物理空间压成一维。全局补全为

$$
H^\sharp=NI-X_0(1-R),\qquad
R=\prod_{k\ne0}p_k .
\tag{4}
$$

任意j≠0有

$$
\|[[H,Z_0],X_j]\|=0,\qquad
\|[[H^\sharp,Z_0],X_j]\|=2 .
\tag{5}
$$

还可直接给固定时间的实际接收差，避免单靠短时导数推断。初态为|0…0〉，j端可选择I或交换|0〉、|2〉且固定|1〉的局部酉。其余端点均不操作。无翻转时R=1，0端保持|0〉；翻转后R=0，0端变为cos(t)|0〉＋i sin(t)|2〉。因此

$$
\frac12\|\rho_0^{(0)}(t)-\rho_0^{(1)}(t)\|_1=|\sin t| ,
\qquad j\ne0 ,
\tag{6}
$$

对任意N、任意距离均成立。t=.275时为.271546936956。故全P＋Q算符／操作合同下，不能有固定t、随距离趋零且与规模无关的传播界。

**范围必须保留：** j端翻转把态带出P。若模型只允许P内态及保P操作，这不是它的超距通信协议。607同时讨论了整个补全系统的Gauss热态与P部门映射；两者不能混成同一个权限范围。本反例要求放弃的是从“分块、正性、规范相容”直接推“整个补全局域”的推断，不是证明所有投影模型非局域，更不是统一研究的反证。

## 4. 在位约束的共同修正

若实际每个p_i只作用一个独立因子，取E_i=C_{p_i}。这些映射两两对易，令

$$
\mathcal E=\prod_i\mathcal E_i,\qquad
H_{\rm loc}=\mathcal E(H)
=\sum_{\boldsymbol s}\Pi_{\boldsymbol s}H\Pi_{\boldsymbol s},
\quad
\Pi_{\boldsymbol s}=\bigotimes_i
\begin{cases}p_i&s_i=0\\I-p_i&s_i=1.\end{cases}
\tag{7}
$$

这里保留每一处约束是否满足的信息，禁止不同联合部门间的跃迁。给定这些部门和各自原形式，该选择被确定；并非唯一可能的局域物理定律。

若H=Σ_S h_S、h_S仅作用S，S外的E_i作用为恒等。每个E_i是正、保单位的分块映射，算符范数收缩。因此

$$
H_{\rm loc}=\sum_S\mathcal E_S(h_S),\qquad
\operatorname{supp}\mathcal E_S(h_S)\subseteq S,\qquad
\|\mathcal E_S(h_S)\|\le\|h_S\| .
\tag{8}
$$

原有界相互作用的范围、加权相互作用范数和规模一致LR充分条件由此继承；无须另选传播常数。对独立系统，E_A⊗E_B作用于张量和，严格得到各自补全的张量和。对严格在位原H，本轮修正仍严格在位。

这项证明不能无条件用于重叠约束、非对易约束、准局域谱投影或无界场。特别是Gauss区域的独立因子合同须先按360的边界规则核实，不能由形式上的张量写法宣布已经成立。

## 5. 已经接通的共同过程没有丢失

令P=∏p_i，则P与Hloc、H♯都对易，且

$$
PH_{\rm loc}P=PHP=PH^\sharp P,\qquad
H_{\rm loc}J=Jh=H^\sharp J,\quad h=J^\dagger HJ .
\tag{9}
$$

所以两种补全在P内有相同等待过程；只要每个指定L_r与所有p_i对易，便保P。对任意有限自适应历史与任意被动参考R，

$$
K_{\boldsymbol r}^{\rm loc}J
=K_{\boldsymbol r}^{\sharp}J=Jk_{\boldsymbol r},
\qquad
(K_{\boldsymbol r}^{\rm loc}\otimes I_R)J_R
=(K_{\boldsymbol r}^{\sharp}\otimes I_R)J_R .
\tag{10}
$$

它保每个未归一化结果后态，故同时保成功率、全部记录和参考关联。不是只比较某个概率。补空间动力学、全空间热态和谱一般改变；P内受限Gibbs态则相同。没有保证全部CP权限或免费准备。

若各p_i和J不依赖所微分的几何参数γ，原来源G=∂γH，则

$$
G_{\rm loc}=\mathcal E(G),\qquad
J^\dagger G_{\rm loc}J=\partial_\gamma h
=J^\dagger G^\sharp J .
\tag{11}
$$

同样，由L_r与p_i对易，仪器的Heisenberg映射I*与E对易，因而

$$
\mathcal I^*(H_{\rm loc})-H_{\rm loc}
=\mathcal E\!\left(\mathcal I^*(H)-H\right),\qquad
\mathcal I^*(G_{\rm loc})-G_{\rm loc}
=\mathcal E\!\left(\mathcal I^*(G)-G\right).
\tag{12}
$$

这把区域组合、局域性、P内完整过程及同一来源合在同一个规则下；没有各自另拟合能源或记录。代价明确：放弃保整块QHQ，增加实际独立在位约束的前提。当前GW场依赖子空间不自动满足它们。

## 6. 复算和验收范围

三组均通过：

1. N=2、3、4的(3)—(6)，距离增大时接收差严格为同一值；逐约束补全的对应双对易子为零。全局补全保QHQ，局部补全确实改变QHQ。
2. 加入允许部门内的非对角演化、近邻耦合及带间项；逐项核支持、范数、独立组合和(9)。P部门非平凡且可交互；全空间两种热配分函数不同，不误报热态等价。
3. N=3、κ=e^(−2γ)、γ=.12，两次有反馈的四分支仪器历史，二维纠缠参考。历史算符误差5.93×10⁻¹⁵，带记录联合态迹范数差6.61×10⁻¹⁵；能源误差1.29×10⁻¹⁴，来源差分误差9.36×10⁻¹¹，仪器注能交换图残差为零。

qutrit诊断检验607规则的逻辑性质，不冒充标准模型物种、真实Gauss准备或当前GW投影。解析的任意N与任意参考结论由上述恒等式承担；小矩阵测试只是复算。

**下一步：[609](609/drafts/STATUS.md)返回606实际一粒子手征投影。** 核能否先构造一粒子协变导数，再提升到Fock空间，并以原Φ项保P部门能源；联查质量项、局部系数和同一来源。不能把本轮在位乘积分解假定到真实GW核，也不能只给系数衰减就宣布无界电动能满足LR界。认知装置设计继续后置，统一目标未结项。
