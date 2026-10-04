# 第660轮：原手征权重、局部配对与观测变换的共同表示

日期：2026-10-02。接[659](../../research_note_659.md)。[代码](../joint_local_mirror_process.py)、[结果](../joint_local_mirror_process_results.json)、[核验](../research_round_660_checks.json)、[条件账](../unified_physics_condition_ledger_660.md)。

## 1. 问题、历史与本轮结论

当前优先整合条件，认知系统设计后置。上一交流轮只确认执行顺序，没有新增科学结果。本轮接659真实边界连接及[660起步材料](physical_process_entry.md)，检验能否在同一个表示中保留原权重、物理时间反射与原观测。

**结果：**存在一个保留原overlap动能、仅含逐点配对的有限格候选，经过显式三角变量变换，精确还原原手征权重及指定保留观测。在单位规范背景下，其局部变量代数的反射正性可由成熟自由overlap定理与逐点相互作用延拓证明。但原左手bar观测一般必须变换；该变换含全局逆算子，尚不能把新表示的局部正性直接转给全部原物理观测。

与659不同，这一映射在指定匹配参数处是有限矩阵精确身份，不再取第五方向极限。它也不把659的体模型配分函数改写为同一对象：若使用659体表示，其体减除仍须保留。两种表示分别标明参考因子。

回查605硬允许域、612谱类型、613朴素手征筛选、615原holonomy、655原完整图过程、656—658空间辅助测度。不重做这些反例。2026报告在612—613已出现；本轮只补原全文和实际对象映射，不将其当作新发现。

## 2. 输入、约定及文献接口

保留原四维Wilson结构、m₀=1、内部16维T矩阵、B=iγ₅γ₂γ₄及单位S⁹场E。自由正性结论采用单位规范链路、周期空间、偶数反周期时间；数值空间例为两条平凡空间方向和一条非平凡方向，并未推出维数。另用615的带电平坦holonomy检验代数及来源，不对它额外宣称一般规范反射正性。

认知动机是“权重、态与观测须共同输送”；模型输入包括上述格、表示、球面处方，以及本轮选择共用E的两种配对。原独立bar-E的所有观测没有因此被保留。没有新增认知公理。

成熟输入采用[Kikukawa–Usui，自由overlap反射正性，§III、V](https://arxiv.org/html/1005.3751v3)：自由核具有正的跨反射分解，逐点且反射配对的相互作用可以纳入同一证明。本轮给下面具体Majorana配对的映射，并不声称其等于该文式51的模型。另核[Kikukawa §6.2—6.3](https://arxiv.org/html/1710.11618v3#S6.SS2)的投影配对及Schur结构；其作用中配对未带1/2，式206也显示动能／配对比的因子2。下文明确用指数+ξᵀNξ/2定义系数，不能把本轮y=1误读为该文相同数值耦合。

令n=64V，r=n/2=32V，固定手征帧J±满足J±J±†=Q±。原负／正H帧为u、v，取共同index-zero且Kₗ可逆的局部背景片：

$$
D=\frac{I+\gamma_5\operatorname{sign}H}{2},\qquad
Du=Q_-u,\quad Dv=Q_+v,\qquad
K=J_-^\dagger Du,\quad K_\ell=J_+^\dagger Dv.
\tag{1}
$$

M(E)逐点为B⊗T(E)，bar-M(E)为B⊗T(E)†，且T(E)ᵀ=T(E)、T(E)T(E)†=I。615代码用内部优先次序，须同步交换张量因子。候选定义为

$$
\xi=\binom{\psi}{\bar\psi^{\mathsf T}},\qquad
\mathcal N(E)=\begin{pmatrix}M(E)Q_+&-D^{\mathsf T}\\D&\bar M(E)Q_-\end{pmatrix},
\qquad \mathcal W=\exp\!\left(\tfrac12\xi^{\mathsf T}\mathcal N(E)\xi\right).
\tag{2}
$$

即+barψDψ、两项配对各系数1/2。仅配对严格逐点；D仍是指数局域而非有限邻接核。本文标题“局部配对”不把整个overlap作用称作超局域。

## 3. 原权重的精确连接

改写ξ=Sχ，χ=(b,c,d,e)，其中ψ=ub+vc，barψᵀ=J−*d+J+*e。记

$$
S=\operatorname{diag}([u,v],[J_-^*,J_+^*]),\quad
A_\pm=u^{\mathsf T}M Q_\pm u,\quad A=A_++A_-,\quad
\bar B=J_-^\dagger\bar M J_-^*,\quad C=u^{\mathsf T}M Q_+v,\quad F=v^{\mathsf T}M Q_+v.
\tag{3}
$$

原B的手征限制为单位反对称矩阵，T(E)单位，故bar-B可逆。固定原方向下Pf(bar-B)=1；由式(1)及原内部矩阵直接得到

$$
K^{\mathsf T}\bar B^{-1}K=A_-,\qquad
S^{\mathsf T}\mathcal N S=
\begin{pmatrix}
A_+&C&-K^{\mathsf T}&0\\
-C^{\mathsf T}&F&0&-K_\ell^{\mathsf T}\\
K&0&\bar B&0\\
0&K_\ell&0&0
\end{pmatrix}=N_c.
\tag{4}
$$

这里不是仅把镜像子块消元：C及F会影响左手bar观测，不能删去。采用行列式为1的三角变换η=Tχ：

$$
b'=b,\quad c'=c,\quad
d'=d+\bar B^{-1}Kb,\qquad
e'=e+K_\ell^{-\mathsf T}C^{\mathsf T}b-\tfrac12K_\ell^{-\mathsf T}Fc.
\tag{5}
$$

代入Grassmann二次型即可逐项消去b-d、b-c及c-c项：

$$
T^{-\mathsf T}N_cT^{-1}=
\begin{pmatrix}
A&0&0&0\\0&0&0&-K_\ell^{\mathsf T}\\
0&0&\bar B&0\\0&K_\ell&0&0
\end{pmatrix}=N_0.
\tag{6}
$$

r为32的倍数，重排符号与自由手征块Pf符号均为正。因此

$$
\operatorname{Pf}\mathcal N(E)=
\frac{\det K_\ell}{\det S}\operatorname{Pf}A(E),\qquad
Z_{\rm loc}=\frac{\det K_\ell}{\det S}\int d\mu(E)\operatorname{Pf}A(E).
\tag{7}
$$

这保留真实相位及帧Jacobian。A可能有零模；式(7)不要求A可逆。Kₗ不可逆处，权重身份可在同一index-zero帧片内由连续性延伸，含逆Kₗ的观测身份不能直接延伸，更不覆盖不同index的拓扑部门。

起步镜像子块的一般系数给yA++z²/y A−。恢复原A要求两系数共同匹配；y=z=1是本文明确选点。动能与两种配对不能各自独立任选。这是C01/C04/C16的一个实际交叉限制，尚不是认知原则唯一选出的耦合。

## 4. 必须同时输送观测与来源

对线性Grassmann来源jᵀξ，变换后来源为

$$
j_\eta=T^{-\mathsf T}S^{\mathsf T}j,\qquad
G_\eta=T S^{-1}G_\xi S^{-\mathsf T}T^{\mathsf T}=-N_0^{-1},
\quad G_\xi=-\mathcal N^{-1}.
\tag{8}
$$

符号按指数+二次型，G为有序二点Wick核；不用逆矩阵时可直接以有限Grassmann多项式定义来源泛函。保留的b、c、e′观测精确对应原辅助配对和自由左手Weyl部门；d′只给本候选的bar配对条件分布。原bar-E独立时的bar观测没有一并识别。

例如在A可逆的背景片，令Kₗ⁻ᵀCᵀ=L，−Kₗ⁻ᵀF/2=R，则e=e′−Lb−Rc。其裸e-e逆核一般非零，目标自由左手bar-bar核却为零。实际2×2自由盒的复算给

$$
\|G_{ee}\|_F\simeq8.63016854219,\qquad
\|G_{e'e'}\|_F<5.4\times10^{-15},\qquad
\|G_\eta+N_0^{-1}\|_{\max}<1.8\times10^{-15}.
\tag{9}
$$

这是固定E的条件核差，不据此声称球面平均后的无来源两bar期望也非零；对称性可能使后者消失。作为包含E观测的泛函，两者确有差别：在该非零差异的开邻域选择有界连续来源系数，与未归一差异的复共轭相乘，积分得到非零模平方。故相同权重仍不能替代实际观测映射。

来源若随背景λ变化，式(8)的T及S也必须随之求导；式(7)同时要求

$$
\partial_\lambda\log Z_{\rm loc}
=\partial_\lambda\log\det K_\ell-
\partial_\lambda\log\det S+
\partial_\lambda\log\int d\mu(E)\operatorname{Pf}A(E).
\tag{10}
$$

615原holonomy θ=.23的完整S⁹积分为0.100662735736129，真实log响应为−20.680469012824。代码对完整Nambu权重作球面径向精确求积再差分，误差1.50×10⁻⁸。若仅将u和v的首列分别重定相位7θ和11θ，却漏掉detS，便产生纯假的18i来源；代码复现。这将帧选择、原测度及C22共同来源真正接在一起。

## 5. 自由背景的物理时间反射正性

采用标准反线性、逆Grassmann乘积次序的Θ：ψ映到(barψγ₄)ᵀ，barψ映到(γ₄ψ)ᵀ，实E映到反射位置的E。记Γ为时反射置换乘γ₄及内部单位，RΘ为其Nambu形式，原B/T矩阵给

$$
R_\Theta=\begin{pmatrix}0&\Gamma^{\mathsf T}\\\Gamma&0\end{pmatrix},\qquad
\mathcal N(E)=-R_\Theta^{\mathsf T}\mathcal N(\theta E)^*R_\Theta.
\tag{11}
$$

式(11)只表示对称性，本身不是正性证明。实际证明利用配对严格逐点且偶，故它的作用V=V++ΘV+，跨时间半区的部分完全来自成熟自由D。将独立正S⁹乘积测度加入自由反射泛函，仍为正；再令F+为正半区局部多项式，有

$$
\mathcal I_{\rm loc}(\Theta F_+F_+)
=\mathcal I_{\rm free\otimes\mu}
\left[\Theta(e^{V_+}F_+)\,e^{V_+}F_+\right]\ge0.
\tag{12}
$$

自由泛函与独立两半E乘积在右式分别因子化为正型核；有限Grassmann次数及紧S⁹保证所有系数可积。也可直接使用成熟跨反射正锥分解。这里是对本候选的解析延拓，数值只核式(11)和V的半区分解，没有用抽样特征值代替证明。

归一也须严格检查。取常E₀，将两配对同时乘λ≥0，λ>0时Schur配对为

$$
A_\lambda=u^{\mathsf T}M(E_0)
\left(\lambda Q_++\lambda^{-1}Q_-\right)u.
\tag{13}
$$

自由B反单位对称性保持原u空间，故它等于一个可逆单位配对乘压缩的正矩阵，再乘内部T(E₀)；无零模。λ=0为无零模的自由D，AP偶数时间使各时间sin分量非零。常E₀反射不变，所以Pf沿λ为非零实数，其在0处为正detD，进而Pf N(E₀)>0。由式(7)及658的原归一球面平均严格正：

$$
Z_{\rm loc}=\operatorname{Pf}\mathcal N(E_0)
\int d\mu(E)\frac{\operatorname{Pf}A(E)}{\operatorname{Pf}A(E_0)}>0.
\tag{14}
$$

这样可归一式(12)。没有声称原权重对每个E都正，也没有为一般规范背景给正性证明。[2026报告](https://indico.yukawa.kyoto-u.ac.jp/event/78/contributions/1982/)的相关正性结论同样限定弱规范极限；本轮保存[原PDF及文本来源](source_receipt.json)，不读图、不凭缺字补公式。

## 6. 为什么还没有得到原共同物理过程

式(12)控制的是局部ξ、E在正半区的代数。式(5)使用Kₗ逆和整体手征帧，可能把正半区原观测映成跨两半的变量。要把正性转给指定原观测代数，还须实际满足

$$
\Phi(\mathcal A^{\rm original}_+)\subseteq\mathcal A^{\rm loc}_+,
\qquad \Phi\Theta_{\rm original}=\Theta_{\rm loc}\Phi,
\qquad \Phi\circ\tau_{\rm original}=\tau_{\rm loc}\circ\Phi.
\tag{15}
$$

前两项为时间支撑和反射兼容，最后一项为时间平移兼容。本轮精确代数变换没有证明它们。尤其不能把一个有限盒反射泛函等同655原固定图Hamiltonian，也不能越过612对原固定有限图全时间强合同的谱限制。

整个依赖链是

$$
\underbrace{D,B,T,S^9,\text{匹配配对}}_{\text{模型输入}}
\ \Longrightarrow\ \underbrace{(7),(8),(10)}_{\text{精确权重／观测／来源映射}},
\qquad
\underbrace{\text{自由RP},\ V_++\Theta V_+}_{\text{既有定理及本轮核对}}
\ \Longrightarrow\ \underbrace{(12),(14)}_{\text{新表示的归一正泛函}}.
\tag{16}
$$

两条链目前由式(15)这个明确条件隔开。共同连续／无限时间极限、一般规范测度、动态量子几何、维数和群的选择仍独立开放。

## 7. 核验与下一步

四组实质检查：完整带相位Nambu消元（512维空间例及两种128维带电背景）；二点与四点来源输送及裸bar反例；两个偶数时间长度的实际反射／逐点分解；615全球面平均与帧相位来源。结果全部通过。代码复用既有Python和NumPy；起步probe重复核对但不计新组；无独立代理、图像检验或目标改动。

本轮同时压缩C01/C04/C16/C19/C22：原权重、有限自由正表示、保留观测与来源可以共同规定，匹配系数和帧Jacobian不是独立自由选择。没有关闭C03/C09的实际记录、C11量子引力约束、C20共同物理极限。

接[661观测代数与共同物理时间接口](../../661/drafts/STATUS.md)：先检验成熟Weyl正性所用的实际局部观测是否能接到本轮保留观测，以及失败是否只属于直接映射。只做到能判断接口相容性，不展开新存储或控制器设计；若必须保留额外条件，写回总账并返回其它可整合接口。
