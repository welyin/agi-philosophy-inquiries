# 第617轮：原规范区域拼接、完整演化与共同能源来源

日期：2026-10-01。接[616](research_note_616.md)、[617入口](round617_drafts/STATUS.md)及[去重记录](round617_drafts/region_source_entry.md)。[代码](joint_region_energy_gluing.py)、[结果](joint_region_energy_gluing_results.json)、[核验](research_round_617_checks.json)、[条件账](unified_physics_condition_ledger_617.md)。三组复算、十二式；主代理审查，无新增独立代理审查。

## 1. 问题、结论与范围

继续先合并成立条件，认知系统的实现设计后置。本轮将**原规范群、区域匹配、完整Hamiltonian、物理态、记录和几何来源**接在同一个映射上。

**条件性结论：** 对589／598／603既有有限图模型，形式上切开原链路，再施加切口Gauss匹配，可以与原完整物理过程精确等价；同一映射保持原费米模式、正常热态和几何变分。原电动能必须按同一系数矩阵分配，包括非对角剪切项。若两边各复制一份原动能，虽有正确边界匹配，仍把能量及其来源加倍。

这减少了“区域拼接正确、动力学正确、来源正确”各自独立挑选处方的自由。**没有**把整体Gauss空间变成各区域物理空间的普通张量积，也没有建立新物理节点、真实空间细化或新的认知装置。

## 2. 历史、成熟结果与额外输入

- [360](research_note_360.md)、[362](research_note_362.md)、[365](research_note_365.md)已经给共享中心、独立准备障碍与Z₂边界协议；本轮不重做这些结论。
- [574](research_note_574.md)、[589](research_note_589.md)给原配置空间和完整正空间度规的闭形式；[598](research_note_598.md)给原有限CAR和Gauss质量；[603](research_note_603.md)给固定图热态。
- [608](research_note_608.md)已经区分整体编码与独立组合；[616](research_note_616.md)确定普通热迹不能换成双宇称周期超迹。

采用成熟的规范区域扩展与匹配构造，而非宣称发现新的因子化理论：[Donnelly，2012](https://arxiv.org/pdf/1109.0036)讨论切开边、表示指标及边界熵；[Donnelly—Freidel，2016，§§2.5、4.2](https://arxiv.org/html/1601.04744v2)区分区域扩展与边界匹配后的整体空间。本轮自己的工作是核其与本项目原群、**非对角动能**、CAR及几何来源的共同字典。上述文献中的引力构造不作为本项目已推出GR的依据。

明确保留：固定有限图、原三维邻接处方、严格正背景几何、原紧商群和物种、耦合与量子处方、选择的态。切口只改变表示；不改变格距、场数、原顶点体积或原链路几何权重。J在旧背景Hilbert空间固定识别下与几何参数无关。

## 3. 原商群上的精确匹配空间

采用原紧群及归一Haar测度：

$$
G=(SU(3)\times SU(2)\times U(1)_Q)/\mathbb Z_6,\qquad
J:L^2(G)\longrightarrow L^2(G\times G),\qquad
(Jf)(A,B)=f(AB).
\tag{1}
$$

这是直接定义在商群上的群乘法，不把商群换成直积。表示中的整数Q约定及所有原物种不变。切口作用为(A,B)变成(Ah⁻¹,hB)。Haar不变性给

$$
\|Jf\|^2=\int dA\,dB\,|f(AB)|^2=\|f\|^2,\qquad
\operatorname{Ran}J=L^2(G\times G)^{G_{\rm cut}},\qquad
JJ^\dagger=P_{\rm cut}.
\tag{2}
$$

满射到不变子空间也须证明：作测度保持换元(U,A)=(AB,A)，切口作用只右移A。紧群上右平移不变的L²函数关于A几乎处处常数，故只依赖U。P_cut是切口作用的Haar平均。

对商群允许的任意不可约酉表示R，使用归一基函数√d_R R(U)_ij，乘法法则给

$$
J|R,i,j\rangle
=\frac1{\sqrt{d_R}}\sum_{k=1}^{d_R}
 |R,i,k\rangle_A\otimes|R,k,j\rangle_B .
\tag{3}
$$

三个或更多切段由群乘法结合律保证匹配一致，n段系数为d_R^{-(n-1)/2}。表示指标k是匹配索引，不是新增的可独立准备物质。每条切边应用J，其他变量和原CAR空间应用身份，得到全图等距J。原端点规范变换与J交织，所以原顶点Gauss限制与新增切口限制可同时施加，且匹配后物理空间与原物理空间酉等价。

## 4. 非对角电动能的共同分配

对每个简单因子取正交Hermitian生成元t_a，U(1)按原整数Q归一。令P_a为−iℏ乘以沿exp(it t_a)U的微分。第一段使用同方向左乘；第二段必须先用A把其方向运回原源点。定义O(A)满足A t_b A⁻¹=Σ_a O_ab(A)t_a。则

$$
P^A_a=-i\hbar X^A_a,\qquad
\widetilde P^B_a=\sum_b O_{ab}(A)(-i\hbar X^B_b),
\qquad
P^A_aJ=JP_a=\widetilde P^B_aJ .
\tag{4}
$$

证明第二个等式：A exp(it t_b)B=exp(it Ad_A t_b)AB，再用O的正交性。O依赖A，B方向微分不作用于O；两个微分算符均关于乘积Haar形式对称，并保切口匹配空间。第一段方向与第二段的裸方向不能在非Abel情形直接认作相同。U(1)的O为1。

继承589：在同一源点i，三条出边以及同一群因子的方向μ、ν共同出现，

$$
q_E[f]=\sum_{i,c,a,\mu,\nu}
 \langle P_{i\mu,ca}f,K^{(c)}_{i,\mu\nu}P_{i\nu,ca}f\rangle,\qquad
K^{(c)}_{i,\mu\nu}
=\frac{b_c}{\epsilon}\,
\frac{\bar\gamma_{i,\mu\nu}}{\psi_{e_\mu}\psi_{e_\nu}},
\qquad K^{(c)}_i>0 .
\tag{5}
$$

原几何为γ=ψ⁴barγ，det barγ=1。不能只保存每条边的Casimir而删掉μ≠ν项。选择任意平滑的半正定矩阵分配K_A+K_B=K；例如K_A=K^{1/2} L K^{1/2}，0≤L≤I，K_B=K−K_A。切口动能形式使用

$$
q_{E,\rm cut}[\Psi]
=\sum\langle P^A\Psi,K_A P^A\Psi\rangle
 +\sum\langle\widetilde P^B\Psi,K_B\widetilde P^B\Psi\rangle,
\qquad q_{E,\rm cut}[Jf]=q_E[f].
\tag{6}
$$

求和包含式(5)的全部指标；未切边的两个提升都取原P。正分配只是同一匹配空间上同一形式的表示自由，不是新的物理参数。若L或K_A、K_B随几何变动，求导时必须连同它们一起微分；总和身份仍保证原来源。没有把两段设成新长度后重新套格距处方。

这里B提升含A，完整跨区域作用仍保留，故式(6)不表示两个可独立执行的区域Hamiltonian。对任意切图，最简选择K_A=K、K_B=0也足够表达相同匹配过程。

## 5. 原物质、相互作用及算符域

所有原乘法算符按U→AB拉回；其他变量不变。特别是原带荷跳跃

$$
c_v^\dagger\,t_e(\gamma)R(U_e)c_w+\mathrm{h.c.}
\ \longmapsto\
c_v^\dagger\,t_e(\gamma)R(A_e)R(B_e)c_w+\mathrm{h.c.}
\tag{7}
$$

中仍是同一t_e、同一端点费米算符。原Higgs的R_H(U)=z³W、磁圈的原乘积、Yukawa／Majorana与F依赖全部保留。跨边项没有变成两个互不关联的局部项。

CAR空间没有增加切口费米子。使用原固定次序或等价的分次张量积，J只作用于规范变量，与总费米宇称交换；不能把跨区域的奇算符当成普通交换变量。

完整闭形式和Hamiltonian不能仅靠逐项符号一致来宣称。记原Gauss物理空间上的下有界闭形式为q_γ，定义匹配空间上的形式及其域为

$$
D(q_{\gamma,\rm cut})=J D(q_\gamma),\qquad
q_{\gamma,\rm cut}[Jf,Jg]=q_\gamma[f,g].
\tag{8}
$$

J是到匹配空间的酉映射，因此闭性、下界和稠密性全部继承。式(4)—(7)在原公共光滑核上给出其显式局部表达，闭包取式(8)。不需要证明未约束的整个扩展空间有良好热迹，也不允许把它的冗余态加入原配分函数。由闭形式表示定理，

$$
H_{\gamma,\rm cut}=JH_\gamma J^\dagger
\quad\hbox{在匹配物理空间上},\qquad
D(H_{\gamma,\rm cut})=J D(H_\gamma),\qquad
e^{-itH_{\gamma,\rm cut}}J=Je^{-itH_\gamma}.
\tag{9}
$$

这包含原无界标量／规范动能和有限CAR矩阵势；不是只在数值截断成立的声明。598的相对形式界和603的固定图热迹仍按原前提继承。几何若另行量子化，或J本身依赖几何，则本证明需要扩充，不能漏掉联络项。

## 6. 来源、正常热态与完整记录

在旧固定Hilbert识别、共同形式域及可微背景族范围内，J与γ无关，因此形式导数满足

$$
(\partial_\lambda q_{\gamma,\rm cut})[Jf,Jg]
=(\partial_\lambda q_\gamma)[f,g],\qquad
J^\dagger(\partial_\lambda H_{\gamma,\rm cut})J
=\partial_\lambda H_\gamma
\quad\hbox{按相应形式意义}.
\tag{10}
$$

不能把第二式误读为任意态上均有有界应力算符；期望与响应沿用603已核的能源域。式(10)保留同一几何输入，不证明它已满足量子Einstein约束。

原正常热态ρ_β=e^{-βH}/Z满足

$$
Z_{\beta,\rm cut}=Z_\beta,\qquad
\rho_{\beta,\rm cut}=J\rho_\beta J^\dagger,\qquad
\operatorname{Tr}(\rho_{\beta,\rm cut}\,\partial_\lambda H_{\rm cut})
=\operatorname{Tr}(\rho_\beta\,\partial_\lambda H).
\tag{11}
$$

迹取匹配后的Gauss物理空间。616的普通迹与固定宇称分支各自照原定义输送，不能在切口偷偷换成超迹。任何已有、保物理空间的仪器K_r输送为J K_r J†，其多时结果概率、条件态和与任意参考系统的关联均保持；这是完整过程等价，不是独立地区任意操作权限的证明。

## 7. 两种具体错误接法与复算

若两个提升各配一份完整K，式(4)立即给

$$
q^{\rm duplicate}_{E,\rm cut}[Jf]=2q_E[f],\qquad
\partial_\lambda q^{\rm duplicate}_{E,\rm cut}[Jf]
=2\,\partial_\lambda q_E[f].
\tag{12}
$$

其差异不是切口Gauss投影能消除的。若改为只留K的对角部分，则一般删去原剪切耦合。两项是明确错误处方的排除，不是一般区域构造不可能。

代码三组实际复算：

|检验|对象与结果|
|---|---|
|原商群匹配|原Q、u、d、L、e、ν表示维数6、3、3、2、1、1；商群误差8.62×10⁻¹⁶，切口非Abel匹配误差1.66×10⁻¹⁵；三段结合与参考保持通过|
|非对角动能／来源|原冻结b_w=0.2101856811；三条弱链路的4096×8切分映射；动能差3.34×10⁻¹⁶，剪切／共形导数误差小于3.50×10⁻¹¹；删除剪切的算符差0.101476，复制动能差0.565006，共形来源差1.130012|
|原物质乘法／CAR|非交换原群链路、原Higgs距离、原磁势及4模式端点跳跃；最大误差1.12×10⁻¹⁵；丢掉半边链路产生3.363913的跳跃算符差，新增切口物质模式为0|

第二组为三个出边的原弱Peter–Weyl系数测试，固定右指标，**在端点Gauss限制之前**；对应有限扇区热来源0.04733234008，差分误差小于3.47×10⁻¹¹。第三组用原L表示的一个固定旋量分量测4个端点CAR模式，不计算全图费米谱。全空间、物理Gauss与正常热态的继承依靠第3—6节证明，不由这些小矩阵代替。

## 8. 哪些条件真正合并，哪些仍开放

- C01／C02／C14：原量子对象、同一商群及匹配后的区域组合共用J，不用增加切口可独立物理态。
- C03／C19／C22：原记录、正常热态及来源由同一完整过程映射共同保留，不能另外选能源账。
- C10的一部分：给定同一几何时，匹配来源不依赖这次形式切分；尚未证明所有物种自动共有物理度规。
- C20：只获得表示切分的精确一致性，**真实细化、跨尺度有效匹配仍开放**。
- C13：式(3)可出现标准的表示边界熵，但其本身不决定物理面积、面积系数或Einstein方程。本轮不重算旧面积熵路线。

同一整体的不同区域表示已经相容；独立主体怎样形成整体、全手征局域过程、实际认知操作、物理连续、维数选择、动态几何与完整GR约束仍未完成。把实现设计后置，没有删除这些最终验收要求。

下一步见[618入口](round618_drafts/STATUS.md)：转到给定引力分支的区域拼接与共同变分，回查旧共同作用和约束结果，先辨认物质匹配与几何边界需要共享哪些数据。不得把规范边切分直接当作引力约束闭合，亦不为模拟切口再设计认知硬件。

