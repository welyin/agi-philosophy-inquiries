# 第648轮：原物质参考、引力角点与区域量子拼接的共同条件

日期：2026-10-02。接[647](research_note_647.md)，复用[573](../archive_554_584/research_note_573.md)、[617](../archive_585_628/research_note_617.md)、[618](../archive_585_628/research_note_618.md)。[代码](648/joint_relational_corner_gluing.py)、[结果](648/joint_relational_corner_gluing_results.json)、[核验](648/research_round_648_checks.json)、[条件账](648/unified_physics_condition_ledger_648.md)、[起步去重](648/drafts/joint_interface_entry_audit.md)、[范围审计](648/drafts/corner_gluing_scope.md)。三组检查、十四式；主代理审查，无新增独立代理审查。

## 1. 本轮合并什么

647完成经典关系定位的来源，但把固定图量子区域接到动态几何还需要边界约束。617的内部紧群切口不能直接代表引力切口；618明确未覆盖角点。这里补这个真实缺口，而不再做554的联合参考精度或606的一般编码导数。

**结果：** 原F同时固定引力角点的boost荷、面积密度和标量来源；573同一物质参考确实可在同一解上确定一个局部角点及其法向夹角。另一方面，把617的普通Hilbert空间切口构造直接照搬到非紧boost参考的指定正则表示，得不到非零正常不变态；有限规范区间的近似等距映射也没有强极限。这限定了一个明确接法，不是对量子引力或全部区域组合的否定。

|层次|本轮地位|
|---|---|
|认知动机|整体切成区域后，内部规范、物质参考和几何边界必须共同匹配|
|继承输入|原3+1、正F、给定二导数Einstein—物质作用及573初值|
|成熟工具|非类光角点作用、Noether荷、非紧群正则表示|
|解析连接|原F的角点来源及框架一致性；指定boost拼接的正常态障碍|
|数值|原573物质／几何上的角点；五场曲率缩并；明确正则表示的积分核验|
|没有完成|引力量子Hilbert空间、关系instrument、统计面积律、全尺度映射|

## 2. 同一原作用固定角点系数

沿618的记号，五实场φ、Jordan度规和Einstein度规满足

$$
F=2-\phi^2/6>0,\qquad g_E=F g_J,\qquad
S_J^{\rm grav}=\tfrac12\int\sqrt{-g_J}FR_J
 +\sum_i\epsilon_i\int_{\Sigma_i}\sqrt{|h_J|}F K_J .
\tag{1}
$$

取非类光边界相交的类空二维角点S，不跨越法向因果类型变化，固定整体取向与角度加法约定。下面选择角点符号为正；反向取向同时反号。由Einstein框架标准角点项及二维诱导面积的Weyl权重，得到

$$
I_S=\int_S c\,\eta,\qquad
c=F\sqrt{q_J}=\sqrt{q_E},\qquad \eta_E=\eta_J .
\tag{2}
$$

η是法向平面的有向boost参数。这里固定使用618的Dirichlet作用，不另加独立角点反项；角度原点、取向或额外边界作用若改变，需要一起记账。[Jubb—Samuel—Sorkin—Surya，式(30)](https://arxiv.org/html/1612.00149)给同一R/2归一的面积乘角度项。本轮的内容是其与原五场F及物质参考的映射，不把一般角点项算作新定理。

还可从原曲率项独立验证。以单位反对称法向二形式ε_ab满足ε_ab ε^ab=−2，去掉体积密度的曲率导数为

$$
E^{abcd}=\frac{\partial(FR/2)}{\partial R_{abcd}}
=\frac F4(g^{ac}g^{bd}-g^{ad}g^{bc}),\quad
E^{abcd}\epsilon_{ab}\epsilon_{cd}=-F,\quad
-2\pi\sqrt{q_J}E^{abcd}\epsilon_{ab}\epsilon_{cd}=2\pi c .
\tag{3}
$$

对在S上为零、但法向一阶导数给定单位boost的向量场，Noether荷的系数为c；含ξ∂F的项在S消失。原Yang–Mills对这种ξ的项也因ξ=0消失，内部规范变换的电荷仍须另行匹配。当前原经典解的费米背景为零，不据此签收所有旋量角点及挠率问题。

[Iyer—Wald的Noether荷框架](https://arxiv.org/pdf/gr-qc/9403028)在合适的平稳黑洞条件下赋予式(3)熵解释。这里S只是普通局部角点，式(3)是几何荷密度；既不是热熵，也未证明637—644的量子区域熵等于它。

## 3. 标量来源与面积匹配不能独立指定

固定嵌入的角点数据共同变分给

$$
\delta I_S=\int_S\left[c\,\delta\eta+
\eta\sqrt{q_J}\left(-\frac{\phi_A}{3}\delta\phi_A
+\frac F2q_J^{ij}\delta q_{J,ij}\right)\right].
\tag{4}
$$

因此角点标量来源为−η√q_J φ_A/3。改变框架时，δq_E=Fδq_J+q_JδF，所以同一项正好进入ηδ√q_E，不能删去。若角点由物质标签定义而随场变化，还要加647已经给出的嵌入输送；式(4)的固定嵌入部分并不取代它。

两侧拼接至少须在共同法向参考和相反边界取向下匹配几何荷：

$$
Q_L[\zeta]+Q_R[\zeta]=0,\qquad
Q_{L,R}[\zeta]=\pm\int_S\zeta\,c_{L,R},\qquad
c_L=c_R\quad\text{在相同角点识别下} .
\tag{5}
$$

它不是充分拼接条件；诱导数据、体约束和切向对称性仍要匹配。原内部Gauss匹配没有约束q或F，故不能单独代替(5)。这增加的是给定引力作用要求的接口，不是新认知公理。

## 4. 573原物质参考提供实际局部角点

在573的类时h片，取两个标签面h=H₀和k=h+λs=K₀，λ为固定的边界选择，未进入动力学。λ≠0时其交集等价于固定h、s。令

$$
A=(dh)^2,\quad B=(ds)^2,\quad C=dh\cdot ds,\quad
D=A+2\lambda C+\lambda^2B<0,\qquad
\sinh\eta=\frac{\lambda\sqrt{C^2-AB}}{\sqrt{(-A)(-D)}} .
\tag{6}
$$

要求A<0、C²−AB>0且两个时间法向同向，方向符号固定。C来自同一字段，不是新独立钟尺；所有缩并共同换框架后比值不变。由原h类时及dh、ds独立，足够小λ确实给这样的局部角点，不能推广为全局图册。

在原点(x,y,z)=(0,π/2,π/4)，h的空间梯度为零，ds只有y分量，法向速度必须使用573的完整曲目标逆度量。原48³复算为

$$
F=1.86897167022480,\quad \psi=1.01826915188282,\quad
v_h=0.00777915012825,\quad v_s=-0.23867620330440,
\quad |\nabla s|_\gamma=0.03158713837300 .
\tag{7}
$$

选择λ=.005，仅作为同一解上的边界诊断。由法向分量直接算出的角度与(6)的Gram缩并一致：

$$
\eta=\operatorname{artanh}\frac{\lambda|\nabla s|_\gamma}{v_h+\lambda v_s}
=0.02398595763107,\quad
c=\psi^4=1.07510368057955,\quad
\sqrt{q_J}=c/F=0.57523808290269 .
\tag{8}
$$

该点交集切向基为∂x、∂z，故上述面积密度不是任意换一张平面所得。物质／几何仍是原573受约束解；代码中独立的变分诊断可以离壳，不冒称新变分都解了全部约束。

若用裸常数M=2乘Jordan面积替代同一原F，或在固定Jordan面积下忽略标量来源，原点就有

$$
\frac{M\sqrt{q_J}-c}{c}=0.07010717811439,\qquad
\delta c\big|_{q_J}=\sqrt{q_J}\delta F=-0.00728496544259
\quad(\delta h=.07,\ \delta s=-.02).
\tag{9}
$$

这排除“保持原体作用，却独立指定常系数角点面积荷”的接法。它不否定同时改变体作用的另一模型。

## 5. 617的紧群Hilbert拼接不能直接照搬

[Donnelly—Freidel，§3—4](https://arxiv.org/html/1601.04744v2)给引力表面对称性及其区域约束框架；其文并未完成引力Hilbert空间构造。这里只用固定法向几何后的单一boost子群R，不借SL(2,R)与372的SL(2,C)名字相似宣称同构。

617的证明依赖内部紧群G的归一Haar体积。对比明确的候选

$$
\mathcal H_G=L^2(G),\quad (J_G f)(A,B)=f(AB),\quad\int_GdA=1;
\qquad
\mathcal H_{\rm boost}=L^2(\mathbb R,da),\quad
(J_{\rm formal}f)(a,b)=f(a+b) .
\tag{10}
$$

后一个是检验中的**正则法向参考表示假设**，不是从认知原则或原引力约束推出的Hilbert空间。其匹配作用为(a,b)→(a+t,b−t)。换元u=a+b后，任何不变函数都沿a方向几乎处处常数，于是

$$
\|J_{\rm formal}f\|^2
=\int_{\mathbb R}da\int_{\mathbb R}du\,|f(u)|^2=\infty
\quad(f\ne0),\qquad
L^2(\mathbb R^2)^{\mathbb R}=\{0\} .
\tag{11}
$$

可用平移不变L²函数的分布导数为零证明最后一式；并非仅检验一个Gaussian。若另限制到这个空间的某个谱子空间，也不会突然产生非零正常不变态。但其他表示或重新定义约化内积不在该结论中。

## 6. 有限规范区间的近似不提供强极限

最直接的归一补救是引入一个计算区间：

$$
(J_Lf)(a,b)=\frac{\mathbf1_{[-L,L]}(a)}{\sqrt{2L}}f(a+b),\qquad
\|J_Lf\|=\|f\|,\quad
\|U_tJ_Lf-J_Lf\|^2=\min\{|t|/L,2\}\,\|f\|^2 .
\tag{12}
$$

这来自两个平移区间的对称差长度。固定t时缺陷趋零，但对归一f，两个越来越大区间的嵌入满足

$$
\langle J_Lf,J_{4L}f\rangle=\tfrac12,\qquad
\|J_Lf-J_{4L}f\|^2=1 .
\tag{13}
$$

故J_Lf不构成强Cauchy族，不能把“小匹配误差＋始终归一”当作存在非零正常物理极限的证明。L只是诊断错误接法的规范区间，不被设为宇宙资源阈值。这里的Gaussian测试态也不是已经解了引力全部约束的物理态。

可以在适当稠密测试域上考虑未归一群平均的分布配对：

$$
\langle\eta(\Psi),\eta(\Phi)\rangle_{\rm phys}
\ stackrel{?}{=}\ \int_{\mathbb R}dt\,\langle\Psi,U_t\Phi\rangle .
\tag{14}
$$

这不是有界Haar投影，必须另核正性、零范数商、内积、算符和过程的下降。它是保留的替代路线，不以(11)否定所有非紧约束量子化，也不在本轮直接宣布(14)已完成原引力拼接。

## 7. 核验与条件压缩

三组检查均通过：

- 原五场F、一般Lorentz基下的曲率缩并、二维Weyl面积和共同角点变分；独立作用差分末步误差≤6.31×10⁻¹²。
- 同一573初值24³、32³、48³的角点；原48³角度两种算法差小于10⁻¹²，角点荷缩并与Einstein面积一致，标量来源差分误差小于10⁻¹¹。
- 正则boost表示的独立积分，48／64点求积相符；L=2、8、32，t=.7时匹配缺陷平方分别为.35、.0875、.021875，而J_L和J_4L距离平方始终为1。无限区间和无强极限由(11)—(13)解析证明。

本轮将C02区域拼接、C08面积、C09物质参考、C11边界约束、C12原引力耦合及C22来源接到同一角点；原F确定的角点系数与标量来源不能再各自任意选。C13只增加几何荷与熵候选的归一接口，没有签收统计熵或Newton系数的生成。

完整区域量子化仍需要对非紧匹配、内部Gauss、物质过程和动态几何采用相容的物理内积。下一项[649](649/drafts/STATUS.md)检验替代约化能否保持原区域过程及来源，不继续增大L做精度扫描，也不转入认知装置设计。统一目标保持开放。
