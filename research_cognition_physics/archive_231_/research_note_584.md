# 第584轮：共同量子测度、原物质记录与几何响应

日期：2026-10-01。接[583冻结条件](research_note_583.md)、[574原量子模型](research_note_574.md)、[577原物质记录](research_note_577.md)及[578能源界](research_note_578.md)。链接按archive_231_解析；审查来源见[核验记录](round584_drafts/final_review.txt)，完成状态见[冻结核验](research_round_584_checks.json)。独立实质审查已收到；末稿文字修正与最终哈希由主代理核验，未取得独立最终哈希签名。

## 1. 本轮接通了什么

583证明固定条件本身也须随框架变换，并未确定完整引力量子测度。本轮保留574已经声明的有限图Laplace–Beltrami处方，把其曲测度表示转换到平直坐标测度，检查演化、Gauss、记录及背景能源响应能否一起保存。

**可以，但必须保留由原算符确定的非恒定量子补项。** 本轮给出完整有限图的酉闭域映射；它保留原准备、记录、能源谱及所列背景权偏导。删掉这个补项仍可定义合法、保持Gauss的量子模型，却会使同一个正规Gauss源、同一个577读口的二阶时间概率不同，也会改变背景几何响应。

这不是从认知原则唯一选出了LB排序；是在既有分支内排除一种错误的“等价简化”。它也不把平直坐标测度解释成平直目标动力学或完整Jordan引力路径积分测度。

## 2. 输入、成熟来源与去重

|层次|本轮地位|
|---|---|
|认知动机|同一来源的状态、演化、记录与资源响应应使用同一表示字典|
|继承输入|固定有限图、原规范群、正空间几何权、五维H⁵目标及LB排序、原非负现场／边／磁势|
|解析增量|原五维Cartesian半密度的明确q；完整图Gauss与闭域映射；删q的同源记录和背景响应差|
|数值范围|五维局部微分算符、Gauss紧支撑波函数积分、原节点二次型与现场势的几何偏导|
|保留边界|无图连续极限、无完整引力测度、无协变量子应力；有限图各节点Hilbert空间仍无限维|

574已使用半密度说明Gaussian移植的边界失败；556／564已有径向测度产生逆平方项。本轮不重复它们。本轮q来自完整五维Cartesian测度，连接的是原全图模型、同一读口及背景权，不是另选二维径向量子化。

[Shubin，math/0007019v2，定理1.1及第1节](https://arxiv.org/pdf/math/0007019)处理完备主符号度量、任意正光滑测度及半有界Schrödinger算符。574的有限乘积完备性、光滑非负势已满足其条件，本轮直接继承闭实现。本文的配置磁一形式为零，不表示省去原动态规范链路。酉共轭与有界扰动是成熟算符工具，不作为新基础理论计数。

## 3. 完整有限图上的测度字典

记F=M−|φ|²/6>0、M=2，G=K⁻¹=F(I−φφᵀ/(6M))，mu=√M F⁻³。节点集为𝒱、链路集为𝓔，紧规范群仍记G_gauge。保持每条链路的Haar测度，定义

$$
\mathcal H_K=L^2\!\left(\prod_{v\in\mathcal V}\mu_vd^5\phi_v\prod_{e\in\mathcal E}d\mathrm{Haar}_e\right),\quad
\mathcal H_0=L^2\!\left(\prod_vd^5\phi_v\prod_ed\mathrm{Haar}_e\right),\qquad
(T\Psi)(\phi,U)=\left(\prod_v\mu_v^{1/2}\right)\Psi(\phi,U).
\tag{1}
$$

T是两个不同Hilbert空间之间的满酉映射。F→0时乘子发散，因此不能说它在同一个平直测度空间上是有界点乘；这不妨碍式(1)的酉性。T及其逆都把球内部的光滑紧支撑核映到自身。

令ell=log mu。直接对任意球内紧支撑光滑函数求导有

$$
T(-\Delta_K)T^{-1}=-\partial_a(G^{ab}\partial_b)+q,\qquad
q=\frac12\partial_a(G^{ab}\partial_b\ell)+\frac14(\partial_a\ell)G^{ab}(\partial_b\ell),\qquad
\partial\ell=\frac\phi F,\quad G\partial\ell=\frac F M\phi.
\tag{2}
$$

五维散度为5−7|φ|²/(6M)，平方项为|φ|²/M，所以

$$
\boxed{q(\phi)=\frac52-\frac{|\phi|^2}{3M}},\qquad
\frac12<q\le\frac52,
\qquad
H_c:=TH_KT^{-1}
=\sum_v\frac{\hbar^2}{2w_v}\big[-\partial G_v\partial+q_v\big]+H_{\mathrm{link}}+V.
\tag{3}
$$

w_v=epsilon³ psi_v⁶为原固定正节点体积。T与原所有乘法势交换，原测地边势及磁势未变；T不依赖链路，所以原正链路动能亦保持。目标逆度量G仍然是变系数，未换成平直I。

由原自伴算符定义

$$
\mathcal D(H_c)=T\mathcal D(H_K),\qquad
\mathcal D(H_c^{1/2})=T\mathcal D(H_K^{1/2}),\qquad
TP_K=P_0T.
\tag{4}
$$

前两式分别是闭算符域与闭形式域；P为两表示中的Gauss投影。最后一式由每个mu_v只依赖规范不变量|φ_v|²得到。故严格物理空间也酉对应，无需假定它按物理区域张量分解。

F=0是原完整目标的无穷距离端。换测度不能在有限坐标球面任意添加新的墙边界；这里的域是式(4)的继承域。

## 4. 原准备、读口及能源共同保留

令rho_c=T rho T⁻¹、E_c=TET⁻¹，并同样转换准备仪器。谱定理给任意时间的精确字典

$$
e^{-itH_c/\hbar}=T e^{-itH_K/\hbar}T^{-1},\qquad
\operatorname{tr}\!\left(\rho_c e^{itH_c/\hbar}E_c e^{-itH_c/\hbar}\right)
=\operatorname{tr}\!\left(\rho e^{itH_K/\hbar}E e^{-itH_K/\hbar}\right).
\tag{5}
$$

特别地，577的相位准备和正式二元读口都是乘法：

$$
U_\theta=\exp\!\left(\frac{i\theta|X_A|^2}{2\hbar}\right),\qquad
b=\sin s_B,\qquad E_+=\frac12+\frac b4,\qquad
TU_\theta T^{-1}=U_\theta,\quad TE_+T^{-1}=E_+.
\tag{6}
$$

整个读者代数仍须共同共轭；含微分或动量的算符不能一概保持旧表达式。局域乘积T及其规范不变性保证577的读者支撑与Gauss定义对应。

同一映射保持全部谱和能量矩，因此578的全态界也直接继承：

$$
\langle H_c\rangle_{\rho_c}=\langle H_K\rangle_\rho
\ge\frac{\hbar^2}{3}\sum_v\frac1{w_v}
\ge\frac{\hbar^2|\mathcal V|^2}{3\sum_vw_v}.
\tag{7}
$$

换测度没有消除固定ℏ、有限体积细化的能源障碍。574未截断Gaussian的有限能失败仍适用。这里图有限不等于全部量子引力已正则化。

## 5. 删除q为何是一个可检验的模型改变

在同一个H_0中定义

$$
Q=\sum_v\frac{\hbar^2q_v}{2w_v},\qquad
H_b=H_c-Q,\qquad
\|Q\|\le\frac{5\hbar^2}{4}\sum_vw_v^{-1},\qquad
\mathcal D(H_b)=\mathcal D(H_c).
\tag{8}
$$

固定有限图、正w下Q为有界实规范不变乘法，故这是同域的自伴有界扰动。H_b核上的散度动能及原势均非负，闭包仍半有界，并保持Gauss。界依赖图与体积，不能据此宣称网格一致。

这说明仅要求“合法量子模型＋保持Gauss”不能选定原排序。比较必须固定同一flat来源和读口字典。令D_j O=(i/ℏ)[H_j,O]，在共同核上

$$
D_cb=D_bb,\qquad
(D_c^2-D_b^2)b=-\frac1{\hbar^2}[Q,[H_b,b]]
=-\frac1{w_B}G_B(db,d_BQ)
=\frac{\hbar^2F_B^2s_B\cos s_B}{3M^2w_B^2}=:A_B.
\tag{9}
$$

D²表示两次生成元作用。原势、其它节点及链路的项在这个差中消去；没有把完整H替换成自由模型。

## 6. 严格Gauss来源与有限时间概率差

取每节点波函数只依赖h=|X|及s，链路部分取Haar常函数。B端选光滑紧支撑幅度

$$
u_B(h,s)=C\,f\!\left(\frac{h-.52}{.10}\right)f\!\left(\frac{s-.40}{.15}\right),\qquad
f(z)=\begin{cases}e^{-1/(1-z^2)},&|z|<1,\\0,&|z|\ge1.\end{cases}
\tag{10}
$$

按完整Cartesian体积归一，径向积分是h³ dh ds乘共同S³面积。支撑h∈[.42,.62]、s∈[.25,.55]，F≥Fmin=M−(.62²+.55²)/6，Fmin≈1.88551666667。函数在h=0邻域恒零，所以径向表达式延拓为完整五维光滑函数。其它节点任选规范不变紧支撑因子；全图源已严格Gauss，不需要丢弃角动能或假设物理空间先分解。

该源在两个H的任意所需幂域内，b也保持共同紧支撑核，故二阶时间Taylor展开成立。支持内A_B严格正。hbar=.7、w_B=.8时有明确界

$$
\langle A_B\rangle\ge
\frac{\hbar^2F_{\min}^2(.25)\cos(.55)}{3M^2w_B^2}
\approx0.0483439909373>0,\qquad
p_c(t)-p_b(t)=\frac{t^2}{8}\langle A_B\rangle+o(t^2)>0
\tag{11}
$$

对所有充分小正t成立。最终采用式(6)的577原归一，故系数为1/8；入口候选另举的(1+b)/2也是合法效果，其1/4不用于本轮正式比较。没有给统一可观测时间窗或无资源制备，只证明这两个明确处方在同源、同读口下可区分，也未排除另造完全不同的抽象酉字典。

## 7. 同一表示必须同时保留几何响应

让原正空间几何权按光滑参数lambda变化，固定目标K、M及ℏ，使用同一Hilbert平凡化。T不依赖这些空间权。因此在共同紧支撑核矩阵元意义下

$$
\partial_\lambda H_c=T(\partial_\lambda H_K)T^{-1},\qquad
\langle\partial_\lambda H_c\rangle_{T\Psi}
=\langle\partial_\lambda H_K\rangle_\Psi.
\tag{12}
$$

这是固定对应态的背景偏导，不是态随背景改变时的总能源导数，不需宣称参数族全部闭算符域相同。若省Q，则遗漏

$$
\partial_\lambda(H_c-H_b)
=-\sum_v\frac{\hbar^2q_v}{2w_v^2}\partial_\lambda w_v,
\qquad w_v=\epsilon^3\psi_v^6\ \Longrightarrow\
\partial_{\psi_v}(H_c-H_b)=-\frac{3\hbar^2q_v}{\epsilon^3\psi_v^7}.
\tag{13}
$$

原链接权和场势随几何变化时也必须共同变换；它们在两表示中完全保留，不影响式(13)给出的差。该有限图Hamiltonian响应不是已经构造的协变重整化应力或量子Einstein约束。

## 8. 四组实际复算

[代码](joint_quantum_measure_records.py)及[结果](joint_quantum_measure_records_results.json)复用原Python与NumPy，默认重新计算并比对完整JSON。入口两组首次纳入正式计数。新增紧支撑源直接核两种二次型：

$$
\int d\mu\,\big|\nabla_K(u/\sqrt\mu)\big|^2
=\int d^5\phi\left[(\partial u)^*G(\partial u)+q|u|^2\right].
\tag{14}
$$

这里没有再用u→r^(3/2)u的径向半密度，两侧都保留原h³体积；因此仍是五维Cartesian公式。

|检查|保存结果与边界|
|---|---|
|完整五坐标算符|独立差分原曲测度LB和flat散度；步长.004→.001，误差9.38×10⁻⁶降至5.85×10⁻⁷；省q的差趋2.5|
|局部双对易子|直接四项计算[Q,[H_b,b]]，解析A=.0896388397486，最细差8.19×10⁻⁸；不是时间演化模拟|
|紧支撑Gauss源|48、80、128点每轴；曲测度范数1，式(14)最细差6.17×10⁻¹⁴；〈q〉=2.42650266800、〈A〉=.0868278572612，原E初值.597228774058，概率t²系数.0108534821576，大于解析下界.00604299886717；原规范变换误差2.78×10⁻¹⁶|
|几何权及原现场势|epsilon=.8、psi=1.13；同源节点能源曲表示197.749351631、正确flat表示197.749351631，省q为197.191644451；缺失能源.557707179952。缺失背景偏导−2.96127706169，中心差分误差按步长²收敛；正确两表示的偏导差≤7.05×10⁻¹²|

最后一项只数值核原节点二次型与现场势；完整图所有项的对应由前述解析酉字典承担，不把节点数值当完整引力来源。紧支撑波函数积分给具体有限资源量，局部算符诊断和解析有限时间结论各自标明，没有虚构全图传播实验。

## 9. 联合条件账与后继

[584条件账](unified_physics_condition_ledger_584.md)把同一量子处方的测度、闭域、Gauss、记录、能源及所列背景响应接在一起。它减少的是分别指定这些表示的自由度；没有取消原排序、空间图、几何、维数或领先作用的独立输入。

583的受限标量作用现在有一套明确继承的有限图量子表示字典；这不等于553全部物种的固定Jordan背景运行已和574相同，也不是由完整引力积分推得条件测度。578—579的连续能源障碍完整保留。

下一步转向**同一有限尺度的量子来源与动态几何约束**：先回查既有均值反馈、半经典态与Einstein初值结果，明确量子能源、流和压力是否能同时进入同一约束。对仅有平均H的拼接给精确限制；优先使用现有共同物质，不继续枚举任意排序项或优化读口。
