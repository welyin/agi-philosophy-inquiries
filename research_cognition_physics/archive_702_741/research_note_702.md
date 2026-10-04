# 第702轮：原完整过程的正转移、共同热态与记录的联合极限

日期：2026-10-02。接[701](../archive_653_701/research_note_701.md)、[702入口](702/drafts/source_scale_entry.md)。[代码](702/joint_gibbs_preserving_transfer.py)、[结果](702/joint_gibbs_preserving_transfer_results.json)、[核验](702/research_round_702_checks.json)、[全条件账](702/unified_physics_condition_ledger_702.md)。两组检查、十六式；主代理审查，无独立代理审查。

## 1. 原对象中的真实缺口

上一目标轮完成701并执行702入口，属于有效进展。本轮先读取导航、最新报告／结果及入口，未发现运行中的Python研究。目标和先整合条件的顺序保持不变。

[655](../archive_653_701/research_note_655.md)已经把正分割接到原全部相互作用的H_F，证明强预解、真实时间和**固定原态**上的记录收敛；其§4明确没有得到热迹的迹范数收敛。[625](../archive_585_628/research_note_625.md)的有限CP压缩也保留指定初态，不允许随意重选热态。[603](../archive_585_628/research_note_603.md)证明原Gibbs态存在，但这不自动证明一列近似Hamiltonian自己的Gibbs态收敛。

**本轮补齐的连接：** 对同一个原有限图H_F重新分配已有束缚势，构造一族正转移，使其自身的Gauss Gibbs态在迹范数中趋于原态；同一族的真实时间、原记录联合后态、全部热能量矩及热熵共同收敛。没有删除Dirac／Majorana、空间跳跃、曲目标或规范电动能。

这是新的明确分割，有限步长时不等于655的质量因子分割，更不等于673的overlap／S9候选。699的反例继续成立；本轮提供一个可用于后续比较的正基准，尚未证明它与手征候选或连续场论相同。

|层次|地位|
|---|---|
|认知动机|参考态、演化、记录和能源必须来自同一个近似过程|
|继承输入|原固定有限图、H⁵场目标、商群、完整CAR、原耦合及严格正几何；β>0平衡态|
|新增表示选择|0<θ<1，把θW移到热生成元内；原总H不变|
|解析新增|有效生成元的统一束缚下界、热尾紧性、迹范数收敛、共同Gibbs／记录／能源与熵|
|数值验证|原64 CAR条件纤维下界；原径向形式与中性Majorana量子耦合的256维校准|
|仍开放|空间细化一致性、一般几何来源导数、热化来源、手征映射及量子引力|

### 1.1 成熟结果与历史复用

迹理想中的乘积收敛是成熟研究；[Cachia–Zagrebnov原论文](https://doi.org/10.1112/S0024610701002526)给一类以Gibbs半群为因子的迹范数定理。不能直接把655的裸动能K当成Gibbs生成元：非紧目标上的裸热迹一般无穷。本轮选取原束缚势的一部分，明确补齐这个输入，并另给下文的直接紧性证明。

强预解／log重建沿655已有证明；复时间导数使用655已核的[Arendt–Nikolski向量值Vitali定理](https://www.uni-ulm.de/fileadmin/website_uni_ulm/mawi.inst.020/arendt/downloads/pubbib/short/2000-AreNik-VctVldHlmFncRvs.pdf)。成熟定理不计原创，新增的是它们在原完整模型及**变化中的自身热态**上同时成立。

空间旧结论按[701入口§5](../archive_653_701/701/drafts/joint_scale_entry.md)直接复用：382／383上界，384删除额外Lipschitz；386真实邻域桥或425半幅／成本桥；522—523热参考和实际方向仪器共同实现。524局域探针已完成，亦不重做。当前h、s、CAR／Gauss到这些接口的映射仍未补齐，辅助E不是记录s。

## 2. 同一个H，给一份已有束缚势两个用途

在623／655的固定有限图及同一Gauss Hilbert空间上，写H=K+V+J+C+C†。K≥0包含原标量／非对角电动能，V≥W≥0含原场势、边势、磁势；J为规范跳跃，C为614字典下的完整质量配对。定义

$$
A_\theta=K+\theta W,\qquad
D_\theta=V-\theta W+J+C+C^\dagger,\qquad
H=A_\theta\dotplus D_\theta,\quad 0<\theta<1.
\tag{1}
$$

点表示同一闭形式和；在原共同算符核上就是原微分算符。θ是分割处方，不是新物理耦合。由655已有‖J‖≤j₀、‖C‖≤c₀+c₁W^(1/4)，

$$
D_\theta\ge-c_\theta I,\qquad
c_\theta=j_0+2c_0+
\frac{3(2c_1)^{4/3}}{4^{4/3}(1-\theta)^{1/3}}<\infty.
\tag{2}
$$

直接用sup_{x≥0}[b x^(1/4)−r x]=3b^(4/3)/(4^(4/3)r^(1/3))。原非线性和所有质量仍保留。Aθ按603的相同势尾和热核比较，对每个正β有有限热迹；紧群电动能未删除，有限Fock身份因子也完整保留。Dθ是自伴矩阵乘法，Gauss协变。原共同核保证形式和仍为623的H。

## 3. 正转移及关键的统一算符下界

对任意a>0定义

$$
S_a=e^{-aA_\theta/2}e^{-aD_\theta}e^{-aA_\theta/2},
\qquad 0<S_a\le e^{ac_\theta}e^{-aA_\theta},\qquad
H_a=-a^{-1}\log S_a.
\tag{3}
$$

严格大于0意为正且核为零，不指有统一严格正谱隙。中间指数有界且单射；两侧热算符同样单射，所以log由谱定理良定义。S_a紧，H_a自伴且具有紧预解。该处方需要Aθ的禁闭热算符，没有声称比655裸动能分割更易计算。

关键是log的算符单调性。对有界正可逆算符由log的预解积分和逆算符序直接得到，非可逆下界情形以ε正则化再取闭形式极限：

$$
\log S_a\le ac_\theta I-aA_\theta
\quad\Longrightarrow\quad
H_a+c_\theta\ge A_\theta\ge0.
\tag{4}
$$

这里没有使用错误的“算符指数保持一般算符序”。下文的热迹比较使用min–max本征值序。cθ只是证明移位，原H、H_a与未移位配分函数均不改；尤其不将依赖几何的cθ扣成某个免费真空项。

沿原共同核微分S_a得−H；稳定界由(3)给出。完全复用655对(I−S_a)/a及−log S_a/a的预解比较：

$$
H_a\xrightarrow[a\downarrow0]{\rm strong\ resolvent}H,
\qquad e^{-\beta H_a}\xrightarrow{\rm strong}e^{-\beta H}
\quad(\beta>0).
\tag{5}
$$

本轮不把该旧强收敛再次当热迹结论；需要下面新的统一紧性。

## 4. 高能尾部不能随a逃走

全部迹先限制在同一个Gauss物理空间。记A=Aθ、L_a=H_a+cθ、R_a(β)=e^(−βL_a)，极限L=H+cθ；由(4)及H=A+D≥A−cθ，二者均≥A。设Z_A(t)=Tr_phys e^(−tA)<∞。谱比较及x e^(−βx)≤2e^(−βx/2)/(eβ)给

$$
\operatorname{Tr}R_a(\beta)\le Z_A(\beta),\qquad
\operatorname{Tr}[A R_a(\beta)]\le
\operatorname{Tr}[L_aR_a(\beta)]
\le M_\beta:=\frac{2}{e\beta}Z_A(\beta/2).
\tag{6}
$$

A和L_a不必对易。中间式由A≤L_a在L_a本征向量上的形式期望逐项相加得到；正迹用递增谱截断定义。式(6)对a统一，也对极限R成立。

取P_R=1_{A≤R}，它有限秩且保Gauss。正算符的块分解及Hilbert–Schmidt Cauchy–Schwarz给

$$
\operatorname{Tr}[(I-P_R)R_a]\le\frac{M_\beta}{R},\qquad
\|R_a-P_RR_aP_R\|_1
\le2\sqrt{\operatorname{Tr}R_a\,
\operatorname{Tr}[(I-P_R)R_a]}
\le2\sqrt{\frac{Z_A(\beta)M_\beta}{R}}.
\tag{7}
$$

因此

$$
\|R_a-R\|_1\le
\|P_R(R_a-R)P_R\|_1
+4\sqrt{\frac{Z_A(\beta)M_\beta}{R}}
\longrightarrow0.
\tag{8}
$$

先固定R，由(5)和有限秩使第一项趋零，再令R趋无穷。**这正是不能仅由强收敛替代的一步。** R是证明用谱截断，不删除原物理高能态；所有状态和演化仍由完整S_a／H_a定义。没有给出跨空间细化统一的Z_A，也没有实际数值计算原完整谱投影。

## 5. 近似过程自己的热态确实恢复原热态

由S_a正且Gauss空间非零，Z_a(β)=Tr_phys e^(−βH_a)严格正且有限。603原Z同样严格正。恢复证明移位后，

$$
Z_a\to Z,\qquad
\rho_{a,\beta}:=\frac{e^{-\beta H_a}}{Z_a}
\xrightarrow{\|\cdot\|_1}\rho_\beta,
\qquad
\|\rho_{a,\beta}-\rho_\beta\|_1
\le\frac{2}{Z}\|e^{-\beta H_a}-e^{-\beta H}\|_1.
\tag{9}
$$

不再借用原ρβ作为所有近似的固定初态，也不先假定Z_a有共同正下界。其下界最终由Z>0及收敛产生。这里证明的是平衡态近似，没有证明封闭系统自动热化或宇宙为何选择这个β。

## 6. 同一H_a的真实时间和实际记录

强预解及共同下界沿655给

$$
U_a(t)=e^{-itH_a}\xrightarrow{\rm strong}U(t)=e^{-itH},
\quad\text{在有限实时间区间对每个固定向量一致。}
\tag{10}
$$

ℏ取原单位，恢复时两边同除以ℏ。保留624的原有界、保Gauss读取L_r(s)，同一有限记录词为C_{a,r}=L_{r_m}U_a(t_{m-1})⋯L_{r_1}。定义含经典记录与未归一条件后态的联合输出Ω_a。有限乘积的强收敛、(9)及迹理想连续性给

$$
\Omega_a:=\sum_{\boldsymbol r}|\boldsymbol r\rangle\langle\boldsymbol r|
\otimes C_{a,\boldsymbol r}\rho_{a,\beta}C_{a,\boldsymbol r}^\dagger,
\qquad \|\Omega_a-\Omega\|_1\longrightarrow0.
\tag{11}
$$

具体将差拆成初态误差及固定ρβ下的词误差；前者由整体记录CPTP映射不增迹距离控制，后者逐支受2‖(C_a−C)ρβ^(1/2)‖₂控制并趋零。故记录概率和非零概率分支的归一后态都收敛；概率趋零的后选择分支不领取统一相对误差。

附加被动参考时，同样结论要求给定的联合初态本身迹范数收敛；仅系统边缘收敛不足以指定未知系统—参考关联。不把Euclidean热权当成真实测量概率，也不把读取权限说成已经由H自主设计出来。

## 7. 能源、常lapse来源与热熵也使用这一族

由L_a≥A的谱比较，在复时间右半平面有

$$
\|e^{-zL_a}\|_1
=\operatorname{Tr}e^{-(\operatorname{Re}z)L_a}
\le Z_A(\operatorname{Re}z),\qquad \operatorname{Re}z>0.
\tag{12}
$$

各族在迹类Banach空间中全纯。式(8)给正实轴上的迹范数收敛，局部一致界(12)允许向量值Vitali及Cauchy导数公式。还原同一固定cθ：

$$
H_a^j e^{-\beta H_a}\xrightarrow{\|\cdot\|_1}
H^j e^{-\beta H},\qquad j=0,1,2,\ldots,\quad\beta>0.
\tag{13}
$$

这里是**固定a的同一有效H_a**的能源插入，不把随β改变分割步长的导数冒充H_a。常lapse η作用于完整H_a；它的一、二阶热来源直接由(13)给出，包含原动能、质量、跳跃和势的共同极限。

于是每个固定j的热能量矩及von Neumann热熵均有

$$
\operatorname{Tr}(\rho_{a,\beta}H_a^j)\to
\operatorname{Tr}(\rho_\beta H^j),\qquad
S(\rho_{a,\beta})=\beta\operatorname{Tr}(\rho_{a,\beta}H_a)+\log Z_a
\longrightarrow S(\rho_\beta).
\tag{14}
$$

不是仅由一般迹距离推出无限维熵连续，而是用Gibbs身份和已证明的能源矩。尚未证明任意空间度规／局部lapse变化的全部来源导数，更没有由热熵收敛推出面积律或Newton常数。

## 8. 正性、表示参数与原手征候选的区别

每个a都有同一个正Gauss量子过程。对于有界、偶、保Gauss的正半区Euclidean观测词，按伴随和反向时间次序反射，在热圆上把两个半区拼成

$$
Q_a(F,F)=Z_a^{-1}\operatorname{Tr}_{\rm phys}(X_{a,F}^\dagger X_{a,F})\ge0.
\tag{15}
$$

正热时间间隔保证所需迹理想条件；这是原Hilbert过程的正性，不为673的另一来源映射自动签收费米反射。涉及奇费米观测时还须保持原分级约定，不能跳过691的物理奇Gauss接口。

不同θ只重新分配原势。对任意两份固定θ、θ′∈(0,1)，同一β下有

$$
\lim_{a,b\downarrow0}
\|\rho^{(\theta)}_{a,\beta}-\rho^{(\theta')}_{b,\beta}\|_1=0,
\tag{16}
$$

真实记录及上述能源／熵同样回到原对象。不要求θ依赖尺度趋近0或1；那将使本证明的热迹或下界常数退化，须另查。

这一正族不能在699已否定的固定参数处精确保留其全部反射来源。即使两边都叫“路径积分”或“手征表示”，也仍需具体观测、尺度和相互作用映射。原overlap物种控制不能免费转移到此族；本轮也不证明H_a在有限a有空间一致局域性。

## 9. 实际核验、适用范围与下一项

第一组沿655原字典保留两节点全部64 CAR模式、原Dirac／Majorana和规范跳跃，用实际BdG最低Fock能源核(2)的全配置Young下界。接近F=0的样点仅校准下界实现；一般非紧结论由原势控制承担，不靠扫描。

第二组复用625的64点(h,s)径向形式及原读取，把中性Majorana Fock四维因子作为量子耦合纳入，得到256维非对易校准。此诊断明确取Dirac Yukawa为零、带荷Fock真空；在单节点、无边模型中中性singlet与径向部门是合法不变部门。有限径向盒和差分仍只是计算校准，**不是原全图热迹或Gauss积分**。解析定理不删任何原耦合。

θ=1/4、1/2、3/4，a从0.2减至0.025，log算符下界、热迹／态、记录联合后态、一二阶能源及热熵均核验。最后一档Gibbs迹范数差分别约8.76×10⁻⁶、7.34×10⁻⁶、6.10×10⁻⁶。不同θ指向同一个校准H；没有用有限矩阵结果证明无限维收敛或空间连续极限。

**本轮合并：** C01态、C03记录、C04演化、C09热参考、C13热熵和C20时间近似，现在在原完整有限图中共用一族自身归一的正过程。这削减了“还需另外挑选近似热态及归一”的独立输入；未削减温度、热化、群、维数、耦合或正背景的输入。

所有热迹常数、Fock维数、cθ和A的谱均可依赖图及几何下界；不抵消578—579的空间细化问题。给定作用的经典几何、指定连续手征物质、原正H_F和辅助候选四分支尚未统一。接[703](703/drafts/STATUS.md)：核这个共同热过程对实际物理来源变化的响应，尤其区分有界记录／规范来源和会改变动能的一般几何来源，不把态收敛自动升级为全部反作用收敛。
