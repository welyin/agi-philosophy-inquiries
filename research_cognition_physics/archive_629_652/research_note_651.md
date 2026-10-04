# 第651轮：原物质参考、因果区域与完整边界变分的共同连接

日期：2026-10-02。接[650](research_note_650.md)，回查[573](../archive_554_584/research_note_573.md)、[618](../archive_585_628/research_note_618.md)、[619](../archive_585_628/research_note_619.md)、[647](research_note_647.md)、[648](research_note_648.md)及[651起步](651/drafts/material_wall_joint_entry.md)。[代码](651/joint_material_boundary_transport.py)、[结果](651/joint_material_boundary_transport_results.json)、[核验](651/research_round_651_checks.json)、[条件账](651/unified_physics_condition_ledger_651.md)。三组复算、十六式；主代理审查，无新增独立代理审查。

## 1. 共同条件与实际增量

650的边界动量和区域通量需要类时切口。573的原物质参考虽已满秩，固定其中一个标签不自动给类时边界。本轮在同一受约束解上证明固定s面类空，并给保持原作用不变的局部类时选面。648已有原参考的法向Gram资料和角点；此处不把这些资料再报为发现，新增的是其与类时区域条件及完整变化的连接。

采用同一物质图册后，参考、关系度规和边界来源不能分别自由指定：两个逆度规组合由参考标签固定，原F也成为标签的既定函数；但法向、外曲率和法向物质导数仍随度规变化。将这些资料共同拉回，原体作用、边界项和角点项有同一个经典移动区域表示，原完整透射匹配得以共同输送。

本轮完成的是限定局部经典接口。未证明所有候选边界辛结构等价，未将经典图册直接提升为量子参考或instrument，也没有改动统一目标。

|层次|采用或保留|
|---|---|
|认知动机|同一物质兼任参考和来源时，区域、几何与变化规则应共同确定|
|继承输入|573同源Einstein—Yang–Mills—sigma解，3+1、正F、给定领先作用；费米背景为零|
|局部输入|同一可逆参考片、固定关系标签域、光滑且不改变因果类型的边界|
|新增原对象连接|固定s面失败、局部类时替代；原参考约束关系度规及F；完整边界几何共同输送|
|成熟工具|场依赖坐标、固定域拉回、标准边界变分、Reynolds与Noether身份|
|数值范围|原573初片的因果诊断；另一个自洽关系图册上的离壳几何差分|
|开放|量子物理态、区域量子过程、跨尺度映射、全局参考、全EFT边界及引力生成|

[Brunetti—Fredenhagen—Rejzner，§2.5—2.6](https://arxiv.org/html/1306.1058v5)提供字段依赖定位方法；其增加四个标量的例子不是本模型的原物质菜单。这里使用原h、s及其导数，不借该文签收本项目的量子化。[Harlow—Wu，§3.5、式(102)—(111)](https://arxiv.org/html/1906.08616v3)明确其法向变分将选面函数固定；本轮先把实际物质选面共同拉回，才使用这种固定面公式。

## 2. 原参考图册不能自动代替区域的因果条件

沿573，X=(h,s,A,B)，A=(dh)²、B=(ds)²，C=dh·ds，均用同一Einstein度规。取λ为固定选面参数而非动力学耦合，定义

$$
Y=(h,s-\lambda h,A,B),\qquad f=Y^1-b_0,\qquad
N^2=(df)^2=B-2\lambda C+\lambda^2A,\qquad \det(dY)=\det(dX).
\tag{1}
$$

因此Y仍是原局部图册。N²>0才使f=0成为类时侧壁；N²<0给类空面。图册满秩本身没有指定每张坐标面的因果类型。

在573原初片点p=(0,π/2,π/4)，空间dh=0、ds=−b dy。令H=1.05h★、S=s★、b=.06S、d=HS/12、β=1−S²/12、p_s=.025−.0056/b，则原速度为

$$
v_h=C_h\psi^{-6},\quad C_h=-Fdp_s,\qquad
v_s=C_s\psi^{-6},\quad C_s=F\beta p_s,\qquad
A=-C_h^2\psi^{-12},\quad C=-C_hC_s\psi^{-12},\quad
B=b^2\psi^{-4}-C_s^2\psi^{-12}.
\tag{2}
$$

原连续界F>1.8、β>.97、|p_s|>.14、b<.033、ψ<3/2给

$$
\frac{C_s^2}{b^2}>\left(\frac{1.8\cdot .97\cdot .14}{.033}\right)^2
>54.86>\left(\frac32\right)^8>\psi^8
\quad\Longrightarrow\quad B<0.
\tag{3}
$$

所以固定s面不是618／650所需的类时边界。这是同一原解的解析结论，不由网格外推，也不表示原物质坐标失效。

同一资料给出一项相容选择：在p计算λ★，随后将它作为常量固定；

$$
\lambda_\star=\frac CA=\frac{C_s}{C_h}=-\frac{12-S^2}{HS},\qquad
N^2(\lambda_\star)=b^2\psi^{-4}>0,\qquad
N^2(\lambda)>0\ \Longleftrightarrow
|\lambda-\lambda_\star|<\frac{b\psi^4}{|C_h|}\quad\text{在p}.
\tag{4}
$$

连续性提供p附近的类时面及小场邻域。λ★与数值ψ无关；它是一次区域选择，不是每个场配置中重新调优的参数。32³诊断给λ★≈−30.6815268210，N²≈.000997747310595，允许区间约(−34.7420138732,−26.6210397688)。没有全局或长期侧壁保证。

## 3. 原参考与关系度规不再是两套独立资料

在同一可逆片取e_Φ=Y_Φ⁻¹，全部场统一拉回，记上横线。原定义强制

$$
\bar h=Y^0,\qquad \bar s=Y^1+\lambda_\star Y^0,\qquad
\bar g_E^{00}=Y^2,\qquad
\bar g_E^{11}+2\lambda_\star\bar g_E^{01}+\lambda_\star^2\bar g_E^{00}=Y^3,\qquad
\bar F=2-\frac{(Y^0)^2+(Y^1+\lambda_\star Y^0)^2}{6}.
\tag{5}
$$

最后一式使用h²=|Higgs|²，保留全部五场；并未要求Higgs的内部取向恒定。A、B是同一物质的一阶导数收缩，故前两个逆度规组合不能另配任意函数。固定Y的允许变化必须满足

$$
\delta\bar h=\delta\bar s=\delta\bar F=0,\qquad
\delta\bar g_E^{00}=0,\qquad
\delta\bar g_E^{11}+2\lambda_\star\delta\bar g_E^{01}=0,\qquad
\delta\bar g_E=\bar F\,\delta\bar g_J.
\tag{6}
$$

这是图册自洽身份，不是场方程，也没有冻结物理径向变化：物质变化同时改变e_Φ和关系度规。更不能从两项代数限制直接计算物理自由度，完整规范约束、边界及导数依赖仍在。647的点读函数退化成标签函数与这里相容；非平凡区域内容仍可在几何及其权重中。

## 4. 同一移动区域的几何变化

设u=δe∘e⁻¹。对每个变化，原逆参考给u而非任意边界驱动；在侧壁上

$$
u^\mu=-(dY^{-1})^\mu{}_I\delta Y^I,\qquad
\delta\bar g=e^*(\delta g+\mathcal L_u g)=:\kappa,\qquad
\delta\bar\phi=e^*(\delta\phi+\mathcal L_u\phi)=:\sigma,\qquad
n\cdot u=-\frac{\delta s-\lambda_\star\delta h}{\sqrt{N^2}}.
\tag{7}
$$

同样输送全部规范连接；内部规范改变须配同一bundle识别。本式是成熟链式法则，不独立算研究成果。实际要求是将它共同用于法向、外曲率和原动量。

下文在固定Y域计算，省略横线，n²=+1，e_a为侧壁固定切向基，q_ab为诱导Lorentz度规，κ_nn=n^μn^νκ_μν。单位法向随度规变化为

$$
\delta q_{ab}=e_a^\mu e_b^\nu\kappa_{\mu\nu},\qquad
\delta n_\mu=\tfrac12\kappa_{nn}n_\mu,\qquad
\delta n^\mu=-g^{\mu\rho}\kappa_{\rho\nu}n^\nu+\tfrac12\kappa_{nn}n^\mu.
\tag{8}
$$

由K_ab=e_a^μe_b^ν∇_μn_ν直接得到

$$
\delta K_{ab}=\tfrac12\kappa_{nn}K_{ab}-e_a^\mu e_b^\nu n_\rho\delta\Gamma^\rho{}_{\mu\nu},\qquad
\delta\Gamma^\rho{}_{\mu\nu}=\tfrac12g^{\rho\sigma}
(\nabla_\mu\kappa_{\nu\sigma}+\nabla_\nu\kappa_{\mu\sigma}-\nabla_\sigma\kappa_{\mu\nu}).
\tag{9}
$$

式中∂κ_nn乘n_ν的项因切向投影消失。物质法向导数v^A=n^μD_μφ^A的变化还包括同一规范连接a的变化：

$$
\delta v^A=\delta n^\mu D_\mu\phi^A+
n^\mu\big[D_\mu\sigma^A+(\delta a_\mu)^A{}_B\phi^B\big].
\tag{10}
$$

即使σ=0，δn也能使δv≠0。因而固定参考标签绝不允许顺便冻结法向来源。原Jordan边界动量继续使用同一F、f_A=−φ_A/3：

$$
\Pi_J^{ab}=\frac{\sqrt{|q_J|}}2\big[F(K_Jq_J^{ab}-K_J^{ab})+q_J^{ab}f_Av_J^A\big],\qquad
\pi_{J,A}=\sqrt{|q_J|}(f_AK_J-v_{J,A}),\qquad
\delta\Pi_J,\delta\pi_J=\delta\big[\Pi_J(q_J,K_J,\phi,v_J),\pi_J(q_J,K_J,\phi,v_J)\big].
\tag{11}
$$

最后一项明确对所有参数取同一变化。不得把式(6)的δF=0误读成δK=δv=0。Einstein／Jordan两框架的边界一形式沿618相同；当前f·σ=δF=0只是限制了其中一项混合来源，没有消灭物质—几何耦合。

## 5. 体作用、边界与角点必须使用同一嵌入

取相对紧标签域D₀，全部配置的e_Φ定义在同一D₀邻域，边界光滑分片且因果类型固定。继承原领先体作用L、带方向的GHY边界项ℓ及648角点作用I_S。逐个配置的协变性给精确身份

$$
S_{\rm rel}[\Phi]=S[\Phi;e_\Phi(D_0)]
=\int_{D_0}L[e_\Phi^*\Phi]+\int_{\partial D_0}\ell[e_\Phi^*\Phi]
+\sum_{S_0}I_{S_0}[e_\Phi^*\Phi].
\tag{12}
$$

这是同一作用的局部表示，包含度规、标量、规范连接；原573费米背景为零。没有把场依赖变换当作一个字段无关线性算符，也不据此证明所有约化量子理论等价。

在原移动域，体积分的变化为Reynolds公式；在固定标签域，变化为δΦ+Lie_uΦ。二者的边界辛势相差在壳Noether全导数：

$$
\delta\!\int_{D_\Phi}L
=\int_{D_\Phi}E\cdot\delta\Phi+\int_{\partial D_\Phi}[\Theta(\delta\Phi)+i_uL],\qquad
\Theta(\mathcal L_u\Phi)=i_uL+dQ_u\quad\text{在壳}.
\tag{13}
$$

第二式在每一条变化上使用当时的u；再次变化必须保留δu。单个面上的dQ_u会给角点项，不能删掉后声称两种区域辛势相同。完整闭合边界中的取向取消、GHY变分的边界全导数以及显式角点作用均须共同保留。这里没有选择或证明唯一的扩展边界相空间。

在固定标签域，远离角点的在壳侧壁一形式，连同未丢弃的规范项，具有

$$
\vartheta_{B_0}=\Pi^{ab}\delta q_{ab}+\pi_A\sigma^A+
\sum_I\mathcal E_I^{a}\delta\bar a_{I,a},\qquad
\mathcal E_I^a=-\sqrt{|q|}\,\kappa_I n_\mu\mathcal F_I^{\mu a},\qquad
\delta I_{S_0}=\int_{S_0}(c\,\delta\eta+\eta\,\delta c).
\tag{14}
$$

κ_I为原规范动能−κ_I F_I²/4的系数，内部Lie指标按原不变内积收缩；没有删除573非零规范场。式(14)还应与GHY变分产生的角点全导数按相同约定合并，或对支撑远离角点的变化使用侧壁部分。c、η由同一拉回几何计算，不能换成固定嵌入的旧变化。

## 6. 完整透射连接可以共同输送，但不保证单区自治

在同一光滑解的两侧使用**同一个e_Φ**、相反向外法向与共同内部规范识别。原诱导字段识别及动量匹配共同拉回后，仍有

$$
\Pi_++\Pi_-=0,\qquad \pi_++\pi_-=0,\qquad
\mathcal E_++\mathcal E_-=0
\quad\Longrightarrow\quad
\vartheta_{B_0,+}+\vartheta_{B_0,-}=0,\qquad
\delta\vartheta_{B_0,+}+\delta\vartheta_{B_0,-}=0.
\tag{15}
$$

须在完整共同匹配子空间取外微分，包括嵌入变化和允许数据变化。两侧物理体域互补、共享几何边界，形状位移的体项也按相反方向抵消；同一角点取向约定下内部角点数据一起取消或合成整体所需项。此处复用光滑解的共同限制，不证明任意两组独立区域初边值都可拼成解。

如果两边各冻结或各选一个独立e，式(15)没有得到保证。正确输送只说明整体组合相容，不能由δθ总和为零推出每边δθ为零，也不撤销650的单侧通量障碍。650常真空的原X退化，其数值反例没有被搬到当前573图册上。

## 7. 可复算检查及其证明角色

第一组复用原573解的16³、24³、32³诊断；另用精确有理数核式(3)的54.8676892562>25.62890625。原类时选面与λ★独立于ψ的身份单独核对。

第二、三组使用一个**离壳诊断族**，并非伪造原573未知的后两个参考时间导数。令Y基点取原点的(H,S−λH,A,B)，设q(Y)为明确仿射函数，构造

$$
g^{00}=Y^2,\quad g^{01}=\varepsilon q(Y),\quad
g^{11}=Y^3-2\lambda\varepsilon q(Y)-\lambda^2Y^2,\quad
\phi=(0,Y^0,0,0,Y^1+\lambda Y^0),\qquad
Y_\Phi=(h,s-\lambda h,(dh)^2,(ds)^2)=Y.
\tag{16}
$$

另两对角分量为代码中给出的光滑正函数，其余混合分量为零。小ε内签名及N²>0保持。这个族的参考确实由同一场重建，故可以独立核式(5)—(11)，而非先指定四标签再假设它们自洽；它未满足Einstein方程，不能用作新在壳例子。

代码从完整度规及其导数重建Christoffel、单位法向和K，再以对称差分与解析变分独立比较。三个步长2×10⁻⁶、10⁻⁶、5×10⁻⁷的误差按约四倍下降，最终归一误差小于6.53×10⁻⁶。在Einstein诊断中δφ=0而δv_h≈.949753651，冻结法向会完全漏掉它。两框架原边界一形式分别为.024362922313288796与.02436292231328881，差约1.39×10⁻¹⁷。它们是诊断族的参数导数，不是原宇宙的观测数值。

实现审查中修正了旧参数函数不返回d的读取错误，并因小g⁰⁰引起的差分截断误差缩小步长；没有放宽核验阈值。正式结果在这些修正后生成。无图像检验。

## 8. 减少的独立输入与下一项

本轮合并了C02区域组合、C05因果类型、C08测度／几何、C09参考和C22来源：原物质图册确定两项关系度规身份及F；同一嵌入确定允许边界变化，并保持既有整体透射连接。剩余选面参数、全局分支、边界辛势选择与量子实现仍须单列，不能宣称已唯一确定。

不继续优化类时选面的宽度或参考条件数。[652入口](652/drafts/STATUS.md)回到跨分支共同量子对象：原h、s和含速度／几何的A、B怎样与完整Gauss过程、区域态及来源共用一个可定义的操作表示。先回查旧量子参考和约束结果；一般不对易公式、给算符加帽子或单独设计仪器均不计接通。
