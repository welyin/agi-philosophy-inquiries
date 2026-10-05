# 第802轮：原完整费米顶点的物理切片表达与连通响应中的抵消

日期：2026-10-05。接[801](research_note_801.md)、[冻结入口](802/drafts/STATUS.md)及[工作稿](802/drafts/research_note_802_working.md)。

配套：[完整系数校准](802/original_bff_vertex.py) · [结果](802/original_bff_vertex_results.json) · [联合复算](802/physical_bff_response.py) · [结果](802/physical_bff_response_results.json) · [核验](802/research_round_802_checks.json) · [复算入口](802/verify_round802.py)。

## 1. 本轮完成了什么

**将原连续作用的一玻色、两费米部门，具体写入785的物理切片和801的共同输出关联，并算出两项确切抵消。** 这里的一玻色、两费米顶点下称BFF顶点。

结果不是一般Kubo公式的再次陈述：原64维Nambu主部、全部Dirac/Majorana质量、度规平方根标架、Levi-Civita和内部连接均已进入同一显式系数。原参考固定使显式singlet质量导数消失，但几何顶点保留；单独非零的体积变化也不构成信号，其乘着的费米拉氏密度在本轮连通交叉收缩中严格为零。

**尚未得到原连续态与指定记录模式的非零总响应。** 本轮签收的是完整系数、抵消及对象映射，不把顶点的非零主部当成输出关联非零，也不把有限jet当成原时空解。803将代入同一原态的有序费米收缩，检查剩余项的总和。

|层次|地位|
|---|---|
|认知动机|同一材料同时承担参考、物质、几何与记录，不能分别计算后直接拼接|
|继承输入|753原四维背景、作用、群谱/参数、正F；730/767原自由核；正则片U和原实紧支完成|
|直接复用|598/604完整矩阵；664无独立挠率；773/785切片；780旋量字典；791/793—796共同量子分支；801输出菜单|
|本轮增量|完整BFF系数；在当前切片中的标量限制；费米在壳密度和对称标架联络项的具体抵消；切换后投影的领圈公式|
|验证|原参数/全部费米矩阵上的系数差分，原方程jet，参考和共同输出的精确校准|
|未增加|新耦合、新记录材料、新认知公理、平稳真空或瞬时谱重置|

历史去重：600/664已有几何来源，735已有联合Ward，791已有一圈均值来源，743—751已有有限图实际信号。此处补的是它们到当前受约束连续输出量的具体连接，不重新证明这些一般结果。382—386、425、522—523及384已消去的旧条件保持。

## 2. 必须澄清：线性切片，不是把参考非线性地冻结

785实际取固定背景线性算符L_R，对完整规范变换后的字段差施加切片条件。由该轮解的唯一性，

$$
L_Rw=0\quad\Longrightarrow\quad\chi(w)=0,\qquad
\mathcal W(w)=w,\qquad S_{\rm red}(w)=S_{\rm cl}(B_0+w).
\tag{1}
$$

最后一式在已经采用780标架截面、去掉纯规范坐标的变量中理解。**在q=0的物理切片上，不应再人为加一套非线性规范位移。** 从原任意坐标v回到w时，旧非线性变换及源伙伴仍须按789—796保留。

773参考为X=(h,s,(∇h)²,(∇s)²)，所以沿线性物理玻色方向w有

$$
X'w=0,\qquad \delta h=\delta s=0,\qquad
\delta F=0,\quad F=2-|\phi|^2/6.
\tag{2}
$$

这是在背景处的切向条件。它不等于X(B₀+w)=X(B₀)的非线性条件，也不把高阶径向变化全部删掉。内部切片使用H、DμH和颜色曲率的联合特征；不能将它擅自改为δH=0的单位规范。Higgs角向分量仍可能进入顶点。

早期工作稿的“固定h/s”仅应按式(2)理解。其参考jet结果保持，但非零度规密度不是响应非零的见证。

## 3. 保原Nambu计数的完整费米二次密度

用Ψ表示780截面中的**未作空间半密度缩放**的原Nambu旋量；原32CAR在这里写成64个自对偶分量，不是64个独立物种。与730的Cauchy半密度核按背景密度共同换写，原W不重选。令v=√|det g|，e=Ê⁻¹，

$$
\mathsf C^\mu=e_0{}^\mu I+\sum_{a=1}^3e_a{}^\mu\Gamma_a,\qquad
\Omega_\mu=\rho_{\rm spin}(\omega_{{\rm LC},\mu})+
\rho_{\rm int}(A_\mu),\qquad
\mathsf B(\phi)=\sum_{a=1}^5 N_a\frac{\phi^a}{\sqrt F}.
\tag{3}
$$

Γₐ和Nₐ是原604/598矩阵；所有复Y参数保留。内部连接取反Hermitian表示，原带色、电弱和sterile块均按原表示。实际Cᵘ与内部表示对易；spin连接不要求是普通Euclidean意义的反Hermitian矩阵。

在固定局部旋量平凡化中，原Hermitian费米密度可以写为

$$
S_F=\frac12\int d^4x\,v\,\ell_F,\qquad
\ell_F=\frac i2\left[\Psi^\dagger\mathsf C^\mu\partial_\mu\Psi
-(\partial_\mu\Psi)^\dagger\mathsf C^\mu\Psi\right]
 +\Psi^\dagger\mathsf Z\Psi,\qquad
\mathsf Z=\frac i2\left(\mathsf C^\mu\Omega_\mu-
\Omega_\mu^\dagger\mathsf C^\mu\right)-\mathsf B.
\tag{4}
$$

整体1/2只补偿Nambu重复，不能删去；所有双线性按原自对偶关系解释。相应普通四分量Dirac或原Weyl表达只是同一二次型的字典，不新增物种。这里没有加入独立扭率或四费米项。

成熟文献接口仍用[Zahn，1210.4031，§2—3](https://arxiv.org/html/1210.4031)：带规范/Yukawa背景的局部Dirac/Wick构造。其§4的特定应力守恒结论不直接覆盖本模型全部变化质量；本轮方程由原完整二次作用求导，联合Ward复用735，而非扩张该文定理的前提。

## 4. 显式一阶变分：全部同阶项先保留

设物理玻色方向为w=(kμν,aμ,ϕ)，其中k=δg。所有系数在B₀评价。780的平方根截面给

$$
\nu=\frac{\delta v}{v}=\frac12g^{\mu\nu}k_{\mu\nu},\qquad
\delta E=\frac12 E_0g_0^{-1}k,\qquad
\delta\mathsf C^\mu=-\frac12 k^\mu{}_{\nu}\mathsf C^\nu,\qquad
\delta\mathsf B=\sum_aN_a\left(
\frac{\varphi^a}{\sqrt F}+\frac{\phi^a(\phi\cdot\varphi)}{6F^{3/2}}\right).
\tag{5}
$$

在BFF阶只需这份标架的一阶导数，但它与旧动态spin字典相连，不能先冻结标架。一般背景的联络导数可直接从

$$
\omega_\mu=E\Gamma_\mu E^{-1}-(\partial_\mu E)E^{-1},\qquad
\delta\Gamma^\rho{}_{\mu\nu}
=\frac12g^{\rho\sigma}
(\nabla_\mu k_{\sigma\nu}+\nabla_\nu k_{\sigma\mu}-\nabla_\sigma k_{\mu\nu}),\qquad
\delta\Omega_\mu=\rho_{\rm spin}(\delta\omega_\mu)+\rho_{\rm int}(a_\mu)
\tag{6}
$$

求得。Γμ是Levi-Civita坐标连接矩阵，不是式(3)的空间Clifford矩阵Γₐ。

取原实紧支开关η，不对以下密度做未经配套的分部积分。原物理三价BFF项为

$$
V_{BFF}=\frac12\int d^4x\,\eta v\left\{
\nu\ell_F+
\frac i2\left[\Psi^\dagger\delta\mathsf C^\mu\partial_\mu\Psi
-(\partial_\mu\Psi)^\dagger\delta\mathsf C^\mu\Psi\right]
 +\Psi^\dagger\delta\mathsf Z\Psi\right\},
\qquad
\delta\mathsf Z=\frac i2\left(
\delta\mathsf C^\mu\Omega_\mu-\Omega_\mu^\dagger\delta\mathsf C^\mu
+\mathsf C^\mu\delta\Omega_\mu-\delta\Omega_\mu^\dagger\mathsf C^\mu
\right)-\delta\mathsf B .
\tag{7}
$$

式(7)是原对象的完整局部系数。794—795中与这份物理代表恰当等价的原实紧支完成、全部旧源伙伴仍保留，不能把任意原坐标下的裸截断改称该完成。该阶的其他原量子项按798计数只有线性来源及Wick线性项；它们未从作用删除，只是在第6节的特定居中混合量中抵消。

式(2)令显式singlet质量变化为零，且δB=Nₐϕᵃ/√F；Higgs角向项未一般消失。反之，不能把式(7)的完整度规响应等同于单独ν times 质量。

## 5. 两个确切抵消，及一个保留的非零主部

### 5.1 对称标架中的spin联络变分

在一点取背景正规标架，δωμab由∇k的反对称组合给出。Clifford反对易式使对称费米主部中的这项只取δωcab的完全反对称部分；k对称令它为零。换回任意固定背景标架，身份保持：

$$
\sum_\mu\left(\mathsf C^\mu\delta\Omega_\mu^{\rm spin}
 -(\delta\Omega_\mu^{\rm spin})^\dagger\mathsf C^\mu\right)=0.
\tag{8}
$$

这没有令原背景spin连接为零：式(7)中的δC·Ω和背景协变导数仍保留。若换成带相对Lorentz转动的截面，还须共同运输旋量，不能只拿式(8)删去它的补偿项。

### 5.2 体积项在所求交叉收缩中消失

将式(4)极化为ℓ_F(u,z)。原自由费米方程及其伴随给

$$
i\mathsf C^\mu\partial_\mu z+
\frac i{2v}\partial_\mu(v\mathsf C^\mu)z+\mathsf Zz=0,
\qquad
\ell_F(u,z)=0\quad\text{当u、z均为原自由解}.
\tag{9}
$$

证明是将方程左乘u†、伴随方程右乘z后相加除以2；密度/主部导数项抵消。这适用于全部原变化质量和规范连接，不依赖常质量或平稳背景。

有限光滑Cauchy记录与一个局部费米二次式的连通Wick收缩，只留下两条被记录模式涂抹的自由双解。故式(9)直接用于这两个交叉腿；局部减除所加的c数在协方差中消掉。因此Cov_F(p,ℓ_F(x))及反序协方差为零。**这不允许在离壳作用、其他插入或更高阶中删除νℓ_F。** 时间序场方程接触也没有被宣布为零；这里使用的是原W的普通有序跨收缩。

### 5.3 仍保留的物理耦合

在上述交叉收缩的意义下，式(7)化成

$$
V_{BFF}\ \simeq_F\ \frac12\int d^4x\,\eta v
\left[-\frac i4 k^\mu{}_{\nu}\,
\Psi^\dagger\mathsf C^\nu\overleftrightarrow{D_\mu}\Psi
+i\Psi^\dagger\mathsf C^\mu\rho_{\rm int}(a_\mu)\Psi
-\Psi^\dagger\delta\mathsf B\Psi\right].
\tag{10}
$$

D含原背景spin及内部连接；≃_F只指费米在壳交叉系数，不是完整相互作用泛函相等。符号和1/2沿式(4)固定。

可明确排除“当前切片让全部BFF项为零”。dh类时，局部正交标架可使dh沿时间轴，ds位于时间—x平面。取只有k_yz=k_zy=a的剪切，标量和内部连接变化为零。两条参考范数的一阶变化及全部原内部特征变化均为零，故这是ker L_R中的非零方向。它无迹，且δCʸ=−aΓ_z/2。选原Γ_z的+1本征旋量，在该点令∂yΨ=iκΨ；其时间导数由完整原费米方程确定。式(10)的密度中随κ的斜率为a/4（含Nambu权），不能被零阶质量、规范项抵消为恒零。

这只证明顶点在物理切片上不是零泛函，即使费米腿满足原方程也如此。**该任意玻色剪切尚未满足全部耦合玻色方程；也未同原W、η和固定记录作完整配对。** 不据此宣布连续读数非零。

## 6. 投影、切换和共同输出：现在可以实际代入的系数

原费米作用自身的联合规范不变性给一个有用的检查。取任意两条自由费米解，令j_F为其全部玻色背景变分的Euler来源，密度/参数伴随沿785。原Noether身份及费米方程给R₀†j_F=0，故

$$
\Pi^\dagger j_F=j_F,\qquad
\Pi^\dagger(\eta j_F)=\eta j_F-
L_R^\dagger[R_0^\dagger,\eta]j_F,
\qquad\Pi=1-R_0L_R.
\tag{11}
$$

第一式说明完整补偿没有删除原物理来源；第二式说明有限切换不能不经处理就移过投影。领圈项须按原紧支代表保留。若先将含∂w的密度作分部积分，η的导数也必须一并计入；不能把式(7)中的裸截断Euler来源直接误写为ηj_F。两种写法的边界字典应共同使用。

将式(7)，或其交叉收缩等价式(10)，展开为原受约束玻色jet Bα与原费米二次式Jα之和。对801的同一输出Φ(f)和p=e₀₀，第一相互作用系数为

$$
C_1(f,p)=i\int d\mu\sum_\alpha\left[
W_B(\Phi(f),B_\alpha)\,\operatorname{Cov}_F(p,J_\alpha)
-W_B(B_\alpha,\Phi(f))\,\operatorname{Cov}_F(J_\alpha,p)
\right].
\tag{12}
$$

此处Bα/Jα的系数已经由式(3)—(10)固定，不能另挑方便质量顶点。原首阶均值位移、纯玻色三价和原线性来源在本式居中相减后抵消；费米奇线性项由宇称排除。原完整BV恰当项由既有Ward处理。更高ε阶没有这些简化的自动保证。

式(12)同时变换两个输出。只计算δΦ而冻结δp会漏掉同阶项，甚至留下不应有的虚部。已有精确校准给单边13/750+4i/125，记录修正−4i/125，共同结果13/750；这是符号/共轭检查，不是原连续模型的数值预测。

## 7. 复算与准确验收

两组联合诊断，包含工作稿中此前未计数的校准，累计由3575到3577；子身份数量不算独立研究发现。

1. 原32CAR/64Nambu全部质量和主部：12个局部jet校准完整作用差分与式(7)，最大差约1.04×10⁻¹¹；质量导数约5.79×10⁻¹²；式(8)约2.38×10⁻¹⁶。漏掉密度、metric或质量的负对照均有明确非零误差。数值连接采样原电弱表示；对颜色表示的推广由同一线性变分公式承担，不声称对所有背景做了数值扫描。
2. 10个极化原费米方程jet的式(9)残差约1.16×10⁻¹⁶；剪切斜率0.05与含Nambu权的解析值一致。并重现两个既有工作probe：参考切片严格零残差、双输出关联的有理数抵消。

所有jet仅检验局部系数；没有把它们称为753的已求解时空，没有替换原W，也没有在该数值层计算一个完整物理三点分布。局部微分身份的解析证明与模拟范围分列。

**完成项：** 原完整BFF到当前物理切片/记录关联的显式映射，特定交叉密度和spin变分的抵消，以及不能漏掉的几何/内部/角向质量系数。

**未完成项：** 同一连续W与固定模式下的C₁总非零/总零判定；相应可区分信号；自主准备/末读/能源；有限耦合和图连续/跨尺度共同实现。它们不由本轮校准通过率解决。原四维、群谱、Einstein作用和参数仍是输入，统一目标保持。

接[803入口](803/drafts/STATUS.md)：把式(10)—(12)变成原自对偶CAR协方差上的具体有序配对；先核原纯准自由参考的选择规则及可探测条件，再决定最低非零阶，不重做一般Kubo公式或有限图配对。
