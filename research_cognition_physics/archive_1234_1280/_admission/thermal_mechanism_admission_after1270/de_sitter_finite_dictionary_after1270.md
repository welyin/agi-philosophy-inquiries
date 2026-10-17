# de Sitter有限类光段：物理面积、完整改进源与模流权的共同字典

日期：2026-10-16。接[上一方向筛查](geometric_thermal_identity_candidate_after1270.md)。这是一次非编号、有限解析核对；使用线性Einstein方程，不是从认知原则生成GR。没有修改原阶段、目标、文献或冻结证据。

## 1. 准入问题及复用边界

[318](../../../archive_301_341/research_note_318.md)已证明平直仿射类光段上改进应力、Wald接触与有限端点须共同保留。[1250，§3—5](../../1250/proof.md)已在另一实际FLRW材料上区分物理仿射参数、面积权、热权及几何余项。这些工具直接复用，不重报一般分部积分。

本次只问：在给定de Sitter物理度规上，原文的共形模矢量与物理应力/面元配对后，究竟对应哪个面积泛函？这一对象和1250的实际材料、窗口不同。

核对来源为[D’Angelo等，PTEP 2024, 021A01，§2、§4](https://doi.org/10.1093/ptep/ptae003)。原文式(16)—(21)涉及共形重标定、来源和面积符号；下面明确固定物理定义。结果只排除**未经运输便把这些符号全部认作同一物理量**的用法，不据此判定整篇论文错误。

所有方程取相干激发振幅 $\eta$ 的 $O(\eta^2)$ 项：物质场为 $O(\eta)$、面积/曲率响应为 $O(\eta^2)$。背景/参考的几何独立异常项按同一背景参考相减。没有证明非线性、改变背景后真空极化或完整重整化余项的有限数值界。

## 2. 固定物理度规、面元和参数

取 $H>0$，背景

$$
ds^2=\Omega^2\left[-du\,dv+\frac{(v-u)^2}{4}d\Omega_2^2\right],
\qquad \Omega=\frac{-2}{H(u+v)} .
\tag{1}
$$

在 $v=0$、$u\in[a,b]\Subset(-2\ell,0)$ 上，

$$
g_{uv}=-\frac{2}{H^2u^2},\quad
\Gamma^u{}_{uu}=-\frac2u,\quad
dA_0=H^{-2}d\Omega_2,\quad
(\mathrm dv)^\sharp=-\frac{H^2u^2}{2}\partial_u .
\tag{2}
$$

选择使正能流为正的面取向，物理类光面元为

$$
d\Sigma^a=(\partial_u)^a dA_0\,du .
\tag{3}
$$

例如直接用 $\sqrt{-g}=\Omega^4r^2\sin\theta/2$ 与 $g^{uv}=-2/\Omega^2$，便得到式(3)的绝对归一，$r=-u/2$。参数重标定时切向量和参数微元抵消，不能单改升指标后的法向量、丢掉横截面积因子。

$u$不是物理仿射参数。一份方便的仿射参数为

$$
\lambda=-1/u,\quad
\frac{du}{d\lambda}=u^2,\quad
k_\lambda=u^2\partial_u,\quad
\lambda_0=\frac1{2\ell}.
\tag{4}
$$

原文给出的模流矢量在该面上的限制，记作 $\zeta$，为

$$
\zeta=f(u)\partial_u,\qquad
f(u)=-\frac{\pi}{\ell}w(u),\quad
w(u)=u(u+2\ell).
\tag{5}
$$

同一矢量在物理仿射参数中恰为

$$
\zeta=2\pi(\lambda-\lambda_0)\partial_\lambda .
\tag{6}
$$

因此这个具体物理面上的权是线性仿射boost权；它不是把二次 $u$ 权原封移到物理仿射聚焦方程。

## 3. 从独立物理面积到正确有限能流

令 $q=\delta\log\sqrt{\gamma}$；在保持所选类光面/生成线识别的线性光学规范中，背景物理膨胀和剪切为零。以同一有限角片定义

$$
\delta A(u)=\int q(u,\Omega)\,dA_0,\qquad
J_{\rm H}(u)=\int\Delta T^{\rm phys}_{uu}(u,\Omega)\,dA_0,\qquad
L_u=\partial_u^2+\frac2u\partial_u .
\tag{7}
$$

非仿射Raychaudhuri与线性Einstein输入给

$$
L_uq=-8\pi G\,\Delta T^{\rm phys}_{uu},\qquad
L_u\delta A=-8\pi GJ_{\rm H}.
\tag{8}
$$

背景膨胀为零，使非仿射系数的一阶变化乘背景膨胀而消失；剪切平方和膨胀平方在此阶不贡献。这不授予非线性有限振幅恒等式。

按物理流定义的有限量是

$$
K_{\rm H}[a,b]
=\int_a^b f(u)J_{\rm H}(u)\,du
=\frac1{8G\ell}\int_a^b w(u)L_u\delta A(u)\,du .
\tag{9}
$$

关键恒等式是

$$
\int_a^b wL_u z\,du
=[wz'+2\ell z]_a^b .
\tag{10}
$$

于是得到完全保留有限端点的共同关系

$$
\boxed{K_{\rm H}[a,b]
=\frac{[w\delta A'+2\ell\delta A]_a^b}{8G\ell}.}
\tag{11}
$$

定义

$$
\mathfrak B_f[s]=\left[s-\frac{f}{2\pi}s'\right]_a^b,
\qquad S_A=\frac{\delta A}{4G},
\tag{12}
$$

则 $K_{\rm H}=\mathfrak B_f[S_A]$。在仿射图同一式为

$$
K_{\rm H}
=\frac1{4G}
\left[\delta A-(\lambda-\lambda_0)\partial_\lambda\delta A\right]_{\lambda(a)}^{\lambda(b)} .
\tag{13}
$$

它保留面积和端点膨胀，不是均匀平均面积。有限源离开并不能自动令两个端点面积或膨胀为零；本式不需要这种条件。

这里 $K_{\rm H}$首先是给定矢量与完整相减应力的物理能流。将它进一步等同某个正常相干态的相对熵，需要同一区域代数、参考、场与边界约定；有限子段本身不自动是一个完整模Hamiltonian。

## 4. 共形来源运输及原式(18)—(21)的比较

对四维共形标量，同背景参考相减后，按本次给定的字典

$$
\Delta T^{\rm phys}_{uu}=\Omega^{-2}\Delta T^{\rm flat}_{uu}
=\frac{H^2u^2}{4}\Delta T^{\rm flat}_{uu}.
\tag{14}
$$

结合物理面元，式(9)等于

$$
K_{\rm H}
=-\frac{\pi}{4\ell}\int_a^b\!\int
u^3(u+2\ell)\Delta T^{\rm flat}_{uu}\,d\Omega_2\,du .
\tag{15}
$$

这与把原文式(18)的 $T_{uu}$ **逐字解释为这里的平直协变分量**所得的
$-2\pi\ell^{-1}\int w\,\Delta T^{\rm flat}_{uu}\,d\Omega_2du$
不是同一权。若原符号已吸收共形/横截面密度或别的归一，必须先给该转换；不能将其视作无关常数。

为明确区分，而非争论记号，定义后一种逐字截断泛函 $\widetilde K_{18}$。在本物理面积图里，它满足

$$
\widetilde K_{18}
=\frac1{G\ell}\int_a^b
\left(1+\frac{2\ell}{u}\right)L_u\delta A\,du
\tag{16}
$$

以及

$$
\widetilde K_{18}
=\frac1{G\ell}\left\{
\left[\left(1+\frac{2\ell}{u}\right)\delta A'
+\left(\frac2u+\frac{6\ell}{u^2}\right)\delta A\right]_a^b
+\int_a^b\left(\frac2{u^2}+\frac{12\ell}{u^3}\right)\delta A\,du
\right\}.
\tag{17}
$$

因此连保留这份字面权的选择，也须用不同的面积体权和端点，不能直接得均匀面积平均。

原文式(20)采用的形式核在本有限段若直接代物理面积，会给

$$
\widetilde K_{20}^{\,\rm naive}
=\frac1{16G\ell}\int_a^b w\,\delta A''\,du
=\frac{[w\delta A'-w'\delta A]_a^b+2\int_a^b\delta A\,du}{16G\ell}.
\tag{18}
$$

原式(21)的均匀平均正来自这一形式积分。它不等于式(11)，也没有纳入式(2)的非仿射项。下节固定波形给出有限不等的精确见证。所否决的是这个**额外的直接物理字典**；不声称已审完原文所有可能的重标定解释或其余结论。

## 5. 同一非最小耦合的Wald接触

令 $\Phi$为物理共形场、$\xi_c=1/6$。在本阶，完整Hilbert类光来源为

$$
\Delta T^{\rm phys}_{uu}
=(\partial_u\Phi)^2-\xi_cL_u(\Phi^2).
\tag{19}
$$

一般曲率耦合还有 $\xi_cG_{uu}\Phi^2$；本背景 $G_{uu}=0$，其几何扰动贡献为更高阶，不是任意删源。令

$$
F_A(u)=\int\Phi^2dA_0,\quad
J_{\rm can}(u)=\int(\partial_u\Phi)^2dA_0,\quad
S_\xi=-2\pi\xi_cF_A .
\tag{20}
$$

因此

$$
J_{\rm H}=J_{\rm can}-\xi_cL_uF_A,\qquad
L_u\left(\delta A-8\pi G\xi_cF_A\right)=-8\pi GJ_{\rm can}.
\tag{21}
$$

对同一个有限流权定义 $K_{\rm can}=\int fJ_{\rm can}du$，得到

$$
\boxed{K_{\rm can}=\mathfrak B_f[S_A+S_\xi],\qquad
K_{\rm H}=K_{\rm can}-\mathfrak B_f[S_\xi].}
\tag{22}
$$

等价的接触差是

$$
K_{\rm H}-K_{\rm can}
=\frac{\pi\xi_c}{\ell}[wF_A'+2\ell F_A]_a^b .
\tag{23}
$$

这使完整来源与同耦合的Wald截面荷共同运输。Wald接触不是另一份非负von Neumann熵，也不由此确定实际热参考。对一般有限端点不能删式(23)。

与318的区别明确：318输入平直仿射线性权；本次从给定物理de Sitter度规和共形模矢量证明了它在**物理仿射图**中的线性形式，并给非仿射图的正确端点字典。

## 6. 一份有限精确函数见证

取 $H=G=1,\ \ell=2,\ a=-2,\ b=-1$，按单位背景立体角计算。取

$$
\Phi(u)=\eta g(u),\qquad
g(u)=
\begin{cases}
(u+2)^3(u+1)^3,&-2\le u\le-1,\\
0,&\text{其余}.
\end{cases}
\tag{24}
$$

这是非零紧支 $C^2$类光数据；$g,g',g''$两端均零。一般推导对更光滑波形同样成立。当前只验证有限经典来源/线性聚焦字典；没有借该 $C^2$函数直接授予指定diamond的正常Hadamard相干态、完整紧支制备或实际端口。

令 $\delta A=8\pi\eta^2Q$，由

$$
L_uQ=-\left[(g')^2-\tfrac16L_u(g^2)\right],\qquad Q(a)=Q(b)=0
\tag{25}
$$

确定一份响应。代码用有理Laurent多项式精确构造，全部系数保存在结果。端点导数为

$$
Q'(a)=\frac3{5720},\qquad Q'(b)=-\frac9{10010}.
\tag{26}
$$

因而两端面积为零并没有使能流消失：

$$
\boxed{\frac{K_{\rm H}}{\pi\eta^2}
=\frac{K_{\rm can}}{\pi\eta^2}
=\frac{\mathfrak B_f[S_A]}{\pi\eta^2}
=\frac{12}{5005}>0.}
\tag{27}
$$

非零响应由端点膨胀承担，未假设未来事件视界条件。面积和完整曲率耦合系数在足够小的非零 $\eta$下保持正域；本次不以此宣称已控制真实非线性余项。

为检验非零接触，在同一波形、同一 $Q$上只取内部段 $[-7/4,-5/4]$，没有重新制备或改源：

$$
\frac{K_{\rm H}}{\pi\eta^2}=\frac{319119267}{167939932160},\quad
\frac{K_{\rm can}}{\pi\eta^2}=\frac{123681291}{83969966080},\quad
\frac{\mathfrak B_f[S_\xi]}{\pi\eta^2}=-\frac{14337}{33554432}.
\tag{28}
$$

三者精确满足式(22)。遗漏接触会留下这份严格非零有限差。

在原有限段上，式(18)的直接物理套用还满足

$$
-\frac3{2500}
<
\frac{\widetilde K_{20}^{\,\rm naive}-K_{\rm H}}{\pi\eta^2}
<
-\frac1{1000}.
\tag{29}
$$

这是严格有限删除：有理积分只含 $\log2$，代码用
$\log2=2\sum_{j\ge0}3^{-(2j+1)}/(2j+1)$及正尾上界认证，不以浮点数代证明。式(17)对另一字面泛函的恒等式也被精确核验。没有把有限窗差解释成整个PTEP理论、GR或认知纲领失败。

## 7. 复算与停止线

[代码](de_sitter_finite_dictionary_after1270.py)仅使用Python标准库Fraction，核验一份同源波形的：

- 共形完整应力运输、非仿射聚焦及正有限能流；
- 物理面积端点和非零Wald接触；
- 原文形式权若未经运输套入物理量时的有限失配。

[结果](de_sitter_finite_dictionary_after1270.json)首次排他创建；默认运行只读匹配。代码绑定本稿、旧318/1250和上一筛查，不新增研究计数。

这项净连接足以决定方向：**物理面积与完整来源可以有共同有限恒等式，但其权和端点须由物理参数运输得到；不能直接用均匀面积平均替代。** 尚未给实际热参考、测温/来源仪器或从认知关系生成几何。到此停止，不扩展装置、全部场论或非线性引力。

