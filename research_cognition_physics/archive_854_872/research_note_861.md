# 第861轮：共享材料坐标的经典约化Hamiltonian与尺度来源

日期：2026-10-06。接[860](research_note_860.md)、[入口](861/drafts/STATUS.md)及[前期推导](861/drafts/magnetic_reference_reduction_working.md)。配套：[代码](861/magnetic_reduced_hamiltonian.py) · [结果](861/magnetic_reduced_hamiltonian_results.json) · [核验](861/research_round_861_checks.json) · [复算入口](861/verify_round861.py)。

## 1. 本轮问题、结论和范围

859—860已让同一探针兼任局部参考和通信材料。本轮继续问：这些参考能否进入**原Einstein约束的实际经典约化**，同时保住Hamiltonian、材料体积和完整来源？这与858另加物理约束的体积匹配不同，不能混称为同一约化。

**结论：** 在既有正F、非零色磁场、类时h钟和满秩材料参考的局部片上，原两导数经典零费米部门具有显式的配置型规范固定及canonical变换。空间约束解出三份标签动量，原Hamiltonian约束再给一条实际正时钟分支。由此获得约化Hamiltonian和材料体积的共同字典；源变化必须同时运输被消去的动量与实际时钟速率。856型材料逆体积候选的一阶约化系数也可明确写出，但这不是原量子图—连续理论已经等价。

|层次|本轮采用或得到的内容|
|---|---|
|认知动机|承担参考、记录和资源的材料应属于同一过程|
|继承的建模输入|572/753的Einstein—规范—曲目标作用；852既有探针p；859的真实约束背景|
|解析增量|复合磁标签的完整辛运输、原四引力约束的局部约化、材料体积及来源拉回|
|另行选择|h时间与Y空间图册；材料计数密度η及胞元求积；若加入减除项，其系数/有限项仍是输入|
|校准|精确颜色与五个无迹度规方向；原859背景上的Hamiltonian根；来源链式变化|
|未宣称|完整高导数EFT的精确Hamiltonian、约化量子化、原图等价、579发散消除或全物理统一|

## 2. 复用与文献定位

直接复用[574](../archive_554_584/research_note_574.md)的曲目标逆度量，[753](../archive_742_763/research_note_753.md)及[859](research_note_859.md)的实际约束背景和参考秩，[856](research_note_856.md)的材料密度/完整变分合同。578—579的谱与尺度下界不重证。382—386、425、522—523空间接口保持原结论，本轮不增加维数证明，不重新引入384已消去的Lipschitz条件。

以材料时钟约化广义相对论是成熟路线；[Giesel与Thiemann，§1—2](https://arxiv.org/html/1206.3807v2)讨论不同材料模型的约化空间与物理Hamiltonian，并强调空间约束也必须处理。该文不能直接替我们证明由**原规范场与度规共同构成的磁标签**可用，更没有证明本项目的量子图等价。以下新增的是这些原字段的具体共同连接；不添加独立尘埃物种，也不引用尘埃模型的特殊Abelian密度结论。

860冻结稿式(12)中`epsilon_*Longrightarrow`的排版应读作epsilon_*后接逻辑推出符号；本轮登记勘误，不改冻结历史。

## 3. 配置型参考及规范保持

限于860的两导数S_base、经典零费米部门，保留p和所有原玻色场。取dh类时且向未来，h=τ。空间磁标签及规范条件为

$$
M=\frac12\gamma^{ik}\gamma^{jl}\langle F^c_{ij},F^c_{kl}\rangle,
\qquad Y=(s,p,M),\qquad
\chi=(h-\tau,Y^a-y^a),\qquad J^a{}_i=D_iY^a .
\tag{1}
$$

M在任意ADM相点是配置函数；只有h切片上才与859协变M_h相同。四参考满秩与dh类时给det J≠0。所有χ彼此对易且内部Gauss不变；内部Gauss仍保留，不冒称已把全部内部规范固定。

若C_perp[N]+C_i[β^i]为原约束，空间生成元选与普通Lie作用相差至多Gauss项的版本，则规范面上的矩阵及局部逆是

$$
\delta\chi^0=v_hN,\quad
\delta\chi^a=\mathcal B^a[N]+J^a{}_i\beta^i,\qquad
N=\frac{f^0}{v_h},\quad
\beta^i=(J^{-1})^i{}_a\left(f^a-\mathcal B^a[f^0/v_h]\right).
\tag{2}
$$

v_h为原C_perp[1]的实际h速度，正分支v_h>0。磁标签含空间导数，因此不能把全部B当乘法函数：

$$
\mathcal B^M[N]=N\mathcal B^M[1]
+2\langle F^{ij},v_{Aj}\rangle\partial_iN,
\qquad v_{Aj}=\frac{\delta C_\perp[1]}{\delta\Pi_A^j} .
\tag{3}
$$

度规正常变化已含在B^M[1]。保持规范用f⁰=1、f^a=0，得到实际lapse和shift；保紧支测试的局部微分逆不需要另做椭圆求逆。第一类约束在完整约束面弱对易，χχ块为零，所以任意两配置读数的Dirac括号为零。此结论不涉及它们与全部物质动量的括号，也不提供全局无Gribov图册。

## 4. 磁坐标必须共同运输电动量与度规动量

固定局部坐标体积，把空间度规写为

$$
\gamma_{ij}=e^{2r/3}\bar\gamma_{ij},\quad\det\bar\gamma=1,
\quad \sqrt\gamma=e^r,\quad
M=e^{-4r/3}\bar M,\quad
r=\frac34\log\frac{\bar M}{M},\quad
\bar M=\tfrac12\bar\gamma^{ik}\bar\gamma^{jl}\langle F_{ij},F_{kl}\rangle .
\tag{4}
$$

在M>0片，用M替换r，保留五个单位行列式度规坐标、原A及其余物质配置b。r与barγ含坐标密度权，不是普通标量/普通度规。由π_r=(2/3)π^{ij}γ_ij和紧支分部积分直接得

$$
\pi_M=-\frac{3\pi_r}{4M},\qquad
P_b^{\rm new}=P_b+\frac34(D_b\bar M)^*\left(\frac{\pi_r}{\bar M}\right),
\qquad
\Theta_{\rm old}=\int\left(\pi_M\delta M+P_b^{\rm new}\delta b+\cdots\right).
\tag{5}
$$

D_b为配置Fréchet导数，星号含空间分部积分；barγ变化限制无迹切空间。特别地

$$
\Pi_A^{{\rm new},j}=\Pi_A^j-\frac32D_i\left(\frac{\pi_r}{\bar M}\bar F^{ij}\right)
=\Pi_A^j+2D_i\left(\frac{M\pi_M}{\bar M}\bar F^{ij}\right).
\tag{6}
$$

因此新Π_A不是旧电场。全部原Hamiltonian、约束、物质读数与来源均须经式(4)—(6)逆运输。变换及生成函数内部规范不变，原Gauss结构随之保持；有限差分校准不代替连续规范论证。

## 5. 解空间约束后，原时钟约束仍可显式求根

记余下配置z为barγ、所有连接和Higgs角变量等，余下共轭动量为P_z；把h,s,p,M的动量记π_h,P_a。由于h和Y都为空间标量，而z的空间变换不依赖h或Y，canonical空间动量映射给

$$
C_i=\pi_hD_ih+P_aD_iY^a+C_i^{z},\qquad
C_i^{z}[\beta]=\int P_z\,\mathcal L_\beta z,
\qquad
P_a=-(J^{-1})^i{}_aC_i^{z}\quad (D_ih=0).
\tag{7}
$$

式中C_i^z按紧支分部积分定义密度；例如barγ的Lie导数要含密度权−2/3，连接可用普通Lie版本并保Gauss。C_i^z不含π_h或P_a。由式(5)等变辛势可证明式(7)，不能把新Π_A当原电场来猜这份约束。在Y=y规范内J=I，故P_a=−C_a^z；其空间导数也须保留。

原Higgs在h>0局部极坐标中，径向(h,s)与角动能正交；原正F目标的径向逆块为

$$
\mathcal F=2-\frac{h^2+s^2}{6}>0,\qquad
K^{hh}=\mathcal F(1-h^2/12),\quad
K^{hs}=-\mathcal F hs/12,\quad
K^{ss}=\mathcal F(1-s^2/12).
\tag{8}
$$

p为852独立中性标量。磁canonical变换不依赖h，且式(7)消去的P_a不依赖π_h。于是完整原C_perp先按式(4)—(7)拉回、再取h=τ,Y=y，得到

$$
\widetilde C_\perp=\alpha\pi_h^2+\beta\pi_h+\zeta,\qquad
\alpha=\frac{K^{hh}}{2\sqrt\gamma}>0,\quad
\beta=\frac{K^{hs}P_s}{\sqrt\gamma},\qquad
v_h=2\alpha\pi_h+\beta .
\tag{9}
$$

ζ**定义为完整原约束的上述拉回减掉已显示两项**：它保原引力动量与曲率、完整规范电磁能、P_s和P_p动能、Higgs角动能、全部空间梯度和原势，不是可任意选取的自由函数。旧Π_A逆运输含∂π_M，因此ζ可含余下动量的空间导数；但它不含π_h或∂π_h。这正是这里可以点式求时钟根的原因。

在继承的约束解附近，Δ=v_h²>0。选原向未来分支：

$$
\Delta=\beta^2-4\alpha\zeta>0,\qquad
\pi_h^*=\frac{-\beta+\sqrt\Delta}{2\alpha},\qquad
\mathcal H_{\rm phys}=-\pi_h^*,\qquad N=\Delta^{-1/2}.
\tag{10}
$$

证明物理Hamiltonian的符号：原canonical作用限制到χ=C=0后，标签不动而dot h=1，辛项成为∫(P_z dot z+π_h*)，因此H_phys=−π_h*。余下Gauss及其乘子继续保留。式(5)没有时间显含，故没有遗漏另一个时间生成项；H_phys本身可经h=τ显含τ。

这给局部约化作用及原解的规范呈现。它不声称物理Hamiltonian半有界、自伴或全局存在，也不把这个复合图册的高空间导数形式当作新PDE适定性证明。原局部解的存在继承既有协变方程；尚未证明新的相空间范数下任意数据均可全局推进。

## 6. 来源必须在消去的动量内继续求导

对仍在同一两导数正则类别、满足约束兼容性的作用参数/来源变化j，先运输全部作用，再解空间约束。固定剩余canonical变量时，隐函数求导给

$$
\delta\mathcal H_{\rm phys}
=\frac{1}{v_h}\left[
(\delta\alpha)(\pi_h^*)^2+(\delta\beta)\pi_h^*+\delta\zeta
\right],\qquad
\delta\widetilde C_\perp
=(\delta C_\perp)_{P_a}+D_{P_a}C_\perp[\delta P_a].
\tag{11}
$$

第二式在π_h固定时用，δP_a来自式(7)，也保源依赖的J、C_i^z和canonical变换。D可含空间导数；积分来源须连1/v_h一起作正确伴随分部积分。体积密度、参数及材料定位变化均包含在此链式导数，不能只保显式δ势。

精确校准中，沿一份声明的空间动量变化P_s'=2/9，漏掉隐式反力会使dH_phys/dj少−2/9。它只检查来源运输代数，不把随意外加j称为新协变理论。任意外源若破坏第一类约束，不满足本节前提。

## 7. 材料体积与逆体积候选的一阶共同连接

复用856定义，η(Y)d³Y为材料计数三形式。在h切片以及Y坐标内分别有

$$
n=\eta\frac{|\det D_iY^a|}{\sqrt\gamma},\quad
V_R=\frac{\sqrt\gamma}{\eta|\det D_iY^a|},\qquad
V_R\big|_{Y=y}=\frac{1}{\eta}\left(\frac{\bar M}{y^3}\right)^{3/4}.
\tag{12}
$$

因此这个经典约化片的体积确是配置读数；它仍随颜色配置和无迹度规变化。η及离散胞元选择未被坐标定理决定，原590独立控制代数亦未被证明等价。

为明确尺度接口，**另声明**856型协变候选I₁=−∫√−g n²，但把Y换成当前(s,p,M_h)。在dh类时、h=τ规范中，其拉回不含额外时间导数：M_h=M_Σ且n只含空间Jacobian。进一步Y=y，得到

$$
I_1\big|_{h=\tau,Y=y}=-\int d\tau d^3y\,N\mathcal U,
\qquad \mathcal U=\frac{\eta^2}{\sqrt\gamma}
=\eta^2\left(\frac{y^3}{\bar M}\right)^{3/4},\qquad
\delta\log\mathcal U=-\frac34\delta\log\bar M\quad(\delta Y=0).
\tag{13}
$$

这里的Y不是新独立物种。只取**沿原局部解的第一阶有效作用系数**：辅助变量的零阶驻值使其一阶改变量不贡献S₀的一阶变化，所以可将原实际N₀代入I₁。紧支变化下得

$$
\left.\frac{d\mathcal H_{\rm phys}}{d\lambda}\right|_0
=N_0\mathcal U,\qquad
\delta(N_0\mathcal U)=\mathcal U\,\delta N_0+N_0\,\delta\mathcal U.
\tag{14}
$$

这包含式(6)的电动量运输及式(10)时钟根的全部依赖。不是把N固定成1，也不是只将式(13)当一份没有材料来源的体积函数。若从原能量减去指定项，λ取相反号；沿856约定a=ℏ²/3，并非本轮推出其值。

式(14)是规范固定后的一阶变分系数：未证明将任意U直接添加到离规范面的C_perp还能保持原约束代数，也未给860全部高导数通信项的精确全阶canonical形式。高阶约束完成、量子测度和有限反项仍须共同检查。

还有一个有用的尺度判据。若**单独研究**保持α、β和空间约束不变、ζ→ζ+λU的约化根候选，则

$$
\mathcal H_{\rm phys}(\lambda)-\mathcal H_{\rm phys}(0)
=\lambda\frac{\mathcal U}{v_h}
+\lambda^2\frac{\alpha\mathcal U^2}{v_h^3}+O(\lambda^3),\qquad
\lambda_*=\frac{v_h^2}{4\alpha\mathcal U}>0.
\tag{15}
$$

λ=λ_*处两根并合、所选h钟速率为零；这仅是该固定剩余数据候选的钟图册端点，不是时空或整个理论奇点。负λ的减除方向不因此被禁止；其零点Taylor展开的最近分支点仍在正λ_*。

同一光滑材料几何上选η_δ=δ⁻³时，U_δ按δ⁻⁶增长，式(15)的解析半径按δ⁶缩小。这个**非一致展开判据**说明固定δ的一阶匹配不能直接交换为固定物理ℏ的δ→0结论。它不重证579，不否定所有重新组织的减除或耦合随尺度变化方案，也不证明这些方案已经成功。

## 8. 可复算核验

[前期颜色校准](861/magnetic_coordinate_canonical_probe.py)及其[结果](861/magnetic_coordinate_canonical_probe_results.json)保当时工作稿元数据：五点周期差分的四份变化严格保辛势，漏电动量平移则都失败；逆运输严格保电能。该离散泛函不是精确非阿贝尔格点理论。

本轮[新增代码](861/magnetic_reduced_hamiltonian.py)复算它，再补：

|核验|结果|
|---|---|
|五个无迹度规方向|精确有理数辛势差全部为0；材料体积与反项变化等大反号|
|完整来源链式运输|微分约束残差严格为0；漏空间动量反力缺陷−2/9|
|原859背景，16³/24³，ε=.02|Hamiltonian密度残差约−5.91×10⁻¹⁴、5.62×10⁻¹⁴|
|用式(10)恢复原π_h|误差8.05×10⁻¹²、7.65×10⁻¹²；判别式均正|
|实际钟速率及lapse|v_h约.00733879，N约136.262；不是单位lapse|
|源一阶有限差分|相对误差1.25×10⁻⁹；有限正/负候选变化的根残差小于2×10⁻²¹|

背景数值仅在旧初片与h法向一致的指定点核验；**没有声称将整个16³/24³网格都变换为Y图册**。全局网格不是本轮解析局部定理的前提。16³与24³是两个既有背景配点检查，不当作量子连续收敛证据。单位η的源数值是该局部坐标密度校准，未擅自改称最终材料网格的物理减除常数。

本轮计一组新增科学验证，正式861、累计3646。历史稿、旧代码和冻结文件全部保留；无图像检查、应用任务、目标或Git改动。

## 9. 下一项与失败判据

C09/C20/C22得到的是同一经典分支的配置参考—Hamiltonian—体积—来源连接。下一轮[862入口](862/drafts/STATUS.md)直接检验这份字典能否接574/598原图过程的**共同量子对象**，优先审计Gauss、材料边及受约束钟引起的混合动量和场依赖权；不把固定图单胞元下界原样移植为约化谱结论。

若两边采用的代数/态域、动量或计时不同，则必须构造运输或明确不匹配；不得靠替换变量名称签收。候选减除的有限耦合/尺度一致控制以及579剩余来源继续开放。整体认知—物理统一目标保持。
