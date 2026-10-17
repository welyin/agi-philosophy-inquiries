# 第1238轮证明：原固定记录顶点的度规速度与CT继承资格

日期：2026-10-15。配套：[准入](admission.md)、[代码](check.py)、[结果](results.json)、[下一步](NEXT.md)。

## 1. 结论、量词与历史复用

本轮证明的是：**原852/860允许并已固定校准的空间应力读口，在573真实正则背景的开邻域内，完整离壳Dirac读写作用有非零的度规速度二次系数。** 因而不能原样把纯引力动量当作完整读写作用的动量，再直接调用1188的最小物质前提。

该命题只针对固定的旧作用表示和直接Legendre继承。它不排除分部积分、含边界的正则变换、辅助场表示或按有效阶重组；不排除CT、旧共同模型或整个GR目标。所用局部费米jet是检查离壳泛函的标签，不是一份新增物理准备。

主要复用：

- [573§5](../../archive_554_584/research_note_573.md)的同一真实约束解、类时h及满秩材料参考开集。
- [847§4](../../archive_819_853/research_note_847.md)原纯空间度规读口；[852§4、8](../../archive_819_853/research_note_852.md)固定校准及完整来源。
- [860§2、3、7](../../archive_854_872/research_note_860.md)的固定内部补偿及完整混合来源。
- [802§3—5](../../archive_797_818/research_note_802.md)的原32 CAR／64自对偶Nambu矩阵、密度、标架与Levi-Civita变分。
- [861§1、7](../../archive_854_872/research_note_861.md)只签收两导数零费米正则部门；[920§3](../../archive_894_922/research_note_920.md)已明确一般的完整动量警告。

因此，非零“原固定顶点系数”及其CT资格用途是净增；速度依赖的一般警告、共同来源原则、正则变换本身均不重新计为发现。

## 2. 原读口，及本轮固定的东西

用γ表示ADM空间度规，h仍为原Higgs径向材料。原材料四标签为

$$
X^A=(h,s,\mathscr A,\mathscr B),\qquad
\mathscr A=g^{\mu\nu}\partial_\mu h\partial_\nu h,
\quad \mathscr B=g^{\mu\nu}\partial_\mu s\partial_\nu s,
\quad Q^A{}_{\mu}=\partial_\mu X^A.
\tag{1}
$$

在一次选定的旧正则片上，记背景Jacobian为J、其逆为E。取847原z方向测试，其局部正交表达为正倍数的两份z协向量。852的材料校准是

$$
B=\frac{E^{T}b_z E}{\kappa v_i},\qquad
C_{\mu\nu}[\Phi]=B_{AB}(X)Q^A{}_{\mu}Q^B{}_{\nu},
\qquad
S_{\rm read}=-\frac{\lambda_i}{2}\int d^4x\,v\,q\,C_{\mu\nu}T_F^{\mu\nu}.
\tag{2}
$$

这里v为真实体积密度，v_i为852的平台解，两者不同；在选定内部窗口v_i>0，κ有限正。852的q=p；860的q=p−fε(X)，其中fε只校准一次。以下固定一个耦合标签λ_i，其他项与全部来源继续存在，不用不同标签的顶点作人为抵消。

“固定B”指固定材料函数，不指每次变更背景后再用新Jacobian重拟合它。下文只变更一点的度规时间jet，保g、所有材料场值及其一阶jet，所以X、B(X)、q、v在这个变化下都固定。该变化是离壳资格测试，不声称改变后的jet满足全套场方程或约束。

## 3. 原真实背景上确有非退化的速度方向

在初片取原N=1、shift=0，保持度规值固定，仅令

$$
\partial_0\gamma_{ij}\longmapsto
\partial_0\gamma_{ij}+\tau V_{ij},\qquad V_{ij}=V_{ji}.
\tag{3}
$$

物质一阶及二阶jet固定。由逆度规微分，只有Q的时间列发生下述变化：

$$
\delta Q^1{}_0=\delta Q^2{}_0=0,\quad
\delta Q^3{}_0=-H^iH^jV_{ij},\quad
\delta Q^4{}_0=-S^iS^jV_{ij},\qquad
H^i=\gamma^{ij}\partial_jh,\quad S^i=\gamma^{ij}\partial_js.
\tag{4}
$$

原573初片h、s均与z无关。在其已证满秩特殊线的开邻域内，可选cos(x+y)≠0且sin y≠0的点，于是空间dh沿(1,1,0)、空间ds沿(0,1,0)，两者非零且独立；h仍类时、J仍可逆。此处没有要求特殊线本身的某一导数额外非零，也没有借有限网格证明连续开集存在。

设w为E的z行，即w_A dX^A=dz。若w₃=w₄=0，则dz=w₁dh+w₂ds；右侧没有dz分量，矛盾。因此(w₃,w₄)≠0。两个独立空间向量的对称张量H⊗H和S⊗S线性独立，所以式(4)作为Sym²到两维空间的线性映射满射。存在一个有限V使

$$
a:=w_A\delta Q^A{}_0\ne0.
\tag{5}
$$

847的局部正交z可沿原坐标z方向选取；该方向与初片h、s梯度均正交。原读口的正系数连同固定归一记为c>0，于是B=cwwᵀ。固定一次该B后，式(3)—(5)精确给

$$
C(\tau)=c\,(dz+\tau a\,dt)\otimes(dz+\tau a\,dt),
\quad
C_{00}=ca^2\tau^2,\quad C_{0z}=ca\tau,\quad C_{zz}=c.
\tag{6}
$$

这不是仅说明存在任意新B：B仍由原847读口按852校准生成。选择正则点与局部窗口属于847原存在性分支；本轮不声称历史曾给每个波包固定数值中心，也不声称所有方向/所有点均非零。

## 4. 仅有C00非零还不够：保完整Dirac应力

原802固定对称标架截面上的未作空间半密度缩放作用为

$$
S_F=\frac12\int d^4x\,v\,\ell_F,\qquad
\ell_F=\frac i2\Psi^\dagger\mathsf C^\mu
\overleftrightarrow{\partial_\mu}\Psi
+\Psi^\dagger\left\{
\frac i2(\mathsf C^\mu\Omega_\mu-\Omega_\mu^\dagger\mathsf C^\mu)
-\mathsf B(\phi)\right\}\Psi.
\tag{7}
$$

Ω含全部Levi-Civita和原内部连接，B含原Dirac/Majorana质量及复Yukawa参数。整体1/2是原Nambu重复计数，不是新增64个独立物种。下面的极化读取沿802的合法费米二次型字典；不把一般复64列都解释为独立经典物质。

在同一点固定旋量值u及其余jet，选一个空间方向x并用实标签k表示

$$
\partial_x\Psi=(\partial_x\Psi)_{\rm fixed}+ik u.
\tag{8}
$$

取原Γx的+1极化，使局部正交坐标中的u†Cx u非零。在保存的共形坐标中该值为ψ⁻²。**不令费米方程随k重新决定时间导数**；这是离壳作用资格，不使用802§5.2的在壳密度抵消。原自对偶伴侣与半权沿802运输，未引入额外物种。

直接从完整Hilbert定义δS_F=(1/2)∫v bμν Tμν及802的 $\delta\mathsf C^\mu=-\tfrac12 b^\mu{}_{\nu}\mathsf C^\nu$（用同一g升第一指标）得到应力的费米导数部分：

$$
T^{\rho\sigma}_{\rm der}
=\frac12 g^{\rho\sigma}\ell_{\rm der}
-\frac i8\left(
\Psi^\dagger\mathsf C^\sigma\overleftrightarrow{\partial^{\rho}}\Psi
+\Psi^\dagger\mathsf C^\rho\overleftrightarrow{\partial^{\sigma}}\Psi
\right),
\qquad
\ell_{\rm der}=\frac i2\Psi^\dagger\mathsf C^\mu
\overleftrightarrow{\partial_\mu}\Psi .
\tag{9}
$$

所有质量、spin、内部连接项仍保留在T中，只是它们不含独立的k标签。对式(3)的变化，它们可依赖度规速度，因此不把完整T设成常数。由于g、v、Cμ的值固定，式(9)中的k系数不依赖τ。在N=1、shift=0下，

$$
\partial_k\ell_{\rm der}=-u^\dagger\mathsf C^x u,
\qquad
\partial_k T_F^{00}=\frac12u^\dagger\mathsf C^x u.
\tag{10}
$$

完整spin联络的变分也没有被遗漏：802式(8)在所用对称截面给

$$
\sum_\mu\left(\mathsf C^\mu\delta\Omega_\mu^{\rm spin}
-(\delta\Omega_\mu^{\rm spin})^\dagger\mathsf C^\mu\right)=0.
\tag{11}
$$

因此此最高费米导数系数可直接由Hilbert应力取得，不需将q乘到一个未经处理的积分变分后丢弃∂q领圈。原背景连接及δC·Ω继续存在；换标架时仍须原spin运输。式(11)也不允许把原连接一概设零。

## 5. 不可抵消的完整作用系数

式(6)中C00二次、C0z一次、Czz零次，而完整T的k系数均不含τ。故原读写密度中kτ²的唯一系数为

$$
[k\tau^2]\,\mathcal L_{\rm read}
=-\frac{\lambda_i q v}{4}\,ca^2\,
u^\dagger\mathsf C^x u\ne0
\qquad(\lambda_i q\ne0).
\tag{12}
$$

这里[kτ²]是有限jet多项式系数，不是概率或相关函数。质量与连接即便另含τ，也没有k；其它应力分量没有同一次数的τ²k，故不能抵消式(12)。这已经检查完整应力 contraction，不能简化成“C00非零所以作用非零”。

设v_V是沿V的度规速度坐标。在该固定有限jet表示中，读写项贡献的速度导数满足

$$
\partial_k\partial_{v_V}^2\mathcal L_{\rm read}
=-\frac{\lambda_i q v}{2}\,ca^2u^\dagger\mathsf C^x u\ne0.
\tag{13}
$$

于是直接的速度—动量关系含记录项。若仅写p_g=∂L/∂dotγ，则它不能仍等于纯引力的同一表达。这里没有构造完整高阶Ostrogradsky约化；原参考还含其它高阶jet，真正的全正则运输需另给。

当860背景q̄=0且Ψ̄=0时，式(12)所显示的qΨ²(dotδγ)²属于第五字段阶；它不改变原背景或二次自由Hessian。C0z还允许较低阶的线性速度项，所以不能说所有动量修正都从第五阶开始。更不能因自由Hessian不变就删掉实际保留阶的记录来源。

## 6. 对CT承接的准确影响

1188的候选类及λ、ab选择以同一独立π、最小物质项和速度字典为前提。本轮说明：**旧852/860完整作用在当前固定变量中，不能未经证明直接满足该前提。** 这不是证明其在一切变量中都不满足。

- 573的最小YM—sigma底座及861限定的零费米两导数部门仍可作为低阶承接对象。
- 803的半密度运输必须随来源共同运输，不能冻结它；本轮没有完成Dirac部门的完整正则消去。
- 741的高阶局部修正、852/860/869的记录顶点可留为声明有效阶的完整作用项；若如此，应保持旧899—917的同源/Ward/有限反作用预算，不把它们重新命名成度规速度无关物质。
- CT只选候选类内两个比率，不选择κ、Λ，也不凭此提供基片约束、全部保持或完整共同源解。CT尚未采用。

这关闭的是一项承接资格误认；未另证明某个有效截断的全部源误差可忽略，也未把新证明负担升级为整个认知纲领必须完成的高能模型。

## 7. 两组可复算核验

### 7.1 原保存背景jet的精确代数

直接读取907 `material_reference_jets_results.json` 的N=24、t=0、(0,π/2,π/4)旧X Jacobian，及900同N、同859扩展初片的ψ。该点数值资料不是573连续证明的替代。按原bz=2ψ⁴dz²及固定κv_i=1归一，由B=J⁻ᵀbzJ⁻¹算出

$$
B_{44}\simeq8.191865615568015\times10^7,\quad
\frac{\partial\dot X^4}{\partial\dot\gamma_{yy}}
\simeq-0.0008586713931688316,\quad
[\tau^2]C_{00}\simeq60.39998187502437.
\tag{14}
$$

将保存小数作为精确有理数输入，对五个非零/零速度标签核完整4×4张量，身份严格成立。保留保存的空间dh舍入小数，不擅自将其改成零。大逆Jacobian不是资源估计，也不据此宣称临界尺度。

预检中double矩阵回乘出现约7.5×10⁻¹¹抵消误差，未通过内部1×10⁻¹²阈值；正式检查改用精确Fraction验证同一输入代数，没有降低阈值或修改旧数据。它仍不证明原PDE近似的连续误差。

### 7.2 原完整64分量Dirac contraction

读取802原完整矩阵及同h、s值，保全部质量、非零内部连接、非零局部LC诊断jet、密度和Nambu半权。沿完整作用的实际度规变分作两种步长中心差分，再提取k和τ²混合系数。非零LC jet只是系数诊断，不冒充原背景发展资料。

用λ=2/7、q=3/11作为离壳系数标签，预期为

$$
-\frac{\lambda q\,\psi^6}{4\psi^2}\,
60.39998187502437 .
\tag{15}
$$

两个步长的完整计算均与式(15)相符；具体值和误差保存在results。删除密度变分会使本混合系数消失，提供同对象删弱检查。式(11)使用非零方向jet另核；不是先把spin项设零。原32 CAR和64自对偶记号保持，没有量子态、操作权限或物理准备的新断言。

默认check只读重放并逐项匹配结果；旧源码、报告、两份输入JSON在执行前后SHA256一致。正式报告与验收由主代理另行完成。
