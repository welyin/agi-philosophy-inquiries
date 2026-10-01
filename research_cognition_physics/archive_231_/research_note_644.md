# 第644轮：同一Gauss区域态的费米复制、二阶熵与几何来源

日期：2026-10-01。接[643](research_note_643.md)。[代码](joint_region_replica_source.py)、[结果](joint_region_replica_source_results.json)、[核验](research_round_644_checks.json)、[条件账](unified_physics_condition_ledger_644.md)。三组检查、十四式；主代理审查，无新增独立代理审查。

## 1. 本轮共同接口

643给原完整Gauss热态的有序物质消元，但尚未把这份历史接到区域熵。**本轮将617同一规范切分、原全部CAR、643共同热历史及原几何来源接到二阶Rényi熵上。** 其复制闭合必须同时保规范边界与费米符号；来源须插入同一完整H，不能给熵另外配一份区域热Hamiltonian。

该结论适用于既定有限图、固定区域切法、正背景几何和原完整正常Gibbs态。解析式覆盖全部非线性玻色、规范和CAR；数值热态是明确冻结玻色的两节点64模式诊断，另有642的真实Gauss波包见证，二者不混称。没有得到面积律、Newton常数或GR。

## 2. 历史与成熟工具

[307](research_note_307.md)已有受限自由费米区域计算；[616](research_note_616.md)已区分普通热迹与宇称超迹；[617](research_note_617.md)给原商群切分J；[636](research_note_636.md)区分边界态、电来源和熵；[637](research_note_637.md)已证明原全态的区域熵有限。此处不将这些结论重新计作发现。

[Casini—Huerta的综述](https://arxiv.org/html/0905.2562)给自由场区域熵及复制方法的标准背景；[Klich的配对迹公式](https://arxiv.org/html/1403.7824)用于有限CAR权重。下面先从普通Hilbert迹自行确定复制符号，再接入原Gauss与来源，不凭“graded swap”名称省略符号核对。成熟复制恒等式本身不是本轮增量；增量是原完整区域、共同参考和来源的具体相容连接。

|层次|地位|
|---|---|
|认知动机|同一区域的物质、关联信息及几何反作用须来自同一态|
|继承输入|原完整H、Gauss、既定图、区域J、原CAR次序、β及正γ|
|解析|完整物理热态的带符号复制、桥可积与同一来源导数|
|数值|全部原两节点64模式的条件热态；真实Gauss波包的边界／宇称见证|
|边界|不把Rényi指数2当指数1；未建立连续尺度或几何面积熵等价|

## 3. 必须先固定原区域对象

以617的J切开跨区域原群链路，(Jf)(U_A,U_B)=f(U_AU_B)，原节点CAR不增加。设Kβ为原Gauss热核在扩展区域空间中的零延拓：

$$
K_\beta=J P_Ge^{-\beta H}J^\dagger,\qquad
Z=\operatorname{Tr}K_\beta,\quad
\rho=K_\beta/Z,\quad\rho_A=\operatorname{Tr}_B\rho,\qquad
R_2(A)=\operatorname{Tr}\rho_A^2,\quad S_2(A)=-\log R_2(A).
\tag{1}
$$

这里Kβ支撑于原顶点Gauss及切口匹配空间；未把扩展空间上的冗余态加入配分函数。J与γ无关。ρ_A是这个明确的扩展区域约化态，不是另选电中心／磁中心后不加说明的熵。

普通区域Hilbert交换Σ_A给

$$
Z_{2,A}=\operatorname{Tr}[(K_\beta\otimes K_\beta)\Sigma_A],\qquad
R_2(A)=Z_{2,A}/Z^2,\qquad 0<R_2(A)\le1.
\tag{2}
$$

纯度正性来自真实密度算符，不来自逐路径非负。617没有使Σ_A与两个副本的Gauss投影对易；两个P_G与全部切口条件必须保留，不能在交换前改成区域独立准备。

## 4. 将Hilbert交换转换成正确的CAR插入

原H、P_G、J都保总费米宇称，因此ρ_A与P_A=(−1)^N_A对易。对两个副本的单粒子空间，在A上定义旋转，在B上取身份：

$$
r_A\big|_{A_1\oplus A_2}=
\begin{pmatrix}0&-I_A\\I_A&0\end{pmatrix},\qquad
r_A\big|_{B_1\oplus B_2}=I,\qquad
\det r_A=1,\qquad \mathcal T_A=\Sigma_{A,\rm bos}\,\Gamma(r_A).
\tag{3}
$$

Γ为同一外代数表示。对同一区域的齐次宇称基|i〉、|j〉，Γ(r_A)的交换系数为(−1)^{N_j+N_iN_j}。ρ_A的非零矩阵元只连接同宇称，因此该系数在参与普通平方迹的项中为+1。玻色指标照常交换。这给

$$
\operatorname{Tr}[(\rho\otimes\rho)\mathcal T_A]
=\operatorname{Tr}\rho_A^2,
\qquad
\operatorname{Tr}[(\rho\otimes\rho)\Gamma(r_A^{\rm unsigned})]
=\operatorname{Tr}(P_A\rho_A^2)
\quad\text{（纯CAR情形）}.
\tag{4}
$$

unsigned把式(3)的−I改成+I。含玻色时第二式也需同样的玻色交换插入。**普通Hilbert交换Σ_A没有错；错的是把它直接认作无额外符号的单粒子换位之第二量子化。** 这是算符表示的区别，不是两种物理熵。

式(4)对任意偶正常态成立，不要求Gaussian。两个副本按分次CAR组合；其各自偶密度的乘积限制到A后，恰是ρ_A的偶乘积态，故跨B的Jordan–Wigner排序不另添可调符号。这也适用于核的偶矩阵块，足以接入643逐路径因子。

## 5. 接入同一有序物质权重和Gauss闭合

每个副本都有643的原Gauss闭合g₁、g₂和费米路径矩阵𝒰₁、𝒰₂。玻色区域交换重新连接两个桥的A端点，B端点仍留在各副本；切边先按J的原乘法拉回。条件费米因子成为

$$
\mathcal W_{2,A}=
\operatorname{Tr}_{\mathcal F\widehat\otimes\mathcal F}
\left[\Gamma(r_A)\Gamma(R_F(g_1)\oplus R_F(g_2))
(\mathcal U_1\widehat\otimes\mathcal U_2)\right].
\tag{5}
$$

单粒子式按643的Nambu加倍和正确模式重排，令𝕄₁₂为两条原有序历史及各自g的直和，𝕣_A为r_A的Nambu作用，则

$$
\mathcal W_{2,A}^{\,2}=\det(I_{4n}+\mathsf r_A\mathsf M_{12}).
\tag{6}
$$

与643一样，真正的迹由原Fock提升确定，不能用平方关系任意选支。两个副本各自保原商群约束；r_A通常不与其局部规范闭合交换，式(5)的相乘次序与桥端点须一起保留。这里不存在独立拟合的“区域熵权重”。

绝对可积也能从643继承，而不是只写形式路径积分。以f(x)=k^(1/2)_cut,β(x,x)为正标量支配核对角，x=(x_A,x_B)，J保证∫f=Tr exp[−β(T+W/2)]有限。交换π_A保持两个扩展区域的乘积测度，标量核Cauchy–Schwarz给

$$
\int dx\,dy\,
\sqrt{f(x)f(y)f(y_A,x_B)f(x_A,y_B)}
\le\left(\int f(x)dx\right)^2<\infty.
\tag{7}
$$

再乘d_F²e^{2βC}即控制带Γ(r_A)与两个Γ(g)的绝对路径权。规范平均的Haar质量为1，归一切口冗余亦为1。因而643热桥、区域复制及Gauss平均可在此界下共同使用，未假设区域独立动力学。

## 6. 区域熵变分必须使用同一个完整来源

保持区域节点集合、J、β固定，a为623允许的几何参数。令G_a=∂aH，A₀=H+c≥1。623给G_a A₀⁻¹及其伴随有界；603给全部热能源矩，故Duhamel插入在迹范数中可积：

$$
K_a'=-\int_0^\beta d\tau\,
J P_G e^{-(\beta-\tau)H}G_ae^{-\tau H}J^\dagger.
\tag{8}
$$

例如τ≤β/2时，将e^{−(β−τ)H}G写为e^{−(β−τ)H}A₀(A₀⁻¹G)，前者迹范数统一受热能源矩控制；另半区间从右侧同理。这比仅有形式期望有限更强，保证以下有界交换插入可微。

$$
Z_a'=\operatorname{Tr}K_a'=-\beta Z\langle G_a\rangle,
\qquad
Z_{2,A,a}'=\operatorname{Tr}[(K_a'\otimes K+K\otimes K_a')\mathcal T_A],
\qquad
\partial_aS_2=-\frac{Z_{2,A,a}'}{Z_{2,A}}-2\beta\langle G_a\rangle.
\tag{9}
$$

最后一项的负号来自+2 log Z；不能把归一化来源丢掉。一般G与H、𝒯_A均不对易，式(8)的虚时间插入不能擅自改成区域独立热平均。完整G含标量、规范非对角电项、CAR和几何系数来源，沿643共同保留。

若区域边界或J随几何移动，需另加其导数，当前固定区域结论不覆盖它。式(9)是二阶Rényi熵的来源合同，不是Einstein方程。

## 7. 原两个节点的全部64模式验证

数值冻结φ₁=(0,.7,0,0,.45)、φ₂=(0,.5,0,0,−.35)，链路取原群身份。原Y、F和Majorana均不变；在598允许的跳跃型中**指定**有界系数t=.17，未声称这一数值已被认知原则决定。条件单粒子矩阵是

$$
h=\begin{pmatrix}h(\phi_1)&tI_{32}\\tI_{32}&h(\phi_2)\end{pmatrix},\qquad
\Delta=\operatorname{diag}(\Delta(\phi_1),\Delta(\phi_2)),\qquad
\rho_F=\frac{e^{-\beta B(h,\Delta)}}{\operatorname{Tr}e^{-\beta B(h,\Delta)}},\quad\beta=1.7.
\tag{10}
$$

A为第一个节点的全部32模式。该冻结配置下有准确的上夸克6重、下夸克6重、带电轻子2重及中性配对1重因子；逐节点计数6×2+6×2+2×2+4=32。两节点中性因子使用256维Fock矩阵，其余使用16维；没有只取一个中性因子代表全物质。

独立Nambu相关矩阵C=(I+e^{βℍ})⁻¹限制到A的湮灭与产生模式，给

$$
S_{2,A}=-\tfrac12\operatorname{tr}
\log[C_A^2+(I-C_A)^2]
=2\log Z_F-\tfrac12\log\det(I+\mathsf r_A\mathsf M_{12}).
\tag{11}
$$

这里最后一式选真实正纯度对应支；数值交叉比较实际Fock约化、Nambu相关矩阵及256维复制行列式。三路同得S₂=18.965143611975，最大差低于4×10⁻¹⁵，纯度约5.80153319×10⁻⁹。

原上夸克两节点实际Fock因子的区域纯度为.380045716092；Γ(r_A)给同值，而无符号单粒子交换的Fock提升给−.232956959658，恰等于Tr(P_Aρ_A²)。该数值负值不是密度不正，而是计算了另一个算符。

对共同条件lapse B→eᵃB，因同一条件H与其自身对易，原归一态的导数为

$$
\rho_F'=-\beta(B-\langle B\rangle)\rho_F,\qquad
\rho_{F,A}'=\operatorname{Tr}_B\rho_F',\qquad
S_{2,A}'=-\frac{2\operatorname{Tr}(\rho_{F,A}\rho_{F,A}')}
{\operatorname{Tr}\rho_{F,A}^2}.
\tag{12}
$$

所有因子共同给S₂′=−4.947048015442。独立差分步长4×10⁻⁴、2×10⁻⁴的误差为4.47×10⁻⁹、1.14×10⁻⁹。删除节点间跳跃、另用各点热态时，区域S₂变为18.870098007864，相差.095045604111。这是所列条件物质模型的差别，不是完整Gauss热积分的物理数值。

## 8. 原真实Gauss态的符号见证

为避免全部符号检查仅停留在冻结条件态，沿642的正常有限能源Ψ_R：节点1一个R粒子、节点2满海中的R空穴，标量为原不变紧支撑波包，链路系数R(A)†R(B)/√d。将连接节点1的边按617切开。区域A的节点1及半边可组成d个正交向量a_l，对侧为正交b_l，因单位表示矩阵逐点满足行／列正交，

$$
J\Psi_R=\frac1{\sqrt{d_R}}\sum_{l=1}^{d_R}a_l\otimes b_l,\qquad
P_Aa_l=-a_l,\qquad
R_2(A)=1/d_R,\quad S_2(A)=\log d_R.
\tag{13}
$$

其他恒定链路及共同标量因子不增加Schmidt秩。每个a_l的费米数为1，因此无符号CAR交换给

$$
\langle\Gamma(r_A^{\rm unsigned})\rangle_{\rho^{\otimes2}}
=-1/d_R,\qquad
\langle\mathcal T_A\rangle_{\rho^{\otimes2}}=1/d_R
\quad\text{（含共同玻色交换）}.
\tag{14}
$$

R=ν、L、u、Q时d=1、2、3、6。原ν的d=1已经给−1与+1之别；原Q同时保d=6规范边界权和奇物质宇称。它们是完整Gauss态，但不是Gibbs态或H不变部门。使用旧态来验证新复制接口，不把旧Gauss存在性另算新发现。

## 9. 关闭的接口与仍开放的物理问题

C02区域切法、C13二阶区域熵、C14规范边界、C15原CAR、C19同一热参考和C22来源可以共用式(1)—(9)，无需另拟合熵的物质组成或来源。普通热参考、费米符号、规范拼接和全部能源必须同步保留。

这没有把[635](research_note_635.md)的连续局部UV几何熵等同于当前有限图S₂：635研究复制参数α=1的一阶变化，本轮指数为2；还缺同一态／尺度的连续匹配、复制数延拓和有限项。没有从本轮推出面积律、固定面积熵系数或Newton匹配。

三组复算均通过，历史冻结证据保留。接[645入口](round645_drafts/STATUS.md)：核真正连接有限图区域态与连续几何熵的条件；若只能重述复制技巧，转向总账的图—连续手征或动态几何共同输入，不继续扫描小型热因子的参数。
