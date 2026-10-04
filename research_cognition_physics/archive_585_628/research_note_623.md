# 第623轮：原完整量子过程的共同算符域与真实热响应

日期：2026-10-01。接[622](research_note_622.md)，回查[574](../archive_554_584/research_note_574.md)、[579](../archive_554_584/research_note_579.md)、[589](research_note_589.md)、[598](research_note_598.md)、[603](research_note_603.md)。[代码](623/joint_operator_domain_completion.py)、[结果](623/joint_operator_domain_completion_results.json)、[核验](623/research_round_623_checks.json)、[条件账](623/unified_physics_condition_ledger_623.md)。三组复算、十八式；主代理审查，无新增独立代理审查。

## 1. 新接通的条件

622证明了规整二阶来源系数收敛，保留“系数是否来自原无限维完整过程”的缺口。本轮针对原具体Hamiltonian证明：**正几何系数及原物质耦合的指定参数族具有共同算符域；同一全Gauss热态中的实际双历史过程具有二阶展开，其系数正是622的噪声、迟致及接触项。** 不必另加一个可解环境，也不必把动态标量冻结成质量矩阵。

这里的“完整”指598的固定有限图玻色、规范链路及有限CAR总系统，包含所有目标坐标的量子变化；不含动态引力度规。二阶结论是沿指定有限条源历史的标量过程展开，未证明传播子在Hilbert空间算符范数中二阶可微。

|层次|地位与范围|
|---|---|
|认知动机|态、演化、物质、几何来源和涨落必须共同成立|
|继承输入|原有限图、商群、目标K、正束缚势、CAR及原量子排序|
|本轮参数范围|固定M、原位L与u；外部正空间γ及原保Gauss耦合在紧参数集内C²变化；时间历史C¹|
|解析增量|共同算符域、来源算符界、瞬时热二点对象、实际CTP二阶展开|
|数值核对|原目标势的Laplace公式、含边界的图范数身份、原边势及32模式质量的增长|
|未消除输入|背景几何、参考热态、有限图、原物种与参数；未完成连续、动态几何或GR|

## 2. 对接已有定理及去重

[Shubin定理1.1](https://arxiv.org/pdf/math/0007019)用于完整配置流形上的光滑非负**标量**Schrödinger算子。先只处理动能加原位势，矩阵费米项随后用本轮算符相对界加入，不直接将标量定理套在矩阵势上。

[Schmid—Griesemer对Kato演化定理的说明及命题2.1](https://arxiv.org/html/1203.4700)给共同稠密域上C¹生成元的演化框架。下面验证原模型的共同域与参数正则性，再使用该定理；不是仅因“H自伴”就断言含时问题已经解决。

589只证明共同形式域，598只给费米项的形式扰动，603的热矩还未保证来源作用后是Hilbert向量，621—622则控制规整与涂抹二点对象。本轮补的就是这些接口；没有把成熟定理自身计为新发现。

## 3. 原势的一个全域表达

在每个H⁵节点取双曲面空间坐标x，而非引入新的物理空间维数：

$$
x=\phi/\sqrt F,\qquad z=\sqrt{M/F}=\sqrt{1+|x|^2/6},\qquad
a=\sum_{j=1}^4x_j^2,\quad b=x_5^2,\quad t=a+b,\qquad
\phi=x\sqrt{M/(1+t/6)}.
\tag{1}
$$

原U恰好成为关于a、b的二次多项式。令1=(1,1)ᵀ、y=(a,b)ᵀ：

$$
P=I-\frac{u1^T}{6M},\qquad \delta=Py-u/M,\qquad
U=\frac14\delta^TL\delta,\qquad
Q=\nabla_y^2U=\frac12P^TLP,\quad \ell=\nabla_yU(0)=-\frac12P^TL(u/M).
\tag{2}
$$

由双曲面诱导度量，grad x_i·grad x_j=δ_ij+x_ix_j/6且Δx_i=5x_i/6。因此

$$
\Delta_{\mathcal K}U=(8+2a)U_a+(2+2b)U_b
+(4a+2a^2/3)U_{aa}+(4b+2b^2/3)U_{bb}+(4ab/3)U_{ab}.
\tag{3}
$$

这包括原Higgs角分布，未把原势换成径向势。设l=max|ℓ_i|、h=max|Q_ij|，逐项得|ΔU|≤10l+(2l+14h)t+(8h/3)t²≤(11l+29h/3)(1+t²)。又因原L正定、u分量正且Σu<6M，取λ₀=detL/trL、η=1−Σu/(6M)>0、b₀=Σu/M：

$$
U\ge\frac{\lambda_0}{8}(\eta t-b_0)^2,\qquad
t\ge t_0:=2b_0/\eta\Longrightarrow U\ge\kappa t^2,\quad
\kappa=\lambda_0\eta^2/32,
\qquad |\Delta_{\mathcal K}U|\le C_\Delta(1+U).
\tag{4}
$$

显式可取C_Δ=(11l+29h/3)max(1+t₀²,1,κ⁻¹)。这是全域解析界；有限点采样只检查公式实现。原参数给C_Δ≈603.269，未优化常数。

## 4. 动能与原位势可以分别控制

记T_γ为589全部节点和电场动能，W_γ=Σ_i w_iU_i，先定义H₀,γ=T_γ+W_γ。令D_γ=−T_γ为配置Laplace算子；它不同于物理空间上的Laplace算子。系数w_i不依赖配置坐标，群导数不作用U，故

$$
\mathscr D_\gamma W_\gamma=\frac{\hbar^2}{2}\sum_i\Delta_{\mathcal K_i}U_i,
\qquad |\mathscr D_\gamma W_\gamma|\le C(1+W_\gamma).
\tag{5}
$$

这里及下文的常数在固定图、正几何紧集内统一。对光滑紧支撑复波函数，分部积分给

$$
\|H_{0,\gamma}\Psi\|^2=\|T_\gamma\Psi\|^2+\|W_\gamma\Psi\|^2
+2\int W_\gamma|\nabla\Psi|_\gamma^2d\mu
-\int(\mathscr D_\gamma W_\gamma)|\Psi|^2d\mu.
\tag{6}
$$

这里梯度范数包括动能中的系数。由式(5)及Young不等式，得到真正的算符图范数估计：

$$
\|T_\gamma\Psi\|^2+\tfrac12\|W_\gamma\Psi\|^2
\le\|H_{0,\gamma}\Psi\|^2+C\|\Psi\|^2,
\qquad D(H_{0,\gamma})=D(T_\gamma)\cap D(W_\gamma).
\tag{7}
$$

域等号不能仅由形式不等式推出。具体证明：配置乘积度量完整、W光滑非负，Shubin给H₀在C_c^∞上的本质自伴。沿其图范数核序列，式(7)不等式使TΨ与WΨ分别收敛，闭性给一个包含关系；反方向，对D(T)∩D(W)中的向量与C_c^∞配对，得到它属于最小算子的伴随域，本质自伴遂给等号。

## 5. 为什么变化γ仍是同一个算符域

取参考T_*=Σ_i(−Δ_Ki)+L_G，其中L_G是所有紧群链路的标准双不变Casimir之和。L_G与每个不变电动量对易，因而与589含非对角方向项的T_el也对易。把L_G分成有限维本征块，再对各节点Laplace作共同谱分解：T_*在每个纤维是标量，T_γ是与它一致可比的正矩阵。

因此不仅有形式比较，还有范数比较。几何和正电系数在紧正集合内时

$$
m\|T_*\Psi\|\le\|T_\gamma\Psi\|\le M_*\|T_*\Psi\|,
\quad D(T_\gamma)=D(T_*),\qquad
D(W_\gamma)=D(W_*),\quad W_*:=\sum_iU_i.
\tag{8}
$$

对参数导数同理：每个纤维上的一、二阶动能导数受C T_*控制。不要求不同γ的电动能彼此对易；所用的是它们共同与中心Casimir对易。此处也没有把节点Laplace的连续谱误作离散谱。

## 6. 边势与原费米物质不改变该域

设ρ_i=d_K(0,φ_i)。所有原规范变换固定目标原点，故d_K(φ_i,R(g)φ_j)≤ρ_i+ρ_j。589的交叉边势用正系数矩阵最大特征值控制；原磁势在紧群上有界。598则已给全有限Fock矩阵范数界。令R_λ=W_grad+W_mag+B，包括原复Dirac、Majorana及有界规范跳跃，则

$$
W_{\rm grad}\le C\sum_i\rho_i^2,\qquad
\|B(q)\|_{\mathcal F}\le C+cW_*^{1/4},\qquad
\|R_\lambda(q)\|\le\varepsilon W_\gamma(q)+C_\varepsilon
\quad\text{对任意 }\varepsilon>0.
\tag{9}
$$

最后一步来自式(4)及t=6sinh²(ρ/√6)，所以原位势在大ρ时至少指数增长。无须对目标Log作高阶配置导数。参数变化只进入边势的有限系数，不改变这一估计；一、二阶参数导数同样受控。有限Fock维数本身并未使质量算符有界。

式(9)给‖RΨ‖≤ε‖WΨ‖+C_ε‖Ψ‖，结合式(7)，R对H₀是无穷小算符相对有界。为说明所用自伴结论，可取z=±iy足够大，使‖R(H₀−z)⁻¹‖<1；Neumann逆保证两个半平面的满射，闭对称和遂自伴。图范数等价也保留原核。其闭形式与598同一形式一致，因此不是重新选择另一个自伴扩张。

得到共同域及参数正则性：

$$
H(\lambda)=H(\lambda)^*,\qquad
\mathcal D=\big(D(T_*)\cap D(W_*)\big)\otimes\mathcal F,
\qquad D(H(\lambda))=\mathcal D,\qquad
\lambda\longmapsto H(\lambda)\in C^2(\mathcal D,\mathcal H_F).
\tag{10}
$$

有限Fock张量指分量逐一处于同一域，D赋任一H的图范数。C²要求原系数族C²；这不是让任意不光滑源自动变得可微。Gauss投影与全部族交换，故限制到D∩H_phys仍成立。固定M、L、u是本轮范围；尤其没有把动态s当成可任意冻结的外部质量参数。正γ趋退化、图细化、动态量子几何和任意新增带荷源不在此定理中。

## 7. 同一热态的瞬时联合来源

基点H=H(0)，令A=H+c≥1。实际来源G_a=∂_aH、C_ab=∂_a∂_bH在共同域上对称，并满足

$$
\|G_a\Psi\|\le k_a\|A\Psi\|,\qquad
\|C_{ab}\Psi\|\le k_{ab}\|A\Psi\|,
\qquad \|G_a\rho^{1/2}\|_{\rm HS}^2\le k_a^2\langle A^2\rangle_\rho<\infty.
\tag{11}
$$

ρ是603原Gauss Gibbs态。来源二点对象定义为这些Hilbert–Schmidt向量的Gram，或等价的Σ_n p_n〈G_a(t)e_n,G_b(u)e_n〉；不要求未经证明的乘积算符G_a(t)G_b(u)有处处定义。其绝对值≤k_ak_b〈A²〉，在有限时间区间连续。

这将621—622在抽象形式类中的涂抹结果加强到**本轮原具体来源族**的瞬时二点对象。621的抽象反例缺少式(11)，所以没有被推翻。也未证明每个对称来源都本质自伴，或任意有限平均能量态都有有限噪声，更没有声称连续应力的点态UV奇性消失。

## 8. 实际演化的一阶与弱二阶展开

取ℏ=1书写演化，恢复时每个来源插入带1/ℏ。沿有限实C¹历史λ(t)，H_ε(t)=H(ελ(t))。式(10)和Kato共同域定理给唯一酉传播子，正反时间都保持D并在其图范数上有有限时间界。在B(D,H)中一致展开：

$$
H_\epsilon(t)=H+\epsilon V(t)+\frac{\epsilon^2}{2}W(t)+o(\epsilon^2),
\qquad V(t)=\lambda_a(t)G_a,\quad W(t)=\lambda_a(t)\lambda_b(t)C_{ab}.
\tag{12}
$$

固定初态，不随历史重新热化。对Ψ∈D，原实际传播子直接满足Duhamel身份

$$
U_\epsilon(T,0)-U_0(T,0)
=-i\int_0^T U_\epsilon(T,t)[H_\epsilon(t)-H]U_0(t,0)dt,
\qquad \|(U_\epsilon-U_0)\Psi\|\le C|\epsilon|\|A\Psi\|.
\tag{13}
$$

积分先作用在D向量上；U₀与A交换，外侧酉性即可给界，不需要先控制带源过程的二阶能量。稠密性将强收敛扩至全部Hilbert空间。对差商作强支配收敛得到

$$
U_1\Psi=-i\int_0^T U_0(T,t)V(t)U_0(t,0)\Psi\,dt,
\qquad (U_\epsilon-U_0)\Psi/\epsilon\longrightarrow U_1\Psi.
\tag{14}
$$

二阶不非法要求V U₀Ψ仍属于D。取另一个Χ∈D，将式(13)中外传播子的差移到bra侧；其伴随是反向演化，对Χ应用同一一阶结论即可。采用内积第二槽线性约定，二阶系数为

$$
\begin{aligned}
\mathcal U_2(\Xi,\Psi)
={}&-\frac i2\int_0^T\langle U_0(t,T)\Xi,W(t)U_0(t,0)\Psi\rangle dt\\
&-\int_{0<u<t<T}\langle V(t)U_0(t,T)\Xi,
U_0(t,u)V(u)U_0(u,0)\Psi\rangle\,du\,dt.
\end{aligned}
\tag{15}
$$

所有内积都只需每侧一个来源作用于D。由式(13)的正反向界和式(12)余项，有

$$
\frac{\langle\Xi,[U_\epsilon-U_0-\epsilon U_1]\Psi\rangle}{\epsilon^2}
\longrightarrow\mathcal U_2(\Xi,\Psi),\qquad
\left|\frac{\langle\Xi,[U_\epsilon-U_0-\epsilon U_1]\Psi\rangle}{\epsilon^2}\right|
\le C\|A\Xi\|\|A\Psi\|.
\tag{16}
$$

证明中的排序积分来自对U_ε(T,u)的一阶展开，插入时间t在u之后，故式(15)的次序与符号固定。这是共同域上的弱二阶Peano展开，不冒充D上的强二阶展开，也没有假定共同D(H²)。

## 9. 从弱展开到原热CTP

分别对两条历史应用上节。设He_n=E_ne_n、Ae_n=a_ne_n。实际Z_ε=Σ_n p_n〈U_−,εe_n,U_+,εe_n〉。两条一阶差的乘积强收敛；各单分支二阶弱项取Χ=U₀e_n∈D，式(16)给Ca_n²支配，603保证Σp_na_n²有限。因此热迹与二阶极限合法交换，得到

$$
\log\operatorname{Tr}_{\rm phys}(U_{-,\epsilon}^{\dagger}U_{+,\epsilon}\rho)
=\epsilon L_1+\epsilon^2L_2+o(\epsilon^2),\qquad
L_1=-i\int d_a(t)\langle G_a\rangle dt.
\tag{17}
$$

取Z=1附近log分支，d=λ_+−λ_-、c=(λ_++λ_-)/2。用式(11)的Gram定义N_ab=½〈{δG_a(t),δG_b(u)}〉及χᴿ_ab=iθ(t−u)〈[G_a(t),G_b(u)]〉，展开式(15)及跨分支项，得到

$$
L_2=-\frac12\iint d_a(t)N_{ab}(t,u)d_b(u)dtdu
+i\iint d_a(t)\chi^R_{ab}(t,u)c_b(u)dtdu
-i\int d_a(t)c_b(t)\langle C_{ab}\rangle dt.
\tag{18}
$$

式(11)使本轮二点积分直接绝对有界；同时，它们就是622固定来源规整系数的极限。故本轮在声明族内接通了“先规整求系数”和“实际完整过程展开”的两条路。不是给所有形式族的统一余项定理，也不是完整过程在所有源函数空间上的Fréchet C²定理；足以得到此处固定历史的实际二阶响应系数。

归一、交换分支共轭和Re L₂≤0继续使用同一态及同一过程。接触项必须是实际参数族的Hessian，不能独立补设。

## 10. 可复算证据与范围

三组检查均直接复用原参数，没有重新拟合：

1. 原球坐标势与式(2)在25点的最大相对差5.78×10⁻¹⁴；沿双曲测地线的独立二阶差分核式(3)，步长.01、.005、.0025的相对误差为6.08×10⁻⁵、1.52×10⁻⁵、3.80×10⁻⁶。全域界由式(4)证明，非由这些点证明。
2. 在原H⁵单节点因子保留全部角依赖，用复径向波函数积分式(6)。有限半径3.7的边界项为.5227358454，明确保留；48和96点求积的身份残差均≤5.62×10⁻¹⁶，两次总范数相对差2.25×10⁻¹⁵。删掉边界项会留下上述非零差额。证明用紧支撑核，此数值则用含显式边界的独立检验。
3. 两节点反向射线ρ从2至10，原边距离平方与4ρ²相符，最大相对差7.38×10⁻¹⁵。边势／原位势比8.82067降至.000167507；原32模式复质量的CAR范数上界／原位势比69.6001降至.00172134。这里的范数上界逐项用‖c†c‖≤1及配对项三角界，并非直接对2³²维Fock矩阵对角化。

这些检查核新恒等式和增长关系；全图共同域及真实热过程结论由第4—9节证明，不由单节点／两节点数值替代。未重复运行622的小Fock动力学并将其计为新轮次。

## 11. 联合条件的实际压缩

在本轮原固定图族中，原位势的同一个约束同时支持：物质稳定性、共同演化域、几何及物质来源、热态二点对象，以及实际二阶过程。C01、C08、C10、C17、C19、C22不再需要分别配一个演化或噪声环境。但该约束和模型输入没有因此成为从认知原则推出的新定理。

下一项[624](624/drafts/STATUS.md)回查实际记录、区域及尺度条件：检查已经存在的记录操作能否与这一真实来源过程共同使用，并把缺失的尺度／几何接口列清；不新设计认知硬件，也不继续单独优化本轮常数。只有新的共同身份或明确冲突才编号。连续手征、物理时空、动态引力、观测匹配与内部资源条件继续开放。
