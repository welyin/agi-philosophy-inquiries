# 第653轮：原手征辅助测度、规范状态与时间拼接的共同表示

日期：2026-10-02。接[652](research_note_652.md)，回查[614](research_note_614.md)—[616](research_note_616.md)、[628](research_note_628.md)—[629](research_note_629.md)、[646](research_note_646.md)及[653入口](round653_drafts/common_chiral_state_entry.md)。[代码](joint_chiral_character_state.py)、[结果](joint_chiral_character_state_results.json)、[核验](research_round_653_checks.json)、[条件账](unified_physics_condition_ledger_653.md)。四组检查、二十式；主代理审查，无新增独立代理审查。

## 1. 实际整合的条件与保留范围

原32模式CAR和手征候选虽有614的内部表示字典，尚不保证共享量子权重。本轮从615实际辅助积分出发，得到三项相互连接的结果：

1. 保持原电荷算符时，任意改变原32CAR状态都不能精确复现该候选的完整单点holonomy权重；不只是不允许直接删辅助项。
2. 同一单点配置扩展到原**整个非Abel商群**后，辅助权重是正状态的群特征函数，且有精确有限表示。
3. 在原m₀=1、零标量、空间平凡的时间链分支上，实际overlap辅助积分给出同一个正转移算符。任意有限时间长度的权重、规范投影及响应可沿这份表示拼接。

这把C14—C16规范／测度与C19状态、C04演化、C22来源在一个明确分支内接通。**没有**证明该分支Hamiltonian等于598原完整图H；空间传播、非零质量、一般配置及连续极限仍待接回。没有把辅助表示认作新增基本粒子，也不重新设计认知装置。

|层次|采用或证明|
|---|---|
|认知动机|同一物质的表示、状态、量子权重和拼接规则应共同相容|
|继承输入|614原G和物种；615的四维Euclidean处方、m₀=1、S⁹辅助测度、反周期时间和零标量|
|明确范围|空间各一格且空间链路为单位；任意有限时间格数及原G时间链路|
|解析新增|完整群的辅助积分；原CAR直接迹障碍；实际核的正时间拼接|
|精确代数|整数Laurent／Weyl证书、Beta矩及有限球谐谱；不靠浮点判断正性|
|数值角色|原手征矩阵与Pfaffian、非交换群元素和实际1—3格时间链独立核验|
|开放|原全图H匹配、空间局域性、非零质量、一般手征场及动态引力|

成熟候选来自[Kikukawa，§3.2、§3.6、§5.2](https://arxiv.org/html/1710.11618v3)：使用其辅助积分及子群限制，保留一般测度局域性的开放范围。[Weyl特征标公式](https://ocw.mit.edu/courses/18-745-lie-groups-and-lie-algebras-i-fall-2020/mit18_745_f20_lec26.pdf)和[Funk–Hecke公式，附录B](https://arxiv.org/pdf/2105.04504)是成熟数学工具；以下须对本项目实际核计算，不能以工具名称代替连接证明。

## 2. 原物种、辅助权重与严格受限的迹障碍

保留614的全局群和一代全左手表示：

$$
G=\frac{SU(3)\times SU(2)\times U(1)_Q}{\mathbb Z_6},\quad
V(g)=\operatorname{diag}(z^{-2}C,z^3W)\in SU(5),\quad
R(g)=\Lambda^{\rm even}V(g),\quad \dim R=16,\qquad
\mathcal F=\Lambda(R\otimes\mathbb C^2).
\tag{1}
$$

这里C、W分别为颜色和弱群矩阵，Q=6Y；内部旋量容器不被规范化为更大的物理群。沿615单点原超荷g(θ)，16个荷为1×6、−4×3、2×3、−3×2、6、0。记物理、辅助及完整条件因子为

$$
W(\theta)=\prod_{i=1}^{16}\cos^2(q_i\theta/2),\qquad
M(\theta)=\sum_{k=0}^{8}\frac{(k+1)(k+2)(9-k)}{990}
\cos^{2k}\theta\,\cos^{16-2k}(3\theta/2),\quad Z_f=WM.
\tag{2}
$$

对原32模式电荷Q₀，各内部q有两个spin占据。Σq=0，因此614粒子—空穴改写没有净常量荷偏移。精确有

$$
\operatorname{spec}Q_0\subset[-36,36],\qquad
2^{-32}\operatorname{Tr}_{\mathcal F}e^{i\theta Q_0}=W(\theta),\qquad
\operatorname{Tr}(\sigma e^{i\theta Q_0})
=\sum_{q=-36}^{36}\operatorname{Tr}(\sigma\Pi_q)e^{iq\theta}.
\tag{3}
$$

最后一式对任何θ无关正常正态σ成立，不要求[σ,Q₀]=0。可是式(2)的Fourier支撑分别扩至±24和±60，并且

$$
[e^{60i\theta}]Z_f=\frac1{55\,4^{23}}>0,\qquad
\sum_{|q|>36}[e^{iq\theta}]Z_f
=\frac{82745376001}{2902710697328640}>0.
\tag{4}
$$

所以保持原电荷、只改32模式状态或θ无关Hamiltonian的直接迹字典不可能精确成立。相同物种名称和反常抵消不弥补这个差别。

范围须严格：这是固定玻色／规范背景后的条件费米因子，尚未积分时间holonomy以施加Gauss。它不是原完整Gauss物理态上的外部化学势，也不排除扩展表示、保留有效测度或受控低能等价。615已有“删M会变来源”；本轮新增的是排除整个固定原电荷的状态重选类，并给下面的实际替代。

## 3. 原整个商群上的辅助积分

把V(g)对角化为diag(e^{iα₁},…,e^{iα₅})，选择Σα=0的提升。原16维偶外代数基A对应相位σ_A=Σ_{j∈A}α_j。复用615实际γ矩阵、固定γ₅基V_±、B=iγ₅C_D及对称内部矩阵T_a，定义

$$
u_A=e^{i\gamma_4(\pi+\sigma_A)/2}V_+,\qquad
v_A=e^{i\gamma_4(\pi+\sigma_A)/2}V_-,\qquad
\mathsf A_g(E)=u^T[T(E)\otimes B]u.
\tag{5}
$$

T_a非零配对的总内部权为±α_j；原反周期π产生同一符号。由原Clifford身份及固定Pfaffian定相得

$$
\mathsf A_g(E)=-T(C_\alpha E)\otimes\varepsilon,\quad
C_\alpha=\operatorname{diag}_{j=1}^{5}\big(\cos(\alpha_j/2)I_2\big),\qquad
\operatorname{pf}\mathsf A_g(E)=|C_\alpha E|^{16}.
\tag{6}
$$

均匀E∈S⁹在五个二维平面上的平方长度r_j服从Dirichlet(1,1,1,1,1)。其多项式矩使积分成为

$$
M(g)=\int_{S^9}\operatorname{pf}\mathsf A_g(E)dE
=\frac{h_8(a_1,\ldots,a_5)}{495},\qquad
a_j=\cos^2(\alpha_j/2),\qquad
h_8(a)=\sum_{k_1+\cdots+k_5=8}\prod_j a_j^{k_j}.
\tag{7}
$$

每项的积分系数为8!4!/12!=1/495。所有a_j≥0；若全部为零，V=−I₅，其行列式为−1，矛盾。因此该**整个单点平坦群族**上M>0且M(e)=1。不同相位提升改变余弦符号却不改变a_j；原测度的共轭不变性将对角证明扩展到G全部元素。并非一般多维规范配置上的非零性证明。

## 4. 正群表示的精确证书

令z_j=e^{iα_j}、Πz_j=1，取四个独立变量。式(7)可写为

$$
M(g)=\frac{N(z)}{D},\qquad
N(z)=h_8(2+z_1+z_1^{-1},\ldots,2+z_5+z_5^{-1}),\qquad
D=495\,4^8=32440320.
\tag{8}
$$

N有6661个非零整数Laurent单项式。对原子群Weyl群S₃×S₂取ρ=(2,1,0,1,0)，权以同时加(1,…,1)等价。代码按整数系数直接核

$$
N(z)\sum_w\operatorname{sgn}(w)z^{w\rho}
=\sum_\lambda c_\lambda\sum_w\operatorname{sgn}(w)z^{w(\lambda+\rho)},
\qquad c_\lambda\in\mathbb Z_{\ge0}.
\tag{9}
$$

不是只在几个群元素上拟合特征标。完整多项式两边相等，495个非零不可约类型的重数均严格正。若类型标为SU(3)的(a,b)、SU(2)最高权j和整数Q，则每个都满足中心下降条件，维数和为

$$
2a+4b+3j+Q=0\pmod6,\qquad
\sum_\lambda c_\lambda\frac{(a+1)(b+1)(a+b+2)(j+1)}2=D,
\qquad M(g)=D^{-1}\chi_{\mathcal K_{\rm int}}(g).
\tag{10}
$$

因此M是原G的正状态特征函数；这比M逐点为正更强。完整整数表在结果文件中。D是一个表示证书的维数，不是空间维数、基本粒子数量或最小实现规模。下一节会从原时间链得到更直接的加权表示，而不任意另造器件。

## 5. 从原overlap核得到真实的时间拼接

保留空间各一格、空间链路为单位和φ=0，但允许时间有任意N≥1格。定义原16维内部空间上的反周期协变移位S：S_{x,x+1}=R(g_x)，越过时间端点取负号。S严格酉。对m₀=1，原Wilson核满足

$$
X=-S^\dagger\otimes P_4^+-S\otimes P_4^-,\qquad
X^\dagger X=I,\qquad D_{\rm ov}=\frac{I+X}{2},\qquad
P_4^\pm=\frac{I\pm\gamma_4}{2}.
\tag{11}
$$

这一步使用原核，不重新指定一个方便的传播Hamiltonian。选γ₄正基W_+、W_−=γ₅W_+；相应手征基可直接取

$$
u=\frac{I\otimes W_++S^\dagger\otimes W_-}{\sqrt2},\qquad
v=\frac{I\otimes W_+-S^\dagger\otimes W_-}{\sqrt2},\qquad
u u^\dagger=\frac{I-\gamma_5X}{2}.
\tag{12}
$$

在相同固定定相下，辅助配对按时间点分块。定义原10维实表示O(g)通过T(O(g)E)=R(g)^*T(E)R(g)†；它是614容器向量表示在G的限制。辅助核精确给

$$
\operatorname{pf}\mathsf A[E,g]
=\prod_{x=0}^{N-1}
\left|\frac{E_x+O(g_x)E_{x+1}}2\right|^{16}
=\prod_x k(E_x,O(g_x)E_{x+1}),\qquad
k(E,F)=\left(\frac{1+E\cdot F}{2}\right)^8.
\tag{13}
$$

E_N=E_0；反周期负号在此双费米辅助项中平方消去，但物理行列式的反周期仍保留。证明使用W_+ᵀBW_+=W_-ᵀBW_-、交叉块零，及式(6)同一Clifford身份。不是将E固定在鞍点，也未以平均矩阵替代Pfaffian积分。

对归一S⁹测度定义转移算符及原群作用

$$
(\mathbb T f)(E)=\int k(E,F)f(F)dF,\qquad
(\mathscr U(g)f)(E)=f(O(g)^{-1}E),\qquad
M_N(g_0,\ldots,g_{N-1})
=\operatorname{Tr}\!\left(\mathbb T^N\mathscr U(g_0\cdots g_{N-1})\right).
\tag{14}
$$

k是内积非负整数幂的正系数和，因此为正定核；旋转不变性给[T,U]=0。式(14)直接来自实际多重积分的核乘法，对非交换时间链路仍成立。

## 6. 有限正支撑、同一Hamiltonian及规范投影

球谐次数ℓ≤8构成T的非零支撑。Funk–Hecke归约中t=E·F、u=(1+t)/2服从Beta(9/2,9/2)。用Gegenbauer Rodrigues式分部积分，或逐项精确Beta矩，得到

$$
\lambda_\ell=\frac{(9/2)_8}{(9)_8}
\frac{8!}{(8-\ell)!\,(17)_\ell}>0\quad(0\le\ell\le8),\qquad
\lambda_\ell=0\ (\ell>8),\qquad
d_\ell=\binom{\ell+9}{9}-\binom{\ell+7}{9}.
\tag{15}
$$

约定下标小于9时二项式为零，(a)_n为上升阶乘。边界因子保证分部积分无端点贡献，ℓ>8时8次多项式的高阶导数为零。核的正性及这些严格正本征值均有解析依据。

在非零转移支撑K_tr=⊕_{ℓ=0}⁸H_ℓ上，有

$$
\dim\mathcal K_{\rm tr}=35750,\qquad
\sum_{\ell=0}^8d_\ell\lambda_\ell=1,\qquad
\rho_{\rm aux}=\mathbb T|_{\mathcal K_{\rm tr}},\qquad
H_{\rm aux}=-a_\tau^{-1}\log\mathbb T\quad\text{在该支撑上}.
\tag{16}
$$

a_τ>0是给定时间格单位，采用ℏ=1的Euclidean约定。H_aux为有限自伴算符；球谐无穷尾是转移零空间，不能在完整L²(S⁹)上直接写有限log。此处的支撑约化由原积分核确定，未任意设置认知阈值，也未证明所有原辅助插入的可访问代数都已经重建。

这也给式(10)的第二份解析证书：Dλ_ℓ依次为735471、346104、134596、42504、10626、2024、276、24、1，全部为正整数。各球谐本就是原G表示，故DM=ΣDλ_ℓχ_{H_ℓ}是实际表示的特征标。代码再核它与式(8)的完整整数多项式相等。

物理负手征块在式(12)基中为(I−S†)/2，按每个spin取行列式；标准循环移位行列式给

$$
W_N^{\rm raw}(g_0,\ldots,g_{N-1})
=2^{-32N}\det(I+R(g)^{\dagger})^2
=2^{32-32N}W(g),\quad g=g_0\cdots g_{N-1},\quad
W(g)=2^{-32}\operatorname{Tr}_{\mathcal F}\Gamma(R(g)\otimes I_2).
\tag{17}
$$

detR=1，所以共轭／方向翻转在这个平方特征标中无影响。所有g无关归一与基定相均保留：N=1回到615的W。不能对每个N各自把完整热迹设成1而丢掉自由能。

于是这整条时间链的完整费米条件因子是同一个正Hamiltonian的扭曲热迹：

$$
Z_N^{\rm raw}(g)=\operatorname{Tr}_{\mathcal F\otimes\mathcal K_{\rm tr}}
\left[e^{-Na_\tau H_*}\,\Gamma(R(g)\otimes I_2)\otimes\mathscr U(g)\right],\qquad
H_*=\frac{32\log2}{a_\tau}I+I_{\mathcal F}\otimes H_{\rm aux}.
\tag{18}
$$

不是为每个温度重选一个状态：N仅改变同一个T的幂。H_*可给有限支撑上的真实酉群，但它只重建这条零标量、空间平凡的Euclidean分支；**不是598全图量子动力学的证明**。

由于两个表示都下降到原G，规范闭合可用同一个Haar投影，且与H_*对易：

$$
P_G=\int_G dg\,\Gamma(R(g)\otimes I_2)\otimes\mathscr U(g),\qquad
\int_G Z_N^{\rm raw}(g)dg=\operatorname{Tr}(P_Ge^{-Na_\tau H_*})>0.
\tag{19}
$$

严格正性因原Fock真空与辅助常数球谐给非零不变态。可进一步归一为匹配空间上的正常正态。这个投影作用于扩展表示，不能冒充原32CAR物理空间未变；式(3)—(4)已证明直接替换失败。此处群积分是时间规范闭合，不是把条件荷当作物理外部可调荷。

同一权重决定同一来源，例如沿原超荷族

$$
\partial_\theta\log Z_N^{\rm raw}(g(\theta))
=\partial_\theta\log W(g(\theta))
+\frac{\sum_{\ell=0}^{8}\lambda_\ell^N\partial_\theta\chi_{H_\ell}(g(\theta))}
{\sum_{\ell=0}^{8}\lambda_\ell^N\chi_{H_\ell}(g(\theta))}.
\tag{20}
$$

只在Z非零处使用对数。N=1精确回到615原完整来源；改变N后的响应也由同一谱确定，没有独立拟合。它是该分支的holonomy响应，不被命名为完整GR应力。

## 7. 四组核验及证据等级

1. **实际原矩阵：** 10组一般对角G配置、40个S⁹向量；式(6)配对误差≤3.16×10⁻¹⁶、Pfaffian误差≤2.56×10⁻¹⁵，物理行列式误差≤3.89×10⁻¹⁶；对615原族误差≤1.12×10⁻¹⁶。
2. **精确表示：** 6661项Laurent多项式与5940项交错式完全按整数核等；495个G类型全部正重数并过Z₆条件。独立32模式占据计数核式(3)—(4)。超出原电荷区间的精确总权约2.85062428292×10⁻⁵，不靠近机器精度的单个系数判断。
3. **共同状态／来源：** 7个非交换群元素的M、WM矩阵正定，最小本征值分别约.697697、.992987；这是实现检查，正性由式(14)—(16)证明。θ=.17的完整来源−14.7903978622与独立差分差约2.12×10⁻⁹。
4. **实际时间拼接：** 直接构造1、2、3格原overlap矩阵、手征基和Pfaffian；投影误差≤1.59×10⁻¹⁵，Pfaffian相对误差≤2.12×10⁻¹⁴，物理行列式相对误差≤8.24×10⁻¹⁴。9个本征值由独立Gegenbauer多项式与有理Beta矩逐项核对；球谐特征标与原N(z)按整数完全相等。平凡holonomy的M₂=6486285/2147483648≈.003020411823、M₃=2451339870343/87081320919859200≈2.8150007883×10⁻⁵，不能将M₁=1套给所有长度。

未数值积分高维配置来猜测正性，未建立巨大Fock矩阵，未做图像检查。有限检查不承担任意N的证明，后者来自式(11)—(18)的矩阵与核身份。

## 8. 减少的自由与下一项

原辅助测度一经固定，就同时决定正状态、允许的表示支撑、时间转移谱及holonomy来源；这几项不必各自增加假设。由此获得一条可以明确接入的测度—状态—时间接口，而非仅有群名称的一致。

代价和边界同样明确：辅助支撑不等于原CAR，原空间传播和非零标量质量被当前分支条件排除；恢复它们时，X†X=I及式(13)的分块未必保持。这个结果既不强迫宇宙新增35750种粒子，也不证明整个手征路线只有此实现。

[654入口](round654_drafts/STATUS.md)接回原空间传播与质量菜单，检查本次正转移表示的哪一部分可与既有全图H及其来源共用，哪一部分必须改变。优先核真实交叉条件，不继续优化辅助维数，不以这条限定分支替代完整统一目标。
