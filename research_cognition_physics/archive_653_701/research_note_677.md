# 第677轮：全谱有理表示与完整物理来源平均的共同极限

日期：2026-10-02。接[676](research_note_676.md)、[已执行入口](677/drafts/rational_support_entry.md)。[代码](677/joint_rational_physical_limit.py)、[结果](677/joint_rational_physical_limit_results.json)、[核验](677/research_round_677_checks.json)、[条件总账](677/unified_physics_condition_ledger_677.md)。两组新数值核验、十六式；主代理审查，无独立代理审查。

## 1. 本轮补什么，哪些旧结果直接使用

676完成并开始677入口，属于有效进展。本轮核对README、direction、state、最新结果及发布后保护记录，未发现运行中的Python研究。当前目标仍是把既有物理分支放进同一可审核过程；先整合共同条件，认知系统设计后置。

用户提醒的旧空间结果已逐项回查：

|旧轮次|直接继承的结论|当前还须连接的对象|
|---|---|---|
|382、383|完整实际方向接口给三维上界；实际反向只需连续自由对合|同一物质过程的真实位置壳与仪器|
|384|可逆重定向与全族协变给下界，已消去额外Lipschitz／Hausdorff输入|不把这些已消去条件重新列为缺口|
|386|渐近位移连接真实邻域；一个有限尺度的严格误差裕量即可给维数证书|实际端点与当前过程的映射；不要求固定绝对可见度|
|425|局部紧端点群、一致半幅与成本收缩给NSS／Lie坐标|与386是替代连接路线，不要求同时满足；不加半幅同态或全局收缩自同构|
|522、523|明确精度任务下的热参考下界；实际方向仪器与条件性三维已有共同实现|把其参考态、仪器及尺度映射接到同一原物质过程|

**本轮不增加空间定理。** 新增连接是：659的符号极限不再只作用于某个固定初始帧，而是作用于673的全部固定尺寸原物理来源；再用669热矩和673零测度结论，把这个极限接入完整配置、辅助及双Haar平均。

截断overlap及其domain-wall关系是成熟方法，分别见[Neuberger原研究](https://arxiv.org/abs/hep-lat/9710089)、[Boriçi原研究](https://arxiv.org/abs/hep-lat/9912040)。本轮不把有理符号近似当新理论，也不从这些文献借来本项目原质量、来源与体减除的身份。

|层次|地位|
|---|---|
|认知动机|改变表示时，原观测、来源、参考及代价应同时保留|
|继承输入|原群／16通道／质量／S⁹测度／H_b；673指定两时间片未归一候选|
|额外表示选择|下面的全谱有理调节；不是新增认知公理|
|解析新增|固定有限盒完整来源泛函的加权L¹收敛，含有限质量来源导数|
|数值新增|原非平坦非零支撑上的二／四来源、质量响应、全矩阵收敛和有限调节的规范输送|
|未证明|有限调节RP、局部体测度身份、原H_F物理时间、连续量子几何|

## 2. 以全谱调节替代固定初始空间

沿676，X=D_W−I、H=γ₅X。0<a<1时第五方向转移T_a严格正，Shamir核为

$$
H_a=\gamma_5aX(2I+aX)^{-1}
=(I-T_a)(I+T_a)^{-1},\qquad
\sigma(H_a)\subset(-1,1).
\tag{1}
$$

对正整数L定义

$$
\varepsilon_{a,L}
=(I-T_a^L)(I+T_a^L)^{-1}
=\tanh(L\,\operatorname{artanh}H_a)
=\frac{(I+H_a)^L-(I-H_a)^L}{(I+H_a)^L+(I-H_a)^L}.
\tag{2}
$$

分式指右乘分母逆；此处各因子同为H_a函数，彼此交换。分母严格正。公式定义在全谱，包括H零模，不选初始J₋、不固定某个spin分区的负谱维数。第五坐标仍是辅助坐标。

记原H非零点的谱隙g=min|σ(H)|、h=‖H‖。m₀=1、d条单位规范移位时，X=(d−1)I−Σμ[(I−γμ)Tμ+(I+γμ)Tμ†]/2，每个括号除以2为酉，故h≤2d−1。当前保留两传播方向，h≤3；这不是空间维数结论。若ah/2<1，直接展开逆得

$$
\left\|\frac{2H_a}{a}-H\right\|
\le\delta_a:=\frac{ah^2}{2(1-ah/2)}.
\tag{3}
$$

当δ_a<g，B=2H_a/a的隙至少g−δ_a。对sign的对称实轴resolvent积分应用恒等式，差的被积范数被δ_a/[(g−δ_a)²+t²]控制，积分给‖sign B−sign H‖≤δ_a/(g−δ_a)。再接式(2)标量误差：

$$
\|\varepsilon_{a,L}-\operatorname{sign}H\|
\le 2e^{-2L b_a}+\frac{\delta_a}{g-\delta_a},
\qquad b_a=\min|\operatorname{artanh}\sigma(H_a)|
\ge\frac a2(g-\delta_a).
\tag{4}
$$

所以对于每个固定非零H背景，

$$
a\longrightarrow0,\quad aL\longrightarrow\infty
\quad\Longrightarrow\quad
\varepsilon_{a,L}\longrightarrow\operatorname{sign}H.
\tag{5}
$$

这是659已有共同极限条件的全谱使用，不重计为新极限规律。单独L→∞只得到sign H_a。入口中保持aL约20的两点不满足式(5)，其误差平台不是反例。

## 3. 同一固定尺寸来源及统一界

原单粒子维数n=2r，原Berezin变量ξ有2n个。定义软矩阵P_u=(I−ε)/2、P_v=(I+ε)/2。它们是互补正收缩；有限a,L一般不幂等。原固定J±、M(E)、bar-M和原反对称质量Pφ不改：

$$
D_\varepsilon=\tfrac12(I+\gamma_5\varepsilon),\qquad
\Phi_\varepsilon=
\begin{pmatrix}
J_-^\dagger P_v&0\\
-J_+^{\mathsf T}M(P_u+\tfrac12P_v)&J_+^{\mathsf T}
\end{pmatrix}.
\tag{6}
$$

$$
\mathcal N_{\varepsilon,\lambda}
=\begin{pmatrix}MQ_+&-D_\varepsilon^{\mathsf T}\\D_\varepsilon&\bar M Q_-\end{pmatrix}
+\lambda\Phi_\varepsilon^{\mathsf T}P_\phi\Phi_\varepsilon,\qquad
F_\varepsilon[j]=\int d\xi\,
e^{\xi^{\mathsf T}\mathcal N_{\varepsilon,\lambda}\xi/2+j^{\mathsf T}\Phi_\varepsilon\xi}.
\tag{7}
$$

这里沿用673的Berezin顺序及来源顺序。该定义在有限参数是一个正则多项式，**不直接套用673精确投影的三角解耦、675辅助因子分離或676零模计数**。仅在ε→sign H后恢复那些原身份。

因‖ε‖≤1、M酉、‖P_u+P_v/2‖≤1，用块范数矩阵[[1,0],[1,1]]，得与673相同的常数c=(1+√5)/2：

$$
\|D_\varepsilon\|\le1,\qquad
\|\Phi_\varepsilon\|\le c,\qquad
\|\mathcal N_{\varepsilon,\lambda}\|
\le2+|\lambda|c^2\|P_\phi\|.
\tag{8}
$$

这个界对全部配置、a及L成立，无谱隙分母。若ε与ε′都是Hermitian收缩，η=‖ε−ε′‖，由Φ的两个变化块分别为J₋†(ε−ε′)/2与J₊ᵀM(ε−ε′)/4可得

$$
\|\Phi_\varepsilon-\Phi_{\varepsilon'}\|\le\frac{\sqrt5}{4}\eta,\qquad
\|\mathcal N_{\varepsilon,\lambda}-\mathcal N_{\varepsilon',\lambda}\|
\le\left(\frac12+\frac{\sqrt5}{2}|\lambda|c\|P_\phi\|\right)\eta.
\tag{9}
$$

任意指定2k个物理来源列组成Z，按673同一增广Pfaffian约定，其系数为

$$
C_{\varepsilon,Z}=
\operatorname{Pf}\begin{pmatrix}
\mathcal N_{\varepsilon,\lambda}&\Phi_\varepsilon^{\mathsf T}Z\\
-Z^{\mathsf T}\Phi_\varepsilon&0
\end{pmatrix},\qquad
|C_{\varepsilon,Z}|
\le\bigl(2+|\lambda|c^2\|P_\phi\|+c\|Z\|\bigr)^{n+k}.
\tag{10}
$$

奇数来源系数为零。Pfaffian平方等于行列式，给出此范数界；Pfaffian的有限多项式性和式(9)给逐来源收敛。无需质量逆、传播子逆或除以配分函数，因此包括原零权重支撑。

## 4. 新连接：整个原平均的加权L¹极限

固定当前有限盒和原H_b热时τ>0，沿673全部原测度，记

$$
d\varpi_\tau=
k_\tau(q_i,a_0q_j)k_\tau(q_j,b_0q_i)\,
d\mu(q_i)d\mu(q_j)d\nu(E_i)d\nu(E_j)\,da_0\,db_0.
\tag{11}
$$

a₀、b₀为两边界群元，区别于调节参数a；没有删除任何一份规范输送。原热核正，因而dϖτ是正的有限控制测度，并非声称费米权重是概率。

673已证明：对每对空间配置，两个时间边界中的H零集Haar零测度。Fubini给H≠0为dϖτ几乎处处。669全部热矩与673的核Cauchy–Schwarz给

$$
\|P_\phi\|\le C[w(q_i)+w(q_j)],\qquad
\int [w(q_i)+w(q_j)]^p\,d\varpi_\tau<\infty
\quad(p<\infty).
\tag{12}
$$

对固定Z，式(10)提供不依赖a,L的可积控制函数；Z若是原场多项式增长且规范协变的有限来源，同理增大有限p即可。式(5)、(9)、(10)及支配收敛因此给本轮主要结论：

$$
\int|C_{\varepsilon_{a,L},Z}-C_{\operatorname{sign}H,Z}|\,d\varpi_\tau
\longrightarrow0,\qquad
\int C_{\varepsilon_{a,L},Z}\,d\varpi_\tau
\longrightarrow\int C_{\operatorname{sign}H,Z}\,d\varpi_\tau.
\tag{13}
$$

这是**包含全部原S⁹、双Haar与配置平均的未归一物理来源泛函**，不是对完整积分做了蒙特卡洛估计。有限源空间意味着所有来源系数同时有此结论。固定qᵢ、qⱼ的完整辅助／Haar核也逐点收敛；式(13)进一步给全配置加权L¹控制。没有要求排除近零谱配置，也没有要求每个扇区有相同秩。

允许有限质量来源θ，使Pφ(θ)=Pφ,0+ΣθαPφ,α，每个Pφ,α反对称且有式(12)的多项式增长界；其作用不改变H、Φ或控制热测度。对每个有限多重指标β和紧θ集合，有限多项式展开同样给

$$
\partial_\theta^\beta\int C_{\varepsilon_{a,L},Z}(\theta)\,d\varpi_\tau
\longrightarrow
\partial_\theta^\beta\int C_{\operatorname{sign}H,Z}(\theta)\,d\varpi_\tau,
\quad\text{在紧θ集合上一致。}
\tag{14}
$$

这不是对任意几何或规范来源交换微分：那些来源可能移动H零集、改变热测度及观测本身。也没有同时取τ→0、体积→∞或物理格距→0。其一致热矩、谱／局域性与源控制仍须单独建立。

676的72／56扇区在极限中的全部纯物理来源为零。式(13)保留这一支撑；不要求每个软矩阵提前精确为零，也不从近零浮点数推断其相位。

## 5. 规范、相位与正性的边界

全谱函数在规范变换R下满足ε′=RεR†。对原ξ和物理观测的单位矩阵Gξ、G_Y，与673同样有

$$
\Phi'_\varepsilon G_\xi=G_Y\Phi_\varepsilon,\quad
G_\xi^{\mathsf T}\mathcal N'_\varepsilon G_\xi=\mathcal N_\varepsilon,\quad
j'=G_Y^{-\mathsf T}j,\quad
\det G_\xi=1.
\tag{15}
$$

这些身份在有限a,L成立，保留完整复相位和原质量输送。规范物理插入仍需正确收缩与边界输送；不是把任意带电j直接叫作物理可观测量。

674的原极限标量／实质量来源反射身份继续可用；式(13)使其在积分极限中恢复。但本轮没有证明每个有限软调节的同一反射身份，更没有证明反射正性。完整平均收敛也不意味着其极限为非零正态：

$$
\text{完整未归一来源的受控极限}
\;\not\Longrightarrow\;
\text{反射正性、非零归一、原 }H_F\text{ 身份或物理时间。}
\tag{16}
$$

局部体动作、体行列式减除以及来源怎样共同给出式(7)，仍须实际检验。659已证明原来源依赖的参考因子不能随意丢掉，这一要求保留。

## 6. 可复算结果与限度

复用既有Python／NumPy。第一组采用673已经指定的原全群非平坦配置、四份不同原φ、原16通道及四份S⁹辅助场。新增检验全部512维N、物理Φ、二／四来源及质量响应；不是重新搜索空间维数或优化层数。固定L=16/a²：

|a|L|sign误差|Φ误差|N误差|
|---:|---:|---:|---:|---:|
|.125|1024|.0733696|.0366150|.0380686|
|.03125|16384|.0185112|.00923094|.00964596|
|.0078125|262144|.00465045|.00231833|.00242537|
|.001953125|4194304|.00116507|.000580807|.000607424|

二来源绝对误差从3.79×10⁻²¹降到5.08×10⁻²⁴，四来源从9.39×10⁻²³降到1.32×10⁻²⁵，原质量导数误差从4.00×10⁻²⁰降到5.44×10⁻²³。这一背景的原Pf约1.10×10⁻²¹，最末权重相对误差仍约6.70%；没有把矩阵的小误差误当全部物理精度已足够。目标N的最小奇异值明确非零，此处的相对诊断可用；676精确零权重背景仍禁止做比值。质量导数的数值检查在本可逆样本用trace公式，解析结论不依赖逆矩阵。

第二组在a=.23、L=7核有限软矩阵，P_u的幂等缺陷约.24098，明确不是精确投影。直接有理幂公式与谱计算误差4.29×10⁻¹⁵，规范下全质量合同误差4.66×10⁻¹⁵；原带相位权重／二来源输送相对误差分别1.75×10⁻¹⁴、3.19×10⁻¹⁴。未取绝对Pf替代原权重，未做正半定投影。

L的谱函数计算不等于实际构造数百万层局部体。完整平均的结论来自以上解析支配，不来自这两组有限样本。既有[障碍背景全谱入口结果](677/drafts/rational_limit_probe_results.json)及[首次平台诊断](677/drafts/rational_support_probe_results.json)原样保留；没有重复登记为新发现。

## 7. 下一项

接[678入口](678/drafts/STATUS.md)：以相同Wilson邻接算子构造全谱有理函数的稀疏块边界身份，检查其是否可提升为式(7)包含原质量、全部来源与体减除的完整Grassmann表示。先核非交换矩阵次序与相位，避免只验证一个传播子便宣称过程等价。

这属于同一候选的表示整合；若仅能得到数学线性系统，明确停在该层。原H_F、真实时间及共同尺度映射仍是不同验收项。空间382—386、425、522—523按§1直接复用，统一目标不改。

