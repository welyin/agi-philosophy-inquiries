# 第602轮：量子来源、参考态与实时响应的共同匹配

日期：2026-10-01。接[601](research_note_601.md)，回查[529](../archive_467_530/research_note_529.md)、[558](../archive_554_584/research_note_558.md)、[586方法桥](586/closed_time_path_bridge_review_586.md)、[590](research_note_590.md)—[591](research_note_591.md)、[598](research_note_598.md)—[600](research_note_600.md)。[代码](602/joint_quantum_response_matching.py)、[结果](602/joint_quantum_response_matching_results.json)、[核验](602/research_round_602_checks.json)、[条件账](602/unified_physics_condition_ledger_602.md)。三组检查、十式，主代理审查；无新增独立代理审查。

## 1. 整合问题与新结果

601给共同经典有效阶，但没有给同一量子态的实时反作用。局部Euclidean圈系数、静态自由能与实时间响应是不同对象，不能只凭“来自同一作用”直接互换。

本轮将598原全部32个费米模式的冻结质量部门接入**同一参考态、平均来源、因果响应、涨落**。计算表明：固定逆温的静态Hessian与隔离演化的零频响应一般不同，差额由同一来源的守恒涨落固定。原参数下差额矩阵严格非零，并有实际酉演化验证。所需补充不只是一个有效势数值，还包括参考态怎样随源变化以及保留哪份双时相关。

一般谱恒等式是成熟结果，新增的是原Dirac＋Majorana全质量部门的对象映射、共同来源矩阵、接触项和实际演化检查。没有引入新的认知装置。

|层次|地位|
|---|---|
|认知动机|共同物质的几何来源和标量力，应与其真实量子态及动态响应一致|
|模型输入|598全部一代32模式CAR质量矩阵、原Higgs/singlet真空值、冻结经典背景；给定β>0|
|准确部门|单节点质量作用，未纳入链路跳跃、动态玻色场和Gauss投影后的全相互作用热态|
|解析连接|静态—动态差额、保守来源涨落及统一温度／lapse参考|
|数值验证|64维BdG全质量块、独立16维中性Fock空间、精确酉淬火|
|未完成|完整Gauss热态与来源域、连续重整化、热化机制、全量子约束|

## 2. 成熟工具及与旧轮次的区别

[Watzenböck等，SciPost Phys. 12,184，§2](https://arxiv.org/html/2112.02903v2)区分孤立Kubo响应与等温静态响应，并将差额写成零频守恒相关。[Hu—Verdaguer，§4](https://arxiv.org/html/0802.0658v1)将平均应力、响应与噪声放在共同闭合时间路径中。本轮采用这些成熟方法，不导入其Hubbard动力学或把给定引力作用当作已生成。

529已区分虚时间响应与等时方差；558已区分能量Hessian与瞬时噪声；591在原玻色非简并真空带上建立亚隙响应。它们没有计算598新增质量部门的混合参考态、lapse与Yukawa来源矩阵。本轮不重做578—579细化障碍，也不重新定义CTP。

## 3. 同一量子态下的精确匹配条件

先在有限维物理空间或声明的有限谱规整上，令H(λ)为C² Hermitian族，G_a=∂_aH、C_ab=∂_a∂_bH。在λ=0取

$$
\rho_\beta=Z^{-1}e^{-\beta H},\qquad
\mathcal F_\beta=-\beta^{-1}\log Z,\qquad
\partial_a\mathcal F_\beta=\langle G_a\rangle,\qquad
\partial_a\partial_b\mathcal F_\beta=\langle C_{ab}\rangle-\chi^E_{ab}.
\tag{1}
$$

χᴱ为连接的虚时间积分协方差，不是等时方差。正号源扰动+λG下，隔离演化的响应写作δ〈G〉=(〈C〉−χᴿ)λ。用H的本征值E_n、G矩阵元及p_n=e^(−βE_n)/Z展开：非等能项在χᴱ及χᴿ(0)中均带(p_n−p_m)/(E_m−E_n)；相同能量块只在χᴱ中留下βp_n项。

用全部简并谱投影而非任取一个对角基底，定义

$$
\mathcal P_0(G_a)=\sum_E P_EG_aP_E,\qquad
D_{ab}=\operatorname{Re}\operatorname{Tr}\rho_\beta
[\mathcal P_0(G_a)-\langle G_a\rangle]
[\mathcal P_0(G_b)-\langle G_b\rangle],\qquad
\boxed{\chi^E_{ab}-\chi^R_{ab}(0)=\beta D_{ab}.}
\tag{2}
$$

D是半正定实对称矩阵。对实方向v，vᵀDv是守恒投影来源的方差，所以“同一态、同一来源、相同接触项”下，两响应相等当且仅当该方向D为零。有限温忠实态中，这要求投影来源在该支持上为常数。

$$
K^R(0)-K^E=\beta D,\qquad
K^R=\langle C\rangle-\chi^R,\quad K^E=\mathcal F_\beta'',\qquad
D_{ab}=\lim_{T\to\infty}\frac1T\int_0^T
\frac12\langle\{\delta G_a(t),\delta G_b(0)\}\rangle\,dt .
\tag{3}
$$

最后一式是有限谱的Cesàro时间平均；不能说有限系统的相关函数逐点收敛。来源加参数依赖的身份项c(λ)I会同步改变平均和接触项，但连接涨落D不变，不能用这种能量零点修复此差额。

该恒等式不否定有完整谱资料的Euclidean理论与实时理论之间的正确解析延拓；被排除的是只取静态零频值、丢掉零频特殊项或改动参考态后仍宣布因果过程一致。

## 4. 原质量部门的具体对象

保持598的Y参数和质量归一，取原真空h=√u_h、s=√u_s，φ=(0,h,0,0,s)。32模式二次CAR算符B(s)含全部Dirac块与右手中微子Majorana配对。冻结的玻色背景可规范协变地变换，但这一条件费米态不是完整Gauss物理态。

加入一份均匀lapse源n，并把s作为原标量源：

$$
H(n,s)=nB(s),\qquad
G_n=B,\quad G_s=nB',\qquad
C_{nn}=0,\quad C_{ns}=B',\quad C_{ss}=nB'' .
\tag{4}
$$

G_s正对应600不能省略的质量标量力；均匀n是当前条件部门的时间权，不是已经构造所有局部lapse及引力约束。这里先明确一份最小共同来源，不从H[1]自动推出全局部耦合。

在Nambu表示Ψ=(c,c†)中

$$
B=\frac12\Psi^\dagger\mathcal B\Psi+\frac12\operatorname{tr}h,\qquad
\mathcal B=\begin{pmatrix}h&\Delta\\\Delta^\dagger&-h^T\end{pmatrix},\qquad
\mathcal F_\beta=-\frac1{2\beta}\operatorname{Tr}
\log(1+e^{-\beta n\mathcal B}).
\tag{5}
$$

原h及其所有s导数的迹为零，因此最后一式无需另加常数；一般BdG模型必须保留tr h/2。64维是Nambu加倍，不是额外物种或64个独立Weyl场。

令u=F^(−1/2)，固定h时F=M−(h²+s²)/6；原分子矩阵随s仿射。实际质量导数使用

$$
u'=\frac{s}{6F^{3/2}},\qquad
u''=\frac1{6F^{3/2}}+\frac{s^2}{12F^{5/2}},\qquad
\mathcal B'=u'\mathcal B_{\rm num}+u\mathcal B_{\rm num}',\quad
\mathcal B''=u''\mathcal B_{\rm num}+2u'\mathcal B_{\rm num}'.
\tag{6}
$$

这保留了F的变化，不能只对Majorana分子的s求导而冻结Dirac块。各矩阵来自598原函数，没有另挑方便的二能级模型。

## 5. 全32模式的来源、噪声及响应

令ε_i为n𝓑的本征值，f_i=(1+e^(βε_i))^(−1)，g_a为相应Nambu来源。Wick代数给

$$
\langle G_a\rangle=\frac12\sum_i f_i(g_a)_{ii},\qquad
\chi^R_{ab}(0)=\frac12\sum_{\epsilon_i\ne\epsilon_j}
\frac{f_i-f_j}{\epsilon_j-\epsilon_i}
\operatorname{Re}[(g_a)_{ij}(g_b)_{ji}],\qquad
D_{ab}=\frac12\sum_{\epsilon_i=\epsilon_j}
f_i(1-f_i)\operatorname{Re}[(g_a)_{ij}(g_b)_{ji}].
\tag{7}
$$

包括全部等能块，避免因spin／color简并任取基底而漏项。代码区分等能的容差为10⁻¹⁰，物理身份来自解析简并块公式，不由浮点阈值定义新物理原则。

同一来源的对称相关以及非零频率的平衡关系为

$$
N_{ss}(t)=\frac12\sum_{ij}f_i(1-f_j)|(g_s)_{ij}|^2
\cos[(\epsilon_i-\epsilon_j)t],\qquad
f_i(1-f_j)+f_j(1-f_i)
=(f_i-f_j)\coth\frac{\beta(\epsilon_j-\epsilon_i)}2 .
\tag{8}
$$

取ℏ=1；第二式只写在非零能差上。零频D另保留，不能把ω=0直接塞入coth后把D丢掉。响应与噪声受到同一态和谱的约束，不是两份任意可调输入；本轮未假设Markov浴或白噪声。

## 6. lapse变化同时牵涉参考温度

H=nB时，B与演化对易，所以χᴿ_nn=0；有限温态的βVar(B)却通常非零。固定坐标β计算静态导数时，ρ∝exp(−βnB)的能级占据随n改变，与固定初态的隔离演化不同。

$$
\beta_{\rm proper}=n\beta_{\rm coordinate},\qquad
n\mapsto n',\quad
\beta_{\rm coordinate}\mapsto\beta_{\rm proper}/n'
\quad\Longrightarrow\quad
\rho\ \hbox{不变}.
\tag{9}
$$

式(9)只处理静态均匀lapse和既定时间单位；不是任意弯曲时空的Tolman定律推导。它说明一部分所谓“不相容”来自比较了不同参考态。s改变时，单一温度重标定一般不能保全部谱占据；需要另给实际状态转换或热化条件，不能免费假设。

Gibbs态在此是来源选择，不表示宇宙外有热浴；从一个共同内部系统实际得到、更新该态仍未证明。等温、固定占据与实际开放演化作为不同分支保留。

## 7. 复算与真实演化

原h≈.665443、s≈.545864、F≈1.876536，β=2时：

$$
D=\begin{pmatrix}
1.0000879248&.0622210600\\
.0622210600&.0178992185
\end{pmatrix},\qquad
\chi^R_{ss}(0)\simeq.0247561108,\quad
\chi^E_{ss}\simeq.0605545478 .
\tag{10}
$$

D的本征值约.0139732、1.0040139，两个来源方向均存在非零静态／动态差额。数值负的某些自由能曲率不能直接判引力不稳；这里只求条件物质来源，未含几何裸作用和全部相互作用。

三组核验：

1. 全32模式质量导数和自由能来源：静态Hessian与直接来源有限差分最大差2.95×10⁻¹¹；一阶质量导数差5.40×10⁻¹²，二阶差5.38×10⁻⁷。原函数与解析导数分别计算。
2. 中性四模式在该背景下为封闭二次因子。独立构造16维Fock矩阵，不用BdG谱直接生成目标值；平均、孤立响应及四个时刻噪声的最大差6.94×10⁻¹⁷。保持初始ρ不变，直接用B(s±10⁻⁵)作酉演化，再读B′(s±10⁻⁵)，t=.5、2、7的线性变化与因果核最大差3.46×10⁻¹²。读数包含接触项，未用重新热化的末态替代真实演化。
3. β=.7、2、5时逐跃迁检验式(8)，误差≤8.89×10⁻¹⁶；n从1变1.23，按式(9)调整β后B均值差≤8.89×10⁻¹⁶。固定坐标β则B均值分别改变−.204048、−.432093、−.349285，确实换了态。

有限矩阵精确成立不等于全图Gauss热态已经存在，也不证明599零温UV系数不正确。零温非简并有隙且参考态一致时D趋零，恢复591型分支；有简并、零隙或非平衡态时须另核极限。

## 8. 合并后的条件与下一步

平均来源、静态有效系数、因果响应和噪声现在有一个共同匹配合同：必须使用相同变量、同一参考态族、同一接触项及完整所需频谱。只匹配一个静态势或共同质量矩阵不够。

该合同减少了将噪声和响应分别任选的自由；没有选出唯一宇宙态、β或热化率。局部有限EFT项只描述某个频率窗口，不能自动承担全部多时过程。

[603候选](603/drafts/STATUS.md)检查原完整非线性有限图、Gauss与费米扩展是否具有所需热态及可控响应，避免把本轮冻结质量因子当成全模型。先解决状态存在和来源域，不设计热化装置；有限截止与量子连续问题继续分开。
