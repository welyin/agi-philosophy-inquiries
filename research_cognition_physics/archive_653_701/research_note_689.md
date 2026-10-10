# 第689轮：原物理Weyl观测与CAR规范投影的同一插入

日期：2026-10-02。接[688](research_note_688.md)，直接继承[653](research_note_653.md)时间链及[670](research_note_670.md)全部物理Y来源正则字典；遵从[归属补正](687/drafts/attribution_correction_653.md)。[代码](689/joint_physical_car_projector.py)、[结果](689/joint_physical_car_projector_results.json)、[核验](689/research_round_689_checks.json)、[条件账](689/unified_physics_condition_ledger_689.md)。两组新核验、十六式；主代理审查，无独立代理审查。

## 1. 新增接口及旧结论复用

688给出扩展Gauss投影后CAR单态概率p_N，但没有把相应投影量写入原物理Y观测。本轮补上该字典：由原矩阵固定等时接触项，再构造同时间片的有限Grassmann多项式，其原插入恰等于重建CAR中的群变换；对原群Haar平均即得688的规范单态投影。

这不是由配分函数相同猜测观测相同。实际插入、全部独立时间链路、多个时间片的有序群词以及零Weyl权重处的未归一来源，都使用同一矩阵身份。

|项目|来源／本轮范围|
|---|---|
|认知动机|同一状态识别须保持实际可读观测，不能只匹配配分函数|
|额外建模输入|653纯时间、空间平凡、m₀=1、φ=0、原G／32CAR／S⁹／反周期闭合|
|旧解析结果|653全时间转移；654受控辅助脱耦；661原Weyl代数；670含零模的完整来源；679反射合同|
|本轮解析增量|原等时接触字典、局部群插入多项式、多个时刻同一CAR群词及原Gauss投影读数|
|验证|完整原32模式、非交换G链路、K或Cayley分母奇异、独立逐时规范变换|
|未完成|原空间动态Q₀正性、完整H_F身份、实际s仪器、共同连续及量子GR|

空间按[679逐项表](research_note_679.md)直接复用382—386、425、522—523：384额外Lipschitz不恢复，386与425是替代路线，523共同条件实现已存在；仍需与当前物质、h、s、Gauss的实际对象匹配。本轮没有新增或撤回空间假设，也不把下文的“局部接触项”当作光滑坐标证明。

## 2. 原Weyl字段自带确定的接触项

用653的反周期移位S与γ₄适配帧。令Q=S†⊗I₂作用于全部32N个负手征系数；J±仅为固定spin基，和661固定手征基相差常数基变换。原矩阵直接给

$$
K=J_+^\dagger Dv=\frac{I-Q}{2},\qquad
A=J_-^\dagger v=\frac{I+Q}{2},\qquad
w=Ac,\quad \bar w^{\mathsf T}=e,\quad S_W=e^{\mathsf T}Kc.
\tag{1}
$$

正指数约定继承661、670。归一有序收缩C=〈w barw〉在K可逆处为

$$
C=-AK^{-1}=I-2(I-Q)^{-1},\qquad
C_{\rm coh}=-(I-Q)^{-1},\qquad C=2C_{\rm coh}+I.
\tag{2}
$$

接触I严格为同时间、同模式项。固定背景下物理部门为Gaussian，所以不只二点函数，所有Wick来源系数都有同一关系。为避免Grassmann来源顺序歧义，定义特征多项式的系数直接为按指标升序排列的Wick矩：

$$
\mathcal G_C=\begin{pmatrix}0&C\\-C^{\mathsf T}&0\end{pmatrix},\quad
J=\begin{pmatrix}0&I\\-I&0\end{pmatrix},\quad
\mathcal W_C(j):=\exp\!\left(\tfrac12j^{\mathsf T}\mathcal G_Cj\right)
=e^{\frac12j^{\mathsf T}Jj}\mathcal W_{\rm coh}(\sqrt2j).
\tag{3}
$$

这只是选定来源顺序后的代数定义；转回670的源指数须保持其原顺序和符号。零权重处不使用归一逆核，下面直接给未归一多项式。

在时间片t，反周期条件使(Qᴺ)tt=−U_t，其中U_t为32维原闭合群词、含两份spin。于是

$$
C_{tt}=(U_t-I)(U_t+I)^{-1},\qquad
\langle a_i^\dagger a_j\rangle_{U_t}
=\left[U_t(I+U_t)^{-1}\right]_{ji}
=\tfrac12\delta_{ij}+\tfrac12(C_{tt})_{ji}.
\tag{4}
$$

这是扭曲迹的条件关联，单个群背景不必为正态。U_t=I时，原Ctt=0而CAR占据为1/2；把等时原w直接当作缩放后的正规CAR会漏掉这个确定的接触项。646重合极限边界和679接触抵消是继承结论，本轮新增是当前实际纯时间CAR插入所需的明确系数。

## 3. 原物理变量中的群插入

令H为作用于时间片t的32维单粒子矩阵，实际投影时取H=R(h)⊗I₂、h∈G。定义

$$
D_H=(I+H)/2,\qquad B_H=(H-I)/2,\qquad
\mathcal O_H(w_t,\bar w_t)
=\det D_H\ \exp\!\left(-\bar w_tB_HD_H^{-1}w_t\right).
\tag{5}
$$

式(5)先写在D_H可逆片。它不是新增动作或外部控制：只是原物理观测代数中对应Γ(H)的符号。Grassmann指数有限截断。其每阶系数为detD_H乘B_HD_H⁻¹的子式；Cauchy–Binet及Jacobi补子式把它改成D_H、B_H的多项式。因此**O_H在所有H上有唯一多项式延伸**，无需删去−1本征值背景。

令P_t选择该时间片的32个分量。将式(5)直接插入原eᵀKc的Berezin积分，而非另指定一个态，得到未归一物理因子

$$
I_W(\mathcal O_H)
=\det\begin{pmatrix}K&P_t^{\mathsf T}B_H\\P_tA&D_H\end{pmatrix}
=2^{-32N}\det(I+U_tH).
\tag{6}
$$

在D_H可逆处第一项是Schur行列式；两端均为多项式，故对K或D_H奇异仍成立，不除原detK。原辅助因子、S⁹积分及其相位不因只插物理Y而改动。670保证该原Weyl插入也是原完整辅助表示中的同一Y多项式。

另一种核对是在非零权重片用式(2)：

$$
\langle\mathcal O_H\rangle_{U_t}
=\det(D_H+B_HC_{tt})
=\frac{\det(I+U_tH)}{\det(I+U_t)}
=\frac{\operatorname{Tr}_{\mathcal F}\{\Gamma(U_t)\Gamma(H)\}}
{\operatorname{Tr}_{\mathcal F}\Gamma(U_t)}.
\tag{7}
$$

D_H、B_H彼此对易，和Ctt一般不对易；这里只用行列式的Sylvester身份，不能任意交换矩阵。式(6)而非式(7)负责零权重延伸。

## 4. 多个时间片共用同一有序CAR过程

不能仅靠单插入签收时间组合。令每片有独立H_t，组成块对角H；未插入片设H_t=I。D=(I+H)/2、B=(H−I)/2，则

$$
I_W\!\left(\prod_{t=0}^{N-1}\mathcal O_{H_t}\right)
=2^{-32N}\det\begin{pmatrix}I-Q&B\\I+Q&D\end{pmatrix}
=2^{-32N}\det(I-HQ).
\tag{8}
$$

证明在D可逆处左乘D给D(I−Q)−B(I+Q)=I−HQ；再作多项式延伸。没有要求不同时间群元可交换。反周期块移位的循环行列式进一步给

$$
\det(I-HQ)=\det(I+\mathcal U_t),\qquad
\mathcal U_t=-(HQ)^N_{tt},\qquad
\det(I+\mathcal U_t)=\operatorname{Tr}_{\mathcal F}\Gamma(\mathcal U_t).
\tag{9}
$$

这明确给出局部H_t与原链路交替相乘的实际顺序。辅助转移保持653的原规范链路，H_t只插在CAR上，不能误把它也旋转进辅助部门。

独立逐时原规范变换g_t使Q→gQg†、H→gHg†，式(8)不变；局部O_H相应协变。对h的Haar平均给规范不变偶物理多项式。

同一时间片的普通Grassmann乘法不自动等于算符乘法：

$$
\mathcal O_{H_1}\mathcal O_{H_2}\ \not\equiv\ \mathcal O_{H_1H_2}
\quad\text{（同一时间片的普通符号积）}.
\tag{10}
$$

必须先按算符顺序合并H或使用相应接触／符号乘法。不同时间片的式(8)已经具有正确时间序；没有把符号字典错误宣称为普通外代数的乘法同构。

## 5. 688的状态差异成为原Y插入的精确读数

定义局部原物理多项式与其匹配的CAR算符：

$$
\mathcal O_{\Pi,t}=\int_G\mathcal O_{R(h)\otimes I_2}(w_t,\bar w_t)\,dh,
\qquad \Pi_{\mathcal F}=\int_G\Gamma(R(h)\otimes I_2)\,dh.
\tag{11}
$$

O_H系数连续且次数有限，Haar积分逐系数存在。它无E插入，作用在原物理Y子代数；没有识别为真实标量s的仪器或新认知装置。

对每个原holonomy g，先做物理插入的h积分。Haar不变性给

$$
\int_G\chi_{\mathcal F}(gh)\,dh=n_0,\qquad
Z_N[\mathcal O_{\Pi,t}]=2^{-32N}n_0\sum_{\ell=0}^{8}s_\ell\mu_\ell^N.
\tag{12}
$$

再对原g施Gauss以及全S⁹积分，分母为688的2⁻³²ᴺΣm_ℓμ_ℓᴺ。因此

$$
\langle\mathcal O_{\Pi,t}\rangle_{\rm original}
=\frac{n_0\sum_\ell s_\ell\mu_\ell^N}{\sum_\ell m_\ell\mu_\ell^N}
=p_N<1\quad(N<\infty).
\tag{13}
$$

这是原未归一插入的多项式证明再除653已证严格正的完整条件分支配分函数；没有删去逐背景零点。实际值N=1为0.0960279096992471，N=2为0.305869658315581，N=4为0.741658077324403。

对两个不同时间片，Π_F与原群演化对易且Π_F²=Π_F。由式(8)精确得到

$$
\langle\mathcal O_{\Pi,t}\mathcal O_{\Pi,u}\rangle_{\rm original}=p_N
\quad(t\ne u),\qquad
\langle\Pi_{\mathcal F}\rangle_{\sigma_{\mathcal F}}=1.
\tag{14}
$$

故688的有限尺度差别不是仅停留在一份配分函数表示里的辅助记号。现在有原Y观测可读到它。零Hamiltonian、无辅助的CAR单独参考σ_F仍不是原完整H_F态；后者的玻色和规范部门也可以补偿CAR电荷。

## 6. 检验、范围和后续

原矩阵核验：N=1—4、全部16内部通道及两spin，最大接触身份误差5.67×10⁻¹⁴；0／2／4／6／8来源的实现与解析Wick身份相符。原群产生的K零模、D_H奇异及二者同时奇异均检验未归一块行列式，最大误差6.31×10⁻¹⁴。至少一个K奇异例中插入后矩阵可逆，说明零权重不能被直接删去。

非交换多时间片及独立逐时规范变换核验最大行列式误差4.89×10⁻¹⁴。Haar后的p_N按688冻结有理数据逐项完全相同；不以群抽样近似完成投影。数值承担实现检查，式(6)、(8)、(12)承担所有有限N的证明。入口代码初稿曾错误要求零模插入的未归一行列式绝对值大于1；改为检验插入矩阵的非零最小奇异值，这不改变解析身份。

本轮减少了一个独立字典要求：纯时间分支的等时修正、群插入、规范投影和原状态读数已由同一K、A确定，无需为每个背景另拟合。但必须保留

$$
\text{本纯时间分支的原Y／CAR投影身份}
\quad\not\Longrightarrow\quad
\text{一般空间动态 }Q_0\text{ 的反射正性或原完整 }H_F\text{ 身份}.
\tag{15}
$$

C01、C14、C19与C22的此具体观测连接完成；C03真实s记录、C07—C09空间对象接口以及量子几何仍原样保留。654旧辅助脱耦不新计结论，四分支继续分别记账。

下一项[690入口](690/drafts/STATUS.md)将这份明确的**原物理子代数**带回673—675完整动态泛函：先核空间非平凡时同一插入是否仍具有原CAR算符含义，再核完整Gauss物理型。不得把静态或纯时间结论直接推广，也不再重做已闭合的球谐时间链。

$$
\text{先识别同一物理观测}\ \longrightarrow\quad
\text{保留完整 }H_b\text{、S}^9\text{ 及Gauss平均}\ \longrightarrow\quad
\text{验收该原物理子代数上的正性与过程身份}.
\tag{16}
$$

