# 第654轮：原完整质量、共同量子态与手征辅助测度的同一时间极限

日期：2026-10-02。接[653](../../research_note_653.md)、[原质量字典614](../../../archive_585_628/research_note_614.md)及[起步探查](spatial_mass_common_entry.md)。[代码](../joint_mass_state_time_limit.py)、[结果](../joint_mass_state_time_limit_results.json)、[核验](../research_round_654_checks.json)、[条件账](../unified_physics_condition_ledger_654.md)。

## 1. 问题与结论层级

按用户当前顺序，先整合条件，认知系统设计后置。本轮连接C04演化、C15质量、C16测度、C19状态、C20时间极限、C22来源。问题是：653的正时间表示能否接回598的原质量，而不是给辅助结构再设计一套装置。

**解析结果：**在明确声明的、空间平凡的物理质量插入下，全部原32个CAR模式的Pfaffian等于同一个正转移的Fock迹。有限步能谱误差从二阶开始，状态却有一阶变化；同一受控极限可以同时恢复原质量演化、热态、相关函数和标量来源。若653原辅助核保持不变，这个极限使其所有非平凡模式退耦，只留下规范单态。

这不是原完整四维手征理论的构造。空间传播、量子标量全配置积分、几何来源与真空项仍未接通；下文严格限定范围。原群、物种、复Yukawa系数和时间离散处方都是已有或声明的建模输入，不称其为认知原则推出的新定理。

## 2. 历史与新增输入审计

612已证明非零空间动量的overlap割线与有限Hamiltonian精确重建的边界，646已处理自由正时间窗口；本轮不重复它们。643已给完整CAR历史的Feynman–Kac表示和排序来源，本轮不把一般路径积分公式重新算成发现。614的固定粒子—空穴字典保证质量、群和来源相容，但未选择完整Euclidean Yukawa作用。

[Kikukawa §5.3](https://arxiv.org/html/1710.11618v3#S5.SS3)指出可加入Higgs和Yukawa耦合并适当使用手征投影，并未在该节完整指定这里要用的物理质量项；其一般局域性问题也仍单列。因此本轮新增输入是：在653的物理负手征时间坐标中加入下述局部配对项，正手征辅助积分沿用653。该分离是待连接的候选处方，不冒称文献已证明的完整标准模型构造。

## 3. 保留原全质量与有限CAR空间

用614与标量无关的固定Bogoliubov变换W，原单点质量变为纯配对。记f=32，F(φ)=M−|φ|²/6>0，保留原五个复Y系数，不重新拟合：

$$
H_m(\phi)=C(\phi)+C(\phi)^\dagger,\quad
C=\frac12\sum_{i,j=1}^{32}\Delta_{ij}(\phi)c_i^\dagger c_j^\dagger,
\quad\Delta^T=-\Delta,\quad\Delta=M_{16}(\phi)\otimes\varepsilon.
\tag{1}
$$

此M₁₆是614的对称全左手质量矩阵；不与标量势中的常数M混同。Hilbert空间仍是Λ(C³²)，没有增加物理粒子。常标量仅用于检验原条件背景因子，不把动态标量永远冻结成外部质量参数。

## 4. 实际候选Pfaffian与正转移的精确身份

N个反周期时间格，a>0。S的超对角元为1，末端回绕为−1。允许各格不同的原φ_t。以每格f个模式分块，定义：

$$
d=\tfrac12(I-S^\dagger)\otimes I_f,\quad
P=\operatorname{diag}_t(\tfrac a2\Delta_t),\qquad
\mathcal A=\begin{pmatrix}P&-d^T\\d&-P^*\end{pmatrix}.
\tag{2}
$$

这个处方的Berezin次序由质量为零时的正迹固定。令σ=(-1)^{fN(fN+1)/2}，则有任意N、任意反对称Δ_t的身份：

$$
T_a(\phi)=e^{-aC(\phi)}e^{-aC(\phi)^\dagger}>0,\qquad
\sigma\operatorname{Pf}\mathcal A
=2^{-fN}\operatorname{Tr}_{\mathcal F}
\big[T_a(\phi_N)\cdots T_a(\phi_1)\big].
\tag{3}
$$

**证明。**先将式(2)整体乘2。其齐次方程给ψ_t=ψ_{t−1}−aΔ_t†barψ_t、barψ_{t+1}=barψ_t−aΔ_tψ_t。因此每格递推矩阵恰为：

$$
\binom{\bar\psi_{t+1}}{\psi_t}
=U_a(\Delta_t)\binom{\bar\psi_t}{\psi_{t-1}},\qquad
U_a=\begin{pmatrix}I+a^2\Delta\Delta^\dagger&-a\Delta\\-a\Delta^\dagger&I\end{pmatrix}>0.
\tag{4}
$$

逐格消去变量只用单位对角块，结合反周期闭合给det(2A)=det(I+U_N⋯U₁)。右边也是式(3)Fock迹的平方：C与C†的Nambu生成矩阵分别是严格上、下三角块，指数相乘正是U；没有正常序迹常数遗漏。这里复用[配对费米子的迹身份](https://arxiv.org/html/1403.7824)，但保留实际Fock提升，不任取行列式主平方根。

最后，两边都是各Δ_t及Δ_t*元素的有限多项式；平方相同且在全部Δ_t=0时同为2^f，故多项式整体相同。穿过迹为零的位置不造成任意改号。恢复A的1/2比例得到式(3)。这同时说明身份不限常质量或质量矩阵彼此对易。

T_a正，并不意味着不同背景T的乘积有逐路径正迹。三格原复质量例子的归一迹有虚部.0022853381；这是原条件历史的相位，不是负概率，也不宣称完整反射正性已经建立。

## 5. 有限步状态为何不等于“只改质量数值”

常φ时取正算符的唯一自伴log，H_a=−log(T_a)/a。BCH在有限CAR空间给：

$$
H_a=H_m+\frac a2[C^\dagger,C]+O(a^2),\qquad
[C^\dagger,C]=\tfrac12\operatorname{Tr}(\Delta\Delta^\dagger)I
-c^\dagger\Delta\Delta^\dagger c.
\tag{5}
$$

真空常数与数算符项必须一起保留。对应Nambu矩阵B_a=−log(U_a)/a满足：

$$
B_a=
\begin{pmatrix}-\tfrac a2\Delta\Delta^\dagger&\Delta\\
\Delta^\dagger&\tfrac a2\Delta^\dagger\Delta\end{pmatrix}+O(a^2),
\qquad B_0=\begin{pmatrix}0&\Delta\\\Delta^\dagger&0\end{pmatrix}.
\tag{6}
$$

对Δ作奇异值分解，U分解成正2×2块。若σ_j是Δ的奇异值，其正准粒子能为：

$$
E_{a,j}=\frac2a\operatorname{arsinh}\frac{a\sigma_j}{2}
=\sigma_j-\frac{a^2\sigma_j^3}{24}+O(a^4).
\tag{7}
$$

能谱二阶接近不能取消式(5)的一阶本征向量变化。在β=Na、Ξ=(c,c†)列向量约定下，同一转移规定热态和等时相关矩阵：

$$
\rho_a=\frac{T_a^N}{\operatorname{Tr}T_a^N},\qquad
\mathcal G_a:=\langle\Xi\Xi^\dagger\rangle_{\rho_a}
=(I+e^{-\beta B_a})^{-1}.
\tag{8}
$$

中性原Majorana模式s=1.2的例子中，即使把原H的能量完全校准为E_a，a=1时所得态与实际ρ_a的迹距离仍为.0204960581，a=.25时仍为.0051783042。标量依赖的另一个Bogoliubov旋转可以改变比较字典，但其导数会进入来源；614的固定字典没有许可把这些导数丢掉。

## 6. 共同来源及同一受控极限

所有纯产生二次式相互对易，故对任一原实标量分量有精确插入：

$$
\partial_jT_a=-a\big[(\partial_jC)T_a+T_a(\partial_jC)^\dagger\big],
\qquad
\partial_j\log\operatorname{Tr}T_a^N
=N\frac{\operatorname{Tr}(T_a^{N-1}\partial_jT_a)}{\operatorname{Tr}T_a^N}.
\tag{9}
$$

原质量分子对φ线性，分母是√F。代码独立使用其解析导数，而非另拟一种来源：

$$
\partial_j\Delta(\phi)=
\Delta(e_j)\sqrt{\frac{F(e_j)}{F(\phi)}}
+\frac{\phi_j}{6F(\phi)}\Delta(\phi).
\tag{10}
$$

**受控共同极限定理。**令φ及所比较的C²参数族处在F>0内部的任意固定紧集，β在正紧区间。有限维矩阵指数和正log在a=0附近解析，且H_a在a=0可延拓为H_m。因此对算符及其一阶来源，统一有：

$$
\|H_a-H_m\|+\|\partial_jH_a-\partial_jH_m\|\le C a,
\quad T_a^{\beta/a}=e^{-\beta H_a}\longrightarrow e^{-\beta H_m},
\quad\rho_a\longrightarrow\rho_m\quad\text{以迹范数}.
\tag{11}
$$

式(11)在N=β/a为整数的格点族上成立；中间自伴log也定义其他a的函数。Duhamel身份和统一有界性同时给归一态的一阶参数导数收敛。任意固定有限个有界观测量及有序时间插入共享此极限，例如：

$$
\frac{\operatorname{Tr}\big(e^{-(\beta-t)H_a}A e^{-tH_a}B\big)}
{\operatorname{Tr}e^{-\beta H_a}}
\longrightarrow
\frac{\operatorname{Tr}\big(e^{-(\beta-t)H_m}A e^{-tH_m}B\big)}
{\operatorname{Tr}e^{-\beta H_m}},\quad 0\le t\le\beta.
\tag{12}
$$

因此不是分别配一个态、一个响应和一个时间规律。平滑有界背景历史也由每步T=I−aH_m+O(a²)及望远镜估计收敛到原排序指数；这不使任意条件历史迹自动正。

紧集假设很重要：F→0时原质量无界，本轮不以有限维CAR为由忽略配置域。623的原完整H共同域和643的热积分是已有工具，但没有自动证明本候选质量分割在全配置积分下可交换极限。

## 7. 现在把653的辅助测度放进同一个极限

记653实际核的转移为T_aux。其支撑为S⁹的ℓ=0,…,8球谐空间，不重做谱推导。原核不依赖a；最大特征值及比值为：

$$
\lambda_0=\frac{7429}{327680},\quad
r_\ell:=\frac{\lambda_\ell}{\lambda_0}
=\prod_{j=0}^{\ell-1}\frac{8-j}{17+j},\quad
r_1=\frac8{17},\quad \sum_{\ell=0}^{8}d_\ell=35750.
\tag{13}
$$

ℓ=0是唯一常函数，是原整个规范群G的单态。记P₀为其正交投影，归一核Q=T_aux/λ₀。对任何N：

$$
\|Q^N-P_0\|_1=\sum_{\ell=1}^{8}d_\ell r_\ell^N
\le35749(8/17)^N,\qquad
\frac12\left\|\frac{Q^N}{\operatorname{Tr}Q^N}-P_0\right\|_1
=\frac{\sum_{\ell>0}d_\ell r_\ell^N}{1+\sum_{\ell>0}d_\ell r_\ell^N}.
\tag{14}
$$

这是真正统一上界，不是有限时间长度外推。辅助相对真空的首个能隙则是：

$$
\Delta E_{\rm aux}(a)=\frac1a\log\frac{17}{8}longrightarrow+\infty.
\tag{15}
$$

因此“保留原固定核＋a趋零＋全部辅助模式保有有限非零能隙”三项不相容。允许辅助非单态退耦，就有下面的共同极限；若要保留非平凡辅助动力学，必须另给随a变化的处方及其测度来源。不是全部手征构造的反证。

## 8. 共同时间过程、真空项与适用边界

在本轮声明的物理／辅助分離候选中，常φ、单位时间holonomy时原始条件因子为：

$$
Z^{\rm raw}_{a,N}=2^{-32N}\operatorname{Tr}(T_a(\phi)^N)
\operatorname{Tr}(T_{\rm aux}^N).
\tag{16}
$$

移除明确的共同真空标量后，正转移的联合极限是：

$$
\mathbb T_a=T_a(\phi)\otimes Q,\qquad
\mathbb T_a^{\beta/a}\longrightarrow e^{-\beta H_m(\phi)}\otimes P_0
\quad\text{以迹范数；}\quad
\rho_a^{\rm joint}\longrightarrow\rho_m(\phi)\otimes P_0.
\tag{17}
$$

有限物理Fock空间、式(11)和式(14)直接给证明。在全部原扩展空间上，极限a→0后t>0的辅助过程是P₀，t=0仍是I，**并非全扩展空间上的强连续半群**；在存活子空间Fock⊗RanP₀才得到原质量的连续半群。不能把此退耦误称有限步全Hilbert同构。

原群协变仍共同成立：Δ(gφ)=R(g)Δ(φ)R(g)^T，故T_a(gφ)=Γ(R(g))T_a(φ)Γ(R(g))†，P₀对原群不变。对任何有定义的共同规范投影P_G，未归一极限可以受压缩：

$$
\|P_G(X_a-X)P_G\|_1\le\|X_a-X\|_1.
\tag{18}
$$

但非零带荷φ固定时T_a(φ)通常不与整个G对易；式(18)不把固定Higgs背景变成完整Gauss动力学。必须让标量同时变换或量子化才谈原全系统的Gauss约束。投影后的归一态还要求极限归一因子非零。

辅助来源也不能独立选择。对原holonomy的任一有限维生成元J，原表示U(g)满足：

$$
\left|\operatorname{Tr}(Q^N\mathscr U(g))-1\right|
\le\sum_{\ell>0}d_\ell r_\ell^N,\qquad
\left|\partial_\theta\operatorname{Tr}(Q^N e^{i\theta J})\right|
\le\|J\|\sum_{\ell>0}d_\ell r_\ell^N.
\tag{19}
$$

因为J在P₀上为零。原U(1)_Q在该辅助空间最大|Q|≤24，代码用实际球谐字符检查来源一同消失。它不会消去物理CAR自身的来源；本轮常标量响应保留式(9)。

原始整体标量并未消失。其每格对应的真空能为：

$$
E_{\rm vac}(a)=\frac{32\log2-\log\lambda_0}{a},\qquad
Z^{\rm raw}_{a,N}=e^{-\beta E_{\rm vac}(a)}
\operatorname{Tr}\mathbb T_a^N.
\tag{20}
$$

对固定a的归一物质态它可抵消，但对几何／lapse来源不能未经记账丢弃。因此本轮没有解决宇宙学常数、物理真空能或GR应力。这是下一轮应共同审查的接口，不能用“归一化”代替证明。

## 9. 可复算检查

四组全部通过，复用原Python与NumPy，无图像检查。

1. **全原质量身份：**保留32个模式、全部复Y、变化的五标量背景；1、2、3格实际Pfaffian与独立Fock分块迹的相对差≤1.63×10⁻¹⁵。原色块完全重建全矩阵，未删物种；原群协变误差≤1.96×10⁻¹⁵。
2. **同一态与五来源：**a从1减至.125，全32模式相关矩阵误差.0706692降至.00913256，呈一阶；能量误差.00739187降至.000119583，呈二阶。原全部五个标量来源的转移插入与独立谱导数差≤1.01×10⁻¹⁰；原中性Dirac／Majorana四模式显式密度矩阵与Nambu相关矩阵差≤5.18×10⁻¹⁶。
3. **有限步状态反例：**原中性Majorana配对的Pfaffian与精确能量相符，但只校准能量仍有上述非零迹距离，证明不能据相同配分谱签收状态。
4. **联合辅助极限：**用653冻结有理谱及原群字符，N=2、4、8、16、32、64检验相同上界与来源；N=32辅助态距单态约3.35×10⁻¹⁰。全N结论由式(13)—(15)证明，数值不代替证明。

## 10. 这轮实际合并了什么

有限步质量的状态、演化、相关函数和来源现在由同一候选转移决定；原辅助核在相同时间极限中的去向也被固定。不能再自由要求“质量变连续而辅助谱完全不变”。这是一条共同条件的相容性结果和一条限定要求的不相容结果。

尚未整合的重点是：原空间传播与该质量处方的实际兼容、原标量／规范量子过程下的共同极限、绝对真空项与同一几何来源。继续[655](../../655/drafts/STATUS.md)检查这些跨分支条件，优先核真空归一与已有几何来源的冲突，再明确空间极限的最低输入。保持完整统一目标，不把本时间分支称为标准模型或引力的推导，不转入认知设备设计。
