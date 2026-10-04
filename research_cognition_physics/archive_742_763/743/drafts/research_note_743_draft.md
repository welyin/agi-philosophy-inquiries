# 第743轮：原量子标量、严格Gauss准备与非高斯记录关联

日期：2026-10-04。接[742](../../research_note_742.md)，保留[入口报告](research_note_743_working.md)。[代码](../joint_native_quantum_record.py)、[结果](../joint_native_quantum_record_results.json)、[核验](../research_round_743_checks.json)、[条件账](../unified_physics_condition_ledger_743.md)。两组检查、十八式；主代理审查，无独立代理或图像检验。

## 1. 结果与去重

**在598／623原完整固定图Hamiltonian内，原量子标量已经能够使费米子从高斯态发展出非高斯四点关联。可以选正规、有限能源且严格满足Gauss约束的初态，无需增加新物种或新耦合。** 原Yukawa参数与同一标量方差同时固定四点缺陷、标量—配对关联和初始纠缠；不是分别拟合的三笔资源。

这解除742限制在原完整模型中的适用前提，未证明633理想占据记录已实现。入口的128维条件参考与下面的真实Gauss准备不同：前者沿730的数值校准，后者是原有限图内的另一个明确准备。没有把有限图CAR空真空当成730连续Hadamard参考。

上轮完成742并保存743入口，属于进展。本轮核导航、最新笔记、结果和进程，未发现活跃Python研究；重查577、598、623、712及717—720。**717式(13)—(16)已经给原Majorana到sin s读口的力；719式(18)已经给配对的一阶生成。直接复用，不把它们重计为新发现。** 新增是原完整Gauss量子准备中的四点方差缺陷与共同关联，补入口尚未验证的物理约束范围。

|层次|地位|
|---|---|
|认知动机|实际相互作用应同时承担关联、读取和来源|
|继承输入|原有限图、原群／物种／Y、曲目标、正给定几何、598完整自伴H、623共同域|
|准备输入|下述正常Gauss波包和偶费米态；不宣称由H自动制备|
|解析增量|原完整H下的非高斯四点缺陷、同一标量—配对关联及纠缠首项|
|数值范围|原条件128维资料；原64模式稀疏CAR初始导数和原正常波包求积|
|保留缺口|原理想仪器、末读自治、永久记录、连续手征／相互作用及动态量子引力|

## 2. 入口的条件性结论保留

复用623、735的质量坐标，不重新推导：

$$
x=\phi/\sqrt F,\qquad B_{\rm mass}=\sum_{v,a}x_v^aO_{v,a},\qquad
\phi=\sqrt{\frac{2}{1+|x|^2/6}}\,x.
\tag{1}
$$

若只取标量—费米必要部门，初态为σ_b⊗ρ_f、C为实际x协方差，扣除平均二次演化后，入口直接计算得到

$$
\rho_{f,I}(t)=\rho_f-\frac{t^2}{2\hbar^2}
\sum_{ab}C_{ab}[O_a,[O_b,\rho_f]]+O(t^3).
\tag{2}
$$

完整Hb给平均二次漂移，不能在总方程中省略；它不改变式(2)该阶非酉项。此式不是Markov极限，也不凭一个随机矩阵求积便认定原量子过程已模拟完成。入口只验证给定连接的条件准备；下面专门补严格Gauss例。

742原记录伙伴的协方差切向量Jₐ给

$$
\mathcal W(t)=-\frac{t^2}{\hbar^2}\sum_aC_{aa}\operatorname{pf}(J_a)+O(t^3),
\qquad-\operatorname{pf}(J_a)=\tfrac14\|J_a\|_{\mathrm F}^2\geq0.
\tag{3}
$$

该Pfaffian身份由纯四维块的旋转切空间得到；代码及独立复核误差约1.39×10⁻¹⁷。高斯／Wick工具参考[Bravyi，§IV，式(17)](https://arxiv.org/html/quant-ph/0404180)，不将其作为原创性主张。

原实际数字中，给定C=diag(10⁻⁴,10⁻⁴,10⁻⁴,10⁻⁴,.04)得到系数1.20010263405×10⁻⁵。该例主要来自第四Higgs方向，不能说只有singlet在起作用。五个原直接质量矩阵均与占据读口矩阵正交且不全对易；由正协方差的双交换子，当前短时噪声一般扰动占据数，因此不等于633的无选择Lüders通道。保留入口全部限制，不继续扫描方差。

## 3. 改用原完整物理模型的正常准备

实际演化使用同一个

$$
H_F=H_b+B_D+B_M+B_{\rm hop},\qquad
B_M=\sum_vx_v^5\left(Y_s a_v^\dagger b_v^\dagger+\bar Y_s b_va_v\right),
\quad a_v=c_{v,30},\quad b_v=c_{v,31}.
\tag{4}
$$

Hb包含全部原玻色动能、现场势、测地边、电／磁项；Dirac和原跨边跳跃没有删去。固定几何权是本阶段输入。原32模式的其它物种也没有删去。

取原CAR空真空Ω，只作为有限图试验态，并取归一化光滑紧支撑的Gauss不变ψ：

$$
\Psi_0=\psi\otimes\Omega,\qquad
\psi\in C_c^\infty(\mathcal K^V\times G^E)^G,
\quad\|psi\|=1,
\qquad P_G\Psi_0=\Psi_0.
\tag{5}
$$

原真空规范不变；每个sterile配对也规范平凡。ψ可用717的Higgs径向／singlet紧包、其它节点同类包和Haar常数链路组成。全H保持Gauss。共同光滑核保证任意有限次H幂域，相关有界CAR期望有足够时间导数；涉及x的矩阵元再用623的势及算符域控制。

所有Dirac和跳跃在Ω上湮灭，所以精确的第一步是

$$
H_F\Psi_0=(H_b\psi)\otimes\Omega
+\sum_vY_s(x_v^5\psi)\otimes|v\rangle,
\qquad |v\rangle=a_v^\dagger b_v^\dagger\Omega.
\tag{6}
$$

不同v的配对向量正交。式(6)不等于宣称整个演化只在真空和单对空间闭合；后续Dirac、传播、多配对和Hb都会出现。

## 4. 四点缺陷的完整量子推导

记μ_v=〈x_v⁵〉ψ、v_v=Varψ(x_v⁵)，n_a=a†a、n_b=b†b。由于这些占据算符在Ω上为零，Taylor二阶项只需式(6)，得

$$
\langle n_a\rangle_t=\langle n_b\rangle_t
=\frac{t^2}{\hbar^2}|Y_s|^2\langle(x_v^5)^2\rangle_\psi+O(t^3),
\qquad
\langle n_an_b\rangle_t
=\frac{t^2}{\hbar^2}|Y_s|^2\langle(x_v^5)^2\rangle_\psi+O(t^3).
\tag{7}
$$

719的一阶配对身份直接复用，而正常自旋交叉两点最早为二阶：

$$
\langle a_vb_v\rangle_t=\frac{it}{\hbar}Y_s\mu_v+O(t^2),
\qquad\langle a_v^\dagger b_v\rangle_t=O(t^2).
\tag{8}
$$

偶高斯态必须满足〈n_an_b〉=〈n_a〉〈n_b〉−|〈a†b〉|²+|〈ab〉|²。真实演化与此关系的差为

$$
\mathcal K_v(t):=\langle n_an_b\rangle_t
-\langle n_a\rangle_t\langle n_b\rangle_t
+|\langle a_v^\dagger b_v\rangle_t|^2-|\langle a_vb_v\rangle_t|^2
=\frac{|Y_s|^2v_v}{\hbar^2}t^2+O(t^3).
\tag{9}
$$

只要Y_s≠0且v_v>0，所有充分小的非零t都有严格非高斯约化费米态。任何非零正常配置L²包不可能全部支撑在一个x_v⁵常值的零测度超面，因此该准备类的v_v严格正。这里的存在范围是原固定图，不是连续场真空结论。

尤其μ_v=0时，一阶配对均值可以消失，但四点缺陷仍严格存在。仅看平均标量系数或平均配对，会漏掉这个过程。

## 5. 与内部关联、纠缠共用同一资源

在同一初态与H下，标量—配对连通关联的一阶项为

$$
\langle x_v^5 a_vb_v\rangle_t
-\langle x_v^5\rangle_t\langle a_vb_v\rangle_t
=\frac{it}{\hbar}Y_s v_v+O(t^2).
\tag{10}
$$

Hb在这一导数中乘上初始〈ab〉=0，不贡献该项。它不是读出另一份额外参考的关联。

更强地，将HΨ₀投影到玻色和费米两边均正交于初态的空间，得到

$$
\chi_\perp=\sum_vY_s[(x_v^5-\mu_v)\psi]\otimes|v\rangle,
\qquad\|\chi_\perp\|^2=|Y_s|^2\sum_vv_v.
\tag{11}
$$

对纯整体态作短时Schmidt展开，费米约化态的线性熵满足

$$
1-\operatorname{Tr}\rho_f(t)^2
=\frac{2t^2}{\hbar^2}\|\chi_\perp\|^2+O(t^3)
=\frac{2|Y_s|^2t^2}{\hbar^2}\sum_vv_v+O(t^3).
\tag{12}
$$

不同节点的玻色涨落可以相关；式(11)的配对向量正交仍消去此处交叉项。式(9)—(12)不是对原相互作用另加的强度，也不意味着任意目标记录或无限稳定记忆已经形成。

## 6. 与旧实际读口和来源连接，避免重做717

717已证明原e₊(s)=1/2+sin s/4只在相应二阶质量力中接到Majorana：

$$
\mathscr F_{e_+}=
\frac{\sqrt F\cos s}{4w_v}
\left(Y_s a_v^\dagger b_v^\dagger+\bar Y_s b_va_v\right).
\tag{13}
$$

Dirac力消去来自原K，不是删掉相互作用。为了说明共同模型已有读出方向，可直接在相同ψ上使用717同类的偶态

$$
|\eta_\pm\rangle=\frac{\Omega\pm e^{i\arg Y_s}|v\rangle}{\sqrt2},
\qquad
\rho_\pm=|\psi\eta_\pm\rangle\langle\psi\eta_\pm|.
\tag{14}
$$

它们严格Gauss，玻色初始边缘相同；取ψ在x⁵上对称，μ=0，原平均总能源也相同。由旧力身份和实际Taylor域得到

$$
p_{+|\eta_+}(t)-p_{+|\eta_-}(t)
=-\frac{|Y_s|}{4w_v}\langle\sqrt F\cos s\rangle_\psi\,t^2+O(t^3).
\tag{15}
$$

这项作为717的直接复用，不另计新测试。它读取配对相干，不能替换为已经读出了n_f；该末端效果仍是声明的操作输入，完整未知态仪器需下一轮核。选择支持|s|<π/2，使该系数严格非零；原R=0紧包已经满足。

所有这些态仍由同一H演化，故

$$
\langle H_F\rangle_t=\langle H_F\rangle_0,
\qquad G_\lambda=\partial_\lambda H_F,
\qquad\text{实际来源为同一态上的 }\langle G_\lambda\rangle_t.
\tag{16}
$$

全能源守恒复用自伴自治H，来源域复用623。不是说各子系统能源不变，也不表示该量子准备已解Einstein约束。741的条件性半经典连接保持，但不能将其730参考的数值原样挪给此不同准备。

## 7. 复算与物理范围

第一组复核冻结入口的128维原准备与质量矩阵，首次计入正式检查；没有把入口另算一轮。第二组在原两节点64模式CAR中保全部Dirac、Majorana及非平凡群链路跳跃，直接核式(6)—(8)，最大误差1.56×10⁻¹⁷。

正常波包完全复用717的R=0、半径.7与宽度1，原测度为d⁵x/√(1+|x|²/6)，不采用配置点态：

$$
\mu=0,\qquad v=0.11362789527824811,\qquad
\langle\sqrt F\cos s\rangle=1.2440968689246032.
\tag{17}
$$

原Y_s=.31+.09i，在ℏ=1时得到

$$
|Y_s|^2v=0.011840026687993452,\qquad
|Y_s|v=0.036679112756376524,\qquad
2|Y_s|^2v=0.023680053375986904.
\tag{18}
$$

这三数依次校准四点缺陷、标量—配对连通导数的模、每个同型节点的纠缠首项。48到80阶求积最大变化3.79×10⁻¹²，只是收敛诊断，不是区间认证。严格正性来自方差和支持证明，不靠浮点小数。没有数值模拟全Hb演化、全Gauss热参考、实际连续或引力。

## 8. 条件合并与下一项

C03／C15／C19／C22现在有一份原完整固定图上的共同资源见证：实际量子标量、Gauss、原物种、非高斯关联及继承的来源／读口可以共存。它比单纯声明“量子玻色会超出高斯类”更强，也明确何处换了准备。

接[744原相互作用诱导的实际仪器](../../744/drafts/STATUS.md)：以原sin s读口和同一H，核未知偶sterile编码输入的完整效果、真实后态及来源，而不是只给两个均值；保准备和末读输入，先不加新物种或优化寿命。复用577、717—718、623—625及743共同域。空间382—386、425、522—523、604及649／699保持。统一目标未改、未完成。
