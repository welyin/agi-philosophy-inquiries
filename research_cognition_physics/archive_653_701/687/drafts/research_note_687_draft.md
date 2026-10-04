# 第687轮：完整非阿贝尔holonomy、球面正核与精确群积分

日期：2026-10-02。接[686](../../research_note_686.md)、[687入口](STATUS.md)。[代码](../joint_full_holonomy_positive_measure.py)、[结果](../joint_full_holonomy_positive_measure_results.json)、[核验](../research_round_687_checks.json)、[全条件账](../unified_physics_condition_ledger_687.md)。两组检查、十八式；主代理审查，无新增独立代理审查。

## 1. 实际推进及范围

上一目标轮完成685—686并更新正式证据，属于进展。本轮读取最新导航、686报告／结果／结项检查，未发现运行中的Python。回查615、657、659、672—675、679—680：固定E、固定边界负数及675少量群抽样都不能判断完整物理平均；不重复这些实验。

本轮利用615**每方向一格、空间链路单位、反周期时间、φ=0**的原有限诊断，将其阿贝尔holonomy结果推广到**整个原非阿贝尔商群**。新结果包括：

- 无需选群平方根的全局谱帧，保原Weyl行列式、辅助配对和相位。
- 全S⁹积分是群上的正型函数；在原群上有显式统一正下界。
- 完整原物理权重同样是群正型函数；整个群Haar积分可精确计算，非零性有解析证书。
- 群正核与物理时间反射正性严格区分。没有把此平坦一格部门替代673的多空间点、两时间片、完整H_b对象。

认知动机是群关系、辅助积分和物理权重要能共用正内积。群、表示、格点、边界、m₀=1及球面测度仍为继承输入。没有添加新认知原则或预先取绝对Pfaffian。

成熟接口：[Kikukawa，1710.11618v3，§4.1、§4.5—4.6](https://arxiv.org/html/1710.11618v3)区分特殊背景的非零／正性论证与一般复杂权重。615已指出原群不是其固定向量Spin(9)分支。本轮不援引一般正性结论，下面针对实际一格对象独立推导。只查文字，不做图像检验。

旧空间382—386、425、522—523按[679复用表](../../research_note_679.md)继承；本轮处理的是内部规范群，不把其参数当空间位置、真实方向壳或物理时间。

## 2. 原整个群与全局谱帧

复用614的精确字典：

$$
G=\frac{SU(3)\times SU(2)\times U(1)}{\mathbb Z_6}
\simeq S(U(3)\times U(2)),\quad
V(g)=\operatorname{diag}(z^{-2}C,z^3W),\quad
R(g)=\Lambda^{\mathrm{even}}V(g)\in SU(16).
\tag{1}
$$

这是原群的另一表示，不扩成新的Spin(10)规范理论。沿615内部优先的矩阵顺序，令C_R=(R+R†)/2、S_R=(R−R†)/(2i)：

$$
X=C_R\otimes I_4-iS_R\otimes\gamma_4,\qquad
H=(I_{16}\otimes\gamma_5)X,\quad H^2=I_{64},\qquad
D=(I+X)/2 .
\tag{2}
$$

式(2)来自原一格Wilson核及反周期边界，并非指定一个替代投影。避免对R取平方根，定义

$$
a=(I-R^\dagger)/2,\quad b=(I+R^\dagger)/2,\qquad
u=i(a\otimes V_++b\otimes\gamma_4V_+),\quad
v=i(a\otimes V_-+b\otimes\gamma_4V_-).
\tag{3}
$$

V±是615固定γ₅手征帧。因为a†a+b†b=I、a†b+b†a=0，且R与其函数可交换，直接相乘给

$$
u^\dagger u=v^\dagger v=I,\quad u^\dagger v=0,\quad
uu^\dagger=(I-H)/2,\quad vv^\dagger=(I+H)/2,\quad
\det[u,v]=\text{常数}.
\tag{4}
$$

最后一项：完整旋转i(a⊗I+b⊗γ₄)在γ₄=+1块是iI，在γ₄=−1块是−iR†，其行列式为(det R†)²=1。故没有随g变化的Jacobian或隐藏相位。R含−1本征值也不使该帧失效；物理Weyl权重可在此处为零。

## 3. 实际辅助配对及完整S⁹平均

复用615的T(E)=ΣE_aT_a、B、ε。定义原实10维表示O(g)：

$$
R(g)^*T(E)R(g)^\dagger=T(O(g)E),\qquad
u^T[T(E)\otimes B]u
=-T\!\left(\frac{E+O(g)E}{2}\right)\otimes\varepsilon .
\tag{5}
$$

原Clifford字典保证O(g)实正交且O(gh)=O(g)O(h)。式(5)由(3)展开：混合手征配对为零，剩下aᵀTa+bᵀTb=(T+R*TR†)/2。615已证固定定相后Pf[−T(e)⊗ε]=|e|¹⁶，故

$$
\operatorname{Pf}A_g(E)
=\left[\frac{1+E^TO(g)E}{2}\right]^8,\qquad
M(g)=\int_{S^9}\left[\frac{1+E^TO(g)E}{2}\right]^8d\nu(E).
\tag{6}
$$

每个原E的配对权重都在[0,1]，这里是实际Pf的恒等式，不是用其绝对值构造另一个模型。一般多格、非平坦背景不具有式(5)的简化。

原O是V的实化表示，或与其共轭等价的实表示。设V的五个单位圆本征值为z_i，∏z_i=1，A=(2I+O+Oᵀ)/4。由于五为奇数，

$$
2=\left|1-\prod_{i=1}^5(-z_i)\right|
\le\sum_i|1+z_i|
\le\sqrt{5\sum_i|1+z_i|^2},\qquad
\operatorname{Tr}A=\frac12\sum_i|1+z_i|^2\ge\frac25.
\tag{7}
$$

球面二阶矩与Jensen不等式因此给整个原群上一致的明确界：

$$
1\ge M(g)=\mathbb E(E^TAE)^8
\ge\left(\frac{\operatorname{Tr}A}{10}\right)^8
\ge25^{-8}>0 .
\tag{8}
$$

这不要求O有不变向量，也不把该界推广到任意Spin(10)元素或一般格点规范场；后者可以有O=−I。物理Weyl零点仍然保留。

若z_i=e^{iα_i}，A的五个成对本征值为λ_i=cos²(α_i/2)。球面在五个复平面的平方长度服从Dirichlet(1,1,1,1,1)。因此完整积分还可写成

$$
M(g)=\frac{h_8(\lambda_1,\ldots,\lambda_5)}{\binom{12}{8}}
=\frac{h_8(\lambda)}{495},\qquad
\sum_{k\ge0}h_k(\lambda)t^k=\prod_{i=1}^5(1-\lambda_i t)^{-1}.
\tag{9}
$$

对原超荷对角族，三个λ等于cos²θ、两个等于cos²(3θ/2)，(9)准确回到615的Beta(3,2)多项式。本轮不重新领取该阿贝尔积分。

## 4. 全物理权重的群正型证书

从同一v和原bar手征帧计算Weyl行列式，固定参考相位后，

$$
W(g)=\det(\bar vDv)
=\left|\det\frac{I+R(g)}2\right|^2,\qquad
Z(g)=W(g)M(g),\qquad Z(e)=1 .
\tag{10}
$$

det R=1保证实际行列式身份；不是事后删相位。虽然M处处正，W可为零，如原超荷θ=π/6。

正数函数不一定是正型函数，需要以下独立证明。令D₀(g)=1⊕O(g)，v_E=(1,E)/√2。则

$$
M(g)=\operatorname{Tr}[\rho_M D_0(g)^{\otimes8}],\qquad
\rho_M=\int|v_E^{\otimes8}\rangle\langle v_E^{\otimes8}|\,d\nu(E)
\ge0,\quad\operatorname{Tr}\rho_M=1 .
\tag{11}
$$

令Γ(R)=⊕_{k=0}^{16}ΛᵏR为有限外代数上的酉表示，dim ℱ=2¹⁶。标准外代数迹恒等式det(I+R)=Tr Γ(R)给

$$
W(g)=2^{-32}\operatorname{Tr}
[\Gamma(R(g))\otimes\overline{\Gamma(R(g))}].
\tag{12}
$$

所以Z是某个明确有限酉表示U在正密度ρ=2⁻³²I⊗ρ_M中的矩阵系数。对任意有限g_i和复c_i，

$$
\sum_{ij}\bar c_i c_j Z(g_i^{-1}g_j)
=\operatorname{Tr}\!\left[
\rho\left(\sum_i c_iU(g_i)\right)^\dagger
\left(\sum_j c_jU(g_j)\right)\right]\ge0 .
\tag{13}
$$

这是解析Gram证书，不依赖数值特征值是否恰好为正。特征Hilbert空间只是该实际群函数的表示证书，没有把它自动认作原物理CAR或H_F时间空间。

归一Haar平均P_G=∫U(g)dg是正交不变投影。外代数真空⊗真空及D₀的纯平凡分量组成不变单位向量；ρ在它上的权重是2⁻³²×2⁻⁸。由此，

$$
\int_G Z(g)\,dg=\operatorname{Tr}(\rho P_G)\ge2^{-40}>0 .
\tag{14}
$$

没有用有限群样本代替完整Gauss群平均。这仍是单holonomy群积分；673的独立两边界和完整H_b热密度不是式(14)。

## 5. 原整个群积分的精确值

Z为共轭不变函数，可在S(U3×U2)的秩四极大环面取z₁,…,z₄独立、z₅=(z₁z₂z₃z₄)⁻¹。包含颜色和弱部门的完整Weyl密度后，

$$
\int_G Z\,dg=\frac1{12}\operatorname{CT}
\left[
\prod_{1\le i<j\le3}(2-z_i/z_j-z_j/z_i)
(2-z_4/z_5-z_5/z_4)\,Z(z)
\right].
\tag{15}
$$

CT是四变量Laurent常数项；12=3!2!。不是把非阿贝尔群积分偷换成不带Weyl密度的相位平均。原Z₆商已在(1)处理，归一Haar不额外乘覆盖次数。

令r_S=∏_{i∈S}z_i，S遍历原16个偶子集。W=4⁻¹⁶∏_S(2+r_S+r_S⁻¹)，M=(4⁸·495)⁻¹h₈(2+z_i+z_i⁻¹)。因此

$$
\int_G Z\,dg=\frac{C}{12D},\qquad
D=4^{24}\cdot495,\qquad
C=\operatorname{CT}P\in\mathbb Z,\quad
0\le C\le12D ,
\quad
\deg_{\pm}P\le(19,19,19,18).
\tag{16}
$$

P是将(15)乘D后的整数Laurent多项式。W每坐标次数≤8，M≤8，Weyl密度分别再加(3,3,3,2)。所以每个变量20阶单位根平均**恰好**提取CT：次数范围内没有非零20的倍数。这只用于提取常数项，不需恢复全部Fourier系数。

代码在两个素数2013265921、1000000021上执行同一有限和。素性用遍历至整数平方根的试除核实，20阶根用阶的素因子2、5核实；全部乘加检查不超过64位整数范围。两素数乘积大于12D，故中国剩余定理结合(16)的界唯一恢复C，无舍入误差。得到

$$
\boxed{\displaystyle
\int_G Z(g)\,dg
=\frac{261746352167}{17416264183971840}
\simeq1.5028845991431662\times10^{-5}>0 .}
\tag{17}
$$

另外用20⁴及23⁴浮点Weyl求积独立核数值，两者包含全部球面矩；浮点结果不是精确性的依据，整数CT才是。这里没有拟合或蒙特卡洛误差。

## 6. 核验及真正补上的条件

第一组复用已保存入口：七份原完整非交换群配置、实际64维H和32维辅助配对，各四个E检查全局帧、Pf与实际Weyl行列式；最大Pf绝对误差约5.91×10⁻¹⁶。另检验靠近单位元的相关群Gram、表示乘法、原Z₆商和物理W零点。群Gram正性由(13)证明，数值只核对象映射。

第二组完成秩四完整群Haar与全S⁹的精确积分及独立浮点检查。未使用675八模式截断、固定E或不受控抽样。对特殊诊断的精确结果不签收一般动态规范模型。

本轮减少的自由：在此整个群部门上，辅助平均、Weyl权重、群正核和Gauss投影必须来自同一个已固定对象，不能分别另配。旧615的“原群不属Spin(9)充分条件”不再阻碍这个**平坦一格部门**的非零性；一般多格问题仍在。

$$
\text{整个原群上的正型函数与非零Haar平均}
\quad\not\Longrightarrow\quad
\text{673完整动态物理反射正性、原 }H_F\text{ 身份或连续时空}.
\tag{18}
$$

## 7. 接续下一项

[688入口](../../688/drafts/STATUS.md)：恢复两个时间片及两个独立原时间链路，先在仅时间传播、空间平凡的原Wilson部门核实际配对是否成为相邻球面之间的同一核k(E,F)=((1+E·F)/2)⁸。若成立，检验完整球面积分与群holonomy的转移表达；实际双边界身份及所有时间归一因子必须保留。

这比另造一个群正核更严格：下一轮必须用原矩阵推出来。一般空间传播、非零质量、原H_b动态配置、真正物理观测与反射仍分别验收。没有更改联合目标，四个原分支继续区分。

