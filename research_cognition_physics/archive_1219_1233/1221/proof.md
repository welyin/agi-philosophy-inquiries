# 1221候选证明：原Gauss匹配比较态与边界关联能源

2026-10-15。已通过[独立审阅](independent_review.md)，最终范围以[验收](acceptance.json)为准；不据此完成TG3。解析对象是617、637、642—644的原完整固定图Hamiltonian及同一切分，不是1219的移动装置。

## 1. 继承对象与准确增量

固定原有限三维图、严格正背景几何、原商群
$G=(SU(3)\times SU(2)\times U(1)_Q)/\mathbb Z_6$、全部H⁵标量、CAR物种及原耦合。保留589的非对角动能/梯度/磁势、598的全部质量与跳跃；采用603完整正常Gibbs态$\sigma_\beta=e^{-\beta H}/Z$、$\beta>0$。几何可在固定图上的紧正光滑参数域内变化。$J$为617事前固定的真实切边等距，且与这些背景参数无关。

复用[Wolf—Verstraete—Hastings—Cirac，式(3)](https://arxiv.org/html/0704.3906v2)：普通张量积Gibbs态的区域互信息受边界比较能源控制。其自由能变分方法不是本轮新定理。原Gauss空间不能任取两边缘的直积；本轮补的是**保原匹配的比较态、原无界材料的有限形式能源以及同一来源的边界消项**。无界边项不能自动套用该文的有界自旋面积上界。

[Casini—Huerta—Rosabal，§3.2](https://arxiv.org/html/1312.1183)已有带中心的区域代数和分块熵。表示索引及中心项沿636、637保留，不把它们都解释成可消去的物理关联。

## 2. 精确匹配与非线性比较态

对有限条切边的紧群作用作Peter–Weyl分解。区域内部Gauss约束已施加，切口表示标签合记$\lambda$，$d_\lambda$为所有切边表示维数之积。固定原CAR次序，可写

$$
\mathcal H_{\rm phys}\simeq
\bigoplus_\lambda(\mathcal H_{a,\lambda}\widehat\otimes\mathcal H_{b,\lambda}),
\qquad
J\big|_\lambda=I_{ab}\widehat\otimes|\Omega_\lambda\rangle,
\qquad
|\Omega_\lambda\rangle=d_\lambda^{-1/2}\sum_{i=1}^{d_\lambda}|i\rangle|\bar i\rangle.
\tag{1}
$$

这里$\mathcal H_a,\mathcal H_b$是重数空间，含区域物质和未切自由度；式(1)不是两区域可独立准备的权限。帽号提醒原CAR分次组合。每个$P_\lambda$为偶投影，与总宇称交换；不把$\lambda$擅自扩成新的局域宇称超选规则。

对任意正常、总偶、有限完整形式能源的$\rho$，定义

$$
\bar\rho=\sum_\lambda P_\lambda\rho P_\lambda,
\quad p_\lambda=\operatorname{Tr}(P_\lambda\rho),
\quad \rho_\lambda=P_\lambda\rho P_\lambda/p_\lambda,
\quad
\tau[\rho]=\bigoplus_{p_\lambda>0}p_\lambda
   (\rho_{a,\lambda}\widehat\otimes\rho_{b,\lambda}).
\tag{2}
$$

零概率块取零。每个条件边缘都是偶态：由$[\rho_\lambda,\Pi_a\Pi_b]=0$，对另一边部分迹并用其宇称共轭的迹不变性，得到$[\rho_{a,\lambda},\Pi_a]=[\rho_{b,\lambda},\Pi_b]=0$。故其分次乘积是正、正常、总偶的态；在固定排序下可按普通密度矩阵张量积计算。没有要求原$\rho_\lambda$分别与两边宇称对易，原跨区奇×奇偶关联允许存在，比较态会删除其中一部分。

式(2)仍在同一个物理匹配空间中，且

$$
\operatorname{Tr}_B J\tau J^\dagger=
\operatorname{Tr}_B J\rho J^\dagger=\rho_A,
\qquad
\operatorname{Tr}_A J\tau J^\dagger=\rho_B.
\tag{3}
$$

证明：跨$\lambda$项因另一侧正交表示空间而消失；对角块给
$p_\lambda\rho_{a,\lambda}\otimes I_{d_\lambda}/d_\lambda$。原切口概率和强制配对都未改变。

**$\rho\mapsto\tau[\rho]$是非线性的变分比较处方，不是未知输入上的普遍CPTP去关联操作。** 不声称装置已经实现它，不登记为免费抹除、实际可提取功或新认知公理。原完整H通常改变切口表示；不假设$[H,P_\lambda]=0$，也不假设$\sigma_\beta=\bar\sigma_\beta$。

## 3. 原完整H的双向形式控制

637已经证明

$$
0\le K_0=c_T\left(\sum_v-\Delta_{H^5,v}+\sum_e-\Delta_{G,e}\right)
       +u\sum_v U_v\le H+C_*,\qquad c_T,u>0,
\tag{4}
$$

以及$J^\dagger(K_A\otimes I+I\otimes K_B)J=K_0$。$K_A,K_B$为该轮比较能源，切边Casimir各分一半，未复制完整物理动能。本轮需要反向上界，不能仅凭式(4)宣布$\tau$有有限H能量。

令$r_v=d_{H^5}(\phi_v,0)$。原群保持原点，三角不等式给

$$
d(\phi_v,R(U_e)\phi_w)^2\le2(r_v^2+r_w^2).
\tag{5}
$$

574的原势为$U=V/F^2$。其边界球上$V$有严格正下界，且
$F=M\operatorname{sech}^2(r/\sqrt6)$，故$U\ge c\exp(4r/\sqrt6)$在足够大$r$成立。紧内区补一个常数便有$r^2\le C(1+U)$，无需修改原势。

逐项比较原H：

- 节点/规范动能：固定图紧正几何的最大特征值控制完整二次形式，包括非对角电项。
- 原梯度：589的$\bar\gamma^{\mu\nu}\langle u_\mu,u_\nu\rangle$由其最大特征值乘$\sum\|u_\mu\|^2$控制，再用式(5)。不是删除交叉项。
- 原磁势：有限图紧群上的光滑函数有界，含原面间交叉项。
- 节点势：正权有固定上下界。
- 598原全部CAR质量：$\|B_{\rm mass}\|\le\eta\sum w_vU_v+C_\eta$；跳跃是有限CAR与紧群酉矩阵乘有界系数，故有界。

所以同一光滑核上、再取共同闭形式延拓，存在有限常数$a,b$使

$$
H+C_*\le aK_0+bI,
\qquad D(q_H)=D(q_{K_0}).
\tag{6}
$$

常数可依赖固定图、温度之外的原系数和紧正几何域；**未证明随图规模或格距一致**。没有从物理应力中扣掉$C_*$或$b$。

由式(3)和637的区域比较算符，先用非负谱截断再取单调极限：

$$
\operatorname{Tr}\tau K_0=\operatorname{Tr}\rho K_0<\infty,
\qquad
\operatorname{Tr}\tau(H+C_*)\le a\operatorname{Tr}\rho K_0+b<\infty.
\tag{7}
$$

这在原无限维H⁵/紧群材料及全部CAR上成立。637进一步给原两个扩展区域熵有限，故$H(p)$、$\sum p\log d_\lambda$、$\sum p S(\rho_{a,\lambda})$与B侧对应项均有限。$\tau$全局熵由这些非负项控制，$\rho$全局熵由原全热迹/有限能量控制。下节不存在$\infty-\infty$。

## 4. 可以删除的关联、不可偷删的中心与配对

固定这些边缘与概率时，$\log\tau$逐块为
$\log p_\lambda+\log\rho_{a,\lambda}+\log\rho_{b,\lambda}$；因此

$$
\mathcal C(\rho):=D(\rho\Vert\tau)
=S(\tau)-S(\rho)
=D(\rho\Vert\bar\rho)
  +\sum_\lambda p_\lambda I(a:b)_{\rho_\lambda}\ge0.
\tag{8}
$$

支持包含由正块及条件边缘保证；有限熵使各项可按有限块/谱近似取极限。$D(\rho\Vert\bar\rho)\le H(p)<\infty$。式(8)还说明$\tau$是保持这些资料的唯一最大熵比较态，而不是某种真实动力学的平衡态。

扩展空间的普通互信息则是

$$
I(A:B)_{J\rho J^\dagger}
=H(p)+2\sum_\lambda p_\lambda\log d_\lambda+\mathcal C(\rho).
\tag{9}
$$

因此不能把左边全部计成式(8)的可消去关联。中心概率与表示配对基线被保留；跨sector相干通常被删除，必须计入式(8)第一项。636两次切口时单边表示熵是$2\sum p\log\dim R$，与这里$d_\lambda=\prod_{e\in\partial A}\dim R_e$一致。

## 5. 同一完整Gibbs态的边界能源身份

在原共同形式域按实际支撑拆分

$$
H=H_A+H_B+H_\partial.
\tag{10}
$$

$H_A,H_B$为各自区域偶代数中的完整项；跨边/跨面/跨梯度簇项全部留在$H_\partial$。切边的电导数可完整放在原源节点相邻的半边，原同源非对角星形项仍在同侧，不能把617的运输公式误称两边独立Hamiltonian。即使采用另一合法分账，只要真正局域项保留式(3)的边缘期望，下式不变。

区域内部、切口不变的算符由紧群Schur分解逐块作用为$X_{a,\lambda}\otimes I_b$或相反，故$\tau$保留其期望。无界项的结论由式(6)—(7)及其形式控制延拓，而非假设所有区域算符有界。

取$\rho=\sigma_\beta$，$\tau=\tau[\sigma_\beta]$。Gibbs对数给

$$
\boxed{
 D(\sigma_\beta\Vert\tau)+D(\tau\Vert\sigma_\beta)
 =\beta\operatorname{Tr}[(\tau-\sigma_\beta)H]
 =\beta\operatorname{Tr}[(\tau-\sigma_\beta)H_\partial].}
\tag{11}
$$

两相对熵有限，右端因真实H域保证有限，且非负。于是
$\mathcal C(\sigma_\beta)\le\beta\Delta E_\partial$。这不是假设跨区相干不可见，而是把其损失连同块内关联一起收费。$\Delta E_\partial$是比较能源差，不是已经实现的准备功。

若再有逐边有界相互作用和独立面积/密度控制，才可进一步讨论面积上界。原H⁵跨边势无界，式(11)本身不允许直接写$\|H_\partial\|$或认定$\Delta E_\partial\propto A$。

## 6. 完整几何来源：固定态变分与实际变态分开

令$\xi$为原背景参数，固定区域与$J$。589、623的来源及上述逐项界给在紧正域内
$|\partial_\xi q_H|\le C_\xi(q_{K_0}+\|\cdot\|^2)$，所以两态的完整来源期望都有限。局域项及其系数导数仍在同一区域代数，故

$$
\operatorname{Tr}[(\tau-\sigma)\partial_\xi H]
=\operatorname{Tr}[(\tau-\sigma)\partial_\xi H_\partial].
\tag{12}
$$

这是同一两份固定物理态的虚变分来源差，不是冻结参考后冒称真实热响应。若比较对本身随$\xi$变化，在相应带能量权迹可微条件下必须使用

$$
\frac{d}{d\xi}\operatorname{Tr}[(\tau_\xi-\sigma_\xi)H_\xi]
=\operatorname{Tr}[(\tau_\xi-\sigma_\xi)\partial_\xi H_\xi]
 +\operatorname{Tr}[\partial_\xi(\tau_\xi-\sigma_\xi)H_\xi].
\tag{13}
$$

式(11)逐个参数成立；本轮不从有限一阶能量单独推出原无限sector族$\tau_\xi$的带权迹可微，也不把式(12)替代式(13)。移动区域/$J$、动态量子几何及内部面积仪器另须验收。相同参考/来源的要求复用639、644、645、843，未追加新自由应力。

## 7. 真实原材料中的有限诊断，不冒充完整热谱

取原至少$3^3$周期图上的一个(12)方环。每节点使用同一非零光滑紧支径向H⁵包
$f(r)\propto\exp[-1/(1-r^2)]$（$r<1$，外部零）；所有其它链路取常数函数。环函数取平凡角色与原弱伴随角色$\chi_1=|\operatorname{Tr}W_\square|^2-1$。它们严格Gauss不变，且弱伴随标签$(0,0,2,0)$满足原Z₆商条件。

保留环上一对跨区相邻节点的原sterile $\nu_R$两spin模式，即旧64模式字典的30、31、62、63；其它CAR取真空。这些模式均是原群singlet。四个模式全部Fock空间维数16，含奇/偶sector；诊断密度与总宇称对易，并未增加局域宇称超选。两环标签乘16维给32个真实物理向量。每个属于原全H所有固定能量矩域。

两个区域各含方环同侧的两个顶点。原$J$给平凡sector Schmidt秩1、伴随sector Schmidt秩$3^2=9$。有限重数之外的边界配对保持；不是另装辅助qubit。

令$P$为该真实物理子空间投影。沿617均匀形状$\bar\gamma=e^{tS}$、$\epsilon=.73$、$t=.8$、$\psi=e^{.12}$，569原冻结弱系数$b_w$及磁伴随系数$\alpha_w$不变，644声明的跳跃$t_\nu=.17$不变。矩阵元给

$$
PHP=C(\gamma)I+
\left[\begin{pmatrix}0&0\\0&2\kappa\end{pmatrix}
-g\begin{pmatrix}0&1\\1&1\end{pmatrix}\right]\otimes I_{16}
+I_2\otimes t_\nu\sum_{s=1}^2(c_{A,s}^\dagger c_{B,s}+\mathrm{h.c.}),
\tag{14}
$$

$$
\kappa=\frac{b_w e^{-2\sigma}}\epsilon
 (2\bar\gamma_{11}+2\bar\gamma_{22}-2\bar\gamma_{12}),
\qquad
g=\frac{\alpha_w e^{-2\sigma}}\epsilon\bar\gamma_{33}.
\tag{15}
$$

各项依据：

1. 环原点两条出边的Gauss关系给$P_1\cdot P_2=-C_2$，其余出边配常数函数；伴随$C_2=2$。故式(15)保原剪切交叉项。
2. 原磁势的伴随项是$\alpha_w(3-\chi_1)$。Haar角色正交及$1\otimes1=0\oplus1\oplus2$给式(14)融合矩阵；代码直接用SU(2) Haar类积分复核。其余表示的非零颜色/超荷使相关非平凡矩阵元消失，其余面有未出现于该环的独立边，Haar平均给共同常数；不同面的O交叉项同理为零。
3. 标量动能和节点势在同一个$f$上给共同常数。每个梯度星的其它节点径向分布可独立吸收各原紧群转动，故其积分与环holonomy无关，也只给共同常数；非对角梯度没有从全H删除。
4. 598全部Dirac/Majorana矩阵元是$\phi^A/\sqrt F$的奇函数，径向包使均值为零。外部带荷CAR真空及所选中性模式给原sterile跳跃；其它跳跃的矩阵元为零。未把完整H的质量项删掉，作用出P的分量一般非零。

$C(\gamma)$有限但本轮不计算；它在同归一化两态的能源差、来源差及压缩Gibbs态中严格相消，**不等于物理真空能源可忽略**。压缩Gibbs态$e^{-\beta PHP}/Z_P$只是同原H矩阵元的Galerkin诊断；$P$不是H不变空间，没有把它当原603完整Gibbs或其误差受控近似。完整Gibbs结论由第3—6节承担。

## 8. 验收范围与后续门槛

[check.py](check.py)核四组：原商群/角色/完整CAR矩阵元；匹配边缘及熵能源分账；原完整来源的bulk消项与变态项；比较处方非线性。数据见[results.json](results.json)。只验证声明的有限诊断；解析域证明不能由这些数值替代。

本轮接起的是原完整Gauss材料中的**有限关联能源—边界来源**，没有独立面积仪器，没有$\delta S=s\delta A$，没有修复原候选的全部连续极限，也不与1219—1220机械材料拼成同一实现。最有价值的下一问题是：当同材料的实际空间区域和面积读口已经就绪，式(11)的边界能差能否由事前声明的有限边界资源/面积读数共同控制；若只能事后拟合比例或追加任意边界态，停止此支。
