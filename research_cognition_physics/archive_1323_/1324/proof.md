# 1324：有限相干参考、有效转换与不可省略的反作用

2026-10-19。前件见[准入](../_admission/round1324_admission.md)。本轮沿1323同一模式父作用建立有限有效接口；未把它升级为局域Higgs场、全部规范物质或现实质量。新增公理0。

## 1. 分层输入与实际共同对象

**认知动机。** 共同参考可以使参与者执行相容操作，但使用后的参考、记录与来源仍在整体内。此动机本身不规定参考必为相干态或下列Hamiltonian。

**候选操作合同。** 继承1323的全未知参考单份交接、无隐藏量子旁路及同钟短时余额一致性。这些尚非已采用公理。一般N的交织结论直接复用；本轮另登记在所选有限窗口中实施该合同。

**物理类别。** 两个玻色参考模式 $a_0,a_1$，单费米部门

$$
\mathcal S=(R\oplus L_0\oplus L_1)\otimes K_n,\qquad
T_\alpha=|L_\alpha\rangle\langle R|\otimes M,\qquad N_b=\sum_\alpha a_\alpha^\dagger a_\alpha.
$$

左右味资料及任意未参与伙伴Z可以未知。单费米式是1323偶CAR算符 $l^\dagger r$ 的压缩，不凭此宣称全部多费米物种已实现。令真空上的 $N_b^{-1/2}=0$，取

$$
B_\alpha=a_\alpha N_b^{-1/2},\quad
B_\alpha^\dagger=N_b^{-1/2}a_\alpha^\dagger
=a_\alpha^\dagger(N_b+1)^{-1/2},\quad
V=\sum_\alpha(T_\alpha\otimes B_\alpha+T_\alpha^\dagger\otimes B_\alpha^\dagger).
\tag{1}
$$

守恒量 $K=N_b+n_L$，其中 $n_L$ 为两个左角色的总投影。实际完整空间为 $K=0,\ldots,L+1$ 的全部部门，不只截断 $N_b\le L$；从左端和参考L的输入，可以到右端和参考L+1。暗方向、失败及反向过程都保留。每个K≥1块的亮转换是 $W_K\otimes M$，故 $\|V\|\le\|M\|$；K=0只有右真空。

自由Hamiltonian先独立登记

$$
H_0=E_R n_R+E_L n_L+\omega N_b.
\tag{2}
$$

在明确工作片 $E_L=E_R+\omega$，有 $[H_0,V]=0$。以下转换在该共同相互作用绘景中比较；回到物理绘景对两种输出施同一 $e^{-itH_0}$，参考相位旋转不被删除。固定有限窗口中可选正的自由能基底，不要求由相互作用项单独定义正能量。

**初态与权限。** 对单位 $u\in\mathbb C^2$、$\alpha=\sqrt\nu e^{i\phi}$，实际参考准备为

$$
|\chi\rangle=Z_L^{-1/2}\sum_{k=0}^{L}
\frac{\alpha^k}{\sqrt{k!}}|k,u\rangle,\qquad
Z_L=\sum_{k=0}^{L}\frac{\nu^k}{k!}.
\tag{3}
$$

$|k,u\rangle$ 是沿u的归一化对称k粒子态。准备的数目、纯度、方向与相位是输入资源；相位相对已登记内部参考解释，不由固定数目态免费取得。初态为 $\rho_{\mathcal SZ}\otimes|\chi\rangle\langle\chi|$。若准备器仍与参考相关，须另外扩大初态与误差模型，本定理不代签。有限钟、开关、相干读口为已声明权限；未构造自治制备工厂，也未赋予任意规范带荷操作。

## 2. 有效矩阵与精确生成元误差

记 $p_k=\nu^k/(k!Z_L)$、$p_0=1/Z_L$，以及

$$
s_L=\sqrt\nu\sum_{k=0}^{L-1}\frac{p_k}{\sqrt{k+1}},\qquad
c_L=\sum_{k=0}^{L}\frac{p_k}{k+1},\qquad
b_\alpha=\langle B_\alpha\rangle=e^{i\phi}s_Lu_\alpha.
\tag{4}
$$

定义物质有效生成元及实际输入嵌入

$$
h=\sum_\alpha(T_\alpha b_\alpha+T_\alpha^\dagger b_\alpha^*),\qquad
E\psi=\psi\otimes|\chi\rangle.
\tag{5}
$$

这里 $h=E^\dagger VE$ 是精确恒等式；并且
$\partial_t\rho_{\mathcal S}|_{0}=-i[h,\rho_{\mathcal S}]$。
它提供一次短时相干校准的操作口径，不将有限时概率先拟合后再宣称预测。

**命题1。** 对L≥1、ν>0和任意M，

$$
\|VE-Eh\|=\eta_L,\qquad
\eta_L=\|M\|\sqrt{\max(1-s_L^2,c_L)}.
\tag{6}
$$

证明：SU(2)协变性允许选u=(1,0)计算，不限制未知物质方向。令 $\delta B_\alpha=B_\alpha-b_\alpha$。右输入的误差输出在左部门，误差Gram为

$$
\sum_\alpha\langle\delta B_\alpha^\dagger\delta B_\alpha\rangle M^\dagger M
=(1-p_0-s_L^2)M^\dagger M.
$$

左输入的误差输出在右部门，方向Gram的两本征值是 $1-s_L^2$ 与 $c_L$，分别对应沿u和垂直u的左态，并乘 $MM^\dagger$。交叉项为零；左右输出部门正交。最大本征值即(6)。一般u时方向矩阵为

$$
(1-s_L^2)uu^\dagger+c_L(I-uu^\dagger).
$$

不能删除横向左态：在大参考区它往往决定最大的全输入误差。式(6)是所有输入及伙伴的统一范数，不是只对右初态的均值检验。

## 3. 有限时间、完整后态与任意伙伴

**命题2。** 对任意实t，

$$
\|e^{-itV}E-Ee^{-ith}\|\le |t|\eta_L.
\tag{7}
$$

矩形Duhamel恒等式直接给

$$
e^{-itV}E-Ee^{-ith}
=-i\int_0^t e^{-i(t-s)V}(VE-Eh)e^{-ish}\,ds.
$$

它避免错误假定 $h\otimes I$ 保持整个K截止空间：最高右态再受h可能越界。实际E的像是合法子空间，矩形公式无需这种假定。也可在完整K空间压缩 $h\otimes I$，其在E的参考支撑上精确实现同一比较过程。

张量任意伙伴不会改变算符范数；由纯态距离和凸性，实际完整物质—参考—伙伴输出与

$$
(e^{-ith}\otimes I_Z)\rho_{\mathcal SZ}(e^{ith}\otimes I_Z)
\otimes|\chi\rangle\langle\chi|
$$

的半迹距离，乃至这两个从物质到物质—参考的通道的半diamond距离，不超过

$$
\epsilon_L(t)=\min(1,|t|\eta_L).
\tag{8}
$$

这证明有限工作窗中的有效转换与共同后态接口。它不是无限次无误差催化：重复使用会累积误差并改变参考，不能每次无来源地重新准备χ。

所有和与空间都是有限的。无须先规定微观连续、物理最小尺度或完成无截止理论。L和ν改变时必须重新登记实际资源与(8)，不能从本轮一个数值窗口签任意高能外推。

## 4. 来源、反作用与控制

精确父守恒给

$$
\Delta\langle N_b\rangle=-\Delta\langle n_L\rangle,\qquad
\Delta\langle n_R\rangle=-\Delta\langle n_L\rangle.
\tag{9}
$$

参考初始数目为 $\mu_L=\nu(1-p_L)$。同一自由能账为

$$
\Delta\langle H_0\rangle
=(E_L-E_R-\omega)\Delta\langle n_L\rangle=0.
\tag{10}
$$

完整H与V的期望也守恒。真实独立来源为
$\partial_{E_R}H=n_R$、$\partial_{E_L}H=n_L$、$\partial_\omega H=N_b$，M实虚部给同一转换算符来源。先保持这些源独立，再施加共振工作条件；不能一边变ω一边暗调E_L而漏记源。

比较过程把参考态冻结，数目期望一般不等于实际值。(8)不意味着无资源消耗；即便有限空间可用 $\|N_b\|\epsilon_L$ 作粗界，这个界也不替代(9)。正确的有效记账是由有界投影n_L预测实际转移，再按(9)计参考反作用：

$$
|\Delta\langle N_b\rangle+\Delta\langle n_L\rangle_{\rm eff}|
\le\epsilon_L(t).
\tag{11}
$$

此界对同一初态成立，不给比较密度矩阵额外修改数目均值。物质从右到左近乎完全转换时，参考实际失去近一个量子，仍可能在迹距离上接近原态。

开关V是公开控制权限。瞬时开关若采用该理想描述，功项分别为 $\langle V\rangle_{\rm in}$ 与 $-\langle V\rangle_{\rm out}$；因V在工作段守恒，净值抵消。对事前右输入二者均为零，不能因此宣布制备、时钟或相干参考无需资源。若研究有限时控制器，须补其具体Hamiltonian；不把未建控制器当作当前代数定理已实现，也不把所有设备制造列为本轮有限比较的门槛。

## 5. 普通线性菜单：同一比较口径

将式(1)的B换为a，定义 $V_{\rm bare}$。它不满足1323跨余额同动作合同；本节检验采用该合同的有限后果。保同一完整K窗口、参考准备、时钟、物质和未知伙伴。记

$$
a_L=\sqrt\nu(1-p_L),\qquad
h_{\rm bare}=\sum_\alpha(T_\alpha e^{i\phi}a_Lu_\alpha+\mathrm{h.c.}).
\tag{12}
$$

对应精确残差为

$$
\eta_{\rm bare}
=\|M_{\rm bare}\|\sqrt{1+\nu p_L(1-p_L)}.
\tag{13}
$$

证明：右残差为 $(\mu_L-a_L^2)M^\dagger M$；沿u左残差为 $(1+\mu_L-a_L^2)MM^\dagger$，横向为 $MM^\dagger$。代入 $\mu_L=\nu(1-p_L)$ 即得。Duhamel与(9)—(11)逐字成立，换用对应η。

不能用固定绝对t的裸残差不随ν缩小，否认普通相干驱动有有效极限；实际转移时间可随 $1/(\|M\|\sqrt\nu)$ 缩短。本轮只比较已固定的同钟任务。

## 6. 固定校准与未拟合的有限预测

事前取 $L=400$、$\nu_0=64$、$\nu_1=256$、$m_0=0.1$、$T=\pi/(2m_0)$。用ν_0的短时完整相干生成元校准，选择

$$
M_{\rm norm}=g_{\rm norm}I_n,\quad g_{\rm norm}=m_0/s_L(\nu_0),\qquad
M_{\rm bare}=g_{\rm bare}I_n,\quad g_{\rm bare}=m_0/a_L(\nu_0).
\tag{14}
$$

两父模型系数不同，各自只校准一次；随后不改L、T、g或准备规则。所有味接受相同倍数，不以观测经验质量设M。这里n=3、4、5均可张量接续；不把n消去再当成代际新定理。

对初始右端，实际概率可从完整父块精确求出：

$$
P_{\rm norm}(\nu)=(1-p_0)\sin^2(g_{\rm norm}T),\qquad
P_{\rm bare}(\nu)=\sum_{k=0}^L p_k\sin^2(g_{\rm bare}T\sqrt k).
\tag{15}
$$

有效概率分别为 $\sin^2(g_{\rm norm}s_LT)$ 与
$\sin^2(g_{\rm bare}a_LT)$，误差由各自(8)控制。两者在校准点拥有同一精确短时生成元，不保证同一有限时概率。

数值见[作者结果](results.json)及[独立结果](independent_results.json)。筛查预计ν_1时实际概率约0.999990与0.009562；有效概率约0.999995与0，完整输出误差上界约0.09837与0.19635。因此即从保守界也能留下大于0.70的概率间隔。正式数值由两套保存结果签核，不靠均值近似冒充实际演化。

区分的是有限模式实验而非现有Higgs观测。若理想读口实现为独立重复Bernoulli记录，在实际可登记的独立制备、无系统误差前件下，每菜单1000次、每个频率误差0.05的联合失败界为 $4e^{-2(1000)(0.05)^2}<0.027$；需支付2000次准备及读口资源。此成熟集中界仅解释有限可辨性，不授予免费独立重复或真实实验设备。

## 7. 相位、味响应与尚未跨过的物理字典

有效转换块为

$$
C_{\rm eff}(\nu,u)=e^{i\phi}s_L(\nu)\,u\otimes M.
\tag{16}
$$

同一参考变化时，所有非零奇异值乘同一 $s_L(\nu_1)/s_L(\nu_0)$，左右味混合方向不变；一般M仍自由。这是有限准备响应的共同关系。它不选择现实参数，也不是独立Higgs源 $\delta S/\delta H$：ν和u首先为已登记准备参数，没有空间动量、极点或场归一字典。误差界本身不能直接微分以制造导数响应定理。

若相位平均，$\langle B_\alpha\rangle=0$而式(15)的右到左概率不变；各相位下h的平均为零不等于平均演化恒等。用零均值代替真实过程会丢掉参考相关。此诊断说明相干准备是有效幺正接口的实质资源，不是可随手略去的数学选择。

本轮完成的是同一1323父作用的“实际有限量子参考—有效转换—完整联合误差—精确资源反作用—固定校准后的区别”。参考低退化、量子相干和Rabi本身已有成熟研究，参见[Åberg](https://arxiv.org/abs/1304.1060)与[van Luijk等](https://quantum-journal.org/papers/q-2023-11-06-1166/)。本轮不声称首创这些一般现象。

下一步应优先判断相反参考荷的多个右物种能否在同一左端和共同父作用中实现，及其产生什么真正跨部门限制；不能停在更多ν、更多L或单模式仪器修补。单下降端口容量限制仍不是全部右物种的分类，完整局域场源与现实质量仍待独立连接。
