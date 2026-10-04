# 第680轮：完整物理代数中的质量约化与严格归一边界

日期：2026-10-02。接[679](../../research_note_679.md)、[680已执行入口](mass_congruence_entry.md)。[代码](../joint_mass_reflection_congruence.py)、[结果](../joint_mass_reflection_congruence_results.json)、[核验](../research_round_680_checks.json)、[联合账](../unified_physics_condition_ledger_680.md)。两组新核验、十四式；主代理审查，无独立代理审查。

## 1. 去重、问题与实际增量

上一目标轮完成679、更新权威索引并执行680入口，属于有效进展。本次核对README、research_direction、RESEARCH_STATE、最新报告、结果及发布核验；未发现运行中的Python研究。回查668、669、673、677、679及相关历史检索。

668已经用半区质量指数证明自由背景下的正性继承；669已有动态标量热矩、有限质量指数与可能的归一零点。这些不算本轮新定理。本轮完成的是：**这份操作在673完整动态规范／双Haar候选的物理代数上保支撑、保Gauss且可逆，从而把全部原局部质量下的RP统一约化到同一个无质量泛函。**

另给精确Grassmann反例：即便原泛函反射正、Z(0)=1，质量变换始终可逆，也可能在有限λ处Z(λ)=0。因此不能顺便删掉严格归一条件。反例只否定这一一般逻辑蕴涵，不是原候选或认知原则的反例。

|层次|内容|
|---|---|
|认知动机|原质量、可观测操作、规范约束与概率归一须在同一过程中相容|
|继承输入|原有限盒、原16通道、实际物理Y、Pφ、H_b、S⁹、完整双Haar与679反射|
|解析增量|完整物理代数上的可逆质量合同；全多项式RP双向等价；严格归一不被该等价保证|
|数值验证|全部264项原质量逐项加入／逆序去除与0／2／4／6来源；精确有理Grassmann反例|
|保留边界|无质量完整候选RP、指定归一、原H_F身份、共同连续及量子GR|

空间382—386、425、522—523直接按679 §1复用。本轮不重证旧空间合同，不恢复384已消去输入，不叠加386／425，也不重开523已完成的条件性共同实现。

反射配对作用的正性插入方法是成熟工具，见[Kikukawa—Usui，§II、III、V](https://arxiv.org/html/1005.3751v3#S5)。原文自由／非规范范围不替代这里的动态规范正性。我们只检验它在本项目完整对象上的双向合同。

## 2. 哪一个代数，哪些参数保持不变

固定原有限盒、正热时τ、全部原H_b参数及几何、原Wilson核和辅助测度。λ是当前候选中的质量插入参数；改变λ时不偷偷重选H_b或参考态。若现实建模要求质量参数连同H_b势一起改变，须重新匹配那个联合族，不能直接套本轮“同一基线”的等价。

取正半区物理Grassmann多项式，系数允许为原正半区配置的多项式增长函数；规范变换下整个表达式为标量。记

$$
\mathcal A_+^{\rm phys}
=\left\{f(q_+,E_+,Y_+):\ 
\text{有限Grassmann次数、多项式增长系数、原局部Gauss不变}\right\}.
\tag{1}
$$

可取其中偶子代数；下面变换保费米奇偶，结论分别成立。源列本身带荷并不自动属于(1)。使用原物理Y的时间标签，不能换成661已证明一般不保时间支撑的裸Φ变量。

令Π±是物理字段正负时间片的选择矩阵。原质量逐格点，故

$$
P_\phi=P_{\phi,+}+P_{\phi,-},\qquad
P_{\phi,\pm}=\Pi_\pm P_\phi\Pi_\pm,\qquad
V_\lambda=V_{\lambda,+}+\Theta V_{\lambda,+},\quad
V_{\lambda,+}=\tfrac\lambda2Y^{\mathsf T}P_{\phi,+}Y .
\tag{2}
$$

link反射没有固定时间片；原Pφ=diag(p,−p*)的反射负号由679承担。实际独立Gauss变换G_Y保时间片，且G_YᵀP'φ,+G_Y=Pφ,+，所以V₊为真正半区规范标量，不依赖跨半区裸谱投影的支撑。

## 3. 可逆性及全部平均的可积性

正半区含m₊个物理Grassmann生成元。V₊偶、二次，指数有限：

$$
R_\lambda f=e^{V_{\lambda,+}}f,\qquad
e^{V_{\lambda,+}}=\sum_{k=0}^{\lfloor m_+/2\rfloor}
\frac{V_{\lambda,+}^k}{k!},\qquad
R_\lambda^{-1}=R_{-\lambda}.
\tag{3}
$$

因此正负λ都不会扩大时间支撑，且不产生无限阶标量指数尾。原质量系数在669的x=φ/√F变量中至多线性；任一原多项式增长系数乘上(3)后仍为有限阶增长。于是

$$
R_{\pm\lambda}:\mathcal A_+^{\rm phys}\longrightarrow
\mathcal A_+^{\rm phys},\qquad
R_{-\lambda}R_\lambda=I,\qquad
R_\lambda(g\!\cdot\! f)=g\!\cdot\!R_\lambda f .
\tag{4}
$$

这也可扩展到各处不同、但满足原质量反射与Gauss合同的有限Dirac／Majorana系数。它不选择真实耦合，也不把所有质量参数当作可以与H_b独立变化的自然常数。

使用673全部热核、S⁹、双Haar平均dϖτ。若w(q)为其原多项式增长控制函数，对每个固定f、g、λ，存在有限p和C，使质量插入后每个Berezin系数受C[w(q₊)+w(q₋)]^p控制。故669／673直接给

$$
\int d\varpi_\tau\,
\left|I_{0,b}\!\left(\Theta(R_\lambda f)R_\lambda g\right)\right|
<\infty,\qquad
M_p(\tau)=\int w(q)^p k_\tau(q,q)\,d\mu(q)<\infty .
\tag{5}
$$

这里没有假定未积分Grassmann权重为正。有限指数逐项可积，因此可合法保留全部原平均再重排。也不需要物理N、质量或配分函数可逆。

## 4. 完整RP的双向约化

在同一原平均上令𝔐λ(f)=∫dϖτ I₀,b(e^Vλ f)。因V₊偶，Θ反线性且反转乘积次序，得到

$$
Q_\lambda(f,g):=\mathfrak I_\lambda(\Theta f\,g)
=\mathfrak I_0\!\left(\Theta(R_\lambda f)R_\lambda g\right)
=Q_0(R_\lambda f,R_\lambda g).
\tag{6}
$$

前向是668—669熟悉的插入机制；当前新增的是(1)—(5)保证其逆也在完整动态Gauss物理代数内。于是

$$
Q_\lambda\ge0\ \text{于全部 }\mathcal A_+^{\rm phys}
\quad\Longleftrightarrow\quad
Q_0\ge0\ \text{于全部 }\mathcal A_+^{\rm phys};
\qquad
Q_\lambda(R_{-\lambda}f_0,R_{-\lambda}f_0)=Q_0(f_0,f_0).
\tag{7}
$$

因此，若原无质量完整泛函有负方向，它在全部有限原质量下都有对应的物理负方向。该测试一般含更高次来源，只看某个二点Gram或无来源权重可能漏掉它。反之，一旦同一无质量完整对象的RP成立，全部允许质量无需逐一重证RP。

672固定边界的负向量尚未经过完整Gauss和S⁹平均，不能用它直接触发(7)来宣布候选失败。当前无质量Q₀的真实正负性仍开放。

## 5. 这不是“质量没有物理作用”

Rλ是可逆线性乘法，不是保单位的代数自同构，也不自动与物理时间演化相容：

$$
R_\lambda1=e^{V_{\lambda,+}},\qquad
(R_\lambda f)(R_\lambda g)=e^{2V_{\lambda,+}}fg
\ne R_\lambda(fg)\ \text{一般成立}.
\tag{8}
$$

如果Q₀确实正，则(6)可在零范数商空间上给等距对应；但真空向量、被表示的可观测量及时间操作还要同时输送。它不证明不同质量具有相同粒子谱，更不证明候选已等于原H_F。这里压缩的是RP检验条件，不是删掉质量部门。

## 6. 严格归一不能由可逆性领取

指定物理总权重为

$$
Z_\lambda=\mathfrak I_\lambda(1)
=Q_0(R_\lambda1,R_\lambda1).
\tag{9}
$$

Rλ1在多项式代数中非零，不代表它在Q₀的零范数商中非零。给两个正半区Grassmann元x₁、x₂及其反射元y₁、y₂；Θxᵢ=yᵢ，并反转乘积。设a=x₁x₂、b=Θa=y₂y₁，规定∫ba=1。a、b为偶且平方为零，取

$$
\Omega_0(F)=\int(1-a)(1-b)F,\qquad
V_\lambda=\lambda(a+b),\qquad
R_\lambda=1+\lambda a,\quad R_\lambda^{-1}=1-\lambda a .
\tag{10}
$$

在正半区基(1,x₁,x₂,a)上，直接Berezin积分给完整Gram：

$$
Q_\lambda=v_\lambda v_\lambda^\dagger,\qquad
v_\lambda=(1-\lambda,0,0,-1)^{\mathsf T},\qquad
Z_\lambda=(1-\lambda)^2 .
\tag{11}
$$

所以对全部实λ及全部复系数多项式反射正，Z₀=1；但λ=1时Z₁=0，尽管R₁仍可逆。此时Q₁(a,a)=1，泛函不是恒为零。代码以有理数进行完整外代数乘法、反射及积分，不以浮点特征值猜正性。

这个例子满足所讨论的一般质量插入与反射条件，但不是原16通道／H_b模型。它只严格反驳“RP＋原归一＋可逆半区插入⇒所有质量处严格归一”的泛化；原候选Zλ仍需自己的证据。

作为反向测试，另取密度1+2a+2b+ba，其偶基Gram为[[1,2],[2,1]]。若f₀=1−a，则

$$
Z_\lambda=1+4\lambda+\lambda^2,\qquad
f_\lambda=R_{-\lambda}f_0=1-(1+\lambda)a,\qquad
Q_\lambda(f_\lambda,f_\lambda)=-2 .
\tag{12}
$$

例如λ=0、1、3的标量权重均为正，负方向却保留。这个有限控制例再次说明低次子空间或正标量权重不足以代替(7)的全代数条件。

## 7. 原完整质量的复算

正式第一组不再只挑入口中的两个质量单项：在原2×2非平坦全群配置、完整16通道、变化E和φ上，实际Pφ有264个非零上三角项，正负半区各132项，全部加入并逆序去除。只在这个明确非奇异的数值样本使用逆核。

对物理收缩C=ΦN⁻¹Φᵀ及权重w=Pf N，加入原αYᵢYⱼ对应的核秩二更新。沿当前Pfaffian排序：

$$
d=1-\alpha C_{ij},\qquad
w'=wd,\qquad
C'=C+\frac{\alpha}{d}
\left(C_{:i}C_{:j}^{\mathsf T}-C_{:j}C_{:i}^{\mathsf T}\right).
\tag{13}
$$

这是Pfaffian／逆矩阵更新恒等式，不是新代数定理。样本全路径最小|d|约0.98837。重建的质量矩阵与原λPφ逐项一致，没有改写质量模型。原完整增广Pfaffian作为独立对照：

$$
C_Z=\operatorname{Pf}
\begin{pmatrix}N&\Phi^{\mathsf T}Z\\-Z^{\mathsf T}\Phi&0\end{pmatrix}
=w\,\operatorname{Pf}(Z^{\mathsf T}CZ)
\quad\text{在该可逆样本上}.
\tag{14}
$$

完整权重相对误差3.27×10⁻¹⁵，物理收缩最大分量误差1.64×10⁻¹⁵；0／2／4／6来源相对误差均小于3.24×10⁻¹⁵。逆序撤回全部264项后，权重相对误差2.67×10⁻¹⁵、收缩误差1.12×10⁻¹⁶。

代码第二组为上述精确Grassmann正例与逻辑反例。未数值计算完整Haar／S⁹积分，也不从固定配置数值推断RP。一般奇异、变秩及零权重范围由有限多项式证明承担，不把数值逆核的前提升级为物理假设。

## 8. 合并后的条件与下一步

本轮将C17原质量与C14 Gauss、C16完整RP、C19归一及C01物理代数合并：在固定共同基线及上述封闭代数内，RP只需判定一次无质量完整对象；严格归一保留为独立条件。对固定有限软调节同样适用，677的来源极限继续有效；不据此宣称τ、体积、几何来源或连续重整化的联合控制。

原有限图H_F、指定连续物质、给定作用的经典几何、手征辅助候选四分支依然未完全识别。没有新增认知公理，没有完成整个统一目标。

接[681](../../681/drafts/STATUS.md)：优先寻找可与原无质量物理泛函及全部来源直接匹配的动态规范成熟构造，核清其反射／时间和辅助减除条件；同时审查完整Gauss正性真正需要什么见证。近期文献只按已证范围接入，不将不同流、不同边界或低维背景计算直接改名为原候选。停止质量调参和重复小Gram扫描；认知设计后置，目标不变。
