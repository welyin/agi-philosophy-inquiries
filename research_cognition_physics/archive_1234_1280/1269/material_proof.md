# 1269材料分项：体积型杆读取的内部流作用与来源闭合

2026-10-16。按[准入](admission.md)交付1个作者组。接[1268材料结论](../1268/material_proof.md)，不修改冻结模型，不采用CT，不新增认知公理。

**净结论：**体积型杆耦合可以由原守恒标签流的内部辛结构承担。明确的密度依赖作用给出单位速率内部钟、保持的指针动作、完整变分来源和精确交换抵销；它在齐次部门接回1268的同一个读取Hamiltonian。在均匀共动jet上，新增动态控制不恢复1268的密度—杆双向主部混合，但零速传输块仍可能有Jordan链。该有限材料连接不等于完整PDE、量子局部端口或全部CC资格。

## 1. 复用、归一与新增模型输入

复用[1245 §§1—3](../1245/proof.md)的无焦散正标签流、[1247 §§1—3](../1247/proof.md)的全部内部变量及[1249 §§1—5](../1249/canonical.md)的齐次辛归一。守恒稠密流的Eulerian变分形式是成熟工具，参见[Brown (1993), §§2.1、5.1](https://arxiv.org/pdf/gr-qc/9304026)；这里新增的是1268指定杆读取与原内部钟/指针的具体来源连接，不将一般流体变分重报为发现。

记守恒向量密度为 $j^\mu$，proper数密度为 $n$：

$$
n=\frac{\sqrt{-g_{\mu\nu}j^\mu j^\nu}}{\sqrt{-g}},
\qquad U^\mu=\frac{j^\mu}{\sqrt{-g}\,n},\qquad U^\mu U_\mu=-1.
\tag{1}
$$

$n$不是新增外部控制。正的守恒标签测度可在无焦散图中给出该流；变分中的标签势和乘子保留原标签自由度，而非增加免费材料。取 $n>0$、$M_0>0$、$\gamma+2d>0$ 的明确有效域。保持原杆 $\phi^A$、$\mathcal B=\sum_A(\nabla\phi^A)^2$、$\gamma=1+p/5$ 及其它物质。

下文 $P=P_{\rm phys}$ 是指针的**物理动作强度**，与它配对的相位为 $\varphi$；不要与旧场 $p$、程序 $B$、记忆 $K$ 混同。定义

$$
C_0=3k^2=\frac3{10},\qquad
d(\theta,P)=\frac{\alpha}{C_0}g(\theta)P .
\tag{2}
$$

在1268的有限量子码中 $P_{\rm phys}=\varepsilon\sigma_z$，所以不存在第二个免费小耦合：$d=\varepsilon\alpha g\sigma_z/C_0$。经典作用本身先以连续内部正则图表示，不把这一图自动当成有限qubit的完整实现；有限码另见本轮指针分项。

保持窗、包括读取前置窗内，原两写入脉冲的残留应按原式保留；以下新交换推导先取旧脉冲关闭的保持子域，并在第5节说明其与旧开窗项的相加。原无新读取质量记为 $M_0$，满足 $\partial_JM_0=1$。在旧保持窗可写

$$
M_0=\frac{71}{100}+J+\frac{q^2+r^2+K^2}{2}.
$$

旧振子、程序、记忆、第二流体及其正基线均不删除。

## 2. 同一内部作用及精确变分

使用一个局部Eulerian标签图，令

$$
\mathcal A_\mu=\partial_\mu\ell+\beta_A\partial_\mu Z^A
 +J\partial_\mu\theta+P\partial_\mu\varphi
 +r\partial_\mu q+B\partial_\mu\psi+K\partial_\mu Y .
$$

新的控制流与杆的作用为

$$
\boxed{
S_{\rm c,\phi}=\int d^4x\,
\left[j^\mu\mathcal A_\mu
-\sqrt{-g}\left(nM_0+\frac{\gamma}{2}\mathcal B
+d(\theta,P)\mathcal B\right)\right].
}
\tag{3}
$$

其余旧作用照加。$\ell$保证 $\partial_\mu j^\mu=0$；$\beta_A,Z^A$给标签的沿流运输。内部辛项没有度规依赖。作用(3)等价于沿原守恒标签流使用

$$
M_{\rm per\ label}=M_0+\frac{d\mathcal B}{n},
\tag{4}
$$

但(4)显含密度、杆导数和度规。不能再套1245中“质量独立于 $n$、仅依赖场值”的世界线Legendre和应力模板。

对内部变量变分，记 $\dot f=U^\mu\nabla_\mu f$，旧保持子域得到

$$
\dot\theta=1,\quad \dot P=0,\quad
n\dot J=-\mathcal B\,d_\theta,\quad
n\dot\varphi=\mathcal B\,d_P.
\tag{5}
$$

所以控制不是外加坐标时脉冲，指针相位也不是无反作用的后设标签。旧自由振子仍满足 $\dot q=r,\dot r=-q$，旧记忆仍满足 $\dot Y=K,\dot K=0$；它们在 $\dot M_0$ 中抵消，给 $\dot M_0=\dot J$。对于未关闭的旧脉冲，(5)右边另加旧 $M_{0,\theta}$、$M_{0,P}$ 项及原内部方程；新的体积交换不替代旧场—流交换。

固定度规时 $\delta n/\delta j^\mu=-U_\mu/\sqrt{-g}$。因 $d$不含 $j^\mu$，流变分给

$$
\mathcal A_\mu=-M_0U_\mu .
\tag{6}
$$

机械化学势和惯性是 $M_0$，不是把(4)直接塞入旧尘流公式得到的 $M_0+d\mathcal B/n$。这是密度依赖项必须完整变分的具体区别。

## 3. 独立Hilbert来源与内部交换抵销

作度规变分时固定 $j^\mu$ 及内部标量，有

$$
\delta n=\frac n2(g_{\mu\nu}+U_\mu U_\nu)\delta g^{\mu\nu}.
$$

于是由同一作用直接得到

$$
T_{\rm c}^{\mu\nu}=nM_0U^\mu U^\nu,\qquad
T_{\rm rod}^{\mu\nu}=(\gamma+2d)S^{\mu\nu}
-\left(\frac\gamma2+d\right)\mathcal B g^{\mu\nu},
\quad
S^{\mu\nu}=\sum_A\nabla^\mu\phi^A\nabla^\nu\phi^A.
\tag{7}
$$

相对旧杆的新增来源是 $2dS^{\mu\nu}-d\mathcal B g^{\mu\nu}$。它已计一次交互能源和压力，不能再在流体来源重复加入 $d\mathcal B UU$。

杆方程与散度是

$$
\nabla_\mu[(\gamma+2d)\nabla^\mu\phi^A]=0,\qquad
\nabla_\mu T_{\rm rod}^{\mu\nu}
=-\frac{\mathcal B}{2}\nabla^\nu\gamma-\mathcal B\nabla^\nu d.
\tag{8}
$$

第一项仍由旧 $p$ 部门抵销。对第二项，利用同一内部方程和Cartan恒等式：

$$
\iota_U d\mathcal A=-dM_0-\frac{\mathcal B}{n}\,dd
\quad\text{（此式中的 $d$ 前缀为外微分）}.
$$

为避免符号混淆，把结果写成分量：

$$
n\bigl[M_0a^\nu+\dot M_0U^\nu\bigr]
=\mathcal B\nabla^\nu d,\qquad
\dot M_0=-\frac{\mathcal B}{n}\dot d,\qquad
M_0a^\nu=\frac{\mathcal B}{n}
(g^{\nu\mu}+U^\nu U^\mu)\nabla_\mu d .
\tag{9}
$$

具体推法是将 $\mathcal A=-M_0U^\flat$ 代入
$U^\mu(\partial_\mu\mathcal A_\nu-\partial_\nu\mathcal A_\mu)$：
左边为 $-\dot M_0U_\nu-M_0a_\nu-\partial_\nu M_0$；内部正则方程给右边
$-\partial_\nu M_0-(\mathcal B/n)\partial_\nu d$。使用守恒数流即得

$$
\boxed{\nabla_\mu T_{\rm c}^{\mu\nu}
=+\mathcal B\nabla^\nu d.}
\tag{10}
$$

(8)新增项与(10)严格相消。旧脉冲/场依赖开启时，Cartan右边还含旧外场力；它们按1247原完整来源相消，而新增 $(\mathcal B/n)\,dd$ 保持。这里没有由Einstein方程反定义来源，没有通过端点 $g=0$ 删除门中压力，也没有声称各部门分别守恒。

### 3.1 原杆的ADM动量必须同步更新

记 $\Gamma=\gamma+2d$。杆的ADM项改为

$$
\mathcal H_\phi
=\frac{\sum_A\Pi_A^2}{2\sqrt h\,\Gamma}
+\frac{\sqrt h\,\Gamma}{2}\mathcal B_{\rm sp},\qquad
\mathcal B=\mathcal B_{\rm sp}
-\frac{\sum_A\Pi_A^2}{h\Gamma^2}.
$$

因此在固定该图的独立变量时

$$
\partial_P\mathcal H_\phi=\sqrt h\,d_P\mathcal B,\qquad
\partial_\theta\mathcal H_\phi=\sqrt h\,d_\theta\mathcal B.
\tag{11}
$$

遗漏动量项的导数会错误地只得到空间梯度。下一节用同一标签变量核新正则资格；不是把(4)当作旧的无密度依赖质量。

### 3.2 无焦散标签图中的新物质括号

改回原标签变量 $X(z)$、$\theta(z),J(z),\varphi(z),P(z)$。内部辛项仍是同一个 $d\mu(z)$ 加权正则配对。把标签标量推前为

$$
d(x)=d\bigl(\theta(Z(x)),P(Z(x))\bigr),\qquad Z=X^{-1}.
$$

体积项在这一图中不含 $\dot X$，杆系数 $\Gamma$ 也不含几何速度。因此机械Legendre仍给

$$
E_c=\sqrt{h^{ij}(X)P_{{\rm mech},i}P_{{\rm mech},j}+M_0^2},
\quad H_c[N]=\int d\mu\,N(X)E_c.
$$

这里 $P_{\rm mech}$ 与内部 $P$ 不同，根式绝不可换为 $M_0+d\mathcal B/n$。场的Hamiltonian是原 $p,X$ 项加(11)前的杆项；总动量约束仍为原场动量项与
$-\int d\mu\,\xi^i(X)P_{{\rm mech},i}$，原其它流体保留。

固定标签内部值而变 $X$ 时，

$$
\delta_Xd(x)=-\partial_i d(x)\,\delta X^i(Z(x)).
\tag{11a}
$$

该式没有 $\partial\delta X$，也不引入额外的密度权变分，因为 $d$ 是**标量推前**，不是动作密度的推前。令
$\mathcal W(z)=\det(\partial X/\partial z)\sqrt{h(X)}/n_{\rm label}(z)$；
对 $d\mu$ 定义泛函导数，有

$$
\frac{\delta H_\phi[N]}{\delta X^i(z)}
=-N(X)\mathcal W\,\mathcal B(X)\partial_i d(X),\quad
\frac{\delta H_\phi[N]}{\delta\theta(z)}
=N(X)\mathcal W\,\mathcal B(X)d_\theta,\quad
\frac{\delta H_\phi[N]}{\delta P(z)}
=N(X)\mathcal W\,\mathcal B(X)d_P.
\tag{11b}
$$

其中 $\mathcal B$ 含杆动量部分。全部新增导数含 $N(X)$ 而无 $\nabla N$。所以新场—流机械交叉、内部交叉在交换 $N,M$ 后均是同点 $NM-MN=0$。在本候选 $d(\theta,P)$ 中，杆Hamiltonian不依赖它们的共轭 $J,\varphi$，也不依赖 $P_{\rm mech}$，其自身不另生内部/机械括号。

杆的场括号中，动能逆系数 $1/\Gamma$ 与梯度系数 $\Gamma$ 相乘为1。$\partial_i\Gamma$、原 $p$ 导数和其余无lapse导数项按1245的方法在反对称组合中抵消；剩下的lapse梯度正是原杆动量生成元。加旧场、流及其原交叉项后，在该正则图与边界条件下得到

$$
\boxed{\{H_{\rm mat}[N],H_{\rm mat}[M]\}
=D_{\rm mat}[h^{ij}(N\partial_jM-M\partial_jN)].}
\tag{11c}
$$

新材料不含 $\pi_h,\partial h,\dot h$，与原引力Hamiltonian的反对称交叉也仍是 $NM$ 型。故1245式(18)的同一待选引力括号结构可以接回这份新物质；若采用其CT条件选后的系数，约束保持结论按同样范围成立。**CT仍未采用。** 此处证明新的经典正则身份，未证明该材料所有非均匀初值的存在性、量子约束、规范固定后的任意量子完成或实际当地CC权限。

## 4. 均匀控制的动态主块

在1268同类局部共动jet取 $U=(1,0,0,0)$、$n,M_0>0$、
$\partial_i\phi^A=s\delta_i^A$、$\partial_t\phi^A=0$，
$\theta=t,\ P=P_0$；$\gamma,d$的背景变化作为低阶系数保留。
令 $\Gamma=\gamma+2d>0$，$D_*=\delta d=d_\theta\delta\theta+d_P\delta P$。

线性化的新主连接是

$$
\partial_t\delta\theta=0,\quad \partial_t\delta P=0,
$$

$$
\partial_t^2\varphi^A-\Delta\varphi^A
=\frac{2s}{\Gamma}\partial_A D_*,
\quad
\partial_t\delta v_i=\frac{\mathcal B_0}{nM_0}\partial_iD_*,
\quad
\partial_t\delta n=-n\,\partial_i\delta v_i .
\tag{12}
$$

最后两个内部量 $\delta J,\delta\varphi$ 由 $\delta\mathcal B,\delta n,\delta\theta,\delta P$ 驱动；它们不反向进入 $d(\theta,P)$。特别是 $\partial_JM_0=1$ 保证钟扰动没有新增密度依赖速率。背景 $\dot d\,\delta v$、原场/旧内部频率等为低阶项，不在(12)中冒充已消失。

用沿一个波矢的杆变量 $w=\partial_t\varphi_L$、$z=\partial_\parallel\varphi_L$，以及 $v,\delta n,\delta\theta,\delta P$ 做一阶主约化，$\partial_t{\bf u}=\mathsf A\partial_\parallel{\bf u}$：

$$
\mathsf A=
\begin{pmatrix}
0&1&0&0&2sd_\theta/\Gamma&2sd_P/\Gamma\\
1&0&0&0&0&0\\
0&0&0&0&bd_\theta&bd_P\\
0&0&-n&0&0&0\\
0&0&0&0&0&0\\
0&0&0&0&0&0
\end{pmatrix},\qquad b=\frac{\mathcal B_0}{nM_0}.
$$

故

$$
\det(\lambda I-\mathsf A)=\lambda^4(\lambda^2-1).
\tag{13}
$$

非零传播块仍为原杆锥；没有1268式(8)中非零常数项产生的负速度平方。这里**保留了动态流体、钟和指针**，不是再次冻结 $n$ 的杆测试。

但(13)不证明整个系统强双曲。零特征值一般有Jordan链，原无压流本就允许密度被速度梯度驱动；控制不均匀会再激发该链。对(12)的冻结常系数块、初始 $v=\delta n=\varphi=\dot\varphi=0$、固定Fourier波数 $q\ne0$，有精确强迫解的模界

$$
|\delta v(t)|\le b|q|t|D_*|,\qquad
|\delta n(t)|\le\frac{nbq^2t^2}{2}|D_*|,
$$

$$
|\varphi_L(t)|\le\frac{4s}{\Gamma|q|}|D_*|,
\qquad
|\partial_t\varphi_L(t)|\le\frac{2s}{\Gamma}|D_*|.
\tag{14}
$$

因此有限 $q_{\max},T$ 有具名多项式放大预算；它不是全部真实时变传播、量子尾或因果输出证书。一般非均匀控制背景、几何扰动的完整低阶耦合、焦散和任意高频健康性未在此解决。只因非零锥正确就宣布全部CC成立，与只因零速块有导数损失就否定原有限任务，都是超范围推断。

## 5. 同一齐次读取、钟余额和压力

沿原共动齐次族 $n=a^{-3}$、$\mathcal B=C_0a^{-2}$，统一强度图的辛前因子记为 $\mathcal A$：1249未缩放单位中为 $V$；当前1268物理相似族中为 $\mathcal A=V_0s^2$，$\varepsilon=\hbar/\mathcal A$。它不等于瞬时物理体积 $a^3V_0s^3$，也不是坐标重标。本文局部系数按同一归一有效单位书写；总内部动量为 $\Pi_\varphi=\mathcal A P$，详见[指针分项§2](pointer_proof.md)。同一作用给

$$
S_{\rm int}=-\mathcal A\int dt\,N\,\alpha g(\theta)P\,a,\qquad
S_{\rm ptr,symp}=\mathcal A\int dt\,P\dot\varphi ,
$$

$$
h_{\rm int}=\alpha gPa,\quad
\dot\varphi=\alpha ga,\quad
\dot P=0,\quad
\dot J\big|_{\rm new}=-\alpha g'Pa.
\tag{15}
$$

代入 $P=\varepsilon\sigma_z$，得到1268的
$\varepsilon\alpha g\,a\sigma_z$。原量子模型用 $F(Q)=\xi_Q e^Q$ 的有界完成；(15)只在该平台核内匹配。不能将全局完成区不加说明地当成同一局部杆作用。原七对自由演化、旧写窗末段和新读取同在一个生成元，而非冻结 $Q$。

新体积来源独立给

$$
\rho_{\rm int}=d\mathcal B=\frac{\alpha gP}{a^2},\qquad
p_{\rm int}=-\frac{d\mathcal B}{3},\qquad
\dot{(nM_0)}+3HnM_0=-\mathcal B\dot d.
\tag{16}
$$

因此新总来源满足连续方程，且新增门不是“零压力、只改相位”。在任意满足 $a\le3,|\dot a|\le4$ 的齐次读窗内，取原1268非负归一脉冲、$g(t_a)=0$，积分(15)得到新门对钟方程的累计贡献

$$
\Delta J_{\rm new}(t)
=-\alpha P\,g(t)a(t)
+\alpha P\int_{t_a}^t g(u)\dot a(u)\,du ,
$$

$$
|\Delta J_{\rm new}|\le|\alpha P|(3g_{\max}+4),\qquad
|\Delta J_{\rm new}(t_b)|\le4|\alpha P|.
\tag{17}
$$

旧脉冲的钟交换仍另计，不能把此附加贡献当全部 $J(t)-J(t_a)$。式(17)按实际 $a$ 取值，不预设新增门后的几何等于旧轨道。

例如采用1268的 $\alpha=1/16$、$\tau=10^{-8}$、
$g_{\max}=218750000$，则前式系数为 $164062501/4$。
若 $|P|\le\varepsilon\le2^{-100}$，新贡献小于 $2^{-74}$；
杆系数增量
$2|d|\le(273437500/3)\varepsilon<2^{-73}$。
这些是给定有限管道的资源充分界；原极小 $\varepsilon$满足它们，但不能仅凭它们为所有量子支撑签严格正 $J$ 或为新全空间解签全窗存在。1268已有的有界来源/波包误差保持其原范围。

## 6. 核验与本轮停止线

[material_check.py](material_check.py)一个作者组，精确有理核验：

- 在有非零流速、非对角度规变分及杆时间梯度的有限jet上，自动微分同一个作用密度，逐分量核(7)；并核(9)投影/收支。
- 核(11)的动量贡献和(11c)中新交叉的局部 $NM$ 抵销；从矩阵本身计算特征多项式(13)，展示零速Jordan链而不以有限样本证明全PDE。
- 核齐次源连续式、原动作归一、脉冲积分和(17)有限预算。代码只作解析身份的有限复算，不运行新PDE或假称已实现当地量子仪器。

结果首次 --write-results 排他创建，默认只读复算匹配 [material_results.json](material_results.json)。

**准入收益：**1268尚未闭合的 $-\mathcal B\nabla d$ 现在由同一原内部流明确承担；同齐次读取有了具名的控制/压力字典。**停止线：**尚未证明全部空间量子码、实体指针/环境/输出来源、一般控制非均匀误差和真实合法输出；也未以本轮替代独立场方程验收。有限码分项与本材料共享归一和作用，额外权限仍各自登记。新增模型不是认知公理唯一推论，CT仍未采用；不为修复所有PDE或制造无限设备扩张本轮。
