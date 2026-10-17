# 1265：原紧支静态共形背景的完整来源约束与有限残差

2026-10-16。已预约1265，author_groups=1；新增认知公理0，不改变已验收1264或当前目标。最终验收以本轮acceptance记录为准。本稿只检验一个明确接法：把1264事前给定的静态共形背景，直接要求为同一最小KG材料及其它非负法向来源产生的Einstein几何。结论是**该精确接法不相容**，不是整个认知研究纲领、任意量子引力或有限有效GR的反证。

[作者复算](check.py)默认只读比对[保存结果](results.json)，仅首次 --write-results 使用排他创建。程序用精确多项式系数核验完整16分量Einstein张量，再核同profile的有限有理预算；连续场域、来源正性及分部积分由下文解析证明承担。

## 1. 旧接口与适用前提

- [1194 §2](../../archive_1161_1200/1194/proof.md)已经区分总Hamiltonian正谱与局部来源正性；本稿不从前者推出任意量子场的能量条件。
- [1235](../1235/proof.md)、[1236 §3—4](../1236/proof.md)已有同对象有限约束与完整ADM来源；[1245 §5—6](../1245/proof.md)、[1247](../1247/proof.md)、[1252](../1252/proof.md)已有必须同步选择合法初始几何、不能只贴新来源的实例。那些动态非紧支同类初片没有被本稿否定。
- [1264](../1264/source_proof.md)只签给定几何中的来源—表示—记录连接，没有声称背景支持来源已经实现。因此本稿是后继准入，不是撤销1264。
- 成熟约束及边界定义直接复用[Gourgoulhon，§6.5.2、§7及§8.2.5](https://arxiv.org/pdf/gr-qc/0703035)，不重报正质量定理。下面在原profile上直接计算有限残差。

采用签名 $(-,+,+,+)$、$\kappa=8\pi G>0$，并仅在条件Einstein接口
$G_{\mu\nu}=\kappa\langle T_{\mu\nu}\rangle$
下判断。CT仍未采用。外部精确平直、来源有限且渐近消失时，有效宇宙学常数须为0；若另外保真空常数，先与$\Lambda$作同一账目，不能保平直远区却遗漏它。

原度规、profile和空间尺度保持为

$$
g_{\mu\nu}=e^{2\omega(x)}\eta_{\mu\nu},\qquad
\omega=\lambda f_R,\qquad f_R=F(x/R),\qquad \lambda>0,
$$

其中1264的径向 $F$ 在单位球为1，在半径2以外为0，光滑且介于0和1。

## 2. 全Einstein张量，不只比较迹

共形Christoffel为

$$
\Gamma^\alpha{}_{\mu\nu}
=\delta^\alpha_\mu\partial_\nu\omega
+\delta^\alpha_\nu\partial_\mu\omega
-\eta_{\mu\nu}\eta^{\alpha\beta}\partial_\beta\omega .
$$

代入Ricci定义，四维中得

$$
R_{\mu\nu}=-2\partial_\mu\partial_\nu\omega
-\eta_{\mu\nu}\Box\omega
+2\partial_\mu\omega\partial_\nu\omega
-2\eta_{\mu\nu}(\partial\omega)^2,
$$

$$
R=-6e^{-2\omega}\left[\Box\omega+(\partial\omega)^2\right].
$$

因此静态情形的完整协变分量为

$$
\boxed{
G_{00}=-2\Delta\omega-|\nabla\omega|^2,\qquad
G_{0i}=0,
}
\tag{1}
$$

$$
\boxed{
G_{ij}=-2\partial_i\partial_j\omega
+2\partial_i\omega\partial_j\omega
+\delta_{ij}\left(2\Delta\omega+|\nabla\omega|^2\right).
}
\tag{2}
$$

这里导数与 $\Delta$ 属于该原空间图的Euclidean坐标。法向
$n=e^{-\omega}\partial_t$、lapse $N=e^\omega$、$h_{ij}=e^{2\omega}\delta_{ij}$，shift为0且 $K_{ij}=0$。独立的三维标量曲率检查是

$$
{}^{(3)}R=e^{-2\omega}
\left[-4\Delta\omega-2|\nabla\omega|^2\right],
\qquad G(n,n)=\tfrac12{}^{(3)}R.
\tag{3}
$$

所以式(1)正是该初片Hamiltonian约束的分量，不是仅迹关系。
完整空间来源还必须匹配式(2)：其无迹部分为
$-2(\partial_i\partial_j\omega)^{\rm TF}
+2(\partial_i\omega\partial_j\omega)^{\rm TF}$，
一般不是各向同性压力。

记法向能源 $\rho$ 和正交空间压力迹 $S$，若满足Einstein方程，必须有

$$
\kappa\rho=e^{-2\omega}(-2\Delta\omega-|\nabla\omega|^2),
\qquad
\kappa S=e^{-2\omega}(4\Delta\omega+5|\nabla\omega|^2).
$$

在线性阶 $S=-2\rho$；这不能自动冒充低速无压力材料。该关系本身不排除一切正能标量态；真正无需另外状态估计的阻碍是下一节的法向来源正性。

## 3. 本次来源正性为何适用，哪里不适用

对于1264采用的最小复杂KG一粒子部门，固定同一静态真空相对正常序。纯一粒子态的协变波函数记q，它与L²位置振幅不同。直接从同一最小作用得到

$$
\langle:T_{00}(x):\rangle
=\|\partial_tq(x)\|^2+\|\nabla q(x)\|^2
+e^{2\omega(x)}q(x)^\dagger M^2q(x)\ge0.
\tag{4}
$$

这里保原完整 $M$，没有用平均质量替代。质量部门的来源块对角不等于测量或去相干；任意原输入、旁挂参考取部分迹以及正常混合，都保持式(4)的非负性。有限能源时可按非负局部可积密度或涂抹形式理解，避免要求点值处处存在。任意额外材料若也有非负 $T(n,n)$，只会加强下文阻碍。

**范围不能省略：** 式(4)是这份最小、正质量、一粒子相对来源的直接身份，不是任意重整化量子应力的能量条件。曲背景真空极化、负能量准备、改进耦合或其它新增部门可能贡献负的局部来源；1264没有给出这些部门的完整有限证书。若引入它们，必须重新公开作用、状态、重整化与来源，不能把它们默默视为已经抵消了式(1)。

本稿也不把正经典/半经典来源方程偷换成完整量子几何的算符约束。

## 4. 直接积分障碍与零ADM边界

取任意非负紧支 $\chi$，在 $\operatorname{supp}f_R$ 上等于1。则

$$
\int \chi G_{00}\,d^3x
=-\int|\nabla\omega|^2d^3x<0.
\tag{5}
$$

散度项的边界为0。任何满足式(4)的来源却有
$\int\chi\kappa\langle T_{00}\rangle d^3x\ge0$，
故该精确静态背景不能满足完整Einstein方程，甚至不能满足初片法向约束。增加其它非负来源不能修复。

式(5)是明确的正权测试，**不是声称坐标积分 $G_{00}$ 就是ADM能源或物理总能源**。例如 $\int G(n,n)dV_h$ 有不同权重，不能混用。

边界的另一种说明是：$h_{ij}=\delta_{ij}$ 在有限球以外严格成立，所以ADM表面积分为0。非零profile已经迫使某处法向来源为负；无需调用一般正质量定理再证明一次。该计算使用完整 $\mathbb R^3$、无内边界且光滑的原1264空间；若换成有洞、边界通量或不同拓扑，应重列任务而非继承此结论。

## 5. 一个线性于幅度的有限同profile残差

仅用式(5)得到二阶幅度。为避免弱场下一阶被积分抵消，使用与原profile完全同源的非负权：

$$
\chi=f_{2R},\qquad
w=\chi(1-f_R),\qquad
Z=\int w\,d^3x>0.
$$

$\chi=1$ 覆盖 $f_R$ 支撑，$w$ 的支撑在半径 $4R$ 内。
分部积分给

$$
\int w\Delta f_R\,d^3x=\int|\nabla f_R|^2d^3x=:I_R>0.
$$

于是完整加权约束残差满足

$$
\begin{aligned}
\mathcal C_w
&=\frac1Z\int w\left(G_{00}-\kappa\langle T_{00}\rangle\right)d^3x\\
&=-\frac{2\lambda I_R+
\lambda^2\int(1-f_R)|\nabla f_R|^2d^3x+
\kappa\int w\langle T_{00}\rangle d^3x}{Z}\\
&\le-\frac{2\lambda I_R}{Z}<0 .
\end{aligned}
\tag{6}
$$

径向profile从 $r=R$ 的1降到 $r=2R$ 的0。Cauchy–Schwarz直接给

$$
1=\left|\int_R^{2R}f_R'(r)dr\right|^2
\le\left(\int_R^{2R}r^2f_R'(r)^2dr\right)
\left(\int_R^{2R}\frac{dr}{r^2}\right),
$$

所以

$$
I_R\ge8\pi R,\qquad
Z\le\frac{256\pi R^3}{3},\qquad
\boxed{\mathcal C_w\le-\frac{3\lambda}{16R^2}.}
\tag{7}
$$

这是该固定任务上有量纲 $({\rm length})^{-2}$ 的有限加权曲率—来源残差。可等价写成

$$
\mathcal C_w=\frac1Z\int
w e^{-\omega}\left[G(n,n)-\kappa\langle T(n,n)\rangle\right]dV_h,
$$

所以是声明法向/图下的正权约束测试，而不是对某个不定号坐标分量任意取平均。

1264已经固定
$r=10^{-4}$、$\lambda=\log(1+r)>r/(1+r)$、
$R=38201213950$。因此严格非零下界可只用有理算术表示：

$$
|\mathcal C_w|>
\frac{3r}{16R^2(1+r)}
=\frac3{233516588884543995193640000}
>1.2847\times10^{-26}.
\tag{8}
$$

该有理值由[作者程序](check.py)直接读取1264的既有r、R参数后精确重算；没有积分新的PDE、扫描函数或修改1264结果。16分量核验采用由Christoffel到Ricci再到Einstein的独立系数展开，不以有限抽样冒充全角恒等式。

**有限精度解释：** 若某候选要求这份归一约束测试误差小于式(8)，当前“原背景由非负来源精确维持”的接法不能通过。若实际容差更大，这个下界本身不足以拒绝近似有效模型。不得把一个极小残差直接当作整个阶段统一失败，也不得在未指定误差时宣布近似问题已解决。

若试图保该精确背景而增加支持者，它至少需要负的该加权 $T_{00}$ 来抵消式(6)，同时还需匹配空间应力和动量；仅有负迹或只调温度/记录参数不够。该说明不是要求所有认知模型必须有负能量。

## 6. 后继应释放什么，哪些单项不够

1. **只放开lapse不够。** 保同一静态空间 $h=e^{2\omega}\delta$、shift0和 $K_{ij}=0$ 时，$G(n,n)={}^{(3)}R/2$ 与lapse无关。把N改成另一个函数可以改变钟与空间演化方程，却不能修复上述空间初片约束。

2. **可以放开非紧支空间尾。** 对时间对称共形空间写 $h=\Psi^4\delta$，Hamiltonian约束为
   $$
   \Delta\Psi=-2\pi G\,\rho\,\Psi^5 .
   $$
   非负局部来源一般带出非零边界通量和 $1/r$ 尾；不能要求一个非零正来源的几何仍在有限球外严格恢复平直。下一候选应由合法约束选择该尾，不是给1264平顶补一个任意修正函数。

3. **lapse与空间势应允许独立。** 普通弱静态正质量分支有 $N\simeq1-U$、$h_{ij}\simeq(1+2U)\delta_{ij}$，而1264把两者锁成同一个 $\omega$。这只是成熟弱场例子，说明应放开锁定；不是本稿已证明原量子来源产生该分支。实际压力、控制与原记录仍需同一作用验证。

4. **另一合法选择是非静态 $K_{ij}\ne0$。** 约束会加入 $K^2-K_{ij}K^{ij}$，已有1245—1252正动态基片可复用。但必须同时满足动量约束和演化；本稿不通过任意指定K来宣布新模型成功。

不必同时把所有选择做完。本次选择力在于：**停止把1264紧支、全时空同共形、静态的测试背景继续修成正来源自引力解**。下一项应从合法的完整约束/边界族选择同一原材料的有限任务，再运输1264的表示及来源工具。不要继续调f的导数、包宽或测量器件来绕过一个已经明确的约束不相容。

## 7. 本次地位

本轮是一组限定候选的来源—几何相容性验收；不宣称新的一般正质量定理，也没有新增认知公理或GR完成声明。1264在给定背景中的来源—记录共同结果保持验收。这里的负向结论只针对明确的精确Einstein接法及已声明非负一粒子来源；量子真空部门、非静态几何、不同边界和有限精度任务均按各自条件另判。


后继范围由本轮[NEXT](NEXT.md)登记。只有归一曲率—来源误差预算低于式(8)时，本有限下界才删除相应近似候选；1264的记录概率误差不具有相同单位，不能直接作该预算。
