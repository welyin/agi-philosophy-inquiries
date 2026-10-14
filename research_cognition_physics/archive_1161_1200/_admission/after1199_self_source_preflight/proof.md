# 完整约束的混合来源与原物质钟：预检推导

2026-10-14。未编号方向预检；不把给定ADM动力学称为从认知公理生成GR。本稿只证明所列混合阶系数，以及所选保留Hamiltonian中原钟概率的初始二阶导数。

## 1. 输入与边界

复用条件性ADM成果，选单位横截面的 $T^2\times[0,1]$，场仅依赖 $z$，最大切片 $\pi=0$，边界 $\psi=N=1$。采用固定平面参考的有限边界Hamiltonian，**不是渐近平直ADM质量**。此处单位吸收引力常数；后文的 $\epsilon$ 是所声明弱耦合模型的输入，不是已推导的现实 $G$。

取偏振度规和形状动量

$$
q=\psi^4\operatorname{diag}(e^f,e^{-f},1),\qquad P=P_f.
$$

记 $\Pi_i=q_{ii}\pi^{ii}$、$C=\Pi_z$，则 $\Pi_x-\Pi_y=P$、$\Pi_x+\Pi_y+C=0$，直接给出

$$
\pi^{ij}\pi_{ij}=\frac{P^2}{2}+\frac{3C^2}{2},\qquad
\mathcal D_z=-2C'+Pf',\qquad C'=\frac12Pf'.
$$

因此 $C$ 从TT二阶开始，$C^2$ 是纯TT四阶。它不贡献静止物质密度与TT二阶的混合项，但不能在声称完整四阶能量时丢弃。两端固定 $C=0$ 须满足 $\int Pf'=0$ 或另计边界纵向动量。

物质初片动量取零，$\rho$ 表示每坐标体积的静止质量能源；在此变分下独立于 $\psi$。**初片静止不等于完整历史静止**；维持固定剖面的支撑、压力和控制若存在，后续必须进入来源和资源账。

## 2. 由同一约束求来源与lapse

令 $A=P^2/2$、$B=f'^2/2$、$G=(-\partial_z^2)^{-1}$ 为Dirichlet正对称Green算符。共形底度规的曲率精确为 $\bar R=-f'^2/2$。保留所需阶数，Hamiltonian约束为

$$
\mathcal C=\psi^{-6}A+\psi^2B+8\psi\psi''+\rho+\cdots=0.
$$

令 $\psi=1+u$，得到

$$
u=\frac18G(A+B+\rho)+\cdots.
$$

保留边界项并使用约束，有限边界能量 $-8[\psi\psi']_0^1=-8[\psi']_0^1$ 给出

$$
H_{\rm red}=\int(\psi^{-7}A+\psi B+\rho/\psi)+\cdots
=\int(A+B+\rho)+\frac18\langle A+B+\rho,G(-7A+B-\rho)\rangle+\cdots.
$$

利用 $G$ 的对称性，严格选出混合系数：

$$
H_{\rho,\mathrm{TT}}=-\langle\rho,GA\rangle
=-\frac12\langle\rho,GP^2\rangle,\qquad
H_{\rho,\rho}=-\frac18\langle\rho,G\rho\rangle.
$$

梯度交叉项 $\langle\rho,GB\rangle$ 的系数抵消；这不是另手填入的平均Newton势。

独立用最大切片保持方程，记 $N=1+n$，对含边界项的Hamiltonian作共形变分，线性所需项为

$$
-6A+2B+16u''+8n''=0.
$$

代入约束及相同Dirichlet条件：

$$
n=-GA-\frac14G\rho
=\frac{\delta H_{\rm red}}{\delta\rho}-1.
$$

该相等检验了来源读数和时间率来自同一混合阶规则。纯TT四阶、有限物质运动和量子余项尚未统一处理。

## 3. 采用原941钟而非新增读取

[941](../../../archive_935_955/research_note_941.md)继承[925原模型](../../../archive_923_934/925/record_resource_exchange.py)。其 $q=0$、资源标签3的闭合二态块为

$$
M=M_0 I+g\sigma_x,\qquad M_0=526/5,\quad g=2/5.
$$

初态仍是原 $|A,0,3\rangle$，效果仍是原接收投影 $\Pi_B=|B,0,3\rangle\langle B,0,3|$。没有新增Y方向门或Ramsey读取。当前只用已知 $q=0$ 初态，未签完整未知qutrit与参考 $R$ 的联合通道。

在上述混合系数基础上明确选择保留模型

$$
H_{\rm ret}=M+H_{\rm geo}-\frac\epsilon2 M X-\frac{\epsilon b_w}{8}M^2,
\quad \epsilon=10^{-3},\quad \rho=Mw,\quad b_w=\langle w,Gw\rangle,
$$

其中 $X=\langle w,GP^2\rangle$，$M$ 与几何算符对易。$\epsilon$ 是显式弱耦合输入；**从完整量子ADM到此有限模式模型的误差还没有证明**。所有仅作用于几何的Hamiltonian项均可放入 $H_{\rm geo}$；不能因此省去带物质的同阶项。

令 $\alpha=1-\epsilon b_wM_0/4$，钟翻转矩阵元为 $g(\alpha-\epsilon X/2)$。因为初态在 $\Pi_B$ 核内，严格有

$$
p_B''(0)=2\|\Pi_BH_{\rm ret}|A,\chi\rangle\|^2
=2g^2\langle\chi|(\alpha-\epsilon X/2)^2|\chi\rangle.
$$

故几何专属项不影响这个初始二阶导数，但会影响有限时间预测。此等式不把 $X$ 换成平均值：方差和相干项均保留。

## 4. 两模见证和精确有理界

取 $w=1$，$Gw=z(1-z)/2$，故 $b_w=1/12$。取两个站波 $r_n=\sqrt2\sin(n\pi z)$，$n=1,3$，自由频率 $\omega_1=\pi,\omega_3=3\pi$，$P(z)=P_1r_1+P_3r_3$。

$$
X=P^TWP,\qquad
W=\begin{pmatrix}
1/12+1/(4\pi^2)&-3/(16\pi^2)\\
-3/(16\pi^2)&1/12+1/(36\pi^2)
\end{pmatrix}.
$$

以 $|\chi_\pm\rangle=(|1,0\rangle\pm|0,1\rangle)/\sqrt2$ 比较。它们在 $H_{\rm geo}=\sum_{n=1,3}\omega_n(a_n^\dagger a_n+1/2)$ 下有相同的自由能谱概率分布，均值均为 $4\pi$，并非能量本征态；完整相互作用能源不同，不宣称等完整能量。以下 $\Delta$ 表示正号态减负号态，梯算符直接给

$$
\Delta\langle X\rangle=-\frac{3\sqrt3}{8\pi},\qquad
\Delta\langle X^2\rangle=(\pi+\pi^{-1})\Delta\langle X\rangle.
$$

所以

$$
\Delta p_B''(0)
=2g^2\epsilon(-\Delta\langle X\rangle)
\left[\alpha-\frac\epsilon4(\pi+\pi^{-1})\right]
>\frac{14236293}{220000000000}>64\times10^{-6}.
$$

最后一步仅用 $3<\pi<22/7$、$\sqrt3>17/10$；有理余量为 $156293/220000000000$。浮点诊断约 $6.59572416\times10^{-5}$，不作为严格下界来源。

$W$ 为正且 $\|W\|\leq\operatorname{tr}W<16/81$。对原二态钟 $M\leq528/5$，保留二次Hamiltonian的完整Fock动量矩阵满足

$$
I-\epsilon MW>\frac{16523}{16875}I>0,
\qquad 1-\epsilon b_wM/8\geq9989/10000>0.
$$

这只保证采用上述自由二次 $H_{\rm geo}$ 时保留模型的正性；任意附加几何项不改变初始钟导数，并不意味着它们都保持正性。高占据态也不保证处在弱场有效域。代码中每模5级足以精确容纳 $X|\chi_\pm\rangle$ 所需总占据数至多3的全部中间态；有限矩阵用于这些多项式矩的交叉核验，未用于声称有限时无截断演化收敛。

## 5. 什么尚未交付

初始二阶导数非零只说明此读口对候选自源项并非不可见。平方给出的 $\epsilon^2$ 项在所选保留模型内精确，不表示完整ADM同阶遗漏已受控。要交付有限读数差，仍需给预定有限读时、Taylor或完整演化余项、物质/参考保真误差，并把运动、支撑、边界、被舍模式的贡献接回同一过程。

因此本预检不计第1200轮，不算完整GR、实际有限引力预测或暗部门。其价值是保留这一具体接口继续做有限过程验证，同时停止泛泛重做已知自耦定理。

成熟ADM背景参见[Schäfer与Jaranowski的Hamiltonian综述](https://link.springer.com/article/10.1007/s41114-018-0016-5)。本文系数按以上明确约定自行展开；综述不是本项目实际认知操作实现的证明。
