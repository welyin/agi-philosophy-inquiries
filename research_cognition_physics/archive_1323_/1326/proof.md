# 1326：同一动态参考的低带、来源与反作用

2026-10-19。作者证明稿。[准入](../_admission/round1326_admission.md)。新公理0；阶段目标保持原范围。

## 1. 问题及完整有限父族

1325的认知操作含义是：交出或接回一份职责，通过同一参考和同一接收端保资料与伙伴。本轮不再添加一句认知口号，而是检验已经选出的作用规则能否在参考有自身动力学时，仍给可核的有效能源和来源。自由动力学及物理字典是公开建模输入；一般低能消元直接复用1316和[Schrieffer–Wolff方法](https://arxiv.org/abs/1105.0675)，源随字典变化的要求复用1322。

令 $\mathcal R_B=\mathrm{span}\{|n_0,n_1\rangle:n_0+n_1\le B\}$，所有升降算符均为这个空间上的压缩，最高边界保留。物质角色为D、L_0、L_1、U及有限味重数；只讨论单费米部门。沿1325定义

$$
K=N_b+n_L+2n_U,\quad
V_B=P_B K^{-1/2}(X+X^\dagger)P_B,
$$

$$
X=\sum_{\alpha ij}l_{\alpha i}^\dagger(M_D)_{ij}d_j a_\alpha
+\sum_{\alpha\beta ij}l_{\alpha i}^\dagger(M_U)_{ij}u_j
\epsilon_{\alpha\beta}a_\beta^\dagger.
\tag{1}
$$

零K逆为零。等价地，从D输入的系数是 $a_\alpha N_b^{-1/2}$，从U输入是 $\epsilon_{\alpha\beta}a_\beta^\dagger(N_b+2)^{-1/2}$；不能把两个分母混同。正反向均来自同一Hermitian矩阵。压缩给

$$
\|V_B\|\le g=\max(\|M_D\|,\|M_U\|).
\tag{2}
$$

同一参数化父族可写为

$$
H=I\otimes H_R(f)+D(p)\otimes I+\sum_s E_s n_s+V_B,
\quad H_R(f)=\omega N_b-\sum_\alpha(f_\alpha a_\alpha^\dagger+f_\alpha^*a_\alpha).
\tag{3}
$$

- 旧设置：f=0、p=0、$E_L-E_D=\omega=E_U-E_L$，在安全K窗口内完整承接1325共振交付。
- 新设置：本轮取各E_s=0，$f=(f,0)$实，使用下述有隙低带；原K不再守恒，不能继续沿用K截断空间。

主例B=66使旧K≤66全块均在同一有限父空间内；打开驱动后，N_b≤66是**新明确定义的精确有限模型**，不称为无截断理论的已证近似。两种设置共享作用及承载者，不同时满足同一个实验合同。

f、ω、门户、动量核、控制切换和读口权限是公开输入。f作为已登记有效控制来源计账；没有构造全部内部自主控制器。制备两个f点的状态也不是已经实现绝热运输。若另给f(t)，驱动功须按 $\int\mathrm{Tr}(\rho\partial_fH)\,df$ 计算，不能省略。

驱动一般也改变被驱动子系统的数目与原荷。其交换由同一 $i[H,K]$ 或 $i[H,Q_{\rm charge}]$ 记到控制来源；静态H(f)守恒不等于驱动承载者已经量子化或全认知整体已自治闭合。这里沿已登记有限控制权限工作，不授予额外无代价的荷/相位供应。

## 2. 参考基态与完整低带

以下谱定理适用于任何具有唯一基态χ、基能e_0及隙Δ的有限H_R。设物质维d、$\|D(p)\|\le p_*$，定义

$$
E\psi=\psi\otimes\chi,\quad P=EE^\dagger,\quad Q=I-P,
\quad h_*=p_*+g,\quad \delta=\Delta-2h_*>0.
\tag{4}
$$

新设置按P/Q分块，记 $A=E^\dagger HE$、$B_0=QHE$、$C=QHQ|_{\mathrm{ran}Q}$，有

$$
A=e_0I+D+E^\dagger V_BE,\quad b=\|B_0\|\le g,
\quad C\succeq(e_0+\Delta-h_*)I.
\tag{5}
$$

**命题1（全部有序谱）。** H恰有d个隔离低态。若a_j为A升序本征值，则

$$
a_j-\frac{b^2}{\delta}\le\lambda_j(H)\le a_j\quad(1\le j\le d),
\qquad\lambda_{d+1}(H)-\lambda_d(H)\ge\delta.
\tag{6}
$$

证明：对未扰动的 $I\otimes H_R$ 使用Weyl，得 $\lambda_d(H)\le e_0+h_*$、$\lambda_{d+1}(H)\ge e_0+\Delta-h_*$。对任意低谱参数λ有C−λ≥δ，故

$$
0\preceq B_0^\dagger(C-\lambda)^{-1}B_0\preceq(b^2/\delta)I\preceq(g^2/\delta)I.
$$

在C−λ正定时，H−λ的负惯性恰等于 $A-\lambda-B_0^\dagger(C-\lambda)^{-1}B_0$ 的负惯性。由该计数及压缩min–max逐个得到(6)，包括简并，不仅是最高正支。余项会产生标量、动能及非质量修正，全部保留；不把谱界改称任意坐标下的矩阵范数界。

对于实轴f，本例H_R按n_1=j分为一维有限Jacobi块：$\omega j+\omega n_0-f(a_0+a_0^\dagger)$，0≤n_0≤B−j。它们的全体谱确定实际Δ；不使用无限Poisson基态冒充有限基态。在本轮主例f>0中，j=0块的唯一正系数基态给χ，较小截断的基能按min–max不低于j=0块，故j≥1的基能至少高ω。j=0第一激发仍需一并核。

## 3. 任意未知态与伙伴的低带字典

令Π为H的整个d维低带谱投影。RanΠ不含纯Q非零向量：否则其期望能源既≤e_0+h_*，又≥e_0+Δ−h_*，矛盾。因此 $E^\dagger\Pi E>0$，定义

$$
J=\Pi E(E^\dagger\Pi E)^{-1/2},\quad
h_{\rm dress}=J^\dagger HJ.
\tag{7}
$$

这是固定canonical选择，满足 $J^\dagger J=I$、$HJ=Jh_{\rm dress}$ 及

$$
e^{-itH}J=Je^{-it h_{\rm dress}}.
\tag{8}
$$

张量任意伙伴恒等算符后仍成立，完整低带信息和关联保持。J依赖参考及完整动力学；不得把它当免费制备权限。

如需从实际产品准备E出发，取低带正交基U，Λ为其谱。写X=E†U、Y=QU，有 $CY-Y\Lambda=-B_0X$。C和Λ有有序谱间隔δ，Sylvester半群积分给 $\|Y\|\le b/\delta=\rho$。在另声明ρ<1的域内，主角关系给

$$
\|J-E\|\le\zeta=\sqrt{2-2\sqrt{1-\rho^2}},
\quad
\|e^{-itH}E-Ee^{-it h_{\rm dress}}\|\le2\zeta.
\tag{9}
$$

因此任意伙伴的半迹距离也受min(1,2ζ)控制。右侧保留精确h_dress，并非把第一阶A当完整生成元。准备、读出、记录及容差均须登记；小误差不意味着参考未受影响。

## 4. 真正源矩阵及完整能源

固定门户、B、p和ω后改变实f，完整父源为

$$
S_f=\partial_fH=-I\otimes(a_0+a_0^\dagger).
\tag{10}
$$

若低带随f光滑，设 $\mathcal A_f=J^\dagger\partial_fJ$，由HJ=Jh_dress及J†J=I直接微分得

$$
J^\dagger S_fJ=\partial_fh_{\rm dress}-[h_{\rm dress},\mathcal A_f].
\tag{11}
$$

连接项对一般未知低带输入不能省略。简单能级有Hellmann–Feynman；简并处使用整块来源矩阵，不随意挑矢量定义导数。改变低带坐标时全部矩阵共同运输。

完整本征态ψ_j的相对参考激发能源满足

$$
\lambda_j-e_0=\langle D\otimes I\rangle_j+\langle V_B\rangle_j
+\big(\langle I\otimes H_R\rangle_j-e_0\big),
\quad \langle I\otimes H_R\rangle_j-e_0\ge0.
\tag{12}
$$

比较激发来源时扣e_0′，总账仍保e_0与e_0′。这是真实相互作用引起的参考能源，不是从近似产品态宣布零消耗。若需绝对正能源，在声明有限f域加入一个固定常数即可，不添加随f变化而未记账的抵消项。

## 5. 实际有限读口及它的代价

在有限参考中 $\|S_f\|\le s=2\sqrt B$。给定已登记的源测量权限，取

$$
F_\pm=\frac12(I\pm S_f/s),\quad L_\pm=\sqrt{F_\pm},
\quad W_{\rm read}=|+\rangle L_++|-\rangle L_-.
\tag{13}
$$

$W_{\rm read}^\dagger W_{\rm read}=I$，记录差给 $s(p_+-p_-)=\langle S_f\rangle$；完整分支、未读记录相干和伙伴均保留。忽略记录的系统通道为ΣL_±ρL_±；其能源变化为

$$
\Delta E_{\rm read}=\mathrm{Tr}\!\left[H\left(\sum_\pm L_\pm\rho L_\pm-\rho\right)\right].
\tag{14}
$$

这份变化须由测量控制/记录装置承担；不能把读源当免费被动过程。纯blank记录、允许的有限测量及其能源边界都是公开输入。有限等距可膨胀为酉，但这不证明全部物理仪器或自主控制已由认知原则生成。本轮只认证这个明确合同及其完整有效账。

## 6. 有限有效物理字典与留出检验

继承Dirac自由核，在单味、固定动量螺旋度块取(D,L_0,L_1,U)次序

$$
D(p)=\operatorname{diag}(p,-p,-p,p).
$$

此空间/手性字典是输入，角色名称不能代替它。轴向参考下第一阶 $A-e_0I=D+E^\dagger V_BE$ 的正负支具有 $\pm\sqrt{p^2+|m_D(f)|^2}$ 和 $\pm\sqrt{p^2+|m_U(f)|^2}$ 形式；m由实际χ_B与同一V求得。完整父谱及全部修正另由(6)—(12)控制。故只得到明确有限Dirac能量核的字典，不签完整场真空、全部洛伦兹过程、现实质量极点或Higgs场。

事前固定ω=1、B=66、p=0.005，归一门户M_D=0.013、M_U=0.011。普通线性对照删除K归一化，在同一有限参考上只于f=1匹配第一阶两门户，随后保持系数不变。裸对照范数可取 $g_{\rm bare}=\sqrt{B+2}\max(|M_D^{\rm bare}|,|M_U^{\rm bare}|)$；也可在完整K块求更紧且事前确定的范数。校准不要求两模型在f=1的完整父谱精确相同。

令 $a_j^{\rm exc}(f)=\lambda_j(A)-e_0(f)$、$\epsilon(f)=b(f)^2/\delta(f)$。b由实际完整跨块计算，不是另设的拟合量；全局g仍用于谱隙分离。实际同一有序低支有限变化属于

$$
\big[\,a_j^{\rm exc}(2)-a_j^{\rm exc}(1)-\epsilon(2),\
a_j^{\rm exc}(2)-a_j^{\rm exc}(1)+\epsilon(1)\,\big].
\tag{15}
$$

无需微分误差界；两个区间分离时，即给一次校准后的可辨后果。使用真实跨块b的公式事前固定，不能事后按每个输出调门户。更粗的g²/δ界仍成立，但随截止增长可能失去选择力，不据此误判实际有限模型。正式数字由代码给出，不预签原B=18原型的数值为B=66结果。

裸线性仍是删除余额合同的对照，不是同合同内反模型；源差亦不能直接等同1322的Higgs源差。f是此模型中真实独立驱动源，进一步认作Higgs必须另有规范/场归一化字典。

## 7. 有限计算如何覆盖完整父过程

实轴驱动保持

$$
\mathcal C=N_1+n_{L_1}+n_U.
$$

χ_B在N_1=0；四物质输入及伙伴只落在C=0或1。前者维2(B+1)，后者维4B+2，总维6B+4。V、D、源S_f及其测量函数保持这份子空间，故全部实际后态/记录可精确计算。C≥2完全位于Q，依据(5)的**全参考隙**不进入低带；不能仅因没有枚举就忽略它们。

作者和独立实现应在小B完整父空间交叉核大窗分块公式；主例B=66计算完整低带、状态、来源及读口。对三/四/五味的直和与一般有限矩阵，谱命题保持任意d的量词；不重新做无选择力的味数扫描，不称该字典选出了代际数。

## 8. 净增与剩余问题

本轮新增的是同一已选原语与动态参考的全低带、伙伴、真实源和反作用连接，以及固定校准后的有限来源/能谱后果。Schur、标准Dirac色散和一般测量膨胀只是复用工具。

新增输入包括参考自由核、驱动、继承动量字典、有效控制及读口；它们尚未成为认知公理推论。原1322标准物质分支、右带荷替代与新参考合同之间的比较必须保各自原量词：不能将一个接口不适用改成整份场论被排除。是否达到阶段所有条件，须以最终联合前件账与五项审计决定；本证明不自动结项。
