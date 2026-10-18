# 第1295轮证明：非退化不可约味资料与有限未知记录的共同开族

2026-10-18。依据[事前准入](../_shared/round1295_admission.md)。本轮只证明明确二次背景有效类中的联合实例；不新增公理，不把一般Duhamel、CAR或1011分类重新计为发现。三代与四代分别满足本有限合同，不据此签收完整SM或整个比较类。

## 1. 原字段、固定参数与两个不同的背景问题

沿[1293 §1—2](../1293/proof.md)取标准荷 $k=0$、含原 $L,e^c,N,H$ 的Dirac分支，$n=3,4$，$\hbar=c=1$、度规签名 $(+---)$。保同一原费米二次作用 $S^{(2)}_{F,n}[g,A,H;Y]$，包含全部原物种；先对背景变分，再取

$$
H_0=(0,v/\sqrt2),\quad v=\sqrt2,\quad A_\mu=0,\qquad
\alpha_*={31\over2^{27}},\quad \epsilon={\alpha_*\over128}={31\over2^{34}},
\tag{1}
$$

$$
D_n=\operatorname{diag}(1,\ldots,n),\quad(P_n)_{ij}=\mathbf1_{|i-j|=1},\quad
A_n=D_n+P_n,\qquad Y_\nu=I_n+\epsilon A_n,\quad Y_e=2I_n+\epsilon D_n.
\tag{2}
$$

其余耦合、原字段及参考来源均为各模型固定输入，不把旁观部门的绝对源置零。$Y_\nu,Y_e$ 是同作用中的真实常数，不是认知公理的结论。“全矩阵”表示保留全 $n$ 维作用及其混合，**不表示每个矩阵元都非零**。原 $\widetilde H^\dagger L$／$N^\dagger$ 字段身份、离散 $r$、相同宇称及经典非零Higgs背景的精确归一化字典按1293及[背景着装接续](../_admission/background_dressing_screen_after1293.md)继承；不改称动态量子Higgs复合场的精确CAR。

在本背景，中微子一粒子生成元及平方为

$$
h_\nu=\begin{pmatrix}-\boldsymbol\sigma\cdot p\otimes I_n&I_2\otimes Y_\nu\\
 I_2\otimes Y_\nu^\dagger&\boldsymbol\sigma\cdot p\otimes I_n\end{pmatrix},\qquad
h_\nu^2=\operatorname{diag}\bigl(-\Delta+Y_\nu Y_\nu^\dagger,-\Delta+Y_\nu^\dagger Y_\nu\bigr)\otimes I_{\rm spin}.
\tag{3}
$$

排列张量指标不改变式(3)。带电块用 $Y_e$。与 $Y_{\nu,0}=I_n,Y_{e,0}=2I_n$ 的同主部基点相比，差为有界零阶矩阵。因此整个保留二次生成元满足

$$
\|h-h_0\|\le\epsilon\max(\|A_n\|,n)\le\epsilon(n+2),\qquad
 d:=\sup_{0\le t\le1/2}\|U_h(t)-U_0(t)\|\le3\epsilon.
\tag{4}
$$

最后一步是原Cauchy空间的Duhamel界。此处未把变化的Yukawa假作同一模型中的免费控制，也不认为任意度规主符号变化都受此算符范数控制。一般背景的主部／来源估计应回用[803 §4—6](../../archive_797_818/research_note_803.md)、[805 §4—5](../../archive_797_818/research_note_805.md)。

## 2. 正性、简单谱及共同不可约性

$\|P_n\|\le2$、$D_n\ge I$，故 $Y_\nu\ge(1-\epsilon)I>0$；$\|Y_\nu\|\le1+6\epsilon<2$。$Y_e$ 的正本征值 $2+i\epsilon$ 两两不同。$Y_\nu$ 是非零副对角的实对称Jacobi矩阵：任一本征向量的第一分量不能为零，且后续分量由首分量及三项递推唯一决定。因此每个本征空间至多一维；Hermitian性给全部谱简单。正性保证奇异值也简单。

令 $H_\nu=Y_\nu Y_\nu^\dagger$、$H_e=Y_eY_e^\dagger$。$H_e$ 简单且对角，故与它对易的Hermitian $Q_L$ 必为对角。另一方面

$$
(H_\nu)_{i,i+1}=2\epsilon+\epsilon^2(2i+1)>0,\quad 1\le i<n.
\tag{5}
$$

于是 $[Q_L,H_\nu]=0$ 迫使相邻对角元相同，继而 $Q_L=qI_n$。这直接调用[1011 §3—4](../../archive_1009_1043/research_note_1011.md)的共同Yukawa分类；左右交织和Hermitian性均按该原定理保留。不存在共同保留的孤立两代块，也不存在1293那份精确两代味 $SU(2)$。此为精确参数合同；微小非零边并不意味着已有仪器能以该精度测出不可约性，1011 §7的误差限制仍在。

## 3. 每模型自身参考、共同未知输入及局部盲性

每份 $h$ 使用**自身**正常Hadamard参考 $\omega_h$；平直常质量块取 $P_h=\mathbf1_{(-\infty,0)}(h)$。通过质量奇异值分解得到有限多个通常质量Dirac真空，故此参考合法。不同质量方程的参考不能直接互换；不要求不同模型在一个预定全局Fock表示中酉等价。旧[730 §3—5](../../archive_702_741/research_note_730.md)提供一般矩阵质量工具。

取1293原 $a=1/32,L=1/4,I=[7/16,1/2]$、正光滑紧支bump $b_A$、平移 $b_B$ 和 $g_{A,B}=b_{A,B}/\|b_{A,B}\|_2$；不作正频投影后再称紧支。固定单位二旋量 $\zeta$，

$$
f_i=(0,g_A\zeta)\otimes e_i,\quad h_i=(g_B\zeta,0)\otimes e_i\quad(i=1,2),\qquad
 V=(f_1,f_2),\quad F=VV^\dagger,\quad Q_\perp=1-F.
\tag{6}
$$

对每模型取 $\sigma_{h,\perp}=\omega_h|_{\operatorname{CAR}(Q_\perp\mathcal K_n)}$ 的整个限制。令 $W$ 将逻辑二维空间等距嵌入两模的单占据子空间，定义同型准备

$$
\mathcal E_h(\rho_{QR})=(W\otimes I_R)\rho_{QR}(W^\dagger\otimes I_R)
\ \widehat\otimes\ \sigma_{h,\perp}.
\tag{7}
$$

所有码密度矩阵为偶，分级乘积合法。它公开了一份正常、有限模式CP替换初态合同，未声称原参考本已因子化或准备器已自主制造。共同其余态可随模型改变，但不随未知输入改变。对于初远端完整偶代数 $\mathfrak B_0\subset\operatorname{CAR}(Q_\perp\mathcal K_n)$，有

$$
\omega_{\mathcal E_h(\rho)}(B B_R)=\omega_h(B)\operatorname{Tr}(\rho_R B_R).
\tag{8}
$$

此式由完整态(7)成立，不靠二点相同推断。任意伙伴与全部后态按同一演化保留。二点函数仅用于二次量：

$$
C_{h,\rho}=Q_\perp P_hQ_\perp+V\rho_QV^\dagger,\qquad
\Delta C_{h,\rho}=-FP_h-P_hF+FP_hF+V\rho_QV^\dagger.
\tag{9}
$$

一般混合单占据态不是准自由态；例如 $\langle n_1n_2\rangle=0$ 而错误Wick替代给 $\det\rho_Q$。有限光滑CAR操作与传播保持相同短距及固定测试的来源域，复用730、[813 §5—7](../../archive_797_818/research_note_813.md)、[818 §3、8](../../archive_797_818/research_note_818.md)。

## 4. 全窗有限记录与未知相干

令二维接收幅矩阵 $A_h(t)_{ji}=\langle h_j,U_h(t)f_i\rangle$。基点同质量给 $A_0(t)=\alpha_0(t)I_2$。1293已用严格时样KG尾核和同相位正bump证明全窗 $|\alpha_0(t)|>\alpha_*$，该连续证书直接复用；不重新求传播积分。

压缩是收缩，故由(4)

$$
\|A_h(t)-\alpha_0(t)I_2\|\le d\le3\epsilon={3\alpha_*\over128},\qquad
 s_{\min}(A_h(t))>{125\alpha_*\over128}>0.
\tag{10}
$$

以同一局部偶效果 $E_1(t)=a_t^\dagger(h_1)a_t(h_1)$ 比较逻辑 $|0\rangle,|1\rangle$，每模型两输入的共同参考基线相消：

$$
\Delta p_h(t)=|(A_h)_{11}|^2-|(A_h)_{12}|^2
\ge(|\alpha_0|-d)^2-d^2
>\alpha_*^2-2\alpha_*d\ge{61\over64}\alpha_*^2.
\tag{11}
$$

这里使用 $d<\alpha_*$ 及右式对 $|\alpha_0|$ 单调。归一非负光滑时间权的平均保持该下界。有限矩阵指数只校准Duhamel，不能替代这一连续时空证明或领取实际读数精度。

接收双线性的输入依赖矩阵是 $A_h\rho_QA_h^\dagger$（指标按二点约定转置），故(10)使全部未知二维密度矩阵的这份有限信息映射为单射，并线性延拓至伙伴。即使 $\rho_Q$ 含相干也不丢成经典标签。此处不声称由这些效果无噪恢复整qubit、已有层析仪器，或接收端与其余系统无关联。

初远端盲性由(8)保持；早于两球因果接触的接收代数没有输入差。当前比较的各模型均为自己的自治常质量演化，没有质量quench。

## 5. 完整背景来源是矩阵身份，而非味盲标量

固定全部物理耦合，以同一局部协变、保背景Ward的重整化规则定义每模型的

$$
J_\zeta=\delta_\zeta S_{F,n}^{(2)},\qquad
\delta S_{F,n}^{(2)}=\int\sqrt{-g}\left[-\tfrac12T^{\mu\nu}\delta g_{\mu\nu}
+j_a^\mu\delta A_\mu^a+\mathcal F_I\delta H^I\right].
\tag{12}
$$

$H^I$ 为全部四个实Higgs方向。原带电块、$W^\pm$、Higgs相位／径向与spin连接先保留再代背景；其它物种及参考绝对源、固定局部反项均保留。不同 $P_h$ 不赋予逐模型任意重调基线的权限。任意Yukawa spurion变分不属本次来源菜单。

把同一固定紧支源测试运输到初片，记其二次算符核为 $\mathfrak j_{h,\zeta}$。在所用光滑形式域中，准备的拉回准确为

$$
\mathcal E_h^*(J_\zeta)=c_{h,\zeta}I_2+B_{h,\zeta},\quad
B_{h,\zeta}=V^\dagger\mathfrak j_{h,\zeta}V,\quad
c_{h,\zeta}=\omega_h(J_\zeta)+\operatorname{Tr}\{\mathfrak j_{h,\zeta}(-FP_h-P_hF+FP_hF)\}.
\tag{13}
$$

迹只涉及光滑有限秩腿，不是未减除的无限海迹。对伙伴效果 $B_R$，联合期望为
$c_{h,\zeta}\operatorname{Tr}(\rho_RB_R)+\operatorname{Tr}[\rho_{QR}(B_{h,\zeta}\otimes B_R)]$。$B_{h,\zeta}$ 一般非标量，允许输入依赖及与伙伴的来源关联；不再借1293的简并Schur结论。带电混合腿的压缩在当前无交叉参考／中性码上可为零，但该结果不能代替先变分全部顶点。

在壳保持原背景Ward：紧支规范变化及微分同胚给 $J_{\delta_\eta(g,A,H)}=0$、$J_{\mathcal L_X(g,A,H)}=0$，其中规范Ward含
$D_\mu j_a^\mu+\mathcal F_I(t_aH)^I=0$。单独左弱流一般不守恒；同Yukawa的Higgs相位项配平。仅在所取平直、零规范势、常H背景上，完整费米应力化为普通平移守恒。背景量子反常范围继承原完整物种证书，不由两味码重算或删除其余场。

有限性来自每份 $\zeta\in C_c^\infty$ 测试及(9)的光滑有限秩腿。式(4)的有界传播差**不控制**无界应力、源高阶矩或跨模型来源误差；若比较来源，须使用805的实际Sobolev／导数域。这里给精确二次源身份和有限固定测试，不签全部噪声预算或自引力解。

## 6. 正有限相对能源

对自身平直真空参考 $P_h=P_-(h)$，实径向右手模式使 $FhF=0$：动量平均为零，质量块将右手送入正交左手。对任意未知输入均有

$$
E_{\rm rel}=\operatorname{Tr}(h\Delta C_{h,\rho})
=-2\operatorname{Tr}(FhP_h)=\operatorname{Tr}(F|h|)>0.
\tag{14}
$$

$P_h$ 的其它部门保持原参考，不新增相对能；绝对基线仍保留。由(3)、$\|Y_\nu\|<2$、1293的 $\|\nabla g_A\|^2<512/a^2$ 和平方根凹性，

$$
E_{\rm rel}\le2\sqrt{\|\nabla g_A\|^2+\|Y_\nu\|^2}
<2\sqrt{524292}<1450.
\tag{15}
$$

这是该正常初态相对自身参考的正总能源，非准备器成本、局部正能、全部源矩阵恒定或弱引力预算。此能量期望碰巧与未知输入无关，不能据此把(13)全部替换为标量。全等待过程的总相对能守恒；原光滑紧支模式保证每份固定能量矩有限。

## 7. 开族、删除及准入结论

(2)点上的质量谱间隔、最小奇异值和共同对易映射在去标量矩阵空间上的最小非零奇异值均严格正。有限维连续性使三项在某邻域保持。一个明确的非空开域为

$$
\mathcal O_n=\{(Y_\nu,Y_e):\ \|Y_\nu-I_n\|<\alpha_*/2,\quad
\|Y_e-2I_n\|<1,\quad
Y_\nu,Y_e\text{ 的奇异谱简单},\quad
\{Y_\nu Y_\nu^\dagger,Y_eY_e^\dagger\}'=\mathbb CI_n\}.
\tag{16}
$$

中心(2)属于该域；质量均非零。共同核为标量等价于去标量矩阵空间上的对易映射满列秩，是开条件。允许一般复矩阵；“正”指健康动能及正的非零质量奇异值，不要求Yukawa矩阵保持Hermitian。

整个域上，中微子记录只使用 $h_\nu$，故 $d_\nu\le\tfrac12\|Y_\nu-I_n\|<\alpha_*/4$，并给
$\Delta p_h>\alpha_*^2/2$、$s_{\min}(A_h)>3\alpha_*/4$。这里**不是**断言整个费米生成元的差小于 $\alpha_*/4$：带电块在(16)中可有更大差。所取 $A=0$、常H背景的二次等待把中性／带电传播分开，故此限制合法；源仍在变分前保全部混合腿。中心保较强(10)—(11)。另可取更小邻域保持中心附近任意严格弱化的强界，不优化半径。每模型参考仍为自身 $P_h$；微小谱隙不授现实分辨率。

删除路径混合而保同字段与准备合同，两个Yukawa均对角，共同对易核升为所有对角矩阵；记录可仍存在。因此一般记录能力不等于不可约。删除全部中微子Dirac质量、使用删除模型自身参考，右→左块严格为零；这个特定任务失败，不否定其它质量机制或记录。

本轮新内容是同一实际二次过程同时具有：全 $n$ 非简并不可约Yukawa、正常未知输入及初远端盲性、全窗左手量子记录、同一完整矩阵来源与正有限能源。它把1011的结构前件接到1293的有限过程，说明共同资格不依赖孤立两代或精确味简并。三个与四个代际都成立，但未成为全部认知角色的代际自由反证；完整非零相互作用SM、量子规范/Higgs复合场高阶、末读装置及其误差没有被赋零。

[代码](check.py)默认及 `--check` 只读重算，`--write` 仅首次排他创建[结果](results.json)。严格结论来自上述代数和连续证明；一组有限矩阵／CAR／Ward诊断只作交叉校准。待主线独审验收，不自授最终签收。
