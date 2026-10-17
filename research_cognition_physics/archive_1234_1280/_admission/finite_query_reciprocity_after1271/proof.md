# 原有限查询帧的定量反作用桥

2026-10-17。**非编号准入引理**；不预约1272，不新增或采用认知公理，不领取实际几何实例。本文证明一个有限维条件桥：已认证的前向概率差，在完整闭合酉和原许可联合查询帧的条件下，给出反向原查询的正下界。它不再另设未知的反向扰动强度。

一般“闭合双系统酉若单向无影响便为乘积”直接复用[176轮§§3—4](../../../archive_171_189/research_note_176.md)，不是本次新定理。[400轮§§2—5](../../../archive_370_428/research_note_400.md)的同过程、原输入域和有限误差范围，以及[1157筛查§§2—4、6](../../../archive_1120_1158/1157/candidate_screen.md)的完整代数条件继续适用。净承接仅为：**前向实际见证—归一Choi恒等式—原许可有限查询帧—反向概率预算**。算符纠缠是成熟概念，见[Zanardi, Entanglement of Quantum Evolutions](https://arxiv.org/abs/quant-ph/0010074)；下文给出本引理所需全部代数证明，不作一般理论优先权声明。

## 1. 对象、顺序与许可量词

设$G,M$为固定有限矩阵因子，维数分别为$g,m\ge2$。本段完整过程为

$$
U:\mathcal H_G\otimes\mathcal H_M
 \longrightarrow\mathcal H_G\otimes\mathcal H_M,\qquad U^\dagger U=UU^\dagger=I .
\tag{1}
$$

控制、库存、程序和环境若在本段参与作用，必须纳入这两个因子的明确划分，不能先固定并删去一个环境，再对所得开放通道使用下面的纯Choi恒等式。把控制归入$G$或$M$也不自动授予控制的任意准备。

输入输出各自选定基。Choi腿顺序始终为

$$
(G_o,M_o,G',M'),\qquad
|U\rangle=\frac1{\sqrt{gm}}
 \sum_{a,\mu,i,\nu}U_{a\mu,i\nu}|a,\mu,i,\nu\rangle,\qquad
\rho=|U\rangle\langle U| .
\tag{2}
$$

撇号是输入副本；它只用于通道表示，不要求实际制备Choi纠缠态。所有转置均相对于式(2)的固定输入基。

**原许可查询要求。** 反向查询的效果$F_j$、$G$准备$\tau_a$以及$M$的两份准备$\sigma_b^\pm$，须为事前登记、资源和到场范围清楚、能够联合执行的有限菜单。个别操作各自许可，不足以证明它们的任意组合许可。未来几何查询只有在已证明它由同一输出$G$经固定许可过程产生、其效果可拉回为$F_j$时，才可进入本菜单；不得偷偷让改变的$M$另经直接通道写入结果。

拉回效果可以来自不同的原许可未来设置，不要求它们共同可交换或在单份未知态上同时读取。每项比较可以使用一份独立、合法准备的运行。后续查询自己的承担者、来源和实际权限仍须记账，不能由有限矩阵效果的存在授予。

## 2. 双向Choi差与精确HS身份

记

$$
B=\operatorname{Tr}_{M_o}\rho,\qquad
A=\operatorname{Tr}_{G_o}\rho,\qquad
S_G=\rho_{G_oG'},\qquad S_M=\rho_{M_oM'} .
\tag{3}
$$

$B$的腿顺序为$(G_o,G',M')$；$A$的腿顺序为$(M_o,G',M')$。定义

$$
R_b=B-S_G\otimes\frac{I_{M'}}m .
\tag{4}
$$

为避免把$R_f$的中间腿写错，以矩阵元定义

$$
(R_f)_{\mu i\nu;\,\tilde\mu\tilde i\tilde\nu}
=A_{\mu i\nu;\,\tilde\mu\tilde i\tilde\nu}
-\frac{\delta_{i\tilde i}}g
 (S_M)_{\mu\nu;\,\tilde\mu\tilde\nu}.
\tag{5}
$$

下标$b$表示$M$输入对$G$输出的反向影响，$f$表示$G$输入对$M$输出的前向影响。HS范数记为$\|X\|_2^2=\operatorname{Tr}(X^\dagger X)$，不是迹范数。

纯Choi态及酉性给出

$$
\rho_{G_o}=\frac{I_g}g,\quad
\rho_{M_o}=\frac{I_m}m,\quad
\rho_{G'M'}=\frac{I_{gm}}{gm},\quad
\operatorname{Tr}B^2=\frac1m,\quad
\operatorname{Tr}A^2=\frac1g .
\tag{6}
$$

$G_oG'$与$M_oM'$互为纯态的互补约化，所以其纯度相同，记

$$
p=\operatorname{Tr}S_G^2=\operatorname{Tr}S_M^2.
$$

此外，

$$
\operatorname{Tr}\!\left[B\left(S_G\otimes I/m\right)\right]
=\frac pm,\qquad
\operatorname{Tr}(S_G\otimes I/m)^2=\frac pm .
$$

展开式(4)的平方，式(5)同理，得到

$$
\boxed{
m\|R_b\|_2^2=g\|R_f\|_2^2=1-p .
}
\tag{7}
$$

这是本次使用的精确身份。它不说任意指定准备下的两向信号一样大，也不说各种因果影响度量相等。

## 3. 反向差所在的有限实向量空间

式(4)显然满足$\operatorname{Tr}_{M'}R_b=0$。又由式(6)，

$$
\operatorname{Tr}_{G_o}B=\frac{I_{G'M'}}{gm},\qquad
\operatorname{Tr}_{G_o}S_G=\frac{I_{G'}}g,
$$

故$\operatorname{Tr}_{G_o}R_b=0$。在各因子的Hermitian HS正交基中，将单位矩阵列为第一项；这两个偏迹条件恰好去掉$G_o$或$M'$为单位的基分量。因此

$$
R_b\in W:=
\operatorname{Herm}_0(G_o)\otimes_{\mathbb R}
\operatorname{Herm}(G')\otimes_{\mathbb R}
\operatorname{Herm}_0(M'),
\qquad
D:=\dim_{\mathbb R}W=(g^2-1)g^2(m^2-1).
\tag{8}
$$

$\operatorname{Herm}_0$表示无迹Hermitian矩阵。选定$W$的HS正交归一实基$\{W_\alpha\}_{\alpha=1}^D$。本引理认证的是整个$W$上的帧；若拟使用更小空间，必须另证该预定空间包含所用实际$R_b$，不能看过结果后删掉不可见方向。

## 4. 原前向概率见证给出反向HS下界

前向菜单采用两份实际许可$G$准备$\rho_G^\pm$、同一$M$准备$\sigma_M$和输出效果$0\le E_M\le I$。令

$$
\Delta\rho_G=\rho_G^+-\rho_G^-,
\quad
\delta_f=
\left|\operatorname{Tr}
 [(I\otimes E_M)U(\Delta\rho_G\otimes\sigma_M)U^\dagger]\right|,
\quad
A_f=\|E_M\|_2\|\Delta\rho_G\|_2\|\sigma_M\|_2 .
\tag{9}
$$

下面只讨论$\delta_f>0$，因而$A_f>0$。归一Choi公式含因子$gm$：

$$
\delta_f
=gm\left|\operatorname{Tr}
 [R_f(E_M\otimes\Delta\rho_G^T\otimes\sigma_M^T)]\right|
\le gm A_f\|R_f\|_2 .
\tag{10}
$$

式(5)中减去的乘积项与$\Delta\rho_G^T$的无迹性配对为零；未使用免费后选择。结合式(7)，

$$
\|R_b\|_2\ge
 \sqrt{\frac gm}\frac{\delta_f}{gm A_f}.
\tag{11}
$$

这是由同一个前向见证产生的下界，没有再假定另一项反向扰动$\eta$。

## 5. 完整原许可联合查询帧

记$\Delta\sigma_b=\sigma_b^+-\sigma_b^-$。把实际允许的联合测试编号为$\ell=1,\ldots,J$，每行对应一个许可三元组$(j,a,b)$，定义

$$
T_\ell=F_j\otimes\tau_a^T\otimes\Delta\sigma_b^T,
\qquad
\mathsf S_{\ell\alpha}=\operatorname{Tr}(T_\ell W_\alpha),
\qquad
s_*=\sigma_{\min}(\mathsf S)>0 .
\tag{12}
$$

$\mathsf S$为实矩阵；这里$s_*$表示满列秩意义下的最小奇异值，故要求$J\ge D$。效果的单位分量不影响它在$W$上的投影。几何效果本身信息完备仍不够：实际$G$准备和$M$准备差也必须共同使式(12)成立。

每行真正比较的两个概率之差为

$$
q_\ell=
\operatorname{Tr}[(F_j\otimes I)
 U(\tau_a\otimes\Delta\sigma_b)U^\dagger]
=gm\operatorname{Tr}(R_bT_\ell).
\tag{13}
$$

令$x_\alpha=\operatorname{Tr}(R_bW_\alpha)$，则$\|x\|_{\ell^2}=\|R_b\|_2$。由最大分量、最小奇异值和式(11)依次得

$$
\begin{aligned}
\max_\ell|q_\ell|
&\ge \frac{gm}{\sqrt J}\|\mathsf Sx\|_{\ell^2}\\
&\ge\frac{gm\,s_*}{\sqrt J}\|R_b\|_2\\
&\ge
\boxed{\frac{s_*}{\sqrt J\,A_f}\sqrt{\frac gm}\,\delta_f}.
\end{aligned}
\tag{14}
$$

因此至少一项**原菜单内**的反向记录有正差。它可能使用与原前向见证不同的$G$准备；本式不保证指定初态、指定单个几何效果或任意未知旁观态都出现同一下界。

本证明不要求实际生成未知态的副本。$U\otimes I_R$当然保留任意被动参考的完整联合演化；式(14)则是一项声明的乘积准备菜单存在性结论，不能将其改写成对所有未知参考态的统一信号下界。

## 6. 每支误差和认证常数

先将真实概率与同一理想$U$及菜单之间的误差公开：前向两支每支至多$e_f$，全部反向测试的两支每支至多$e_b$。这些预算可来自准备、过程、末效果和统计估计的已认证界，须保持同一对象；它们不自动控制能源或应力误差。

若已认证的前向概率估计差为$\delta_{\rm obs}$，则模型前向差至少为$(\delta_{\rm obs}-2e_f)_+$。取帧的认证下界$\underline s\le s_*$及$\overline A\ge A_f$，有

$$
\boxed{
\max_\ell|\Delta p_{\ell,\mathrm{actual}}|
\ge
\left[
\frac{\underline s}{\sqrt J\,\overline A}\sqrt{\frac gm}
 (\delta_{\rm obs}-2e_f)_+
 -2e_b
\right]_+ .
}
\tag{15}
$$

若这些界只在统计置信事件上成立，结论亦只在同一联合事件上成立，必须保留失败概率；单次二元样本不是$\delta_{\rm obs}$的概率认证。式(15)严格为正才有有限可分辨证书。

校准帧若满足$\|\mathsf S-\widehat{\mathsf S}\|_{\rm op}\le e_S$，可用

$$
\underline s=\sigma_{\min}(\widehat{\mathsf S})-e_S>0 .
\tag{16}
$$

不能把浮点“非零秩”或仅定性的可辨性当作严格正的$\underline s$。各项误差按其实际定义计入，不能既把某项校准误差计入$e_b$又把它当成完全准确的物理帧。反之，满足精确模型输入时可取相应误差为零，但不等于已经实验认证。

## 7. 精确非空例：36项Pauli菜单

这一例只证明有限数学条件非空，**不把Pauli效果认作实际几何，不授予CNOT的原物质实现或来源**。

取$g=m=2$，令$X,Y,Z$为Pauli矩阵。使用

$$
F_i=\frac{I+\sigma_i}{2}\quad(i=x,y,z),\qquad
\{\tau_a\}_{a=0}^3=\{I/2,F_x,F_y,F_z\},\qquad
\Delta\sigma_b=\sigma_b
=\frac{I+\sigma_b}{2}-\frac{I-\sigma_b}{2}.
\tag{17}
$$

全部$3\cdot4\cdot3=36$项组合均在本例菜单中。取$W$的归一Pauli张量基。$F_i$的无迹系数矩阵是$I_3/\sqrt2$，$\Delta\sigma_b$的系数矩阵是$\sqrt2 I_3$；转置只给$Y$列一个符号，不改奇异值。中间准备帧为

$$
\mathsf A=\frac1{\sqrt2}
\begin{pmatrix}
1&0&0&0\\
1&1&0&0\\
1&0&1&0\\
1&0&0&1
\end{pmatrix},\qquad
\mathsf A^T\mathsf A=\frac12
\begin{pmatrix}
4&1&1&1\\
1&1&0&0\\
1&0&1&0\\
1&0&0&1
\end{pmatrix}.
\tag{18}
$$

完整帧是上述三个帧的张量积，前后常数相消。式(18)的特征值为$1/2,1/2,(5+\sqrt{21})/4,(5-\sqrt{21})/4$，所以

$$
J=D=36,\qquad s_*^2=\frac{5-\sqrt{21}}4.
\tag{19}
$$

取$U$为$G$控制、$M$目标的CNOT。前向选$G$的$|0\rangle,|1\rangle$，$M$准备$|0\rangle$，效果$E_M=|1\rangle\langle1|$，则$\delta_f=1$、$A_f=\sqrt2$。式(14)给出

$$
\max_\ell|q_\ell|
\ge\frac{\sqrt{(5-\sqrt{21})/4}}{6\sqrt2}
>\frac1{30}.
\tag{20}
$$

最后的严格有理护栏可用$\sqrt{21}<23/5$，从而$s_*>3/10$，以及$\sqrt2<3/2$。本例实际最大反向差为1：$G$准备$|+\rangle$、$M$分别准备$|+\rangle,|-\rangle$，用$F_x$读取即可。所有状态正常归一，比较没有后选择；实际准备的纯度、库存和相干参考仍是具体物理模型必须承担的资源。

## 8. 删除检验、物理范围与停止线

1. **只保经典标签不足。** 对同一CNOT，若原几何菜单只含$Z_G$的函数，所有几何标签在交互中保持不变，尽管前向差可为1。此时帧在$W$上有核，式(12)不成立。不能将缺失的$X/Y$读口临时补入后称原几何已反应。
2. **效果完备不补准备权限。** 即使允许全部Pauli读取，若$G$只能准备$Z$对角态，CNOT仍不给这份实际准备族反向信号。联合帧仍退化。完整资源要求可能恰好固定某个程序或时钟准备，使实际帧不满秩；不能用“它属于完整量子母空间”解除该限制。
3. **开放过程不能套纯度身份。** 固定并迹掉第三方后的通道通常有混合Choi态，式(6)—(7)不再成立。把第三方收入$M$后，也必须重新核实际准备帧，不能免费改动原ready状态。
4. **未来读取不等于免费层析。** [391轮§§3、7](../../../archive_370_428/research_note_391.md)已区分共同尖锐标签、信息完备统计和单份未知态恢复。这里的有限菜单仍需原访问合同；重复准备不等于克隆当前未知态。
5. **非零记录不等于曲率。** [962轮§§4—6](../../../archive_956_989/research_note_962.md)已有原材料关系反馈，仍未成为几何。本桥若仅影响一个内部钟、尺或器件，也不能推出共同度规、曲率、普适响应系数或Einstein方程。
6. **不得靠定义几何排除反例。** [热边界角色准入§2](../thermal_boundary_role_screen.md)已禁止将原几何资料事后扩充成完整状态层析。这里只给一个条件结论；尚未证明目前的公共关系承担者及原有限未来几何菜单满足式(12)。

本次可关闭的数学差额是：不再先假设未知的反向强度，再从可辨性推出几何记录。待验的真实角色仍是**原许可联合查询帧及其资源、未来支持和标定**。在它未通过之前，本引理不否决1271的固定平直共同模型，不宣告认知公理不足或广义相对论已完成。没有实际共同实例时，不据本引理继续制造仪器或展开一轮候选修补。
