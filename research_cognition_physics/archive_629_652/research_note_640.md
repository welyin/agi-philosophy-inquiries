# 第640轮：一次真实空间合并的共同热态、时间核与几何来源

日期：2026-10-01。接[639](research_note_639.md)及[640入口](640/drafts/STATUS.md)。[代码](640/joint_spatial_block_reference.py)、[结果](640/joint_spatial_block_reference_results.json)、[核验](640/research_round_640_checks.json)、[条件账](640/unified_physics_condition_ledger_640.md)。三组检查、十四式；主代理审查，无新增独立代理审查。

## 1. 合并问题、结论和边界

639有尺度一致的信息上界，却没有证明相邻节点合并后，参考态、演化和来源能用同一份朴素局部Hamiltonian。本轮进行真实相邻格点的平均／差分变换。

**结果：** 原两中性标量的明确二次近似中，空间平均后的热态可精确写为有效热Hamiltonian的Gibbs态；但该Hamiltonian依赖温度，不能直接作为原时间过程的生成元。原过程保留两个频率，而一个正则变量的二次自主Hamiltonian只有一个频率。保留消去模产生的频率核可恢复传播；同一消元还须保几何依赖的配分函数归一化，才能恢复平均几何来源。

需要同一个空间消元规则共同交付**态、时间核和来源归一化**。不能分别拟合后直接签收，也不将此要求命名为新认知公理。

|层次|范围|
|---|---|
|认知动机|改变尺度仍描述同一物质、同一参考和同一几何响应|
|原对象|574正F目标、测地邻接与原势；576不含新增指针的两径向正规模|
|近似输入|平直常真空、给定三维周期图、树级中性二次涨落的正则量子化；ℏ=1|
|比较规则|相邻两点作正则Haar平均／差分，随后迹掉差分变量|
|解析及数值|原空间谱、精确Gaussian约化热态、时间核、共同共形来源|
|未覆盖|完整相互作用Gauss Gibbs、带荷边界、多时非线性仪器、连续与GR|

本轮二次Gibbs不是603完整模型的Gibbs。正规模来自原非线性势Hessian；这里采用树级二次涨落的正则量子化，不包括曲目标的圈或排序修正。尚无有限温度下忽略高阶相互作用的统一误差界。结论严格限于本二次分支，不排除完整模型的所有尺度方案。

## 2. 历史去重与成熟方法

[567](../archive_554_584/research_note_567.md)已有传播及质量消元记忆，[558](../archive_554_584/research_note_558.md)区分有效能源与仪器，不重复一般Schur补或压缩不保乘积。[608](../archive_585_628/research_note_608.md)、[617](../archive_585_628/research_note_617.md)区分独立区域和表示切分，[625](../archive_585_628/research_note_625.md)是同图谱近似。本轮真实移除一半空间坐标。

[612](../archive_585_628/research_note_612.md)是overlap连续谱／热迹的全时障碍；本轮来自原标量空间平均的两个有限原频率，二者不能混同。

[Yeo—Shim，§II](https://arxiv.org/html/2412.02074v1)给约化热态和Gaussian协方差方法；这里被消去的是原格点差分模，不引入其外部热浴。[Morinelli等，§§3.2、5.1](https://arxiv.org/pdf/2010.11121)区分算符细化与状态粗化，指出阶梯剖面的连续动量正则性问题。这里仅用有限图Haar变换，不套用其小波连续极限定理。

新增量是原曲目标、质量和格图中的共同字典及明确失败数值；Gaussian热态与小波理论本身不计作发现。

## 3. 原两标量的共同二次模型

原x=(h,s)、x★=(√u_h,√u_s)、F★=M−|x★|²/6、M=2。用原目标度量和势Hessian选正规矩阵B：

$$
K_\star=\frac I{F_\star}+\frac{x_\star x_\star^\top}{6F_\star^2},
\qquad
V_\star''=\frac{2\operatorname{diag}(x_\star)\mathsf L
\operatorname{diag}(x_\star)}{F_\star^2},
\qquad
B^\top K_\star B=I,\quad
B^\top V_\star''B=\operatorname{diag}(m_L^2,m_H^2).
\tag{1}
$$

复算得m_L²=.0468385104353、m_H²=.138056726486，未增加576的圆指针。令δx=ε^(-3/2)Bq，用原格点积分动量定义p。空间γ_ij=ψ⁴δ_ij、均匀ψ=e^g时，每一正规模：

$$
H_g^{(2)}=\frac12[A_gp^\top p+q^\top V_gq],\qquad
A_g=e^{-6g},\quad
V_g=e^{2g}\epsilon^{-2}L_{\rm lat}+e^{6g}m^2I,\qquad
A_gV_g=m^2I+e^{-4g}\epsilon^{-2}L_{\rm lat}.
\tag{2}
$$

L_lat是给定三维周期图Laplacian。原测地边距离二阶展开与K★相同，两模共用梯度；不是把567缺singlet梯度的早期模型重新当当前574模型。只量子化二次涨落，未声称已经积分掉规范／费米相互作用。横向零动量可抽出二次模型的一个空间方向，不改变给定维数，也不证明完整Gauss空间因子化。

## 4. 实际相邻合并留下两个频率

沿x每两点定义：

$$
Q_j=\frac{q_{2j}+q_{2j+1}}{\sqrt2},\quad
R_j=\frac{q_{2j}-q_{2j+1}}{\sqrt2},\qquad
P_j=\frac{p_{2j}+p_{2j+1}}{\sqrt2},\quad
\Pi_j=\frac{p_{2j}-p_{2j+1}}{\sqrt2}.
\tag{3}
$$

有限正则正交变换有酉实现；随后对差分变量偏迹是CPTP状态映射，原参考态也沿它输送。均匀体积下原singlet两点平均为s★＋Σ_a B_sa Q_a/√(2ε³)。故保两正规模Q就能表示原singlet线性平均；完整非线性读量和多时仪器尚未签收。

粗格动量2k混合细格k及k+π。横向零动量、ε=1：

$$
a=\cos^2(k/2),\quad b=\sin^2(k/2),\quad a+b=1,\qquad
d_0=m^2+4e^{-4g}\sin^2(k/2),\quad
d_1=m^2+4e^{-4g}\cos^2(k/2),\qquad\omega_i=\sqrt{d_i}.
\tag{4}
$$

取N_x=8允许的k=π/4。从8节点原Laplacian直接构造粗实余弦平均及其差分伙伴，它们构成二维不变空间：

$$
L_{\rm pair}=\begin{pmatrix}1&1\\1&3\end{pmatrix},\qquad
\lambda_{0,1}=2\mp\sqrt2,\qquad
(a,b)=\left(\frac{2+\sqrt2}{4},\frac{2-\sqrt2}{4}\right).
\tag{5}
$$

这是空间变换给的权重，不是任意选两个振子拟合质量。其余动量块同法处理；本轮数值仅此真实块的两个中性物种，块配分函数不是全宇宙配分函数。

## 5. 同一状态映射决定有效热态

原二次整体Gibbs中，保留Q、P模为零均值Gaussian态，交叉对称协方差零：

$$
X:=\langle Q^2\rangle=\frac{A_g}{2}\sum_{i=0}^1
\frac{a_i}{\omega_i}\coth\frac{\beta\omega_i}{2},\qquad
P:=\langle P_Q^2\rangle=\frac1{2A_g}\sum_{i=0}^1
a_i\omega_i\coth\frac{\beta\omega_i}{2},\qquad(a_0,a_1)=(a,b).
\tag{6}
$$

不另选Gibbs参考。设ν=√(XP)>1/2，唯一到加性常数的二次有效热Hamiltonian：

$$
\rho_Q=\frac{e^{-\beta h_\beta}}{\operatorname{Tr}e^{-\beta h_\beta}},
\quad h_\beta=\frac12(A_\beta P_Q^2+V_\beta Q^2),\qquad
\Omega_\beta=\frac1\beta\log\frac{2\nu+1}{2\nu-1},\quad
A_\beta=\Omega_\beta\sqrt{X/P},\quad
V_\beta=\Omega_\beta\sqrt{P/X}.
\tag{7}
$$

这是本块平均力Hamiltonian。g=0时：

|原模|β|A_β|V_β|Ω_β|
|---|---:|---:|---:|---:|
|轻|.5|.993340192|.718572533|.844859147|
|轻|2|.921032236|.708907513|.808038781|
|轻|30|.138566512|.114444817|.125929421|
|重|.5|.993539718|.819316453|.902232474|
|重|2|.923676862|.807337249|.863550078|
|重|30|.135411023|.125547218|.130385878|

温度依赖不只是数值现象。原基态在平均／差分划分下纠缠：

$$
\nu_\infty^2=\frac14
\left(\sum_i\frac{a_i}{\omega_i}\right)
\left(\sum_i a_i\omega_i\right)>\frac14
\qquad(ab>0,\ \omega_0\ne\omega_1).
\tag{8}
$$

严格性由Cauchy—Schwarz等号条件给出。两物种ν∞=.523406395647、.520417467115。固定正定单模二次Hamiltonian在β→∞给纯Gaussian基态，不能匹配这份整个约化温度族。若提议更一般的β无关算符，在任一β取满秩Gaussian态的对数已迫使它为相应二次算符加常数，仍受此限。

不排除一个固定温度的匹配，也不要求宇宙实际处于热平衡；这是对当前声明热参考分支的审计。

## 6. 热Hamiltonian不能直接替代原时间过程

原保留Q的Euclidean传播核：

$$
\mathcal C_Q(z)=A_g\left(\frac a{z+d_0}+\frac b{z+d_1}\right)
=\frac{A_g}{S(z)},\qquad
S(z)=z+u-\frac{c^2}{z+v},
\quad u=ad_0+bd_1,\ v=bd_0+ad_1,\ c^2=ab(d_1-d_0)^2,
\quad z=\nu_E^2.
\tag{9}
$$

热Euclidean过程的ν_E取Matsubara频率；同一有理函数也定义解析核。ν_E不是式(7)的辛本征值ν。Schur分母保留差分变量的记忆，S(0)不能替代全部频率。

单模二次自主Hamiltonian的线性Q只有一频，此处两正权、两频不同。用M_n=ad₀ⁿ+bd₁ⁿ可严格区分：

$$
M_0M_2-M_1^2=c^2>0,\qquad
\mathcal K_Q(t):=i[Q(t),Q(0)]
=A_g\sum_i a_i\frac{\sin(\omega_it)}{\omega_i}.
\tag{10}
$$

单频对应矩行列式为零。原块g=0两物种均有c²=1；因此不是靠三次数值采样证明全时失败。β=2时用式(7)演化，在t=.3、.8、1.4的最大对易核误差分别.029568175804、.031122197740；保留原Schur核误差不超过5.6×10⁻¹⁷。

排除的是这个块内“精确热态Hamiltonian＋同一个线性Q”直接充当原全时二次过程。不排除保留差分部门、有记忆过程、非局部时间核或受控有限窗口近似。

## 7. 同一消元还固定几何归一化

ρ_Q不包含所有被消去自由能。对这个实际二模块：

$$
Z_{\rm pair}=\prod_{i=0}^1\frac1{2\sinh(\beta\omega_i/2)},\qquad
Z_\beta^{\rm osc}=\frac1{2\sinh(\beta\Omega_\beta/2)},\qquad
C_\beta(g)=\frac{\log Z_\beta^{\rm osc}-\log Z_{\rm pair}}{\beta},
\quad h_\beta^{\rm full}=h_\beta+C_\beta(g)I.
\tag{11}
$$

Tr_Q exp(−βh_full)=Z_pair且归一态不变。C从原块计算，不是任意新反项。β及原正则坐标固定，原ψ=e^g几何来源为

$$
G_g=\partial_gH_g^{(2)}
=\frac12\left[-6A_gp^\top p+
q^\top(2e^{2g}L_{\rm lat}+6e^{6g}m^2I)q\right]
\qquad(\epsilon=1).
\tag{12}
$$

同一原热迹给

$$
\langle G_g\rangle_{\rm pair}
=-\frac1\beta\partial_g\log Z_{\rm pair}
=\frac12[(\partial_g A_\beta)P+(\partial_g V_\beta)X]
+\partial_gC_\beta .
\tag{13}
$$

仅为热平衡平均来源身份，不代替非平衡响应与噪声。β=2、g=0：

|原模|原完整块来源|仅二次h_β的来源|∂g C_β|恢复总来源误差|
|---|---:|---:|---:|---:|
|轻|−3.039754108|−1.015113166|−2.024640945|3.81×10⁻⁹以内|
|重|−2.892747659|−.897000214|−1.995747449|3.37×10⁻⁹以内|

频率平方矩阵同时满足

$$
(z+d_0)(z+d_1)=(z+v)S(z),\qquad
\partial_g\log[(z+d_0)(z+d_1)]
=\partial_g\log(z+v)+\partial_g\log S(z).
\tag{14}
$$

若使用未归一时间动能的路径积分，还需同时处理A_g和相空间测度。本轮热迹直接从原正则振子谱计算，不靠漏掉测度得式(13)，也不直接相乘无限Matsubara裸行列式。

639有归一相对信息对几何标量能源移位失明的旧限制；本轮进一步计算实际空间消元中不能省掉的C_β和来源。约化态单独不能恢复全部几何账。

## 8. 核验及真正合并的自由

三组首次通过：

1. 原目标、势、8节点空间平均共同决定原质量与空间谱；原非线性势和测地边的二阶差分误差最终8.85×10⁻⁹以内，随步长减小。
2. 同一偏迹的四温度Gaussian态重建误差4.45×10⁻¹⁶以内；直接区分原二频核与热Hamiltonian的一频核。
3. 原几何变分、配分函数、约化归一化复现同一来源，并核频率行列式身份。

C01、C03的一阶读量、C19、C20、C22在这个块上共用一个空间消元：A_β、V_β、记忆核、C_β由原数据确定，不能分别挑选。它们不是新增基础自由参数，完整非线性Gauss闭合仍开放。

[641入口](641/drafts/STATUS.md)检查同一个有记忆区域过程如何交付记录与来源；也保留正则性更好的尺度函数分支，不把有限Haar一步外推为连续动量存在。认知硬件设计后置。
