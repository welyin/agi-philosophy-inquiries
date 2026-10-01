# 第641轮：同一空间合并过程的非线性记录与几何响应

日期：2026-10-01。接[640](research_note_640.md)。[代码](joint_blocked_record_history.py)、[结果](joint_blocked_record_history_results.json)、[核验](research_round_641_checks.json)、[条件账](unified_physics_condition_ledger_641.md)。两组检查、十二式；主代理审查，无新增独立代理审查。

## 1. 从传播核到实际记录合同

640已给原中性二次场的空间合并、约化热态、时间核和几何归一化。本轮不再重做二频论证，而是把原二值平方根正弦读取函数直接作用于同一个相邻singlet平均，计算两次记录及其几何响应。

**新增正向连接：** 原同一Gaussian过程的等时方差、对称时间关联和对易核，可以共同确定这份非线性仪器的完整两次记录概率。虽然第一次读取后态通常不是Gaussian，以下计算仍保留其真实反作用，无须另拟合一个“读后热态”。

**具体排除：** 用640精确复现等时热态的有效Hamiltonian代替原演化，两次记录分布仍会改变。只保对称关联、删去对易核也会改变记录。均匀几何变化下，仅更新等时参考而固定两种时间核，会漏掉可复算的记录响应。

本轮仍限于原树级二次中性场的正则量子化。原函数的参数未改，但它作用于线性化singlet平均；不声称这是原完整非线性Gauss热态中的精确仪器结果。没有给高阶相互作用误差或自主装置实现。

## 2. 历史和本轮增量

[586](research_note_586.md)已区分真实记录合并与另做平方根粗读；本轮固定同一原函数，不将细记录合并偷换为新仪器。[624](research_note_624.md)已有条件概率及得分，[625](research_note_625.md)已有同图来源共同近似。本轮新对象是640的**实际空间平均**，由原同一三维图及参考计算其非线性两时记录。

Weyl乘法和Gaussian特征函数是成熟量子数学；本轮不将它们另称为理论发现。新增的是原曲目标质量、空间合并、二值操作和几何参数的同一字典。所需原函数、质量与空间对象均从本地冻结结果复用，无新增物种／外部浴。

|层次|范围|
|---|---|
|继承输入|原两中性质量、B、N³周期图、ψ=e^g几何、β热参考，ℏ=ε=1|
|操作|原L±(s)=√(.5±.25 sin s)作用于线性化相邻singlet平均|
|解析|同一Gaussian过程给精确非线性两记录公式，含反作用和几何导数|
|数值|N=8、β=2；完整细格求和与粗格别名求和交叉验证|
|边界|完整Gauss、相互作用、实际实现、连续尺度及GR仍开放|

## 3. 真实空间平均与一份共同过程

在原3维格点上取相邻0、e_x，保留两正规模a=L,H。沿640原B，

$$
S(t)=s_\star+\sum_a\frac{B_{sa}}2
 [q_{a,0}(t)+q_{a,e_x}(t)],\qquad
\mu=s_\star,\quad
V=\langle(S-\mu)^2\rangle,\quad
C(t)=\tfrac12\langle\{S(t)-\mu,S(0)-\mu\}\rangle,\quad
[S(0),S(t)]=iD(t)I.
\tag{1}
$$

这是一个相邻平均量，非640单个空间Fourier模式；代码把所有允许动量重新组合回来。给定均匀几何，
 
$$
A=e^{-6g},\qquad
\omega_a(\mathbf k)^2=m_a^2+
4e^{-4g}\sum_{j=1}^3\sin^2(k_j/2),\qquad
f(\mathbf k)=\frac{1+\cos k_x}{2}.
\tag{2}
$$

原整体二次Gibbs直接给

$$
V=\frac1{N^3}\sum_{a,\mathbf k}B_{sa}^2f(\mathbf k)
 \frac A{2\omega_a}\coth\frac{\beta\omega_a}{2},\qquad
C(t)=\frac1{N^3}\sum_{a,\mathbf k}B_{sa}^2f(\mathbf k)
 \frac A{2\omega_a}\coth\frac{\beta\omega_a}{2}\cos(\omega_at),\qquad
D(t)=\frac1{N^3}\sum_{a,\mathbf k}B_{sa}^2f(\mathbf k)
 \frac A{\omega_a}\sin(\omega_at).
\tag{3}
$$

与640粗动量、两别名频率及平均权的计算一致，最大差1.12×10⁻¹⁶。两套求和避免把一个随意配置当热平均。式(3)只是本明确二次分支的热平均，不是原完整Gauss热迹。

## 4. 原仪器的两记录概率

令r,s∈{+1,−1}，使用同一个原函数，第一记录在0、第二在t：

$$
L_r(S)=\sqrt{\tfrac12+\tfrac r4\sin S},\quad
E_s(S)=L_s(S)^2,\qquad
p_{rs}(t)=\langle L_r(S(0))E_s(S(t))L_r(S(0))\rangle,\qquad
p_r=\tfrac12+\tfrac r4e^{-V/2}\sin\mu .
\tag{4}
$$

L_r顺序不能忽略；并未把第一次条件态重新热化。设L_r(x)=Σ_k c_{r,k}e^(ikx)，c_{r,−k}=c*_{r,k}。同一过程的Weyl身份为

$$
\left\langle e^{ikS(0)}e^{iS(t)}e^{ilS(0)}\right\rangle
=\exp\left\{i[(k+l+1)\mu-\tfrac12(k-l)D(t)]
-\tfrac12[((k+l)^2+1)V+2(k+l)C(t)]\right\}.
\tag{5}
$$

负号由[S(0),S(t)]=iD给出，不独立拟合相位。Gaussian只用于初始整体及线性Heisenberg场，不要求仪器后的状态仍Gaussian。于是

$$
F_r(t)=\sum_{k,l}c_{r,k}c_{r,l}
\exp\left\{i[(k+l+1)\mu-\tfrac12(k-l)D]
-\tfrac12[((k+l)^2+1)V+2(k+l)C]\right\},\qquad
p_{rs}(t)=\tfrac12p_r+\tfrac s4\operatorname{Im}F_r(t).
\tag{6}
$$

此处两个Fourier系数直接相乘；共轭已体现在L自伴及c_−k=c*_k，不能多加一个共轭而改变仪器顺序。

有限K并非新增认知阈值，而是已定义操作的计算误差参数。在|Im z|≤a₀=1.2上，cosh a₀<2，平方根无零点且解析；M₀=√(.5+.25 cosh a₀)控制其模。移围道给

$$
|c_{r,k}|\le M_0e^{-a_0|k|},\qquad
\|L_r-L_{r,K}\|_\infty\le
\delta_K=\frac{2M_0e^{-a_0(K+1)}}{1-e^{-a_0}},\qquad
|p_{rs}-p_{rs}^{(K)}|
\le\frac14\delta_K(2\sqrt{3/4}+\delta_K).
\tag{7}
$$

最后一界用于式(6)中第一概率p_r精确保留、只截断F_r的计算。K=16时δ≤3.860×10⁻⁹，每事件尾界≤1.672×10⁻⁹；离散Fourier网格2048的解析别名尾远小于此，浮点舍入另由数值交叉核验。K=16与24的概率差2.78×10⁻¹⁷，不用这个小差单独冒充严格截断界。

## 5. 同一等时态仍给不同的实际历史

β=2、g=0时μ=.545863690227、V=.272616120626。第一记录p+=(约).613250579888，原模型和有效热Hamiltonian模型精确相同。

对照模型使用每个粗动量上640算出的h_β，保留相同粗变量与仪器，并以其生成时间演化；它精确复现整个粗等时Gaussian态，不只是一个方差。另一对照保留原V、C却设D=0，把过程当经典Gaussian变量。两者的两记录总变差分别为：

|t|改用有效热Hamiltonian|只删除对易核D|
|---:|---:|---:|
|0|0|0|
|.4|.000238135132|.000453885974|
|1|.000401998445|.000291292637|
|2|.000862375348|.000082765384|

这些是声明模型内可区分的结果，不是实验预测的精度承诺。差异远大于式(7)计算尾界。t=0还可独立验证：

$$
p_{rs}(0)=\frac14+\frac{r+s}{8}e^{-V/2}\sin\mu
+\frac{rs}{32}(1-e^{-2V}\cos2\mu).
\tag{8}
$$

每行求和都等于p_r；非选择第一测量仍可能改变第二次边缘，因为它扰动了原态。不以边缘变化宣称超距通信。

## 6. 几何响应也由同一过程决定

沿g变化，均匀体积权归一后相邻平均定义不变，μ不变；639一般非均匀权重的额外操作导数仍适用于另一种比较，此处不能遗漏后又叫通用。

记n=k+l，式(5)指数的几何导数为

$$
\partial_g\log W_{kl}
=-\tfrac12(n^2+1)V'-nC'-\tfrac i2(k-l)D',
\qquad
p_r'=-\frac r8 e^{-V/2}\sin\mu\,V'.
\tag{9}
$$

因此同一Fourier和给

$$
p_{rs}'=\tfrac12p_r'
+\frac s4\operatorname{Im}
\sum_{k,l}c_{r,k}c_{r,l}W_{kl}
[-\tfrac12(n^2+1)V'-nC'-\tfrac i2(k-l)D'].
\tag{10}
$$

本轮V'、C'、D'从式(3)的实际g族差分取得，再用不同步长对完整p_rs差分交叉核对。解析公式是精确链式法则；所报导数值是数值验证，不假称所有导数已获严格区间证书。

由于每个原效应在1/4与3/4之间，所有两记录概率至少1/16，得分可正常定义：

$$
\ell_{rs}=\partial_g\log p_{rs},\qquad
\sum_{rs}p_{rs}\ell_{rs}=0,\qquad
I_g=\sum_{rs}p_{rs}\ell_{rs}^2.
\tag{11}
$$

|t|记录Fisher信息I_g|只更新V、漏C'及D'的最大概率导数差|
|---:|---:|---:|
|.4|.032063713770|.017077352003|
|1|.043073764908|.028335908468|
|2|.051420065702|.032310546270|

完整差分最终误差低于9.35×10⁻¹⁰。几何作用于参考与等待演化，不能仅改初态方差。这里的得分是记录对几何参数的统计响应，**不是应力张量或Einstein方程**。640的几何归一化常数C_β仍须另保，记录概率不能反推出其绝对来源。

## 7. 共同核条件而非独立记录模型

若两份声明的Gaussian过程具有相同均值及全部相关时间上的对称核／对易核，其有序Weyl积相同：

$$
\left\langle\prod_{j=1}^m e^{iu_jS(t_j)}\right\rangle
=\exp\left[
i\sum_j u_j\mu_j-\frac12\sum_{ij}u_iu_jC_{ij}
-\frac i2\sum_{i<j}u_iu_jD_{ij}\right],
\qquad [S(t_i),S(t_j)]=iD_{ij}I.
\tag{12}
$$

对原L_r的绝对收敛Fourier级数，任何固定有限历史可按同一身份输送；原初始关联及操作反作用都包含在有序积内。此为Gaussian分支的条件性充分接口，本轮数值只算两次历史，未宣称任意长度的统一误差界。

这给C03记录、C19参考、C20尺度与C22记录响应一个共用入口。其来源为同一原物质过程，不是新增一个独立噪声浴或记忆控制器。

## 8. 核验和下一步

两组检查覆盖：原非线性函数Fourier恢复、解析尾界、细／粗三维求和、t=0身份、两次真实记录、两种错误替代及几何得分。首版结果保留在[草稿](round641_drafts/first_results_before_tail_bound.json)，随后增加独立细格核对与解析截断界；未覆盖旧结果。

[642入口](round642_drafts/STATUS.md)：返回完整非线性／Gauss空间的尺度条件，核上述共同过程接口能否扩展，而非继续优化Gaussian记录精度。原二次分支不是统一模型的最终验收范围。
