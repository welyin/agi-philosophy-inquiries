# 相关A/B差值的格点界：混合通道导数独立核验

2026-10-13。未编号数学准入补充；AI内部独立复核。只审同一剪切对照的空间离散误差，不改物理参数、准备、原读口或冻结文件。已读[候选](tensor_candidate/candidate.md)及[诊断预览](tensor_candidate/results_preview.json)。

## 1. 结论

主线提出的因子正确。在第2节固定背景、固定48维码及同一末端效果的量词下，

$$
\boxed{
\left|[p_A^{\rm cont}-p_B^{\rm cont}]
      -[p_A^{\rm lat}-p_B^{\rm lat}]\right|
\le 2dJ+4pdJV,
}
\tag{1}
$$

其中

$$
d=p-\frac{\sin(\ell p)}\ell\ge0,\qquad
J=\int_0^T|g z_t(t)|dt,\qquad
V=\int_0^T\bigl(|\bar v(t)|+|g z_t(t)|\bigr)dt.
\tag{2}
$$

当前$\bar v=1+\nu z_s\ge1-\nu>0$，因此V可写成主线的$\int(\bar v+|gz_t|)$。无需假定剪切不变号。

这是一份**同一A/B差值的界**。它没有把各支格点误差分别缩小，也没有认证插入任意中途背景操作后的历史。它可用于原移动钟×位置效果及同一H的质量零末端位置效果，保任意未知码态、跨质量／动量相干和任意参考R。

诊断预览的带符号积分代入给约$1.3768\times10^{-4}$；该数值尚须第5节的绝对值积分和轨道包络认证。仅见采样点上$z_t\le0$，不能把该数值签成严格界。

## 2. 精确比较对象与量词

令$q(t)=g z_t(t)$。A/B的共同标量背景为N(t)、$\bar v(t)$，两者的剪切符号相反。在所有插值参数中固定这些函数，不随u或λ重算source或G；它们可取待证真实Hartree轨道，或一条明确定义的Hermitian近似轨道，后一方式须另付真实轨道误差。

原准备码为四模式$0,+e_x,+e_y,+e_z$乘Dirac四维及质量三标签$0,.3,.5$，总48维。对非零轴b的块，定义

$$
h_{u,\lambda}^{(m,b)}(t)
=N(t)m\beta+(p-\lambda d)
   [\bar v(t)+u\sigma_bq(t)]\alpha_b,
\quad \sigma=(1,-1,0),
\tag{3}
$$

$$
u\in[-1,1],\qquad \lambda\in[0,1].
$$

零动量块是$Nm\beta$，对u、λ均无导数。λ=0为连续动量比较器，λ=1为实际格点动量；u=±1为A/B，互换A/B名称只翻转所比较差值，不改界。

H是这些块的一个直和算符，**不对输入作块间退相干**。初始密度矩阵可以含任意块间相干，并可与任意R纠缠。以下范数取在整个48维空间上；直和导数范数是各块最大值，故不乘模式数或质量数。质量项虽影响实际U和原钟相位，但其u／λ导数为零；证明不删质量，也不要求质量和动量项对易。

假定N、$\bar v$、q在有限窗内可积且实，H Hermitian。有限维酉演化存在，参数依赖的Duhamel导数可逐次成立。参数插值不是新的物理准备或额外装置，只是对同一两个端点的数学比较。

这一步只离散比较probe；完整高动量source仍按原格点协方差产生背景。它没有将source也作连续低动量近似。

## 3. 混合通道导数：时间排序与全部因子

记$U_{u,\lambda}(t,s)$为式(3)的完整酉传播，$\mathcal C_{u,\lambda}(t,s)=\operatorname{Ad}_{U(t,s)}$。对Hamiltonian变分定义Liouvillian

$$
\mathcal L_r(t)=-i[H_r(t),\,\cdot\,],\quad
\mathcal L_{u\lambda}(t)=-i[H_{u\lambda}(t),\,\cdot\,].
$$

基本完全有界迹范数估计为

$$
\|\mathcal C(t,s)\|_\diamond=1,\quad
\|\mathcal L_r(t)\|_\diamond\le2\|H_r(t)\|,
\quad
\|\mathcal L_{u\lambda}(t)\|_\diamond
\le2\|H_{u\lambda}(t)\|.
\tag{4}
$$

式(4)由任意R上$\|AX-XA\|_1\le2\|A\|\|X\|_1$及酉不变性成立；不是只对无参考的纯态估计。

对$\mathcal C(T,0)$作两次Duhamel微分，得到一项混合插入和两项时间排序插入：

$$
\begin{aligned}
\partial_u\partial_\lambda\mathcal C(T,0)
={}&\int_0^T\mathcal C(T,t)\mathcal L_{u\lambda}(t)
                         \mathcal C(t,0)dt\\
&+\int_{0<s<t<T}\mathcal C(T,t)
 \bigl[\mathcal L_u(t)\mathcal C(t,s)\mathcal L_\lambda(s)
       +\mathcal L_\lambda(t)\mathcal C(t,s)\mathcal L_u(s)\bigr]
 \mathcal C(s,0)\,dsdt.
\end{aligned}
\tag{5}
$$

各中间传播均用当前相同(u,λ)。定义

$$
A_u=\int_0^T\|H_u(t)\|dt,\quad
A_\lambda=\int_0^T\|H_\lambda(t)\|dt,\quad
A_{\rm mix}=\int_0^T\|H_{u\lambda}(t)\|dt.
$$

两种排序在取范数后给

$$
\int_{s<t}\![a(t)b(s)+b(t)a(s)]dsdt
=\left(\int a\right)\left(\int b\right).
$$

这是一个完整矩形，**不再额外乘2**。由(4)—(5)因此有

$$
\boxed{\|\partial_u\partial_\lambda\mathcal C(T,0)\|_\diamond
\le2A_{\rm mix}+4A_uA_\lambda.}
\tag{6}
$$

对式(3)，因$0\le p-\lambda d\le p$且$|\sigma_b|\le1$，

$$
\begin{aligned}
H_u&=(p-\lambda d)\sigma_bq\alpha_b,& A_u&\le pJ,\\
H_\lambda&=-d(\bar v+u\sigma_bq)\alpha_b,& A_\lambda&\le dV,\\
H_{u\lambda}&=-d\sigma_bq\alpha_b,& A_{\rm mix}&\le dJ.
\end{aligned}
\tag{7}
$$

这些上界对整个参数矩形统一成立。

## 4. 从矩形积分到原末端概率

令

$$
\mathcal D=
\mathcal C_{+,0}-\mathcal C_{-,0}
-\mathcal C_{+,1}+\mathcal C_{-,1}.
$$

微积分基本定理给$\mathcal D=-\int_{-1}^{1}du\int_0^1d\lambda\,
\partial_u\partial_\lambda\mathcal C$。u区间长度为2，所以

$$
\|\mathcal D\|_\diamond\le4dJ+8pdJV.
\tag{8}
$$

对同一未知$\rho_{PR}$，$(\mathcal D\otimes\mathrm{id}_R)(\rho)$是Hermitian且迹为零。任意同一效果$0\le E\le I$因此满足

$$
|\operatorname{tr}[E(\mathcal D\otimes\mathrm{id}_R)(\rho)]|
\le\tfrac12\|(\mathcal D\otimes\mathrm{id}_R)(\rho)\|_1
\le\tfrac12\|\mathcal D\|_\diamond.
\tag{9}
$$

(8)的参数区间因子2与(9)的half-diamond恰相抵，得到(1)。不是对一般非零迹误差任意除以2；这里四个酉通道的迹精确抵消是必要依据。

可在$\mathcal D$后接同一终端CP、迹保持记录仪器；其diamond范数不增。保留所有细记录的未条件化输出及其同一汇总因此可使用这份相关差值界。不允许先选成功支再免费归一化；条件概率需另计成功率。

在原连续平滑菜单中，四模式的$f_b=[1+\sin(2\pi x_b/16)]/2$及$f_b\otimes P_{-Y}$经压缩具有与65格点完全相同的Fourier矩阵系数，故可用同一E。这个事实不是“离散位置PVM就是连续精确位置PVM”；若改读口，须重新证明效果／仪器字典。原格点细记录仍是实际采样菜单。

## 5. 变号、数值轨道与绝对值积分

定义$\xi_x=\int(\bar v+q)$、$\xi_y=\int(\bar v-q)$。一般仅有

$$
\begin{aligned}
J&=(\xi_y-\xi_x)/2+2\int_0^Tq_+(t)dt,\\
V&=\xi_y+2\int_0^Tq_+(t)dt
\quad(\bar v\ge0),
\end{aligned}
\tag{10}
$$

其中$q_+=\max(q,0)$。只有认证q全窗非正，才能直接令$J=(\xi_y-\xi_x)/2,V=\xi_y$。采样点不变号不证明点间不变号。

**不必为此增加变号分析工程。** 设某步长度h、数值轨道的cubic Hermite多项式$\hat q$在该步转成Bernstein系数$B_0,\ldots,B_3$。Bernstein基非负，且每个基的积分为h/4，因此直接有

$$
\int_{\rm step}|\hat q(t)|dt
\le\frac h4\sum_{k=0}^3|B_k|.
\tag{11}
$$

若四系数同号，(11)就是该多项式绝对值的精确积分；混号时仍为保守上界。各步求和并为系数与求和舍入向外付预算即可，不需解三次根。更一般次数n时系数为h/(n+1)。

若已证$\int|q-\hat q|\le\delta_q$、$\int|\bar v-\hat{\bar v}|\le\delta_v$，可取

$$
J\le J_{\rm Ber}+\delta_q,\qquad
V\le V_{\rm Ber}+\delta_q+\delta_v.
\tag{12}
$$

例如对q=gz_t，$\delta_q\le|g|\int|z_t-\hat z_t|$。函数$2dJ+4pdJV$对非负d,p,J,V单调，故上界直接代入。

另一合法组织方式是：把(1)应用于**精确定义的Hermitian数值背景**，用(11)认证它的J、V，再在总概率预算另加真实背景与该背景的误差。两种组织任选其一，避免把同一项重复计入。若probe本身用非严格酉的RK输出矩阵，还需单独支付残差／归一化误差；(9)不能对该未规范化误差直接照搬half-diamond。

## 6. 当前预览的独立算术检查

从当前保存预览读取，而非重新积分背景：

$$
p=\pi/8,\quad\ell=16/65,\quad
\xi_x=1.7775632797963896,\quad
\xi_y=1.8687891391179694.
$$

独立普通浮点计算得到

|量|诊断值|
|---|---:|
|精确sine差的浮点值d|0.0006112779835064641|
|解析三次上界$d_* = \ell^2p^3/6$|0.0006115636426094639|
|带符号$(\xi_y-\xi_x)/2$|0.04561292966078989|
|若q≤0，使用d的相关界|0.0001376121384819067|
|若q≤0，使用d_*的相关界|0.0001376764466380996|

在此低动量范围，$0\le d\le d_*$由$0\le x-\sin x\le x^3/6$得到；正式数值可直接向外界定$d_*$，不需另给sine库精度假设。上表仅交叉核量级，**没有认证q的符号、实际背景误差或全数值概率包络**。因此目前可核准定理及约1.38×10^-4的可行量级，不能把上表单独当作已经完成的严格差值证书。

## 7. 允许使用与禁止外推

本界可与各自exact BG的finite-K准备误差、原完整末读实现误差及严格数值包络相加，支撑同一个A/B终端概率差的正余量。两份BG本来可以不同；它们的实际读取与反馈并未被此处抵消。

本界不支持以下推论：

- 每支单独的连续／格点全通道误差均降到(1)。当q=0时(1)为零，真实含义只是A/B差值同时为零，而两支自身都可有非零格点偏差。
- 不同A/B输入态、不同效果或另选轨迹也自动享有相同抵消。若改变其中任何项，须重建同一个参数矩形。
- source全Brillouin离散误差、任意BG干预或中途把probe送出四模式码后的未来也被覆盖。
- 非酉浮点矩阵、两个不同归一化的后选择结果仍可以免费使用(9)。

**审核结论：解析界通过，因子无缺失；正式应用等待J/V及已有数值轨道责任的外包络。** 它是利用已固定对照结构收紧比较误差，没有调整物理参数、改见证或另开科研轮。
