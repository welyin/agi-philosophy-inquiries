# 1176解析证明：公共尺度核中的来源如何进入原共同无迹运输

2026-10-13。本轮采用明确新增的有限接续规则，证明其条件结论。没有新增采用公理。以下“完整”指给定有限CAR父模型的全部材料、关系自旋、交互及保留记录，不是完整现实物理。

## 1. 同一有限局部父过程

复用[1175](../1175/proof.md)。格点 $\Lambda=(\mathbb Z/65\mathbb Z)^3$、$V=65^3$、周期 $L=16$、格距 $\ell=16/65$。源每站4个CAR模式、质量 $m_S=2$；探针每站12个CAR模式，质量标签 $M_P=\operatorname{diag}(0,.3,.5)$。取 $\beta=\sigma_z\otimes I$、$\alpha_b=\sigma_x\otimes\sigma_b$。正常空CAR真空约定固定，不作Dirac海近似。

完整二次量子化的每站质量算符为 $m_i=c_i^\dagger\beta M c_i$；正向边跳跃为 $-ic_i^\dagger\alpha_b c_{i+e_b}/(2\ell)+\mathrm{h.c.}$，每轴 $t_{bi}$ 将两邻边各半分。原公共尺度与新增无迹路径的来源为

$$
F_{si}=\epsilon m_i+\nu\sum_bt_{bi},\qquad
F_{ti}=g(t_{xi}-t_{yi}),\qquad
(\epsilon,\nu,g,\Omega,T)=(.15,.35,.2,.5,2).
\tag{1}
$$

每复制保原关系自旋 $G_s$，增加 $G_t$。定义

$$
\begin{aligned}
h_r&=H_{S,0,r}-\Omega\sum_i(X_{si,r}+X_{ti,r}),\\
W_{rs}&=\sum_{i,c=s,t}(Z_{ci,r}F^S_{ci,s}+F^S_{ci,r}Z_{ci,s}),\\
H_{BG,K}&=\sum_rh_r+\frac1K\sum_{r<s}W_{rs},\\
H_{\mathrm{act}}&=C_KI+H_{BG,K}+H_{P,0}
+\sum_{i,c}F^P_{ci}\overline Z_{ci},\qquad
\overline Z_{ci}=K^{-1}\sum_rZ_{ci,r}.
\end{aligned}
\tag{2}
$$

不含 $r=s$；背景有限复制的 $(K-1)/K$ 与探针的系数差在有限K界内。所有CAR耦合为偶算符，空间支持至多一条邻边；同站不同复制是内部自由度。式(2)是一份有限、自伴、自治Hamiltonian，完整总能量精确守恒。它没有宣称全局连续时间邻边Hamiltonian给出严格离散光锥或已实现所有有限控制装置。

完整Fock范数而非仅单粒子范数给

$$
\begin{aligned}
b_S&=2\epsilon m_S+6\nu/\ell+4g/\ell=12.38125,\\
b_P&=2\epsilon(.3+.5)+18\nu/\ell+12g/\ell=35.58375,\\
\|W\|&\le w=2Vb_S=6800401.5625,\qquad
\mathcal L\le Vb_P=9772187.34375.
\end{aligned}
\tag{3}
$$

理由是源每轴半边和范数不超过 $2/\ell$，三质量探针为 $6/\ell$；质量Fock范数分别为 $2m_S$ 和 $2(.3+.5)$。例如

$$
C_K=1+KV(2m_S+6/\ell+2\Omega+b_S)
+V[2(.3+.5)+18/\ell+b_P]
\tag{4}
$$

足以使全谱正。该常数也必须放入完整比较器、保相同相位。它不是现实真空能推导；关系自能尚未由普适几何方程选择。

## 2. 同能谱来源和原公共尺度核

每个格点动量令 $k_b=\sin(\ell p_b)/\ell$。来源A准备完整正常Slater态，协方差

$$
C_A(p)=
\begin{cases}
[I+\operatorname{sgn}(k_x)\alpha_x]/2,&k_x\ne0,\\
(I+\beta)/2,&k_x=0 .
\end{cases}
\tag{5}
$$

每动量两个模式占据。每份复制背景相同且已知，两关系自旋先取 $+X$；背景准备跨空间可以有关联，未由式(5)免费获得局域制备器。

来源B经 $x/y$ 交换、$z$ 不变的反射运输。令

$$
u=-i(\sigma_x-\sigma_y)/\sqrt2,\qquad U_D=\beta(I\otimes u).
\tag{6}
$$

直接乘法得 $U_D\beta U_D^\dagger=\beta$，$\alpha_x,\alpha_y$ 互换、$\alpha_z$ 保持。因此原四动量 $0,+e_x,+e_y,+e_z$ 及质量钟码闭合。此处使用给定Dirac模型的parity，不要求所有认知原理具有该对称。若改用交换x/y同时反转z的正旋转，原正轴码不闭合，不能混用。

将整个背景空间反射、旋量运输并共轭所有 $G_t$ 的 $X_t$，式(2)的 $H_{BG,K}$ 不变且两指定背景准备互换。因此它们具有相同完整背景总能量谱分布，所有演化时刻亦然。此断言不适用于“背景改变但探针输入固定”的整个实验总Hamiltonian谱分布。

**原尺度核。** 令 $g=0$。平移均匀Hartree演化只通过来源的 $E_m$ 和总 $E_t$ 写入原 $G_s$。式(6)的对称保持两标量，故A/B的 $z_s,N,v$ 历史及同一探针整个准备等距相同。背景完整态本来可以不同；读取背景可以辨认来源。本命题仅指原探针/R终端菜单，不排除全部局域标量理论。

## 3. 自治协方差：同一物质来源推动两条公共关系

取 $\Gamma=(\beta,\alpha_x,\alpha_y,\alpha_z)$，$L_{ab}=i\Gamma_a\Gamma_b$，$s_{ab}=-s_{ba}$，写

$$
C_p=\frac12\left[I+\sum_ar_a\Gamma_a+\sum_{a<b}s_{ab}L_{ab}\right],
\qquad h_p=\sum_ad_a\Gamma_a .
\tag{7}
$$

由 $\dot C=-i[h,C]$ 精确得到

$$
\dot r_a=-2\sum_bs_{ab}d_b,\qquad
\dot s_{ab}=2(d_b r_a-d_a r_b).
\tag{8}
$$

这是四向量加六双向量的闭合，不丢两个频率支。正负动量反射及相应Clifford运输将 $65^3$ 模精确归并成 $33^3$ 类，保全部权重；轴置换不再是可用的6545类归并。

每份Hartree cluster仍分成源Slater与各站两种关系自旋：其有效Hamiltonian中来源和关系项分离。用各自Bloch向量 $\mathbf g_s,\mathbf g_t$，

$$
\begin{aligned}
N&=1+\epsilon z_s,\quad \bar v=1+\nu z_s,\\
v_x&=\bar v+gz_t,\quad v_y=\bar v-gz_t,\quad v_z=\bar v,\\
d&=(Nm_S,v_xk_x,v_yk_y,v_zk_z),\\
E_m&=2m_S\langle r_0\rangle,\quad E_{tb}=2\langle k_br_b\rangle,\\
b_s&=\epsilon E_m+\nu\sum_bE_{tb},\quad b_t=g(E_{tx}-E_{ty}),\\
\dot{\mathbf g}_c&=2(-\Omega,0,b_c)\times\mathbf g_c .
\end{aligned}
\tag{9}
$$

来源方向确实写入自治关系后态，$z_t(t)$ 不是外给函数。每站总平均能量

$$
NE_m+\sum_bv_bE_{tb}-\Omega(x_s+x_t)
\tag{10}
$$

守恒。A/B的 $z_s,E_m,\sum E_t$ 相同，$z_t$ 反号、$E_{tx},E_{ty}$ 互换；两个背景各自保完整资源。

新 $F_t$ 不是任意命名来源。共同连续比较器在暂时独立的 $N,a_b$ 下有

$$
h=N\beta M+\sum_b(N/a_b)\alpha_bp_b,\qquad a_b=N/v_b .
\tag{11}
$$

沿 $(z_s,z_t)$ 拉回，导数分别为 $\epsilon\beta M+\nu\sum\alpha_bp_b$ 与 $g(\alpha_xp_x-\alpha_yp_y)$；格点父过程正是其完整二次量子化。只匹配这两条路径，未生成任意局部四维应力或关系自能的GR自耦。

“无迹”在此精确指速度系数的方向 $\delta v=(g,-g,0)$，在等方初态也是空间度规的无迹切向。有限 $z_t$ 下 $a_b=N/v_b$ 非线性包含高阶标量/体积变化，不声称整条路径处处是纯无迹度规变分；A/B的标量和体积仍因两轴互换而完全相同。

## 4. 原钟、实际三维采样和非共形比较器

探针准备码是四动量乘Dirac4乘三质量，总48维；所有幅度、跨质量/动量相干及任意参考R允许未知。共同生成元为

$$
h_{P,p}=N\beta M_P+\sum_bv_b\alpha_bk_b .
\tag{12}
$$

终端仪器保留全部实际站点细记录和内部钟结果。有限汇总 $f_b(j)=[1+\sin(2\pi j_b/65)]/2$，移动钟使用后两质量标签的 $P_{-Y}$。固定见证为四模式等幅、$\chi=(1,0,0,1)/\sqrt2$、钟 $|+\rangle$；质量零实验在同一父材料的零质量标签，非更换Hamiltonian。报告所有三轴，不按最大特征态重新选输入。

连续低带字典是对角Bianchi I

$$
ds^2=N^2dt^2-\sum_ba_b^2dx_b^2,\qquad a_b=N/v_b.
\tag{13}
$$

用场重标 $\psi=(a_xa_ya_z)^{1/2}\Psi$，共同概率内积和时间自旋联络的体积项吸收后，Dirac生成元为式(11)。这是同一完整质量和三轴导数的成熟几何比较，不是由一个末钟读数层析完整曲线。离散位置PVM通过完整有限Fourier映射运输为采样仪器，不能叫连续精确点位置PVM；声明的平滑效果在原码内精确匹配，细记录不被删掉。

初时全部尺度为1、导数为0。式(5)给

$$
E_m(0)=4/65,\quad E_{tx}(0)=Q=\frac2L\cot\frac{\pi}{130}>0,\quad E_{ty}(0)=E_{tz}(0)=0,
$$

$$
z_s''=-4\Omega(\epsilon E_m+\nu Q),\qquad
z_t''=-4\Omega gQ,\qquad
(a_x'',a_y'',a_z'')=(u-r,u+r,u),
\quad u=(\epsilon-\nu)z_s'',\ r=gz_t'' .
\tag{14}
$$

因初始 $N=1,\dot a_b=0$，坐标时间与正交固有时间的二阶尺度导数相同。此时电Weyl的平方给

$$
C_{\mu\nu\rho\sigma}C^{\mu\nu\rho\sigma}
=2\sum_b(a_b''-\overline{a''})^2=4r^2>0.
\tag{15}
$$

故这一**已声明比较器**不属于任何共形平直族，不能由共同坐标换表变成原FLRW。然而有限钟/站点菜单并没有层析全部Weyl或排除一切不同模型的观测等价；式(15)是比较器解析性质，实际有限证书是下一节的方向留出。

## 5. 全未知输入和精确背景的有限K桥

复用[有限副本证明](../_admission/after1174_direction/mean_field_bound.md)，在新的 $h,W,F$ 上使用式(3)。对每份固定已知背景准备及整个未知探针/R，

$$
\left\|U_{\rm act}(t)E-
\bigl[U_{BG,K}(t)\phi^{\otimes K}\bigr]\otimes U_P(t)\right\|
\le \eta_K=\frac{4\mathcal LT e^{6wT}}{\sqrt K}.
\tag{16}
$$

两边保同一 $C_K$ 相位；比较背景是精确BG，不是整体乘积态。选 $K=10^{70881071}$ 足以令每支 $\eta_K\le10^{-5}$。这是巨大但有限的数学资源，不称已模拟或现实可行。后接同一完整终端记录CP不增加误差；不允许中途改变BG后继续免费沿原轨道。

式(16)中的 $E$ 是固定背景乘任意探针输入的准备等距，右侧同样按输入向量作用，算符范数对全48码及任意R统一。实际记录的分布式制备、当地采样和后来的因果汇集仍由声明的有限资源合同承担，不授予瞬时全局读取或免费重置活钟。

在 $g=0$ 的原尺度核中，A/B探针/R半diamond距离不超过两支 $\eta_K$ 之和。整个背景不相同且仍保留；这个量词是结论的一部分。

## 6. 同一A/B差的相关格点界和数值责任

普通两支独立格点界会覆盖固定钟差。使用已独立证明的[相关通道混合导数界](../_admission/after1175_selection/correlated_lattice_bound.md)：固定同一真实Hartree背景、同一未知码态和同一效果，

$$
|(\Delta p)_{\rm cont}-(\Delta p)_{\rm lat}|
\le 2dJ+4pdJV,\qquad
p=\pi/8,\quad 0\le d=p-\sin(\ell p)/\ell\le d_*=\ell^2p^3/6,
\tag{17}
$$

$$
J=\int_0^T|gz_t|dt,\qquad V=\int_0^T(\bar v+|gz_t|)dt .
$$

证明使用参数矩形 $p\mapsto p-\lambda d$、$gz_t\mapsto ugz_t$ 的两次Duhamel微分；两种时间排序合成完整矩形，最后四通道差迹为零，转概率时可取half-diamond。质量项保留，直和码不退相干。该界只收紧A/B差，单支误差没有缩小。也不离散近似来源的高动量模式。

新十分量残差不能沿用1175旧数值。独立[数值包络](../_admission/after1175_selection/tensor_candidate/numerical_bound.py)用2048步RK轨道、逐步cubic Hermite残差、六次Bernstein凸包及显式基本算术包络认证。反对称旋转项消去，取

$$
E^2=|e_{g_s}|^2+|e_{g_t}|^2+2\langle|e_q|^2\rangle,\quad
E'\le2\sqrt2F_{\rm op}E+R,\quad
F_{\rm op}^2=\lambda_{\max}\langle f_cf_{c'}\rangle,
\tag{18}
$$

其中 $f_s=(\epsilon m_S,\nu\mathbf k)$、$f_t=(0,gk_x,-gk_y,0)$。完整立方权重下交叉Gram为0，$F_{\rm op}=F_s$。探针误差由同一背景积分界和自身多项式残差估计，保全部48码及R；用右极分解处理非严格酉的浮点传播矩阵，终端算术也支付预算。说明与原始包络保存在同目录，默认检查不冒充重新执行此积分。

[数值说明](../_admission/after1175_selection/tensor_candidate/numerical_bound.md)公开残差和基本算术假设；[独立极分解后验补证](../_admission/after1175_selection/tensor_candidate/polar_endpoint_guard.md)重放同一2048步终点，将存储矩阵逐元素转为精确有理数，不假定特征值程序精度。存储矩阵、数学精确极分解、原见证及解析Fourier效果的单概率责任小于 $2.961\times10^{-11}$，已在每支 $10^{-10}$ 末端预算内。报告中的SVD效果距离只作诊断，不参与概率下界。

固定证书给两支数值和小于 $9\times10^{-5}$；$J,V$ 使用绝对值积分和真实轨道误差，不靠采样点不变号。仪器合同每支 $\eta_{\rm inst}=2\times10^{-5}$ 是额外输入，覆盖准备、全部读出记录及允许后态；不声称已制造仪器。

因此本轮采用总预算

$$
B=2d_*J_{\rm up}+4p_{\rm up}d_*J_{\rm up}V_{\rm up}
+2\times10^{-5}+4\times10^{-5}+9\times10^{-5}.
\tag{19}
$$

前三类分别为空间A/B差、有限K两支、仪器两支、数值两支。舍入均向外。结果中保留精确数值，固定 $y$ 移动钟联合及质量零 $x$ 差减去 $B$ 后仍正。此保守共同预算同时支持实际父任务可分和声明连续比较器的方向响应，不签任意历史/高能/曲率仪器。

## 7. 控制、去重及结论边界

- $g=0$：原规则加旁观自旋，已证明A/B整个Hartree探针等距相同；有限父过程只在式(16)范围内相同。来源初态不同于1175旧 $\beta$ 填充，不复用旧数值轨道。
- 共同反射：在全部48码上，式(6)连同动量置换交织A/B生成元；将输入、所有端口、钟和细记录一起运输时概率不变。只换来源而固定原任务是不同物理实验。
- 质量零和有质量钟由同一固定父H的不同准备承担。原固定 $\chi$ 有方向性，零质量的 $y,z$ 定向均值为零不应隐去或换态优化。
- 1167已有公共尺度核，1175已有主动共同尺度；本轮净增是同一材料的真实无迹来源方向、自治承担者和原共同读口接续。形式投影、一般HDA、PF分类和GR结论均未重证。

成功不意味着 $g,\Omega$、允许方向、关系自身作用或其量子约束已被认知原则唯一确定。没有取得普适自旋2、Lorentz全物种、非线性Einstein或暗部门经验解释。
