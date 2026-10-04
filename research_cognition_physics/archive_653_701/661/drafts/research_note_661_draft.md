# 第661轮：原物理Weyl观测与辅助测度的共同反射泛函

日期：2026-10-02。接[660](../../research_note_660.md)与[实际入口](observable_interface_entry.md)。[代码](../joint_physical_auxiliary_state.py)、[结果](../joint_physical_auxiliary_state_results.json)、[核验](../research_round_661_checks.json)、[条件账](../unified_physics_condition_ledger_661.md)。

## 1. 本轮合并的对象与边界

660已经精确匹配原权重、变换后的观测和来源，但观测变换是否保持正时间支撑尚未解决。本轮完成两件相互关联的核对：

1. 原物理Weyl局部观测若通过660变换实现，确实跨越时间半区；直接使用候选裸局部分量也不保持原条件关联。这两个具体接法不能承担原局部观测的正性转移。
2. 对**原保留物理Weyl代数与原E辅助代数**，可以直接从原积分建立同一背景下的归一反射泛函，不必经过该非局部变量映射。两部门的原来源必须共同变化，分别固定它们的背景会改变实际联合响应。

这是自由有限格、指定原观测代数的真实连接。乘积保持正性的数学定理本身是既有工具，不作为新发现。没有改换原辅助核，没有引入新的认知装置，也没有将所有镜像场、原标量记录或动态引力纳入已证范围。

上一轮完成660并保存本轮实际探查，属于有效进展。本轮先读导航、660报告及结果、完整659条件账和660增量，回查612、624—625、647、649—650、652—658。入口无运行中的Python进程。初始读取657代码时旧文件名不正确，检索后使用实际的joint_auxiliary_reflection_gluing.py；没有由不存在的文件推断结果。

## 2. 原积分，而非另配两种状态

沿原自由单位规范背景、m₀=1、周期空间及2L格反周期时间；原D、16内部T、B与均匀S⁹处方不变。四维格、内部表示和零物理标量依旧是输入。将物理Weyl投影记P_v=vv†，避免与镜像u投影混淆：

$$
\widehat\gamma_5=\gamma_5(1-2D),\quad
P_v=\frac{I-\widehat\gamma_5}{2}=Q_-+\gamma_5D,\quad
w=J_-^\dagger vc,\quad \bar w^{\mathsf T}=e,\quad K_\ell=J_+^\dagger Dv.
\tag{1}
$$

w是原左手场的负γ₅自旋分量，e是原左手bar场的正γ₅分量。时间“局部”按这些原Weyl场分量的位置定义，不把手征帧系数c本身当格点局部字段。

原饱和镜像配对对b积分得Pf A_D(E)，原独立bar-E配对积分为1。对不含b和镜像bar字段的保留观测，原积分因此精确为

$$
\mathcal I_D(Ff)=C_D\!int d\mu(E)\operatorname{Pf}A_D(E)f(E)
\int dc\,de\,e^{e^{\mathsf T}K_\ell c}F(w,e),
\tag{2}
$$

C_D固定所有原基方向和Grassmann次序；不可在背景求导时随意重设。式(2)也可由660的式(6)得到，但本轮先把它当原变量的积分身份。归一后

$$
\Omega_D(Ff)=\omega_{{\rm W},D}(F)\,\omega_{{\rm aux},D}(f),\qquad
\mathcal A_{+,D}=\mathcal A^{\rm W}_{+,D}\otimes\mathcal A^{E}_+.
\tag{3}
$$

所有部门共用D、格、时反射和切口。代数取Weyl局部分量的有限Grassmann多项式与E半区连续函数的有限张量和。有限格和紧球面使原积分定义良好；不将欧氏Grassmann变量当作经典随机变量。

## 3. 对照成熟Weyl代数，定位660映射的实际限制

[Kikukawa—Usui §IV式(46)—(50)](https://arxiv.org/html/1005.3751v3#S4)证明的是指定Weyl自旋分量代数，不是全部投影后分量。其反射交换w与bar-w，并反线性地反转乘积次序。原核满足

$$
J_-^\dagger vK_\ell^{-1}
=J_-^\dagger D^{-1}J_+,\qquad
G_{w\bar w}=-J_-^\dagger D^{-1}J_+.
\tag{4}
$$

负号按原指数+barψDψ的有序Wick约定。证明从P_vD⁻¹=Q−D⁻¹+γ₅出发，两边选负／正自旋块后接触项消失；Wick展开将完整多项式泛函与成熟自由Dirac块相连。612已经使用此逆核身份，本轮不重算为新定理。

[入口代码](observable_interface_probe.py)在原2空间×4时间盒核对式(4)，然后比较660候选的裸局部分量、准确变换后的分量及时间支撑：

$$
\|G_{\rm bare}-G_{\rm original}\|_F=19.5362350126,\quad
\|G_{\rm mapped}-G_{\rm original}\|_{\max}<2.2\times10^{-15},\quad
\|\Phi_{+,-}\|_F=3.40879388268>0.
\tag{5}
$$

最后一项是原正半区观测到候选负半区变量的映射块。因而本次明确排除“660映射保时间支撑”，并未排除其它正实现。前两项是固定E条件核，不能说无E来源的球面平均二点函数也必有相同差异。

支撑失败还有无需数值阈值的解析见证：在原单空间格、两时间格，D₁₀=−γ₄/2；变换后w=J−†P_vψ=J−†(I−D)ψ，所以正区t=1的w对负区t=0的ψ系数为J−†γ₄/2，其算符范数精确为1/2。这个已经非零的子块足以否定整套变换保正时间支撑，不需依赖式(5)的小数。

## 4. 原物理代数的共同正内积

取正半区t=L,…,2L−1和同一反射θt=2L−1−t。657按另一半命名其正区；对标量E交换两半使Gram核转为其共轭，正性不变，因此可与这里的Weyl反射使用同一切口。658给原辅助积分严格非零，成熟Weyl自由归一也非零。

对任意有限和A=Σ_i F_i f_i，由式(3)得到

$$
\Omega_D(\Theta A A)=\sum_{ij}W_{ij}E_{ij},\qquad
W_{ij}=\omega_{{\rm W},D}(\Theta F_iF_j),\quad
E_{ij}=\omega_{{\rm aux},D}(\Theta f_i f_j),\quad
W\succeq0,\ E\succeq0,\ W\circ E\succeq0.
\tag{6}
$$

将求和理解为全1向量对W∘E的二次型即可非负；一般复系数可以并入F_i。E函数为偶，故不引入额外费米交换符号。这是原物理代数的直接反射证明，绕过了式(5)失败的变量映射。

取零范数商并完备化，还得到这一个有限盒的等距识别

$$
\mathcal H_{{\rm OS},D}^{\rm joint}\simeq
\mathcal H_{{\rm OS},D}^{\rm W}\otimes\mathcal H_{{\rm OS},D}^{\rm aux},\qquad
[Ff]\longmapsto[F]\otimes[f],\qquad
[1]\longmapsto[1]\otimes[1].
\tag{7}
$$

这里OS空间由反射配对构造。尚未证明各时间盒的状态、平移或区域代数有统一极限；不能据此宣称共同Hamiltonian或原正常Gibbs态已经找到。也没有把E定义为624的原物理标量s。

## 5. 用实际完整球面权重检验联合接口

单空间格、两个时间格仅作可完全求积的交叉检查，不改变上节的一般自由有限格论证。653已给原闭链权重k(E,F)k(F,E)=((1+E·F)/2)¹⁶。于是

$$
Z_{\rm aux}=\frac{(9/2)_{16}}{(9)_{16}}=\frac{6486285}{2147483648},\qquad
\omega_{\rm aux}(E^a_+E^b_-)=\frac{\delta^{ab}}{10}\frac{16}{25},\qquad
G_E\big|_{\{1,\sqrt{10}E^a\}}=\operatorname{diag}(1,\tfrac{16}{25}I_{10}).
\tag{8}
$$

球面对积t的权为(1−t²)⁷ᐟ²。代码以10点Gauss–Jacobi对原64维Pfaffian精确径向求积；旋转协变性给全部向量分量，并非将原球面替换为圆。原Weyl正半区64个线性分量给单位Gram，联合704维Gram最小本征值为16/25。

另在2空间×4时间格保留全部内部矩阵，四组半区配置的原辅助核与四个Weyl波包的反射矩阵共用同一D；联合Gram最小本征值0.01938717028。这个空间核检查是有限样本，不是该大球面配置空间的完整积分。

为避免只对手写乘积自证，代码将一个原反射双线性观测输送进660完整1024维Nambu矩阵，直接扰动Pfaffian作来源插入。其复关联为0.01114798507−0.02084096414i，与原Weyl核差4.2×10⁻¹²，真实相位保留。全代数正性依赖式(2)—(6)解析证明，不依赖这几个样本。

## 6. 共同背景与来源不因因子化而独立

原未归一对象在任一光滑无零点帧片为

$$
Z(\lambda)=\frac{\det K_\ell(\lambda)}{\det S(\lambda)}
\int d\mu(E)\operatorname{Pf}A_{D(\lambda)}(E),\qquad
\partial_\lambda\log Z=\partial_\lambda\log\det K_\ell-
\partial_\lambda\log\det S+\partial_\lambda\log Z_A.
\tag{9}
$$

660已经核帧相位的共同抵消；本轮向真正联合观测推进。采用615原holonomy θ，r(E)=Σ_{a=0}⁵ E_a²，r服从Beta(3,2)。在原16表示选择实际存在的电荷q=−4，定义固定局部Weyl双线性O，使其条件均值为tan(qθ/2)。则

$$
\Omega_\theta(Or)=o(\theta)\mu(\theta),\qquad
o=\tan(-2\theta),\quad
\mu=\frac{\int_0^1 12r^3(1-r)[a r+b(1-r)]^8dr}{Z_A},\quad
a=\cos^2\theta,\quad b=\cos^2(3\theta/2).
\tag{10}
$$

O按代码中i·tr(G_wbar U)/2定义，U=VP†γ₄VM；它是所选两个自旋分量的双线性，非新增测量设备。这个带电一格背景用于来源代数核验，不扩张式(6)的自由单位规范正性范围。

记未归一r矩为m，原辅助积分为z。完整条件归一要求

$$
\mu'=\frac{m'}z-\mu\frac{z'}z,\qquad
\mu''=\frac{m''}z-\mu\frac{z''}z-2\mu'\frac{z'}z.
\tag{11}
$$

于是同一θ的联合响应是

$$
\partial_\theta\Omega_\theta(Or)=o'\mu+o\mu',\qquad
\partial_\theta^2\Omega_\theta(Or)=o''\mu+2o'\mu'+o\mu''.
\tag{12}
$$

代码以7点Gauss–Legendre精确积分原权重及r来源，再从完整Nambu逆核和原Pfaffian计算联合观测；这一径向求积次数覆盖原多项式及Beta密度。θ=.23时

$$
\Omega(Or)=-0.30782086419309,\qquad
\partial_\theta\Omega(Or)=-1.63973680134112,\qquad
o\mu'=-0.09212374253905\ne0.
\tag{13}
$$

完整一／二阶响应与有限差分分别相差5.1×10⁻⁸、1.6×10⁻⁷。漏掉辅助部门的共同来源会遗漏最后这个非零项。这里没有新发现一般乘积求导律，新增的是原带电测度、物理Weyl观测及同一来源的实际匹配与定量见证。

## 7. 明确返回其它未整合条件

本轮关闭自由有限盒上“原物理Weyl局部观测与原E测度能否共享归一反射泛函”这个接口。660全部镜像观测的非局部映射不必继续作为此接口的前置条件。剩余工作须针对真实共同对象：

$$
\text{有限自由反射泛函}\ \not\Rightarrow\
\left\{\text{原完整H／Gibbs},\ \text{原s记录instrument},\
\text{动态几何与量子约束},\ \text{共同物理连续极限}\right\}.
\tag{14}
$$

624的读口是sqrt(1/2±sin s/4)，647的关系定位使用原h、s及其梯度，652给同一原H下的参考算符；它们不能由“有一个E场”替代。612的全时谱限制继续有效。没有宣称一般规范反射正性、标准模型或量子GR完成。

两组新核验覆盖原联合反射对象及共同来源；入口的映射反例单独复算但不增计组数。代码初稿曾按5维向量电荷误选q=−2，断言前索引失败；读取实际16维电荷后改为其真实q=−4，再完整运行通过。这是代码开发修正，不是物理排除。

接[662同一来源、关系区域与引力约束](../../662/drafts/STATUS.md)。先回查585局部lapse缺口、647关系定位、649—651区域配对以及652原量子参考，检验共同局部来源能否承担约束接口。不要继续优化镜像装置或把本轮子代数当成全目标；条件账保留全部四个模型分支。
