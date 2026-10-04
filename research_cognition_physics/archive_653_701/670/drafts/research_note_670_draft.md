# 第670轮：非平坦规范场中的共同观测，以及无质量零模处的正则延伸

日期：2026-10-02。接[669](../../research_note_669.md)、[已执行入口](nonflat_gauge_entry.md)。[代码](../joint_nonflat_mass_measure.py)、[结果](../joint_nonflat_mass_measure_results.json)、[核验](../research_round_670_checks.json)、[全条件账](../unified_physics_condition_ledger_670.md)。两组新核验、十八式；主代理审查，无独立代理审查。

## 1. 接续、去重及结论层次

本轮读取最新导航、669报告／结果、670入口，并回查660—661、668、628及667的历史空间接口。当前工作仍是整合共同条件，认知系统设计后置，研究目标不变。

669已经证明：只积分原玻色变量而保留自由费米传播，不能恢复原局部Gauss。670入口已将同一质量和观测接到真实非零超荷曲率。本轮新增两项：

1. 以全部原颜色、弱、超荷表示和非交换链路核共同测度、物理二点函数及质量来源，明确一般代数合同。
2. **消去660观测变换中的逆无质量传播块。** 原配对与手征投影自身就给正则公式；无质量零模不再是这个字典的必需排除条件。

|层次|本轮地位|
|---|---|
|认知动机|表示改变应保留同一可测过程，不因某种坐标写法奇异而误判过程不存在|
|继承的建模输入|598原质量、614全左手字典、615原辅助配对、660局部候选、668质量拉回；有限给定格点和overlap核|
|解析结论|原配对下无逆K的三角／观测公式，完整带来源多项式身份和局部规范协变|
|数值验证|非交换非平坦原SM背景；原holonomy中无质量零模与全部原非零质量共存|
|未证明|动态规范反射正性、所有背景严格归一、原完整Hamiltonian同一性、连续极限及量子GR|

### 1.1 旧空间结论直接继承

按用户提醒，复用[382](../../../archive_370_428/research_note_382.md)完整反向接口上界、[383](../../../archive_370_428/research_note_383.md)自由连续对合、[384](../../../archive_370_428/research_note_384.md)已经消去额外Lipschitz的下界、[386](../../../archive_370_428/research_note_386.md)真实邻域桥与有限尺度证书、[425](../../../archive_370_428/research_note_425.md)半幅／成本收缩的替代坐标桥，以及[522](../../../archive_467_530/research_note_522.md)—[523](../../../archive_467_530/research_note_523.md)热参考与实际方向instrument的共同实现。具体输入对照保持[667 §1.1](../../research_note_667.md)和全条件账C07—C09、C20的记录。

386与425不叠加为必要条件，384已消去的要求不重新引入，523的共同存在性不重证。本轮未修改这些定理，也未新增空间假设；当前仍需把旧实际端点、参考、仪器和尺度映射识别为同一物质过程的对象。下面的正则性是有限费米字典的正则性，**不是**空间光滑流形的另一证明。

## 2. 原核给出的两种投影

在同一有限、index-zero、Wilson核H不穿零的背景片，u、v分别为H的负、正谱正交帧。采用660记号，固定bar手征帧J₋、J₊以及其投影Q₋、Q₊。有

$$
P_u=uu^\dagger,\quad P_v=vv^\dagger,\quad P_u+P_v=I,
\quad Q_\pm=J_\pm J_\pm^\dagger,
\quad D=\tfrac12(I+\gamma_5\operatorname{sign}H).
\tag{1}
$$

直接作用于H的正负谱子空间，而非使用自由动量，得到

$$
Du=Q_-u,\qquad Dv=Q_+v,\qquad
K_\ell=J_+^\dagger Dv=J_+^\dagger v,\qquad K=J_-^\dagger u.
\tag{2}
$$

注意H有隙不蕴涵Kₗ可逆。原S⁹配对M为逐点B⊗T(E)（一站holonomy按内部优先次序写为T(E)⊗B）；这些只是同一矩阵的排列。由原B、T可直接核

$$
M^{\mathsf T}=-M,\qquad MQ_+=Q_+^{\mathsf T}M,
\qquad M^\dagger M=I\quad(|E_x|=1).
\tag{3}
$$

这不是向一般配对模型添加新条件，而是615／660已使用配对的恒等式。若另换不保此手征分块的配对，下面消去式不能自动复用。

## 3. 消去旧三角字典中的逆K

660定义C=uᵀMQ₊v、F=vᵀMQ₊v。式(2)—(3)逐项给

$$
C=u^{\mathsf T}MJ_+K_\ell,
\quad C^{\mathsf T}=-K_\ell^{\mathsf T}J_+^{\mathsf T}Mu,
\quad F=K_\ell^{\mathsf T}J_+^{\mathsf T}Mv.
\tag{4}
$$

因此在旧逆式可用处，其两块恰好约掉；更一般地直接定义

$$
d'=d+\bar B^{-1}Kb,\qquad
e'=e-J_+^{\mathsf T}Mu\,b-\tfrac12J_+^{\mathsf T}Mv\,c,
\quad b'=b,\ c'=c.
\tag{5}
$$

原bar辅助配对的\(\bar B\)为可逆单位配对，Pfaffian为1；这与Kₗ无关。式(5)仍为det T=1的三角变换，**包括Kₗ奇异处**。直接代入原660二次型，各交叉项按式(4)相消，得

$$
T^{-\mathsf T}S^{\mathsf T}\mathcal N S T^{-1}=
\begin{pmatrix}
A&0&0&0\\0&0&0&-K_\ell^{\mathsf T}\\
0&0&\bar B&0\\0&K_\ell&0&0
\end{pmatrix},\qquad A=u^{\mathsf T}Mu.
\tag{6}
$$

无需诉诸“可逆点稠密”的额外条件。原物理Y=(w,e′)中w=J₋†vc，转回ξ=(ψ,barψᵀ)，得到**不依赖帧相位、不含逆K的实际观测映射**：

$$
Y=\Phi\xi,\qquad
\Phi=
\begin{pmatrix}
J_-^\dagger P_v&0\\
-J_+^{\mathsf T}M(P_u+\tfrac12P_v)&J_+^{\mathsf T}
\end{pmatrix}.
\tag{7}
$$

投影的范数不超过1，\(P_u+P_v/2\)的范数为1。对任意两分量输入，以a、b表示其范数，有输出范数平方≤a²+(a+b)²。因此

$$
\|\Phi\|\le\frac{1+\sqrt5}{2},\qquad
\|\lambda\Phi^{\mathsf T}P_\phi\Phi\|
\le |\lambda|\frac{3+\sqrt5}{2}\|P_\phi\|.
\tag{8}
$$

这是算子范数界，不是结果表中随模式数增长的Frobenius范数。界不含最小奇异值σmin(Kₗ)。H不穿零时投影随背景光滑，故Φ也光滑；接近H的谱隙关闭处，其导数仍可能失控。本轮没有移除H及定秩条件。

**对660的精确补充：**660已延伸无质量权重，并指出其含逆K的观测写法不能直接延伸。本轮利用原配对恒等式完成观测本身的延伸；这是旧限制的缩减，不覆盖或改写被冻结报告。661指出的时间支撑问题仍在：谱投影一般跨越两个时间半区，光滑有界不等于严格局部。

## 4. 同一质量和全部物理来源

使用668原Dirac／Majorana质量Pφ，在同一ξ上定义

$$
\mathcal N_\lambda=\mathcal N+\lambda\Phi^{\mathsf T}P_\phi\Phi,
\quad L=\operatorname{diag}(J_-^\dagger v,I),\quad
N_{W,\lambda}=\begin{pmatrix}0&-K_\ell^{\mathsf T}\\K_\ell&0\end{pmatrix}
+\lambda L^{\mathsf T}P_\phi L.
\tag{9}
$$

取独立Grassmann来源j。三角合同、bar配对Pf=1和原方向约定给

$$
\int d\xi\ e^{\frac12\xi^{\mathsf T}\mathcal N_\lambda\xi+j^{\mathsf T}\Phi\xi}
=\frac{\operatorname{Pf}A}{\det S}
\int d\eta_W\ e^{\frac12\eta_W^{\mathsf T}N_{W,\lambda}\eta_W+j^{\mathsf T}L\eta_W}.
\tag{10}
$$

两边均为有限Grassmann多项式；未要求A、Kₗ或N_W,λ可逆。故包括零模处的所有**未归一物理来源系数**，不只一个无来源权重。对原S⁹积分仍成立。无来源为

$$
\operatorname{Pf}\mathcal N_\lambda=
\frac{\operatorname{Pf}A\operatorname{Pf}N_{W,\lambda}}{\det S},
\qquad
G_Y=-L N_{W,\lambda}^{-1}L^{\mathsf T}\quad
\text{仅在该归一高斯核可逆时}.
\tag{11}
$$

零权重时不使用右侧逆核；相应插入由Pfaffian余子式给出。质量或背景来源也由同一多项式求导；在非零片可以写为

$$
\partial_s\log Z=
\partial_s\log\!\left(\frac{\operatorname{Pf}A}{\det S}\right)
+\tfrac12\operatorname{tr}\left(N_{W,\lambda}^{-1}\partial_sN_{W,\lambda}\right).
\tag{12}
$$

这里Z指固定背景、固定E权重；整体积分后需对整体Z求导，不互换两个log。式(9)含L及Pφ的共同导数；局部表示同样含dΦ，不能冻结观测映射。

## 5. 非平坦规范协变与帧选择

对原SM商群的真实局部表示R，链路为U′xy=RxUxyRy†，投影随R共轭，M′=R*MR†。令r为物理两分量手征表示、Gξ=diag(R,R*)、GY=diag(r,r*)。原质量同步变换，则

$$
\Phi'G_\xi=G_Y\Phi,\qquad
G_Y^{\mathsf T}P_{\phi'}G_Y=P_\phi,\qquad
G_\xi^{\mathsf T}\mathcal N_\lambda'G_\xi=\mathcal N_\lambda.
\tag{13}
$$

这是式(7)的代数结果，不需自由谱，也不需Kₗ可逆。det Gξ=1保证完整带相位Pfaffian不变。u→uVᵤ、v→vVᵥ的独立帧选择还给

$$
\operatorname{Pf}A\mapsto\det V_u\operatorname{Pf}A,
\quad\operatorname{Pf}N_{W,\lambda}\mapsto\det V_v\operatorname{Pf}N_{W,\lambda},
\quad\det S\mapsto\det V_u\det V_v\det S.
\tag{14}
$$

整体相消。局部帧片中的这项身份不是全拓扑部门的光滑测度或相互作用正性定理。原615借用的[SO(10)辅助测度文献](https://arxiv.org/abs/1710.11618)有自己的全拓扑定义及局域性问题；本轮原质量／观测候选的范围须独立验收，不能只靠引用升级。

第一组使用2空间×2时间盒、原完整16维内部表示，两条非交换链路的plaquette满足

$$
\|U_\square-I\|_F=0.206418365813,
\quad\|U_aU_b-U_bU_a\|_F=0.0110756718210,
\quad\min|\sigma(H)|=0.936248437462.
\tag{15}
$$

分别在四个时空点实施原SU(3)、SU(2)、U(1)变换，同步改变E、φ和链路。原权重、物理二点核及质量来源协变；最大误差<2.6×10⁻¹⁴。帧相位独立改变时仍相消；正则局部式与物理式的权重／二点核也相同。比较完整复数，不把绝对值当测度。

## 6. 原无质量零模处的实际见证

复用615一站holonomy，仅作为代数边界检验，不把平坦holonomy重算成非平坦证据。原H隙恒为1；取θ₀=π/6，超荷q=6通道使Kₗ有精确零模。保原φ=PHI[0]、全部原Y和λ=.37，得到

$$
\sigma_{\min}(K_\ell)=0\ \text{解析上},\qquad
\sigma_{\min}(N_{W,\lambda})\simeq0.00129614093723,
\quad\operatorname{Pf}\mathcal N_\lambda\simeq2.05917074181\times10^{-10}.
\tag{16}
$$

有限精度Kₗ最小奇异值为1.23×10⁻¹⁶。零点处三角身份误差2.23×10⁻¹⁶，完整权重相对误差2.45×10⁻¹⁵，归一二点核绝对误差2.88×10⁻¹³。四组明确余子式与未归一二点系数也吻合。源响应跨该点有限：

$$
\partial_\theta\log Z\simeq-4654.46764291,
\quad |\partial_\theta\log Z-\text{中心差分}|<8.8\times10^{-6},
\quad |\partial_\lambda\log Z-\text{中心差分}|<3.1\times10^{-10}.
\tag{17}
$$

响应较大，不能将“有限”说成“对所有背景一致稳定”。此处存在一个可归一例和其连续小邻域；没有证明所有背景或指定动态积分均非零。

### 首次失败检查及修正

开发时曾预期含逆K的旧观测矩阵随零模逼近而发散，增长断言失败。实际Frobenius范数约7.70→7.675。已保留[首次代码](first_attempt_mass_measure.py)和[诊断](first_attempt_diagnostic.json)。没有放宽断言凑通过，而是由式(2)—(4)找到确切相消，改验正则身份。该失败否定一个坐标发散预期，不是模型失败或认知原则反证。

## 7. 真正关闭的条件及下一轮

本轮关闭的是原配对下“共同观测／质量字典必须排除无质量Kₗ零模”的要求；非平坦、非交换规范链路与原质量和来源共同相容。保留清楚的边界：

$$
\text{有界正则观测＋规范协变＋完整来源身份}
\ \not\Rightarrow\ 
\text{动态规范反射正性＋严格归一＋原实时间过程同一性}.
\tag{18}
$$

[成熟自由overlap反射正性结果](https://arxiv.org/abs/1005.3751)证明自由及所述非规范相互作用的情形，不能直接充当本轮一般动态规范正性证明。605—610原热态全支撑、变化投影和电动能限制继续有效；643原Gauss有序历史不重建为新成果；221—222、230、506的完整过程合同不降格成相同权重。

接[671入口](../../671/drafts/STATUS.md)：用同一正则字典检查带真实空间曲率、反射兼容背景上的物理半区泛函，再定位动态规范积分所需的正时间组合合同。禁止用任意不反射对称背景的Gram非厄米直接反证总体正性，也不能用有限二点矩阵正性代替全部多项式和动态测度。四分支尚未统一，认知设计后置，目标不改。
