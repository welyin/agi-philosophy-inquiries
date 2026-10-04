# 第703轮：共同几何来源、热响应与原记录的同一极限

日期：2026-10-02。接[702](research_note_702.md)及[703入口](703/drafts/bounded_source_entry.md)。[代码](703/joint_geometry_thermal_limit.py)、[结果](703/joint_geometry_thermal_limit_results.json)、[核验](703/research_round_703_checks.json)、[全条件账](703/unified_physics_condition_ledger_703.md)。两组检查、十六式；主代理推导和审查，无独立代理审查。

## 1. 这次关闭哪个缺口

上一目标轮完成702并执行703入口，属于有效进展。本次读取导航、702报告／结果及入口；详细进程查询被系统拒绝，普通进程查询未发现Python研究，未把草稿状态当活跃进程。目标不变。

702把原完整固定图的近似自身热态、真实时间、记录及热熵接到同一个H。703入口覆盖有界读口来源；一般几何改变原动能，无界来源不能直接套入口的有界扰动界。

**本轮新增：** 在原正几何的局部参数域内，从原动能和束缚势中取出一个固定的小份额作为共同热参照。余下部分仍含原全部几何变化，保持稳定和共同域。由同一正分割得到热算符、归一态及几何来源导数的共同极限；原热应力响应、Hessian接触项、Euclidean有序来源和立即实施的原记录共用这份极限。

这是新的有限步表示，不改变任何背景处的原总H；不是在702旧处方上未加说明地微分。尚未得到近似真实时间传播的全部来源导数，也未产生动态量子几何、连续手征物质或Einstein约束。

|层次|本轮地位|
|---|---|
|认知动机|共同几何不能为态、来源、记录分别配一个过程|
|继承输入|589正几何、603热态、623共同域、624读取、原完整CAR及Gauss|
|新增表示|固定正参照A，原H(x)分成A和D(x)，有限步处方另行声明|
|成熟工具|闭扇形形式的解析半群、Schatten乘积界、向量值Vitali|
|解析连接|无界几何来源的热极限与导数交换；同一归一、接触和记录|
|数值范围|原非对角几何系数；原径向—中性Majorana量子校准|
|仍开放|近似实时间来源导数、空间细化、几何量子化、四分支身份|

### 1.1 直接继承的结果

[589](../archive_585_628/research_note_589.md)—[590](../archive_585_628/research_note_590.md)已有正几何系数和形式域；[603](../archive_585_628/research_note_603.md)、[623](../archive_585_628/research_note_623.md)已有热迹、共同算符域及原真实来源响应；[624](../archive_585_628/research_note_624.md)、[643](../archive_629_652/research_note_643.md)、[662](../archive_653_701/research_note_662.md)、[665](../archive_653_701/research_note_665.md)已有记录选择项、有序热历史、局部lapse分配及接触来源。这些不重新列为空白。

核读[Vogt–Voigt定理0.1、3.1](https://user.math.uni-bremen.de/hvogt/papers/vovo17.pdf)：共同域解析闭扇形形式给解析预解，统一半群界再给参数解析性。后面逐项验证其输入。继续使用655／702已核的[Arendt–Nikolski定理2.1](https://www.uni-ulm.de/fileadmin/website_uni_ulm/mawi.inst.020/arendt/downloads/pubbib/short/2000-AreNik-VctVldHlmFncRvs.pdf)，不将成熟定理计为新发现。

空间382—386、425、522—524直接复用：384额外Lipschitz已消去；386与425替代；523共同实现及524局域探针不重做。当前h、s、CAR／Gauss与旧空间端点的身份仍待接，辅助E不是记录s。

## 2. 原几何族中的固定热参照

固定623的M、L、u、配置测度和有限图。正几何、局部lapse及原耦合进入有限个系数，用同一Hilbert空间上的固定形式Qν写为有限和。取623参考动能及原位势：

$$
H(x)=\sum_{\nu=1}^{d}f_\nu(x)Q_\nu,\qquad
R=T_*+W_*\ge0,\quad W_*=\sum_iU_i,\qquad
H(x)\ge mR-c_0,\quad m>0 .
\tag{1}
$$

最后的界在足够小实参数邻域一致。动能系数与原位势权重有正下界；边势、磁势非负；原完整质量和跳跃用623无穷小W相对界吸收入部分原位势。各Qν在V=D(R^(1/2))上相对R+I形式有界。中心Casimir、全部Fock和Gauss保留；R按603具有所有正热时的有限热迹。

固定0<η<m，并小于原动能／原位势系数的保留余量，定义

$$
A=\eta R,\qquad D(x)=H(x)-A,\qquad H(x)=A\dotplus D(x),
\qquad D(x)\ge(m-\eta)R-c_0 .
\tag{2}
$$

D不再是702的矩阵乘法，它也含剩余原动能。扣除小份参考动能后，节点与群簇配置动能仍严格椭圆，原位势权重仍正。623的图范数和原核论证适用：D与原H具有相同共同核及域，没有对非紧标量加墙或删去费米相互作用。η是表示参数，A不是新环境；需要计算相互作用D的热算符，未声称它更便宜。

## 3. 实际来源的小复邻域

先沿有限系数的一条实解析路径，令z为复延拓。实、虚Yukawa分量作为独立实坐标延拓，不把参数共轭当全纯函数。有限系数在V形式范数中连续，缩小复圆盘使扰动足够小：

$$
\operatorname{Re}d_z[v]\ge\delta\,r[v]-c\|v\|^2,\qquad
|\operatorname{Im}d_z[v]|\le C\bigl(r[v]+\|v\|^2\bigr),
\qquad \delta>0 .
\tag{3}
$$

Re d_z+(c+1)‖·‖²与r+‖·‖²范数等价，故形式闭、扇角局部一致，D(z)+c最大增生。由上述成熟定理，

$$
z\longmapsto e^{-tD(z)}\ \text{为有界算符值全纯函数},\qquad
\|e^{-tD(z)}\|\le e^{ct}\quad(t\ge0).
\tag{4}
$$

这一步处理了无界几何变化；没有假定‖D′‖有限，也没有凭逐点自伴性跳过共同形式条件。

## 4. 原模型与同一正分割

实x时定义

$$
S_a(x)=e^{-aA/2}e^{-aD(x)}e^{-aA/2},\qquad
H_a(x)=-a^{-1}\log S_a(x),\qquad
0<S_a(x)\le e^{ac}e^{-aA},\quad H_a(x)+c\ge A.
\tag{5}
$$

因子单射，log按谱定理定义；最后使用log算符单调性，不用错误的一般指数单调性。原共同核上的切向生成元是H(x)，655强预解论证适用。702固定A低能投影和统一热尾证明遂逐背景给

$$
T_n(x;\beta):=S_{\beta/n}(x)^n
=e^{-\beta H_{\beta/n}(x)}
\longrightarrow e^{-\beta H(x)}
\quad\text{以迹范数},\qquad\beta>0.
\tag{6}
$$

每个固定实背景处，真实酉演化、原记录、能量矩和热熵按702继承。逐背景收敛还不是来源导数收敛。

## 5. 统一迹界及来源导数

复z时T_n用同一乘积定义。2n个A热因子分别属于Schatten 2n类，乘积Hölder与(4)给

$$
\|T_n(z;\beta)\|_1
\le e^{c\beta}\operatorname{Tr}_{\rm phys}e^{-\beta A}
=:B_\beta<\infty .
\tag{7}
$$

界与n及圆盘内z无关。T_n在迹类中全纯；实轴上(6)的迹范数收敛经Vitali升级为复紧集一致收敛。极限在有界算符空间中由解析唯一性识别为原热半群。Cauchy给出误差合同：

$$
\left\|\partial_z^j T_n(z_0;\beta)
-\partial_z^j e^{-\beta H(z_0)}\right\|_1
\le\frac{j!}{r^j}
\sup_{|z-z_0|=r}
\|T_n(z;\beta)-e^{-\beta H(z)}\|_1\longrightarrow0 .
\tag{8}
$$

圆周须在共同扇形域内；未给n的通用数值速率。有限多个来源用同一系数域的界，逐变量解析唯一性及Cauchy处理；混合二阶亦可由方向二阶极化求得。

实Z>0，局部一致收敛允许在实背景附近选共同无零点归一区域。因此

$$
\partial_x^\alpha\rho_n(x)\longrightarrow\partial_x^\alpha\rho(x)
\quad\text{以迹范数},\qquad
\partial_x^\alpha\log Z_n\longrightarrow\partial_x^\alpha\log Z,
\qquad|\alpha|\le2 .
\tag{9}
$$

解析坐标中任意固定阶也成立，本轮主要验收原几何一二阶来源。β固定；不把改变β时随之改变a的总导数认作有限a的能源插入。

## 6. 几何来源与Hessian接触不可分开

沿623共同域记G_a=∂_aH、C_ab=∂_a∂_bH。原热响应为

$$
\begin{aligned}
\partial_a\log Z&=-\beta\langle G_a\rangle_\beta,\\
\partial_a\partial_b\log Z
&=\beta^2\int_0^1\operatorname{Tr}
\bigl(\rho_\beta^{1-t}\delta G_a\rho_\beta^t\delta G_b\bigr)\,dt
-\beta\langle C_{ab}\rangle_\beta .
\end{aligned}
\tag{10}
$$

域和迹配对沿603／623，不将非交换Duhamel项替换为普通等时方差。式(8)—(9)把整个表达接到同一近似热泛函。有限a的来源必须微分整个分割；不能将裸G任意放到某个切口便声称精确相等。703入口已经给这种替换的有限步差异。

**无需新增“几何路径必须解析”的物理要求。** 先以有限独立系数构造局部解析热泛函，再与623允许的任意C²系数路径复合。前两阶链式法则保留

$$
\frac{d^2}{d\lambda^2}F(f(\lambda))
=\sum_{\mu,\nu}f_\mu'f_\nu'\partial_\mu\partial_\nu F
+\sum_\nu f_\nu''\partial_\nu F .
\tag{11}
$$

亦可用共同域的二阶系数jet计算基点导数；有限n与极限使用相同jet。原C²路径不必具有全邻域解析延拓。度规坐标变化、指定耦合及665接触项自身的二阶导数不能丢。固定M、L、u、配置测度和Gauss识别仍是范围，改变这些对象须另核。

## 7. 原记录及有序热来源共用这一族

对624的固定有界Kraus立即读取，记录—量子输出是同一有界CP映射作用于ρ_n：

$$
\Omega_n(x)=\sum_{\boldsymbol r}|\boldsymbol r\rangle\langle\boldsymbol r|
\otimes L_{\boldsymbol r}\rho_n(x)L_{\boldsymbol r}^\dagger,\qquad
\partial_x^\alpha\Omega_n\longrightarrow\partial_x^\alpha\Omega
\quad\text{以迹范数}.
\tag{12}
$$

参数依赖的仪器／字典须加其受控导数。各固定背景中含真实等待的记录已有702零阶收敛；**其来源导数未由本轮证明**。624原连续过程的响应存在不撤销，缺的是与当前近似族共同取导数极限。

对正热时段τ_j及固定有界插入B_j，令β=Στ_j，按每份热因子占β的比例分配Schatten指数：

$$
\left\|B_mT_{n_m}(z;\tau_m)\cdots B_1T_{n_1}(z;\tau_1)\right\|_1
\le e^{c\beta}\operatorname{Tr}e^{-\beta A}\prod_j\|B_j\|.
\tag{13}
$$

实轴上的有限乘积收敛，经相同解析论证可共同插入固定阶几何来源。保留原次序，不把Euclidean权叫真实测量概率，也不将正热时间推到未控的实时边界。奇费米反射仍需691字典；新正族不自动修复699候选。

## 8. 可复算校准与不可省略的来源项

第一组使用589原电场、梯度和磁二形式系数，包含非对角形状。|logψ|≤r_v、‖log barγ‖≤r_s、|log N|≤r_N时，边／面均值仍在同一ψ区间，各正系数相对单位参考的共同下界至少是exp[−r_N−max(6r_v,2r_v+2r_s)]。取r_v=.12、r_s=.35、r_N=.18及下界的1/4为参照份额。全参数盒界来自谱比较，18个配置只校准实现，不求全图Gauss热谱。

第二组复用702的256维原径向—中性Majorana诊断，真正改变体积动能与原位势：

$$
H(u)=e^{-6u}K+e^{6u}W+B,\qquad
G=-6e^{-6u}K+6e^{6u}W,\qquad
C=36e^{-6u}K+36e^{6u}W.
\tag{14}
$$

canonical现场质量B的体积因子按原字典抵消。单节点诊断无空间边，且沿702声明Dirac零耦合、有限径向盒与差分；原完整解析定理不删跳跃、Dirac／Majorana或Gauss。

取A=.15(K+W)、β=1.2、u₀=.03。在|Re u|≤.1、|Im u|≤.06内，

$$
\operatorname{Re}D(u)\ge
\bigl(e^{-0.6}\cos(0.36)-0.15\bigr)(K+W)+B,\qquad
e^{-0.6}\cos(0.36)-0.15>0 .
\tag{15}
$$

指数及一二阶jet用缩放平方的Taylor乘积计算，保全部非对易交叉项。原热Hessian另以本征基Duhamel积分独立核验，差约1.95×10⁻¹⁴；中央差分仅作第三份检查。

log Z二阶来源约10.8984146941，Hessian接触贡献−26.8366720040；漏掉它会得到37.7350866981。它不是可随意忽略的小修正。

n=64时，热态一、二阶导数迹范数误差约1.575×10⁻⁵和1.307×10⁻⁴，log Z二阶误差约1.140×10⁻⁶。原连续两次“+”读取的二阶概率误差约1.631×10⁻⁷，完整记录联合后态二阶误差约1.276×10⁻⁴。

单次“+”在该对称诊断恒为1/2，故正式校准用两次原读口、效果L_+⁴，避免以零响应代替检验。初版单次诊断保存在round703_drafts。有限矩阵不证明原全图热积分、仪器自主实现或连续应力理论。

## 9. 合并的条件和下一项

两份合格固定参照A、A′回到同一原模型与来源坐标：

$$
\partial_x^\alpha[\rho_n^A(x)-\rho_m^{A'}(x)]
\longrightarrow0\quad(n,m\to\infty),\qquad|\alpha|\le2 .
\tag{16}
$$

归一热来源及固定原记录同样一致。参照份额、时间分割和热来源无需各拟合一种不同响应；但没有选择背景、温度、耦合、维数、lapse分配或Einstein作用。

C01／C03／C04／C09／C20与C22之间新接通的是：**原无界几何来源及接触项与近似自身热态和立即原读取共同收敛。** 623原真实来源并非新发现，新的是这份正近似族也保它。

四分支仍分开；所有界依赖固定图及几何正余量，没有空间细化一致性、退化几何、量子约束或GR生成。接[704](704/drafts/STATUS.md)：回查621—624、643、655的热规整和实际等待，核真实记录响应与同一近似族的连接；停止本轮静态热导数和格点参数优化。
