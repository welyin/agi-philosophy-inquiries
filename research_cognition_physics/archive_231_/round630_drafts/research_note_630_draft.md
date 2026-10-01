# 第630轮：原物质谱、真实时间响应与几何来源的共同连接

日期：2026-10-01。接[629](research_note_629.md)，复用[599](research_note_599.md)、[600](research_note_600.md)、[620](research_note_620.md)及[625—627](research_note_625.md)。[代码](joint_continuum_source_spectrum.py)、[结果](joint_continuum_source_spectrum_results.json)、[核验](research_round_630_checks.json)、[全账回填](unified_physics_condition_ledger_630.md)。三组检查、十九式；主代理审查，无新增独立代理审查。

## 1. 问题、继承范围与实质新增

现在优先整合条件，认知实现继续后置。599的连续费米有效作用与623—625的原图量子过程，不能仅凭质量相同便认定已相连。本轮先在599已经声明的连续候选内，补出同一整代质量的正实时谱，并检查它是否同时约束局部作用、量子噪声和一类几何来源。

**新增连接：原完整一代的径向质量来源和共形度规来源，共用一个正的成对产生谱。它的两个对数UV系数准确复现599的质量和二梯度极点；谱还共同固定非局部响应、真空噪声及低频展开的适用窗口。** 来源交叉项被这个共同谱约束，不能分别配置。

这不是再次计算反常表、热核导数次数或一般涨落—响应身份。这里给原16-Weyl质量的明确Lorentz谱、与旧热核的归一连接、具体几何来源矩阵和解析误差界。

|层次|内容|
|---|---|
|认知动机|物质、测量扰动与几何来源共用同一个量子过程|
|继承输入|598诊断Y、原F和质量函数；599同一径向方向及canonical连续费米动能|
|本轮明确分支|给定3+1平直Einstein背景、零规范背景、常量非零标量、真空；一代全部费米子量子化|
|实际计算|自由费米真空的连通二点函数，等价于外部来源的一个费米圈；任意原非零诊断质量，不作重质量近似计算谱|
|附加几何比较|光滑共形度规及同一径向来源；变换在远过去／未来回到原背景，使用相容的初始真空|
|未接通|原图热态到本真空的映射、完整相互作用、一般背景张量应力、实际记录仪器、动态引力|
|未减少的输入|维数、真空选择、canonical场论、参数值及局部有限匹配常数|

这条连续候选是总目标的一段明确接口。既未替换掉原图候选，也没有把局部成功视为统一模型完成。

## 2. 从原完整质量取得共同方向

复用599的单位径向方向，M表示原常数2，矩阵质量写作 \(\mathsf M\) 以免混淆：

$$
\phi(q)=\sqrt{6M}\tanh(q/\sqrt6)\,\hat n,\quad
\hat n=(2,0,0,0,1)/\sqrt5,\quad
F=M-\phi^2/6,\quad
\mathsf M(q)=\sqrt6\sinh(q/\sqrt6)\,\mathsf M_0 .
\tag{1}
$$

\(\mathsf M_0\) 是原一代16维复对称Weyl质量：三个颜色的u、d Dirac块、一个e Dirac块，以及原ν与singlet Majorana的二阶块。所有Y继承598，没有再拟合。q★=0.27是599既有诊断幅度，不是物理真空预测。质量沿这条射线只改变共同实幅度；一次常量Takagi变换足够，不产生随位置变化的质量本征基联络。

令m_a为 \(\mathsf M(q_\star)\) 的16个正奇异值，每项取半权：

$$
d_a=\frac12,\quad
\sum_a d_a=8,\quad
S_2=\sum_a d_am_a^2=\frac12\operatorname{tr}\mathsf M^\dagger\mathsf M,\quad
S_4=\sum_a d_am_a^4=\frac12\operatorname{tr}(\mathsf M^\dagger\mathsf M)^2 .
\tag{2}
$$

这相当于七个Dirac和两个Majorana质量；Dirac重复两次后恢复权1。不可再乘两次spin或Nambu数。代码另从598原64维BdG质量谱核对全部奇异值，误差2.78×10⁻¹⁶。

本点F=1.97589548335，S₂=0.328052986277，S₄=0.0249188434630；最小质量0.05127563160来自原e诊断块，最低阈值ω₀=0.10255126320。本轮没有使用真实粒子质量排序。

## 3. 同一质量来源的正实时谱

先用线性外部来源η定义质量 \(\mathsf M_\star(1+\eta)\)。在平直背景中，O是Hamiltonian对η的一阶导数，即完整的质量双线性算符；各Majorana项含标准半因子。采用 \(\hbar=c=1\)，空间零动量及单位体积归一。定义连通谱和迟致核：

$$
\rho(\omega)=\int dt\,d^3x\,e^{i\omega t}
 \langle[O(t,\mathbf x),O(0)]\rangle,\quad
\Pi_R(t)=i\theta(t)\!\int d^3x\,\langle[O(t,\mathbf x),O(0)]\rangle,\quad
N(\omega)=\frac12\rho(|\omega|).
\tag{3}
$$

最后一式限真空，且ρ奇延拓。此约定 \(\operatorname{Im}\Pi_R(\omega+i0)=\rho(\omega)/2\)；因为Hamiltonian扰动为＋ηO，实际一阶期望变化是−Π_R η。局部接触项与来源的非线性二阶导数另计，避免符号或“静态曲率＝正协方差”混淆。

对ω>0，原完整物质给

$$
\rho(\omega)=\sum_a\frac{d_am_a^2\omega^2}{4\pi}
 \left(1-\frac{4m_a^2}{\omega^2}\right)^{3/2}
 \theta(\omega-2m_a).
\tag{4}
$$

这是成熟自由费米成对阈值的应用；标量耦合的三次速度因子亦见[Djouadi原综述§2.1.1](https://arxiv.org/pdf/hep-ph/0503172)。没有使用该综述的旧实验数值。为固定本项目归一，直接在canonical Hamiltonian中计算：

$$
H_{\mathbf k}=\boldsymbol\alpha\cdot\mathbf k+\beta m,\quad
P_\pm=\frac12(1\pm H_{\mathbf k}/E_{\mathbf k}),\quad
\operatorname{tr}(P_-m\beta P_+m\beta)
 =\frac{2m^2|\mathbf k|^2}{E_{\mathbf k}^2},\quad
E_{\mathbf k}=\sqrt{\mathbf k^2+m^2}.
\tag{5}
$$

将(5)乘 \(2\pi\delta(\omega-2E_{\mathbf k})\,d^3k/(2\pi)^3\) 积分便得Dirac权1的(4)。Majorana取半权。解析归一不依赖数值拟合；原全部质量的投影矩阵检查误差1.39×10⁻¹⁷。

真空正性同时给谱和噪声非负。这里没有把热态的零频占据记忆设为零后还称之为同一热态；本轮明确选择真空。602和626对其它状态的限制继续有效。

## 4. 实时谱与599局部有效作用的共同系数

形式未减除的色散关系及Euclidean连通核为

$$
\Pi(z)=\frac1\pi\int_0^\infty
 \frac{\omega\rho(\omega)}{\omega^2-z^2}\,d\omega,\quad
\operatorname{Im}z>0,\qquad
K_E(Q)=-\Pi(iQ),\quad K_E=\Gamma_E^{(2)} .
\tag{6}
$$

(6)的裸积分发散，不能当有限答案。高频展开是

$$
\rho(\omega)=\frac{S_2}{4\pi}\omega^2-\frac{3S_4}{2\pi}
 +O(\omega^{-2}).
\tag{7}
$$

因而同一正谱同时决定Π的质量二次发散、质量对数项和z²对数项。二次发散依赖调节方案；两个对数项可与599核对。该轮在相同canonical连续处方、Majorana计数及给定度规下，质量有关的局部系数为

$$
b_{4,\mathrm f}
 =2S_4(1+\eta)^4+2S_2(\partial\eta)^2
 +\frac{R}{3}S_2(1+\eta)^2,\qquad
\Gamma_{\mathrm{div}}=\frac1{32\pi^2\epsilon}\int\!\sqrt g\,b_{4,\mathrm f}.
\tag{8}
$$

纯曲率项没有被删除，只未写在(8)；也未把本轮当成规范、标量或引力圈总和。热核公式来源已由599核实为[Vassilevich式4.28](https://arxiv.org/pdf/hep-th/0306138)，此处不重新计为新定理。

对(8)求两次η导数，以 \(1/\epsilon\leftrightarrow2\log\Lambda\) 比较对数项，并使用(6)的符号，得到

$$
\Pi_\Lambda(z)=\frac{S_2}{8\pi^2}\Lambda^2
+\left(-\frac{3S_4}{2\pi^2}+\frac{S_2}{4\pi^2}z^2\right)\log\Lambda
+\text{有限项},\qquad
(\ell_0,\ell_2)=(-0.00378721007200,\ 0.00830967921674).
\tag{9}
$$

Λ在(9)是频率截止；它不是Lorentz协变正则化的完整定义，也不用于预测物理真空能。对数比较不宣称两种方案的有限常数相同。曲率—质量项仍须随一般度规共同变分，不能由平直零动量切割独立确定其全部有限系数。

代码从谱积分的截止端点导数求(9)，不作拟合；在Λ=500m_max时两系数分别为−0.00378719535和0.00830951013，余差符合高频幂次。另直接算原同一自由海能的二阶质量导数，与其零频谱积分核对，和的误差5.59×10⁻¹⁵。这给有效势曲率、动能对数系数与实际实时谱的具体共同连接。

## 5. 有限的共同核与适用尺度

两次减除后写

$$
\Pi_{\mathrm{ren}}(z)=c_0+c_2z^2+\mathcal R(z),\qquad
\mathcal R(z)=\frac{z^4}{\pi}\int_0^\infty
 \frac{\rho(\omega)}{\omega^3(\omega^2-z^2)}\,d\omega .
\tag{10}
$$

c₀、c₂为有限匹配资料，未由谱或认知原则选出；不同有限高阶局部来源项若另行允许，也必须共同列账。以下“余项”专指(10)约定的色散部分，不把原标量／几何高阶作用偷设为零。非局部切割和真空噪声由(4)固定，实局部多项式不能修改它们。

对原全部m_a>0：

$$
L_4=\frac1\pi\int_0^\infty\frac{\rho(\omega)}{\omega^5}\,d\omega
 =\frac{\sum_ad_a}{80\pi^2},\qquad
L_6=\frac1\pi\int_0^\infty\frac{\rho(\omega)}{\omega^7}\,d\omega
 =\frac1{1120\pi^2}\sum_a\frac{d_a}{m_a^2}.
\tag{11}
$$

代换x=2m_a/ω把积分化为[0,1]上的Beta积分，分别为1/5与2/35。对实频率 \(|z|<\omega_0=2m_{\min}\)：

$$
0\le\mathcal R(z)\le
 \frac{z^4L_4}{1-z^2/\omega_0^2},\qquad
0\le\mathcal R(z)-z^4L_4\le
 \frac{z^6L_6}{1-z^2/\omega_0^2}.
\tag{12}
$$

本点L₄=0.0101321183642，L₆=0.0647834323545。用192／384点两份独立阶数的Gauss积分核对四个亚阈值频率与解析界，最大积分变化5.17×10⁻¹⁸；这是数值复算精度，不是全理论误差。式(12)才是所声明谱部分的解析展开界。

最低质量同时限制所有来源共用的低频展开域。沿原射线q→0，质量和阈值趋零，L₆发散，不能把本窗口推广为整个原场空间的统一窗口。质量严格为零时对应O也变零，不能交换“先低频展开”和“先取零质量”两个极限。627已经排除全态统一重质量隙，本轮没有重复该证明。

严格带限来源不能同时严格紧支于有限时间；实际有限时间来源还要计超阈值频谱尾。不能仅因脉冲变化慢，就把(12)当全部有限窗操作的误差证书。

## 6. 同一谱怎样约束几何来源

取 \(g_{\mu\nu}=e^{2\sigma}\eta_{\mu\nu}\)。在四维canonical费米动力学中，使用 \(\psi_g=e^{-3\sigma/2}\psi_{\rm flat}\) 后，质量成为

$$
\mathsf M_{\mathrm{flat}}=e^\sigma(1+\eta)\mathsf M_\star,\qquad
1+\zeta=e^\sigma(1+\eta).
\tag{13}
$$

这是给定维数和动能的共形协变，不是推导维数；量子变量变换的局部Jacobian／迹反常、重整化项及一次来源造成的接触项必须保留。它们不改变平直真空中非重合点的切割谱。对于原径向坐标，

$$
\eta(q)=\frac{\sinh(q/\sqrt6)}{\sinh(q_\star/\sqrt6)}-1,\quad
\delta\zeta=b\,\delta q+\delta\sigma,\qquad
b=\frac1{\sqrt6}\coth(q_\star/\sqrt6)=3.71869156775 .
\tag{14}
$$

在作用的负来源导数约定下，q和σ的质量来源为bO与O；几何项相应为负的应力迹，方程运动项只引入接触差别。故切割、非局部迟致部分与非重合真空噪声共同有

$$
\rho_{AB}=v_Av_B\rho,\quad
N_{AB}=v_Av_BN,\quad
\Pi^{\mathrm{nonlocal}}_{AB}=v_Av_B\Pi^{\mathrm{nonlocal}},\qquad
v=(b,1),\quad A,B\in\{q,\sigma\}.
\tag{15}
$$

**这不是完整Hessian秩1，更不是所有应力分量秩1。** 变量变换二阶导数乘一次期望、局部几何项和迹反常给接触项；横向无迹度规变化、其它标量方向及其它量子内线尚未包含。620证明过一般交叉来源必须保留；这里新增实际连续物质谱对交叉项的定量固定。

在本分支，补偿路径可直接写为

$$
q(\sigma)=\sqrt6\operatorname{arsinh}
 \left[e^{-\sigma}\sinh(q_\star/\sqrt6)\right],
\qquad e^\sigma\mathsf M(q(\sigma))=\mathsf M_\star .
\tag{16}
$$

它消去这个费米部门的非局部质量扰动；局部Jacobian和原标量／引力作用仍可能变化，故不是全理论的新规范对称性。原完整BdG质量在此路径上的数值误差≤5.73×10⁻¹⁷。

例如ω=4m_max时，(15)的谱矩阵为

$$
\rho_{AB}\simeq
\begin{pmatrix}
0.363411217717&0.097725560482\\
0.097725560482&0.026279555242
\end{pmatrix}.
\tag{17}
$$

向量(1,−b)的二次型为舍入量−1.33×10⁻¹⁷；若保两个对角谱而删除交叉谱，得到0.726822435435。此见证直接排除“物质和几何各配一份独立噪声，再相加”的指定接法，未排除其它完整共同过程。

## 7. 局部有效作用何时不够

在本真空分支，实局部来源作用的有限多项式近似不能单独表示开放的成对通道：

$$
\omega>2m_{\min}:\quad
\rho(\omega)>0,\quad N(\omega)>0,\qquad
\operatorname{Im}P_{\rm real}(\omega)=0 .
\tag{18}
$$

这里P_real指消去费米子后只留下的实局部二次来源核，不是含显式费米子的局部场论。保留显式物质及同一真空，或保留完整因果影响核，均可表示这些通道。少数局部系数相同不证明完整量子过程相同。

作为实际非绝热响应诊断，取外部η脉冲 \(\epsilon e^{-t^2/(2\tau^2)}\)，按原自由场真空计算最低非平凡阶的吸收能量密度：

$$
\frac{\Delta E}{V}
 =\int_0^\infty\frac{d\omega}{2\pi}\,
   \omega\rho(\omega)|\widetilde\eta(\omega)|^2+O(\epsilon^3),
\qquad
\widetilde\eta=\epsilon\sqrt{2\pi}\tau e^{-\tau^2\omega^2/2}.
\tag{19}
$$

τ=1、3、8时，二阶系数分别为0.00876001780、0.000162094298、0.00000128743233。计算截在变换变量y=τ(ω−2m)≤12，指数尾受e⁻¹⁴⁴乘多项式控制；这些数字是诊断积分，并非已构造内部实验装置。脉冲无限时间高斯尾也没有被冒充严格有限支集。真空涨落、吸收与(9)局部对数来源使用同一个谱，没有三套独立系数。

## 8. 对共同条件账的实际贡献及下一步

- C01／C15／C17／C19：原整代质量在指定连续真空分支具有明确正实时谱；不把该真空当原Gauss热态。
- C10／C12／C20／C22：质量、局部量子系数、噪声和共形几何来源共同受(4)、(9)、(15)约束；有限匹配常数仍独立。
- C20／C25：同一最低阈值控制亚阈值导数展开；越过阈值须保因果非局部核或显式物质，不能以实局部作用单独签收。
- C02／C03／C21：本轮没有连续区域仪器或内部准备设计；不把二点谱当完整记录过程。
- 628的反常结论与629的拓扑角依然成立于各自范围，本平直平凡束计算不选择拓扑角，也不补完一般局域手征测度。

完整C01—C27及623—629的更新见[条件账](unified_physics_condition_ledger_630.md)，整理部分不另计轮次。接[631](round631_drafts/STATUS.md)：将这项明确的连续响应放回共同几何条件，优先核一般度规变化中新增的应力张量通道能否与同一原物质和旧几何有效作用匹配。共形迹方向不能代表全部度规；不延长脉冲／谱积分精度优化，不进入认知装置设计。原图到连续理论、实际态与尺度映射仍是最终必须补的连接。

