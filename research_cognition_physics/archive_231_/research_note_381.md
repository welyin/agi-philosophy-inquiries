# 第381轮：方向层析不能认证标签忠实——Hopf纤维与相干源恢复的精确代价

日期：2026-09-23。科学基线冻结至380；不依赖任何同时进行的后继轮。14项检查通过；24个连续编号公式。

## 1. 本轮问题与结论

380轮已经给出条件性三维判据：完整qubit的方向差值能层析状态，加上实际位置边界、连续单射和全族传递协变，才推出该边界为S²。本轮只检查一个明确问题：**单射能否被其余条件及有限概率证书替代？**

答案是否定的，且不只来自多加一个离散标签。

- Hopf映射把连续S³标签压到完整S²效果轨道。其标签边界、全族传递SU(2)协变、对比完整性都满足；唯独单射失败。四设置证书与S²模型完全相同。
- 若纤维没有任何实际可读区别，就应取操作商S²；不能把纯表示相位算成额外物理方向。
- 若另有明示的相干准备源，纤维可以成为相对于源参考的可读关系。它不能从已经压成qubit密度矩阵的数据自动补回。可见度为v时，任意副本、任意合法恢复的最坏迹距误差，精确最优值为**v／2**。
- 北、南局部截面的过渡具有非零绕数，说明这一族不能写成一个全局连续的“方向加独立相位标签”。截面切换若用于相干源，必须同步调整参考相位。

最后两项针对具体输入、目标及误差量词；不把已有Hopf几何、不可测整体相位或相干控制限制当作本轮原创。

### 1.1 去重与来源

| 已有成果 | 本轮继承的边界 | 本轮增加什么 |
|---|---|---|
| [85轮](../archive_001_222/research_process/research_note_85.md)、[171轮§3](../archive_001_222/research_process/research_note_171.md)、[212轮](../archive_001_222/research_process/research_note_212.md) | 共同J、整体相位与内部参考不能混同；把参考纳入整体不等于省去参考 | 参考源是否提供被方向效果压掉的实际标签，给恢复误差量词 |
| [225轮§5](../archive_223_230/research_note_225.md)、[242轮§2](research_note_242.md) | 相同普通通道不固定相干实现；不能自动控制未知黑箱 | 区分同一S³族的密度态、普通门、已标定相干源三个接口 |
| [371轮§6](research_note_371.md) | 额外可读标签与最小性不同 | 同一个连通球面、全族传递作用及全部380其他条件同时成立的反例 |
| [374轮§5](research_note_374.md) | 已计算Berry曲率与Chern数 | 只用局部截面过渡作范围证明，不重复曲率积分或宣称新拓扑定理 |
| [378轮](research_note_378.md)、[379轮](research_note_379.md) | 构形、姿态、位置、预测状态分别识别 | 不把Hopf纤维预先认作位置；明确两种不同的定位／源合同 |
| [380轮§4、§8](research_note_380.md) | 方向差值层析与有限误差证书 | 保留其余条件后删除单射的精确失败，以及不可恢复性的最优误差 |

直接核对的原始文献：

1. [Mosseri–Dandoloff，*Geometry of entangled states, Bloch spheres and Hopf fibrations*，2001，§2.1—2.2、式(2)—(4)](https://arxiv.org/pdf/quant-ph/0108137)：一比特归一化向量的S³、Bloch球面及圆纤维。本文用标准Pauli约定自行展开矩阵与局部截面。
2. [Araújo–Feix–Costa–Brukner，*Quantum circuits cannot control unknown operations*，2014，普通黑箱反证及式(7)后的物理直和接口](https://arxiv.org/pdf/1309.7976)：仅给未知酉通道并不自动提供其相干控制；已知平凡支路的物理实现是更丰富的接口。242轮已经使用这条限制，本轮延续其资源账。

数学Hopf结构是既有知识。本地新增在于380前提审计、指定相干目标的尖锐恢复界及其可实现上界，而非再提出一个额外标签故事。

## 2. 一个满足其余条件的完整连通边界

**额外模型输入：**先给一个四维局部位置候选邻域U，取其完整边界X＝∂U≅S³。这里d＝4是用于检验条件的输入，不是由量子态球推导宇宙维数。把z视为边界点的坐标，定义一个qubit二值方向效果：

$$
X=S^3=\{z=(z_0,z_1)\in\mathbb C^2:\ z^\dagger z=1\},\qquad
E_z=zz^\dagger=\frac{I+n(z)\cdot\sigma}{2},\qquad
n(z)=\bigl(2\Re\bar z_0z_1,\ 2\Im\bar z_0z_1,\ |z_0|^2-|z_1|^2\bigr).
\tag{1}
$$

每个E_z是合法秩一投影，n(z)∈S²；所有S²向量都出现。秩一投影相同当且仅当它们张成同一复直线，所以

$$
E_w=E_z\ \Longleftrightarrow\ w=e^{i\phi}z,
\qquad
\operatorname{Tr}(E_w\rho)=\operatorname{Tr}(E_z\rho)
\quad\text{对每个qubit态 }\rho.
\tag{2}
$$

因此同一纤维中的所有标签，不仅在有限校准表中相同，而且在这个完整qubit效果接口的**全部**Born读数中相同。

定义下面的SU(2)矩阵。它的第一列恰为z，因而给出了任意两个标签之间的群元：

$$
V_z=\begin{pmatrix}z_0&-\bar z_1\\z_1&\bar z_0\end{pmatrix}\in SU(2),
\qquad (V_wV_z^\dagger)z=w,\qquad
E_{gz}=gE_zg^\dagger\quad(g\in SU(2),\ \forall z\in X).
\tag{3}
$$

这是同一个g搬运整个效果族，不是380已经排除的逐点酉等价。SU(2)在X上传递，并给四维实坐标的正交变换；不要求它是整个SO(4)，因为380并无该附加前提。

完整S²效果像当然满足方向差值层析。尤其取380的四面体nᵢ，各选任意一个提升zᵢ；可以再给每个提升任意相位而不改变概率表：

$$
E_{z_i}=\frac{I+n_i\cdot\sigma}{2},\qquad
D_{i\cdot}=\frac{(n_i-n_0)^{\mathsf T}}2,\qquad
\operatorname{sv}(D)=\frac{(2,1,1)}{\sqrt3}.
\tag{4}
$$

于是相同概率表与相同校准误差预算下，两种标签模型得到完全相同证书；特别地任何在S²表上已为正的认证裕量，在这个S³提升上也同样为正。认证的对象是效果差值对qubit的完整性，不是实际标签到效果的单射。这里得到的是

$$
\left.
\begin{gathered}
X\cong S^{d-1}\text{ 为完整球面边界},\quad d=4,\quad
z\mapsto E_z\text{ 连续},\\
\text{全族传递酉协变},\quad C_\Delta
\end{gathered}\right\}
\quad\not\Longrightarrow\quad d=3
\quad\text{在删除单射前提之后}.
\tag{5}
$$

有限概率证书不能补上这个缺口；即便效果接口的全部精确数据也不能补上。这里没有证明所有字典都失败，更没有构造满足所有未形式化认知要求的四维宇宙。

## 3. 纤维是规范，还是被指定接口遗漏的可读关系

### 3.1 只保留qubit数据时，应取商

若两个标签的区别在全部获准实验中都不可见，它们不是两个不同的完整操作状态。对式(2)的纯效果接口取商，得到S³／U(1)≅S²是完全一致的做法。

但从“选定接口只能看到S²”到“这个商就是完整实际位置边界”，仍须独立证明位置身份的含义。如果从未有别的定位依据，就不能坚持把无读出意义的相位叫第四个空间方向；如果已有独立位置或源接口识别了它，也不能以qubit压缩为理由删除已可读的区别。

本轮接下来的相干源只是一个**额外物理接口的可实现对照**。只有另有定位证据把源的标签差与实际位置对应，它才说明该qubit效果遗漏了位置身份。单纯测得相位本身，仍不足以把相位认作空间。

### 3.2 不能全局连续地选择一个相位代表

在S²北、南两片上分别取

$$
z_N(\theta,\varphi)=
\begin{pmatrix}\cos(\theta/2)\\e^{i\varphi}\sin(\theta/2)\end{pmatrix},\qquad
z_S(\theta,\varphi)=
\begin{pmatrix}e^{-i\varphi}\cos(\theta/2)\\\sin(\theta/2)\end{pmatrix}.
\tag{6}
$$

z_N在南极之外连续且北极无歧义，z_S在北极之外连续且南极无歧义。重叠部分与赤道的过渡满足

$$
z_S=e^{-i\varphi}z_N,\qquad
g_{NS}(\varphi)=z_N^\dagger z_S=e^{-i\varphi},\qquad
\operatorname{wind}(g_{NS}|_{\rm equator})=-1.
\tag{7}
$$

若存在全局连续截面，可写成z_N f_N＝z_S f_S。f_N、f_S分别延拓到两个半球圆盘，所以各自赤道绕数都为零；但f_N＝e^(-iφ)f_S要求绕数相差−1，矛盾。这是已知Hopf非平凡性的一个直接证明；数值只核对显式过渡符号，不代替一般拓扑论证。

代码对解析g_NS用64段核对绕数−1。没有声称只凭有限节点就认证任意未知连续过渡函数的绕数；若段间没有额外规则，有限样本不能排除隐藏绕行。本轮也不重复374的Berry曲率积分。

## 4. 合法相干源能读相对相位，但不是给未知态补全局相位

增加控制寄存器C和新准备的qubit Q，给出**已经标定相位的物理联合操作**：

$$
W_z=|0\rangle\langle0|_C\otimes I_Q+
|1\rangle\langle1|_C\otimes V_z,\qquad
|\Phi_z\rangle=W_z|+\rangle_C|0\rangle_Q
=\frac{|0,0\rangle+|1\rangle\otimes z}{\sqrt2}.
\tag{8}
$$

W_z是明确的四维酉矩阵。只要这个联合源已被给定，标准有限维量子操作可以实现它；并没有从未知E_z或未知Ad(V_z)黑箱中推导W_z。若通过已知控制参数准备z，参数、源和控制器本来就携带所需信息；若把z当未知位置，则必须另给源如何与位置相耦合的物理说明。

用0≤v≤1表示控制两臂的可见度，定义CQ态

$$
\Xi_{v,z}=\frac12
\begin{pmatrix}
|0\rangle\langle0|&v|0\rangle\langle z|\\
v|z\rangle\langle0|&E_z
\end{pmatrix}_{C}.
\tag{9}
$$

这是合法去相干的结果，不是任意删矩阵元。让新环境E初始为|0〉，C＝0支路不作用，C＝1支路用实酉旋转，使

$$
|e_0\rangle=|0\rangle,\qquad
|e_1\rangle=v|0\rangle+\sqrt{1-v^2}|1\rangle,
\qquad\langle e_0|e_1\rangle=v.
\tag{10}
$$

对环境作已声明的丢弃得到式(9)。环境未被说成免费重置；若保留它，可逆描述要包含其记录。

在同一纤维上，两个态的对角块相同，差只有非对角块。这个差的两个非零本征值为±v|sin(φ／2)|，故采用D(ρ,σ)＝||ρ−σ||₁／2的约定：

$$
D(\Xi_{v,z},\Xi_{v,e^{i\phi}z})
=v\left|\sin\frac\phi2\right|,
\qquad D(\Xi_{v,z},\Xi_{v,-z})=v,
\qquad D(E_z,E_{e^{i\phi}z})=0.
\tag{11}
$$

v＝1、φ＝π时两扩展态正交。一般情况下，用差算子的正谱投影Π作二值效果，概率差恰为D；Kraus算子√Π、√(I−Π)给完整仪器，没有后选择。最优二元测量假定两个候选源已知；它不是从单份未知态读出连续相位的方案。

比较式(6)的两个局部准备，未补偿的截面切换会产生可读差别；同步调整控制相位才是同一个相干源描述：

$$
D(\Xi_{v,z_N},\Xi_{v,z_S})=v|\sin(\varphi/2)|,
\qquad
\bigl[\operatorname{diag}(1,e^{i\varphi})_C\otimes I_Q\bigr]
\Xi_{v,z_S}
\bigl[\operatorname{diag}(1,e^{-i\varphi})_C\otimes I_Q\bigr]
=\Xi_{v,z_N}.
\tag{12}
$$

qubit效果层没有改变；改变的是两臂相对于参考的关系。若连参考和装置一起作相应变换，全部事实可以不变；不能把纯坐标换相位与固定参考下另一次准备混同。

## 5. 一个固定读取接口及其有限精度

令a＝|0,0〉、b_j＝|1,j〉，j＝0、1。在四维CQ上给四个固定Hermitian可观测量，其范数均为1：

$$
B_{j,R}=|a\rangle\langle b_j|+|b_j\rangle\langle a|,
\qquad
B_{j,I}=-i|a\rangle\langle b_j|+i|b_j\rangle\langle a|,
\qquad F_{j,R/I}=\frac{I+B_{j,R/I}}2.
\tag{13}
$$

每个F都是一个二值设置的合法效果。对式(9)有

$$
p_{j,R}=\frac{1+v\Re z_j}{2},\qquad
p_{j,I}=\frac{1+v\Im z_j}{2};\qquad
|\widehat p-p|\le\epsilon\text{ 逐项},\ v>0
\ \Longrightarrow\
\|\widehat z-z\|_2\le\frac{4\epsilon}{v}.
\tag{14}
$$

估计按四个线性反演分量组成；它未必恰好归一化，这不影响所报误差界。v必须已标定，若它也未知或有误差，应另增预算。四项概率一般需要许多次**同一源的重新准备**，不是四次测量，更不是把一个未知态克隆为四份。

这里确有新读出的方向。把复向量改写为四个实分量，在单位S³的切空间上，Hopf微分与四个源期望值的微分满足

$$
\ker(dn_z|_{T_zS^3})=\operatorname{span}_{\mathbb R}\{iz\},\qquad
\operatorname{sv}(dn_z|_{T_zS^3})=(2,2,0),\qquad
\operatorname{sv}\bigl(d(v\Re z_0,v\Im z_0,v\Re z_1,v\Im z_1)|_{T_zS^3}\bigr)=(v,v,v).
\tag{15}
$$

可在z＝(1,0)直接计算第一项，再用SU(2)及其正交作用推广到全部z；后项就是实四维恒等嵌入限制到三维切空间。它区分两个接口实际保留的局部参数数目，不将源参数维数自动解释成物理空间维数。

## 6. 任意副本、任意合法恢复的精确最优误差

给定v及任意整数N≥1。恢复器只接收E_z的N份副本，以及与z无关的固定辅助态τ_A。允许任意CPTP通道Λ_N；不能再把真实z作为免费经典参数传给它。目标是式(9)完整CQ态。则

$$
\boxed{
\inf_{\Lambda_N\ {\rm CPTP}}
\sup_{z\in S^3}
D\!\left(\Lambda_N(E_z^{\otimes N}\otimes\tau_A),\Xi_{v,z}\right)
=\frac v2.}
\tag{16}
$$

这里的输入接口甚至比有限概率表更强：直接允许量子操作处理副本。因而同一下界也适用于只拿有限测量记录作恢复的较弱方案。

**下界。** z和−z给完全相同输入，因此同一恢复器输出同一个σ_z。三角不等式与式(11)给

$$
v=D(\Xi_{v,z},\Xi_{v,-z})
\le D(\Xi_{v,z},\sigma_z)+D(\sigma_z,\Xi_{v,-z}),
\qquad
\max\{D(\sigma_z,\Xi_{v,z}),D(\sigma_z,\Xi_{v,-z})\}\ge v/2.
\tag{17}
$$

这个论证不随N改变，不需要假设恢复器只作测量，也不需要有限维辅助寄存器的特别结构。辅助参考若本身预装了z或与相干源已有标签相关性，就不满足输入合同，不能再引用本界。

**可达到的上界。** 丢弃N−1份输入和辅助态，保留一份ρ，执行

$$
\mathcal R(\rho)=\frac12|0\rangle\langle0|_C\otimes|0\rangle\langle0|_Q
\operatorname{Tr}\rho
+\frac12|1\rangle\langle1|_C\otimes\rho.
\tag{18}
$$

它不是非线性地选择一个相位代表，而是一个固定线性CPTP过程。显式Kraus表示为

$$
K_0=\frac{|0,0\rangle\langle0|}{\sqrt2},\qquad
K_1=\frac{|0,0\rangle\langle1|}{\sqrt2},\qquad
K_2=\frac{|1\rangle_C\otimes I_Q}{\sqrt2},\qquad
\sum_{j=0}^2K_j^\dagger K_j=I_Q.
\tag{19}
$$

对ρ＝E_z，输出就是相位平均态Ξ_0,z。与目标之差的非零本征值为±v／2，所以对**全部**z有

$$
\mathcal R(E_z)=\Xi_{0,z}
=\frac1{2\pi}\int_0^{2\pi}\Xi_{v,e^{i\phi}z}\,d\phi,
\qquad D(\mathcal R(E_z),\Xi_{v,z})=v/2.
\tag{20}
$$

积分仅用于识别输出；实施式(19)不需要无成本连续Haar随机化。这证明式(16)的精确最优值，而非仅给一个松下界或数值优化猜测。v＝0时目标本来就不含纤维关系，误差为零；v＞0时增加qubit副本不能补回已经删除的关系。

### 6.1 若真实输入泄露了一点标签

精确Hopf相同输入与有误差、可能带标签泄露的源不是同一合同。若用于±两种准备的实际输入σ₊、σ₋满足D(σ₊,σ₋)≤δ，迹距收缩给

$$
\max_{s=\pm1}D(\Lambda(\sigma_s),\Xi_{v,sz})
\ge\max\{0,(v-\delta)/2\}.
\qquad
D(\sigma_+^{\otimes N},\sigma_-^{\otimes N})\le\min\{1,N\delta\}.
\tag{21}
$$

后式由张量积差的逐项展开和迹范数性质得到。因而不能把理想模型的N无关下界不加条件地外推到每份都泄露非零标签的实际源；许多副本可能积累这种新信息。反之，只有概率估计噪声、而真实输入仍精确相同，并不提供任何标签信息。

## 7. 状态源、普通酉门、相干门是三个不同接口

为避免把同一纤维的全部门也误判为相同，直接展开式(3)：

$$
V_{e^{i\phi}z}=V_z\operatorname{diag}(e^{i\phi},e^{-i\phi}),\qquad
V_{-z}=-V_z,\qquad
D(J(\operatorname{Ad}V_z),J(\operatorname{Ad}V_{e^{i\phi}z}))=|\sin\phi|.
\tag{22}
$$

J为用正规化Bell态生成的Choi态。因此一般φ的普通门已经可辨；只有相差±时它们给完全相同普通通道。三个接口分别识别

$$
z\mapsto E_z:\ S^3/U(1)\cong S^2;\qquad
z\mapsto\operatorname{Ad}V_z:\ S^3/\{\pm1\}\cong\mathbb{RP}^3;\qquad
z\mapsto\Xi_{v,z}:\ S^3\hookrightarrow\mathcal D(\mathbb C^4)\quad(v>0).
\tag{23}
$$

最后的单射可由式(14)直接反演证明。普通门Ad(V)＝Ad(W)迫使W是V的相位倍，二者都属于SU(2)时只能为±，给第二项。

即使允许查询普通未知门任意多次、带固定辅助与适应性处理，±两个标签的黑箱说明仍完全相同，所以式(17)下界继续成立。一次普通门作用在已准备的|0〉上就能生成E_z，再用式(19)，故在至少允许一次查询的这类恢复任务中，v／2上界也可达到。这里的“查询”只有普通通道接口；若同时提供已知旁路或相干的1⊕V，实现说明就已经改变。

对任意可能与参考R关联的未知ρ_QR，**已给定**的W_z仍是合法联合酉过程。取C初态|+〉，定义Ω_±为W_±z作用后的完整CQR态，则

$$
\Omega_\pm=(W_{\pm z}\otimes I_R)
(|+\rangle\langle+|_C\otimes\rho_{QR})
(W_{\pm z}^\dagger\otimes I_R),\qquad
D(\Omega_+,\Omega_-)=1,\qquad
\operatorname{Tr}_{CQ}\Omega_\pm=\rho_R.
\tag{24}
$$

证明：W_−z＝Z_CW_z，且二者与Z_C对易；共同反演W_z之后，两输出分别变成|+〉〈+|⊗ρ_QR和|−〉〈−|⊗ρ_QR，支撑正交。酉过程当然可反演，R的边缘保持；没有声称原QR联合态在操作期间不受改变。式(8)—(14)的源层析则使用新准备的|0〉输入，不能把其结论直接套到每个未知QR输入的同一目标态上。

所以这里没有“测得孤立态的全局相位”。可测的是一个已声明更大系统中、相对于旁路的相位；把来源删掉再说它由密度矩阵自己产生，会违反225、242已经记下的接口边界。

## 8. 复算与资源账

| 对象 | 保存结果 |
|---|---:|
| S³提升后的四方向对比奇异值 | 1.15470054，0.57735027，0.57735027 |
| Hopf切向微分奇异值 | 2，2，约1.53×10⁻¹⁶ |
| v＝0.65相干源期望值2p−1的切向奇异值 | 0.65，0.65，0.65 |
| 北／南过渡绕数 | −1 |
| 赤道φ＝π切换：qubit／相干源迹距 | 0／1 |
| v＝0.65、φ＝1.2：同纤维相干源迹距 | 0.3670176077 |
| 同φ的普通门Choi迹距 | 0.9320390860 |
| v＝0.65的恢复最优误差 | 0.325 |
| 显式CPTP恢复在65个相位上的误差范围 | 0.325上下浮点误差 |
| 四源概率各误差≤10⁻⁴时的标签误差／上界 | 0.000500887／0.000615385 |

14项检查覆盖：投影与精确纤维；全族传递协变；S³提升仍通过同一有限对比证书；局部截面过渡；受控源及其协变；截面切换补偿；含环境的可见度；完整Helstrom仪器与距离；固定四读取及精度；普通门与受控门在未知参考上的区别；达到最优界的Kraus通道；多副本与通用Stinespring通道的下界；实际标签泄露的收缩界；被源接口补回的切向方向。一般CPTP、任意副本与拓扑结论均由上面解析证明承担，有限样本不承担这些量词。

本轮加入的资源包括C、空白Q、可见度环境、相位标定、源控制器、完整读取及重复准备。一个已知四维联合酉矩阵在现有操作合同中可实现，不等于一个固定有限自主处理器能精确编程整个连续族；[345轮](research_note_345.md)的资源限制仍然适用。没有证明参考永不耗散、可免费复制、可从未知E_z恢复，或能够不扰动地读取未知标签。

四个二值设置共给四项概率，不是一次四结果测量；它们无需彼此对易。给未知关联输入执行测量仍须按完整仪器记录扰动。有限概率估计的置信条件沿用既有统计接口，本轮只验证明示的确定性误差传播。

文件：[代码](direction_hopf_faithfulness_audit.py)、[结果](direction_hopf_faithfulness_audit_results.json)、[主代理核验](research_round_381_checks.json)。使用既有Python 3.12.14和NumPy 2.3.5，14项检查通过后首次保存结果；没有覆盖旧研究、安装依赖或生成图像。

~~~powershell
& 'C:/Users/admin/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B -X utf8 'research_cognition_physics/archive_231_/direction_hopf_faithfulness_audit.py'
~~~

## 9. 对三维主线的实际推进

380的对比证书、传递协变和球面边界都不能替代标签忠实性。这里真正待解释的是：**为什么选定方向载体恰好保留实际定位任务中所有已可读的方向区别？** 这不是继续增加同类qubit样本就能回答的问题。

后继研究应给该忠实接口独立的定位／传播依据，或把它坦诚保留为额外桥梁。若仅有qubit效果，应取其操作商而不虚构隐藏物理相位；若有额外相干源关系，应在整体中记录它及参考资源，不能用压缩后的S²重新定义掉已经确认的定位区别。

本轮不选择宇宙四维，也不否定380含单射前提的三维条件定理；它确定了这条前提不能由其余条件及有限概率认证自动消去。未推出自然三维、SR、GR或标准模型。
