# 固定5³格点准入诊断：独立复核

2026-10-16。审阅者：`/root/spatial_parent_next`。只核[作者说明](README.md)、[完整代码](diagnostic.py)及[保存结果](diagnostic_results.json)；不修改这三份文件、导航、轮次或目标。沿工作目录祖先及本目录检查，未发现适用的 `AGENTS.md`。本页是非编号准入审阅，不启动量子改造或预约1257。

**结论：PASS_WITH_SCOPE。** 固定非齐次初片的完整来源、正域、代数约束及局部钟能源律缺陷成立。作者的严格界与完整4500维正则诊断相容。本结果只否定这份有限离散化在该初片精确保留原逐标签钟能源律的声明；不否定小扰动有限误差合同、另一离散方案或研究纲领。当前实际任务尚未准入。

## 1. 核查对象和独立性

使用原1255第一互动钟及1256第二记录流字典。逐项检查了密度动量权重、度规非对角因子2、全部Ricci项、直接差分的密度协变散度、第二流的复合密度与相位动量、完整root以及其逆度规变分。代码用完整Hamiltonian求全部4500个正则偏导，并未将各格当作七对齐次模型，也未借用连续HDA作为离散恒等式。

以下第2节的解析重建独立完成于完整代码落盘之前。其后实际读取全部代码，并于2026-10-16运行默认只读复算；用时约9.2秒，返回：

```text
READ_ONLY_RECOMPUTATION_MATCH
INITIAL_DEFECT_CERTIFIED_TASK_NOT_ADMITTED
FULL_GRID [5, 5, 5] DIMENSION 4500
delta 0.0 Jdot [0.0, 0.0, 0.0, 0.0, 0.0]
delta 0.1 Jdot [0.0, -0.00284306680762, 0.0, 0.00142153340381, 0.00142153340381]
```

这属于独立解析审查加原代码的独立执行，不伪称另写了第二份完整4500维实现。复杂步值仍是浮点诊断；严格非零界由解析式及精确根外包承担。

绑定的作者文件SHA-256：

|文件|SHA-256|
|---|---|
|README.md|`bd7641b9ddbdca7439ffa565024287e4eccecafdeacb55e51719b1c5b77d90c1`|
|diagnostic.py|`e9bbcda917d9b9cebcffbc7c474320c418f952e50fd2e48034c3da3f97a4a759`|
|diagnostic_results.json|`27cbe825d80e2dda8ef9fc834055bcd50cdeb0aa5a2a8ffc79174e7968cd9035`|

## 2. 完整正则贡献与lapse Hessian抵消

以下只针对作者指定初片，沿另两方向重复，令 $D=D_1$、$\beta=\beta^1$，并写

$$
h_{ij}=\delta_{ij},\qquad \pi^{ij}=P\delta^{ij},\qquad
P=-\frac{\sqrt{111/100}}3,
$$
$$
A=A_0-5d^2,\quad A_0=\frac{187}{200},\quad
N=\frac A{\mathcal R},\quad\beta=\frac d{\mathcal R}.
$$

在常度规、常 $P$ 的初片，作者直接离散 $D_g$ 的一阶变分为

$$
\delta D_{g,1}=-2D\delta\pi^{11}
-PD\delta h_{11}+PD\delta h_{22}+PD\delta h_{33}.
$$

这只是一阶初片身份；不将其与另一种离散Lie生成元在一般数据上等同。由完整 $H_\Delta$ 得

$$
\dot h_{11}=-6PN+2D\beta,\qquad
\dot h_{22}=\dot h_{33}=-6PN.
\tag{R1}
$$

令 $G$ 表示总 $C_0$ 对一个对角度规分量的无空间导数部分，在本初片有

$$
G=-\frac32P^2+\frac{11}{400}-\frac52d^2.
\tag{R2}
$$

这里包含引力动能、原 $X$ 势能、三杆空间梯度以及新增第一杆动能对度规的变分。不能把它当作纯引力贡献；第二流在机械动量为零的此初片，其能源 $W$ 的度规变分为零。

完整root的显式逆度规变分在 $\dot\pi^{11}$ 中给 $+d^2/(2\mathcal R)$。其余对角式为

$$
\dot\pi^{11}=-NG-PD\beta+\frac{d^2}{2\mathcal R},
$$
$$
\dot\pi^{22}=\dot\pi^{33}=-NG+PD\beta-\frac16D^2N.
\tag{R3}
$$

线性曲率变分为

$$
\delta C_{g,\mathrm{curv}}=\frac16D^2(\delta h_{22}+\delta h_{33}).
$$

所以其对 $\dot C_0$ 的贡献为 $-2P D^2N$，而(R3)经 $\partial C_0/\partial\pi^{aa}=-6P$ 给出的lapse Hessian贡献为 $+2P D^2N$，两者准确抵消。这里没有省略完整几何变分，也没有先假定离散HDA。

对作者见证中的 $d_i=0$ 格点，$(h,\pi_h)$ canonical配对对总 $\dot C_0$ 的贡献因而为

$$
(2G-6P^2)D\beta
=\left(-9P^2+\frac{11}{200}\right)D\beta
=-\frac{211}{200}D\beta.
\tag{R4}
$$

第一杆的完整变化是

$$
\dot\phi^1=\frac{Nd}{k}+k\beta,
\qquad
\dot\Pi_{\phi1}=kDN+\frac1kD(\beta d),
$$
$$
\left.\dot C_0\right|_{\phi^1,\Pi_{\phi1}}
=dDN+10dD(\beta d)+D(Nd)+\frac1{10}D\beta.
\tag{R5}
$$

在 $d_i=0$ 点前两项为零。其余两杆在本一方向初片的对应正则贡献为零，但其度规变分已在(R2)中保留。$p,X$ 的初始时间动量和梯度为零，故其正则配对在此刻不另添项，其度规来源同样没有删除。

第二记录流不可省略：初始 $\eta=\Pi_\vartheta=1/10$、机械动量为零，逆标签与相位输运给

$$
\dot\eta=\frac1{10}D\beta,\qquad
\dot\Pi_\vartheta=\frac1{10}D\beta,
\qquad
\left.\dot C_0\right|_r=\frac1{50}D\beta.
\tag{R6}
$$

合并(R4)—(R6)，获得

$$
\dot C_{0i}=D(Nd)_i-
\left(\frac{211}{200}-\frac1{10}-\frac1{50}\right)D\beta_i
=D(Nd)_i-A_0D\beta_i.
$$

利用 $Nd=A\beta$，得到作者式

$$
\boxed{\dot C_{0i}=-5D(d^2\beta)_i.}
\tag{R7}
$$

该点 $D_0=0$，故 $\dot{\mathcal R}=-\dot C_0$；初时内部变量及两窗亦给总 $\dot U=0$。因此

$$
\dot J_i=+5D(d^2\beta)_i.
$$

$\dot C_0$与$\dot J$的号相反。作者当前说明及代码采用此正确号。

## 3. 固定严格界与证据范围

取 $\delta=1/10,\Delta=2$，则

$$
\mathcal R_\delta^2=\frac{30929}{40000},\qquad
\frac1{25}<J_\delta<\frac9{200},
$$
$$
\dot J_1=-\frac{5\delta^3}{\Delta\mathcal R_\delta}
=-\frac1{400\mathcal R_\delta}<-\frac1{400}.
$$

总法向材料能源仍是 $111/100$，杆与钟动量逐格相抵，正root给代数Hamiltonian约束零。新增空间应力包含杆的各向同性动能部分以及钟机械运动的各向异性部分；仅保持 $T_{nn}$不表示完整原状态或空间应力没有改变。

本例证明的是：**指定离散方案在一份约束和正域合格的有限初片，未精确保留原逐标签钟能源律。** 它不是新连续解、固定未来读数差、全场量子反例或新公理非蕴涵；也不证明所有非齐次数据都有同样缺陷。作者未把非零 $d$ 点的浮点近零另作解析签署，这一限制恰当。

936、937已有一般约化、有限量子工具及约束范围警示。这里的净证据是原1255—1256材料与这份固定差分的具体同源失配，不应把一般“离散约束可能失真”重报为新理论。

## 4. 钟余额函数的偏导约定

必须区分总导数与固定剩余canonical资料的显式偏导。设 $z$ 为本离散系统全部剩余canonical变量：

$$
J_i(t,z)=\frac{\mathcal R_i(z)}{n_c}-U_i(t,z),\qquad
\frac{dJ_i}{dt}=\{J_i,H_\Delta\}
+\left.\partial_tJ_i\right|_z.
$$

原钟能源律要求

$$
\frac{dJ_i}{dt}=-\left.\partial_\theta U_i\right|_{\text{其余原资料}}
=-\left.\partial_tU_i(t,z)\right|_z,
$$

后一等号使用当前材料时间规范 $\theta=t$；右边只微分 $w(t),f(t)$ 等显式钟依赖，**不**微分正在变化的 $p,X,q,Y$等变量。因此合适的离散残差是

$$
\boxed{\Gamma_i(t,z)=\{J_i,H_\Delta\}
+\left.\partial_tJ_i\right|_z
+\left.\partial_tU_i\right|_z.}
\tag{R8}
$$

对本目录的未完成母Hamiltonian，$\mathcal R_i(z)$无显式时间依赖，所以两个偏导恰相消，$\Gamma_i=\{J_i,H_\Delta\}$。不能将最后一项替成总 $dU_i/dt$，否则会错误增加 $\{U_i,H_\Delta\}$。本初片两者碰巧都为零，单靠初片数值无法检出这种误读。

若以后加入显含时间的完成项、改写重建 $J$ 或变更母对象，必须重核(R8)与原钟能源身份；不能不经检查继续使用 $\Gamma_i=\{J_i,H\}$。任何量子版本也须说明采用的同一 $J,U,H$ 字典和排序，不能仅把当前经典等号换成交换子便宣称实际量子钟已实现。

## 5. 小波包仍可能成立的条件及停止线

固定 $\delta=1/10$ 的见证不能排除更窄的实际准备域。原均匀中心在该离散Hamiltonian中不变的解析理由与作者第4节一致：原仿射背景精确、常背景差分和曲率总一阶变分消失，三杆应力保持各向同性；原第二流相对资料并未从相空间删去。因此可把旧1247齐次历史作为下一任务的中心，但仍须对实际将使用的完整域登记正余量与误差。

有限小波包路线的必要可检验条件是：

1. 原两条齐次中心在声明全窗确为**同一**离散Hamiltonian的解，所选钟余额和弱Einstein残差在中心为零；不能仅由本代码的 $t=0$检查推断全窗。
2. 在事前固定正根、正 $J$、正标签Jacobian的紧域内，这些残差与来源有可计算导数上界；该域包住中心并有非零距离余量。
3. 实际量子准备保留全部格点相对模，并控制其在该域内的偏离和域外的来源加权尾；不能只给平均态或总能源误差。
4. 将真实cell、末端访问窗口与有限弱几何/来源菜单同时固定，检查所有预算仍保原有限记录分离。当前粗cell与短访问窗的身份差额并未由此审阅消除。

若 $\Gamma_i(z_{\rm cl}(t))=0$，光滑性可给该紧域上的

$$
|\Gamma_i(t,z)|\le L_i\,\operatorname{dist}(z,z_{\rm cl}(t)),
$$

其时间积分控制同一候选的钟余额累计误差。量子应用还必须计入排序、实际状态偏离和坏域尾；这是可检验路线，不是本次已取得的有限时间量子证书。弱场方程残差也应直接作为本候选的菜单检验，不要求为有限有效任务先恢复全部高频连续HDA。

Duhamel只能在已指定两份传播时控制它们的差。例如同一有限Hamiltonian的真实量子传播与其二次近似可以比较；当前没有另一份已合格的“精确GR生成元”，故不能据一般Duhamel形式宣称已修复离散HDA。

**停止判据：**若中心本身存在不消失的余额，或新完成使中心失配，则单纯缩窄准备宽度无效；若实际cell/原语仍无共同身份，也不能用量子误差小补足。反之，上述条件可获得时，允许在有限预算内继续，不把这份离散候选的缺陷升级成整个阶段的失败。

本页不认证阶段结项、Jacobson热—几何身份、完整自治量子钟、实际三结果仪器或全量子GR。无新科学轮次及组数。
