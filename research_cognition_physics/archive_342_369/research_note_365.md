# 第365轮：固定 Gauss 约束下的未知主体加入、桥记录回收与内部执行时钟

日期：2026-09-23。已读项目及研究目录导航、[第362轮最终稿](research_note_362.md)、[第360轮](research_note_360.md)及其结果，检索旧档案，并复读[第346轮](research_note_346.md)自主时钟。科学基线为已冻结第362轮，不依赖第366轮等未冻结工作。本轮21项检查、19个连续编号公式，无图像检查。

## 1. 关闭一个具体的加入缺口

第362轮的末态编码保留了未知状态和参考，但编码器不与同一组固定约束对易。它不能直接当作全过程满足 Gauss 法则的加入动力学。

本轮得到一个**有明确输入条件的正协议**：

- 两个旧主体已经各自处于给定规范物理空间，携带任意未知逻辑态及任意有限参考。
- 一个已计入总空间的中性空白桥，使初态从开始就满足同一组扩展 Gauss 约束。
- 七个有连续生成元的操作，全程保持该组约束，并在两个旧逻辑比特之间真正产生可测交互。
- 桥在中途承载量子关联，最后经过相干逆计算回到空白；没有将其偷偷丢弃。
- 复用346轮的工程化历史时钟，可将七门编入一个时间不变、仍严格保持全部 Gauss 约束的 Hamiltonian，在指定时刻确定性完成。

这不是任意裸未知态可以免费进入任意约束空间的定理；它也没有生成规范群、关系图、自然耦合或引力。桥和时钟的制备、最终读取与长期存档仍需独立资源。

| 层次 | 本轮地位 |
|---|---|
| 认知动机 | 主体在保留未知信息及可操作能力时，参与更大整体的真实交互 |
| 额外输入 | 有限简单图、物质／链路比特、Z₂ Gauss 法则、旧物理编码、空桥、允许生成元及其强度 |
| 解析证明 | 初态相容、全时间约束保持、任意参考交织、真实旧主体相互作用、桥精确回收、两体权限限制 |
| 既有成果复用 | 规范物质带衣操作；346轮内部历史时钟及完美端点传输 |
| 数值核验 | 直接矩阵与解析指数、各中间时刻、完整逻辑仪器、参考、桥重置反例、内部时钟谱与回流 |
| 未得到 | 任意入射裸态的免费准备、自然图／Hamiltonian选择、几何局域实现、永久输出或 Einstein 动力学 |

## 2. 对接的已有结构

[Homeier 等，*Realistic scheme for quantum simulation of Z₂ lattice gauge theories with dynamical matter in (2+1)D* (2023)](https://arxiv.org/html/2205.08541)，式 (1) 使用物质奇偶乘入射电场的 Gauss 生成元。只对链路交换 X/Z 基约定后，正是下文的形式；其式 (4) 中相等强度的 hopping 与 pairing 合成物质—链路—物质的三比特翻转。这是成熟的规范物质结构，本轮不把它称作认知原理的新推导。

该文还明确区分精确有效规范项与实际微观装置的近似保护：高阶违反、适用时间及热力学极限另有条件。本轮直接指定精确允许项，不借该论文证明现实两体装置已经无误实现全部协议。

自主执行继续使用346轮核验过的 [Caha–Landau–Nagaj 历史时钟，§2.1](https://arxiv.org/html/1712.07395)和 [Christandl 等工程化完美传输链，式 (13)—(15)](https://arxiv.org/html/quant-ph/0309131)。新增接口是：**固定 Gauss 对易条件在这一编译中严格保留**，并逐项核对加入协议的未知输入和资源。

## 3. 固定模型与全部物理逻辑自由度

给定有 V 个顶点、E 条无自环、无重边的有限图。每顶点放物质比特 m_v，每边放链路比特 e，取 ℏ=1：

$$
\mathcal H_{\rm kin}=(\mathbb C^2)^{\otimes(V+E)},\qquad
G_v=Z_{m_v}\prod_{e\ni v}Z_e,\qquad
\mathcal H_{\rm phys}=\bigcap_v\ker(G_v-I),\qquad
\dim\mathcal H_{\rm phys}=2^E .
\tag{1}
$$

V 个约束相互对易，并因各自含有独有的物质 Z 而独立。以下所有时间演化使用同一组 G_v；没有沿过程替换约束算符。

设 B 为 GF(2) 上的 V×E 关联矩阵，每列在两个端点为1。用链路比特串 s 标记全部物理态：

$$
\mathcal V_\Gamma|s\rangle
=|Bs\rangle_m|s\rangle_e,\qquad
\mathcal V_\Gamma^\dagger\mathcal V_\Gamma=I_{2^E},\qquad
G_v\mathcal V_\Gamma=\mathcal V_\Gamma,\qquad
\mathcal V_\Gamma\mathcal V_\Gamma^\dagger=P_{\rm phys}.
\tag{2}
$$

这里的等距只用来标明**已经有效的输入对象**及操作表示，不把实现此编码的任意线路冒称为保持固定 G 的初态制造过程。未知叠加在线性编码中整体保留，不是复制未知量子态。

对每条 e=(v,w)，完整逻辑 Pauli 可以选为

$$
\overline X_e=X_{m_v}X_eX_{m_w},\qquad
\overline Z_e=Z_e,\qquad
\overline Y_e=i\overline X_e\overline Z_e,\qquad
[\overline O_e,G_u]=0,\qquad
\overline O_e\mathcal V_\Gamma=\mathcal V_\Gamma O_e .
\tag{3}
$$

每个端点处的两次反对易相消；不同边的逻辑 Pauli 相互对易。同边满足完整 Pauli 关系，因而全部边生成物理空间上的 M_(2^E)。原能力是这些物理逻辑算符的能力，不是裸链路自身边缘的完整量子态。

## 4. 从一开始就满足约束的两个旧主体与空桥

最小具体实例：旧边 A=(0,1)、B=(2,3)，新增桥 C=(1,2)。先把桥寄存器计入总系统并准备 Z_C=+1，即 |0⟩_C。四个固定约束为 Z_m0 Z_A、Z_m1 Z_A Z_C、Z_m2 Z_B Z_C、Z_m3 Z_B。

用 a、b 表示旧逻辑基值，入射等距为

$$
\begin{aligned}
W_0|a,b\rangle
&=\mathcal V_\Gamma|a,b,0\rangle\\
&=|a,a,b,b\rangle_m|a,b,0\rangle_{ABC},\\
G_vW_0&=W_0,\qquad W_0^\dagger W_0=I_{AB}.
\end{aligned}
\tag{4}
$$

这恰好是两个旧物理对象的联合态加一条独立空桥，可取任意 ρ_ABR，包括独立主体分别与外部参考关联的情况，也包括更一般的旧联合关联态。旧主体的逻辑维数均为2，联合为4；加桥后的物理维数为8，真实运动学资源是7个比特。

“加入”在本轮指**已有中性桥资源上启用跨主体交互**，不指从无创建一个 Hilbert 因子。t=0 前的寄存器准备和硬件连接是输入。在 Z_C=+1 的准备中，扩展 Gauss 条件约化为两个旧主体的原条件；这和362轮把不相容裸输入改写成另一约束的编码不同。

同一构造可以接任意两张有限旧图：桥连接一个旧顶点 u 与另一个旧顶点 v，选择分别入射 u、v 的旧逻辑边 A、B；其他旧边作为旁观者全部保留。没有要求旧未知态为计算基态或纯态。

## 5. 三类严格允许的连续生成元

定义桥的逻辑 Hadamard、相邻链路控制投影及桥—右主体相位生成元：

$$
\mathsf H_C=\frac{\overline X_C+\overline Z_C}{\sqrt2},\qquad
D_{AC}=\frac{(I-\overline Z_A)(I-\overline Z_C)}4,\qquad
Q_{CB}=\overline Z_C\overline Z_B,\qquad
\mathsf H_C^2=Q_{CB}^2=I,\quad D_{AC}^2=D_{AC}.
\tag{5}
$$

全部与每个 G_v 对易。取正强度 g_H、g_D、g_P，允许 Hamiltonian 为 g_H(I−ℋ_C)、g_D D_AC、sign(θ)g_P Q_CB：

$$
\begin{aligned}
e^{-i\frac\pi2(I-\mathsf H_C)}&=\mathsf H_C,&
e^{-i\pi D_{AC}}&=\mathsf{CZ}_{AC},&
e^{-i\theta Q_{CB}}&=\mathsf R_{CB}(\theta),\\
\tau_H&=\frac{\pi}{2g_H},&
\tau_D&=\frac{\pi}{g_D},&
\tau_P&=\frac{|\theta|}{g_P},\\
[H(t),G_v]&=0
&\Longrightarrow\quad U(t)P_{\rm phys}&=P_{\rm phys}U(t)
\quad\text{对每个中间时刻}.
\end{aligned}
\tag{6}
$$

最后一句来自传播方程或分段指数，不依赖只在门结束处投影回物理空间。θ=0 时中间相位步骤可以是零时长恒等操作。

第一类生成元的非恒等 Pauli 项最多为三体：X_m1 X_C X_m2 及 Z_C；第二类只用 I、Z_A、Z_C、Z_AZ_C，第三类用 Z_CZ_B。支集分别是桥的两端及桥本身、共享左顶点的两条链路、共享右顶点的两条链路。这是**给定关系图上的访问范围**，没有由此获得物理距离或自然空间局域性。

本轮把所列项作为实际演化生成元；不声称在任意未知旧 Hamiltonian 同时不可控作用时，这段门串仍自动完成。

## 6. 旧主体之间真正发生相互作用

令 ℂ_AC=ℋ_C CZ_AC ℋ_C，它是旧 A 控制、桥 C 为目标的逻辑 CNOT。按时间先后执行 ℋ_C、CZ_AC、ℋ_C、R_CB、ℋ_C、CZ_AC、ℋ_C，共七步：

$$
\begin{aligned}
\mathsf U_\theta
&=\mathsf C_{AC}\mathsf R_{CB}(\theta)\mathsf C_{AC}
=e^{-i\theta\overline Z_A\overline Z_B\overline Z_C},\\
\mathsf U_\theta W_0
&=W_0U_{AB}(\theta),\qquad
U_{AB}(\theta)=e^{-i\theta Z_AZ_B}.
\end{aligned}
\tag{7}
$$

证明只用 CNOT 将 Z_C 共轭为 Z_A Z_C，而 Z_B 不变。第二行在整个旧输入空间上成立；桥回到0，不是仅对某些输入碰巧回零。

对任意有限 R、任意 ρ_ABR，

$$
\begin{aligned}
(\mathsf U_\theta W_0\otimes I_R)\rho
(\mathsf U_\theta W_0\otimes I_R)^\dagger
&=(W_0\otimes I_R)
(U_{AB}\otimes I_R)\rho(U_{AB}^\dagger\otimes I_R)
(W_0^\dagger\otimes I_R),\\
(W_0^\dagger\mathsf U_\theta^\dagger\otimes I_R)
\,\rho_{\rm out}\,
(\mathsf U_\theta W_0\otimes I_R)&=\rho,\qquad
(\rho_{\rm out})_R=\rho_R .
\end{aligned}
\tag{8}
$$

信息保持指已知联合酉作用后的完全可逆性与参考保持；真实交互当然可以改变两个旧边缘，不能同时要求它们逐一完全不变。

例如旧输入 |+⟩_A|+⟩_B 经该作用后的 A 约化纯度及一个可控响应为

$$
\operatorname{tr}[(\rho_A^{\rm out})^2]
=\frac{1+\cos^2(2\theta)}2,\qquad
\langle\overline Y_A\rangle_{\,
|+\rangle_A|b\rangle_B}
=(-1)^b\sin(2\theta).
\tag{9}
$$

θ=π/4 时两旧逻辑比特形成1 bit纠缠熵；B 的0／1输入使 A 的 Y 期望为+1／−1。它不是只准备了桥自身的相干：单独作用桥的逻辑旋转时，旧逻辑联合态完全不变，代码另行对照。

原边的逻辑 X、Y、Z 在加桥后仍与全部 G_v 对易，物理支集也仍是各自原边及原端点。逻辑仪器可通过这些完整矩阵代数表示；代码核对一个振幅阻尼 instrument 的各未归一化分支、完备性及任意参考样本。此处只证明可表示和指定实现，不宣称任意复杂仪器都由固定有限装置零成本执行。

## 7. 中间桥不能随便丢弃，但可以相干回收

三步完成 ℂ_AC 后，对旧输入 |+0⟩，在逻辑表示中状态为

$$
\frac{|0,0,0\rangle+|1,0,1\rangle}{\sqrt2}_{ABC}.
\qquad
S(\rho_C)=1,\qquad
D(\rho_{AB},|+0\rangle\langle+0|)=\frac12 .
\tag{10}
$$

若只把桥逻辑自由度追踪掉，A 的未知相位便失去可恢复性。若实际把裸桥链路直接重置到0而不补偿物质荷，该例的固定物理投影成功权重只剩1/2；后续保持 G 的酉也不能将剩余不相容扇区免费改回。

七步协议没有这样重置。相位步骤后再次执行 CNOT，将暂存的相关信息相干撤回旧逻辑对象，并严格恢复独立空桥。因此不违背362轮的记忆下界：本轮保留了两个旧逻辑因子，输出没有被强制压进一个二维裸匹配空间。

需要区分**裸桥链路**与**带衣桥逻辑因子**的熵。第一步 ℋ_C 后，桥逻辑态是纯 |+⟩，但裸链路与端点物质已关联，其边缘可以有1 bit熵；这两个子系统定义不同。代码和结果分别记录两者。

桥的初态条件有实际作用：若采用已经满足 Gauss 的逻辑桥 |1⟩，同一七门在旧主体上实现 U_AB(−θ)，而不是 U_AB(θ)。因此协议没有宣称任意未准备的桥也可无差别使用。

## 8. 操作时间、作用支集与两体权限边界

上述外部调度版本的资源如下；θ=0 时直接省略相位步骤，下列相位生成元范数针对非零步骤：

$$
\begin{aligned}
\|H_H\|&=2g_H,\qquad
\|H_D\|=g_D,\qquad
\|H_P\|=g_P,\\
T_{\rm pulses}
&=\frac{2\pi}{g_H}+\frac{2\pi}{g_D}+\frac{|\theta|}{g_P},\\
\int_0^{T_{\rm pulses}}\|H(t)\|\,dt&=6\pi+|\theta|.
\end{aligned}
\tag{11}
$$

这是这套明确七脉冲的账，不是最短时间或最低能耗定理。范数控制不等于热力学做功；逐门切换暂时由外部调度提供，下一节把执行次序放入内部时钟。

另有一个精确的控制限制。对任意 Pauli 串，令 f_m、f_e 标明哪些物质／链路位置含 X 或 Y。它与所有固定 G 对易，当且仅当

$$
f_m+Bf_e=0\pmod2 .
\tag{12}
$$

在无自环、无重边的简单图上，若总 Pauli 支集最多为2，该式只能有 f_m=f_e=0：

- 只翻转物质时，至少一个对应顶点约束被翻转。
- 翻转一条链路必须同时翻转其两个端点物质，至少需要3处支集。
- 翻转两条不同链路时，其边界不能为空；两边若要共同端点完全抵消就必须是被排除的平行重边。

Pauli 基是全部 Gauss 共轭作用的共同本征基，所以不同违规串的线性组合不能相互抵消这个条件。故所有**算符层面严格保持约束的至多两体 Hamiltonian**在本模型中只含 I/Z 串，不能产生桥逻辑翻转。

这并不排除两体相位相互作用，更不证明三体权限是所有加入协议的必要条件。例如在当前最小图物理空间内，

$$
Z_{m_1}Z_{m_2}\mathcal V_\Gamma
=\mathcal V_\Gamma Z_AZ_B,\qquad
e^{-i\theta Z_{m_1}Z_{m_2}}W_0=W_0U_{AB}(\theta).
\tag{13}
$$

因此若允许该端点—端点耦合，只用一个两体相位脉冲就能实现旧主体交互。七步构造的价值是额外展示**桥实际参与、临时记录、全部参考保持及记录回收**，不是追求最优门数。采用有虚跃迁的两体微观装置近似生成有效三体项，也不属于这里“每时刻精确对易的两体总 Hamiltonian”前提。

## 9. 复用346轮：将执行次序放入固定内部 Hamiltonian

取上一节七个完整物理门按时间先后为 U_1,…,U_7。加入一个中性的8能级时钟 K，扩展约束为 I_K⊗G_v。对给定 θ，固定门串可以直接编入耦合，不另用程序寄存器：

$$
H_F=\sum_{j=0}^{6}\gamma_j
\left(|j+1\rangle\langle j|_K\otimes U_{j+1}
+|j\rangle\langle j+1|_K\otimes U_{j+1}^\dagger\right),
\qquad
[H_F,I_K\otimes G_v]=0 .
\tag{14}
$$

对易性逐项来自每个完整门与 G_v 对易。因此**整个自主演化**保持同一物理部门，无需在中途施加投影或外部切换。它是另一种实现：时钟态叠加不同已完成前缀，不等于逐脉冲版本在每个中间时刻具有同一数据态。

令 C_0=I、C_j=U_j⋯U_1、𝒯=Σ_j|j⟩⟨j|⊗C_j。直接计算每个相邻块：

$$
H_F=\mathcal T(h_K\otimes I)\mathcal T^\dagger,\qquad
h_K=\sum_{j=0}^{6}\gamma_j
(|j+1\rangle\langle j|+|j\rangle\langle j+1|),\qquad
e^{-itH_F}=\mathcal T(e^{-ith_K}\otimes I)\mathcal T^\dagger .
\tag{15}
$$

这复用346轮完整空间的酉等价，不是重新提出历史时钟。取工程化权重

$$
\gamma_j=\frac{\Omega}{2}\sqrt{(j+1)(7-j)},\qquad
h_K=\Omega J_x\quad(J=7/2),\qquad
T=\frac\pi\Omega,\qquad
e^{-iTH_F}(|0\rangle_K\otimes W_0)
=(-i)^7|7\rangle_K\otimes W_0U_{AB}(\theta).
\tag{16}
$$

裸链的端点幅度是 [−i sin(Ωt/2)]^7，故 T 时模长严格为1。由算子恒等式，对任意旧输入和参考，

$$
|0\rangle\langle0|_K\otimes
(W_0\otimes I_R)\rho(W_0^\dagger\otimes I_R)
\ \longmapsto\
|7\rangle\langle7|_K\otimes
(W_0U_{AB}\otimes I_R)\rho
(U_{AB}^\dagger W_0^\dagger\otimes I_R).
\tag{17}
$$

这是指定时刻的确定性传输，不靠挑选“时钟恰好到达”的成功子样本。改变 θ 通常改变编入的 H_F；本轮没有构造一台有限处理器，免费精确执行任意连续未知程序。

### 9.1 内部时钟的资源不能隐藏

时钟是一个抽象8能级系统。运动学总维数8×128=1024，固定物理部门维数8×8=64；后者可用于无损计算压缩。将8能级系统记成三比特，只说明信息容量，不能免费取得特定空间局域耦合。

若单个时钟跃迁强度限制为 γ_max，则本七步链有

$$
\max_j\gamma_j=2\Omega\le\gamma_{\max},\qquad
T\ge\frac{2\pi}{\gamma_{\max}},\qquad
\|H_F\|=\frac{7\Omega}{2},\qquad
\langle H_F\rangle_{\rm in}=0,\qquad
\langle H_F^2\rangle_{\rm in}=\frac{7\Omega^2}{4}.
\tag{18}
$$

等式资源在 Ω=γ_max/2 时达到。γ_max=1 给 T=2π、‖H_F‖=1.75、初始能量方差0.4375。时间不变使这些总能量矩保持；均值为零不是零制备功或零能耗的证明。

H_F 的每个跃迁还要联合访问时钟与相应数据门。它**不继承逐脉冲版本的至多三体物理 Hamiltonian 条件**，也没有自动成为某个三维邻近布局。两个实现的范数、时间属于不同耦合资源模型，不能把它们的数字直接当成同一硬件的加速比较。

### 9.2 自主完成不是永久停机

同一精确解给出

$$
p_7(t)=\sin^{14}(\Omega t/2),\qquad
p_7(T+\delta)=\cos^{14}(\Omega\delta/2),\qquad
e^{-i(2T)H_F}=-I .
\tag{19}
$$

T/2 时终点概率仅1/128；忽略时钟得到不同计算前缀的混合，不能无条件称作已完成目标操作。T 后还会回流，在2T返回原输入。永久保存输出、停止耦合、重新准备时钟和装置复用都没有由本轮免费提供。

因此内部化的是**给定准备之后的有限执行区间**。它关闭了“保持 Gauss 必须外部逐门切换”的疑问，没有关闭整体资源从准备到存档的全生命周期问题。

## 10. 可复算结果与实现检查

| 检验 | 数值结果／解析值 |
|---|---:|
| 七脉冲对旧全部输入的终点等距误差，θ=π/4 | 9.33×10⁻¹⁶ |
| 三维参考混态，实际与目标输出的迹距离 | 4.10×10⁻¹⁶ 以下 |
| 两旧加态输入的最终逻辑纠缠熵 | 1 bit |
| B 的0／1引起旧 A 的 Y 期望，θ=π/4 | +1／−1 |
| 第一段 CNOT 后桥逻辑记录熵 | 1 bit |
| 同时刻直接重置裸桥后的物理权重 | 1/2 |
| 完成七门后独立空桥的熵 | 约0，浮点约1.3×10⁻¹⁵ |
| 7比特中权重≤2的 Pauli 串／严格允许者 | 211／29；允许者全为 I/Z |
| 内部时钟的物理部门全部输入传输误差 | 1.12×10⁻¹⁴ |
| 2T 回流的64维算子误差 | 6.12×10⁻¹⁴ 以下 |
| γ_max=1 时 T+0.1 的终点概率 | 0.99563410 |

同一参考测试没有通过有限样本证明“所有 R”；任意参考的量词来自式 (8)、(17) 的完整等距关系。中间时刻数值投影泄漏为0来自此矩阵实现的块结构；一般全时间结论由精确对易证明承担。

代码先通过检查，再以 --write-results 保存。21项检查包括：

1. 不同有限图的物理维数、编码与投影。
2. 完整逻辑 Pauli 及全部 Gauss 交织。
3. 独立旧物理输入加空桥的初始相容。
4. 每类完整生成元的固定约束对易。
5. 任意脉冲角的闭式指数与独立谱指数。
6. 多个 θ 下完整运动学空间的七门恒等式。
7. 每段多个中间时刻对全部物理输入的保持。
8. 混态参考、确定性终点与逆恢复。
9. 真正的旧主体纠缠及仅旋转桥的对照。
10. 一个旧主体对另一个旧主体的可测响应。
11. 临时桥记录、丢弃和非法重置反例。
12. 非空逻辑桥的相反相位操作。
13. 加入后的 instrument 各分支、完备性和参考。
14. 两体严格规范不变 Pauli 的完整枚举。
15. 脉冲支集与实际生成元范数。
16. 完整1024维历史 Hamiltonian 的约束对易及物理部门交织。
17. 64维内部时钟的独立酉等价与全谱。
18. 指定时刻的全部输入算子传输和参考。
19. 自主中间状态的 Gauss 保持与能量矩。
20. 计时、不到达概率、固定跃迁预算及回流。
21. 两体权限障碍不排除所有旧主体交互的独立反例。

没有安装依赖，没有将64维部门模拟冒称为未知完整理论的普遍实验验证。

## 11. 本轮完成范围与后继问题

本轮关闭的具体接口是：

**已有相容物理对象＋明确中性桥资源＋允许的规范不变交互，可以全程保持固定约束，进行真实联合量子操作，并在保留任意参考的同时回收桥记录；有限执行次序也可由内部时钟承担。**

这与362轮的末态表示变换不同，也不反驳其裸匹配记忆下界；对象、约束和保留的逻辑维数已全部重新交代。

仍需推进的真正问题是：如何让准备、到达读取、持久记录及后续调用也在一个明示整体中完成；或者检验严格三体逻辑控制怎样在有限误差和资源预算下由更受限的相互作用实现。可以直接对接已有自主计算、规范保护与器件控制成果。若进一步要求特定几何局域性，必须先给时钟和各个寄存器的物理布局，不能将这里的关系图或计算进度图直接当成已生成的时空。

本轮从给定 Z₂ 结构构造相容协议，没有让“可以执行规范不变加入”自动选择规范结构，更没有推出引力约束代数或 Einstein 方程。

## 12. 文件与复算

- [代码：fixed_gauss_admission_audit.py](365/fixed_gauss_admission_audit.py)
- [结果：fixed_gauss_admission_audit_results.json](365/fixed_gauss_admission_audit_results.json)
- [本轮整合检查](365/research_round_365_checks.json)
- [基线362](research_note_362.md)
- [复用346自主时钟](research_note_346.md)
- [346代码](346/autonomous_history_clock_audit.py)

在本目录用已有运行时：

    & 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 fixed_gauss_admission_audit.py

当前保存为 Python 3.12.14、NumPy 2.3.5 下21项全部通过的结果。代码及结果没有修改既有历史证据。
