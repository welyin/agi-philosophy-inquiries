# 第355轮：自主声学几何为何仍不自动成为Einstein引力

日期：2026-09-23。本轮16项检查。基线为已完成的231—352轮；不依赖本批未冻结的353、354轮。

## 1. 本轮的具体问题与结论

已先读项目README、本目录README、research_direction、RESEARCH_STATE，检索旧笔记中的声学／BEC／Gross–Pitaevskii接口，并回读[344轮](research_note_344.md)、[349轮](research_note_349.md)和[350轮](research_note_350.md)的范围。既有结果已区分“给定曲背景上的量子传播”和“几何方程”；尚未完整核验一个由介质自身动力学产生时变几何的实例。

**本轮得到更具体的分界：** 无外力、参数不随时间调制的排斥型Gross–Pitaevskii（GP）模型，允许一个精确的均匀膨胀背景。其低波数扰动感受到真正有曲率的Lorentz声学度规，度规随背景自主变化。但计算所得Einstein张量既不满足任意Λ的真空Einstein方程，也不能由正耦合常数下、满足该声学度规零向能量条件的源产生。

这并不等于“所有物质源下都不能写Einstein方程”。事后把Einstein张量定义成某种应力，总能得到一个形式源；缺失的是**由独立物质变量、作用量及操作解释给出这个源与几何响应的关系**。尤其不能因为未加入声子扰动，就把形成背景的全部原子物质当作不存在。

| 层次 | 本轮内容 |
|---|---|
| 认知动机 | 可测传播关系随内部状态变化，是否足以产生引力动力学 |
| 额外输入 | 三维实验室空间、实验室时间、GP近似、原子质量与接触相互作用、初始膨胀相位 |
| 解析证明 | 精确GP解；声子二次作用量；声学曲率与源兼容性障碍；低波数适用窗口 |
| 数值验证 | 真实复序参量差分；连续性、Euler与局部能量；直接联络收缩；扰动方程与误差检查 |
| 解释范围 | 自主有效几何与Einstein动力学之间仍有条件；不是完整F＋U＋C＋P的反模型 |

## 2. 直接对接已有声学引力路线

[Unruh（1981）原始论文](https://doi.org/10.1103/PhysRevLett.46.1351)提出流体声学视界中的辐射类比；本轮核验其出版社摘要，未取得收费全文。[Visser（1998）原文](https://arxiv.org/html/gr-qc/9712010)§2给出无黏、正压、无旋流中的声学度规，并明确区分传播几何与Einstein动力学。

[Barceló、Liberati、Visser（2001）原文](https://arxiv.org/html/gr-qc/0011026)式(5)—(7)、§IV—VI与式(77)分别核验了接触耦合、GP近似、低波数几何和Bogoliubov色散。其声学度规与原子所在的实验室度规不同；相位扰动的几何方程不能直接当成原子的引力方程。

另读同作者[《Analogue models for FRW cosmologies》（2003）作者版](https://homepages.ecs.vuw.ac.nz/~visser/Articles/Journals/ijmpd-frw.pdf)§2—3。利用外流或改变声速模拟FRW几何已经是成熟路线。本文新增的是**本项目中的精确无外力分支、实际Einstein张量与独立源条件核验**，不把声学时空或模拟宇宙学当作新发现，也不承诺这种无限均匀分支能直接制成实验装置。

下面采用明确的正相位约定，自行从GP方程推导，不依赖不同文献对相位符号和质量密度记号的选择。

## 3. 给定GP模型中的精确自主背景

实验室坐标为t、x，原子质量m＞0，接触耦合g＞0，外势为零，所有参数固定。n是**原子数密度**，不是质量密度；a_s是散射长度。采用

$$
\begin{aligned}
i\hbar\partial_t\Psi
&=-\frac{\hbar^2}{2m}\nabla^2\Psi+g|\Psi|^2\Psi,
&g&=\frac{4\pi\hbar^2a_s}{m},\\
\Psi&=\sqrt n\,e^{i\phi},
&\mathbf u&=\frac{\hbar}{m}\nabla\phi,
&c^2&=\frac{gn}{m}.
\end{aligned}
\tag{1}
$$

这里Ψ是凝聚体**序参量的均值场变量**。GP方程是在弱相互作用、稀薄气体、凝聚耗尽较小等条件下使用的有效近似；它不是把全部未知量子态改成按非线性方程演化。完整原子场与声子量子涨落、耗尽及其反作用未在本轮求解，不与350轮的线性量子操作边界混用。

取实虚部，得到连续性与Euler方程：

$$
\begin{aligned}
\partial_t n+\nabla\!\cdot(n\mathbf u)&=0,\\
\partial_t\mathbf u+(\mathbf u\!\cdot\nabla)\mathbf u
&=-\frac gm\nabla n
+\frac{\hbar^2}{2m^2}\nabla\!\left(\frac{\nabla^2\sqrt n}{\sqrt n}\right).
\end{aligned}
\tag{2}
$$

选择n₀＞0、H₀＞0，令a＝1＋H₀t＞0。以下整个背景无需外力或时间编程：

$$
\begin{aligned}
a(t)&=1+H_0t,&
H(t)&=\frac{H_0}{a},&
n(t)&=\frac{n_0}{a^3},&
\mathbf u(t,\mathbf x)&=H(t)\mathbf x,\\
\phi(t,\mathbf x)
&=\frac{mH_0|\mathbf x|^2}{2\hbar a}
-\frac{gn_0}{2\hbar H_0}(1-a^{-2}),&
c(t)&=c_0a^{-3/2},&
c_0^2&=\frac{gn_0}{m}.
\end{aligned}
\tag{3}
$$

n在空间中均匀，所以背景量子压**精确为零**；并非先忽略它才得到解。H点＝−H²使流体质点无加速度，n点＝−3Hn满足连续性。更直接地，

$$
\begin{aligned}
\frac{\partial_t\Psi}{\Psi}
&=-\frac32H+i\left(-\frac{mH^2|\mathbf x|^2}{2\hbar}
-\frac{gn}{\hbar}\right),\\
\frac{\nabla^2\Psi}{\Psi}
&=\frac{3imH}{\hbar}-\frac{m^2H^2|\mathbf x|^2}{\hbar^2},\\
\frac{i\hbar\partial_t\Psi}{\Psi}
&=-\frac{3i\hbar H}{2}+\frac{mH^2|\mathbf x|^2}{2}+gn
=-\frac{\hbar^2}{2m}\frac{\nabla^2\Psi}{\Psi}+gn.
\end{aligned}
\tag{4}
$$

这证明它是**所给GP方程的精确解**。H₀趋零时，相位的时间项连续变成−gn₀t／ℏ，恢复静态均匀背景。

### 3.1 内部资源与范围

令实验室压强P＝gn²／2，背景量子梯度能为零。局部能量也有明确账目：

$$
\begin{aligned}
\mathcal E&=\frac12mn|\mathbf u|^2+\frac12gn^2,
&\mathbf J_E&=(\mathcal E+P)\mathbf u,\\
\partial_t\mathcal E+\nabla\!\cdot\mathbf J_E&=0,
&n(t)a(t)^3V_0&=n_0V_0.
\end{aligned}
\tag{5}
$$

膨胀的物质单元保持原子数；其内能降低包含对邻接流体的压力功，不能把有限截取区域称作无边界交换的封闭系统。

模型定义在无限均匀空间，总粒子数无限，Ψ不是归一化的单粒子态。速度在空间无穷远无界，因此实际非相对论原子解释只取速度远小于真空光速的有限局部区域及有限时窗。截成有限凝聚体会引入边缘密度梯度和量子压，必须另解边界问题。本轮**不是有限封闭宇宙的实现证明**。

## 4. 从扰动作用量确定声学度规，保留共形因子

写相位扰动为δφ＝mϑ／ℏ，故速度扰动是∇ϑ；密度扰动为δn。D_t＝∂_t＋u·∇。展开GP作用量至二阶，在本背景上得到

$$
S_2=\int dt\,d^3x\left[
-m\,\delta n\,D_t\vartheta
-\frac g2(\delta n)^2
-\frac{mn}{2}|\nabla\vartheta|^2
-\frac{\hbar^2}{8mn}|\nabla\delta n|^2
\right].
\tag{6}
$$

背景量子压为零，不代表式(6)最后一项为零。只有对足够长波的扰动，才可先忽略该项并消去δn：

$$
\delta n=-\frac mgD_t\vartheta,\qquad
\frac{S_{2,\mathrm{hyd}}}{m}
=\frac12\int dt\,d^3x\,\frac n{c^2}
\left[(D_t\vartheta)^2-c^2|\nabla\vartheta|^2\right].
\tag{7}
$$

完整系数而非仅主光锥确定声学度规。以签名−＋＋＋，

$$
\begin{aligned}
ds_{\mathrm{ac}}^2
&=\frac nc\left[-c^2dt^2+
|d\mathbf x-\mathbf u\,dt|^2\right],\\
\frac{S_{2,\mathrm{hyd}}}{m}
&=-\frac12\int d^4x\,\sqrt{-g_{\mathrm{ac}}}\,
g_{\mathrm{ac}}^{\mu\nu}\partial_\mu\vartheta\partial_\nu\vartheta,\\
\partial_\mu\!\left(\sqrt{-g_{\mathrm{ac}}}
g_{\mathrm{ac}}^{\mu\nu}\partial_\nu\vartheta\right)&=0.
\end{aligned}
\tag{8}
$$

程序实际求逆和行列式，验证式(7)、(8)的系数矩阵一致。因n、c＞0，度规恰有一个负特征值。实验室原子具有不同的本底钟尺；这里的c是声速。

为便于比较，只除以固定常数n₀／c₀，定义归一化声学度规。此常数缩放不改无质量波动方程；时间变化的n／c则不能任意删除。取共动坐标X＝x／a，有

$$
\begin{aligned}
d\bar s_{\mathrm{ac}}^2
&=\frac{n/c}{n_0/c_0}
\left[-c^2dt^2+|d\mathbf x-\mathbf u\,dt|^2\right]\\
&=-c_0^2a^{-9/2}dt^2+a^{1/2}d\mathbf X^2.
\end{aligned}
\tag{9}
$$

声学proper time在此归一化下用长度单位表示，初始dτ＝c₀dt。声学尺度因子记A，不能与原子膨胀尺度a混淆：

$$
\begin{aligned}
d\tau&=c_0a^{-9/4}dt,&
\tau&=\tau_*(1-a^{-5/4}),&
\tau_*&=\frac{4c_0}{5H_0},\\
A(\tau)&=a^{1/4}=(1-\tau/\tau_*)^{-1/5},&
d\bar s_{\mathrm{ac}}^2&=-d\tau^2+A(\tau)^2d\mathbf X^2.
\end{aligned}
\tag{10}
$$

## 5. 计算曲率后，Einstein源项问题才出现

撇号表示对τ求导，声学Hubble量记ℋ。采用R_μν＝∂_αΓ^α_μν−∂_νΓ^α_μα＋Γ^α_αβΓ^β_μν−Γ^α_νβΓ^β_μα的曲率约定。由式(10)：

$$
\begin{aligned}
\mathcal H&=\frac{A'}A=\frac{H_0}{4c_0}a^{5/4},
&\mathcal H'&=5\mathcal H^2,\\
R&=6(\mathcal H'+2\mathcal H^2)=42\mathcal H^2,\\
G_{\hat0\hat0}&=3\mathcal H^2,
&G_{\hat i\hat i}&=-13\mathcal H^2\quad(i=1,2,3).
\end{aligned}
\tag{11}
$$

帽指标为声学度规的正交标架。程序不调用FRW曲率公式：从实验室时间下式(9)的度规、一阶及二阶导数构造Γ、∂Γ和全部Ricci／Einstein张量，再与式(11)比较。

取零向量k̂＝(1,1,0,0)，任何Λ项都在零向收缩中消失：

$$
\begin{aligned}
(G_{\hat\mu\hat\nu}+\Lambda\eta_{\hat\mu\hat\nu})
k^{\hat\mu}k^{\hat\nu}
&=G_{\hat0\hat0}+G_{\hat1\hat1}
=-10\mathcal H^2<0\qquad(H_0\ne0),\\
G_{\mu\nu}+\Lambda g_{\mu\nu}=\kappa T_{\mu\nu},\quad
\kappa>0,\quad T_{\mu\nu}k^\mu k^\nu\ge0
&\quad\Longrightarrow\quad\text{与本背景不相容}.
\end{aligned}
\tag{12}
$$

由此得到两个**限定清楚的否定结果**：

1. 该声学度规不是任意Λ下的真空Einstein度规，甚至允许逐点调Λ也消不掉这个零向收缩。
2. 若另行指定一个在声学度规上满足零向能量条件（NEC）的源，并要求κ＞0，则它不能作为这个度规的Einstein源。

第二条的源条件是**新增检验输入**。本轮没有证明真实原子的非相对论应力就是声学度规上的NEC源；将二者直接等同，反而是被审计的缺口。第一条也不能用“没有声子，所以全部T＝0”来解释：原子凝聚背景本身一直存在。

若事后取Λ＝0、定义T_eff＝G／κ，则形式上当然可以写成Einstein方程，并得到

$$
\rho_{\mathrm{eff}}=\frac{3\mathcal H^2}{\kappa},
\qquad
p_{\mathrm{eff}}=-\frac{13\mathcal H^2}{\kappa},
\qquad
w_{\mathrm{eff}}=-\frac{13}{3},
\qquad
\rho_{\mathrm{eff}}'+3\mathcal H
(\rho_{\mathrm{eff}}+p_{\mathrm{eff}})=0.
\tag{13}
$$

最后的守恒只是Bianchi恒等式的结果。它不提供独立物质作用量，不说明哪种原子观测量就是ρ_eff，也不证明这种形式“幻影流体”真实存在或模型有负能激发。排斥GP模型在实验室中的能量密度与压强均为正，这和式(13)不矛盾：两种应力定义及所用度规不同。

因此，**有内部动力学产生传播几何**尚不能替代**独立物质源与几何响应共同满足Einstein关系**。本轮没有排除其他声学背景偶然满足某些Einstein解，或另加一致物质映射后得到某个有效引力扇区。

## 6. 声子窗口与不可外推的未来端点

对一个共动Fourier波数K，物理实验室波数k＝K／a。令δ＝δn／n。保留式(6)中的量子压，线性化GP给出

$$
\begin{aligned}
\dot\delta&=\frac{K^2}{a^2}\vartheta,\\
\dot\vartheta&=-\left(c_0^2a^{-3}
+\frac{\hbar^2K^2}{4m^2a^2}\right)\delta,\\
\ddot\delta+2H\dot\delta+
\left(c_0^2K^2a^{-5}
+\frac{\hbar^2K^4}{4m^2}a^{-4}\right)\delta&=0.
\end{aligned}
\tag{14}
$$

忽略量子压时，由式(14)得到的ϑ方程与声学曲时空方程完全一致：

$$
\ddot\vartheta+3H\dot\vartheta+c_0^2K^2a^{-5}\vartheta=0
\quad\Longleftrightarrow\quad
\vartheta''+3\mathcal H\vartheta'
+\frac{K^2}{A^2}\vartheta=0.
\tag{15}
$$

这里“声学方程成立”不要求几何光学WKB条件；只有进一步把扰动解释为局部准定常频率、射线或粒子时，才需背景变化慢于该模式。冻结系数时的Bogoliubov频率和两个不同窗口为

$$
\begin{aligned}
\omega_{\mathrm{com}}^2
&=c^2k^2+\frac{\hbar^2k^4}{4m^2},&
\eta&=\frac{\hbar k}{2mc}
=\frac{\hbar K}{2mc_0}\sqrt a,\\
\eta&\ll1
\quad\text{（扰动量子压小）},&
\epsilon_{\mathrm{WKB}}
&=\frac{H}{ck}
=\frac{H_0}{c_0K}a^{3/2}\ll1,\\
\frac{\omega_{\mathrm{com}}}{ck}-1
&=\sqrt{1+\eta^2}-1,&
na_s^3&=n_0a_s^3a^{-3}\ll1.
\end{aligned}
\tag{16}
$$

η与通常愈合长度ξ＝ℏ／(√2mc)的关系是η＝kξ／√2。实验室均匀流的Doppler项为u·k；式(16)写在随流局部参考中，时间变化时它是瞬时频率，而非整个演化中的守恒频率。

稀薄条件随膨胀改善，但固定K的量子压比及WKB参数均增大。若同时要求η≤η_max、ε_WKB≤ε_max，则共同窗口需满足

$$
\frac{H_0a^{3/2}}{c_0\epsilon_{\max}}
\le K\le
\frac{2mc_0\eta_{\max}}{\hbar\sqrt a},
\qquad
\frac{\hbar H_0}{2mc_0^2}a^2
\le\eta_{\max}\epsilon_{\max}.
\tag{17}
$$

所以式(10)在t趋无穷时出现有限τ_*和曲率增长，并不能被外推为真实宇宙的未来奇点：每个固定非零K最终离开声学长波窗口，更强的局部频率解释通常更早失效。高波数下GP的色散也不具有同一个严格最大声速；本轮没有产生所有自由度共用的基本光锥。

## 7. 可复算结果与误差

### 7.1 背景及曲率

使用NumPy直接复算，无SciPy。GP差分检查采用m＝0.7、ℏ＝0.8、g＝0.4、n₀＝3、H₀＝0.2的数学参数，在t＝0.3、x＝(1.2,0.6,−0.4)上直接评价复序参量；五点时间导数和三方向五点Laplacian不调用解析残差。

| 差分步长 | GP残差绝对值除以｜Ψ｜ |
|---|---:|
| 0.04 | 1.74817×10⁻⁶ |
| 0.02 | 1.09218×10⁻⁷ |
| 0.01 | 6.82542×10⁻⁹ |

误差按四阶预期缩小。单位及离散步长在此按明示数值单位取值；这不是实际原子实验的参数拟合。

几何与声子表采用m＝ℏ＝1、g＝0.01、n₀＝100、H₀＝0.002，故c₀＝1、初始气体参数约5.0393×10⁻⁸。

| 原子尺度a | 声学尺度A | 声学proper time τ | 直接联络计算的R |
|---|---:|---:|---:|
| 1 | 1 | 0 | 1.05×10⁻⁵ |
| 1.5 | 1.106682 | 159.0395 | 2.89346×10⁻⁵ |
| 2 | 1.189207 | 231.8207 | 5.93970×10⁻⁵ |

全部张量分量与式(11)相符。表中的曲率针对保留共形因子、固定归一化的声学度规；它不是实验室时空曲率。

### 7.2 扰动比较与守恒检查

程序用s＝H₀t和无量纲线性向量y＝(δ,w)，其中w＝Kϑ／c₀，在0≤s≤1积分：

$$
\frac{d}{ds}
\begin{pmatrix}\delta\\w\end{pmatrix}
=\frac{c_0K}{H_0}
\begin{pmatrix}
0&a^{-2}\\
-a^{-3}-\eta_0^2a^{-2}&0
\end{pmatrix}
\begin{pmatrix}\delta\\w\end{pmatrix},
\qquad
\eta_0=\frac{\hbar K}{2mc_0},
\qquad
\frac d{ds}\operatorname{Im}(\delta^*w)=0.
\tag{18}
$$

水动力近似去掉η₀²项。两边采用同一个初始向量(1,−i)；这只是线性解的归一化，物理扰动可整体乘任意足够小的幅度。守恒量是式(18)的辛Wronskian，不是向量Euclidean范数；该范数差仅是明示的数值比较量，不能称量子态保真度或观测概率误差。

K＝0.1时，a从1到2，η从0.05到0.07071，瞬时频率修正不超过0.24969%，但有限时间累计向量差达到约0.0534534。4096与8192步的完整GP线性方程积分之差约1.104×10⁻⁹，Wronskian最大漂移约1.358×10⁻¹¹。**小的瞬时色散修正不等于任意长时演化误差小。**

K＝0.5时，末态η约0.35355，频率修正约6.0660%；累计向量差约2.26986，而积分加倍差约4.018×10⁻⁶，数值误差显著小于近似差。该对照不在η≪1的严格长波窗口，作用是显示遗漏量子压的实际后果。

RK4用256、512、1024、2048步得到相对4096步结果的差，分别约7.710×10⁻⁵、4.820×10⁻⁶、3.002×10⁻⁷、1.766×10⁻⁸，符合四阶收敛。这是明确ODE上的数值验证，不是无限维GP长时误差的解析上界，也不是完整量子声子产生计算。

本轮16项检查全部通过：精确GP与静态极限、真实复函数差分、连续性／Euler／局部能量、物质体积数守恒、作用量系数与Lorentz签名、坐标拉回、声学钟与尺度、完整曲率张量、Λ零向障碍、Bianchi恒等式、声学模式方程、Bogoliubov频率、各近似窗口、RK4收敛、Wronskian守恒、量子压差异。

## 8. 本轮关闭了什么，下一步缺什么

本轮关闭的捷径是：**“只要传播几何能够由内部状态自主变化，就已经获得Einstein动力学。”** GP提供了自主背景及有效Lorentz传播，却仍由原子流体方程规定几何；声学曲率、真实源映射与普适物质耦合不是同一个结论。

准确的逻辑范围是：本例给出非Einstein真空度规，并给出针对独立声学NEC源的严格不相容见证。它没有以“背景没有声子”为前提证明所有原子T为零，没有排除事后构造应力的形式Einstein表示，也没有核验全部认知合同。因而**完整认知原则到GR的持续目标仍未完成或整体证伪**。

后续应检验有实际源项的几何响应：从底层变量同时取得可测物质应力、可变几何自由度与约束，验证其响应系数和约束是否满足[351轮](research_note_351.md)、[352轮](research_note_352.md)明确列出的条件；或对接诱导引力／量子反作用路线，逐项区分几何动力学项与介质背景项。不能仅继续扩充能模拟哪些已知曲时空。

[代码](acoustic_geometry_dynamics_audit.py)、[结果](acoustic_geometry_dynamics_audit_results.json)、[本轮核验](research_round_355_checks.json)。

    python -B -X utf8 research_cognition_physics/archive_231_/acoustic_geometry_dynamics_audit.py

保存结果加 --write-results；先执行Checks，再写入结果，拒绝覆盖不同旧结果。没有图像检查。
