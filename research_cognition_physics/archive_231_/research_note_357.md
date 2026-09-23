# 第357轮：非对易自运动能否支持动态源与受控双向经典极限

日期：2026-09-23。先复读项目 README、研究目录三个导航、第353、354、356轮及已有结果，检索未发现相同的动态双部门证明。本轮接续已冻结[356轮](research_note_356.md)，独占新增代码、结果和本笔记；无图像检查。

## 1. 问题与结论

第356轮给出了小集体涨落与有限双向响应共存的正例，但其两个 Z 源守恒，改变的只是 X/Y 分量。本轮加入**固定、不对易的内部运动**，检验源本身随时间变化时，能否继续从同一联合量子规则得到受控的双向经典动力学。

得到条件性正答案：

- 选定纯 iid 初态族后，两个部门的 Z 源可真正变化，且各自的加速度依赖另一部门。
- 集体均方偏离有显式上界 4 exp(4|g|t)/N，均值偏差还有 O(1/N) 的单独上界；固定有限时间的集中性与双向响应同时成立。
- 证明来自完整量子演化的涨落恒等式及 Grönwall 不等式，不靠数值拟合或直接截断关联层级。
- Hamiltonian、全连接耦合、初态及 1/N 缩放仍是新增输入。所得是内部集体动力学，不是空间局域几何或 Einstein 方程。

本轮证明的是指定态族的集体／单体极限，**不继承第356轮对任意未知观察块和参考的 diamond 通道结论**。后者在加入不对易项之后须另证。

## 2. 对接已有定理，并把前提逐项说清

Carollo–Lesanovsky 的 [*Applicability of Mean-Field Theory for Time-Dependent Open Quantum Systems with Infinite-Range Interactions*，PRL 133, 150401 (2024)](https://arxiv.org/html/2403.17163) 式 (1)、式 (4) 和 Theorem 1 已对有限 d、全连接二次集体生成元给出均方误差的 Grönwall 控制。原文明确要求合法完全正传播、适当系数正则性，以及初始集体涨落趋零。

本模型逐项嵌入：

| 原文前提 | 本轮对应 |
|---|---|
| N 个有限 d 单元 | 把第 a 个 G、M 比特作为一个 d=4 单元 |
| 线性和 1/N 二次集体项 | X_G、X_M 及 Z_G Z_M 集体项 |
| 系数正则、N 无关 | Ω_G、Ω_M、g 是常数 |
| 合法完全正传播 | 无耗散，固定自伴 H 生成酉群 |
| 初始均方误差趋零 | 每个单元为同一纯乘积态，涨落 O(1/N) |

选择 Hilbert–Schmidt 正交基时，可取二比特 Pauli 张量积除以 2，连同单位基。于是原文集体变量是本轮相应总和的一半，线性系数为 2Ω，两个对称交叉二次系数各为 2g，仍为常数。配对只服务于定理记账，不建立物理空间邻接。

因此这里不是在质疑已有平均场定理是否成立，而是将其接入项目，并给出本模型更直接的显式常数、源变化见证和资源范围。下面的证明不依赖引用原文未展开的常数。

## 3. 模型、时间和初态输入

记两个部门分别有 N 个量子比特，ℏ=1。定义

$$
\begin{aligned}
S_{\nu,i}&=\sum_{a=1}^{N}\sigma_{\nu,i}^{(a)},\qquad
\widehat s_{\nu,i}=S_{\nu,i}/N,\qquad \nu=G,M,\\
H_N&=\Omega_G S_{G,x}+\Omega_M S_{M,x}
+\frac gN S_{G,z}S_{M,z},\\
|\Psi_N(0)\rangle&=
|\boldsymbol n_G\rangle^{\otimes N}
\otimes|\boldsymbol n_M\rangle^{\otimes N},
\qquad |\boldsymbol n_G|=|\boldsymbol n_M|=1 .
\end{aligned}
\tag{1}
$$

Ω_G、Ω_M、g 是实常数。Hamiltonian 不依赖输入态，也不随时间调度；加入的 X 项是内部固定生成元。态制备、共同轴的识别和这些耦合常数仍是输入，不由 F+U+C+P 自动选择。

候选经典轨道从相同初始 Bloch 向量出发：

$$
\begin{aligned}
\dot{\boldsymbol m}_G
&=2(\Omega_G\boldsymbol e_x+g m_{M,z}\boldsymbol e_z)
\times\boldsymbol m_G,\\
\dot{\boldsymbol m}_M
&=2(\Omega_M\boldsymbol e_x+g m_{G,z}\boldsymbol e_z)
\times\boldsymbol m_M,\\
|\boldsymbol m_\nu(t)|&=1,\qquad
h(\boldsymbol m_G,\boldsymbol m_M)
=\Omega_Gm_{G,x}+\Omega_Mm_{M,x}+g m_{G,z}m_{M,z}.
\end{aligned}
\tag{2}
$$

该光滑常微分方程在两个单位球面上有全时间唯一解，且 h 守恒；长度不变来自叉积与自身正交。它可以用已有自旋 Poisson 括号表示，但这两个球面是内部相空间，不是物理空间坐标。

## 4. 完整量子方程：源确实不再守恒

以下采用 Heisenberg 图景，并在固定初态中取期望。同一时刻两个部门的算子仍互相对易，因为共同酉共轭保持初始交换关系。精确方程为

$$
\begin{aligned}
\dot{\widehat{\boldsymbol s}}_G
&=2\Omega_G\boldsymbol e_x\times\widehat{\boldsymbol s}_G
+2g\widehat s_{M,z}\boldsymbol e_z\times\widehat{\boldsymbol s}_G,\\
\dot{\widehat{\boldsymbol s}}_M
&=2\Omega_M\boldsymbol e_x\times\widehat{\boldsymbol s}_M
+2g\widehat s_{G,z}\boldsymbol e_z\times\widehat{\boldsymbol s}_M,\\
\dot{\widehat s}_{G,z}&=2\Omega_G\widehat s_{G,y},\qquad
\dot{\widehat s}_{M,z}=2\Omega_M\widehat s_{M,y}.
\end{aligned}
\tag{3}
$$

这里没有将任何乘积期望提前因子化。相互作用先改变不对易分量，再经内部 X 运动改变 Z 源；“源可变”与“任意加入一个态依赖量子生成元”完全不同，底层生成元仍是同一个线性 H_N。

## 5. 关键恒等式：转动项在均方涨落中抵消

令偏离算子与均方误差为

$$
\boldsymbol\delta_\nu(t)
=\widehat{\boldsymbol s}_\nu(t)-\boldsymbol m_\nu(t)I,\qquad
D_\nu(t)=\sum_i\langle\delta_{\nu,i}(t)^2\rangle,\qquad
E_N(t)=D_G(t)+D_M(t).
\tag{4}
$$

从式 (2)、(3) 相减，得到例如 G 部门的精确式：

$$
\dot{\boldsymbol\delta}_G
=2\Omega_G\boldsymbol e_x\times\boldsymbol\delta_G
+2g\widehat s_{M,z}\boldsymbol e_z\times\boldsymbol\delta_G
+2g\delta_{M,z}\boldsymbol e_z\times\boldsymbol m_G .
\tag{5}
$$

前两项都是转动项。对平方和求导时，各项按反对易乘积出现，反对称叉积使其成对抵消。第二项虽然含算符系数，但该系数与所有 G 偏离分量对易，所以同样精确抵消；无需忽略有限 N 的交换子。

因此

$$
\begin{aligned}
\dot D_G
&=4g\left\langle
\delta_{M,z}\,(\boldsymbol e_z\times\boldsymbol m_G)
\cdot\boldsymbol\delta_G
\right\rangle,\qquad
|\dot D_G|\le4|g|\sqrt{D_GD_M},\\
|\dot D_M|&\le4|g|\sqrt{D_GD_M},\\
\dot E_N&\le8|g|\sqrt{D_GD_M}\le4|g|E_N .
\end{aligned}
\tag{6}
$$

第一条界来自态内积的 Cauchy–Schwarz：沿任意实单位方向的偏离平方期望不超过三个方向的平方和，且 |e_z×m_G|≤1。最后一步为 2√(D_GD_M)≤D_G+D_M。同一推导也控制 |dot E_N|。

纯 iid 初态使每个部门的总集体方差为 2/N，初始偏差为零。于是，对全部 0≤t≤T，

$$
E_N(0)=\frac4N,\qquad
E_N(t)\le \frac{4e^{4|g|t}}N
\le\frac{4e^{4|g|T}}N .
\tag{7}
$$

Ω 不进入指数并不是省略了内部运动：那些运动以等距转动的形式精确抵消。它们仍决定具体经典轨道及实际误差。

## 6. 均值偏差还可单独控制到 O(1/N)

记 b_ν=〈s_ν〉−m_ν。把式 (5) 中的 s_M,z 拆为 m_M,z+δ_M,z，可得

$$
\begin{aligned}
\dot{\boldsymbol b}_G
&=2(\Omega_G\boldsymbol e_x+gm_{M,z}\boldsymbol e_z)
\times\boldsymbol b_G
+2g b_{M,z}\boldsymbol e_z\times\boldsymbol m_G
+2g\boldsymbol e_z\times\boldsymbol C_G,\\
\boldsymbol C_G&=\langle\delta_{M,z}\boldsymbol\delta_G\rangle,\qquad
|\boldsymbol C_G|\le\sqrt{D_GD_M},\\
u(t)&=|\boldsymbol b_G|+|\boldsymbol b_M|,\qquad
D^+u(t)\le2|g|u(t)+2|g|E_N(t).
\end{aligned}
\tag{8}
$$

D^+ 表示上右导数，包含 b=0 时普通范数导数可能不定义的情形。转动项与 b 正交，其余项使用三角不等式；两侧关联余项总和最多为 4|g|√(D_GD_M)≤2|g|E_N。

由 u(0)=0 和式 (7) 积分：

$$
u(t)\le
\frac4N\left(e^{4|g|t}-e^{2|g|t}\right).
\tag{9}
$$

g=0 时右侧为零，两个独立自旋旋转确实与经典单体运动完全一致。这个界是保守充分界，未声称最优时间增长率。

由偏差—方差恒等式以及各部门的交换对称性，进一步有

$$
\begin{aligned}
\sum_{\nu,i}\operatorname{Var}(\widehat s_{\nu,i})
&=E_N-\sum_\nu|\boldsymbol b_\nu|^2\le E_N,\\
D(\rho_{G,1}(t),\gamma_G^{\rm cl}(t))
+D(\rho_{M,1}(t),\gamma_M^{\rm cl}(t))
&=\frac{u(t)}2 .
\end{aligned}
\tag{10}
$$

第二行仅用单比特迹距离等于 Bloch 向量距离的一半；ρ_ν,1 是本指定大系统态的单体约化态。这不是对任意未知输入或任意参考的通道范数定理。集体单个分量的测量概率可由方差和 Chebyshev 界控制，也不代表三个不对易分量可以联合无扰动精确读取。

## 7. 双向动态源的可检验证据

单纯观察 Z 随自运动变化还不够，还要看另一部门是否改变它。从式 (2) 得到

$$
\ddot m_{G,z}
=4\Omega_Gg\,m_{M,z}m_{G,x}-4\Omega_G^2m_{G,z},
\qquad
\ddot m_{M,z}
=4\Omega_Mg\,m_{G,z}m_{M,x}-4\Omega_M^2m_{M,z}.
\tag{11}
$$

交叉项给出明确的双向源依赖；当 g、Ω 及相应分量非零时，它一般不会消失。初始纯乘积态上，完整量子双重交换子的期望恰好给出同一初始加速度。代码同时核对量子双重交换子与经典 ODE，避免只在经典方程中人为声明“存在反馈”。

采用 Ω_G=0.7、Ω_M=−0.45、g=0.35，初态分别为 (z_G,φ_G)=(0.3,0.25)、(z_M,φ_M)=(−0.4,−0.3)，在 t=0.9：

| 量 | G 部门 | M 部门 |
|---|---:|---:|
| 初始 Z | 0.3 | −0.4 |
| 耦合经典轨道 Z | 0.2144593154 | −0.1518107671 |
| 去掉 g 后的 Z | 0.3164463291 | −0.0796271022 |
| 相互作用造成的 Z 差 | −0.1019870137 | −0.0721836649 |

因此两侧的源变化既有内部转动，也有有限的另一部门影响。其数学形式仍是两个内部自旋的非线性平均场方程，不是应力张量驱动度规的引力方程。

## 8. 状态族、时间窗口和全局量词

式 (6) 的代数不等式本身对任何初始量子态都成立，只要比较的经典初值已指定；但是收敛需要初始集中性：

$$
E_N(t)\le e^{4|g|t}E_N(0).
\qquad
\text{本纯 iid 族中，若 }
e^{4|g|T_N}/N\longrightarrow0,
\text{ 则区间 }[0,T_N]\text{ 上均方误差趋零。}
\tag{12}
$$

这给出一种充分的增长时间条件，**没有证明到达此条件之外就必然失效**，也没有说对数时间是普适最优界。长程相关 cat 态等可能不满足 E_N(0)→0；第353、354轮已说明为何不能以规模大代替初始集中性。

两个初始 Bloch 向量变化时，比较的是不同的 N 份已准备态。准备 γ→γ^⊗N 不是从一个未知副本免费获得 N 份，也不是单副本的仿射量子操作。因而有效非线性 ODE 与底层线性酉演化相容，不重犯第350轮的准备／分支混淆。

本轮未证明：

- 任意观察块、任意参考及任意初始部门间关联的同一 diamond 近似；
- 完整 2N 比特量子态接近经典轨道的乘积态；
- 任意长时间的统一误差，或热化、不可逆耗散；
- 全部允许态都集中到一个经典点。

有限块或全局的更强结论必须给相应比较对象和范数，不能直接由六个集体均值推断。

## 9. 资源与数值方法的误差账

资源保持为两个实际增长的部门，另加每比特固定内部运动：

$$
\dim\mathcal H_{GM}=4^N,\qquad
\dim\mathcal H_{\rm sym}=(N+1)^2,\qquad
\|H_N\|\le N(|\Omega_G|+|\Omega_M|+|g|).
\tag{13}
$$

交叉相互作用仍有 N² 对，每对 g/N，每个比特承担的交叉耦合绝对值之和为 |g|。每粒子能量尺度有限，但没有固定维数的空间局域性或有限传播光锥。对称部门压缩减少计算量，不减少实际 2N 个量子比特的资源。

数值使用三条互相核对的实现：

1. N≤3 时，构造完整 4^N 维张量积 Hamiltonian 并直接谱分解演化。
2. 较大 N 使用两个最大自旋部门的 (N+1)×(N+1) 状态系数表；不截断该部门内的态。Hamiltonian 对状态的作用由三对角自运动与对角交叉相互作用直接计算。
3. 经典六维 ODE 用 RK4 积分，另以时间步减半核查误差。

第二条的指数作用采用每子步 18 阶 Taylor 多项式。若 x=‖HΔt‖ 的已知上界不超过 1/2，p 为阶数、L 为子步数，则精确算术中

$$
\epsilon_{\rm step}\le e^x\frac{x^{p+1}}{(p+1)!},
\qquad
\|P_p(-iH\Delta t)^L-e^{-iHL\Delta t}\|
\le L\epsilon_{\rm step}e^{L\epsilon_{\rm step}} .
\tag{14}
$$

第二条界由幂的伸缩求和和 ‖P_p‖≤1+ε_step 得到。这只是**截断误差界**，没有把浮点舍入误差也严格认证。实际另查范数、能量、反向演化、子步减半与提高阶数，并用小 N 完整张量计算比对。

N=64、t=0.9 时，173 个子步的精确算术截断上界约 4.37×10^−21；实际范数误差约 1.33×10^−15，能量每 N 的误差约 2.78×10^−16，符合浮点误差占主导的预期。

## 10. 结果与复算

上述参数、固定 t=0.9：

| 每部门 N | 两侧均值偏差和 u | 式 (9) 上界 | 总集体方差 | E_N | 式 (7) 上界 |
|---:|---:|---:|---:|---:|---:|
| 8 | 0.0411536621 | 0.823905454 | 0.578345119 | 0.579191970 | 1.762710744 |
| 16 | 0.0206925246 | 0.411952727 | 0.289621003 | 0.289835102 | 0.881355372 |
| 32 | 0.0103771226 | 0.205976364 | 0.144927037 | 0.144980882 | 0.440677686 |
| 64 | 0.00519652424 | 0.102988182 | 0.0724932264 | 0.0725067288 | 0.220338843 |

理论界保守但明确；1/N 结论由证明给出，不由这四个规模的数值趋势代替。

13 项检查覆盖：完整张量与对称演化；初始均方预算；非守恒源交换子；涨落抵消恒等式及其随机相关态核验；有限时间均方和均值界；经典长度及能量；积分细化和截断界；量子能量及逆演化；两侧 Z 对耦合的响应；g=0 精确解；Ω=0 回到356轮；量子与经典初始源加速度；生成元自伴性及范数预算。

[代码](dynamic_two_sector_limit_audit.py)、[结果](dynamic_two_sector_limit_audit_results.json)、[独立核验](research_round_357_checks.json)。

在仓库根目录运行：

    & 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 research_cognition_physics/archive_231_/dynamic_two_sector_limit_audit.py
    & 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 research_cognition_physics/archive_231_/dynamic_two_sector_limit_audit.py --write-results

复用冻结354轮的相干态及小 N 张量工具；没有修改它们。运行环境为既有 Python 3.12.14 与 NumPy 2.3.5，无新安装。

## 11. 对主线关闭了什么，仍缺什么

**关闭的限定缺口：** 第356轮的“源静止”并不是双向经典极限必不可少的条件。给固定联合 H 加上明示的不对易内部运动，仍能严格控制同一初态族的动态源、有限响应和小涨落。

**未关闭的生成缺口：** 本轮仍没有从认知合同选出这个 H、乘积初态、两部门划分或耦合图，也没有让内部集体量成为空间度规。它不提供 Einstein 约束、普适物质耦合或物理空间维数。

下一项真正有用的接口应转向**有关系／局域结构的多个集体单元**：先明确它们的连接和资源，再检验有效动力学、传播与约束能否相容。需要保留“局域单元内部的大 N”与“物理空间上的位置数”两种不同极限，防止再次把 N 或内部球面当成空间。继续更换全连接耦合或参数本身，已不足以推动几何生成这一缺口。
