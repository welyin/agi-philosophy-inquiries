# 第361轮：给定关系图上的集体经典极限与不随位置数增长的传播控制

日期：2026-09-23。先复读项目 README、研究目录 README、research_direction、RESEARCH_STATE，以及已冻结[第357轮](research_note_357.md)和相关代码、结果。本轮科学基线为第357轮，不依赖第358—360轮的新结论。没有图像检查。

## 1. 本轮补上的接口

第357轮处理两个集体部门，证明动态源与双向经典响应可以在同一量子生成元下共存。本轮把集体单元放在**预先给定的有限加权关系图**上，回答此前两个部门不能回答的问题：位置数增加以后，经典误差是否必然失控？初始改变如何沿关系边传播？

得到三个解析结果：

1. 若每个位置的绝对耦合权重之和一致有界，则逐位置集体均方误差及均值偏差有与总位置数 M 无关的上界；固定时间下分别为 O(1/N)。
2. 两条经典轨道的差异受加权邻接矩阵指数控制；距离较远的响应具有路径阶数带来的阶乘尾。它是小尾巴控制，不是严格光锥。
3. 逐位置误差趋零不等于一次测量中整个网络所有位置都准确。一个初始乘积态的精确二项分布反例给出：令 M 随 N 足够快增长，全网至少一处大偏差的概率仍趋于 1−exp(−1)。

这里 N 是**每个位置内部的量子比特数**，M 是**关系图的位置数**。把它们混成同一个“系统规模”，会误读极限。本轮没有选择关系图、空间维数、Lorentz 度规或引力约束代数。

## 2. 已有研究接口及其适用范围

Carollo–Lesanovsky 的 [*Applicability of Mean-Field Theory for Time-Dependent Open Quantum Systems with Infinite-Range Interactions*，PRL 133, 150401 (2024)](https://arxiv.org/html/2403.17163)，式 (1)、式 (4) 和 Theorem 1，对固定有限单元维数、合法完全正传播、满足正则性条件的集体生成元及初始集中态族，给出平均场误差控制。

对本轮的**固定 M**，可将每个位置的第 a 个比特合成一个单元，得到 N 个维数 d=2^M 的单元，图上的二次集体 Hamiltonian 因而可以接入该框架。但 M 增长时，d 及所需算子基也增长；不能直接把固定 d 的常数宣称为对 M 一致。本轮以下证明直接在位置上估计，显式得到不依赖 M 的常数。

Raz–Sims 的 [*Lieb-Robinson Bounds for Classical Anharmonic Lattice Systems* (2009)](https://arxiv.org/html/0902.0025)，特别是其经典格点振子模型与 Theorem 1，说明经典动力学也可以具有体积一致的传播上界。该文对象是带适当势能控制的振子格点及相应可观测量，并非本轮的紧致自旋球面；本轮不套用其速度或常数，而是直接推导下面的图路径界。

| 分类 | 本轮内容 |
|---|---|
| 认知动机 | 大量可操作单元组成关系网络时，局部经典描述与关系响应能否同时成立 |
| 额外建模输入 | 图、权重、每位置 N 个比特、共同 Pauli 轴、Hamiltonian、时间参数、特定初态族 |
| 解析证明 | 逐位置涨落恒等式、比较矩阵、均值偏差、路径传播界、两种规模极限的区别 |
| 数值验证 | 小图完整张量积核验、置换对称部门计算、负权重与不同图对照、守恒及误差检查 |
| 物理解释 | 一种给定关系结构上的集体经典动力学接口；没有把该关系图推导成物理空间 |

## 3. 模型：每个位置有 N 个真实量子比特

取任意有限无向图，位置 v=1,…,M，实对称权重 J，且 J 的对角元为零。每条边在 Hamiltonian 中只计一次，取 ℏ=1：

$$
\begin{aligned}
S_{v,i}&=\sum_{a=1}^{N}\sigma_i^{(v,a)},\qquad
\widehat s_{v,i}=S_{v,i}/N,\\
H_N&=\sum_v\Omega_v S_{v,x}
+\frac1N\sum_{\{v,w\}\in E}J_{vw}S_{v,z}S_{w,z},\\
|\Psi_N(0)\rangle&=\bigotimes_{v=1}^{M}
|\boldsymbol n_v\rangle^{\otimes N},
\qquad |\boldsymbol n_v|=1 .
\end{aligned}
\tag{1}
$$

不同位置可以有不同纯态。每位置制备 N 份同态是给定的制备资源，不是从未知单份状态免费克隆。整个演化由同一个固定自伴 Hamiltonian 生成，没有按输入态重新选择生成元。

定义绝对权重矩阵及入射预算：

$$
A_{vw}=|J_{vw}|,\qquad
\kappa_v=\sum_w A_{vw}\le\kappa,\qquad
|\Omega_v|\le\Omega_* .
\tag{2}
$$

图族改变 M 时要求 κ、Ω_* 可统一选择。图度数可以变化，真正进入估计的是权重之和；例如星形或完全图增加位置时，适当缩小边权可维持相同预算。该缩放本身也是模型输入。

候选经典轨道从 m_v(0)=n_v 出发：

$$
\begin{aligned}
\dot{\boldsymbol m}_v
&=2\left(\Omega_v\boldsymbol e_x+
\sum_wJ_{vw}m_{w,z}\boldsymbol e_z\right)
\times\boldsymbol m_v,\\
|\boldsymbol m_v(t)|&=1,\qquad
h=\sum_v\Omega_vm_{v,x}
+\sum_{\{v,w\}\in E}J_{vw}m_{v,z}m_{w,z}.
\end{aligned}
\tag{3}
$$

每个有限图上，球面的紧致性及向量场光滑性保证全时间唯一解；叉积给出长度守恒，直接代入给出 h 守恒。这些球面是内部相空间，不能把它们的三分量解释成三维物理空间。

## 4. 逐位置误差：局部转动不会放大平方偏离

使用 Heisenberg 图景，并在固定初态中取期望。不同位置的算子在同一时刻仍然对易，因为共同酉共轭保持原先的交换关系。定义

$$
\boldsymbol\delta_v(t)
=\widehat{\boldsymbol s}_v(t)-\boldsymbol m_v(t)I,
\qquad
D_v(t)=\sum_{i=x,y,z}\langle\delta_{v,i}(t)^2\rangle .
\tag{4}
$$

不对量子关联作因子化，可以直接得到精确恒等式：

$$
\begin{aligned}
\dot{\boldsymbol\delta}_v
={}&2\Omega_v\boldsymbol e_x\times\boldsymbol\delta_v
+2\sum_wJ_{vw}\widehat s_{w,z}
  \boldsymbol e_z\times\boldsymbol\delta_v\\
&+2\sum_wJ_{vw}\delta_{w,z}
  \boldsymbol e_z\times\boldsymbol m_v .
\end{aligned}
\tag{5}
$$

在 D_v 的导数中，第一项按反对称性抵消；第二项也抵消，因为其来自其他位置的系数与 δ_v 各分量对易。剩下

$$
\dot D_v
=4\sum_wJ_{vw}
\left\langle\delta_{w,z}
(\boldsymbol e_z\times\boldsymbol m_v)
\cdot\boldsymbol\delta_v\right\rangle .
\tag{6}
$$

对任意态，Cauchy–Schwarz 给出右侧每项的绝对值不大于 4|J_vw|√(D_vD_w)，因为 |e_z×m_v|≤1。再用 2√(D_vD_w)≤D_v+D_w：

$$
\dot D_v\le 2\kappa_vD_v+2\sum_wA_{vw}D_w,
\qquad
B=2\left(\operatorname{diag}(\kappa_v)+A\right),
\qquad
\boldsymbol D(t)\le e^{Bt}\boldsymbol D(0).
\tag{7}
$$

最后一个不等式按分量成立：B 的非对角元非负，其线性系统保持非负锥，可用微分不等式比较。它保留每个位置的连接信息，比先对所有位置求和再估计更细。

纯 iid 制备给出每个位置 D_v(0)=2/N。又因 B 的第 v 行和为 4κ_v≤4κ，对 t≥0 有

$$
D_v(t)\le
\frac2N\left(e^{Bt}\boldsymbol1\right)_v,
\qquad
\sup_v D_v(t)\le\frac{2e^{4\kappa t}}N .
\tag{8}
$$

式 (5)—(7) 的恒等式与比较适用于其他初态的实际初始 D；但集中性结论需要 D(0) 足够小，不能把纯 iid 的 2/N 初始值套给任意关联态。孤立位置的 κ_v 为零，其 D_v 恰保持 2/N，而非随全网规模增长。

## 5. 集体均值还有单独的 O(1/N) 上界

令 b_v=⟨ŝ_v⟩−m_v，u_v=|b_v|。将精确量子期望展开，保留所有关联：

$$
\begin{aligned}
\dot{\boldsymbol b}_v={}&
2\left(\Omega_v\boldsymbol e_x+
\sum_wJ_{vw}m_{w,z}\boldsymbol e_z\right)
\times\boldsymbol b_v\\
&+2\sum_wJ_{vw}b_{w,z}\boldsymbol e_z\times\boldsymbol m_v
+2\sum_wJ_{vw}\boldsymbol e_z\times\boldsymbol C_{vw},\\
\boldsymbol C_{vw}&=\langle\delta_{w,z}\boldsymbol\delta_v\rangle,
\qquad |\boldsymbol C_{vw}|\le\sqrt{D_wD_v}.
\end{aligned}
\tag{9}
$$

相关向量界是逐分量 Cauchy–Schwarz 后求和的结果。局部转动不改变 b_v 长度。使用上右 Dini 导数处理 u_v=0 及最大位置切换，令 u_∞=max_v u_v，则

$$
\begin{aligned}
D^+u_v&\le2\sum_wA_{vw}u_w
+2\sum_wA_{vw}\sqrt{D_vD_w},\\
D^+u_\infty&\le2\kappa u_\infty+
\frac{4\kappa e^{4\kappa t}}N,\qquad u_\infty(0)=0,\\
u_\infty(t)&\le
\frac2N\left(e^{4\kappa t}-e^{2\kappa t}\right).
\end{aligned}
\tag{10}
$$

κ=0 时取零，和独立转动的精确结果一致。每位置保持内部置换对称，因此任意一个真实比特的约化态与 Bloch 向量 m_v 所代表纯态的迹距离为 u_v/2。

这是指定初态族的集体／单比特结论；没有证明完整联合态接近乘积态，没有处理任意未知参考输入，也不是 diamond 通道误差。

## 6. 关系传播：路径阶数控制，而不是硬边界

考虑同一 J、Ω 下两条经典轨道 m、m̃，只有初值不同。令 r_v=|m_v−m̃_v|。差分方程中的局部转动再次不改变范数，得到

$$
D^+r_v(t)\le2\sum_wA_{vw}r_w(t),
\qquad
\boldsymbol r(t)\le e^{2At}\boldsymbol r(0).
\tag{11}
$$

设初始改变只在位置集合 Y，且每点初始差异不超过 δ≤2。若 v 到 Y 的图距离为 r，则 A 的 n 次幂在 n<r 时没有连接这些位置的路径。因此

$$
\begin{aligned}
r_v(t)
&\le\delta\sum_{w\in Y}(e^{2At})_{vw}\\
&\le\delta\sum_{n=r}^{\infty}\frac{(2\kappa t)^n}{n!}\\
&\le\delta\exp\!\left(2\kappa t e^\mu-\mu r\right),
\qquad \mu>0 .
\end{aligned}
\tag{12}
$$

所有界还可与最大距离 2 取较小者。第二行使用 A^n 的行和≤κ^n，所以不引入 |Y| 或 M 的额外因子。若两个位置不在同一连通分量，相应矩阵指数元恰为零；有限连通距离下则一般有小尾巴，不能称作严格锥外零影响。图距离尚未被赋予米的意义，μ 给出的衰减参数也不是已经导出的光速。

对两套相应的纯乘积量子制备，将各自均值与其经典轨道比较，还可得到

$$
\left|
\langle\widehat{\boldsymbol s}_v(t)\rangle_1-
\langle\widehat{\boldsymbol s}_v(t)\rangle_2
\right|
\le
\delta\sum_{w\in Y}(e^{2At})_{vw}
+\frac4N\left(e^{4\kappa t}-e^{2\kappa t}\right).
\tag{13}
$$

该式是所选制备族的集体均值响应控制；最后一项是把两个经典逼近误差相加，远距离时可能很松。它不是对所有量子干预或算子对易子的普适传播定理，也不声称排除了更强的微观 Lieb–Robinson 控制。

## 7. 为什么 N 与 M 不能随意互换

逐位置界对任意有限 M 成立，因此也适用于预算一致的任意有限序列 M=M(N)。但这与全网同时无大偏差是不同命题。

在某个固定时刻，每位置选择**一个固定单位方向** q_v，联合测量 ŝ_v·q_v。不同位置的读出相互对易，所以这个联合概率有定义；没有要求同时测量同一位置的三个非对易分量。二阶矩界及并集上界给出

$$
\Pr\!\left[
\exists v:
\left|\widehat{\boldsymbol s}_v\cdot\boldsymbol q_v-
\boldsymbol m_v\cdot\boldsymbol q_v\right|>\epsilon
\right]
\le
\min\!\left(1,\frac{\sum_vD_v(t)}{\epsilon^2}\right)
\le
\min\!\left(1,\frac{2M e^{4\kappa t}}{N\epsilon^2}\right).
\tag{14}
$$

因此 M exp(4κT)/N→0 是这个估计所给的一个充分条件，不是最优或必要条件；这里是每个固定时刻的联合测量，不是整个时间区间上所有连续测量历史的概率。

更强的“不保证”见证在 t=0 已经成立。固定 0<ε<1，每位置 N 个比特全制备为 |+x⟩，读出每位置平均 Z。其经典预测为零，各位置的读出独立：

$$
\begin{aligned}
p_N(\epsilon)
&=2^{-N}
\sum_{\substack{0\le k\le N\\|2k/N-1|>\epsilon}}
\binom Nk,\\
P_{\mathrm{any}}(N,M)
&=1-(1-p_N)^M,\\
M_N&=\left\lceil1/p_N\right\rceil,\qquad
p_N\longrightarrow0,\quad
P_{\mathrm{any}}(N,M_N)\longrightarrow1-e^{-1}.
\end{aligned}
\tag{15}
$$

最后一个极限使用 M_Np_N→1。单个位置的偏差概率确实趋零，但足够多位置中至少一个失败仍有有限概率。这并不否定式 (8)—(10)，而是严格区分它们没有包含的全网量词。

## 8. 可复算结果

### 8.1 固定小图增加每位置 N

取三个位置的路径，每边 J=0.25，κ=0.5；数值位置按 v=0,1,2 编号，Ω_v=0.55 cos(0.6v)，初态在各位置不同，按代码的 initial_vectors 固定生成。t=0.8：

| N | 最大均值偏差 u_v | 最大均方误差 D_v | 式 (8) 一致上界 |
|---:|---:|---:|---:|
| 2 | 0.06302882 | 1.12246979 | 4.95303242 |
| 4 | 0.03200551 | 0.56212725 | 2.47651621 |
| 8 | 0.01614689 | 0.28132544 | 1.23825811 |
| 16 | 0.00811232 | 0.14073327 | 0.61912905 |

数值符合解析控制，但本轮的 O(1/N) 结论来自证明，不来自四点拟合。解析上界并不紧；小 N 时超过误差本身的简单绝对上界仍是合法但较松的控制。

N=16 时，位置比较矩阵给出的 D 上界为 (0.31834308, 0.51168615, 0.31834308)，比一致上界 0.61912905 更细；实际值为 (0.13216368, 0.14073327, 0.13363649)。

### 8.2 相同预算不选择同一图或同一响应

取 M=4、N=4、t=0.5，并把四种图归一到相同 κ=0.5：

| 图 | 最大实际均值偏差 |
|---|---:|
| 路径 | 0.01463818 |
| 环 | 0.01466970 |
| 星形 | 0.00924382 |
| 完全图 | 0.01015883 |

它们共同满足 D 的一致上界 1.35914091、均值偏差上界 0.53478028，但响应不相同。这表明预算一致与经典极限成立并未自动决定关系网络。代码另外检查 M=4、16、64 的图预算及同一解析常数，不对这些大 M 做指数成本的完整量子态模拟。

### 8.3 初始扰动沿路径传播

七位置路径、每边 J=0.25、Ω_v=0.6，t=0.8。只在第零位置将初始向量绕 y 轴转动 0.8，初始差异 δ=0.77529673：

| 到改变位置的距离 | 实际经典响应 | 加权矩阵指数界 | 仅使用 κ 的路径尾界 |
|---:|---:|---:|---:|
| 0 | 0.76514701 | 0.83899666 | 1.72545460 |
| 1 | 0.24893380 | 0.32699270 | 0.95015787 |
| 2 | 0.01712817 | 0.06454474 | 0.32992049 |
| 3 | 0.000339903 | 0.00853802 | 0.08182554 |

将经典 RK4 步长从 0.002 减半后，全部七位置响应向量的差异最大约 2.48×10⁻¹³。结果 JSON 保留全部距离数据；远端响应小不构成严格为零的证据。代码还检查矩阵幂在不足图距离的阶数上恰没有连接路径，并对相应小图量子均值验证式 (13)。

### 8.4 全网量词的精确反例

式 (15) 取 ε=0.5，严格使用“>ε”，并选择 M=ceil(1/p_N)：

| N | 单位置失败概率 p_N | M | 全网至少一处失败概率 |
|---:|---:|---:|---:|
| 16 | 0.02127075 | 48 | 0.64371035 |
| 32 | 0.00210240 | 476 | 0.63278090 |
| 64 | 0.00002436457 | 41044 | 0.63213221 |

这是二项系数的确定性计算，没有 Monte Carlo 抽样。固定 N=16 而将 M 从 1、10、100 增至 1000 时，该概率分别为约 0.02127、0.19346、0.88352、0.99999999954。

## 9. 数值方法与资源账

每位置内部置换对称性由 H_N 保持，初态也在该对称部门。因此可在每位置维数 N+1 的部门内做精确模型演化。实际资源与计算压缩分别为

$$
\begin{aligned}
\text{真实比特数}&=MN,\qquad
\dim\mathcal H_{\mathrm{full}}=2^{MN},\qquad
\dim\mathcal H_{\mathrm{sym}}=(N+1)^M,\\
\text{逐边微观配对数}&=|E|N^2,\qquad
\sum_{w:\{v,w\}\in E}\sum_{b=1}^{N}|J_{vw}|/N=\kappa_v,\\
\|H_N\|&\le
N\left(\sum_v|\Omega_v|+\sum_{\{v,w\}\in E}|J_{vw}|\right)
\le MN(\Omega_*+\kappa/2).
\end{aligned}
\tag{16}
$$

第二行对每条相邻位置边有 N 个伙伴求和，再对相邻位置求和。对称压缩降低了计算成本，没有消除真实比特、制备或交互资源；每比特相互作用预算虽可有界，总资源仍随 MN 增加。

数值演化使用 18 阶 Taylor 分步，每步的 |dt|‖H‖上界≤0.5。代码报告精确算术中的截断误差上界，并另做浮点守恒、逆向恢复、步长细化及完整小张量积对照；不把截断界当作浮点严格认证。经典参考使用 RK4，配合步长减半和能量、长度守恒检查。

三位置 N=16 的例子使用 4913 维对称部门，对应完整 48 比特系统。44 个 Taylor 子步的精确算术截断上界约 9.54×10⁻²²；实际范数误差约 1.33×10⁻¹⁵，单位 N 的能量误差约 2.22×10⁻¹⁶。因此这里浮点误差主导于所列截断上界。

**15 项科学检查通过**：完整张量积与对称部门对照；各位置初始涨落；随机关联态和带符号权重下逐位置涨落恒等式；局部比较及均值界；孤立位置；经典能量、长度及细化；经典路径响应；路径阶数；量子响应界；四种图的同预算对照；负 J 的绝对权重控制；量子守恒、逆向及细化；精确二项反例；不同 M 的预算一致性；Hamiltonian 与资源上界。检查先通过，再写入结果。

## 10. 范围结论与后续接口

本轮给出的正接口是：**指定关系网络和相互作用以后，可以同时拥有不随网络位置数失控的逐位置经典集中，以及受关系路径控制的响应。** 不必先把所有位置合并成维数随 M 指数增长的大单元，才能证明误差控制。

仍未解决的生成问题包括：

- 关系图和权重为何出现、如何自主演变；本轮没有生成这些数据。
- 哪个可观测量及何种测量决定物理距离、体积或维数。
- 如何把有限图族接为相容的无限体积量子态、代数及连续极限；一致误差常数本身不能代替这项构造。
- 若把关系权重变成量子动力学变量，能否同时保留图局部控制、资源界及闭合的约束结构。

因此下一步适合研究**关系本身具有动力学时的控制条件**，或对接有明确图族与尺度的连续场极限；不能仅凭本轮图上的经典动力学，把尚未生成的物理空间或 Einstein 方程写入结论。

## 11. 复算与文件

- [代码：relational_collective_limit_audit.py](relational_collective_limit_audit.py)
- [结果：relational_collective_limit_audit_results.json](relational_collective_limit_audit_results.json)
- [本轮整合检查](research_round_361_checks.json)
- [科学基线：第357轮](research_note_357.md)
- [复用的集体表示：第354轮代码](collective_poisson_limit_audit.py)

在本目录使用已有 Python/NumPy 运行时：

    & 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 relational_collective_limit_audit.py

只有需要创建结果且科学检查已通过时，追加 --write-results；运行器拒绝以不同内容覆盖已有结果。当前已保存 Python 3.12.14、NumPy 2.3.5 下的 15 项通过结果。
