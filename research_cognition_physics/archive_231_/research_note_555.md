# 第555轮：同一相互作用能量矩统一控制参考读取——非Gaussian来源与全时间范围

日期：2026-09-30。接续[554](research_note_554.md)、[548](research_note_548.md)与[共同条件合并顺序](joint_condition_consolidation_order.md)。[代码](joint_matter_energy_moment_control.py)、[保存结果](joint_matter_energy_moment_control_results.json)、[核验](research_round_555_checks.json)、[条件账](unified_physics_condition_ledger_555.md)。

## 1. 本轮合并了什么

554证明：保持同一物质H的平均能量和全部二阶协方差，仍不能统一控制导数平方读口的未裁剪风险。本轮给出相应的充分条件：**在同一有限胞元、正定四次势模型中，一个有界的二阶能量矩预算，可以同时控制四个参考读数的实际二阶矩；来源无需Gaussian，并且这个预算在未受测量干预的同一静态H演化下保持。**

具体增量是把来源的高阶矩、相互作用演化和554实际共同POVM接成同一个估计。不再为每个读口分别添加一个高阶尾部条件；但保留一个明确的来源限制及原有仪器输入。这是条件性压缩，不是已经从认知原则推出来源预算，也不声称它必要、最优或足以给清晰坐标。

|层次|地位|
|---|---|
|认知动机|参考的统计可靠性与内部能量账应由共同模型连接|
|继承|544同一双标量势、554全部有限模式、同一导数菜单和独立Gaussian种子POVM|
|来源新条件|同一完整相互作用H的二阶矩存在；统一风险族另给其共同上界|
|解析结论|图范数界、全部菜单二阶矩界、任意非Gaussian来源的实际记录界与未测演化下保持|
|数值见证|一／二／八胞元；含Fock非Gaussian源、非零梯度、独立正测量密度积分及尾反例对照|
|未包含|连续极限、维数选择、Higgs径向约化有效性、内生准备与探针、连续监测反馈和量子引力|

非负Schrödinger算符的自伴性是成熟工具，可对接[Simon 1973](https://link.springer.com/article/10.1007/BF01427943)；下面为本光滑有限模型给出所需截断证明，复用[528定义域处理](research_note_528.md)，不把它算作新物理。相空间协变POVM直接继承[Werner §3](https://arxiv.org/pdf/quant-ph/0405184)及554；这里新增的是该仪器与相互作用能量的具体连接。

## 2. 同一H与可核的势常数

N个胞元，每胞体积v>0，总V=Nv。q=(q_h,q_s)共有d=2N个实坐标，p=−iℏ∇。给定实差分D_r，K=Σ_rD_rᵀD_r≥0；两场使用同一K。设L=[[λ_h,p],[p,λ_s]]正定，u=L⁻¹C的两分量为正。门户p与动量p_a以指标区分。取与554相同的真空常数Cᵀu/4：

$$
H=T+W,\qquad T=-a\Delta,\quad a=\frac{\hbar^2}{2v},\qquad
W=\frac v2\sum_{a=h,s}q_a^{\mathsf T}Kq_a+
\frac v4\sum_{j=1}^N(x_j-u)^{\mathsf T}L(x_j-u)\ge0,
\quad x_j=(q_{h,j}^2,q_{s,j}^2).
\tag{1}
$$

以下固定有限N、v、K及势参数，不把H的能量零点任意平移后继续沿用常数。令ℓ为L最小本征值，α=max(0,3λ_h+p,3λ_s+p)。对任意r>0，逐坐标使用q_i²≤u_i+(q_i²−u_i)²/(2r)+r/2，可得全空间逐点界

$$
\Delta W\le A_r W+B_r,\qquad
A_r=\frac{2\alpha}{\ell r},\qquad
B_r=2v\operatorname{Tr}K+
v\alpha\!\left[N(u_h+u_s)+Nr\right]+vN(|C_h|+|C_s|).
\tag{2}
$$

证明：势Laplacian恰为2vTrK+vΣ_j[(3λ_h+p)q_hj²+(3λ_s+p)q_sj²−C_h−C_s]；正四次项给Σ_i(q_i²−u_i)²≤4W/(vℓ)。代入即得式(2)。常数保守，r是证明参数，不是新物理阈值。

## 3. 图范数估计与定义域

对紧支撑光滑ψ先分部积分，随后可对Schwartz波函数使用同式：

$$
\begin{aligned}
\|H\psi\|^2
&=\|T\psi\|^2+\|W\psi\|^2
+2a\int W|\nabla\psi|^2-a\langle\psi,\Delta W\psi\rangle,\\
\|T\psi\|^2+\|W\psi\|^2
&\le\|H\psi\|^2+aA_r\langle\psi,H\psi\rangle+aB_r\|\psi\|^2.
\end{aligned}
\tag{3}
$$

这里仅用了非负二次型〈W〉≤〈H〉，没有使用通常不成立的非对易平方序H²≥T²。

**自伴与延拓。** W为光滑非负多项式。若ψ∈ker(H*+1)，局部椭圆正则允许以χ_R²ψ测试，其中χ_R为半径R内等于1、梯度O(1/R)的光滑截断。得到a∫|∇(χ_Rψ)|²+∫(W+1)|χ_Rψ|²=a∫|∇χ_R|²|ψ|²→0，故ψ=0。非负对称算符的该核判据给出C_c^∞上的本质自伴，取唯一闭包H。

对任意ψ∈D(H)，取C_c^∞的H图范数逼近ψ_n。把式(3)用于差向量，Tψ_n与Wψ_n分别Cauchy；T与乘法W的闭性给ψ∈D(T)∩D(W)及同界。反向包含由分布意义下Tψ+Wψ∈L²及H*=H得到。因此D(H)=D(T)∩D(W)。这也证明所需二阶导数及四次位置矩不会在闭包步骤被遗漏。

对正常混态ρ，以谱分解逐项求和。记E₁=TrρH、E₂=TrρH²；后者是正算符二次型期望，不要求每个向量属于D(H²)。当E₂有限，定义

$$
\mathcal G=E_2+aA_rE_1+aB_r,\qquad
\operatorname{Tr}\rho(T^2+W^2)\le\mathcal G,\qquad
\left\langle\sum_iq_i^2\right\rangle
\le S:=N(u_h+u_s)+Nr+\frac{2E_1}{v\ell r},\qquad
\left\langle\sum_ip_i^2\right\rangle\le2vE_1.
\tag{4}
$$

若只指定一个预算E₂≤𝓔²，则E₁≤𝓔；将这两个上界代入即可得到同一来源族的统一常数。因此不需再把E₁当作独立预算。

## 4. 全部原菜单的源二阶矩

保留554的Q_a=N⁻¹Σ_jq_aj及A_a=N⁻¹(q_aᵀKq_a−p_aᵀp_a/v²)。令U_a=v q_aᵀKq_a/2，T_a=p_aᵀp_a/(2v)，则A_a=2(U_a−T_a)/V。U_a与W都是乘法且0≤U_a≤W；T_h、T_s强对易且非负。因此

$$
\langle Q_a^2\rangle\le\frac SN,\qquad
\|A_a\psi\|^2\le\frac8{V^2}
\left(\|W\psi\|^2+\|T\psi\|^2\right),\qquad
\langle A_a^2\rangle\le\frac8{V^2}\mathcal G.
\tag{5}
$$

A_a²期望指‖A_aψ‖²及其混态和，不要求ψ∈D(A_a²)。这是同一相互作用H下的界，未用自由模能量替代完整能量，也没有把平均导数平方换成均匀模平方。

## 5. 任意来源的真实POVM二阶矩

沿用554独立、零均值Gaussian种子，其相空间协方差记为T_τ，避免与动能T混淆。源态ρ完全任意。令z=(q,p)，Ω=[[0,I],[−I,0]]。对实对称矩阵M，F=(zᵀMz)^W为Weyl量子化，实际后处理Y=z_outᵀMz_out−Tr(MT_τ)。协变POVM卷积矩恒等式与二次Weyl乘积给

$$
\mathbb E_\rho Y=\langle F\rangle_\rho,\qquad
\mathbb E_\rho Y^2
=\langle F^2\rangle_\rho+
4\langle(z^{\mathsf T}MT_\tau Mz)^W\rangle_\rho+
2\operatorname{Tr}(MT_\tau MT_\tau)
-\frac{\hbar^2}{2}\operatorname{Tr}(M\Omega M\Omega).
\tag{6}
$$

推导可先对Schwartz源用特征函数求导：加Gaussian种子后的二次多项式平方比原Weyl符号平方多4zᵀMT_τMz+2Tr(MT_τMT_τ)；而(F²)的Weyl符号为(zᵀMz)²+ℏ²Tr(MΩMΩ)/2。相减得式(6)。这不假设源Wigner函数为正，不把非Gaussian源当作经典正概率。

式(6)先在Schwartz核上陈述。对下述实际菜单，式(6)在D(H)上的二次型意义下延拓：对应测量等距V_M:ψ↦联合输出振幅，YV_M在Schwartz核上的范数由下面的H图范数界控制；闭算符Y与等距V_M使其沿H图范数逼近收敛。谱分解再给混态。这里不依赖未证的四阶算符强域，也未为任意未受控菜单宣称相同延拓。

取各模式同一T_τ=diag(t_qI,t_pI)，t_qt_p≥ℏ²/4。代入本菜单M_q=K/N、M_p=−I/(Nv²)，得到精确式

$$
\mathbb E Y_{A_a}^2
=\langle A_a^2\rangle+
\frac{4t_q}{N^2}\langle q_a^{\mathsf T}K^2q_a\rangle+
\frac{4t_p}{N^2v^4}\langle p_a^{\mathsf T}p_a\rangle+c_\tau,\qquad
c_\tau=\frac2{N^2}\!\left[t_q^2\operatorname{Tr}K^2+
\frac{t_p^2N}{v^4}\right]
-\frac{\hbar^2\operatorname{Tr}K}{N^2v^2}.
\tag{7}
$$

最后一项是**负号**：本菜单Tr(MΩMΩ)=+2TrK/(N²v²)。独立正测量密度积分检出了初稿具体代入时的符号错误，已修正；不是把原量子方差直接当作记录方差。每个K本征值κ≥0有2(t_qκ−t_p/v²)²+(4t_qt_p−ℏ²)κ/v²≥0，故c_τ≥0。

K²≤‖K‖K及q_aᵀKq_a≤2W/v给最终共同上界：

$$
\begin{aligned}
\mathbb E Y_{Q_a}^2&\le B_Q:=\frac SN+\frac{t_q}{N},\\
\mathbb E Y_{A_a}^2&\le B_A:=
\frac8{V^2}\mathcal G+
\frac{8E_1}{N^2}\left(\frac{t_q\|K\|}{v}+
\frac{t_p}{v^3}\right)+c_\tau .
\end{aligned}
\tag{8}
$$

两个场共同使用B_Q、B_A，既包括全部源涨落，也包括既定仪器噪声。只减去了已知种子偏置，没有免费扣除未知来源的偏置。

## 6. 时间保持的准确量词

对未进行中途测量、未外加控制的静态H演化，谱定理给

$$
\rho(t)=e^{-itH/\hbar}\rho(0)e^{itH/\hbar},\qquad
\operatorname{Tr}\rho(t)H^k=\operatorname{Tr}\rho(0)H^k
\quad(k=1,2).
\tag{9}
$$

所以式(8)在任意选择的末读时刻成立，无需假设状态保持Gaussian。它是“任选时刻的一次读取”的统一界，不是对全部连续时刻同时成功的概率，也不是连续监测、重复制备或测量后反馈的能量守恒。参考均值可以随时间移动；这个矩界不证明原雅可比标定、h正径向片或真实坐标身份一直保持。

对任意预先声明的常数基准c_i与非零固定标定j_i，设e_i=(Y_i−c_i)/j_i。由L²三角不等式与Markov并集界，

$$
\mathbb E e_i^2\le R_i:=
\frac{(\sqrt{B_i}+|c_i|)^2}{j_i^2},\qquad
\Pr\{\max_i|e_i|>r_0\}\le\frac{\sum_iR_i}{r_0^2}.
\tag{10}
$$

取r₀=√(ΣR_i/δ)，可沿554有限分箱形成四个有限符号记录；预算大时箱宽也大，不宣称精确定位。误差基准必须声明，不能把无噪声经典曲线当作自动保持的均值。这里证明有限风险而非有用分辨率；物理定位仍需准备、标定保持及尺度连接。

## 7. 反例对照与条件边界

复用554的正交动量尾混合：以权c/M²施加±M位移。势与位置资料不变，平均完整H增加c/(2v)，但二阶能量和实际平方记录会增长。两胞元的同一冻结势、v=1及c=1给

$$
E_1=2.6435576961,\qquad
\begin{array}{c|rr}
M&E_2&\mathbb E Y_{A_h}^2\\ \hline
2&8.76492024&5.50924556\\
4&11.76492024&8.50924556\\
16&71.76492024&68.50924556\\
64&1031.76492024&1028.50924556
\end{array}
\tag{11}
$$

这解释新增预算排除了什么：固定平均能量族并不落在同一个有界E₂来源族内。不是重做554“不足”的命题，而是核对本轮充分条件恰好控制该障碍。有限E₂不是所有仪器的必要条件，也不是所有可能来源条件中最弱的一个。

**正四次约束不能无说明地删掉。** 取N=1、两自由场、v=ℏ=1及W=0。频率参数ω的两份Gaussian包只平移h位置M，则同一自由H下

$$
E_1=\frac{\omega}{2},\qquad E_2=\frac{\omega^2}{2},
\qquad \mathbb E Y_{Q_h}^2=M^2+\frac1\omega\longrightarrow\infty.
\tag{12}
$$

所以即使控制全部动能矩，也不能单靠它控制绝对场值。这里的充分条件依赖模型已有的约束势，而不是普适地宣布“能量控制所有读数”。ω=1.3、M=64的记录二阶矩约4096.769。仍不排除只读差值、改用相对参考等其它任务。

## 8. 可复算验证与来源存在

代码以多项式乘Gaussian的波函数直接施加完整H及菜单，计算真正算符范数；并用独立的正Husimi输出密度积分核实际二阶矩。Fock数态包含n=1、2，故不是只测试Gaussian源。非零梯度两胞元例保留Moyal修正，八胞元例复用554几何及初值型。

|胞元／来源|E₁|E₂|实际h导数读口二阶矩|共同上界B_A|
|---|---:|---:|---:|---:|
|1，数态(1,2)|3.064431|10.286951|12.252500|289.657095|
|2，数态(1,0,2,1)，v=0.7|7.596794|60.873339|19.105945|641.620625|
|2，Gaussian|2.143558|4.971363|2.943861|56.005971|
|8，位移Gaussian|32.149428|1072.561176|5.812500|187.344549|

这些上界较宽，作用是统一控制而非精度优化。六组检查覆盖势常数、图范数／全菜单、独立实际POVM积分、共同预算与有限置信、同平均能量尾族、无约束势反例。图范数恒等式最大相对残差8.06×10⁻¹⁵，独立记录积分最大差1.60×10⁻¹²。全族及全时间结论由解析证明承担，没有把有限维截断传播模拟冒充无限维演化验证。

这些显式Schwartz态证明来源类非空；它们并非相互作用真空、热平衡态或已实现的制备协议。h仍在整个实轴上量子化，非Gaussian例不能直接解释成严格Higgs径向物理态。仪器种子、作用时机与记录装置的来源仍为独立输入。

## 9. 条件账与下一项

本轮将C19来源矩、C22完整能量和C03／C21实际记录接到同一有限物质H，并给随未测演化保持的共同充分条件。它减少逐读口指定高阶来源资料的需要，未消除制备预算本身。

此接口达到有限模型的充分性结论后，不继续只优化常数。下一项返回统一条件账，优先检查实际仪器、区域组合与共同物质的连接：既定POVM是否能由声明的探针物质和交互实现，以及规范约束下哪些读取依赖外部参考。复用旧有限维仪器和区域拼接结果；只有同一候选的新增约束或共同实现才增加轮次。维数、连续几何、完整引力及统一目标仍未完成。
