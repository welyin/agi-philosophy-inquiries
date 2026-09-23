# 第353轮：联合量子几何的均值极限、时间窗口与反作用尺度

## 0. 本轮问题与结论

本轮承接[第350轮](research_note_350.md)已冻结的联合量子反作用模型。先复读当前 README、research_direction、RESEARCH_STATE，并只读复核第352轮新增笔记；后者的多标量约束结论与 P(X) 声锥反例相容。本轮不依赖第351、352轮的未冻结新推导。

问题是：**当“几何寄存器”越来越大，用其均值替代算符，能否成为保持任意参考关联的受控近似？与此同时，物质对几何的作用是否仍然存在于同一联合演化中？**

得到四个明确结论：

1. 独立乘积几何、固定总耦合强度和固定时间下，均值酉通道的误差有严格的 1/N 上界，包含任意物质—参考关联。
2. 当时间增长至 √N 尺度，误差趋向非零值；同一均值但长程相关的几何态甚至在固定时间也不收敛。均值本身不足以决定是否能够经典化。
3. 同一封闭联合演化会改变几何的不对易观测量，并产生可计算的几何—物质纠缠；约化噪声不能替代完整联合账本。
4. **单个物质探针、固定总耦合范数的本模型，固定时间的平均几何反作用也随 1/N 消失。** 因而“几何趋于确定”和“存在非零宏观反作用”是两个不同待证条件。本轮没有推出经典引力反作用，更没有推出 Einstein 方程。

## 1. 来源、认知动机和额外输入

### 1.1 已有原始成果的接口

本轮模型属于已有 central-spin 纯退相干模型。Cucchietti、Paz、Zurek 的 [Decoherence from Spin Environments，Phys. Rev. A 72, 052113 (2005)](https://arxiv.org/html/quant-ph/0508184) 第 II 节式 (1)、(8)、(12)–(16) 给出联合模型及乘积相干因子，式 (22) 讨论高斯退相干。这里采用等耦合且随 N 归一化的特殊族，直接推导自己的误差界及适用时间，不把这些已知机制宣称为原创，也不借用一般环境态都满足高斯极限的断言。

本轮新增的项目内成果是：把已有机制接入第350轮的参考保持标准，明确固定时间、长时间、相关态和宏观反作用四个量词，给出同一联合量子规则的可复算账目。

### 1.2 认知动机

“整体内部交代清楚”的要求促使我们同时保留作为几何接口的 G 和物质 M，检查约化近似隐藏了哪些关联。认知动机没有直接给出下面的 Hamiltonian、独立初态或耦合缩放。

### 1.3 明示的额外建模输入

- G 为 N 个量子比特，M 为一个量子比特；R 为任意未受操作的参考。
- 选择集体 Z 平均量与 ZZ 型耦合；暂不加入 G、M 的自 Hamiltonian。
- 初态为固定的几何态与任意 MR 态的直积。允许 G 内部相关，但不允许在同一个约化通道声明中任意改变初始 G—MR 关联。
- 主要收敛命题进一步要求几何为同均值的独立乘积态。相关态另作反例，不能混用。
- 随 N 扩大时，每对耦合强度按 1/N 缩小；比较固定 g、固定 t，因此固定 θ=gt。

“几何”在此是接口名称；尚未定义空间位置、度规、曲率或 Einstein 约束。N 也不是空间维数。

## 2. 同一联合规则与精确约化通道

取 ℏ=1，并定义

$$
\mathcal H_G=(\mathbb C^2)^{\otimes N},\qquad
Q_N=\frac1N\sum_{j=1}^{N}Z_j,\qquad
H_N=gQ_N\otimes Z_M,\qquad
U_N(\theta)=e^{-i\theta Q_N\otimes Z_M},\qquad
\theta=gt,\quad \|H_N\|=|g|.
\tag{1}
$$

设几何单比特态为 γ，满足 tr(γZ)=μ。它可以是混态，也可以是具有非零 X 相干的纯态；两者在下面的物质约化通道中只通过 Z 概率分布出现。初态为 γ 的 N 重张量积与 ρ_MR 的直积。

把 ρ_MR 按 M 的 Z 基底写成参考空间上的正块矩阵。精确通道为

$$
\rho_{MR}=
\begin{pmatrix}A&C\\C^\dagger&B\end{pmatrix}
\ \longmapsto\
(\Phi_{N,\mu,\theta}\otimes\mathrm{id}_R)(\rho_{MR})
=
\begin{pmatrix}A&\Gamma_N C\\\Gamma_N^*C^\dagger&B\end{pmatrix},
\qquad
\Gamma_N=
\left[\cos\frac{2\theta}{N}-i\mu\sin\frac{2\theta}{N}\right]^N .
\tag{2}
$$

**证明。** 将 Q_N 谱分解代入式 (1)，迹掉 G 后，每个本征值 q 贡献相位 exp(−2iθq)。独立乘积态使这一特征函数分解为 N 个单比特因子。G 的非对角矩阵元在约化通道中不出现，但仍保留在完整联合态中，不能因此把联合模型宣布为经典模型。

均值近似是固定酉 V=exp(−iθμZ_M)。这里 μ 是预先指定几何初态的参数，不是随未知物质态重新计算的函数。因此精确通道、均值通道及它们的极限均保持混合仿射性，没有重复第350轮所排除的非仿射分支演化规则。

## 3. 任意参考保持的均值误差

定义 ε 为通道差的 half-diamond 范数，即允许任意参考时输出态迹距离的最大值。对更一般的固定几何态 σ_G，记其 Q_N 均值仍为 μ，相干因子为 Γ。则

$$
\begin{aligned}
\varepsilon
&=\frac12\|\Phi-\mathcal V_\mu\|_\diamond
=\frac12\left|\Gamma-e^{-2i\theta\mu}\right|\\
&\leq \min\!\left\{1,\theta^2\operatorname{Var}_{\sigma_G}(Q_N)\right\},\\
\sigma_G=\gamma^{\otimes N}
&\quad\Longrightarrow\quad
\varepsilon_N\leq
\min\!\left\{1,\frac{\theta^2(1-\mu^2)}{N}\right\}.
\end{aligned}
\tag{3}
$$

**精确范数的证明。** 两通道输出之差只有非对角块 δC，其中 δ=Γ−exp(−2iθμ)。其迹距离为 |δ|‖C‖₁。正块矩阵满足 ‖C‖₁≤√(tr A tr B)≤1/2；M 的 |+〉态达到等号，Bell 参考态也达到等号。所以该界已包含任意维 R，而不是只对一个选定输入的拟合。

**方差界的证明。** 对 Q_N 的谱随机变量 q，令 ξ=q−μ，则 Eξ=0。标量 Taylor 积分余项给出 |exp(−ix)−1+ix|≤x²/2。取 x=2θξ 并平均，得到 |E exp(−2iθξ)−1|≤2θ² Var(Q_N)，再除以 2。此论证适用于任意固定几何态；独立性只用于最后的方差缩小。

固定 θ、固定 |μ|<1 时还可得到主项：

$$
\log\!\left(\Gamma_N e^{2i\theta\mu}\right)
=-\frac{2\theta^2(1-\mu^2)}{N}
+O_{\theta,\mu}(N^{-2}),
\qquad
\varepsilon_N=
\frac{\theta^2(1-\mu^2)}{N}+O_{\theta,\mu}(N^{-2}).
\tag{4}
$$

当 μ=±1 时方差为零，均值近似对所有 N、θ 都精确。式 (4) 的余项量词是固定 θ，不能直接拿来作全时间一致结论。

## 4. 长时间窗口：√N 尺度留下有限退相干

取 θ_N=c√N，固定 c≠0 及 |μ|<1。对单因子在 2θ_N/N=2c/√N 附近展开，N 倍三阶余项为 O(N^−1/2)，得到

$$
\Gamma_N(\theta_N)e^{2i\theta_N\mu}
\ \longrightarrow\
e^{-2c^2(1-\mu^2)},\qquad
\varepsilon_N(\theta_N)
\ \longrightarrow\
\frac{1-e^{-2c^2(1-\mu^2)}}2>0.
\tag{5}
$$

这不是对固定时间结论的否定，而是其适用窗口：在固定 g 下，t 也随 √N 增长。特别是 θ=o(√N) 可由式 (3) 保证误差趋零；达到 √N 尺度一般已无法用确定均值酉代替。

对 μ=0.3、c=0.6，极限误差为 0.2403309077。N=100000 的闭式计算给出 0.2403312726。持续增加几何自由度并不自动保证任意长时间的经典近似。

## 5. 同一均值不是同一宏观状态：相关几何反例

因为不同 Z_j 相互对易，任意几何态都满足

$$
\operatorname{Var}(Q_N)
=\frac1{N^2}\sum_{i,j=1}^{N}
\left(\langle Z_iZ_j\rangle-\langle Z_i\rangle\langle Z_j\rangle\right).
\tag{6}
$$

因此独立性是一个充分条件；更一般地，若所有协方差绝对值的总和为 O(N)，也足以得到 O(1/N) 方差。这是额外的相关性条件，不能由“大 N”本身推出。

考虑一个与乘积族具有相同局部 Z 均值 μ 的宏观 cat 态：

$$
\begin{aligned}
|\mathrm{cat}_{N,\mu}\rangle
&=\sqrt{\frac{1+\mu}{2}}|0\rangle^{\otimes N}
+\sqrt{\frac{1-\mu}{2}}|1\rangle^{\otimes N},\\
\langle Z_j\rangle&=\mu,\qquad
\operatorname{Var}(Q_N)=1-\mu^2,\qquad
\Gamma_{\mathrm{cat}}(\theta)=\cos(2\theta)-i\mu\sin(2\theta).
\end{aligned}
\tag{7}
$$

其 Q_N 只取 ±1，方差与相干因子均不随 N 缩小。μ=0、θ=π/4 时，物质 |+〉相对均值预测的迹距离恒为 1/2。这里仅要求同一个局部 Z 均值，并没有声称 cat 态与具有非零局部 X 的纯乘积态具有相同完整局部密度矩阵。

把 cat 的两分支相干删去，会得到同一物质约化通道；但联合态不同。初始物质为 |+〉且采用上述 μ、θ 时，纯 cat 的几何—物质 negativity 为 1/2，经典分支混合的 negativity 为零。前者是纠缠，后者是可分的经典记录，不能仅凭物质退相干曲线判断。

## 6. 几何的共轭响应与宏观反作用的缩放缺口

选择几何单比特纯态，Bloch 向量为 (x,0,μ)，其中 x=√(1−μ²)。设初始物质的 Z 均值为 m。式 (1) 的各对耦合对易，直接计算给出

$$
\begin{aligned}
U_N^\dagger Y_jU_N
&=Y_j\cos\frac{2\theta}{N}
+X_jZ_M\sin\frac{2\theta}{N},\\
\langle Y_j\rangle_{\mathrm{out}}
&=xm\sin\frac{2\theta}{N},\\
\operatorname{Cov}_{\mathrm{out}}(Y_j,Z_M)
&=x(1-m^2)\sin\frac{2\theta}{N}.
\end{aligned}
\tag{8}
$$

Z_j 保持不变不等于没有反作用：其不对易的 Y_j 发生变化。即使 m=0 使平均 Y_j 为零，联合 Y_j—Z_M 关联仍可以非零。cat 态的相应响应发生在逻辑 Y_L=−i|0…0〉〈1…1|+i|1…1〉〈0…0|上，其期望为 xm sin(2θ)，与 N 无关；该观测量是全局逻辑观测量，不是平均局部 Y。

然而规范化的集体几何观测量满足

$$
\overline Y_N=\frac1N\sum_{j=1}^{N}Y_j,\qquad
\|\overline Y_N\|=1,\qquad
\langle\overline Y_N\rangle_{\mathrm{out}}
=xm\sin\frac{2\theta}{N},\qquad
\left|\Delta\langle\overline Y_N\rangle\right|
\leq \frac{2|\theta xm|}{N}.
\tag{9}
$$

本模型的初始平均 Y 为零。由此应区分三个时间尺度：

| 无量纲时间 | 乘积几何的物质均值误差 | 平均几何 Y 响应 |
|---|---|---|
| 固定 θ | O(1/N) | O(1/N)，趋零 |
| θ=c√N | 一般趋于非零值 | O(1/√N)，仍趋零 |
| θ=αN | 对固定 abs(μ)<1、sin(2α)≠0，误差趋于 1/2 | xm sin(2α)，可保持非零 |

第三行排除了相干因子模长等于 1 的特殊复现点；它不保证每个 α 都有非零响应。g 固定时，等待时间分别是常数、√N 和 N 尺度，并没有免费加速。

总和 ΣY_j 的固定时间响应可趋向 2θxm，但其算符范数为 N。这不等于范数为 1 的宏观平均量具有有限变化，也不能直接解释为经典引力反馈。**增加几何规模同时固定单个探针和总耦合范数，得到的是受控的近似背景，而不是非零经典反作用的闭合理论。**

## 7. 完整联合账本、参考与资源

在 μ=0、初始 G 为 |+〉的 N 重乘积且 M 为 |+〉时，完整联合态始终纯。令 D 为迹距离，Neg 为几何—物质 negativity，则

$$
\begin{aligned}
D\!\left(U_N|\!+\rangle\langle+\!|^{\otimes(N+1)}U_N^\dagger,\,
|\!+\rangle\langle+\!|^{\otimes(N+1)}\right)
&=\sqrt{1-\cos^{2N}(\theta/N)},\\
\mathrm{Neg}_{G:M}
&=\frac12\sqrt{1-|\Gamma_N|^2},\\
\operatorname{Cov}\!\left(\sum_jY_j,Z_M\right)
&=N\sin(2\theta/N).
\end{aligned}
\tag{10}
$$

第一式来自初末纯态内积 cos^N(θ/N)。第二式来自纯态 Schmidt 系数，即物质密度矩阵的本征值 (1±|Γ_N|)/2。它们是联合量子账，不是只给物质端添加一个经典噪声。

固定时间时，物质通道误差为 O(1/N)，但联合态相对不变几何加均值酉的误差一般仅为 O(1/√N)。更一般地，对固定任意 σ_G 与任意 ρ_MR 的初始直积，纯化后应用 |exp(−ix)−1|≤|x|，再利用偏迹的迹距离压缩性，可得

$$
D\!\left(
U_N(\sigma_G\otimes\rho_{MR})U_N^\dagger,\,
\sigma_G\otimes(V_\mu\otimes I_R)\rho_{MR}(V_\mu^\dagger\otimes I_R)
\right)
\leq
\min\!\left\{1,|\theta|\sqrt{\operatorname{Var}_{\sigma_G}(Q_N)}\right\}.
\tag{11}
$$

纯化证明中，两演化之差由 exp[−iθ(Q_N−μ)Z_M]−I 给出，其平方范数上界为 θ²Var(Q_N)，因为 Z_M²=I。这再次限定了初始 G 与 MR 独立；不声称任意初始 G—MR 关联都能由固定的物质通道表示。

有限系统没有不可逆地丢失关联。该模型还具有严格复现：

$$
\dim\mathcal H_{GM}=2^{N+1},\qquad
U_N(\pi N)=(-1)^NI_{GM},\qquad
\langle H_N\rangle_{\mathrm{out}}=\langle H_N\rangle_{\mathrm{in}}.
\tag{12}
$$

总联合演化守恒能量；几何存储成本为 N 个量子比特。式 (1) 的归一化是输入：若改成每对固定耦合 g，则 Hamiltonian 范数随 N 增长，不能沿用同一固定资源的结论。初态制备、重新使用时的重置、观测记录和读出没有由本模型自动供给。约化噪声下的后续使用也不能未经检查就当作独立同分布的环境调用。

## 8. 可复算结果与检查

代码：[collective_geometry_limit_audit.py](collective_geometry_limit_audit.py)。

结果：[collective_geometry_limit_audit_results.json](collective_geometry_limit_audit_results.json)。

汇总核验：[research_round_353_checks.json](research_round_353_checks.json)（根代理独立审核后生成）。

固定 θ=0.7、μ=0.3：

| N | 精确 half-diamond 误差 | 方差上界 |
|---:|---:|---:|
| 10 | 0.0427739551 | 0.04459 |
| 100 | 0.00443929909 | 0.004459 |
| 1000 | 0.000445701356 | 0.0004459 |
| 10000 | 0.0000445880122 | 0.00004459 |

N=4、μ=0、θ=0.7 的完整联合矩阵核验：物质误差 0.1106665067，联合态误差 0.3403698736，negativity 为 0.3137187132。两种误差并不等同。

17 项 unittest 覆盖：

- Q_N 谱范数、乘积均值与方差；小 N 全矩阵对闭式通道。
- 纯及混几何、Bell 参考饱和误差、任意参考随机态与全参数方差界。
- 固定时间主项、√N 非一致极限及零方差端点。
- Heisenberg 共轭响应、平均 Y 与联合协方差、纯联合态距离与纠缠。
- 同均值 cat、经典混合的联合区别、cat 逻辑观测响应和协方差总和。
- 三个时间尺度、能量守恒、有限复现及固定几何下的仿射性。

运行环境：既有 Python 3.12.14、NumPy 2.3.5，无新依赖、无图片。

在仓库根目录运行：

    & 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 research_cognition_physics/archive_231_/collective_geometry_limit_audit.py
    & 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 research_cognition_physics/archive_231_/collective_geometry_limit_audit.py --write-results

## 9. 本轮关闭的缺口和下一步

**解析已证明：** 在指定的联合量子模型和乘积/弱相关几何条件下，固定时间的确定均值酉是一个保持任意参考的受控近似；相关态和长时间给出精确边界。完整联合反作用存在，但单探针的平均几何固定时间响应消失。

**数值已验证：** 小 N 的完整联合矩阵与解析式一致；大 N 闭式验证收敛速度及两种非收敛情形。数值没有承担“对所有 N”的证明，所有全称结论均来自前述推导。

**物理解释的限度：** 这不从 FUCP 选择唯一 Hamiltonian，不建立空间、普适物质度规、时空约束或 Einstein 方程。它也没有证伪全部半经典理论；它给出某一确定均值近似成立和失败的可检验条件。

下一独立问题应当是：**能否让物质自由度也随规模增长，在清楚交代相互作用能量、耦合归一化和时间预算的前提下，同时保有小量子涨落与非零的归一化几何反作用？** 如果能，需要证明联合量子演化到相应双向平均场方程的误差；如果不能，应指出冲突来自哪一组具体资源缩放。该问题仍只处理“量子联合动力学怎样产生经典反馈”这一步，尚未选择引力的具体方程。
