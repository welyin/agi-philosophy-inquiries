# 第三十五轮：把数值耦合转成整数证书，完整可行下界达到 0.144

日期：2026-09-15。接续第三十四轮的降强度定理。仍固定 h=1/4、α=sinc(h) 和原 λ=1 fresh 族。

## 1 主要结果

本轮在 η=18/125=0.144 构造了**严格的连续圆周鞅耦合**。构造通过有限矩阵找到，但其非负性和连续边缘由独立整数区间证书及解析修正保证，不以线性规划的成功标记为证明。

由第三十四轮，整个原仪器族在所有 η≤0.144 都有准备非情境实现，包含所有方向、真实准备混合及任意有限自适应协议。因此

\[
\boxed{0.144\le\eta_*\le\eta_B\approx0.144739232095}.
\]

相比上一阶段的 [0.1182,η_B]，未解区间宽度从约 0.0265392 缩为 **0.0007392321**。η_B 是否可达仍未解决。

## 2 数值探索暴露的网格问题

先在均匀圆网格上最大化内圆半径，求满足全部行列质量和条件向量均值的正矩阵。在 η_B：

| 网格数 | 离散问题最大半径 | 实际 α |
|:---|---:|---:|
| 32 | 0.9894883830 | 0.9896158370 |
| 64 | 0.9895842062 | 0.9896158370 |
| 128 | 0.9896079437 | 0.9896158370 |

这些离散问题看似都排除了实际 α，但对偶函数只是带截距的一维折线。将它们放回连续分布积分，得到的凸序差反而为正。这说明网格的支撑位置与积分误差会制造伪排除。

探索代码与结果分别保存在 [half_circle_lp_probe.py](half_circle_lp_probe.py)、[half_circle_lp_probe_results.json](half_circle_lp_probe_results.json)。它们不是连续存在性证书。

下面采用明确的方向：先将真实内圆分布随机送到包含它的有限多边形，再由有限耦合送到外圆弧的均值，最后恢复均匀外圆。每一步都严格保持条件均值。

## 3 连续输入到有限多边形的精确量化

固定

\[
N=256,\quad \Delta=\frac{2\pi}{N},\quad
h=\frac\pi N,\quad \rho=\frac{9897}{10000}.
\]

有限源点为 z_i=ρn_{iΔ}，另加 z_0^*=0。证书验证

\[
\rho\cos h-\alpha>9.6403921\times10^{-6}>0,
\]

因此该正多边形包含完整的 α 内圆。

对 ψ∈[iΔ,(i+1)Δ]，取三个非负重心概率

\[
q_i(\psi)=\frac{\alpha\sin((i+1)\Delta-\psi)}
{\rho\sin\Delta},\qquad
q_{i+1}(\psi)=\frac{\alpha\sin(\psi-i\Delta)}
{\rho\sin\Delta},\qquad
q_*(\psi)=1-q_i-q_{i+1}.
\]

前两项非负，而 q_i+q_{i+1}≤α/(ρcos h)≤1。它们满足

\[
q_i z_i+q_{i+1}z_{i+1}+q_*0=\alpha n_\psi.
\]

所以给定内圆点 Y=αn_ψ，按这些概率抽取有限源点 Z，便有 E[Z|Y]=Y。

令 w_i 为 ψ∼R_ηdψ/(2π) 时的源点质量。它不是简单地在网格上取 R/N，而是完整积分：

\[
w_i=\frac{\alpha}{2\pi\rho\sin\Delta}
\int_{-\Delta}^{\Delta}
R_\eta(i\Delta+u)\sin(\Delta-|u|)\,du.
\]

中心质量 w_*=1−Σ_iw_i。π 周期性保证 Σ_iw_i z_i=0。

为严格计算这些质量，定义

\[
I_0=2(1-\cos\Delta),\quad
I_k=\frac{2[\cos(2k\Delta)-\cos\Delta]}{1-4k^2},
\]

便有

\[
w_i=\frac{\alpha}{2\pi\rho\sin\Delta}
\left[I_0+2\sum_{k\ge1}c_k I_k\cos(2ki\Delta)\right].
\]

|I_k|≤I_0 可由正积分核直接得到。使用前 12 个模式，并用

\[
\sum_{k>K}c_k
=r^{K+1}\left[
\frac{1+2(K+1)q}{1-r}+\frac{2qr}{(1-r)^2}
\right]
\]

严格包住尾项。这里 q=√(1−η²)、r=(η/(1+q))²。中心质量约为 0.0000348407234993；没有将小质量删掉。

## 4 有限目标点到连续外圆

令目标点

\[
x_j=\beta n_{j\Delta},\qquad
\beta=\frac{\sin h}{h}.
\]

x_j 正好是外单位圆弧 [jΔ−h,jΔ+h] 的均匀均值。只要有限目标质量为 1/N，再在各自圆弧上均匀抽样，就恢复整个均匀单位圆 σ，且 E[X|x_j]=x_j。

现在有限问题是寻找 J_ij≥0，使

\[
\sum_jJ_{ij}=w_i,\qquad
\sum_iJ_{ij}=1/N,\qquad
\sum_jJ_{ij}x_j=w_i z_i.
\]

源索引 i 包含 256 个多边形顶点和中心，共 257 行；目标共有 256 列。这个辅助矩阵不是原系统的有限隐藏状态替代品：连续角仍被保留，精确旋转仍作用于连续角。

## 5 线性规划仅生成候选

使用 SciPy 1.18.1 的 HiGHS 接口，最大化共同正余量 ε，约束

\[
J_{ij}=Q_{ij}+\varepsilon w_i/N,\qquad Q_{ij}\ge0.
\]

数值结果 ε≈0.0000306889448。这样即使 Q 稀疏，J 的每个元素都保有小的正余量，可用于修复数值误差。

将候选的每个元素向下取整为 p_ij/2^60，p_ij 是非负整数，得到固定矩阵 P。它保存于 [polygon_coupling_candidate.json](polygon_coupling_candidate.json)，SHA-256 为

~~~
15b0efd2836f6e4495154b2413c66320fb5f04c9e382d92af4c5cf099c1868c8
~~~

实际复核的是这个整数矩阵，而非求解器的浮点残差。线性规划的参数、非负边界及等式接口按 [SciPy 官方 linprog 文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html) 核对；本轮没有把求解器输出当作数学证明。

## 6 精确修复所有约束

根据真实解析 w、z、x，定义 P 的行质量残差和行向量残差

\[
a_i=w_i-\sum_jP_{ij},\qquad
b_i=w_i z_i-\sum_jP_{ij}x_j.
\]

利用均匀多边形的恒等式

\[
\sum_jx_j=0,\qquad
\sum_jx_jx_j^T=\frac{N\beta^2}{2}I,
\]

先作行修正

\[
L_{ij}=\frac{a_i}{N}+
\frac{2b_i\cdot x_j}{N\beta^2}.
\]

P+L 的行质量和条件均值恰好正确。再令

\[
d_j=1/N-\sum_i(P_{ij}+L_{ij}),\qquad
\boxed{J_{ij}=P_{ij}+L_{ij}+w_i d_j.}
\]

由于源和目标的总质量都为 1、总向量均值都为零，

\[
\sum_jd_j=0,\qquad \sum_jd_jx_j=0.
\]

所以第二次修正只修复列质量，不破坏已修好的行质量和均值。这给出了**精确矩阵 J 的明确公式**，没有把残差“视为零”。

## 7 独立整数区间证书

[certified_intervals.py](certified_intervals.py) 使用分母 2^100 的整数区间。加减、乘除、平方根均向外舍入；输入 η、ρ 是精确分数；α=4sin(1/4)，π 用 Machin 恒等式及反正切交错级数包围。

Machin 恒等式可直接核对：若 a=atan(1/5)、b=atan(1/239)，则 tan(4a)=120/119，tan(4a−b)=1，结合角度范围即得 π=16a−4b。正弦、余弦使用 Taylor 多项式并加入明确的余项上界；Fourier 尾项使用第三节的精确求和界。所有最后几位的舍入均进入区间，没有依赖平台浮点 sin/cos 的误差假设。

令

\[
e_i=\frac{|a_i|+2\|b_i\|_1/\beta}{N},
\qquad
|L_{ij}|\le e_i.
\]

则

\[
|J_{ij}-P_{ij}|
\le e_i+w_i\left(
\max_j\left|1/N-\sum_iP_{ij}\right|+\sum_i e_i
\right).
\]

整数验证器重建全部解析区间，得到：

| 证书项目 | 严格界的十进制展示 |
|:---|---:|
| 原整数矩阵最小元素 | ≥4.1766598860×10^−12 |
| 任意元素修正的绝对值 | ≤3.9522339644×10^−16 |
| 修正后全部元素 | ≥4.1762646626×10^−12 |

最后一项的精确有理数下界是

\[
\boxed{
\min_{i,j}J_{ij}\ge
\frac{5294044406263354720}
{1267650600228229401496703205376}>0.
}
\]

因此 J 真正非负，且第六节解析公式保证所有约束严格成立。源质量、中心质量和多边形包含关系也分别有正区间证书。验证程序还会拒绝人为篡改一个矩阵元素的控制例。

数值求解器只负责寻找候选；此证书可以在不安装 SciPy 的环境中复核。

## 8 恢复完整 fresh 仪器

三步随机化依次为

\[
Y=\alpha n_\psi
\longrightarrow Z=z_i
\longrightarrow x_j
\longrightarrow X=n_x.
\]

第一步用重心概率，第二步用 J_ij/w_i，第三步在目标圆弧上均匀抽样。各步保持条件均值，边缘链又恰好连接 ν_η 与 σ，因此 E[X|Y]=Y。

若记第 i 个重心概率为 q_i(ψ)，则反向核在第 j 个输出圆弧上具有密度

\[
p(x\mid\psi)=\frac{N}{2\pi}
\sum_i q_i(\psi)\frac{J_{ij}}{w_i}.
\]

由第二十九轮的联合测度换元，正向结果核可明确写为

\[
\boxed{
K_s(y\mid x)=f_s(y)\,
p(x\mid\psi_s(y)).
}
\]

这不是只求出三个摘要矩。对原完整准备编码，有

\[
K_sF_z=F_{A_sz}
\]

作为整个选择性输出密度恒等式。代码通过精确输入圆弧积分、分段输出角积分和原 fresh 矩阵分别核对它，包含原子/圆弧辅助层的全部质量。

第三十四轮再将这个 η=0.144 的耦合作正旋转混合，得到全部较低强度。结合静默旋转和分支归纳，任意有限自适应、后选择协议及真实准备混合都有同一非情境编码。

## 9 复算、依赖与下一轮

主要文件：

- [polygon_martingale_certificate.py](polygon_martingale_certificate.py)：解析量化、候选生成、精确修正与验证。
- [polygon_coupling_candidate.json](polygon_coupling_candidate.json)：固定二进制分母的整数候选矩阵。
- [polygon_martingale_certificate_results.json](polygon_martingale_certificate_results.json)：严格余量与文件哈希。
- [certified_intervals.py](certified_intervals.py)：整数区间算术和余项控制。

默认复核不调用求解器：

~~~
python -X utf8 research_cognition_physics/strength_degradation.py --write-results
python -X utf8 research_cognition_physics/polygon_martingale_certificate.py --write-results
python -X utf8 -m unittest discover -s research_cognition_physics -p "*.py"
~~~

第三十四轮 7 项检查；第三十五轮 9 项仪器及证书检查、4 项区间算术检查。目录全部 **429 项检查通过**，Python 3.12.14、NumPy 2.3.5。

若重新生成数值候选，才需要可选 SciPy 1.18.1。本机已安装在项目的 .research_runtime 目录，没有修改系统 Python 包。可选依赖见 [requirements_research_optional.txt](requirements_research_optional.txt)：

~~~
python -m pip install --no-deps --target .research_runtime -r research_cognition_physics/requirements_research_optional.txt
python -X utf8 research_cognition_physics/polygon_martingale_certificate.py --solve --write-results
python -X utf8 research_cognition_physics/polygon_martingale_certificate.py --write-results
~~~

重新生成可能得到不同候选和哈希，必须重新通过独立验证。后续更高强度的候选应保存为新轮次文件，保留当前证书。

**下一项任务**：用更细的源多边形与同样的解析修正继续逼近 η_B，并检查证书余量随网格加密的变化；同时研究第三十三轮的正半圆耦合，寻找能在极限中消除辅助中心质量的解析构造。有限次网格成功不能单独证明 η_B 可达，当前严格区间仍为 [0.144,η_B]。

准备非情境性仍是新增条件，未从原认知原则推出；原经典解释及第七轮的独立横向保留问题均保持。上述区间改进不是量子理论或引力的推导。
