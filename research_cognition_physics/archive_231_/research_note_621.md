# 第621轮：原全Gauss热态的有限时间窗联合来源

日期：2026-10-01。接[620](research_note_620.md)，回查[590](research_note_590.md)、[591](research_note_591.md)、[603](research_note_603.md)与[612](research_note_612.md)。[代码](joint_smeared_gauss_noise.py)、[结果](joint_smeared_gauss_noise_results.json)、[核验](research_round_621_checks.json)、[条件账](unified_physics_condition_ledger_621.md)。两组复算、十二式；主代理审查，无新增独立代理审查。

## 1. 关闭哪一处缺口

620在原条件质量部门给出联合噪声的交换约束，但有限CAR矩阵不能证明完整无限维Gauss物质的来源二阶矩存在。603已有全模型热态、全部能量矩与相对形式界，结论只覆盖守恒涨落及零频响应，没有签收瞬时噪声。

**本轮证明：在603原固定有限图的同一Gauss热态中，有限时间窗涂抹后的共同形式来源具有有限二阶矩，多个来源的完整复二点矩阵共同收敛，并有显式能量截断尾界。** 没有另换谐模型、删除Majorana或重设每次截断的热态。

这为620要求的联合量子来源提供了一类全模型数学对象。它仍不是连续重整化应力、实际读取仪器或动态几何解。

|层次|地位|
|---|---|
|认知动机|几何、物质与能源应能共用一份合法的多来源关联|
|继承条件|603原完整玻色—规范—CAR有限图H、Gauss空间、正几何、β>0热态、来源形式界|
|本轮额外选择|固定非零时间窗；来源保持Gauss并属于已有共同形式族|
|解析成果|时间窗二阶矩界、联合GNS向量、截断尾界及热谱条件|
|数值核对|原中性Fock因子；另一个明确标记的抽象无限维形式例子|
|仍缺|连续UV、图细化统一常数、全局动态几何、有限幅响应与实际内部实现|

## 2. 原全模型的适用前提

完全复用603，不重新证明热态存在。固定原图和给定正几何，在Gauss空间取H=H_F、A=H+c≥1。每个方向G都必须是保Gauss的Hermitian形式，例如603已核过的几何参数来源；不能据此宣布任意未核的局部带荷算符也满足前提。

$$
\rho_\beta=Z^{-1}e^{-\beta H},\qquad
\langle A^k\rangle_\beta<\infty\ (k\ge1),\qquad
G=A^{1/2}\mathsf B A^{1/2}\ \hbox{作为形式},\quad
\mathsf B=\mathsf B^\dagger,\quad\|\mathsf B\|\le L.
\tag{1}
$$

H有离散谱，能级简并不影响以下证明。在联合本征基记a_n=E_n+c，p_n为同一原热权。形式G可能不能直接作用于每个能量本征矢，故不能先假定G²存在再积分。

## 3. 时间窗使形式来源成为有限热二阶对象

令实测试函数f∈W¹,¹(R)，可取紧支集；定义Fourier变换及谱涂抹

$$
\widehat f(\omega)=\int f(t)e^{i\omega t}dt,\qquad
(G_f)_{mn}=\widehat f(E_m-E_n)\sqrt{a_ma_n}\,\mathsf B_{mn},\qquad
M_0=\sup_\omega|\widehat f|^2,\quad
M_1=\sup_\omega|\omega|\,|\widehat f|^2.
\tag{2}
$$

M₀≤‖f‖₁²、M₁≤‖f‖₁‖f′‖₁：分部积分给|ω f̂|≤‖f′‖₁，再乘|f̂|≤‖f‖₁。由a_m≤a_n+|E_m−E_n|以及每列Σ_m|𝔅_mn|²≤L²，直接得到

$$
\boxed{\ \sum_{n,m}p_n|(G_f)_{mn}|^2
\le L^2\big[M_0\langle A^2\rangle_\beta+M_1\langle A\rangle_\beta\big]<\infty .\ }
\tag{3}
$$

此处没有假设统一能隙。近简并由〈A²〉控制，远能级由f̂控制；简并块也未删去。实际上任意与H对易且〈A²〉有限的正规态都可用同一证明；本项目使用的是603已存在的明确Gauss Gibbs态。不得推广为任意非平稳有限平均能量态。

一个归一的真实有限窗是三角形，半宽τ>0：

$$
f_\tau(t)=\frac{(1-|t|/\tau)_+}{\tau},\qquad
\widehat f_\tau(\omega)=\operatorname{sinc}^2(\omega\tau/2),\qquad
M_0=1,\quad M_1\le\frac2\tau.
\tag{4}
$$

sinc x=sin x/x。此f属于W¹,¹；若需要C∞窗口，可再与归一紧光滑非负核卷积，Fourier模不增，界仍有效、支撑相应增宽。窗口是声明的分辨率选择，不是设计好的内部时钟。

## 4. 联合对象与算符域的准确地位

用P_R=1_{A≤R}先构造G_f^R=P_RG_fP_R。式(3)说明G_f^Rρβ^(1/2)在Hilbert–Schmidt范数下收敛。各能量本征矢的列范数也有限；实f使有限本征矢核上的G_f对称，从而可闭。

本轮所需的是这份唯一的热GNS向量及其二点资料，而非未经证明的本质自伴性、谱测量或仪器实现。多个来源和各自窗口组成

$$
X_a=\lim_{R\to\infty}G_{a,f_a}^{R}\rho_\beta^{1/2},\qquad
\mu_a=\langle\rho_\beta^{1/2},X_a\rangle_{\rm HS},\qquad
C_{ab}=\langle X_a-\mu_a\rho_\beta^{1/2},X_b-\mu_b\rho_\beta^{1/2}\rangle_{\rm HS}.
\tag{5}
$$

C是Hermitian半正定矩阵；实部给对称联合噪声，虚部保留涂抹算符交换子的期望。所有交叉项由Cauchy–Schwarz共同收敛，不得像620已排除的接法那样保边际却任意删除交叉项。这里不声称这些二阶资料决定全部历史。

G和H都保Gauss，谱投影和涂抹也在同一Gauss部门中进行。有限R只规整来源算符，没有重新归一热权或换成另一个Gibbs态。

## 5. 可审计的截断尾界

记M_j(Ω)为式(2)上确界限制到|ω|≥Ω后的值。将遗漏对分成a_n>R/2，以及a_n≤R/2但a_m>R两组；后一组有|E_m−E_n|>R/2。得到

$$
\|(G_f-G_f^R)\rho_\beta^{1/2}\|_{\rm HS}^2
\le L^2\left[
M_0\langle A^2\mathbf1_{A>R/2}\rangle
+M_1\langle A\mathbf1_{A>R/2}\rangle
+M_0(R/2)\langle A^2\rangle+M_1(R/2)\langle A\rangle\right].
\tag{6}
$$

M₁(Ω)→0：f′∈L¹使ωf̂→0，f̂有界。这同时证明全谱极限，并允许对Gram矩阵误差用两边HS范数估计。三角窗还给显式多项式界

$$
M_0(\Omega)\le\frac{16}{\tau^4\Omega^4},\qquad
M_1(\Omega)\le\frac{16}{\tau^4\Omega^3},\qquad
\langle A^k\mathbf1_{A>R/2}\rangle
\le\left(\frac2R\right)^4\langle A^{k+4}\rangle.
\tag{7}
$$

最后所需第五、第六能量矩由603提供。这些常数通常很松，并依赖τ、β、图、几何和L。τ→0时界恶化，不准把有限窗结论直接当成瞬时有限性或对图细化一致性。

## 6. 与热谱及因果资料的连接范围

对一个Hermitian来源，定义正跃迁测度。即使其总质量可能无限，任何有限频带都有界：

$$
d\mu_G(\omega)=\sum_{n,m}p_n|G_{mn}|^2
\delta(\omega-(E_m-E_n))d\omega,\qquad
\mu_G([-\Omega,\Omega])\le L^2(\langle A^2\rangle+\Omega\langle A\rangle).
\tag{8}
$$

这是增长至多多项式的正分布。换n、m给有限温详细平衡

$$
d\mu_G(-\omega)=e^{-\beta\omega}d\mu_G(\omega)\quad(\omega>0).
\tag{9}
$$

该成熟热谱关系不是新定理；新增的是原完整形式来源的测度存在及可涂抹范围。它保留零频块。有限个涂抹来源的复Gram矩阵同时容纳噪声与交换子，但完整迟致核仍涉及时间排序、接触项和共同背景变分；本轮不从式(9)直接签收全部非线性因果响应。

## 7. 为什么形式界与热态还不等于瞬时噪声有限

以下是**抽象无限维例子，不是原Gauss模型的另一个反例**。取ℓ²(N₀)、H|n〉=n|n〉、A=H+1；对n≥1取v_n=√6/(πn)、v₀=0，‖v‖=1。令

$$
\mathsf B=|0\rangle\langle v|+|v\rangle\langle0|,\quad\|\mathsf B\|=1,
\qquad G=A^{1/2}\mathsf B A^{1/2},\qquad
|G_{n0}|^2=\frac6{\pi^2}\frac{n+1}{n^2}.
\tag{10}
$$

这是有效的相对形式；|λ|<1时A^(1/2)(I+λ𝔅)A^(1/2)正且与A可比，故不是用不存在的扰动来源作弊。对任意β>0，原H的热态所有能量矩有限，但p₀=1−e^(−β)>0，所以

$$
\langle G^2\rangle_{\beta,\rm cutoff}
\ge\frac{6p_0}{\pi^2}\sum_{n=1}^N\frac1n\longrightarrow\infty.
\tag{11}
$$

本例平均来源为零。三角窗后，完整二阶矩有限，N以后的正尾有直接界

$$
0\le V_\tau-V_{\tau,N}
\le\frac{16p_0(6/\pi^2)}{\tau^4N^4}.
\tag{12}
$$

证明使用(p₀+p_n)≤2p₀、(n+1)/n²≤2/n、|f̂τ(n)|²≤16/(τ⁴n⁴)，再用Σ_{n>N}n⁻⁵≤1/(4N⁴)。故仅凭603那组抽象形式／热矩前提不能删除时间窗；原模型是否另有更强瞬时域，需要额外证明。

## 8. 复算结果与下一接口

两组检查：

- 原中性16维Fock因子，τ=.05、.2、.8。直接分段时间积分与Fourier涂抹差≤1.12×10⁻¹⁶；联合协方差正性、热详细平衡和同态截断界均通过。τ=.2时能量—标量噪声交叉项≈.01598400；标量方差≈.02875840，解析二阶矩上界≈.60840014。此小矩阵检验归一、符号及界，**全图结论由第2—6节证明**。
- 抽象无限维例β=.7，瞬时截断二阶矩从N=32的2.12285增长到N=8192的3.82453，发散来自式(11)。τ=.2的有限窗值收敛到约1.68798630，N=8192的解析遗漏尾界6.80×10⁻¹³。该尾界是数学界，不含浮点求和舍入误差，未称区间认证。

本轮把C19同一完整热态、C22多来源涨落、C08／C10几何形式来源在**固定图、有限时间窗**下合并。没有增加新的物质模型，也没有将窗口选择认定为认知公理。该全图结果补强620，但不自动赋予完整有限图一个连续微分同胚Ward身份。

[622入口](round622_drafts/STATUS.md)接着核这份同模型的对称关联、交换子与接触项能否共同组成有限窗二阶影响泛函及因果响应；优先补源和过程的实际连接，不继续优化窗口，不设计读取装置。统一目标继续，动态几何、连续UV、GR和观测仍未完成。
