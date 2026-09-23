# 第386轮：从渐近位移到真实邻域，以及有限尺度的三维读出证书

日期：2026-09-23。独立完整轮次，科学基线冻结至385。12项检查通过。

## 1. 新增结果与旧结果的边界

**认知动机。** 实际有限位移未必服从一个全局齐次群；若结构只在较小尺度上趋于稳定，方向读数还能否约束真实位置维数？

本轮接入已有的strong dilatation结构定理，给出两项连接：

1. 在同一真实位置集合上，非退化切向距离与原距离局部具有相同拓扑。因此所得局部Lie结构提供的是实际位置邻域的完整小球壳，而不只是另一个抽象切空间中的球面。
2. **一个有限尺度就可以给出三维上界。** 实际近似反向的无迹读数分离，若大于位置偏差经过读出变化模放大后的误差，就必然分离切向精确反向。无需读出随尺度收敛，也无需绝对可见度保持非零。

这仍是条件定理。完整位置上的强缩放合同、非退化性、可实施的归壳、误差模及全域覆盖都没有由F＋U＋C＋P推出。只有另满足384的**精确**方向下界合同，才能合成三维，而不是仅有n≤3。

| 已有轮次 | 复用、不重复计成果 | 本轮增量 |
|---|---|---|
| [385](research_note_385.md) | 全局精确缩放群的Lie分类、Lyapunov截面 | 将群运算放到局部极限；证明局部群拓扑确实连接原位置拓扑 |
| [383](research_note_383.md) | 完整球形边界自由对合与qubit上界 | 实际近逆不必本身是对合；用有限误差转移到精确切向逆 |
| [382](research_note_382.md) | 已认证覆盖、概率误差和球面网格 | 同时控制近逆、归壳条件与尺度相关读出；允许无绝对信号极限 |
| [384](research_note_384.md) | 全族协变、对比完整及极小重定向给下界 | 不把近似协变或几个采样点的满秩冒充该定理的精确前提 |

全目录去重未发现此前完整使用strong dilatation定理。数学分类及维数不变性归已有文献；代码的已知三维校准空间不是由数值重新发现的物理维数。

## 2. 成熟定理的准确输入

核读[Buliga，*Infinitesimal affine geometry of metric spaces endowed with a dilatation structure*，Houston J. Math. 36(1)，2010，§3、Theorem 6.2及Corollary 6.3](https://imar.ro/~mbuliga/buliga10.pdf)。本文选缩放群Γ＝(0,∞)。输入是完备、局部紧的真实位置度量空间(X,d)；为全局采用通常的拓扑流形定义，再明确要求X可分。先不假定X是流形或已有维数。

必须使用完整定义：A0保证局部定义域、像与跨基点复合相容；每个基点的缩放是可逆局部同胚，满足同基点尺度复合；它们联合连续并一致收缩。距离和跨基点差操作分别满足紧集上一致极限：

$$
\begin{gathered}
\varepsilon^{-1}d(\delta^x_\varepsilon u,\delta^x_\varepsilon v)
\longrightarrow d^x(u,v),\qquad d^x\ \text{非退化},\\
\Delta^x_\varepsilon(u,v)
:=\delta^{\delta^x_\varepsilon u}_{1/\varepsilon}
       \delta^x_\varepsilon v
\longrightarrow\Delta^x(u,v).
\end{gathered}
\tag{1}
$$

上行是A3，下行是A4；极限的统一性、非退化性和A0不是可省的技术修饰。式(1)不要求有限尺度的差操作本来就是群差，也不保证实际度量在有限位移下严格齐次。

相较385，这换成一份**不同的**几何操作合同，不是只删除一个假设：现在明确有度量、局部紧性、完整跨基点缩放和一致极限。它没有把整个X预设为一个自由传递位置群，但仍不是从认知压缩或单点的小脉冲控制自动得到的结构。

## 3. 切向拓扑为什么确实对应真实位置

先在原d拓扑中选足够小的紧邻域K。由A3，d^x是K×K上的连续函数；非退化使它成为K上的距离。因此恒等映射为连续双射，源紧、靶Hausdorff：

$$
\operatorname{id}:(K,d)\longrightarrow(K,d^x)
\ \text{为同胚}.
\tag{2}
$$

这份局部论证足够核实所需拓扑连接，不依赖把一般一致距离收敛误读成自动的全空间一致同胚。

Buliga的Theorem 6.2在**同一个U(x)**上给出范数化局部锥群；Corollary 6.3给出局部正分次Lie群。结合式(2)，可以把结论写成：

$$
(U(x),\Sigma^x,\delta^x) \text{是局部锥群},\qquad
U(x)\ \text{在x附近同胚于}\ \mathbb R^{n_x}.
\tag{3}
$$

这里只作局部识别，不把U(x)任意全局化成一个位置群。正分次也不等于Carnot分层：后者还要求第一层生成整个Lie代数。更没有由拓扑结论给出原d的欧氏形式或体积增长指数。

在局部Lie图中选一个足够小的、求逆对称的指数椭球闭域D_x，使它和涉及的操作都留在公共定义域。385的Lyapunov截面证明可在这里局部使用：实际缩放的线性生成元取适应正定二次型，缩小时保持椭球内域。其边界L和有限尺度像满足：

$$
L=\partial D_x\simeq S^{n_x-1},\qquad
\iota=\text{局部群求逆}|_L,\qquad
\iota^2=\operatorname{id},\quad\iota u\ne u,\qquad
j_\varepsilon:=\delta^x_\varepsilon|_L.
\tag{4}
$$

小域包含一个完整位置开邻域；δ是其上的同胚，所以jε(L)是实际邻域δᵡε(D_x)的完整拓扑边界。它**不必是物理距离d下的等距球**。在椭球参数中取普通单位球弦长度量ρ，则ι是ρ等距的对径；这仅为误差记账选择方便的壳坐标，不把该弦长宣布为物理距离。

## 4. 实际近逆、归壳和物理误差的归一化

由实际缩放复合可执行候选逆操作：

$$
i^x_\varepsilon(u)=\Delta^x_\varepsilon(u,x)
\longrightarrow\iota u,\qquad
R_\varepsilon(u)=\Pi\big(i^x_\varepsilon(u)\big)\in L,
\qquad
b_\varepsilon=\sup_{u\in L}\rho(R_\varepsilon u,\iota u).
\tag{5}
$$

Π是公共紧环带上的径向归壳。原iε通常不在L，不可直接把它当作完整边界的反向。定性上，A4和局部拓扑一致性使iε一致趋近精确逆；它最终远离基点，连续归壳使bε→0。

**操作上仍须交代Π。** 需要根据已校准的端点选择一个基点缩放比例并实际实施；若端点或比例有误差，要一并计入bε。包含1／ε的复合和归壳不是免费控制。A4本身既不提供有限收敛速度，也不保证已有装置可以知道所需比例。数值证书必须有独立误差界或校准。

若实际完成归壳后，仅掌握物理端点距离误差β_ε，则应先去掉尺度。令αε为式(1)第一行在公共紧集上的上界误差，κ_x为从d^x到ρ的紧壳连续模，便有：

$$
\sup_u d(j_\varepsilon R_\varepsilon u,j_\varepsilon\iota u)\le\beta_\varepsilon
\quad\Longrightarrow\quad
b_\varepsilon\le
\kappa_x\!\left(\frac{\beta_\varepsilon}{\varepsilon}+\alpha_\varepsilon\right).
\tag{6}
$$

κ_x(t)可定义为d^x距离≤t的壳点对之最大ρ距离，紧性保证它趋零。仅β_ε→0不够；整个邻域本来就在缩小。还须控制相对尺度的误差，以及下节中的读出敏感度。

## 5. 主命题：一个有限尺度即可传递到三维上界

固定一个有限ε。实际qubit效果以完整壳端口为标签，连续分解为标量和无迹部分：

$$
E_\varepsilon(u)=a_\varepsilon(u)I+f_\varepsilon(u)\cdot\sigma,
\qquad 0\le E_\varepsilon(u)\le I,\qquad
\omega_\varepsilon(t)=
\sup_{\rho(u,v)\le t}|f_\varepsilon(u)-f_\varepsilon(v)|.
\tag{7}
$$

允许效果标量aε(u)随方向变化，主命题只使用三个无迹实系数。不能把包含迹变化的总效果差代替这份三维向量差。

记实际近逆的全域无迹分离为gε。三角不等式直接给出：

$$
\begin{gathered}
g_\varepsilon:=\inf_{u\in L}
|f_\varepsilon(u)-f_\varepsilon(R_\varepsilon u)|,\\
\inf_{u\in L}|f_\varepsilon(u)-f_\varepsilon(\iota u)|
\ge g_\varepsilon-\omega_\varepsilon(b_\varepsilon)>0
\quad\Longrightarrow\quad n_x\le3.
\end{gathered}
\tag{8}
$$

最后一步复用383的自由对合Borsuk–Ulam结果：连续fε:L→R³分离每对精确逆，因而维数不能大于3。也可把fε写成固定迹合法效果I／2＋fε·σ，因为原效果正性保证|fε|≤1／2。这里只用数学上的三分量映射和概率差，没有宣称存在一个免费物理通道抹去原标量信息。

式(8)没有要求Rε自己两次严格返回，也没有要求所有ε共享某个读出极限。只要**某一个**有限尺度满足严格正裕量即可。

若信号尺度vε→0，足够的相对条件例如：

$$
\liminf_{\varepsilon\to0}\frac{g_\varepsilon}{v_\varepsilon}>0,
\qquad
\frac{\omega_\varepsilon(b_\varepsilon)}{v_\varepsilon}\longrightarrow0.
\tag{9}
$$

它使式(8)最终成立，无需vε有正下界。还可以使用更一般的严格正差liminf。重要的是相对裕量，不能把“信号没有消失”当成惟一可行合同，也不能只让几何误差单独趋零。

## 6. 有限概率证书：覆盖、去迹和成本

用六个已校准的±Pauli探针准备。每个端口的无迹分量可直接由同一端口的概率差恢复：

$$
f_{\varepsilon,j}(u)=\frac{p_{j,+}(u)-p_{j,-}(u)}2,
\qquad
|\widehat p-p|\le\tau
\ \Longrightarrow\
|\widehat f_\varepsilon-f_\varepsilon|\le\sqrt3\,\tau.
\tag{10}
$$

τ应包括读出和已知准备偏差；不能靠重复采样消除未知系统误差。两个端口的向量差因此扣除2√3τ。直接使用单个任意准备的概率差会混入标量差，本轮不那样做。

假定已**独立认证**u_i覆盖整个L，覆盖半径h；实际只需在这些设置执行Rε，并知道各样点到精确ιu_i的ρ误差≤b。因ι为ρ等距，从最近样点得到：

$$
\begin{gathered}
\widehat g=\min_i|\widehat f_\varepsilon(u_i)
                  -\widehat f_\varepsilon(R_\varepsilon u_i)|,\\
\inf_{u\in L}|f_\varepsilon(u)-f_\varepsilon(\iota u)|
\ge\widehat g-2\sqrt3\tau
          -\omega_\varepsilon(h)-\omega_\varepsilon(h+b).
\end{gathered}
\tag{11}
$$

证明只比较u与最近u_i，以及ιu与Rεu_i；因此不需要另猜Rε具有统一Lipschitz常数。若已知读出Lε-Lipschitz，最后两项可用Lε(2h＋b)替代。这个严格正下界认证的是**精确反向分离**，虽只实施了近似反向。

网格只是读取设置；有限qubit数据不能证明它已经覆盖全部未知位置。代码复用382在已给S²上的114点网格及h≤π／8，不将这份模拟网格当作未知现实边界的覆盖证明。

若有m个端口，各自与其实际近逆端口均作六探针测量，则共有12m个Bernoulli设置。对于独立重复读出、准备无额外偏差的特定协议，Hoeffding与并合界给充分预算：

$$
N_{\rm each}\ge
\frac{\log(24m/\alpha)}{2\tau^2},\qquad
\Pr(\text{任一概率误差}>\tau)\le\alpha.
\tag{12}
$$

若τ按vε同比缩小，这个证书的充分采样预算按vε⁻²增加。不是所有量子测量策略的必需下界；本轮只模拟有界校准扰动，没有执行所列的大量采样。

## 7. 从单点上界到整个位置分量：还差什么

式(8)或式(11)只给n_x≤3。要等于3，仍需在同一完整壳上满足384：实际联合连续可逆重定向、每条轨道稠密、全族**精确**酉协变、效果差空间为Herm₀(2)。采样点满秩、近似协变或可控制任意内部酉都不能替代这些条件。

若几何strong dilation合同在同一X的每点成立，式(3)给各点局部Euclidean邻域。邻域交叠中的维数由[Hatcher，*Algebraic Topology*，§2.1，Theorem 2.26，作者全文](https://pi.math.cornell.edu/~hatcher/AT/ATch2.pdf)的维数不变性一致，故n_x局部常数。在连通分量C中：

$$
n_x\ \text{局部常数}\quad\Longrightarrow\quad
n_x\equiv n_C;
\qquad
\big[\text{某一点满足上界证书和384下界}\big]
\Longrightarrow\ C\ \text{为拓扑三维流形}.
\tag{13}
$$

X为可分度量空间确保Hausdorff和第二可数。这里要求每点几何合同及同一真实位置身份，不是把一个点的量子实验外推到任何未知远方。也未构造由操作决定的相容光滑图册、共同物理度量或跨尺度自然传播规律。

## 8. 独立数值校准：有限差不是已经存在的群差

为检查实际复合、归壳和误差传递，给一个完整显式模型。它的X＝R³与欧氏d是**数值校准输入**，不承担一般几何定理的证明，也不借预设曲流形来展示一个已知切空间。

令κ＝0.4，e₂为第二坐标基向量。基点相关的非线性剪切图与实际缩放为：

$$
\begin{gathered}
A(x)=\kappa\sin x_1,\qquad
\chi_x(y)=(y-x)+A(x)(y_1-x_1)^2e_2,\\
\chi_x^{-1}(q)=x+q-A(x)q_1^2e_2,\qquad
\delta^x_\varepsilon y=\chi_x^{-1}(\varepsilon\chi_x(y)).
\end{gathered}
\tag{14}
$$

各基点图都是全局剪切同胚，故同基点缩放复合与逆精确成立；公共小定义域可按A0限制。A(x)和导数一致有界，紧集收缩和跨基点域相容可控。该例的切向距离与实际距离误差可直接算出；记U＝χ_x(u)、V＝χ_x(v)：

$$
d^x(u,v)=|U-V|,\qquad
\left|\varepsilon^{-1}d(\delta^x_\varepsilon u,\delta^x_\varepsilon v)
-d^x(u,v)\right|
\le |A(x)|\varepsilon|U_1^2-V_1^2|.
\tag{15}
$$

极限群差是χ_x⁻¹(V−U)，极限群和是χ_x⁻¹(U＋V)。把真实复合的q＝δᵡεu代入，可以独立给出A4的有限误差：

$$
\left|\Delta^x_\varepsilon(u,v)-\chi_x^{-1}(V-U)\right|
\le |q-x|+|A(x)|\varepsilon|V_1^2-U_1^2|
       +\kappa\varepsilon(1+|U_1|)|V_1-U_1|^2
\quad(0<\varepsilon<1).
\tag{16}
$$

此界来自|A(q)−A(x)|≤κ ε|U₁|，不是数值拟合。有限和Σᵡε＝δᵡ₁/ε δ^(δᵡεu)_ε v在本例一般不结合，极限和才结合；代码同时核验这一区别。没有把有限位移暗中当成385的全局精确群。

在基点x＝(0.7,−0.2,0.1)取完整壳χ_x(u)＝n、|n|＝1。直接复合得到原近逆的坐标：

$$
q_\varepsilon(n):=\chi_x(i^x_\varepsilon(u))
=(\varepsilon-1)
\big[n+\{A(x+\varepsilon n)-A(x)\}n_1^2e_2\big],
\qquad R_\varepsilon n=\frac{q_\varepsilon(n)}{|q_\varepsilon(n)|}.
\tag{17}
$$

归壳通过同基点δ及比例1／|qε|实现。该Rε一般不为精确对合，但已知界足以用于式(8)。有端点图坐标误差e、|e|≤η时：

$$
\begin{gathered}
|q_\varepsilon(n)|\ge(1-\varepsilon)(1-\kappa\varepsilon)=:m_\varepsilon>0,\\
\eta<m_\varepsilon\quad\Longrightarrow\quad q_\varepsilon(n)+e\ne0,\qquad
\sup_n\left|\frac{q_\varepsilon(n)+e}{|q_\varepsilon(n)+e|}+n\right|
\le2\kappa\varepsilon+\frac{2\eta}{m_\varepsilon}.
\end{gathered}
\tag{18}
$$

这里e是端点图坐标误差，不是量子效果的无迹误差。上面的归一化不等式对实际非零向量成立；η<mε则保证误差球内所有实施都非零。若η≥mε，即使公式返回一个大于2的空泛上界，也不保证操作可定义。ε接近1使原近逆逼近基点，小绝对误差可被归壳放大；代码用一个实际非零向量核验这一点。

对这份模型的壳坐标n、m，真实有限尺度距离还有双向界：

$$
\varepsilon(1-2\kappa\varepsilon)|n-m|
\le d(j_\varepsilon n,j_\varepsilon m)
\le\varepsilon(1+2\kappa\varepsilon)|n-m|,
\qquad 2\kappa\varepsilon<1.
\tag{19}
$$

它直接显示物理误差必须除以尺度后才是方向误差。原点附近不同环带的归一化条件与测量精度均有自己的代价，不能只记物理端点越来越接近。

## 9. 无非退化读出极限的正例与高频反例

正例的效果取：

$$
E_\varepsilon(n)=\tfrac12I+\tfrac{v_\varepsilon}2
\big[O(\log\varepsilon)n\big]\cdot\sigma,
\qquad v_\varepsilon=\varepsilon^2,\qquad
\tau=v_\varepsilon/100.
\tag{20}
$$

O是绕第三轴的旋转。效果本身趋于无信息的I／2，而归一化方向图无限旋转，没有固定方向极限。精确逆差为vε；使用真实近逆式(17)、误差界式(18)与有限证书，仍可保持正的相对裕量。此例的具体H作用若用于384下界，仍须作为实际重定向控制另行给出；式(14)的缩放权限不自动提供全部旋转。

必要的反面对照针对的是**读出变化模**。在圆周壳上取偶数k，精确逆ιθ＝θ＋π，实际近逆R_kθ＝θ＋π＋π／k，并取合法的无迹读出：

$$
f_k(\theta)=\tfrac12(\cos k\theta,\sin k\theta,0),\qquad
b_k=\pi/k\to0,\qquad
f_k(\iota\theta)=f_k(\theta),\quad
|f_k(R_k\theta)-f_k(\theta)|=\omega_k(b_k)=1.
\tag{21}
$$

实际位置的反向偏差越来越小，但读出频率同步增大，于是读数偏差一直最大；严格判据正确地不给证书。这是“忽略读出变化模”的反例，不是满足完整三维合同的高维反模型。R_k两次操作还会转过2π／k，体现近似返回也不等于精确对合。没有把此圆例另算一次隐藏标签或消相干研究。

## 10. 复算结果

[代码](asymptotic_displacement_readout_audit.py)、[结果](asymptotic_displacement_readout_audit_results.json)、[本轮核验](research_round_386_checks.json)。核验链接由主代理冻结时生成。既有Python 3.12.14和NumPy 2.3.5，无安装、无图像；12项通过后执行`--write-results`。

```powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'research_cognition_physics/archive_231_/asymptotic_displacement_readout_audit.py'
```

| ε | 真实复合差的误差 | 解析上界 |
|---:|---:|---:|
| 0.2 | 0.14773764 | 0.25607951 |
| 0.1 | 0.07397231 | 0.12766949 |
| 0.05 | 0.03701471 | 0.06374391 |

| ε | vε | 有限证书下界 | 下界／vε | 每设置独立采样的充分预算，α＝0.01 |
|---:|---:|---:|---:|---:|
| 0.25 | 0.0625 | 0.02855467284 | 0.4568747654 | 16,024,861 |
| 0.125 | 0.015625 | 0.00795238670 | 0.5089527486 | 256,397,773 |
| 0.0625 | 0.00390625 | 0.00208846580 | 0.5346472454 | 4,102,364,354 |

这些大数是式(12)的充分预算，未实际采样。代码加入已知界内的模拟概率扰动以检查证书。固定τ＝10⁻⁵、ε＝0.001时，证书下界变为−3.26319×10⁻⁵；并非分离物理上已消失，而是这份精度不能认证它。

高频对照k＝2000时，位置角偏差约0.001570796，近逆读数差仍为1，精确逆读数差解析为0；矩阵／三角计算残差约3.2×10⁻¹³。代码另检查有限和不结合、原近逆非对合、未知迹的六探针消除、双端概率误差、归壳环带裕量与样本计账。21个编号公式连续。

## 11. 推进与剩余接口

这次把“有抽象切向球面”推进为“完整真实位置邻域具有相同局部拓扑”，再给出有限尺度读出的严格操作入口。只需相对误差足够小，而非假定无限尺度下还有固定绝对信号；但定性收敛本身没有免费提供可认证速度。

下一步应集中追问所需几何条件的操作来源，特别是：为什么实际端点具有非退化一致缩放极限；怎样独立验证完整覆盖、归壳及变化模；以及384的完整方向重定向合同如何在实际传播中取得。它们仍是明确新增条件，不能用模型校准代替认知推导，也不能把条件性拓扑三维结论当作度量、SR或GR已经完成。
