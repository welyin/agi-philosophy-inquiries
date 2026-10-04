# 第688轮：完整Gauss支撑与辅助结构留下的CAR边缘态

日期：2026-10-02。接[687](../../research_note_687.md)，以[成果归属补正](../../687/drafts/attribution_correction_653.md)校正其新增范围；直接继承[653](../../research_note_653.md)、[654](../../research_note_654.md)及[661](../../research_note_661.md)。[代码](../joint_gauss_support_marginal.py)、[结果](../joint_gauss_support_marginal_results.json)、[核验](../research_round_688_checks.json)、[条件账](../unified_physics_condition_ledger_688.md)。两组新核验、十六式，主代理审查，无独立代理审查。

## 1. 先校正归属，再回答新的问题

已回查README、research_direction、RESEARCH_STATE及最新保存结果。653已完成完整群正型、实际任意有限纯时间链、球谐谱与扩展Gauss正投影；687及688入口中的重叠部分不再作为新增定理。保留其独立复算，687的显式下界、全局帧和完整Haar精确值仍可使用。

本轮问题是：**扩展系统完成Gauss约束后，忽略辅助自由度，剩余CAR态是否就等于CAR自身的Gauss态？** 653的条件特征函数障碍发生在投影之前，不能直接用来回答投影后的问题。

答案在既有零标量、空间平凡纯时间分支上是否定的。本轮精确算出各辅助球谐层与原32CAR共同满足Gauss条件的重数；据此给出投影后CAR态与限定的CAR单独参考态之间的精确迹距离。没有反证原完整模型，也没有重新设计认知系统。

|层次|本轮地位|
|---|---|
|认知动机|共同约束、共同状态与子系统实际可读内容必须相容|
|继承建模输入|653原G、16通道、m₀=1、反周期时间、空间平凡、φ=0和S⁹测度|
|继承证明|653完整时间转移与谱；654的辅助脱耦极限；661物理来源的独立接口|
|解析新增|各层Gauss重数的精确证书；投影后CAR单态权重和迹距离|
|数值角色|整数模运算与有理数复算；不以浮点拟合判断态是否相同|
|未完成|原Weyl Y等时复合观测字典、完整动态Q₀正性、H_F身份和共同连续极限|

空间接口逐项继承[679表](../../research_note_679.md)：382上界、383自由连续反向对合、384已消去额外Lipschitz的下界；386非退化缩放接真实邻域及有覆盖／误差量词的有限证书；425为局部紧端点群、半幅及成本收缩的替代坐标路线；522全域绝对差坐标热下界；523实际仪器的共同条件性实现。本轮没有补上这些空间接口的任何新假设，也不将已完成的523共同实现重新列为不存在；待接的是当前h、s、CAR／Gauss对象身份。E仍是辅助变量，不是实际s记录。

## 2. 直接复用653的正转移

沿用原实际表示，不以687证书空间替换物理CAR：

$$
G=S(U(3)\times U(2)),\quad
R=\Lambda^{\rm even}V,\quad
\mathcal F=\Lambda(R\otimes\mathbb C^2),\quad \dim\mathcal F=2^{32}.
\tag{1}
$$

这里G与原Z₆商群同构，V为五维载体，detV=detR=1。653已证明的辅助转移及支撑为

$$
k(E,F)=\left(\frac{1+E\cdot F}{2}\right)^8,\quad
\mathcal K=\bigoplus_{\ell=0}^{8}\mathcal H_\ell,\quad
\mathbb T=\bigoplus_{\ell=0}^{8}\mu_\ell I_{\mathcal H_\ell},\quad
\mu_\ell>0,\quad \mu_1/\mu_0=8/17.
\tag{2}
$$

谱与维数逐项读取并核对653保存的有理数结果，不重复证明Funk–Hecke公式。令U_F与U_K分别为原群表示，则

$$
P_G=\int_G U_F(g)\otimes U_K(g)\,dg,\qquad
Z_N^{\rm raw}=2^{-32N}\operatorname{Tr}\{P_G(I_{\mathcal F}\otimes\mathbb T^N)\}.
\tag{3}
$$

式(3)及其正性均属653。N是原时间格数；这一条件分支没有把原全图玻色和规范动态积分进去。

## 3. 新增的各层Gauss支撑

定义原CAR与每层辅助表示共同形成的单态数，以及辅助自身的单态数：

$$
m_\ell=\dim(\mathcal F\otimes\mathcal H_\ell)^G
=\int_G\chi_{\mathcal F}(g)\chi_\ell(g)\,dg,\quad
s_\ell=\dim\mathcal H_\ell^G,\quad n_0=m_0=\dim\mathcal F^G.
\tag{4}
$$

Haar积分取乘积而非任意共轭替代；辅助球谐实表示的特征标为实。五个环面本征值z满足乘积为1。以偶子集A标记R的权r_A，原实际32CAR特征标为

$$
\chi_{\mathcal F}(z)=\prod_{A\ {m even}}(1+r_A)^2
=\prod_{A\ {m even}}(2+r_A+r_A^{-1}),\quad
r_A=\prod_{j\in A}z_j,\qquad
\chi_\ell(z)=h_\ell(z,z^{-1})-h_{\ell-2}(z,z^{-1}).
\tag{5}
$$

第二个等号使用∏r_A=1，不把R暗中换成共轭物种；约定负次数h为0。完整原群Weyl积分给

$$
12m_\ell=\operatorname{CT}\left[\chi_{\mathcal F}\chi_\ell\Delta\right],\quad
12s_\ell=\operatorname{CT}\left[\chi_\ell\Delta\right],\quad
\Delta=\prod_{\substack{i<j\le3\ \text{或}\ i=4,j=5}}
(1-z_i/z_j)(1-z_j/z_i),\quad z_5=(z_1z_2z_3z_4)^{-1}.
\tag{6}
$$

CT为四个独立变量的常数项。对ℓ≤8，χ_F各变量的绝对次数最多8，χ_ℓ最多8，Δ分别最多3、3、3、2，因此完整乘积的次数上界为(19,19,19,18)。20阶单位根的四重平均准确提取常数项，没有混叠。复用687两素数算法，模数乘积超过每个非负常数项的上界12·2³²·dimH_ℓ，故唯一恢复整数。不是从小数四舍五入得到重数。

G在C³和C²球面乘积上传递；实多项式不变量由两块模平方生成。次数2k的对称张量有k+1个不变量，减去乘总半径平方的次数2k−2部分后，谐波只剩一个；奇次没有。因此解析地有

$$
s_\ell=\begin{cases}1,&\ell\text{ 为偶数},\\0,&\ell\text{ 为奇数},\end{cases}
\qquad
(m_0,\ldots,m_8)=(228248,1511416,4996416,11074440,18658304,25771104,30876104,33661096,34779760).
\tag{7}
$$

s_ℓ也由独立常数项计算核对。所有m_ℓ均为4的倍数，符合原中性两个spin模式贡献；整数证书、素数和次数界均保存在结果中。用新重数独立恢复687已有完整Haar值：

$$
2^{-32}\sum_{\ell=0}^{8}m_\ell\mu_\ell
=\frac{261746352167}{17416264183971840}.
\tag{8}
$$

此相等是精确有理数检查；该积分值不再计为688新增结果。

## 4. 投影后，辅助电荷仍能补偿CAR电荷

在653重建的正表示中定义联合态和其CAR边缘态：

$$
\rho_N=\frac{P_G(I_{\mathcal F}\otimes\mathbb T^N)}{\sum_\ell m_\ell\mu_\ell^N},\qquad
\rho_{\mathcal F,N}=\operatorname{Tr}_{\mathcal K}\rho_N,\qquad
\Pi_{\mathcal F}=\int_G U_F(g)\,dg.
\tag{9}
$$

T与原G对易，分子为正；长度归一在正规化时相消，未改原未归一权重。Pi_F是原CAR规范单态子空间上的投影，是有界、规范不变的偶CAR算符。在有限32模式全CAR代数中有相应多项式表示，无需测量E。

若CAR已在单态空间，整体Gauss条件只要求辅助自身为单态。因此精确有

$$
(\Pi_{\mathcal F}\otimes I)P_G
=\Pi_{\mathcal F}\otimes\Pi_{\mathcal K},\qquad
\Pi_{\mathcal K}=\int_G U_K(g)\,dg,
\quad
p_N:=\operatorname{Tr}(\Pi_{\mathcal F}\rho_{\mathcal F,N})
=\frac{n_0\sum_{\ell\ {m even}}\mu_\ell^N}{\sum_\ell m_\ell\mu_\ell^N}.
\tag{10}
$$

比较态明确限定为“无辅助、零CAR Hamiltonian、只施加CAR自身Gauss条件”的正规化态：

$$
\sigma_{\mathcal F}=\Pi_{\mathcal F}/n_0,\qquad
\Pi_{\mathcal F}\rho_{\mathcal F,N}\Pi_{\mathcal F}=p_N\sigma_{\mathcal F},\qquad
[\rho_{\mathcal F,N},\Pi_{\mathcal F}]=0.
\tag{11}
$$

后一对易也可直接从整体投影得到。Pi_F内部的差为−(1−p_N)σ_F，正交补内则是正态块、迹为1−p_N；两个块支撑正交。所以

$$
\frac12\|\rho_{\mathcal F,N}-\sigma_{\mathcal F}\|_1=1-p_N
=\frac{\sum_\ell(m_\ell-n_0s_\ell)\mu_\ell^N}{\sum_\ell m_\ell\mu_\ell^N}>0\qquad(N<\infty).
\tag{12}
$$

严格大于零来自m₁>0、s₁=0、μ₁>0；所有系数m_ℓ−n₀s_ℓ非负，因为单态×单态是联合不变空间的子空间。故这不是若干N样本的猜测，而是所有有限正整数N的结论。

|N|与限定CAR自身Gauss参考态的迹距离|
|---:|---:|
|1|0.903972090300753|
|2|0.694130341684419|
|3|0.451359419152173|
|4|0.258341922675597|
|8|0.0157017490123955|
|16|0.0000383016801159834|

这些小数只展示结果，计算和相等检验全部使用有理数。

## 5. 与旧连续极限兼容，不能越界成物理反例

654已有原辅助谱隙在a→0、N=β/a时脱耦的结论，本轮不重证它。对上式可直接继承谱比：

$$
0<1-p_N\le\frac{\sum_{\ell=1}^{8}(m_\ell-n_0s_\ell)}{n_0}
\left(\frac8{17}\right)^N\longrightarrow0.
\tag{13}
$$

有限尺度不同与受控极限趋同没有冲突。这里没有用N拟合现实参数，也没有新增连续时空证明。源于653核的本征值不允许每轮重新选择。

特别要区分三个层次：

$$
\text{完整联合Gauss单态}\ \not\Longrightarrow\quad
\text{CAR子系统本身为Gauss单态},\qquad
\rho_{\mathcal F,N}\ne\sigma_{\mathcal F}\ \not\Longrightarrow\quad
\rho_{\mathcal F,N}\ne\text{原完整 }H_F\text{ 的CAR边缘态}.
\tag{14}
$$

原H_F还有玻色和规范自由度，它们也可能补偿CAR电荷；不能把比较态σ_F偷换成原H_F的物理态。此结果也不说明整体规范条件失败。

此外，本轮Pi_F可观测是在653重建CAR表示中成立。661、670、679已经给出原Weyl Y的来源字典、自由物理正性、正则延伸及反射合同；尚不能只凭相同配分函数把所有等时CAR复合算符与原Y插入认作同一个操作：

$$
Z_N^{\rm reconstructed}=Z_N^{\rm original}
\quad\not\Longrightarrow\quad
\langle\Pi_{\mathcal F}\rangle_{\rm reconstructed}
=\langle\mathcal O_Y\rangle_{\rm original}
\ \text{对未经证明的 }\mathcal O_Y.
\tag{15}
$$

因此688没有原Q₀反射负方向，没有否证整个手征候选，更不否证认知原则蕴涵命题。它给下一次来源匹配一个明确待保持的投影量。

## 6. 复算、条件合并和下一接口

两组核验分别是完整原群精确重数、CAR边缘态有理数证书；688入口原矩阵复算归入653继承，不新增组数。复用既有Python／NumPy；joint_gauss_support_marginal.py省略参数只读复算，--write-results仅用于不存在结果时首次保存。核验同时复算687、检查旧谱和全部保护文件散列，不做图像检查。

本轮将C01、C14、C16、C19、C20的同一分支状态与约束支撑具体接通；C22中的原完整Y来源身份仍单列。四个分支及全局目标不变：

$$
\text{原完整 }H_F\ ;\quad
\text{指定连续物质}\ ;\quad
\text{给定作用的经典几何}\ ;\quad
\text{原手征辅助候选}
\qquad\text{尚未共同等同}.
\tag{16}
$$

[689入口](../../689/drafts/STATUS.md)先回查646、660—661、670、679的接触项与原物理观测，核纯时间分支的**全部来源及等时复合观测**是否映射到653转移CAR；若能，则把688的投影量接到实际Y。一般场重定义、旧质量拉回和既有时间谱不再独立发表；只有新增的同一对象接口可计进展。
