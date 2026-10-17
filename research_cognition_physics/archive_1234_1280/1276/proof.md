# 1276：同父集体关系的条件参考计数与有效边界面积

2026-10-17。作者证明；验收状态见 [independent_review.md](independent_review.md)。本轮复用 [1175完整过程](../../archive_1161_1200/1175/proof.md) 和其 [保存结果](../../archive_1161_1200/1175/results.json)，不更改原准备、Hamiltonian、普通参数、复制数或探针菜单。新增量是明确声明的条件参考计数，不是新的热学公理。

## 1. 固定同一过程、边界和参考

原空间图为 $65^3$ 周期格点，$\ell=16/65$。取事前固定的 $32^3$ 站点立方区域 $A$；其内部边界站点指有原最近邻边通向补集的站点，每站只计一次：

$$
M=|\partial_{\rm sites}A|=32^3-30^3=5768,
\qquad n_{\rm cut}=6\cdot32^2=6144.
\tag{1}
$$

$M$ 与切边数不同，角、棱站点不重复计入 $M$。沿用1175齐次有效胞元度规的六面面积比较器：

$$
A_g(t)=A_0a(z(t))^2,\qquad
A_0=6(32\ell)^2=\frac{1572864}{4225},\qquad
 a(z)=\frac{1+3z/20}{1+7z/20}.
\tag{2}
$$

面积由旧传播主符号的独立有效几何字典给定，不由下述计数反定义。式(2)不是本轮新造的实际面积仪器，也不把旧1225的材料与读口证书搬入1175。

每个边界站点 $i$ 有原来的 $K$ 个关系qubit。令 $\Pi_{i,m}$ 为恰有 $m$ 个 $+Z$ 的扇区投影，$\bar Z_i=K^{-1}\sum_{r=1}^KZ_{ir}$，于是其本征值为 $z_m=2m/K-1$、维数为 $d_m=\binom Km$。复制标号 $r$ 不是空间方向。声明条件参考及结构算符为

$$
\sigma_{i,m}=\frac{\Pi_{i,m}}{d_m},\qquad
\sigma_{\boldsymbol m}^{\partial A}
 =\bigotimes_{i\in\partial_{\rm sites}A}\sigma_{i,m_i},\qquad
L_K=\sum_{i\in\partial_{\rm sites}A}\sum_{m=0}^K
 \log\binom Km\,\Pi_{i,m}.
\tag{3}
$$

各联合扇区的权重来自同一真实完整后态 $\rho_{\rm act}(t)$；只评价

$$
S_{\rm ref}(t)=\operatorname{Tr}\rho_{\rm act}(t)L_K,
\qquad s_{\rm ref}(t)=\frac{S_{\rm ref}(t)}{MK}.
\tag{4}
$$

$\sigma$ 是额外公开的参考选择，不是原全Hamiltonian的Gibbs边缘。式(4)是条件参考熵的期望，不是实际状态的von Neumann熵、扇区概率的Shannon熵，亦不包含后者。实际纯相干态的固定扇区条件态可以是纯Dicke态，而 $\sigma_{i,m}$ 通常混合。这里不执行测量、去相干、替换、制备或重置；源、关系、未知probe/R及全部交互均保留在旧完整演化中。

## 2. 有限计数与真实关联态的运输

以下对数均为自然对数。记 $h_2(p)=-p\log p-(1-p)\log(1-p)$，端点以连续极限定义。有限类型计数给

$$
\frac{e^{Kh_2(m/K)}}{K+1}\le\binom Km\le e^{Kh_2(m/K)},
\quad
0\le h_2(m/K)-\frac1K\log\binom Km
 \le\frac{\log(K+1)}K.
\tag{5}
$$

一个直接证明是：在成功概率 $m/K$ 的二项分布中，第 $m$ 项是最大项；其概率介于 $1/(K+1)$ 与1。$m=0,K$ 同样成立。

令 $\Phi_K(t)$ 为原**精确、可纠缠**背景演化。旧 [均场界](../../archive_1161_1200/_admission/after1174_direction/mean_field_bound.md) 式(13)作用于单cluster算符 $Z_i$，给每个站点

$$
\| (\bar Z_i-z(t))\Phi_K(t)\|
 \le r(t):=\frac{4e^{6wt}}{\sqrt K}.
\tag{6}
$$

这不是全背景接近Hartree乘积态的声明。对 $\bar Z_i$ 在 $\Phi_K$ 中的谱概率分布，Cauchy给 $\mathbb E|\bar Z_i-z|\le r$。二元熵连续性、凹性和 $h_2$ 在 $[0,1/2]$ 的单调性给出（$r/2\le1/2$）

$$
\begin{aligned}
\left|\left\langle h_2\!\left(\frac{1+\bar Z_i}{2}\right)\right\rangle
 -h_2\!\left(\frac{1+z}{2}\right)\right|
&\le\mathbb E h_2\!\left(\frac{|\bar Z_i-z|}{2}\right)\\
&\le h_2\!\left(\frac{\mathbb E|\bar Z_i-z|}{2}\right)
\le h_2(r/2).
\end{aligned}
\tag{7}
$$

因此没有要求不同站点或复制在末态独立。旧同页式(14)还给完整准备等距的算符范数误差

$$
\eta(t)=\frac{4\mathcal L t e^{6wt}}{\sqrt K}
 =\mathcal L t\,r(t).
\tag{8}
$$

该界对旧指定背景准备、任意原许可probe与任意未知伙伴 $R$ 一致；比较器仍保留精确 $\Phi_K$。其输出半迹距至多 $\eta$。由 $0\le L_K/(MK)\le\log2$，结合式(5)—(8)，

$$
\boxed{\left|s_{\rm ref}(t)-h_2\!\left(\frac{1+z(t)}2\right)\right|
 \le b_K(t):=\frac{\log(K+1)}K+h_2(r(t)/2)+\eta(t)\log2.}
\tag{9}
$$

这里平均后没有额外的 $M$ 因子。范数界保留真实探针反作用；不为任意未知背景或任意后选择分支授予同一保证。原模型的统一正谱基线及全部 $K$ 份资源仍在父Hamiltonian内，计数没有免除这些成本。

## 3. 直接继承原 $T=2$ 终点，不重算动力学

采用1175原参数 $\epsilon=3/20,\nu=7/20,\Omega=1/2,m_S=2$、原初态和

$$
T=2,\qquad K=10^{52275242},\qquad
\mathcal L=7094593.59375>7\times10^6.
\tag{10}
$$

旧保存证书给 $\eta(T)\le10^{-4}$。旧 [数值残差证明](../../archive_1161_1200/_admission/after1174_direction/mean_field_numerical_bound.md) 式(4)定义的加权Bloch误差含 $|g-\widehat g|^2$，故其终点界确实控制 $z=g_z$，不是仅控制探针概率。保存结果为

$$
\widehat z(2)=-0.1140832906148773,\qquad
|z(2)-\widehat z(2)|\le0.0012623025381518402,
\quad -\frac{29}{250}<z(2)<-\frac{14}{125}.
\tag{11}
$$

最后两个有理端点留有严格余量；本轮以保存十进数再向外留 $10^{-15}$ 作有理核验。此继承保留旧证书公开的IEEE基本运算／向外残差包络假设，不提升为逐硬件指令形式验证。

无需构造巨整数 $K$ 或新积分：由式(8)、$T=2$，

$$
\frac{r(T)}2=\frac{\eta(T)}{2\mathcal LT}
 <\frac1{280000000000}<2^{-38},\qquad
h_2(r(T)/2)<\frac{39}{2^{38}}<10^{-8}.
\tag{12}
$$

这里用 $h_2(x)\le x(1-\log x)$ 及 $\log2<1$。对 $n=52275242\ge12$，$\log(K+1)/K<(3n+1)/10^n$，而后一序列严格递减，因此小于 $37/10^{12}<10^{-8}$。同样的界控制 $t=0$（其实际 $\eta(0)=0$）。故可统一保守取

$$
 b_K(0),b_K(T)<\frac1{10000}+\frac2{10^8}
 =\frac{5001}{50000000}<\frac{101}{10^6}.
\tag{13}
$$

令 $f(z)=\log2-h_2((1+z)/2)$，则 $f(0)=f'(0)=0$ 且 $f''(z)=1/(1-z^2)\ge1$，二次积分即给 $f(z)\ge z^2/2$。从 $z(0)=0$ 得到真实完整过程的有限结论

$$
\begin{aligned}
s_{\rm ref}(0)-s_{\rm ref}(2)
&>\frac{98}{15625}-2\frac{5001}{50000000}\\
&=\frac{151799}{25000000}=0.00607196
>\frac3{500}.
\end{aligned}
\tag{14}
$$

这是非零有限两端净差，不声称整段逐点单调。式(14)是原固定 $K$、原指定背景准备和原probe/R域的统一结论；中等 $K$ 诊断不能替代这个大 $K$ 有限证明。

## 4. 同一有效面积与原钟

因为 $a'(z)=(\epsilon-\nu)/(1+\nu z)^2<0$，式(11)给

$$
 a(2)>a(-14/125)=\frac{1229}{1201},\qquad
\frac{A_g(2)-A_g(0)}{A_0}
>\frac{68040}{1442401}>\frac{47}{1000}.
\tag{15}
$$

另有一个受限的正向充分性：给定本参考选择和此齐次同父族，$a=\sqrt{A_g/A_0}$ 与 $z$ 一一对应，

$$
z=\frac{1-a}{\nu a-\epsilon},\qquad
\left|s_{\rm ref}-h_2\!\left(\frac12+\frac{1-a}{2(\nu a-\epsilon)}\right)\right|\le b_K.
\tag{16}
$$

所以这里不是几何与计数毫无联系，而是有限误差下有明确的非线性关系；它只在本齐次族及已声明参考内成立。关系计数与面积在上述同一原过程、相同两端上反向变化。若要求本参考量满足固定 $\lambda>0$ 的 $S_{\rm ref}=\lambda A_g+C$，该两端已不相容；即使仅比较归一量，任何固定正系数的零修正增量关系也不成立。此删除仅限式(3)明确选择的参考候选，不删除用户机制的其它参考／修正，也不把正面积律升级为新公理。

沿用原共同钟 $\tau(t)=\int_0^tN(z(s))ds$，$N=1+3z/20\in[17/20,23/20]$，故两端在同一原钟上有严格顺序。保存值 $\tau(2)\simeq1.8851275756769434$ 的误差是旧独立钟预算；本轮不把 $t$ 暗换成新钟，也不从旧终端移动读数自动取得瞬时面积的新增读口。式(15)是旧有效几何比较器的有限结论，不是对有限 $K$ 量子度规算符或全Lorentz面积测量的认证。

## 5. 核验与停止线

[check.py](check.py) 不导入或重跑旧ODE。它读取并哈希绑定旧证据，精确核终点区间、两种边界计数、幂指数预算和式(14)—(15)，穷举中等 $K$ 的有理类型概率；另外用非乘积向量的谱概率作式(7)诊断。后者仅查实现错误，全域量词由上面证明及旧误差定理承担。结果见 [results.json](results.json)。

净新增是：在已有自主关系—来源—传播父过程中，明确条件参考计数可运输到真实有限 $K$ 后态，并与独立有效面积作有限比较。它既不重报一般类型定理，也不拼接1274的Casimir载体。没有签收该参考的热学熵身份、完整Gibbs温标、穿界热流、局部Hilbert应力、Einstein动力学或强生成猜想的反例。完成此候选筛查即停止，不沿它扩建制备、面积或温度设备。
