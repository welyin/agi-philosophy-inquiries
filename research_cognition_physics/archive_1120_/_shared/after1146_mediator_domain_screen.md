# 1146后资格：正常量子中介乘积交互的共同域与加权来源桥

2026-10-12。独立资格推导，未分配1147。复用[1129正常来源](../1129/proof.md)、[1144连续中介](../1144/proof.md)和[1146来源/CAR资格](../1146/source_car_review.md)，响应[当前方向](after1146_direction_decision.md)。本文件只处理固定平滑模型的共同域、资源与来源比较，不证明局部材料已经形成1144的两份独立逻辑记录。

**资格结论：这项乘积交互的域问题可以关闭，且来源替换只需中介的第一能源矩。** 可以保完整正常Fock与Q/R，不必把中介截断。长时限不必造成来源误差发散；但控制源的正常准备必须覆盖实际因果范围，代价不能仍记作1129的一秒部署。真正剩余接口是局部CAR密度与1144固定独立逻辑源的共同实现。

## 1. 角色、固定输入与记号

材料空间 $\mathcal H_D$ 为完整正常Dirac粒子/空穴表示，$H_D\ge0$。不作正频投影或删掉海。另有两个独立正常场：

- $\phi$：来源/执行场，仍可采用1129的无质量KG场，$H_\phi\ge0$；它在比较准备中提供坡道均值。
- $\chi$：真正量子中介，质量 $\mu>0$，$H_\chi=d\Gamma(\omega)$，$\omega(k)=\sqrt{|k|^2+\mu^2}$；其态及与两份材料的相关性始终保留。

固定 $\epsilon,\ell>0$。$\phi_\epsilon(x)$ 使用旧控制场核；令实光滑紧支 $\rho_\ell$ 为中介核，采用酉Fourier约定，

$$
\chi_\ell(x)=a(v_x)+a^\dagger(v_x),\qquad
v_x(k)=\frac{e^{-ik\cdot x}\widehat\rho_\ell(k)}{\sqrt{2\omega(k)}}.
\tag{1}
$$

所有材料密度 $J_i(x)=J_i(x)^\dagger$ 是有限CAR偶元，可为旧有限旋量/味核。假设强可测且

$$
j_i(x):=\|J_i(x)\|,\qquad
C_i:=\int |\eta_i(x)|j_i(x)\,d^3x<\infty .
\tag{2}
$$

$\eta_i$ 是固定实紧支部署函数；可把实耦合常数并入 $F_i$。$F_i$ 为实有界Lipschitz函数，$M_i=\|F_i\|_\infty$、$L_i=\operatorname{Lip}F_i<\infty$。需写测试力时再取有界连续导数；需要强时间传播时，比较系数在有限时间窗取 $C^2$ 或有限分段 $C^2$ 且匹配 $C^1$。这不要求材料 $J_i$ 彼此对易。

不同因子上的 $F_i(\phi_\epsilon/A)$、$J_i$、$\chi_\ell$ 对易。模型为

$$
H_0=H_D+H_\phi+H_\chi,\qquad
V=\sum_i\int\eta_i(x)F_i(\phi_\epsilon(x)/A)J_i(x)\chi_\ell(x)\,d^3x,
\qquad H=H_0+V.
\tag{3}
$$

这是明确有限分辨率的新材料合同。控制坡道、固定站点、各场类型、核和准备均为输入；没有用一个全局控制算符在每点免费复制。$\phi_\epsilon(x)$ 是同一正常场的不同局部涂抹，近乎均匀的均值必须由整片正常准备支付资源。

## 2. 可直接使用的场界

置

$$
a_\ell:=2\|\omega^{-1/2}v_0\|
=\sqrt2\|\widehat\rho_\ell/\omega\|,\qquad
b_\ell:=\|v_0\|
=\|\widehat\rho_\ell/\sqrt{2\omega}\|,
$$
$$
d_\ell:=\|\omega^{1/2}v_0\|
=\|\rho_\ell\|_2/\sqrt2 .
\tag{4}
$$

固定平滑核和 $\mu>0$ 使它们有限，且

$$
a_\ell\le\frac{\sqrt2}{\mu}\|\rho_\ell\|_2,\qquad
b_\ell\le\frac1{\sqrt{2\mu}}\|\rho_\ell\|_2.
$$

由湮灭算符的标准能源界和CCR：

$$
\|a(v_x)\Psi\|\le
\|\omega^{-1/2}v_0\|\|H_\chi^{1/2}\Psi\|,
$$
$$
\|a^\dagger(v_x)\Psi\|
\le\|\omega^{-1/2}v_0\|\|H_\chi^{1/2}\Psi\|+b_\ell\|\Psi\|.
$$

因此在完整Fock上、对全部其它因子及任意被动参考都有

$$
\boxed{\|\chi_\ell(x)\Psi\|
\le a_\ell\|H_\chi^{1/2}\Psi\|+b_\ell\|\Psi\|.}
\tag{5}
$$

这里 $b_\ell$ 保留真空涨落，不能因中介初始为真空就删掉。该场界与固定平滑玻色耦合的标准相对有界论证一致；例如[Gubinelli–Hiroshima–Lőrinczi的原论文](https://arxiv.org/abs/1304.6662)使用正则化玻色源框架。本文件没有调用其去截断结论。

## 3. 共同算符域、半有界与实际总能源

记

$$
M:=\sum_i C_iM_i,\qquad
\alpha:=a_\ell M,\qquad \beta:=b_\ell M,\qquad
K:=\frac{\alpha^2}{2}+\beta .
\tag{6}
$$

式(5)使(3)先成为从 $D(H_\chi^{1/2})$ 到 $\mathcal H$ 的强向量积分，满足

$$
\|V\Psi\|\le\alpha\|H_\chi^{1/2}\Psi\|+\beta\|\Psi\|.
$$

它不再是1129的有界交互，不可照搬其算符范数结论。但对 $\Psi\in D(H_0)$ 及任意 $\zeta>0$，

$$
\|H_\chi^{1/2}\Psi\|^2
\le\|\Psi\|\|H_0\Psi\|,
$$
$$
\|V\Psi\|\le
\zeta\|H_0\Psi\|+
\left(\frac{\alpha^2}{4\zeta}+\beta\right)\|\Psi\|.
\tag{7}
$$

积分对称，因此Kato–Rellich给

$$
D(H)=D(H_0),\qquad H=H^\dagger.
\tag{8}
$$

各项保完整正常Fock，不要求核去截断。对归一向量另有形式界

$$
|\langle V\rangle|\le\alpha\sqrt{\langle H_\chi\rangle}+\beta
\le\frac12\langle H_\chi\rangle+K,
$$
$$
\boxed{H\ge H_D+H_\phi+\frac12H_\chi-K\ge\frac12H_0-K.}
\tag{9}
$$

上下形式比较也给 $D((H+K+1)^{1/2})=D((H_0+1)^{1/2})$。故实际自治演化保持共同算符域和共同形式域。初始总自由能源有限时，

$$
\boxed{
\langle H_0\rangle_t\le
2\left(\langle H_0\rangle_0+\alpha\sqrt{E_{\chi,0}}+\beta+K\right).
}
\tag{10}
$$

同时同一H保完整能源账

$$
\Delta\langle H_D\rangle+\Delta\langle H_\phi\rangle+
\Delta\langle H_\chi\rangle+\Delta\langle V\rangle=0.
\tag{11}
$$

这是独立资源证明，未用小态距推无界能源。若输入在算符域，则总能源二阶矩有限并由实际H图范数保持；若只需均能和下述误差，可在共同形式域工作，不把二阶总能源矩升级为物理门槛。

## 4. 时间比较过程及严格域依据

设自由来源态 $a_t=e^{-itH_\phi}a$ 在整个部署管内满足

$$
\left\|(\phi_\epsilon(x)/A-\tau(t,x))a_t\right\|
\le\nu(t,x).
\tag{12}
$$

定义外开关参考系数 $f_i(t,x)=F_i(\tau(t,x))$ 和完整材料—中介参考H：

$$
\widehat H(t)=H_D+H_\chi+
\sum_i\int\eta_i(x)f_i(t,x)J_i(x)\chi_\ell(x)\,d^3x.
\tag{13}
$$

这份参考保 $\chi$ 为真实量子场；被替换为均值比较的只有执行场 $\phi$。不得把量子中介也同时换成期望值。

若 $f_i$ 及前两阶时间导数在 $\int|\eta_i|j_i\,dx$ 的加权范数中有界连续，则(7)统一给共同算符域 $D(H_D+H_\chi)$，相应形式在共同形式域 $C^2$ 且均匀半有界。加常数后形式范数等价；常形式域传播定理给唯一酉 $\widehat U(t,s)$、共同形式域及算符域保持。可直接核对[Balmaseda–Lonigro–Pérez-Pardo，Theorem 3.10](https://arxiv.org/pdf/2112.11063)：一阶形式正则性给形式域传播，二阶给强解及算符域保持。有限分段光滑情形在共同域拼接，不在接缝偷偷更换H。

本条件比仅写“每时刻自伴”完整。1144的多项式开关在闭区间 $[0,T]$ 内本来光滑；区间外零延拓仅 $C^1$ 不阻碍该有限窗。一般 $C^1$ 开关可先用形式域传播，并由同一能源界和光滑近似得到下文加权比较。

## 5. 中介第一能源矩：无需高阶 Grönwall 的直接界

令

$$
c(t):=\sum_i\int|\eta_i(x)|j_i(x)|f_i(t,x)|\,d^3x,\qquad
E_\chi(t):=\langle\widehat U(t,0)\Psi,H_\chi\widehat U(t,0)\Psi\rangle .
$$

在共同算符域先计算。材料、来源系数与中介算符分属不同因子，故

$$
[H_\chi,\chi_\ell(x)]
=a^\dagger(\omega v_x)-a(\omega v_x).
$$

对任意有界自伴其它因子 $B$，$[B,a]=0$。CCR使交换子期望的两项互为共轭，因此

$$
\left|i\langle B[\chi_\ell,H_\chi]\rangle\right|
\le2\|B\|\,\|a(\omega v_x)\Psi\|
\le2d_\ell\|B\|\sqrt{E_\chi}.
\tag{14}
$$

此处是期望界：无需另加创造算符的真空常数；向量范数界(5)仍必须保 $b_\ell$。积分给

$$
|E_\chi'(t)|\le2d_\ell c(t)\sqrt{E_\chi(t)}.
$$

以 $\sqrt{E_\chi+\varepsilon}$ 正则化后令 $\varepsilon\downarrow0$，得到

$$
\boxed{
\sqrt{E_\chi(t)}
\le\sqrt{E_{\chi,0}}+
d_\ell\int_0^t c(u)\,du .
}
\tag{15}
$$

此界允许材料内部演化、$\chi$ 与材料纠缠、不同 $J_i$ 不对易；没有把 $\chi$ 当作经典受迫振子。通过共同域的有限能源近似，它扩张到第一能源矩有限的输入。$\Psi$ 可以连同任意R；常数不随R维数变化。

实际自治过程同理，只将 $c(t)$ 换成 $M$，因为实际 $F_i(\phi/A)J_i$ 仍为有界自伴、与 $\chi$ 对易：

$$
\sqrt{E_\chi^{actual}(t)}
\le\sqrt{E_{\chi,0}}+d_\ell Mt.
\tag{16}
$$

不必先估计整个巨大控制源能源才能控制中介均能。

若需一个常规总能源Grönwall证书，令 $M_*=\sup c(t)$、
$K_*=(a_\ell M_*)^2/2+b_\ell M_*$，
$\mathcal E(t)=\langle\widehat H(t)\rangle+K_*+1\ge1+\langle H_D+H_\chi\rangle/2$，以及

$$
c_1(t)=\sum_i\int|\eta_i|j_i|\partial_t f_i|\,dx .
$$

显含时间功身份和(5)给

$$
|\mathcal E'(t)|\le
c_1(t)(a_\ell\sqrt{2\mathcal E(t)}+b_\ell)
\le c_1(t)(\sqrt2a_\ell+b_\ell)\mathcal E(t),
$$
$$
\mathcal E(t)\le\mathcal E(0)
\exp\!\left[(\sqrt2a_\ell+b_\ell)\int_0^tc_1(u)\,du\right].
\tag{17}
$$

它控制参考过程的总第一能源矩；来源误差采用更直接的(15)，避免把不相关的巨大材料或来源能源带入误差。

## 6. 保整个来源后态的加权 Duhamel

谱演算只与标量均值比较，由(12)得

$$
\|[F_i(\phi_\epsilon(x)/A)-f_i(t,x)]a_t\|
\le L_i\nu(t,x).
\tag{18}
$$

定义

$$
\kappa(t):=\sum_i\int|\eta_i(x)|j_i(x)L_i\nu(t,x)\,d^3x .
$$

对完整参考轨道 $\widehat U(t,0)\Psi\otimes a_t$，来源差与材料—中介向量因子分立；用(5)和(15)：

$$
\|(H-\widehat H(t)-H_\phi)(\widehat U(t,0)\Psi\otimes a_t)\|
\le\kappa(t)\left[
a_\ell\left(\sqrt{E_{\chi,0}}+d_\ell\int_0^tc(u)\,du\right)+b_\ell\right].
$$

共同算符域上对实际与参考酉逐时求导，再以形式域近似扩张，即得完整联合输出界

$$
\boxed{
\left\|e^{-iTH}(\Psi\otimes a)
-\widehat U(T,0)\Psi\otimes a_T\right\|
\le\int_0^T\kappa(t)
\left[a_\ell\left(\sqrt{E_{\chi,0}}+
d_\ell\int_0^tc(u)\,du\right)+b_\ell\right]dt .
}
\tag{19}
$$

这是包含材料、量子中介、执行场后态及相位的向量界。来源自由态只是参照，真实后态并未被迹掉或固定。若编码输入满足统一的中介第一能源上界 $E_{\chi,0}\le E_*$，则(19)对所有未知Q及任意被动R一致，给相同的全输出等距上界与半迹距/能量限制通道误差；再接同一仪器和继续也不会放大。

为明确这不是单一输入态的估计，设完整正常编码为
$\mathcal I:\mathcal H_Q\to\mathcal H_D\otimes\mathcal F_\chi$，并要求

$$
\mathcal I^\dagger\mathcal I=I,\qquad
\|H_\chi^{1/2}\mathcal I\|^2\le E_*.
\tag{19a}
$$

将(19)右端中的 $E_{\chi,0}$ 换为 $E_*$，得到同一个常数 $\Delta(T)$，则

$$
\left\|
e^{-iTH}(\mathcal I\otimes|a\rangle)
-\big(\widehat U(T,0)\mathcal I\big)\otimes|a_T\rangle
\right\|\le\Delta(T).
\tag{19b}
$$

这是编码空间上的完整等距算符范数；张量任何 $I_R$ 不改变它。固定中介真空时 $E_*=0$，无论未知Q及其与R的纠缠如何，条件一致成立。固定有限能源中介准备也同理；一般与Q相关的准备则必须先核(19a)，不能只核某个状态的均能。

无需控制参考轨道的 $H_\chi^2$ 或更高矩。证明时在 $D(H_0)$ 计算，是技术密集域；最终误差通过有界能源近似只要求中介第一矩。

**不可去掉能源限定而宣称全部中介输入的普通diamond界。** 最简单的单个非零线性场差已说明：任意 $\delta\ne0$ 时，连续谱场算符满足 $\|\exp(-i\delta\chi)-I\|=2$。小系数不导致无能源限定的算符范数小误差。已有1144从真空中介出发不受这一量词障碍影响。

1129坡道条件 $\tau=t$、$\nu\le3/(\epsilon A)$ 下，若所有 $F_i$ 相同且 $C=\sum_iC_i$，(19)给可直接代入的粗界

$$
\delta_\phi(T)\le
\frac{3CL}{\epsilon A}
\left[(a_\ell\sqrt{E_*}+b_\ell)T+
\frac12a_\ell d_\ell CM_FT^2\right].
\tag{20}
$$

这比无中介的旧公式多出的正是实际中介能源权，不能把旧的有界 $J$ 结论原样乘一个无界 $\chi$。

## 7. 直接对接1144慢开关：长时限并不强迫误差发散

取1144的

$$
f(u)=30u^2(1-u)^2\,1_{0<u<1},\qquad
F_i(z)=\lambda_iT^{-1/2}f(z/T).
\tag{21}
$$

有 $\|f\|_1=1$、$f_{\max}=15/8$、$\operatorname{Lip}f=10/\sqrt3<6$，且对称性给

$$
\int_0^1(1-u)f(u)\,du=\frac12.
$$

这里让实际坡道满足 $\phi_\epsilon/A\simeq t$，因此 $F_i(t)$ 恰是相同的慢开关，而不是把原过程偷偷压回 $T=1$。设

$$
K_\lambda=\sum_i C_i|\lambda_i|,\qquad L_f=6,\qquad
\nu\le\sigma_\phi/A.
$$

则

$$
c(t)=K_\lambda T^{-1/2}f(t/T),\qquad
\kappa(t)\le\frac{\sigma_\phi L_fK_\lambda}{AT^{3/2}}.
$$

把这些精确代入(19)，得到

$$
\boxed{
\delta_\phi(T)\le
\frac{\sigma_\phi L_f}{A}
\left[
\frac{K_\lambda(a_\ell\sqrt{E_*}+b_\ell)}{\sqrt T}
+\frac12a_\ell d_\ell K_\lambda^2
\right].
}
\tag{22}
$$

有限 $\ell,\epsilon,C_i,\lambda_i$ 固定时，上界随 $T$ 不发散，且可通过有限 $A$ 达到任何预定正容差。第一项保真空涨落，第二项付出中介在参考过程中的能源增长；无需假定中介全过程接近真空。该粗界不利用1144较精细的绝热返回，因而即使 $J_i$ 不对易也仍成立。

然而实际控制准备有不可省略的代价。若 $\bigcup_i\operatorname{supp}\eta_i\subset B_R$，要使 $0\le t\le T$ 上自由平滑均值一直为 $At$，初始控制源 $\pi_{cl}=A\beta_0$ 的平台须覆盖

$$
B_{R+T+\epsilon}.
$$

可取支撑在更大半径 $R+T+\epsilon+1$ 的平滑平台；它仍是正常Fock准备，且

$$
E_\phi(0)=\frac{A^2}{2}\|\beta_0\|_2^2
\le\frac{2\pi}{3}A^2(R+T+\epsilon+1)^3.
\tag{23}
$$

不能在 $T=10^9$ 的合同里继续复用1129的 $B_8$ 平台和 $1073A^2$。三维红外可积性及高频平滑性仍成立，但预部署的体积和准备是新的有限资源输入；这不证明经济性、自治准备或SR运动运输。

## 8. 来源和中介反作用各自保留

在测试函数和共同域/形式的意义下，独立CCR给两条不同的真实力：

$$
\frac d{dt}\Pi_{\chi,H}(g)
=\chi_H(\Delta g-\mu^2g)
-\sum_i\int\eta_i(x)(\rho_\ell*g)(x)
[F_i(\phi_\epsilon(x)/A)J_i(x)]_H\,dx ,
\tag{24}
$$

$$
\frac d{dt}\Pi_{\phi,H}(h)
=\phi_H(\Delta h)
-\frac1A\sum_i\int\eta_i(x)(\rho_\epsilon*h)(x)
[F_i'(\phi_\epsilon(x)/A)J_i(x)\chi_\ell(x)]_H\,dx .
\tag{25}
$$

中介力的交互项有界；执行场力的交互项只对 $H_\chi^{1/2}$ 相对有界，不能再称全算符有界。式(5)、(16)可逐时控制其向量范数和期望，再积分给有限反冲。光滑核给 $\omega v_0$ 等所需测试函数的有限范数，域论证可由有限粒子芯上的Weyl正则化及统一形式界完成。它不要求中介态被删除或经典化。

若各核和 $J_i(x)$ 家族共同空间平移协变，总自由空间动量的形式责任式为

$$
\frac d{dt}P_{j,H}
=-\sum_i\int(\partial_j\eta_i)(x)
[F_i(\phi_\epsilon(x)/A)J_i(x)\chi_\ell(x)]_H\,dx .
\tag{26}
$$

右端同样受中介第一能源矩控制，固定 $\eta_i$ 的支撑力未消失。若 $J_i$ 含其它固定空间结构，还须相应增加其显含位置责任；不能只列(26)便宣称全装置四P闭合。

## 9. 仍未证明的材料连接及正式准入门槛

式(19)/(22)把**同一个局部材料—量子中介模型**的外开关替成正常执行场，未把它自动等同1144。

查重已回看[1144后正常X/Y自治中介资格](after1144_autonomous_spatial_source_screen.md)：它已经证明固定独立A/B和正常中介下的控制/位置加权Duhamel及完整Q/R，且明示直接NR钟与既有SR身份不相容。故“把外开关换成内部控制器”本身不是新发现。本文件只补针对局部非对易完整CAR密度、独立执行场与量子中介的乘积交互域及统一能源权；在实际材料匹配通过前，它只是一份联合资格，不够单独登记完整新轮次。

1144的精确线性源相位依赖固定材料算符 $A,B$、独立性和对易关系。这里完整局部密度 $J_i(x)$ 在 $H_D$ 下会演化，不同点的受控模式可能不对易，并可扰动正常海。仅有 $\|J_i\|$ 积分界，不能推出其压缩就是固定 $A_i\rho_i(x)$，也不能套用旧Magnus闭式得到AB记录相位。

下一项真正的最小问题是：给正常材料编码及实际区域读口，沿完整参考轨道证明局部密度耦合与目标独立逻辑源之间的**能源加权等距残差**；保留所有海/材料余部，而非只比较 $J_i$ 的期望。将这一残差积分后，才可与(19)、1144的相位/残余预算相加。中介有限涂抹还会改变有效空间剖面；须明确实际剖面与旧 $\rho_i$ 是否相同，不能不核对就借用旧Yukawa系数。

本资格可支持正式准入的前半项：

1. 固定核下共同正常域、半有界和同H有限能源：已给。
2. 从外开关到正常执行来源、保完整Q/R及两量子场后态：给出可代入的一阶能源加权界。
3. 中介和执行来源各自真实反作用：给出同H测试力及资源界。
4. 实际局部CAR材料能否完成独立的两个不对易问题及继续：本文件未证。
5. 原CO位置/钟、共同boost、动态部署与末读设备的四P、点态材料微因果：本文件未证。

所以不应再把本候选的正常域当成需要连续多轮维修的方向门槛。应先检第4项是否有实际可行的材料身份；若失败，仅是此连接/材料合同失败，不能反证整个认知纲领。

## 10. 最小复算与证明层级

本文为资格解析推导，没有新正式轮次、科学校准或认知公理。未模拟连续场。下列使用既有Python即可复算(22)的两个精确积分：

    from fractions import Fraction as Q
    coefficients = [Q(0), Q(0), Q(30), Q(-60), Q(30)]
    area = sum(c / Q(n + 1) for n, c in enumerate(coefficients))
    weighted = sum(c * (Q(1, n + 1) - Q(1, n + 2))
                   for n, c in enumerate(coefficients))
    assert area == 1
    assert weighted == Q(1, 2)
    assert Q(100, 3) < 36  # (10/sqrt(3))^2 < 6^2
    print(area, weighted)

已只读运行：输出为1与1/2。解析证明承担无限Fock、共同域、完整量词和能源权；标准传播定理承担相应存在性，文献一般定理不替代未完成的逻辑材料匹配。
