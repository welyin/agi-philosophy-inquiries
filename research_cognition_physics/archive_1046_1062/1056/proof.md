# 1056 解析证明：电子二标签与总角能流

2026-10-09。作者稿；待主线代码复核及独立终审。单位取 c=ℏ=1，度规为 (+,−,−,−)。本件证明的是已衰变分支中的领先概率及加权效应；不以概率核替代子粒子的相干后态。

## A. 同一个物理过程与相位空间效应

母粒子为静止 μ⁺，质量 m，输入自旋态 ρ=(I+s·σ)/2。采用低能标准 V−A 作用，忽略电子和中微子质量及辐射修正。末态记为

$$
\mu^+(p)\longrightarrow e^+(e)+\nu_e(k)+\bar\nu_\mu(l),\qquad
p=(m,\mathbf0),\quad p=e+k+l.
\tag{A1}
$$

所有末态能量为正。三体相位空间为

$$
d\Phi_3=(2\pi)^4\delta^4(p-e-k-l)
 \prod_{a=e,k,l}\frac{d^3\mathbf p_a}{(2\pi)^3\,2E_a}.
\tag{A2}
$$

原始依据是 Bróncano–Mena 的 μ⁻ 极化振幅平方、相位空间及单粒子谱，式(1)—(7)；取 CP 共轭得到这里的 μ⁺ 号。文献与归属见 [sources](sources.md)。在 μ⁺ 静止系，其自旋和后的矩阵元平方为

$$
\sum_{\rm final\ spins}|\mathcal M|^2
 =64G_F^2\,[(p+ms)\cdot k]\,(e\cdot l),\qquad
s^\mu=(0,\mathbf s).
\tag{A3}
$$

混合极化按线性延拓解释，不把 $|s|<1$ 当成纯态四自旋。用同阶总率

$$
\Gamma_0=\frac{G_F^2m^5}{192\pi^3}
\tag{A4}
$$

归一化，得到作用于母自旋的算符值测度

$$
F(d\alpha)=\frac{32G_F^2}{\Gamma_0}
 (e\cdot l)\,[E_k I-\mathbf k\cdot\boldsymbol\sigma]\,d\Phi_3,
\qquad \alpha=(e,k,l).
\tag{A5}
$$

因 $e\cdot l\ge0$、$E_k=|\mathbf k|$，该核为正。积分为 I：迹的归一由(A4)给出，自旋项由完整旋转积分消失。这是条件于发生本衰变的 POVM；未衰变分支没有被解释成概率零。

给任意有限末态分箱 C，令 $F_C=\int_C F$。对任意未知 $\rho_{\mu R}$，包括与任意被动 R 纠缠的输入，其经典箱标签与 R 输出为

$$
\mathcal C(\rho_{\mu R})=
 \sum_C |C\rangle\langle C|\otimes
 \operatorname{Tr}_\mu[(\sqrt{F_C}\otimes I_R)\rho_{\mu R}
 (\sqrt{F_C}\otimes I_R)].
\tag{A6}
$$

完整分箱使其 CPTP。R 分支也可写为 $\operatorname{Tr}_\mu[(F_C\otimes I_R)\rho_{\mu R}]$。这最后一种写法不表示未迹掉的乘积算符本身为正。未分箱连续结果可按算符值测度解释，本文判别只需有限窗口和有界加权矩。

**未领取的后态：** (A3)—(A6)足够确定这些经典记录与 R 边缘，却不确定实际 e⁺ν_eν̄_μ 的非对角相干、检测环境或相干重组仪器。由 F 任造 Stinespring 等距，只是相同概率的一种扩张，不能认作真实衰变后态。真实后态须保完整振幅、正常入射包、准备／时间及探测作用；本轮没有完成此项。

## B. 电子读口与完整补集

令 $x=2E_e/m$，$\mathbf n=\mathbf e/E_e$。对另外两个粒子积分，复用已采用 Michel 效应：

$$
E(x,\mathbf n)=f(x)[A(x)I+B(x)\mathbf n\cdot\boldsymbol\sigma],
\quad f=\frac{x^2}{2\pi},\ A=3-2x,\ B=2x-1.
\tag{B1}
$$

本轮只采用一个二标签电子仪器：

$$
b=\{2/5\le x\le1/2,\ n_z>0\},\qquad \bar b=\text{全部其余电子结果}.
\tag{B2}
$$

边界是零测集，不影响概率。球面积分 $\int_{n_z>0}d\Omega=2\pi$、$\int_{n_z>0}\mathbf n,d\Omega=\pi\hat z$ 给出

$$
E_b=\frac{851}{20000}I-\frac{113}{120000}\sigma_z,
\qquad E_{\bar b}=I-E_b.
\tag{B3}
$$

其两本征值为 $4993/120000$、$5219/120000$，均严格位于(0,1)。对

$$
\rho_\pm=(I\pm\sigma_x)/2
\tag{B4}
$$

两份完整二标签分布完全相同，$p_b=851/20000$。没有只保成功而删掉 $\bar b$。这不声称全部电子角谱相同：全部 Michel 方向读口已信息完备，可以区分 $\rho_+$ 与 $\rho_-$。

## C. 所声明的总源读口

对每个末态事件定义

$$
Q_{ij}(\alpha)=\sum_{a=e,k,l}E_a n_{a,i}n_{a,j}
             =\sum_a\frac{p_{a,i}p_{a,j}}{E_a}.
\tag{C1}
$$

它是**全部三个子粒子**的角能流二阶矩。在自由远场，该量也对应空间积分的、真空扣除的出射自由应力 $T^{ij}$ 的动量对角读口。它不是完整局域应力、应力关联函数或完整几何响应。尤其

$$
\operatorname{tr}_{\rm spatial}Q=m,\qquad |Q_{xz}|\le m/2.
\tag{C2}
$$

第一式来自逐事件能量守恒，第二式用 $|n_xn_z|\le1/2$。同一相位空间还逐事件保 $\mathbf e+\mathbf k+\mathbf l=0$。未读中微子没有从源中删除。

对电子标签定义有界加权效应

$$
W_b=\int_{\alpha:\,e\in b}Q_{xz}(\alpha)F(d\alpha).
\tag{C3}
$$

它是 Hermitian 的源矩，不是一个正 POVM 效应。其操作意义是同一完整联合读数中的 $\mathbb E[\mathbf1_b Q_{xz}]$；可用有限 Q 分箱近似，不能将其当作与电子记录无关的另一份准备。

## D. 两体剩余相位空间与固定方向计算

先固定电子方向 $\mathbf n=\hat z$。令 $q=p-e$，并在 ν 对静止系采用

$$
Q=\sqrt{q^2}=m\sqrt{1-x},\quad
\beta=\frac{x}{2-x},\quad \gamma=(1-\beta^2)^{-1/2},
\quad k_*=(Q/2)(1,\mathbf w),\quad l_*=(Q/2)(1,-\mathbf w).
\tag{D1}
$$

记 $c=w_z$、$w_x=\sqrt{1-c^2}\cos\phi$。回到母静止系的速度为 $-\beta\hat z$：

$$
E_k=\gamma Q(1-\beta c)/2,\quad
E_l=\gamma Q(1+\beta c)/2,\quad
k_z=\gamma Q(c-\beta)/2,\quad
l_z=\gamma Q(-c-\beta)/2,\quad k_x=-l_x=Qw_x/2.
\tag{D2}
$$

于是同一事件的总交叉矩（电子项此时为零）为

$$
Q_{xz}=Q\,\frac{(1-\beta^2)w_xc}{1-\beta^2c^2}.
\tag{D3}
$$

在这一 ν 对球面上，相位空间为常数倍 $dc,d\phi$，自旋权重可写成

$$
\mathcal W=(1+c)\{\gamma m(1-\beta c)
 +\gamma m(\beta-c)s_z-m(w_xs_x+w_ys_y)\}.
\tag{D4}
$$

完整原相位空间常数不必另拟合：无极化球积分为 $4\pi\gamma m(1-\beta/3)$，将它归一到(B1)的 $fA$ 即得所有加权效应。此处理没有改变自旋相对系数。

定义

$$
I(\beta)=\int_{-1}^{1}\frac{c^2(1-c^2)}{1-\beta^2c^2}\,dc.
\tag{D5}
$$

对固定 $x,\hat z$，$Q_{xz}$ 加权效应的 $\sigma_x$ 系数为

$$
\mathsf C(x)=
-\frac{fA\,Q(1-\beta^2)}{4\gamma(1-\beta/3)}I(\beta).
\tag{D6}
$$

这使用 $\int_0^{2\pi}w_x^2d\phi=\pi(1-c^2)$，其余奇方位项消失。对横向输入 $s_z=0$，除以电子密度 $fA$ 后，特别在 $x=1/2$ 有

$$
I(1/3)=150-216\ln2,\qquad
\mathbb E[Q_{xz}\mid x=1/2,\mathbf n=\hat z]
=-m(25-36\ln2)s_x.
\tag{D7}
$$

该式是密度层的校准，非零宽实际事件由下一节给出，不能用(D7)自行领取零概率后选择。

## E. 从固定方向到有限半球

旋转协变及(D4)中仅有标量积的结构，使对称源张量的自旋部分具有形式

$$
\mathsf W_{ij}(x,\mathbf n)\big|_{\rm spin}
=\mathsf C(n_i\sigma_j+n_j\sigma_i)
 +\mathsf D\,\delta_{ij}(\mathbf n\cdot\boldsymbol\sigma)
 +\mathsf E\,n_i n_j(\mathbf n\cdot\boldsymbol\sigma).
\tag{E1}
$$

这里是同一树级相位空间积分的张量分解，不是假设物理拥有新增转向控制。无自旋部分只有 $u\delta_{ij}+v n_i n_j$；在半球的 xz 分量积分为零。由于(A5)、(C1)只涉及普通向量积的分量乘法和点积，方位积分不产生 Levi-Civita 赝张量项。

对 $\mathbf n=\hat z$，令 $\mathsf L$ 为 $\mathsf W_{zz}-\mathsf W_{xx}$ 中的 $\sigma_z$ 系数，则

$$
\mathsf L=2\mathsf C+\mathsf E.
\tag{E2}
$$

方位平均给

$$
\overline{Q_{xx}}=\frac{Q}{2\gamma}
 \frac{1-c^2}{1-\beta^2c^2},\qquad
\overline{Q_{zz}-Q_{xx}}=m-3\overline{Q_{xx}}.
\tag{E3}
$$

这里第二式包括电子的 $Q_{zz}=mx/2$。由(D4)的纵向权重得到

$$
\mathsf L=
\frac{fA}{2(1-\beta/3)}\int_{-1}^{1}
(\beta-c^2)\left[m-\frac{3Q}{2\gamma}
\frac{1-c^2}{1-\beta^2c^2}\right]dc.
\tag{E4}
$$

用半球矩

$$
\int_H n_zd\Omega=\pi,\qquad
\int_H n_x^2n_zd\Omega=\pi/4
\tag{E5}
$$

可知总 $Q_{xz}$ 半球加权效应仅有 $\sigma_x$ 分量，其系数为 $\pi(\mathsf C+\mathsf E/4)=\pi(2\mathsf C+\mathsf L)/4$。代入(D6)、(E4)，再用 $Q/\gamma=m(1-\beta)$，得

$$
\boxed{\kappa(x)=
\frac{m x^2(3-2x)}{16(1-\beta/3)}
\int_{-1}^{1}\frac{N_\beta(c)}{1-\beta^2c^2}\,dc,}
\tag{E6}
$$

$$
\boxed{N_\beta(c)=\frac{\beta(3\beta-1)}2
 +\left(-\frac12+\beta-\frac{\beta^2}2-2\beta^3\right)c^2
 +\left(\beta^3+\frac\beta2-\frac12\right)c^4.}
\tag{E7}
$$

因而

$$
\boxed{W_b=C_b\sigma_x,\qquad C_b=\int_{2/5}^{1/2}\kappa(x)dx.}
\tag{E8}
$$

独立的点值化简为 $\kappa(1/2)/m=(31-45\ln2)/16$。它与(D7)不同：前者已经对电子半球旋转积分，不能混用两个系数。

## F. 有理严格界，不靠小数认定非零

在所选窗内 $1/4\le\beta\le1/3$，(E7)的三个系数分别满足

$$
\frac{\beta(3\beta-1)}2\le0,\qquad
-\frac12+\beta-\frac{\beta^2}2-2\beta^3\le-\frac16,
\qquad
\beta^3+\frac\beta2-\frac12\le-\frac8{27}.
\tag{F1}
$$

分母正且不超过1，故

$$
\int_{-1}^{1}\frac{N_\beta(c)}{1-\beta^2c^2}dc
<-\frac16\int_{-1}^{1}c^2dc=-\frac19.
\tag{F2}
$$

严格性来自(F1)中的负 $c^4$ 项，它在正测度集合上非零。正前因子有 $x^2\ge4/25$、$3-2x\ge2$、$1-\beta/3\le1$，所以至少为 $m/50$。乘以负积分须保正确不等号：

$$
\boxed{\kappa(x)<-\frac m{450},\qquad
C_b<-\frac m{4500}<0.}
\tag{F3}
$$

小数只作校准：$C_b/m\simeq-0.00113809757804$。正文严格判别完全由(F1)—(F3)承担。

## G. 联合判别与电子摘要恢复界

由(B3)、(E8)，$\rho_\pm$ 的电子标签概率相同，但

$$
\mathbb E_{\rho_\pm}[\mathbf1_bQ_{xz}]=\pm C_b,
\qquad
|\Delta\mathbb E[\mathbf1_bQ_{xz}]|
=2|C_b|\ge\frac m{2250}.
\tag{G1}
$$

这是未归一分支矩，未将 $p_b$ 除掉来掩盖稀少事件。若需要条件期望，应除以同一已明列的 $p_b=851/20000$。

记两份完整联合经典输出为 $P_\pm$。因为 $\mathbf1_bQ_{xz}\in[-m/2,m/2]$，

$$
\operatorname{TV}(P_+,P_-)\ge
\frac{|\Delta\mathbb E[\mathbf1_bQ_{xz}]|}{m}
\ge\frac1{2250}.
\tag{G2}
$$

有限 Q 分箱的重构误差若每事件至多 $\delta_Q$，两输入矩差的损失至多 $2\delta_Q$。若采用这条分箱路线，恢复下界相应变为 $\epsilon\ge |C_b|/m-\delta_Q/m>1/4500-\delta_Q/m$，不能对任意分箱仍保原常数。例如 $\delta_Q<m/9000$ 给严格大于 $1/9000$ 的恢复下界。因此证书不依赖无限读数分辨率。下面另给不经此误差估计的精确三结果任务。

还可将证书直接写成三个有限结果，不需要连续源读数的无限分辨率。在 b 分支上采用有界源效果 $f_\pm(Q)=(1\pm2Q/m)/2$，定义

$$
G_{b,+}=E_b/2+W_b/m,\qquad
G_{b,-}=E_b/2-W_b/m,\qquad
G_{\bar b}=I-E_b.
\tag{G3a}
$$

因为 $0\le f_\pm(Q_{xz})\le1$，前两式是 $\int_b f_\pm(Q_{xz})F(d\alpha)$，均为正；三项求和为I，且合并前两项恰好恢复原电子 b 结果。对 $\rho_\pm$，这份三结果分布的总变差**恰为** $2|C_b|/m$。对任意未知 R，仍用(A6)即可给出其完整经典结果／R边缘。

这是有界源任务的有限数学效果。把它作为先测能流再作有偏抽签的物理实现，还需声明源读取、随机资源与实际装置；本文不赠送这些权限，不声称三结果效果唯一确定材料测后态或引力检测器。

若仅从电子二标签施加一个与未知自旋无关的后处理，再尝试恢复上述三结果任务，两个输入的预测仍相同。对二者都承诺总变差误差不超过 ε 时，三角不等式强制

$$
\boxed{\epsilon\ge |C_b|/m>1/4500.}
\tag{G3}
$$

相应地，源矩本身的最坏绝对预测误差至少 $|C_b|\ge m/4500$。这是对(B2)这份摘要的限定否决，不否定保完整电子角谱、实际子粒子状态或额外关系变量的字典。任意 R 合同包含平凡 R，因此(G3)也阻止一个声称对全部未知输入及 R 都达更小误差的恢复；它不要求复制未知母态。

## H. 为何必须保联合标签：无条件抵消

由同一(A3)得到 μ⁺ 单粒子能流（积分各自能量）

$$
\frac{d\langle\mathcal E_e\rangle}{d\Omega}
=\frac{d\langle\mathcal E_{\bar\nu_\mu}\rangle}{d\Omega}
=\frac m{80\pi}(7+3\mathbf s\cdot\mathbf n),\qquad
\frac{d\langle\mathcal E_{\nu_e}\rangle}{d\Omega}
=\frac m{80\pi}(6-6\mathbf s\cdot\mathbf n).
\tag{H1}
$$

三者之和为 $m/(4\pi)$，自旋项严格消失，故无条件 $\langle Q_{xz}\rangle=0$。补集满足 $W_{\bar b}=-W_b$。本轮信号是实际电子标签与**总**能流之间的相关，不是没有收到标签也能看到的总平均差；更不是通过丢掉电子能量制造的中微子子部门信号。

## I. 有限衰变机会与停止线

可以另声明一个自旋无关的衰变机会参数 $p\in[0,1]$：衰变分支效应乘 p，另外保 $(1-p)I$ 的未衰变标签。对只给已衰变产物计分的联合任务，(G1)—(G3)相应乘 p。这个扩展仅为声明模型，不证明真实有限窗恰满足指数律，也不把未衰变母粒子的能源置零。

本轮的源是领先、自由远场、有界角能流矩。一般局域 $T^{\mu\nu}(x)$、相干应力、有限频率响应和噪声需要完整振幅与正常波包；近衰变顶点还须共同变分交互作用。线性／非线性GR响应、引力memory仪器、QED软辐射及IR完成、探测装置总源和后续几何反作用均未认证。以上限制不撤销有限窗口的领先联合概率证书，也不把完成全部装置当成本命题的前提。
