# 1258证明：全占据来源预算到有限耦合记录误差

2026-10-16。接[1257](../1257/proof.md)。CT仍为待验证猜想，新增公理0。本轮在原自由质量1的$X$场、固定FLRW保持背景和事先选定的共同参考内，证明一整族有限正激发能源工作态对局部有限探针记录的影响有统一上界。并不要求所有工作态相同，也不重新选择或证明TR共同热参考角色。

## 1. 命题、既有接口与本轮范围

认知动机是：共同参考已确定之后，参与系统仍可有不同的非平衡内容；需要判定这些差异何时小于一个有限任务的误差。这个动机不替代场、背景、准备和端口的模型输入。

直接复用：

- [正式TR采用记录](../../archive_1201_1218/_shared/thermal_roles_adoption_after1218.md)与[1210 §§2—6](../../archive_1201_1218/1210/proof.md)：共同参考角色、Gibbs分类、共同$\beta$及温标已经建立。不将工作态是否唯一当成这些成果的新前提。
- [730](../../archive_702_741/research_note_730.md)、[810](../../archive_797_818/research_note_810.md)、[842](../../archive_819_853/research_note_842.md)：合法参考、同处方来源及正常热参考预算；这些旧结果没有给本轮全占据族的局部尺度误差。
- [1004](../../archive_990_1008/research_note_1004.md)与[1257](../1257/proof.md)：同背景的有限记录可以依赖状态；1257保全平均来源的壳见证继续有效。
- [1247 §7.1](../1247/proof.md)、[1248 §2](../1248/proof.md)：原保持背景存在且有明确正域。此处只用其保守界，不重跑背景ODE。
- [1232 §§3—6](../../archive_1219_1233/1232/proof.md)：有界耦合的有限全概率比较及实际钟、运输、来源合同。其Dirac谱非空证明和数值不能移作本轮$X$场的谱证书。

本轮净连接是

$$
\text{同一正激发能源预算}
\ \Longrightarrow\ 
\text{全占据族的局部协方差统一界}
\ \Longrightarrow\
\text{同一有限耦合后任意探针记录的统一误差}.
\tag{1}
$$

有限正激发预算不是实际制备功、控制能源或绝对真空来源。以下概率界在同一固定背景上比较，不能自动吸收态改变造成的几何反馈。

## 2. 原场、正常占据态与模式能源

取1257的同一环面、场和模式归一：

$$
I=(2,9/4),\qquad \mathcal V=10^3,\qquad
ds^2=-dt^2+a^2(t)d\mathbf x^2,
$$

$$
\widehat X(t,\mathbf x)=X_{\mathrm{cl}}(t)I+
\mathcal V^{-1/2}\sum_{\mathbf n}
\left(a_{\mathbf n}u_{\mathbf n}(t)e^{i\mathbf k_{\mathbf n}\cdot\mathbf x}
+a_{\mathbf n}^{\dagger}u_{\mathbf n}^*(t)e^{-i\mathbf k_{\mathbf n}\cdot\mathbf x}\right),
\qquad \mathbf k_{\mathbf n}=\frac{2\pi}{10}\mathbf n.
\tag{2}
$$

共同纯准自由Hadamard参考$\omega_0$和相同经典均值$X_{\mathrm{cl}}$事前固定；$\mathbf n$与$-\mathbf n$仍是两个行波振子，实性由Hermitian共轭保证。保持期内

$$
\ddot u_{\mathbf n}+3H\dot u_{\mathbf n}+\omega_{\mathbf n}^2u_{\mathbf n}=0,
\quad
\omega_{\mathbf n}^2=1+\frac{|\mathbf k_{\mathbf n}|^2}{a^2},
\quad
a^3(u_{\mathbf n}\dot u_{\mathbf n}^*-u_{\mathbf n}^*\dot u_{\mathbf n})=i.
\tag{3}
$$

令$t_0\in I$固定。工作态$\omega_{\mathbf n_{\mathrm{occ}}}$是在这些模式上的居中、无异常配对准自由占据态，并加与参考相同的经典均值。为避免整数动量与占据混淆，以下用$N_{\mathbf n}\ge0$表示占据均值。每个模式的密度为1257式(5)的几何分布，规定

$$
e_{\mathbf n}(t)=a^3(t)
\left(|\dot u_{\mathbf n}(t)|^2+\omega_{\mathbf n}^2(t)|u_{\mathbf n}(t)|^2\right),
\qquad
\sum_{\mathbf n}N_{\mathbf n}e_{\mathbf n}(t_0)\le E_*,\qquad 0<E_*<\infty.
\tag{4}
$$

这是相对共同模式真空的正平均场能源增量。Wronskian和算术—几何均值不等式给

$$
e_{\mathbf n}(t)\ge
2a^3\omega_{\mathbf n}|u_{\mathbf n}||\dot u_{\mathbf n}|
\ge\omega_{\mathbf n}(t)\ge1.
\tag{5}
$$

因此$\sum N_{\mathbf n}<\infty$。全模式乘积占据态在同一Fock表示中正常：其总粒子数期望有限，且$\prod_{\mathbf n}(1+N_{\mathbf n})^{-1}>0$；有限模式密度的极限由有限总占据控制，不是无穷热浴的形式乘积。无需一个共同最大波数。

**有限能源不自动推出Hadamard。** 若要保留同一Hadamard减除和点态来源，另限定占据差具有相应光滑性，例如有限支持或快速衰减占据；本轮使用该正常Hadamard子族。下面的误差常数不依赖其截止、衰减阶或具体占据。两点与概率估计本身只需式(4)及所列准自由结构，不能据此宣称任意有限能源态均已具有Hadamard资格。也不把该族扩成有任意挤压、异常配对或非Gaussian高阶关联的所有态。

从(3)直接求导：

$$
\dot e_{\mathbf n}
=Ha^3\left[-3|\dot u_{\mathbf n}|^2+
\left(3+\frac{|\mathbf k_{\mathbf n}|^2}{a^2}\right)|u_{\mathbf n}|^2\right],
\qquad
|\dot e_{\mathbf n}|\le3|H|e_{\mathbf n}.
\tag{6}
$$

对包含所有测试支持和$t_0$的固定紧子窗$J\subset I$，设$a\ge a_->0$，并取

$$
\mathcal B=\exp\left(3\sup_{t\in J}
\left|\int_{t_0}^{t}|H(s)|\,ds\right|\right).
\tag{7}
$$

于是所有模式、所有测试时刻共同满足$e_{\mathbf n}(t)\le\mathcal B e_{\mathbf n}(t_0)$。这保留背景的真实时间变化，没有假设瞬时场能源守恒。

## 3. 局部涂抹与全模式统一界

取固定实$\chi,h\in C_c^\infty$，空间支持位于一个有限坐标块，并满足

$$
\int h(\mathbf y)d^3y=0,\qquad
\int y_i h(\mathbf y)d^3y=0.
\tag{8}
$$

令$M_\chi=\|\chi\|_1$。选有限常数$C_0,C_2$使

$$
|\widehat h(\mathbf q)|\le
\min(C_0,C_2|\mathbf q|^2).
\tag{9}
$$

总可取$C_0=\|h\|_1$和$C_2=\frac12\int|\mathbf y|^2|h(\mathbf y)|d^3y$；后面的明确函数用Fourier恒等式得到更小的$C_2$，不声称其绝对二阶矩等于该常数。

对测试中心$(\tau,\mathbf x_\tau)$及固定正尺度$\ell$定义

$$
F_{\tau,\ell}(t,\mathbf x)
=a^{-3}(t)\ell^{-3}
\chi\left(\frac{t-\tau}{\ell}\right)
h\left(\frac{\mathbf x-\mathbf x_\tau}{\ell}\right).
\tag{10}
$$

所有测试支持须留在$J$和不绕行的空间图内。体积测度$a^3dt\,d^3x$与$a^{-3}$相消，时间换元仍留下$\ell$。所以

$$
\Phi(F_{\tau,\ell})=\sum_{\mathbf n}
\left(z_{\mathbf n,\tau,\ell}a_{\mathbf n}
+z_{\mathbf n,\tau,\ell}^*a_{\mathbf n}^{\dagger}\right),
$$

$$
z_{\mathbf n,\tau,\ell}
=\frac{\ell}{\sqrt{\mathcal V}}\,
e^{i\mathbf k_{\mathbf n}\cdot\mathbf x_\tau}
\widehat h(-\ell\mathbf k_{\mathbf n})
\int \chi(s)u_{\mathbf n}(\tau+\ell s)\,ds.
\tag{11}
$$

齐次经典均值被(8)消去。各测试因而有相同的零均值；任意其他相同均值只会在下面Weyl期望中增加一个共同相位。

对非零$\mathbf k$，利用

$$
a^3\omega_{\mathbf n}^2=a^3+a|\mathbf k_{\mathbf n}|^2
\ge a_-|\mathbf k_{\mathbf n}|^2
$$

及

$$
\frac{\min(C_0^2,C_2^2\ell^4|\mathbf k|^4)}{|\mathbf k|^2}
\le C_0C_2\ell^2
$$

得到

$$
|z_{\mathbf n,\tau,\ell}|^2
\le
\frac{\mathcal B M_\chi^2 C_0C_2}
{\mathcal V a_-}\,
e_{\mathbf n}(t_0)\ell^4.
\tag{12}
$$

零模须单独处理：$\widehat h(0)=0$，故其$z$严格为零；没有在$\mathbf k=0$使用除以$|\mathbf k|^2$的式子。

相对参考的额外方差为非负量

$$
d_{\tau,\ell}
=2\sum_{\mathbf n}N_{\mathbf n}|z_{\mathbf n,\tau,\ell}|^2
\le
D_\ell:=
\frac{2\mathcal B M_\chi^2 C_0C_2}
{\mathcal V a_-}\,E_*\ell^4.
\tag{13}
$$

这是同一能源球上、任意允许测试中心及所有模式的统一界。高波数不会因没有统一截止而逃出该估计。

### 3.1 多个测试的共同协方差

额外对称协方差是正半定二次型：

$$
d(F)=2\sum_{\mathbf n}N_{\mathbf n}|z_{\mathbf n}(F)|^2\ge0.
$$

所以对任意$r$个允许测试和$\sigma_j=\pm1$，

$$
d\left(\sum_{j=1}^r\sigma_jF_j\right)
\le r\sum_{j=1}^r d(F_j)\le r^2D_\ell.
\tag{14}
$$

测试可以时空重叠，彼此不必对易。此处的正半定性属于占据增量，不能对任意两态的有符号协方差差直接使用。

对固定有限壳$|\mathbf k|\le K_*$，也可用$a^3\omega^2\ge a_-^3$和(9)得到$d_{\tau,\ell}\le2\mathcal B M_\chi^2C_2^2K_*^4E_*\ell^6/(\mathcal V a_-^3)$。这只是附属强化；本轮主结果是没有固定壳要求的(13)。

## 4. 任意有序有界场乘积的期望差

记$W(F)=e^{i\Phi(F)}$。同一场代数保证两态拥有相同CCR，故

$$
W(F_1)\cdots W(F_r)
=\exp\left[-\frac i2\sum_{j<k}\sigma(F_j,F_k)\right]
W\left(\sum_jF_j\right).
\tag{15}
$$

相位只依赖共同交换子。准自由特征函数和相同均值给

$$
\left|
\omega_{\mathbf n_{\mathrm{occ}}}\!\left(\prod_{j=1}^rW(\sigma_jF_j)\right)
-\omega_0\!\left(\prod_{j=1}^rW(\sigma_jF_j)\right)
\right|
\le1-e^{-d(\sum_j\sigma_jF_j)/2}
\le\frac{r^2D_\ell}{2}.
\tag{16}
$$

令$B_\tau=\cos\Phi(F_{\tau,\ell})=(W(F_{\tau,\ell})+W(-F_{\tau,\ell}))/2$，则$B_\tau$自伴且$\|B_\tau\|\le1$。展开$r$个cos为$2^r$个Weyl积，并保留其原次序，得到

$$
\left|\omega_{\mathbf n_{\mathrm{occ}}}(B_{\tau_1}\cdots B_{\tau_r})
-\omega_0(B_{\tau_1}\cdots B_{\tau_r})\right|
\le\frac{r^2D_\ell}{2}.
\tag{17}
$$

$r=0$时两边之差严格为零。这一步没有将不同时刻的场替换为可交换随机变量。

## 5. 有限耦合、旁挂参考与实际记录误差

设探针$P$有限维，旁挂参考$R$可任意；初态$\sigma_{PR}$任意，允许未知纠缠。**整份$PR$与场初态独立**，分别比较

$$
\sigma_{PR}\otimes\omega_{\mathbf n_{\mathrm{occ}}},
\qquad \sigma_{PR}\otimes\omega_0.
\tag{18}
$$

两过程使用同一自由探针演化及同一相互作用图景耦合

$$
V(\tau)=g(\tau)M(\tau)\otimes B_\tau,
\quad M(\tau)=M(\tau)^\dagger,\quad \|M(\tau)\|\le1,
\quad g(\tau)\in\mathbb R,\quad
\mathcal L=\int|g(\tau)|d\tau<\infty.
\tag{19}
$$

$M$和$B$取有界强连续族，$g$取实可积有限窗函数。由此强积分定义的Dyson级数绝对受控，给同一个有界自伴相互作用的酉过程$U$；不对无界裸$\Phi(F)$直接作未经控制的Dyson展开。

任取探针和参考的末效果$0\le Q_{PR}\le I$。展开$U^\dagger(Q_{PR}\otimes I)U$，左支$p$个、右支$q$个耦合包含一个$r=p+q$项场有序积。探针系数的期望绝对值不超过1；左右有序积分的绝对权不超过$\mathcal L^p/p!$及$\mathcal L^q/q!$。由(17)和

$$
\sum_{p+q=r}\frac{\mathcal L^p\mathcal L^q}{p!q!}
=\frac{(2\mathcal L)^r}{r!},
\quad
\sum_{r=0}^\infty\frac{r^2x^r}{r!}=x(1+x)e^x
$$

得到全阶有限概率界

$$
\boxed{
|P_{\mathbf n_{\mathrm{occ}}}(Q_{PR})-P_0(Q_{PR})|
\le
\min\left\{1,\,
D_\ell\mathcal L(1+2\mathcal L)e^{2\mathcal L}\right\}.
}
\tag{20}
$$

没有丢掉有限耦合的高阶项；这不是只对领先二阶热响应的估计。界不依赖$\sigma_{PR}$、$R$的维数、所选末效果或具体占据分布。因为初场独立，两过程定义探针通道$\Lambda_{\mathbf n_{\mathrm{occ}}},\Lambda_0$，同一证明也给

$$
\frac12\|\Lambda_{\mathbf n_{\mathrm{occ}}}-\Lambda_0\|_\diamond
\le
\min\left\{1,\,
D_\ell\mathcal L(1+2\mathcal L)e^{2\mathcal L}\right\}.
\tag{21}
$$

这里半diamond是允许旁挂参考的最坏输出可辨概率差；(20)已足够作为有限任务结论。若初探针与场相关，便没有本稿所用的共同独立环境通道，不能无条件使用(20)—(21)。

时空涂抹$B_\tau$是同一场代数中的有界比较元；将它写成(19)不自动赋予瞬时实施一个时空区域操作的权限。实际探针若要继承(20)，仍须按1232 §6认证区域运输、固有钟、因果实现和与该比较过程的误差。这是允许有效端口的合同，不是已经制造控制器、原$K$读口或微观局部仪器。

## 6. 原保持窗内的一份明确有限证书

定义一维光滑bump

$$
b_1(s)=
\begin{cases}
C\exp[-1/(1-s^2)],&|s|<1,\\
0,&|s|\ge1,
\end{cases}
\qquad \int b_1=1.
\tag{22}
$$

在$|s|\le1/2$有指数项至少$e^{-4/3}>1/9$，因此$0<C<9$。对$v=(1-s^2)^{-1}\ge1$，未归一指数函数$\varphi$满足

$$
\varphi''=e^{-v}(-2v^2-8s^2v^3+4s^2v^4).
$$

用$v^ke^{-v}\le(k/e)^k$及$e>2$，得到$|\varphi''|<2+27+64=93$，故

$$
\|b_1''\|_1<2\cdot93\cdot9=1674.
\tag{23}
$$

取

$$
b(\mathbf y)=b_1(y_1)b_1(y_2)b_1(y_3),\qquad
h=(\partial_1^2-\partial_2^2)b,\qquad \chi=b_1.
\tag{24}
$$

支持是单位坐标立方体，不称为单位球。积分分部保证零阶及一阶矩为零，并且

$$
M_\chi=1,\quad C_0<3348,\quad
\widehat h(\mathbf q)=-(q_1^2-q_2^2)\widehat b(\mathbf q),\quad
|\widehat b|\le1.
$$

因此(9)可直接采用$C_2=1$。这是Fourier估计，不是绝对矩积分的等式。

原背景取$t_0=17/8$，$a\ge1$、$|H|<11/10$。所有测试时刻留在$I$时，到$t_0$的距离小于$1/8$，所以

$$
\mathcal B\le e^{33/80}<\frac{80}{47}<2.
\tag{25}
$$

指数护栏来自$e^x<1/(1-x)$，$0<x<1$。取$0<\ell\le1/16$，耦合中心与空间中心事前确定；例如

$$
g(\tau)=\frac18\,32b_1\!\left(32(\tau-17/8)\right)
\tag{26}
$$

给$\mathcal L=1/8$。它的中心支持半宽$1/32$，加测试半宽$\ell\le1/16$仍严格位于$I$；空间支持也不绕行原周期10环面。这里(26)只是有界比较过程的明确包络，不是免费实际开关。

将$\mathcal V=1000$代入(13)，得到

$$
D_\ell<14E_*\ell^4.
\tag{27}
$$

另外

$$
\mathcal L(1+2\mathcal L)e^{2\mathcal L}
=\frac5{32}e^{1/4}<\frac5{24},
$$

从而

$$
\boxed{
|P_{\mathbf n_{\mathrm{occ}}}(Q_{PR})-P_0(Q_{PR})|
<3E_*\ell^4.
}
\tag{28}
$$

例如$E_*=1,\ell=1/16$时，

$$
|P_{\mathbf n_{\mathrm{occ}}}-P_0|
<\frac3{65536}<\frac1{10000}.
\tag{29}
$$

这是一个正的、有限的探测尺度及完整耦合强度证书；没有要求微观必须连续到任意小尺度，也没有借无限采样制造精确结论。任意两个同预算工作态还可通过共同参考用三角不等式得到两倍(20)的保守差界。固定参考与固定菜单没有随各工作态结果重调。

## 7. 同一来源预算及仍需另外履行的责任

占据态齐次，且与参考有相同经典均值。保持同一减除和局部反项，新增平均场能源密度恰为

$$
\Delta\rho_X(t)
=\frac1{\mathcal V}\sum_{\mathbf n}N_{\mathbf n}
\left(|\dot u_{\mathbf n}|^2+\omega_{\mathbf n}^2|u_{\mathbf n}|^2\right)
=\frac{\sum N_{\mathbf n}e_{\mathbf n}(t)}
{\mathcal V a^3(t)}
\le\frac{\mathcal B E_*}{\mathcal V a_-^3}.
\tag{30}
$$

同一正占据协方差还给正交标架中各通量、压力分量的绝对值不超过该能源密度上界。无需把各向异性占据删掉；若$N_{\mathbf n}=N_{-\mathbf n}$，平均通量才进一步为零。

对任一类光向量$k^a$，用$|\mathbf k^{\mathrm{orth}}|=|k^t|$和$|A+B|^2\le2|A|^2+2|B|^2$，

$$
0\le\Delta T^X_{ab}k^ak^b
\le4\Delta\rho_X(k^t)^2
\le\frac{4\mathcal B E_*}{\mathcal V a_-^3}(k^t)^2.
\tag{31}
$$

因此明确证书给$\Delta\rho_X<E_*/500$。这不是把来源置零：基准真空、共同经典均值、其余物质及新增占据都必须保留。若要把这个预算换成几何或读数的反作用误差，还须采用已认证的同过程响应界；本轮没有证明原纯经典$g$在加入量子源后仍精确满足半经典方程。

(30)—(31)只给**准备占据增量**的来源界。耦合读出还可能改变场、探针、控制、钟和记录；其全过程来源必须按已有完整资源合同登记。缩小$\ell$可能提高实际空间选择、运输或开关成本，不能由(28)反推这些成本也趋零。

## 8. 对已有热比较合同的用法与停止线

本轮可以为一个已选参考的有限比较合同补上状态误差分量。若该参考在指定区域、时间流和菜单下已有温标预测，且参考比较误差为$\epsilon_{\mathrm{ref}}$、实际端口与比较过程的匹配误差为$\epsilon_{\mathrm{impl}}$，则同一概率任务可使用

$$
\epsilon_{\mathrm{total}}
\le\epsilon_{\mathrm{ref}}+
D_\ell\mathcal L(1+2\mathcal L)e^{2\mathcal L}
+\epsilon_{\mathrm{impl}},
\tag{32}
$$

另计本任务需要的来源反作用误差。它允许有限能源的非平衡工作态与共同参考不同，不需要先证明所有工作态都是TR、Gibbs或KMS。

但(32)不会产生一个尚未给出的参考温度、正热响应、实际加速流或源到几何响应。若目标改成概率比或对数热比，还须保留原1232所需的分母正下界，再把绝对概率误差传递过去；不能因(29)很小便省去稀有事件分母。当前$X$的参考不自动等于Dirac楔参考，其热合同不能按名称移植。

解析证明承担无限模式量词、正常态边界、能源输送、正协方差、全阶Dyson求和及旁挂参考。配套有限核验可检查bump的有理护栏、明确常数、协方差/高频一致性与有限探针比较；有限矩阵不替代以上证明。停止于这项“来源预算—有限耦合记录”的新连接，不扩大为完整GR、Unruh仪器、任意量子态或热参考唯一性的证明。
