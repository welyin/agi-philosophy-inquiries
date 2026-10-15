# 1229：实际场后态的信息、辅助熵与独立楔区 boost 预算

## 0. 本轮只签什么

本轮接[1228 NEXT](../1228/NEXT.md)及[准入](../_admission/after1228_geometric_information_admission.md)，在1126/1228**同一两味有质量自由 Dirac 真空**和已登记的有限区域模式门许可内，构造一个实际非纯酉场准备，再用同一个固定场读口比较它与原真空。

核心结论为不等式，而非预设的熵等式：

$$
0<D_W(\rho_F'\Vert\omega)+I(W:E)_{\rho'}
\le \Delta K_W+S(\rho_E')<\infty. \tag{1}
$$

所有对象属于同一实际有限过程；没有先假设目标相对熵有限。固定后继场读口的信息严格大于 $11$ nat，真实实体读口与它的有限 KL 比较误差小于 $10^{-3}$。独立几何合同的有限误差可再小于 $10^{-2}$。

净增不是再次证明一般链式分解，而是：**完整 Gaussian 门尾的严格楔定位、真实无限维辅助的忠实参考与模域、实际同场后继记录、同源 boost 与独立几何误差，在本菜单里同时满足。**

限制仍保留：本轮区域双局域 CAR 核是1126的有效许可，不签逐点微观因果装置；钟不是加速固有时钟；数学参考 $\tau_E$ 不是实际热浴。几何坐标/校准由独立输入合同提供，没有造出新的几何仪器。所得是信息—几何加权能量预算，**不是面积律、面积熵系数、Unruh测温或Einstein方程**。本轮不证明(1)取等，不需要优化1228的极端钟参数。

## 1. 固定材料、两个门和两个时间截面

复用1228的 $m=1,K=10^6$、实紧支包络、旋量与两味模式：

$$
f_\pm=e^{\pm iKz}\phi\chi,\quad
\|f_\pm\|=1,\quad
f=(f_+,f_-)/\sqrt2,\quad
B^2=\frac{137}{5}.
$$

$f_\pm$ 各自为单味单位模式，完整材料保全部正、负谱，$H_F\ge0$、真空 $\Omega$。包络为 $C^7\cap H^8$，支撑 $x,y\in[-1,1],z\in[1,3]$。第二读口的模式 $u=(f_+,0)$ 嵌入味1。令

$$
G_v=n_v\otimes K_{PE},\quad n_v=c^\dagger(v)c(v),\quad
K_{PE}=\frac\pi2(I_\ell-X_\ell),\quad 0\le G_v\le\pi.
$$

两对完整四维 $P_jE_j$ 保留，空白为 $|00\rangle$。两钟规格不再优化：

$$
\delta=10^{-24},\qquad \sigma=10^{-120},\qquad M=10^{260},
\quad H_{C_j}=P_j^2/(2M)+1.
$$

两钟初态为正常纯 Gaussian，$\langle X_j\rangle=0,\langle P_j\rangle=M$，方差 $\Delta X=\sigma,\Delta P=1/(2\sigma)$。令

$$
g_j(X)=\delta^{-1}b_c((X-t_j)/\delta),\quad
b_c(x)=\frac{15}{16}(1-x^2)^2_+,\quad t_1=1,\ t_2=3,
$$

$$
k_R=e^{-2ih}u,\qquad
H_{\rm all}=H_F+H_{C_1}+H_{C_2}+4
             +g_1(X_1)G_f+g_2(X_2)G_{k_R}. \tag{2}
$$

全窗为 $0\le t\le4$，**准备截面为 $T_p=2$**。初辅助 $E=P_1E_1P_2E_2C_1C_2$ 全纯且与场独立。第二钟在 $[0,2]$ 的前尾真实保留，第一钟在 $[2,4]$ 的后尾也保留。没有重置、关门或丢辅助。

以第一门中心 $t=1$ 为场自由相互作用绘景基准，两门场模式在全窗分别是

$$
v_1(t)=e^{i(t-1)h}f,\qquad
v_2(t)=e^{i(t-3)h}u.
$$

故理想中心仪器是先记录 $n_f$，后记录**固定场读口** $n_u$。精确地，若 $U_A(T,0)$ 是(2)的 Schrödinger 演化，$H_E^0=H_{C_1}+H_{C_2}+4$，本稿使用
$$
U_I(T)=e^{i(T-1)H_F}e^{iT H_E^0}\,U_A(T,0)e^{iH_F}.
$$
$U_I(0)=1$，且初真空被最后的 $e^{iH_F}$ 保持。定义 $U_p=U_I(2)$，其实际场约化态为 $\rho_F'$。辅助同时转到自由绘景不改变其熵。后文场态、boost与楔定位均采用此共同绘景，不能和未经运输的 Schrödinger 态混用。

## 2. 全部时间尾下的严格楔定位与 boost 符号

事前固定

$$
W_8=\{(t,x,y,z):z+8>|t-1|\}.
$$

物理 Schrödinger 模式 $f$ 的支撑纵向在 $[1,3]$；$k_R=e^{-2ih}u$ 的有限传播支撑在 $[-1,5]$，横向在 $[-3,3]^2$。所以整个 $t\in[0,4]$，包括所有 Gaussian 钟尾，均有

$$
z+8-|t-1|\ge4.
$$

这是固定有限时间窗的严格定位，不声称无限未来的同一楔定位。等价地，在首门切面的 Cauchy 代数里，全部 $v_j(t)$ 支撑纵向不低于 $-2$，距离边界 $z=-8$ 至少6。

每个场门为偶的局部 CAR 多项式；辅助系数虽作用在完整无界钟上，其 $g_j$ 是有界算符。因此 Dyson 强积分和强闭性给

$$
U_p,\ U_I(4)\in\mathcal A(W_8)\ \overline\otimes\ \mathcal B(E). \tag{3}
$$

有限区域核许可与楔定位是不同层次：(3) 不把有效核自动升级成逐点可制造装置。实际辅助全被包含，$W_8$ 外场代数的状态仍为原真空。

首门切面的一体 boost 为

$$
b_8=k_z+8h,\quad k_z=\tfrac12(zh+hz),\qquad
[h,k_z]=-ip_z.
$$

由此 $e^{ith}k_ze^{-ith}=k_z+tp_z$。对应实际终截面 $T$ 的 Schrödinger 观测量是

$$
B_{\rm Sch}(T)=B+8H_F-(T-1)P_z. \tag{4}
$$

这里 $B$ 表示未平移边界的空间 boost。在全局初始切面 $t=0$，同一几何荷为 $B+8H_F+P_z$；在首门切面为 $B+8H_F$。后文简写 $\mathsf B$ 为共同首门切面的第二量子化 $b_8$，$\mathsf B\Omega=0$。不得把不同切面的 $B$ 直接互换。

BW 前件由原自由 Dirac 真空满足，故楔参考模生成元差为

$$
\Delta K_W=2\pi\langle\mathsf B\rangle_{\rho_F'}. \tag{5}
$$

沿1228给 $C=16K$ 控制单门 $[H_F,G]$。对全窗，$|t-1|,|t-3|\le3$；由支撑 $|z+8|\le11$ 得

$$
\|b_8 v_j(t)\|\le14(K+B)+\tfrac12<15K.
$$

因此单门 boost 交换子有共同界

$$
\|[\mathsf B,G_{v_j(t)}]\|\le C_B,\qquad C_B=128K. \tag{6}
$$

这是原场的矩界，未用巨大钟总 Hamiltonian 代替它。对 $U_p$，粗界
$\|\mathsf B U_p\Psi_0\|\le2T_p C_B/\delta$ 有限。实际到经典脉冲再到理想记录的两步运输还给更强的期望值界，见第6节。

## 3. 忠实辅助参考：不把初纯态当成热态

在辅助自由绘景定义初始 Gaussian 的振子模

$$
a_j=X_j/(2\sigma)+i\sigma(P_j-M),\qquad
N=N_1+N_2,\quad N_j=a_j^\dagger a_j.
$$

实际初钟是 $N=0$ 的纯态。仅作相对熵证明用的忠实正常参考为

$$
\tau_E=\frac{e^{-N}}{Z}\otimes\frac{I_{PE}}{16},
\qquad Z=(1-e^{-1})^{-2}.
$$

它没有被称为实际热浴或物理 boost 能源。初实际辅助相对于它的交叉熵为 $\log Z+\log16$，有限。后态必须另核 $N$ 矩；不能先宣称所有相对熵有限。

1228的力界沿全窗可取 $L=16/\delta^2$。实际钟在自由绘景的平均占有满足

$$
\langle N_j\rangle
\le \bar n(T)=\left(\sigma LT+\frac{LT^2}{4M\sigma}\right)^2,\qquad T\le4. \tag{7}
$$

这是原真实无限秩后态的界，不假定 Gaussian 或末态独立。因此

$$
S(\rho_E')\le4\log2+2g(\bar n(4))<2.8,
\quad
\operatorname{Tr}\rho_E'(-\log\tau_E)
=\langle N\rangle+\log Z+\log16<\infty. \tag{8}
$$

初始完整辅助熵为0。若换为混合初态，后续预算必须保初熵差，本轮未默默延拓到该情形。

## 4. 有限辅助逼近及目标相对熵的非循环证明

### 4.1 只为证明的辅助截断

取 $P_L=1_{\{N\le L\}}\otimes I_{PE}$，把自由绘景交互替成 $V_L(t)=P_LV(t)P_L$，在补空间取零交互，得到 $U_{p,L}$。这是**辅助证明序列**，不是给实际场或钟增加物理截断。

$V_L(t)$ 有限辅助维数，并在场 boost 导子域；(6)控制其导子。有限辅助版本的 KMS 熵平衡可直接复用 [Jakšić–Pillet，定理3.1、式(3.9)—(3.13)](https://arxiv.org/html/1406.0034)，应用到本轮楔 KMS 参考及有限辅助。初辅助可以纯，参考保持忠实；嵌回固定 $\tau_E$ 时要补回辅助参考归一因子。得到

$$
D(\rho'_{WE,L}\Vert\omega_W\otimes\tau_E)
=2\pi\langle\mathsf B\rangle_L+\langle N\rangle_L
 +\log Z+\log16. \tag{9}
$$

具体消去步骤是：在有限辅助空间上先取论文的正比较泛函 $\eta_L=\omega_W\otimes I_{E,L}$（未归一），初态为 $\omega_W\otimes|0,\mathrm{blank}\rangle\langle0,\mathrm{blank}|$。初相对熵 $D(\rho_0\Vert\eta_L)=-S(\rho_{E,0})=0$；场导子 $\delta_R=i[\mathsf B,\cdot]$ 的有界域已由(6)核定，故定理3.1给
$$
D(\rho'_{WE,L}\Vert\eta_L)
=-i(2\pi)\rho_0(U_{p,L}^\dagger\delta_R(U_{p,L}))
=2\pi\langle\mathsf B\rangle_L.
$$
然后用式(3.10)把辅助比较换成 $\tau_{E,L}=e^{-N}P_L/(Z_L\,16)$，加上交叉熵 $\langle N\rangle_L+\log Z_L+\log16$；嵌回固定的完整 $\tau_E$ 又加 $\log(Z/Z_L)$，恰得(9)。初态无须等于忠实参考，纯辅助不是隐藏的满秩假设。

该成熟熵平衡不计作本轮新定理；这里需要的是本菜单确实进入其域。未归一 $\eta_L$ 只在上面有限级证明中使用，最终所有极限比较均取同一个正常、归一且忠实的 $\tau_E$。

### 4.2 截断一致的辅助二阶矩

令 $d_T=\sigma+T/(2M\sigma)$。自由钟的 $X_j(t)$ 是初始 $a_j,a_j^\dagger$ 的线性组合，系数模不超过 $d_T$。由

$$
\|g_j'\|_\infty\le2/\delta^2,\qquad
\|g_j''\|_\infty\le8/\delta^3
$$

及 $[N,a_j]=-a_j$ 可得，对每个 $V_L$（$P_L$ 与 $N$ 对易），

$$
\|[N,V_L(t)]\psi\|
 \le a\,\|(N+1)^{1/2}\psi\|+b\,\|\psi\|,
\quad
a=32d_4/\delta^2,\quad b=64d_4^2/\delta^3. \tag{10}
$$

推导中用
$[N,g(X_t)]=\bar\xi\,a^\dagger g'(X_t)-\xi g'(X_t)a$；
把 $a^\dagger$ 过到右边产生至多 $|\xi|^2\|g''\|$。两门范数 $\pi$ 与两钟求和均已经计入(10)。

令 $r_L(t)=\|(N+1)U_L(t)\Psi_0\|$。则
$r_L'\le a\sqrt{r_L}+b\le(a+b)\sqrt{r_L}$，初值1，因而

$$
\langle(N+1)^2\rangle_L
\le[1+2(a+b)]^4<4,\qquad 0\le t\le4, \tag{11}
$$

常数与 $L$ 无关。同法用闭图极限给实际过程的有限矩。这里不是从(7)的一阶矩擅自推出二阶矩。

### 4.3 极限与相对熵

有界交互的 Dyson 强极限给 $U_{p,L}\Psi_0\to U_p\Psi_0$。由(11)，

$$
\langle N1_{\{N>R\}}\rangle_L\le4/R,
$$

故 $\langle N\rangle_L\to\langle N\rangle_{\rho_E'}$；只用态范数收敛则没有这一结论。由(6)、$P_L$ 与场 boost 对易及 Duhamel 的交换子积分，还得到
$[\mathsf B,U_{p,L}]\Psi_0\to[\mathsf B,U_p]\Psi_0$，因此 $\langle\mathsf B\rangle_L\to\langle\mathsf B\rangle_{\rho_F'}$。

相对熵下半连续性应用于固定忠实 $\omega_W\otimes\tau_E$，把(9)运输为

$$
D(\rho'_{WE}\Vert\omega_W\otimes\tau_E)
\le\Delta K_W+\langle N\rangle+\log Z+\log16<\infty. \tag{12}
$$

这一步先证明了所需联合相对熵有限；没有使用“无限局域代数的相对熵对小迹距连续”。

现在才能用链式分解：

$$
D(\rho'_{WE}\Vert\omega_W\otimes\tau_E)
=D_W(\rho_F'\Vert\omega_W)+I(W:E)_{\rho'}
 +D(\rho_E'\Vert\tau_E).
$$

由(8)，最后一项等于
$\langle N\rangle+\log Z+\log16-S(\rho_E')$，全为有限量。相减得(1)的上界，三项均有限，不出现 $\infty-\infty$。**下半连续极限只给本轮所需不等式，不能将它写回等号。**

## 5. 后继是同一固定场读口，而非旧辅助记录

### 5.1 理想场信号及真空正概率

固定二结果场 channel $\mathcal M_u$，效应为 $n_u,1-n_u$，属于 $\mathcal A(W_8)$。先理想记录 $n_f$ 并不选择结果。CAR恒等式为

$$
n_f n_u n_f+(1-n_f)n_u(1-n_f)
=\frac12(n_{f_+}+n_{f_-}).
$$

两味真空占据互补，所以理想准备后的 $\mathcal M_u$ 概率为 $p_I=1/2$。这不是1228两个不同第二读口的比较，也不是把第一份辅助指针的数据计为场信息。

原真空同一 channel 的概率为 $s=\langle u,Su\rangle$。由1228高载频上界和逐动量正下界，

$$
s_*\equiv\frac1{4(K^2+B^2)}
\le s\le s^*\equiv\frac{B^2}{K^2}. \tag{13}
$$

下界证明：在 $\alpha_z\chi=\chi$ 上，负谱权重
$(1-p_z/\omega_p)/2\ge1/(4\omega_p^2)$，再用 Jensen 及
$\mathbb E\omega_p^2=K^2+B^2$。这里不是算符序 $P_-\ge cI$；仅是该旋量的矩阵元。$s$ 和 $1-s$ 都有严格正下界，未裁剪小概率。

### 5.2 真正实际过程与同一 channel 的区别

记 $p=\mathcal M_u(\rho_F')_1$。它是**实际准备终态**上的固定场 channel，数据处理直接给

$$
D_{\rm bin}(p\Vert s)\le D_W(\rho_F'\Vert\omega_W). \tag{14}
$$

实际第二实体指针在 $T=4$ 的概率记为 $\widetilde p$。因为第二钟前尾可能已与场相关、第一钟后尾仍在耦合，真实 $[2,4]$ 不是作用于 $\rho_F'$ 的独立固定 CP channel；本轮不把它说成如此。

设 $\epsilon_A(T)=LTd_T$，沿1228二步工具（包括脉冲期间 $H_F$）：

$$
|p-\tfrac12|\le
\alpha=2[\epsilon_A(2)+C\delta],\qquad
|\widetilde p-\tfrac12|\le
\gamma=2[\epsilon_A(4)+2C\delta]. \tag{15}
$$

准备经典比较中第二脉冲在 $T_p=2$ 前严格为零，实际前尾仍由 $\epsilon_A(2)$ 收费；全窗比较则两个窗口都在，旧第一钟后尾保留。这些是完整过程对比，未假定钟与场的中间独立性。因此 $|\widetilde p-p|\le\eta=\alpha+\gamma<10^{-15}$。

真空读出校准采用同一个 $G_{k_R}$、同一个新读出钟和指针，只省去第一准备门；原场准备为 $\Omega$。其实体概率 $\widetilde s$ 满足

$$
|\widetilde s-s|\le\zeta=
2[\epsilon_A(4)+C\delta]<10^{-16}.
$$

这仍是同一理想场 channel 的两份实际校准。它没有把相关旧辅助伪装为相同的真实只作用场的 channel；其差异由上述有限误差明示。

取粗正下界 $q_*=1/(5\cdot10^{12})<s_*$，$1-s>q_*$。二元熵连续性与对数均值估计给

$$
|D_{\rm bin}(\widetilde p\Vert\widetilde s)
 -D_{\rm bin}(p\Vert s)|
\le h_2(\eta)+2\eta\log(1/q_*)+\frac{\zeta}{q_*-\zeta}
<10^{-3}. \tag{16}
$$

这只认证有限记录KL，不声称无限局域相对熵连续。

$p,\widetilde p>0.49$，且 $s,\widetilde s<e^{-24}$，所以

$$
D_{\rm bin}(p\Vert s)>0.49\cdot24-\log2>11,\qquad
D_{\rm bin}(\widetilde p\Vert\widetilde s)>11. \tag{17}
$$

结合(1)、(14)—(17)，实际有限信息证书为

$$
D_{\rm bin}(\widetilde p\Vert\widetilde s)-10^{-3}
 +I(W:E)_{\rho'}
\le\Delta K_W+S(\rho_E').
$$

同时 $\Delta K_W>11-2.8=8.2$，这是同一实际准备过程的正余量，不是给两个无关模型的数值相减。

## 6. 同源 boost/能源预算与独立几何合同

### 6.1 实际到理想的加权误差

准备过程令 $\epsilon_C=\epsilon_A(2),\epsilon_p=C\delta$。经典与理想准备的场范数可取 $R_H=2C,R_B=2C_B$；实际范数取 $R_H^A=4C/\delta,R_B^A=4C_B/\delta$。两步运输给

$$
|\Delta E_F-\Delta E_{F,I}|
\le(R_H^A+R_H)\epsilon_C+2R_H\epsilon_p<10^{-6},
$$

$$
|\langle\mathsf B\rangle-\langle\mathsf B\rangle_I|
\le(R_B^A+R_B)\epsilon_C+2R_B\epsilon_p<10^{-6}. \tag{18}
$$

乘 $2\pi$ 后的 $\Delta K$ 误差也小于 $10^{-6}$。所有量来自实际同终态；不以小概率误差替代能源权重。为几何运输可取 $\Delta E_F<2C+10^{-6}$、
$|\langle\mathsf B\rangle|<2C_B+10^{-6}$，且 $|\langle P_z\rangle|\le\Delta E_F$。

同一场正常应力在自由区域保守；用(3)保证楔外改变为零，将全 boost 差写成首门切面的单侧源分布积分：

$$
\Delta K_W=
2\pi\int_{\mathbb R^2} dA
 \int_0^\infty r\,\Delta\langle T_{00}(r-8,x,y)\rangle\,dr. \tag{19}
$$

(19)是**同截面的重整化期望分布和生成元二次型身份**，没有声称任意尖时片 $T_{00}$ 都是可测自伴算符。这里的积分可用紧支测试序列定义：相互作用绘景的场支持并集为紧集，偶局部酉使其因果补上两点函数改变精确为零；$H^8$ 模式经自由输送保同阶局部正则性。在有限区间 Dyson 展开中，对两点函数差所需的一阶空间导数作用在这些涂抹模式上，每阶仅产生次数多项式因子，仍由有界交互的阶乘余项求和。因此正常序源差在本有限过程具有所需的一阶点分裂期望分布。取包含其紧支改变的测试函数，再按原Noether生成元二次型极限，得到(19)；不能仅从总boost期望有限就推广到任意局部应力算符。

有限区域控制体仍需保侧面通量；门期间来源由(2)交换子产生。它不是1228未来自由源泛函被替换成持续驱动未来历史，也不是全装置 boost 守恒。

### 6.2 独立几何及误差如何登记

基准几何输入为平直 Minkowski 度规、事前边界 $z=-8,t=1$、独立纵距 $r=z+8$ 以及横向面元 $dA=dx\,dy$。**这些是本候选的背景/校准合同，不是本轮从熵造出的测量器。** 没有借1225的Gauss读口，也没有定义 $A=S/s$。

为了给有限误差而非无误差背景，登记一个可独立检验的**均匀校准**合同：纵向原点误差 $|\delta a|\le10^{-12}$，时间原点误差 $|\delta t|\le10^{-12}$，横向面元公共比例误差 $d\widehat A=(1+\xi_A)dA$、$|\xi_A|\le10^{-12}$。不将此延拓到任意位置依赖的形变/密度误差。相应几何加权能源报告为

$$
\widehat{\Delta K}
=(1+\xi_A)\left[
 \Delta K_W+2\pi\delta a\,\Delta E_F
                  -2\pi\delta t\,\Delta P_z\right].
$$

由(18)的同源预算，

$$
|\widehat{\Delta K}-\Delta K_W|
\le10^{-12}|\Delta K_W|
 +(1+10^{-12})2\pi(2\cdot10^{-12})\Delta E_F
<10^{-2}. \tag{20}
$$

这是**若独立校准满足合同**时的误差运输，代码不能充当真实几何标定证据。面元误差之所以可以直接乘总源积分，是因为此处只登记均匀公共比例；任意非均匀误差不能用有符号积分的总值来界定。

最终得到有限、带明确几何前件的操作预算

$$
D_{\rm bin}(\widetilde p\Vert\widetilde s)-10^{-3}
 +I(W:E)_{\rho'}
\le\widehat{\Delta K}+10^{-2}+S(\rho_E').
$$

其信息下界11远大于合计误差0.011。

## 7. 能源、钟反作用与有限性范围

(2)的总 Hamiltonian 自伴正能且守恒。两钟初能与四寄存器基线沿1228完整保留，自由初能为
$M+1/(4M\sigma^2)+6$，严格总初能还须加非零但有界的 $\langle V(0)\rangle$；不能把Gaussian初尾设零。全窗 $T=4$ 的初末交互尾和反作用可用

$$
r_{\rm end}(4)=L\left(3\sigma+\frac{L\,4^2}{2M}\right)
$$

及同一位置跟踪预算界定；第一、第二脉冲的非零 Gaussian 尾从未设为零。实际场注能由两钟与交互能共同承担。几何模能源 $\Delta K_W$ 不是钟的物理总能，也不要求程序钟有物理 boost 荷。

代码检查仅使用旧参数、精确有理护栏、有限 CAR 恒等式以及误差运输，不对完整场相对熵作有限矩阵代算，不数值截断真空来声称无限域结论。辅助截断只用于第4节解析极限，实际(2)仍未截断。

本轮留下的真正接口是：独立几何读口的实际实现、超出有限区域核许可的逐点因果工程、真正共动温度计及来源、信息预算到面积变分的响应。条件性信息—boost关系本身已经在本具体实际菜单中获得有限共同域；不能把这些剩余项删掉后宣布面积律或引力完成。

