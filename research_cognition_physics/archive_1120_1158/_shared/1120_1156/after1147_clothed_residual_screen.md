# 1147后补充：正常窄取向衣着态的质量残差

2026-10-12。继续同一共同运动资格，未分配1148，不改旧资格或导航。保正常窄SO(3)取向包，不以均匀O替代。回查[1147实际证明](../../1147/proof.md)、[材料谱与源](../../1147/material_bound.md)及已冻结[矩与取向资格](after1147_moment_rotation_screen.md)。

**判断：所需两项衣着图范数均可严格给有限常数。** 这使真实有限I/M质量相对于冻结取向四支谱的残差任意小；不将冻结谱冒认实际质量精确谱，也未代替主线待办的裸CAR读口boost估计。

## 1. 冻结对象、完整自旋与统一谱圆

记 $D=H_{D,A}+H_{D,B}$、$H_0=D+H_\chi$。本稿将1147每个味荷对 $ab$ 的最低正模同时保两个 $m_j=\pm1/2$，裸等距为

$$
W_{0,ab}:\mathbb C^2_A\otimes\mathbb C^2_B
\longrightarrow\mathscr H_{ab},\qquad
W_{0,ab}(s_A\otimes s_B)
=b^\dagger(u_{a,s_A})b^\dagger(u_{b,s_B})
\Omega_{D,A}\Omega_{D,B}\Omega_\chi.
\tag{1}
$$

总自旋两份均保留；各站真空为旋转单态。每站最低能仅有这个自旋二重性；材料文件的下一不同正能与成对代价界，给每味荷对的秩4裸最低簇，至其余谱隙仍至少1147定义的g。并未将全部Fock截断到16维：这是完整正常空间中的低谱等距。

整球J为材料旋转标量，故冻结

$$
H_f(O)=D+H_\chi+\lambda\sum_iJ_i\chi(g_*(\cdot-Or_i))
\tag{2}
$$

与两站**分别的** $SU(2)_A\times SU(2)_B$ 强对易。两站相对槽取中心距3，与1147相同；源及海固定不变。

在部门ab记 $e_{ab}^{0}=e_a+e_b$、
$A_{0,ab}=H_0-e_{ab}^{0}\ge0$、
$Z_{ab}=A_{0,ab}+1$，
并取复谱圆 $\Gamma=\{|z|=g/2\}$。沿用1147的

$$
D_*=\frac{10}{\sqrt g}+\frac5g,\quad
R=\frac1{4D_*},\quad
q=|\lambda|D_*<1,\quad p=\frac q{1-q}.
\tag{3}
$$

实际λ满足 $q\le1/40000$，因此p远小于1。对

$$
\mathcal R_O(z)=(H_f(O)-e_{ab}^{0}-z)^{-1}
$$

有统一界

$$
r:=\sup_\Gamma\|\mathcal R_O(z)\|
\le\frac{2/g}{1-q},
$$

$$
s:=\sup_\Gamma\|Z_{ab}^{1/2}\mathcal R_O(z)\|
\le\frac{2/\sqrt g+2/g}{1-q}.
\tag{4}
$$

证明是 $\mathcal R_O=\mathcal R_0(I+\lambda V_O\mathcal R_0)^{-1}$；
$Z_{ab}^{1/2}\le A_{0,ab}^{1/2}+1$，且自由预解式分别由 $2/\sqrt g,2/g$ 控制。场数算符满足 $N_\chi+1\le Z_{ab}$，无需用全空间能源零点补一个未估常数。

设 $\Pi_{ab}(O)$ 是相应秩4 Riesz投影。其秩在小λ解析圆盘中不变，且 $\|\Pi_{ab}-\Pi_{0,ab}\|\le p$。其表示与裸 $(1/2,1/2)$ 表示等价，故Schur引理给

$$
H_f(O)\Pi_{ab}(O)=E_{ab}(\lambda)\Pi_{ab}(O)
\tag{5}
$$

的一条四重自旋简并实际本征值，和1147指定向上自旋的能源同一。不同味部门仍不混合；不存在为了旋转而删除钟频的操作。

## 2. 衣着等距的规范化不会引入未知角导数

材料J的独立旋转对称及场旋转协变给

$$
W_{0,ab}^\dagger\Pi_{ab}(O)W_{0,ab}
=\alpha_{ab}I_4,\qquad
1-p^2\le\alpha_{ab}\le1.
\tag{6}
$$

下界来自正交投影关系
$\|(I-\Pi)W_0s\|\le p\|s\|$。由于 $g_*$径向，

$$
H_f(O)=U_\chi(O)H_f(I)U_\chi(O)^\dagger,\qquad
U_\chi(O)W_{0,ab}=W_{0,ab}.
\tag{7}
$$

因此 $\alpha_{ab}$对O**精确不变**。取正重叠规范

$$
W_{\lambda,ab}(O)=
\alpha_{ab}^{-1/2}\Pi_{ab}(O)W_{0,ab}.
\tag{8}
$$

它是完整正常Hilbert空间中的实际本征等距，保未知自旋/味输入及外部R；满足
$W_{\lambda,ab}(O)=U_\chi(O)W_{\lambda,ab}(I)$。
这里不是按O重新准备真空。

## 3. 衣着向量自身的 $P_\chi^2$ 图范数

令 $X=H_0+1$。先前矩资格§2用热平滑J及光滑g证明了

$$
D((H_f(O)+1)^2)=D(X^2),
$$

$$
\|X^2\psi\|\le
2\|(H_f(O)+1)^2\psi\|+C\|\psi\|,
\tag{9}
$$

其中C可独立定义如下，全部常数仅依固定材料、源、λ而不依I、$M_0$或O：

$$
v_i=\omega^{-1/2}g_i/\sqrt2,\quad
c=2\sum_i\|J_i\|\|v_i\|,
$$

$$
d=2\sum_i
\big(\|[H_{D,i},J_i]\|\|v_i\|
+\|J_i\|\|\omega v_i\|\big),\quad
e=8\left(\sum_i\|J_i\|\|v_i\|\right)^2,
$$

$$
b=2|\lambda|c+|\lambda|d+\lambda^2e,\quad
L=\max(1,4b^2),\qquad C=2bL^{3/2}.
\tag{10}
$$

为明确本次逻辑：(9)不是对静止裸轨道的估计拿来直接替换衣着向量。衣着等距由(5)是真本征向量，故先有
$W_{\lambda,ab}(O)\subset D((H_f(O)+1)^2)$；
再用**算符域身份**(9)才得到其独立X图范数。能源满足
$|E_{ab}-e_{ab}^{0}|<g/2$，所以可取

$$
G_P=2(1+2e_{\max}+g/2)^2+C<\infty,
$$

$$
\|X^2W_{\lambda,ab}(O)\|\le G_P,\qquad
\|P_\chi^2W_{\lambda,ab}(O)\|\le G_P.
\tag{11}
$$

最后一步来自 $|P_\chi|\le H_\chi\le H_0$ 的共同自由谱，故
$P_\chi^2\le H_0^2\le X^2$ 也给相应作用范数控制。
这是全部正常χ及海的作用范数，不是动量均值为零，不是小态距对无界量的错误推广。

(9)的解析链仍保：$V$对 $X^{1/2}$有界，$[X,V]$对 $X^{1/2}$有界，$V^2$对X有界；
$(H_f+1)^2-X^2=2\lambda VX+\lambda[X,V]+\lambda^2V^2$ 对 $X^2$相对界0。该式在 $D(X^2)$上的自伴实现包含于自伴 $(H_f+1)^2$，两者遂相等；所有有限常数来自完整场估计。

## 4. Riesz角导数给正常窄包的 $L^2$ 图范数

令 $\mathcal X_a$为SO(3)左不变单位角导数，$\mathbf L^2=-\sum_a\mathcal X_a^2$。
对全部单位轴及二阶混合角导数，取有限常数

$$
c_1=2\sum_i\|J_i\|
\sup_{O,a}\|\mathcal X_a v_{i,O}\|,
$$

$$
c_2=2\sum_i\|J_i\|
\sup_{O,a,b}\|\mathcal X_a\mathcal X_b v_{i,O}\|.
\tag{12}
$$

球形g的旋转仅改变有限槽位 $Or_i$；其空间一、二阶导数平方可积，SO(3)紧致，因此这两个数明确有限、可由同一g积分求取。定义 $k_1=|\lambda|c_1,k_2=|\lambda|c_2$；因为 $N_\chi+1\le Z_{ab}$，

$$
\|(\mathcal X_aH_f)\xi\|\le k_1\|Z_{ab}^{1/2}\xi\|,\quad
\|(\mathcal X_a\mathcal X_bH_f)\xi\|
\le k_2\|Z_{ab}^{1/2}\xi\|.
\tag{13}
$$

相对共同域上这些系数C²；(4)使 $\mathcal R_O(z)$在算符范数下C²。
预解式导数为

$$
\mathcal X_a\mathcal R=-\mathcal R(\mathcal X_aH_f)\mathcal R,
$$

$$
\mathcal X_a^2\mathcal R=
2\mathcal R(\mathcal X_aH_f)\mathcal R
(\mathcal X_aH_f)\mathcal R
-\mathcal R(\mathcal X_a^2H_f)\mathcal R.
\tag{14}
$$

混合二阶含两个有序插入，同一范数界适用。设 $u_1=k_1s,u_2=k_2s$，用圆周长度 $|\Gamma|/(2\pi)=g/2$ 得

$$
\|\mathcal X_a\Pi_{ab}\|
\le\frac{u_1}{1-q},\quad
\|\mathcal X_a\mathcal X_b\Pi_{ab}\|
\le\frac{2u_1^2+u_2}{1-q}.
\tag{15}
$$

由(6)—(8)，归一因子无需角微分。因此可选

$$
w_1=\frac{u_1}{(1-q)\sqrt{1-p^2}},\qquad
w_2=\frac{2u_1^2+u_2}{(1-q)\sqrt{1-p^2}},
$$

$$
\|\mathcal X_a W_{\lambda,ab}\|\le w_1,\quad
\|\mathcal X_a\mathcal X_b W_{\lambda,ab}\|\le w_2.
\tag{16}
$$

这些是衣着等距的算符导数界；没有假设场角动量受场能源自动控制。

保选定正常窄包 $\psi_a(O)$，不改变它的宽度。复用1137已核的
$\ell_1=40/a,\ell_2=2000/a^2$，定义

$$
\Xi_{ab}s=\psi_a(O)W_{\lambda,ab}(O)s.
\tag{17}
$$

乘积法则、三轴Cauchy–Schwarz及 $2\sqrt3<4$ 给

$$
\|\mathbf L^2\Xi_{ab}\|\le
G_O:=\ell_2+4\ell_1w_1+3w_2<\infty,
\qquad
\|P_\chi^2\Xi_{ab}\|\le G_P.
\tag{18}
$$

其存在性包括完整正常源、全部海与任意自旋和被动R，不仅选定一个已知码向量。光滑窄包边缘平坦，SO(3)上无角δ态或补丁奇点。

## 5. 完整旋转交织与窄包下的真实S矩

将(8)按O直积分，定义 $\mathcal W_\lambda$ 从
$L^2(SO(3))\otimes\mathbb C^4_{\rm flavor}
\otimes\mathbb C^2_A\otimes\mathbb C^2_B$
进入完整实际内部空间。对共同旋转

$$
(U_{\rm int}(R)\Phi)(O)=
U_{D,A}(R)U_{D,B}(R)U_\chi(R)\Phi(R^{-1}O),
$$

有精确交织

$$
U_{\rm int}(R)\mathcal W_\lambda
=\mathcal W_\lambda
\big(U_{\rm left}(R)\otimes I_{\rm flavor}
\otimes D^{1/2}(R)\otimes D^{1/2}(R)\big).
\tag{19}
$$

证明用(7)、各站分别SU(2)对称及裸等距自旋交织；没有预设自旋0。对于固定窄包，得到

$$
\|\sqrt{S^2}\Xi\|
=\|\sqrt{(L+S_A+S_B)^2}(\psi_a\otimes s)\|
\le\ell_1+\sqrt2.
\tag{20}
$$

场衣着的角导数与O坐标变换在**共同S**中配对，不能把这个等式误称 $J_\chi$ 或L分别不变。窄包自身不是旋转不变态；(19)保其全部旋转像，正是后续诱导表示所需的真实旋转身份。

## 6. 有限惯性与材料惯性归入平方根时的实际残差

此节采用主线本次明确选定的父材料约定：将D放入支撑惯性质量。它和先前矩资格的“D在平方根外”是不同配方，不能悄然互换。

令 $m_S=M_0+m_O$，冻结参考为

$$
B=m_S+H_f(O),\qquad
B\Xi_{ab}=\mu_{ab}\Xi_{ab},\quad
\mu_{ab}=m_S+E_{ab}(\lambda).
$$

真实零K质量是

$$
A_{\rm act}=m_O+H_\chi+\frac{\mathbf L^2}{2I}
+\sqrt{(M_0+D)^2+P_\chi^2}+\lambda V(0,O)
=B+\frac{\mathbf L^2}{2I}+\mathcal R_M,
$$

$$
\mathcal R_M=
\sqrt{(M_0+D)^2+P_\chi^2}-(M_0+D).
\tag{21}
$$

D与 $P_\chi$ 强对易；$D\ge0$，故联合谱演算给
$0\le\mathcal R_M\le P_\chi^2/(2M_0)$，且相应作用范数由同式控制。
(11)、(18)使 $\Xi_{ab}\subset D(A_{\rm act})$，于是

$$
\boxed{\;
\|(A_{\rm act}-\mu_{ab})\Xi_{ab}\|
\le\rho:=\frac{G_O}{2I}+\frac{G_P}{2M_0}.
\;}
\tag{22}
$$

$G_P,G_O$不依I、$M_0$，故对任意指定 $\rho_*>0$，
取 $I\ge G_O/\rho_*$、$M_0\ge G_P/\rho_*$
即有 $\rho\le\rho_*$，均为有限正数。为实际质量严格正还需同时满足旧正性下界；可再增静质量，不改变残差常数。

实际质量对共同旋转强对易；因此 $\mathbf1_{|A_{\rm act}-\mu_{ab}|>h}$与S强对易。由谱定理，至少无权谱尾已有

$$
\|\mathbf1_{|A_{\rm act}-\mu_{ab}|>h}\Xi_{ab}\|
\le\rho/h,\qquad h>0.
\tag{23}
$$

(23)不说明实际质量恰有四条离散支；其可有连续谱或混合。也不自行控制带S的谱尾：若后续证明需要
$\|S\,\mathbf1_{\rm tail}\Xi\|$，
须另用足够光滑包的加权残差或仅用全态(20)加尾部有界算符策略，不能将无权(23)直接升级。

## 7. 本稿已关闭与未领取的边界

已关闭：当前完整正常衣着等距的 $P_\chi^2$、窄正常O的 $L^2$、真实完整自旋/旋转交织以及有限I/M实际质量残差。给出的常数为同一材料/剖面和旧g定义的解析量；本稿没有假装数值认证热迹、最短门时或设备可制造性。

仍由主线核验：谱窗有限boost与原裸CAR钟/记录的同一身份、正常位置末读的联合误差，以及全部有限资源反作用的实际合同。本稿没有靠“谱窗存在”领取这些结论，没有把有限尺寸刚体或MIT壁提升为认知新公理。静止矩、衣着矩、质量残差及运动读取是四个不同接口，不能互相替代。
