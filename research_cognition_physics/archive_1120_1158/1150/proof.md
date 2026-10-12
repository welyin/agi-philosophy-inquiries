# 1150证明：原长任务中的非对易干预、完整后态与整站区域收支

2026-10-12。按[准入](admission.md)检验有限有效共同接口。目标保持active。本轮不把有限带访问叫连续局部Dirac场，不签收共同Lorentz运输、自治脉冲装置或全装置四P。

## 1. 同一任务与实现选择

旧过程是1147的固定两站正常MIT材料与连续KG，完整自旋量词和裸读口复用1148§2。新过程在**相同的全正常CAR×KG Hilbert空间**上定义，自由$H_{D,A}+H_{D,B}$及所有未访问模式均保留。未访问模式初始仍为同一正常真空，没有迹掉海或场。

认知任务是：更换内部实现之后，原未知资料、实际钟、记录和声明中途操作能否继续，并由同一过程交代其能源责任。对应方法是保合同及后续能力，而非强制两种内部排班逐项相同；本轮不由社会类比推出物理律。

新增模型选择只有交互顶点的最低带投影和明确有限能源匹配。另声明一次**新增中途整站控制权限**，同一算符作用于新旧过程；它不是由各自H自主制造。四原则、CO、SR、D-P保持，新公理0。

## 2. 投影实际变成什么

复用[有限材料资格](../_shared/1120_1156/after1149_regulated_material_screen.md)中的直接Dirac验证。味1最低正负完整双重态可以取
$$
u_s=N\begin{pmatrix}F\xi_s\\iG\sigma\cdot\widehat x\,\xi_s\end{pmatrix},
\qquad
v_s=N\begin{pmatrix}-iG\sigma\cdot\widehat x\,\xi_s\\F\xi_s\end{pmatrix}.
\tag{1}
$$
$F=j_0(qr),G=bj_1(qr)$，边界$F=G$。两者满足同一MIT条件且能源分别为$\pm E_1$。逐点
$$
u_s^\dagger\beta u_t=\rho_s\delta_{st},\qquad
v_s^\dagger\beta v_t=-\rho_s\delta_{st},\qquad
u_s^\dagger\beta v_t=0,\quad
\rho_s=N^2(F^2-G^2).
\tag{2}
$$
因此同一正常序下
$$
:\bar\psi_{\min,1}(x)\psi_{\min,1}(x):_S
=\rho_s(x-y_i)N_i,\qquad
N_i=\sum_s(b_{i1s}^\dagger b_{i1s}+d_{i1s}^\dagger d_{i1s}),
\quad0\le N_i\le4.
\tag{3}
$$
正负模都在，却没有这一投影内的pair顶点。该模型不能代表旧完整海的同一响应机制；它只是一个可与旧输出比较的实现。$\rho_s$零延拓属于$H^1$，不是原g的$C^\infty$剖面，不能领取任意高能源矩。

取1149已经定义的正幅度$a_*=\kappa Q_g(1)/Q_\rho(1)$，$\rho=\rho_s/s_1$。新Hamiltonian为
$$
H_m=H_{D,A}+H_{D,B}+H_\chi+
\lambda a_*\sum_iN_i\chi(\rho_{s,i})
+\lambda^2[b(N_A+N_B)+2cI].
\tag{4}
$$
$b$是有限带能源参数，不宣称为完整Dirac局部质量反项。所有未选模式只自由演化。数源有界且相互交换；按其有限谱块作正规Weyl位移即可证明(4)在$D(H_D+H_\chi)$同域自伴、有下界。

## 3. 保原参数的匹配和资源

设旧单站二阶自系数为$L_a$，按1147的全海积分精确定义而非截断数值。中介质量1及旧$\|J\|<1,\|g\|_2^2<3$给
$-3/2<L_a\le0$。选择
$$
v=a_*^2\langle\rho_s,\omega^{-2}\rho_s\rangle,\qquad
c=L_0,\qquad b=L_1-L_0+v/2 .
\tag{5}
$$
若用资格稿记号，则$\delta m=b/s_1$，它仅表示同一有限带参数。这保背景能源常数，未按未知分支改变真空。它们是明确匹配输入，不是公理选出的常数。

原热迹包含4个正负最低模，故$\kappa<1/8$。$1\le Q_f(1)\le\sinh1$给$a_*<3/20$。又在$r\le1/\pi$上$F\ge5/6$，从归一积分得到$N^2\le27\pi^2/25<11$；$0\le\rho_s\le N^2,\int\rho_s=s_1<1$，所以
$$
\|\rho_s\|_2^2<11,\quad
|b|<1.624,\quad
\|b(N_A+N_B)+2cI\|<16 .
\tag{6}
$$
原$g\le1$及$R=[4(10/\sqrt g+5/g)]^{-1}$给
$$
0<\lambda\le R/10^4\le1/600000 .
\tag{7}
$$
本轮没有重新缩小λ、换T0或数值虚构热迹/完整谱隙。原选择另有$\lambda^2\le R^2\nu/(10^5g)$，而$R\le\sqrt g/40$，故
$$
\frac{|\lambda^2b|}{\nu}<\frac{1.624}{1.6\times10^8}<1.02\times10^{-8},
\qquad \nu_{\rm eff}=\nu+\lambda^2b>0.
\tag{7a}
$$
这也认证能源附件中输入00的脉冲功$\nu_{\rm eff}/2$为正。再用$C_2<3$得到$C_2\lambda^2/\nu<3/(1.6\times10^8)$，所以实际自由钟实验$\tau_c=\pi/(2\nu)$小于原T0。


独立于态距，场二次型界给
$$
H_m\ge H_{D,A}+H_{D,B}+\tfrac12H_\chi-32\lambda^2.
\tag{8}
$$
码×真空输入的总能至多$2E_{\max}+16\lambda^2$，脉冲前
$\langle H_\chi\rangle\le4E_{\max}+96\lambda^2$。这是实际分项资源有限性，不从后文小范数推断无界能源。

## 4. 全16维码的同一裸比较

在码上$N_i=P_{1,i}$。精确衣着位移为
$$
\alpha_{ab}=-\frac{\lambda a_*}{\sqrt2}\omega^{-3/2}
(a\rho_{s,A}+b'\rho_{s,B}),\qquad a,b'\in\{0,1\}.
\tag{9}
$$
这里$b'$为B味标签，与能源参数b不同。场未被取均值。两材料标签正交，衣着等距与原裸等距$W_0$的差统一未知spin及R：
$$
\|e^{-itH_m}W_0-W_0e^{-itK_m}\|
\le d_m:=2\sqrt2|\lambda a_*|\|\rho_s\|_2
<2.346\times10^{-6}\quad(t\in\mathbb R).
\tag{10}
$$
$K_m$含真实位移自能源。由(5)和1149的跨系数匹配，它**精确**等于旧码能源的二阶多项式：
$$
(K_m)_{ab'}=e_a+e_{b'}+
\lambda^2[L_a+L_{b'}-C_2ab'],\qquad C_2=C_2^{\rm old}(3).
\tag{11}
$$
位移解复用1144，不作新解法发现。

旧模型的每站旋转标量源使全双重态只是四标量谱乘真实spin恒等。1147的解析圆盘及1148§2的完整量词给
$$
\|e^{-itH_o}W_0-W_0e^{-itK_o}\|\le d_o=4/39999,
$$
$$
\|K_o-K_m\|\le r_o=
\frac g2\frac{(\lambda/R)^4}{1-(\lambda/R)^2}.
\tag{12}
$$
没有把一般局部Yukawa模型的高阶自旋混合忽略；(12)只用旧整站标量H。原固定参数满足
$$
T_0r_o\le r_T:=\frac{\pi}{2\cdot10^5(1-10^{-8})},\qquad
T_0=\frac{\pi}{C_2\lambda^2}.
\tag{13}
$$
因此相同全Hilbert空间上的真实输出满足
$$
\|(e^{-iT_0H_m}-e^{-iT_0H_o})W_0\otimes I_R\|
\le d_m+d_o+r_T
<1.19\times10^{-4}.
\tag{14}
$$
未选模式在新侧不激发，旧侧的全部高模/海衣着包含在$d_o$，没有通过删除环境取得小距离。

使用原spinblind裸CAR完整Lüders读取，误差由收缩性保持。连原读取容差，无干预记录误差小于0.00262，同一效果的记录差大于0.99476。实际钟的独立实验仍用原$\tau_c=\pi/(2\nu)$和原读口；由原参数界有$\tau_c<T_0$，故(14)亦给其附加误差。没有将钟改成新理想指针。

## 5. 真实非对易插入后，场和海不重置

在$\tau=T_0/2$，新旧两边实施同一合法全Fock酉
$$
U_A=\exp[-i(\pi/4)G_A],\qquad
G_A=\sum_s(b_{A0s}^\dagger b_{A1s}+b_{A1s}^\dagger b_{A0s}).
\tag{15}
$$
它保持总电荷，在码上$u_A=e^{-i(\pi/4)X_A}\otimes I_{\rm spin}$，
$\|[u_A,P_{1,A}]\|=1/\sqrt2$，且$U_AW_0=W_0u_A$。它不与原实际钟或源对易；不是把未知资料当经典标签。

定义两份实际完整历史
$$
\mathcal U_\#=e^{-iH_\#(T_0-\tau)}U_Ae^{-iH_\#\tau},
\qquad \#=m,o.
\tag{16}
$$
每段均以(10)、(12)替换为裸码比较，尾部经下一段**真实酉**传播，保全部场记忆。码生成元差用Duhamel，两个时长之和仍T0。得到
$$
\boxed{\|(\mathcal U_m-\mathcal U_o)W_0\otimes I_R\|
\le2(d_m+d_o)+r_T
<2.21\times10^{-4}.}
\tag{17}
$$
这对任意码—spin—R相干成立，同一终末完整仪器和同一后续不放大误差。不要求第一段输出精确回到码或真空。

**菜单必须区分：** 原未干预记录与(16)是两次不同实验；本轮不声称插入后仍执行原CZ。若换成完整X，中点可有echo；本轮角度为$\pi/4$，不以echo或其它改任务获得比较。没有证明所有任意多次干预均有相同固定误差。

## 6. 干预后的实际短窗因果

由(4)定义的实际H，而非瞬时Yukawa码H，承担因果命题。每站源支持单位球，中心距3。对任意重启时刻τ之后长度$0<s<1$，在自由相互作用绘景中，所有跨站时间点之差小于1，KG的Pauli–Jordan支持和独立材料代数给交互因子交叉强对易。复用1147的Trotter论证：
$$
U_I(s)=U_{A,I}(s)U_{B,I}(s)=U_{B,I}(s)U_{A,I}(s).
\tag{18}
$$
因而裸B效果的实际Heisenberg回拉只涉及B和KG，与所有A-only Kraus对易。任意已演化正常关联态$\sigma$及A-only迹保持操作$\mathcal E_A$满足
$$
\operatorname{tr}\big(E_B\,U_m(s)\mathcal E_A(\sigma)U_m(s)^\dagger\big)
=\operatorname{tr}\big(E_B\,U_m(s)\sigma U_m(s)^\dagger\big),
\quad0<s<1.
\tag{19}
$$
包括(15)，也包括保结果在A的仪器无条件统计；不提前广播结果，不领取后选归一条件态的无信号。无需把τ时场换回真空或解除已有相关。

这个因果接口以整站为操作分辨率。球内任意小区域权限没有获得；最小投影处处调用同一个$N_i$。有限谱投影的区域测试可覆盖整个截断代数，但逆Gram资源未受控，故不能将无限制子区反例直接扩成每个有限菜单失败。

## 7. 同一过程的做功和动态区域收支

完整推导见本轮[能源证明](energy_balance.md)。它对(4)、(15)及真实有限相干历史成立，并包含未知R。每段的场是有限个正常相干态的叠加；$\rho_s\in H^1$足以给有限经典Cauchy能源及相关交叉矩阵元，不假设点态等时能量算符存在。

以$h_i=a_*\rho_{s,i}$，写$K_i=H_{D,i}+\lambda^2(bN_i+cI)$及$V_i=\lambda N_i\chi(h_i)$。取固定$f\in C_c^\infty$，$f=1$覆盖A完整源，$f=0$覆盖B完整源。正常序场能和流只以有限相干历史上的期望/交叉矩阵元及时间弱形式定义：
$$
\mathcal E_A(t)=\langle K_A+V_A\rangle_t+
\frac12\int f\,\langle:\pi^2+|\nabla\chi|^2+\chi^2:\rangle_t ,
$$
$$
S=-\tfrac12:\{\pi,\nabla\chi\}:,\qquad
\Phi_f=-\int\nabla f\cdot\langle S\rangle .
\tag{20}
$$
段内KG受力项与交互存能变化相消，$\dot{\mathcal E}_A=-\Phi_f$。τ处的实际功为
$$
W_A=\langle U_A^\dagger(K_A+V_A)U_A-(K_A+V_A)\rangle_{\tau^-}
=\langle U_A^\dagger H_mU_A-H_m\rangle_{\tau^-}.
\tag{21}
$$
自由场能不跳，常数c不删掉但其脉冲差为0。对于时间测试$\zeta$：
$$
\boxed{-\int\zeta'(t)\mathcal E_A(t)\,dt
=-\int\zeta(t)\Phi_f(t)\,dt+\zeta(\tau)W_A.}
\tag{22}
$$
这保同一H中的交互存能；不只计算自由场通量后漏掉材料或控制功。没有把$\Phi_f$反向定义为$\dot H$，它由KG的标准流独立给出。材料反作用在源相关与做功中显式保留；空间固定壁仍是有效支撑输入。

资源另证：$U_A$与$H_D$交换子有界、与$H_\chi$对易，保共同能源域。旧侧也独立满足
$$
|W_o|\le2\nu+2\sqrt6|\lambda|\sqrt{4E_{\max}+2K_\lambda}.
$$
新侧
$$
|W_m|\le2\nu+2\sqrt2|\lambda a_*|\|\omega^{-1}\rho_s\|_2
\sqrt{4E_{\max}+96\lambda^2}+2|b|\lambda^2.
\tag{23}
$$
两边功均有限，后续资源由实际总能源及脉冲功接续。**没有从(17)推出无界能源或功的数值精度相等。** 末读需保原能源域权限；统计容差本身不足以保证资源。

## 8. 复算、复用与结果边界

[代码](check.py)/[结果](results.json)有三组：

1. 实际正负最低模、MIT边界与标量pair块消失；
2. 原λ/T0解析预算及投影匹配算符/能源下界常数；
3. 独立有限双振子诊断：非对易中插后的完整相干历史，与直接谱演化比较，保一个外部参考，核完整脉冲功及后段守恒。

第三项仅核位移/历史/做功公式，参数不是原连续场的数值替身，不模拟完整海或认证热迹。连续对象的结论由解析证明承担。三组不是科学实验校准。

1144的交换数源、1147的衣着、1148的spin读口、1149的幅度以及1138/1141的弱能源警戒直接复用。**新增连接**是原固定长任务中，同一真实中途操作之后的完整比较、因果短窗和同H区域做功—存能—通量同时成立。

尚未证明：连续Dirac局部父过程，任意子球权限，共同boost及原正常支撑，局部量子应力/能源仪器，脉冲与壁的自主四P。它们仍属于原阶段未完成要求，不能把本轮的有限整站接口改名为整个局域相对论场已生成。
