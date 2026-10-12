# 1122证明：一个守恒载荷码与两个可寻址载荷的区别

2026-10-11。接[1121](../1121/proof.md)的有限码开放项。继承F＋U＋C＋P、CO、SR及D-P；新增认知公理0、科学校准0。沿用1121明示的3+1平直、共同正质量、自由双味Dirac分支及其正能源真空表示，不再选另一物种或背景。

本轮证明两个不同范围的结论：一个局部CAR载荷码上的混味目标有完整四动量守恒的全局扩张；同一场内两个可独立寻址的载荷码，却不能在固定输入独立正常辅助及完整能量守恒下精确实现只转一个载荷的目标。第二项只排除此精确合同，不排除有限误差或其它码，也不是全认知体系的反模型。

## 1. 同一场、正能表示与明确测试资料

设单味Cauchy空间为$\mathcal K_0=L^2(\mathbb R^3,\mathbb C^4)$，

$$
h_D=-i\alpha^j\partial_j+\beta m,\quad m>0,\quad
\omega=\sqrt{-\Delta+m^2},\quad P_\pm=\frac12(I\pm h_D/\omega),\quad S=P_-.
\tag{1}
$$

两味的真空协方差为$S\otimes I_2$，真空向量记$\Omega$。实际物理Hamiltonian$\mathbb H\ge0$作用在粒子／空穴Fock空间$\mathcal F(\mathcal K_+\oplus\overline{\mathcal K_-})$，分别由$h_D|_+$及$-\overline{h_D}|_-$生成；空间动量及全部四平移与1121相同。这里不以满Cauchy空间上的$d\Gamma(h_D)$替代正能物理Hamiltonian。

取单位旋量$\chi$满足$\alpha_x\chi=\chi$。Dirac反对易关系给

$$
\chi^\dagger\beta\chi=\chi^\dagger\alpha_y\chi
=\chi^\dagger\alpha_z\chi=0.
\tag{2}
$$

取实、非负、非零、归一的径向$\varphi\in C_c^\infty(B_R)$；例如归一化的$\exp[-1/(1-|x|^2/R^2)]$在球内取此值、球外取零。令

$$
f(x)=\varphi(x)\chi,\quad g(x)=\varphi(x-Le_z)\chi,
\qquad L>2R.
\tag{3}
$$

两者的支撑严格分离、$\langle f,g\rangle=0$。在动量表象，$\chi^\dagger h_D(p)\chi=p_x$。对$r,s\in\{f,g\}$，$\widehat\varphi$径向而位移相位只依赖$p_z$，故含$p_x/\omega(p)$的积分关于$p_x$为奇函数。于是

$$
\langle r,P_\pm s\rangle=\frac12\delta_{rs},\qquad
\langle r,P_\pm e^{it\omega}s\rangle
=\frac12\langle\varphi_r,e^{it\omega}\varphi_s\rangle.
\tag{4}
$$

此处的$f,g$始终是紧支Cauchy测试模式，并非紧支纯正能粒子波函数。$P_\pm f,P_\pm g$一般有空间尾。

## 2. 单载荷：守恒扩张存在，但实际局部仪器尚未给出

令$c_i(r)$为第$i$味的局部CAR模式，$n_{ir}=c_i^\dagger(r)c_i(r)$。先仅使用$f$：

$$
e_0=(1-n_{1f})(1-n_{2f}),\quad
\psi_0=2e_0\Omega,\qquad J^{(1)}e_i=c_i^\dagger(f)\psi_0.
\tag{5}
$$

式(4)及准自由真空给$\|e_0\Omega\|^2=1/4$；$c_i(f)\psi_0=0$，CAR遂给$J^{(1)\dagger}J^{(1)}=I_2$。该码在固定奇宇称部门，任意逻辑叠加不混合宇称。

取$X=\begin{pmatrix}0&1\\1&0\end{pmatrix}$、$E=e^{-i\theta X}$，复用1121目标

$$
U_f=\exp[-i\theta(c_1^\dagger(f)c_2(f)+c_2^\dagger(f)c_1(f))].
\tag{6}
$$

空模式态被$U_f$固定，故$U_fJ^{(1)}=J^{(1)}E$。另一方面，对任意$Q\in SU(2)$，全场味旋转$G(Q)$在粒子及空穴部门采用互为共轭的味表示，固定$\Omega$并保全部实际四动量。$e_0$是两味全空投影，亦被$G(Q)$固定，故

$$
G(Q)J^{(1)}=J^{(1)}Q,\qquad
[G(Q),e^{ia\cdot\mathbb P}]=0.
\tag{7}
$$

由不可约双重态的Schur引理，

$$
J^{(1)\dagger}e^{ia\cdot\mathbb P}J^{(1)}=\phi(a)I_2.
\tag{8}
$$

因此$G(E)$是目标在此码上的四动量守恒扩张，且等距身份张量任意参考$R$仍成立。这是一项明确正例，说明1121全CAR范数不能直接移给任意有限码。

该存在性结论不声称$G(E)$是有限区域内实际执行的$U_f$。它作用于全场，且本轮未给$\psi_0$的完整确定性制备、控制器、定位来源或期限证书。式(5)是成功概率1/4的条件投影态，不是保余部原边缘的确定性CPTP重置；实际仪器和失败分支尚未认证。

## 3. 同一场内的双载荷码及其可寻址性

同时使用$f,g$，定义

$$
e_{00}=\prod_{r=f,g}\prod_{i=1,2}(1-n_{ir}),\quad
\psi_{00}=4e_{00}\Omega,\qquad
J(e_i\otimes e_j)=c_i^\dagger(f)c_j^\dagger(g)\psi_{00}.
\tag{9}
$$

四个模式的真空协方差恰为$I_4/2$，故$\|e_{00}\Omega\|^2=1/16$。CAR给$J^\dagger J=I_4$。四维码在固定偶宇称部门；$f$和$g$内各恰有一个局部CAR占据。两区域的偶混味门分别作用于两个逻辑因子：

$$
U_fJ=J(E\otimes I),\qquad U_gJ=J(I\otimes E),\qquad
G(Q)J=J(Q\otimes Q).
\tag{10}
$$

因此这确实保两个独立未知载荷以及任意参考的操作身份。式(9)是同一场的共同编码，不把两个单码的整个场态误当成张量积，也不把真空的其它关联丢掉。其条件投影成功概率1/16，同样没有领取确定性重置或原879／881固定仪器菜单的制备权限。代数上可独立寻址仍不等于实际独立控制器已造出。各码向量均为有限个光滑Cauchy场多项式作用于真空，因而落在任意有限次实际能源幂的定义域；障碍不依赖场码自身的无限能源。

定义有界能量特征压缩

$$
K(t)=J^\dagger e^{it\mathbb H}J.
\tag{11}
$$

因$\mathbb H$与全局味旋转对易，而$\mathbb C^2\otimes\mathbb C^2$分成重数均为一的单态和三重态，必有

$$
K(t)=A(t)I_4+B(t)\mathsf S,\qquad
\mathsf S(u\otimes v)=v\otimes u.
\tag{12}
$$

一般$K(t)$并不是四维酉群，也不等于$e^{itJ^\dagger\mathbb HJ}$。下一节计算实际完整场群酉的压缩，保留码外传播。

## 4. 粒子／空穴字典与交换系数

令$u_r=\sqrt2P_+r$、$v_r=\sqrt2\overline{P_-r}$；各自是两组正交归一物理模式。每味有

$$
c(r)=\frac{b_r+d_r^\dagger}{\sqrt2}.
\tag{13}
$$

按$(b_f,b_g,d_f,d_g)$确定Fock顺序。一味两模式全空的归一态及其占据态为

$$
\begin{aligned}
\xi_0&=\frac12(1-b_f^\dagger d_f^\dagger)(1-b_g^\dagger d_g^\dagger)\Omega,\\
\xi_f&=\frac1{\sqrt2}b_f^\dagger(1-b_g^\dagger d_g^\dagger)\Omega,\\
\xi_g&=\frac1{\sqrt2}b_g^\dagger(1-b_f^\dagger d_f^\dagger)\Omega,\\
\xi_{fg}&=b_f^\dagger b_g^\dagger\Omega.
\end{aligned}
\tag{14}
$$

定义

$$
p_0(t)=\langle\varphi,e^{it\omega}\varphi\rangle,\quad
p_L(t)=\langle\varphi,e^{it\omega}\varphi(\cdot-Le_z)\rangle,\quad
C(t)=\begin{pmatrix}p_0&p_L\\p_L&p_0\end{pmatrix},\quad D(t)=p_0^2-p_L^2.
\tag{15}
$$

由(4)，$C$是物理粒子群酉在$u_f,u_g$的压缩；空穴群酉在$v_f,v_g$的压缩也为$C$。后一个身份需保共轭空间：其矩阵元为$2\langle s,P_-e^{it\omega}r\rangle$，由本模式的反射对称性等于$C_{rs}$。**所以这里出现$p_0^2,p_L^2$，不是$|p_0|^2,|p_L|^2$；同一正向时间的粒子和空穴均携正能源。** $C$通常是收缩而非酉矩阵。

自由Fock传播的等粒子数矩阵元由对应一粒子矩阵元的行列式给出。将(14)按1粒子、3粒子等正交部门相加，得到

$$
\begin{aligned}
F_{ff}=F_{gg}&=\langle\xi_f,e^{it\mathbb H}\xi_f\rangle
=\tfrac12p_0(1+D),\\
F_{fg}=F_{gf}&=\langle\xi_f,e^{it\mathbb H}\xi_g\rangle
=\tfrac12p_L(1-D),\\
F_{00}&=\langle\xi_0,e^{it\mathbb H}\xi_0\rangle
=\tfrac14(1+2p_0^2+2p_L^2+D^2),\\
F_{fg,fg}&=\langle\xi_{fg},e^{it\mathbb H}\xi_{fg}\rangle=D.
\end{aligned}
\tag{16}
$$

特别地，$\xi_f$的3粒子项为$-b_f^\dagger b_g^\dagger d_g^\dagger\Omega/\sqrt2$，$\xi_g$的3粒子项为$+b_f^\dagger b_g^\dagger d_f^\dagger\Omega/\sqrt2$，故交叉项是$p_L-Dp_L$，不是$p_L+Dp_L$。

把两味按味1、味2排序，$J_{12}$对应$+\xi_f\otimes\xi_g$，$J_{21}$对应$-\xi_g\otimes\xi_f$；该负号来自一次费米交换。同味输入则对应$\xi_{fg}\otimes\xi_0$或反序。因此

$$
\boxed{A=F_{ff}^2=\tfrac14p_0^2(1+D)^2,\qquad
B=-F_{fg}^2=-\tfrac14p_L^2(1-D)^2.}
\tag{17}
$$

同味对角也一致，因为$A+B=DF_{00}$。这额外检查了味分解及Fock符号，不能用将$p_0,p_L$取绝对值的概率计算替代。

## 5. 严格分离支撑仍给非零四阶交换项

置

$$
\varepsilon_0=\langle\varphi,\omega\varphi\rangle>0,\qquad
\varepsilon_L=\langle\varphi,\omega\varphi(\cdot-Le_z)\rangle.
\tag{18}
$$

因$\omega^2=-\Delta+m^2$为局部微分算子且两支撑分离，$\langle\varphi,\omega^2\varphi_L\rangle=0$。光滑紧支保证所需任意有限阶能源矩存在，故

$$
p_0(t)=1+i\varepsilon_0t+O(t^2),\qquad
p_L(t)=i\varepsilon_Lt+O(t^3),\qquad
F_{fg}(t)=\varepsilon_0\varepsilon_Lt^2+O(t^3).
\tag{19}
$$

$\varepsilon_L$严格负而不是仅一般非零。利用正自伴算子的平方根热核公式以及分离支撑，

$$
\varepsilon_L=-\frac1{2\sqrt\pi}\int_0^\infty
s^{-3/2}e^{-m^2s}\langle\varphi,e^{s\Delta}\varphi_L\rangle\,ds<0.
\tag{20}
$$

热核为严格正；两个非负非零函数给每个$s>0$严格正的内积。支撑距离$L-2R>0$给$s\downarrow0$端的高斯抑制，$m>0$控制无穷端，积分收敛。公式来自$\sqrt A=(2\sqrt\pi)^{-1}\int_0^\infty(1-e^{-sA})s^{-3/2}ds$，其中恒等项的交叉矩阵元为零。

所以

$$
B(t)=-\varepsilon_0^2\varepsilon_L^2t^4+O(t^5),\qquad
B(t)\ne0\quad\text{对充分小的 }t\ne0.
\tag{21}
$$

这项四阶行为是本条件空模式码的一粒子和三粒子转移抵消后的结果；若改用纯正能双粒子码，通常从二阶开始，那已是另一编码。有限矩阵程序不承担(20)的连续场和严格支撑量词。

由$[\mathsf S,E\otimes I]=(I\otimes E-E\otimes I)\mathsf S$可得

$$
\|[K(t),E\otimes I]\|=2|B(t)||\sin\theta|.
\tag{22}
$$

因此当$\theta\notin\pi\mathbb Z$时，选择性门不保持完整能量特征压缩。$\theta\in\pi\mathbb Z$在本码上只给标量，逻辑通道为恒等，必须排除；这与1121全CAR门仅在$2\pi\mathbb Z$恒等的范围不同。

## 6. 固定正常辅助与精确合同：无需额外矩假设

假设辅助$C$初态$\eta_C$是与全部未知双载荷及参考独立的固定正常态；联合幺正$W$完整保能量，即强意义下

$$
[W,e^{it(\mathbb H+H_C)}]=0\quad(t\in\mathbb R).
\tag{23}
$$

要求对全部码输入$\rho$精确实现

$$
\operatorname{Tr}_C W(J\rho J^\dagger\otimes\eta_C)W^\dagger
=JV\rho V^\dagger J^\dagger,\qquad V=E\otimes I.
\tag{24}
$$

这里输出包含整个场态，而不只是末端测量概率或仅一个逻辑边缘。式(24)是明确较强的精确编码合同；其完全正通道身份自动保任意被动参考。

先将$\eta_C$纯化成$|a\rangle\in C\otimes R_C$，在$R_C$上令守恒生成元为零、$W$取恒等。纯输入输出的场态为纯态，故$W(J|z\rangle\otimes|a\rangle)=JV|z\rangle\otimes|b_z\rangle$。比较一组正交基及其叠加，线性性迫使共同$|b_z\rangle=|b\rangle$（统一吸收相位）。**允许辅助最终改变，不要求复原**。

记$\phi_a(t)=\langle a|e^{itH_C}\otimes I|a\rangle$、$\phi_b(t)$类似。将(23)夹在这些输入向量之间得矩阵恒等式

$$
\phi_a(t)K(t)=\phi_b(t)V^\dagger K(t)V.
\tag{25}
$$

在$t=0$有$K(0)=I_4$、$\phi_a(0)=\phi_b(0)=1$。有限维压缩和正常态的群特征函数连续，故近零$\operatorname{Tr}K(t)\ne0$且$\phi_a(t)\ne0$。对(25)取迹，先得$\phi_a(t)=\phi_b(t)$，再得$K(t)=V^\dagger K(t)V$，与(21)—(22)矛盾。

因此对$\theta\notin\pi\mathbb Z$，(23)—(24)不能同时成立。只保总能量已经足够排除，因此完整四动量守恒也不能挽救这一精确合同。论证全程用有界群酉，不假设辅助初末能量期望或方差有限；其正常性和强连续性是关键。它不排除非正常理想参考，也不把允许初始输入相关辅助的另一合同混入。

一般精确酉实现的守恒限制属于成熟WAY理论；[Kuramochi—Tajima，arXiv:2208.13494，Theorem 2](https://arxiv.org/html/2208.13494)给无界连续守恒量的精确酉限制。本轮不声称发现新的普遍WAY定理，也不把其全系统幺正结论直接套到不被$\mathbb H$保持的码。式(25)给本嵌入码所需的自足有界谱论证；净新增对象是同一连续场内(9)的精确可寻址合同及(17)—(21)的非零证书。

## 7. 结论范围、复算与停止线

单码正例和双码障碍并不矛盾：全局旋转在单码上恰为所需门，在双码上却同时旋转两个载荷。单载荷保四动量的扩张存在性，不自动保证多个空间分离载荷可独立调用同样权限。只核$J^\dagger\mathbb HJ$或少量均值也会漏掉本例四阶才出现的谱差别。

本轮没有证明这些码能由原仪器菜单实际制备，没有给有限区域控制器、定位反作用、原CO端点身份、完整应力或有效精度预算。尤其本排除只针对精确式(24)：有限误差门、不同编码、改变场输出合同、输入相关准备及有资源预算的共同过程均另审。它不是F＋U＋C＋P／CO／SR／D-P的全部共同模型反证。

[check.py](check.py)用有限CAR产生湮灭算子和自由Fock外积传播独立核等距、局部门身份、Fock符号、(16)—(17)、全局守恒扩张及特征压缩障碍；[results.json](results.json)保存结果。有限模型仅校准这些代数恒等式，不称连续Dirac离散逼近、不以格点代替严格紧支或热核证据。默认运行只读比较，不改保存结果。方向去重和下一步分别见[selection_audit.md](selection_audit.md)、[NEXT.md](NEXT.md)。
