# 1129证明：连续材料码的正常自治标量来源

2026-10-11。接[1127](../1127/proof.md)的连续Dirac双材料双轨正常码，直接构造一个固定有限分辨率的自治标量来源。事前要求三个条件共同成立：完整未知Q及任意被动R的交接半diamond误差小于0.01；声明有限观测菜单的锥外概率变化小于0.001；同一个H保留实际来源后态、能源、反作用及固定支撑的空间动量责任。本例全部通过。

新增认知公理0，科学校准0。这里的新增量是连续材料、正常量子来源和这三个合同的共同实例；一般Duhamel、CAR、自由波动方程和有界扰动工具不另计新轮。实际CO/SR钟、动态定位架、准备与读取仪器以及共同Lorentz仍未完成，阶段active。

## 1. 模型输入、证明层级与固定数据

沿用1127的标准正能Dirac Fock表示、完整正负谱及CAR，固定自然单位、质量$m=1$、终时$T=1$。单粒子空间为

$$
\mathcal K=L^2(\mathbb R^3,\mathbb C^4)\otimes\mathbb C^2_{AB}\otimes\mathbb C^2_{rail},\quad
h=(\boldsymbol\alpha\cdot p+\beta)\otimes I,\quad S=1_{h<0}.
\tag{1}
$$

$H_F\ge0$为完整粒子／空穴自由能源，$[H_F,c^\dagger(v)]=c^\dagger(hv)$对$v\in D(h)$成立。材料A、B及双rail是真实模型标签；没有只保一粒子态或把局部真空当空模。另取独立的三维质量零实标量场Fock，$H_\phi=d\Gamma(|k|)\ge0$。它是有效模型输入，未由认知公理唯一推导。

固定$0\le\eta\le1$、$\eta\in C_c^\infty$，在$B_4$上为1、支撑于$B_5$；$0\le\beta_0\le1$、$\beta_0\in C_c^\infty$，在$B_7$上为1、支撑于$B_8$。具体令$R(u)$为$e^{-1/[s(1-s)]}1_{0<s<1}$从0到$u$的归一累积积分（两侧延拓为0、1），取$\eta(x)=1-R(|x|-4)$、$\beta_0(x)=1-R(|x|-7)$。下标避免与Dirac矩阵$\beta$混淆。源与材料准备、共同切片、坐标和静态部署函数均明确外给。

令

$$
q(s)=e^{-1/(1-s^2)}1_{|s|<1},\quad
f(x)=\frac{q(x_1)q(x_2)q(x_3)}{\|q\|_2^3}\xi,\quad
\|\xi\|=1,\quad\alpha_z\xi=\xi.
\tag{2}
$$

则$f\in C_c^\infty(B_2,\mathbb C^4)$、$\|f\|=1$。$\|q\|_2^2\ge e^{-8/3}>1/16$，$|q'|\le8e^{-2}<2$，所以$\|q'\|_2^2<8$、$\|\nabla f\|^2<384<400$，从而$\|hf\|<21$。实空间振幅的动量密度为偶函数，$\xi^\dagger\beta\xi=\xi^\dagger\alpha_{x,y}\xi=0$，故$\langle f,Sf\rangle=1/2$。

取实、偶、非负、光滑径向核

$$
\rho(x)=Z^{-1}e^{-1/(1-|x|^2)}1_{|x|<1},\quad
\rho_\epsilon(x)=\epsilon^{-3}\rho(x/\epsilon),\quad K_\epsilon v=\rho_\epsilon*v,
\quad\epsilon=10^{-5},\quad A=10^{34}.
\tag{3}
$$

$\int\rho=1$、$\mu_1=\int|x|\rho(x)dx\le1$；$Z\ge(\pi/6)e^{-4/3}>1/8$，且原bump平方不超过自身，故$\|\rho\|_2^2\le8$。于是

$$
K_\epsilon=K_\epsilon^\dagger,\quad\|K_\epsilon\|\le1,\quad
\|(K_\epsilon-I)v\|\le\epsilon\|\nabla v\|,\quad
q_\epsilon:=\|\rho_\epsilon\|_2^2\le8\epsilon^{-3}.
\tag{4}
$$

门函数直接采用1128的解析归一函数：$b(s)=e^{-1/[s(1-s)]}1_{0<s<1}$，$I_b=\int_0^1b$，$F(z)=2b(2z-1/2)/I_b$。因此

$$
\operatorname{supp}F=[1/4,3/4],\quad F\ge0,\quad\int F=1,
\quad\|F\|_\infty\le16,\quad L:=\|F'\|_\infty\le512,\quad\lambda=\pi/2.
\tag{5}
$$

所有面积、支撑、范数及后述无限维结论由解析论证承担。代码只复算精确预算、独立核／谱诊断及有限CAR／有界源机制，不把小矩阵、数值积分或离散波动当连续证书。有限$\epsilon$是本有效模型的分辨率，未声明物理最小尺度或$\epsilon\to0$的重整化极限。

## 2. 正常源、实际有界交互与自治H

写$\rho_{\epsilon,x}(y)=\rho_\epsilon(y-x)$，$\phi_\epsilon(x)=\phi(\rho_{\epsilon,x})$。按旋量$s$、rail$i$求和，定义完整CAR偶元

$$
J_\epsilon(x)=\sum_{s=1}^4\sum_{i=1}^2
\left[c^\dagger_{Ais}(\rho_{\epsilon,x})c_{Bis}(\rho_{\epsilon,x})+
c^\dagger_{Bis}(\rho_{\epsilon,x})c_{Ais}(\rho_{\epsilon,x})\right].
\tag{6}
$$

每个相互正交的A/B对的hopping范数为$q_\epsilon$，故$\|J_\epsilon(x)\|\le8q_\epsilon$。完整模型是

$$
H=H_0+V,\quad H_0=H_F+H_\phi,\qquad
V=\lambda\int\eta(x)F(\phi_\epsilon(x)/A)J_\epsilon(x)\,d^3x.
\tag{7}
$$

积分取强算符积分（也等价于其有界弱积分）；平滑核的Fock正则表示给强可测性，统一范数界及$\eta\in L^1$给可积性。这里不假设$x\mapsto F(\phi_\epsilon(x))$范数连续。标量和材料因子对易，$V$自伴且

$$
C:=8q_\epsilon\|\eta\|_1\le8\cdot8\epsilon^{-3}\cdot524=:C_* =3.3536\times10^{19},
\qquad\|V\|\le16\lambda C_*.
\tag{8}
$$

因$H_0\ge0$自伴，$H$在$D(H_0)$自伴、半有界。这是完整Fock结论，不以一粒子能源代替全场能源，也不涉及裸的逐点场乘积。

初始源$|a\rangle$为经典Cauchy均值$\phi_{cl}(0)=0$、$\pi_{cl}(0)=A\beta_0$对应的相干向量。其位移振幅正比$iA\widehat\beta_0(k)/\sqrt{2|k|}$，在三维低频平方可积，高频因$\beta_0$光滑紧支快降；故是正常Fock向量并属于$D(H_\phi)$，而非非正常无穷平台态。其能源

$$
E_\phi(0)=\frac{A^2}{2}\|\beta_0\|_2^2\le1073A^2.
\tag{9}
$$

有限传播使$0\le t\le1$、$x\in\operatorname{supp}\eta$上的平滑均值恰为$At$：所需整个过去球半径不超过$5+1+\epsilon<7$。在自由相干态$a_t=e^{-itH_\phi}a$上，$\phi_\epsilon(x)$的标准差与真空相同。使用酉Fourier约定，分割$|k|\le1/\epsilon$及其补集，得到

$$
v_\epsilon^2=\frac12\int\frac{|\widehat\rho_\epsilon(k)|^2}{|k|}dk
\le\frac{4+1/(8\pi^2)}{\epsilon^2}<\frac9{\epsilon^2},
\quad
\|(\phi_\epsilon(x)/A-t)a_t\|=v_\epsilon/A.
\tag{10}
$$

低频用$|\widehat\rho|\le(2\pi)^{-3/2}$，高频用Plancherel和$\|\rho\|_2^2\le8$；没有删除红外或真空涨落。

## 3. 连续材料比较及完整未知码

仅作比较，令源因子换为$F(t)$，材料Hamiltonian为$H_F+\lambda F(t)B_\epsilon$，其中$B_\epsilon=\int\eta J_\epsilon$本身有界。其一粒子生成元是

$$
h+\lambda F(t)K_\epsilon\eta K_\epsilon\otimes\sigma_x\otimes I_{rail}.
\tag{11}
$$

相应完整材料传播$\widehat U_\epsilon(t)$由同一有界CAR生成元定义，实施该一粒子传播$U_\epsilon(t)$，并固定其相位；无需另挑外场lift或删除其负谱块。

令$u_t=e^{-iht}f$，$F_i=f\otimes e_A\otimes e_i$，$G_i=u_1\otimes e_B\otimes e_i$。$u_t$支撑于$B_{2+t}$，其$K_\epsilon u_t$支撑于$B_{3+\epsilon}\subset B_4$，所以沿理想旋转轨道$K_\epsilon\eta K_\epsilon u_t=K_\epsilon^2u_t$。$h$与梯度对易，(4)和Duhamel给

$$
\|U_\epsilon(1)F_i-(-iG_i)\|\le
2\lambda\epsilon\|\nabla f\|=:d_m\le40\lambda\epsilon.
\tag{12}
$$

$G_i$严格位于B材料的$B_3$区域；实际耦合及读取仍有限。理想全空间旋转只用来估算两指定模式，没有被授予实际控制。

正常源码逐字采用1127结构：

$$
n_i=c^\dagger(F_i)c(F_i),\quad e_0=(1-n_1)(1-n_2),\quad
\Omega_0=2e_0\Omega,\quad J_i=c^\dagger(F_i)\Omega_0.
\tag{13}
$$

$\|e_0\Omega\|^2=1/4$、$c(F_i)\Omega_0=0$、$J^\dagger J=I_2$。两列同奇宇称，全部逻辑偶代数和任意被动参考R允许。准确正常准备仍是输入，不把概率$1/4$当作已实施的准备仪器。

复用1128 §4的有限CAR引理：对单位$v,w$存在支撑于至多两模式的SU(2)旋转$R$，$Rv=w$，其完整CAR实施满足$\|\Gamma(R)-I\|=\|v-w\|$，包括复overlap和纯相位情形。在两个独立rail分别校正(12)，得到偶CAR酉$\mathcal R$，$\|\mathcal R-I\|\le2d_m$。令$\Xi=\mathcal R\widehat U_\epsilon(1)\Omega_0$，则

$$
c(G_i)\Xi=0,\quad
\mathcal R\widehat U_\epsilon(1)J_i=-i c^\dagger(G_i)\Xi.
\tag{14}
$$

目标两模式单占据的全部矩阵单位精确重现Q；有限CAR因子化给一个与Q无关但实际改变的余部。$\mathcal R$是证明比较，不是隐藏的物理装置。B的确定性解码把非单占据送入正交失败标记，故

$$
\tfrac12\|\mathcal N_\epsilon-\operatorname{id}\|_\diamond
\le\delta_m:=2d_m\le80\lambda\epsilon<0.001257143.
\tag{15}
$$

没有后选择删除失败支路，没有要求整个真空输出逼近一个未改变的真空。

## 4. 正常自治来源的完整联合比较

实际输入为$J\rho_{QR}J^\dagger\otimes|a\rangle\langle a|$，未知Q不改变H或来源。参考传播是$\widehat U_\epsilon(t)\otimes e^{-itH_\phi}$。谱演算的标量Lipschitz估计与(10)给

$$
\|[F(\phi_\epsilon(x)/A)-F(t)]a_t\|\le Lv_\epsilon/A.
\tag{16}
$$

在完整参考态上使用有界扰动Duhamel，真实传播因子酉，$J_\epsilon(x)$保持完整CAR范数，遂对任意材料输入和任意被动参考一致有

$$
\left\|\left(e^{-itH}-\widehat U_\epsilon(t)\otimes e^{-itH_\phi}\right)
(I\otimes|a\rangle)\right\|
\le\delta_s(t):=\lambda C_*L\frac{3t}{\epsilon A},\qquad
\delta_s:=\delta_s(1)<8.094633\times10^{-7}.
\tag{17}
$$

左端是从全部材料Hilbert空间到实际联合输出的等距差，保来源相位、后态及其相关性，强于只比较两个源码基输入。参考源自由演化只是比较过程，实际源始终由H运输。

完整等距差控制任意被动R的半迹距；再接固定解码及(15)，得到事前合同

$$
\boxed{\varepsilon_{handoff}:=\tfrac12\|\mathcal N_{actual}-\operatorname{id}\|_\diamond
\le\delta_m+\delta_s<0.001257953<0.01.}
\tag{18}
$$

理想比较输出包含固定$\Xi$及$a_1$，实际联合输出靠近它但不等于它；不据此删除来源或宣称无反作用。

## 5. 本次有限因果观测菜单

合同事前固定任一有限组$(X,Y,W_X,g)$：$X,Y$有界区域，$W_X$为由有限个支撑于X的光滑CAR模式生成的偶酉，$g$为支撑于Y的归一光滑材料模式，测量effect $n_g=c^\dagger(g)c(g)$；可包括X中的rail交换。准备同一正常源$a$，比较任意正常材料态$\varrho$与$W_X\varrho W_X^\dagger$，被动R原样保留。需要有限能源时再要求$\varrho$有限能源；上述光滑有限CAR门保此性质。这里不准任意源干预、选择性过去准备或整个无限局部代数的diamond菜单。

一个具体单项菜单为$X=B_1(0)$、$Y=B_1(4e_z)$，$\kappa(x)$取(2)空间振幅的归一压缩$2^{3/2}\prod_jq(2x_j)/\|q\|_2^3$；$a_i=\kappa\xi\otimes e_A\otimes e_i$，$g=\kappa(x-4e_z)\xi\otimes e_B\otimes e_1$，$W_X=\exp[-i(\pi/2)(c^\dagger(a_1)c(a_2)+c^\dagger(a_2)c(a_1))]$。它在X切换两个未知rail，$\operatorname{dist}(X,Y)=2>1+20\epsilon$；接收菜单不是空的形式量词。

对(11)，$B=K_\epsilon\eta K_\epsilon\otimes\sigma_x$范数不超过1、核空间跳跃半径不超过$2\epsilon$。把传播展开为夹着自由Dirac段的Dyson级数。含n次B的项只能跨越总自由距离$t$加$2n\epsilon$；自由段长度之和恰为$t$，不把相互作用绘景中的膨胀支撑误算成额外传播。非负$F$且$\lambda\int_0^tF\le\lambda$。因此若

$$
\operatorname{dist}(X,Y)>t+20\epsilon,\quad 0\le t\le1,
\qquad
\|1_XU_\epsilon(t)^\dagger g\|\le
r_{10}:=\sum_{n\ge10}\frac{\lambda^n}{n!}
\le\frac{\lambda^{10}}{10!}\frac1{1-\lambda/11}.
\tag{19}
$$

前十项为零；严格距离实际上允许更紧的起始项，此处保守保留$n=10$。令$v=U_\epsilon(t)^\dagger g$，$q=1_{X^c}v$。由于$\|v\|=1$、$\|q\|\le1$，CAR给$\|n_v-n_q\|\le2\|v-q\|\le2r_{10}$。偶$W_X$与$n_q$对易，故参考两支的概率差不超过$4r_{10}$。用(17)对两支各加一次实际联合等距误差，得到

$$
\boxed{|p_{actual}(\varrho)-p_{actual}(W_X\varrho W_X^\dagger)|
\le4r_{10}+2\delta_s<0.000119703<0.001.}
\tag{20}
$$

该界也成立于$n_g\otimes E_R$，$0\le E_R\le I$，不随有限菜单大小或被动参考维数变化；未签全部effect联合统计的总变差相同界。$20\epsilon$是本菜单的有限空间余量，不是全局精确光锥。源自身的自由相对论性不自动消除材料平滑核的有限尾，也不产生任意源控制下的因果证书。

## 6. 同H的有限输入／输出能源

对(13)，全Fock交换子给$\|[H_F,n_i]\|\le2\|hf\|$、$\|[H_F,e_0]\|\le4\|hf\|$。$\Omega\in D(H_F)$且$H_F\Omega=0$，所以$\Omega_0,J_i\in D(H_F)$，

$$
\|H_F\Omega_0\|\le8\|hf\|,\qquad
\|H_FJ_i\|\le9\|hf\|,\qquad
\sup_{\|z\|=1}\|H_FJz\|\le9\sqrt2\cdot21<268.
\tag{21}
$$

这些来自有界CAR在自由生成元域上的交换子，不仅是一粒子Sobolev判据。结合相干源的$D(H_\phi)$正则性，全部编码输入在$D(H_0)=D(H)$。自治传播保$D(H)$，故实际输出也在$D(H_0)$，特别所有时刻材料和源均能有限。更直接地，由正性及有界V，

$$
\langle H_0\rangle_t\le268+1073A^2+2\|V\|,
\quad
\Delta\langle H_F\rangle+\Delta\langle H_\phi\rangle+\Delta\langle V\rangle=0.
\tag{22}
$$

这使用同一完整H的能源守恒，不从普通通道误差推无界能量，也未假设$H_F$或$H_\phi$分别守恒。资源巨大而有限，(9)已明确其量级；不作经济性或制造可行性承诺。

## 7. 实际反作用、端点交互及算符域

记$\Phi_0(\tau,x)=e^{i\tau H_\phi}\phi_\epsilon(x)e^{-i\tau H_\phi}$、$U(t)=e^{-itH}$。CCR给

$$
[\phi_\epsilon(y),\Phi_0(\tau,x)]=iD_\epsilon(\tau,x-y)I,
\quad |D_\epsilon(\tau,x-y)|\le|\tau|q_\epsilon,
\tag{23}
$$

因为$\|\sin(|D|\tau)/|D|\|\le|\tau|$。光滑紧支F的Fourier表示可逐项对中心CCR求交换，给$[F(\phi_\epsilon(y)/A),\Phi_0]=iD_\epsilon F'/A$。由此是有界交换子；无需假设F对两个任意非对易算符是operator-Lipschitz。

域的依据另行说明：三维中光滑紧支$g$满足$\widehat g/|k|\in L^2$及$\widehat g/\sqrt{|k|}\in L^2$，因此Segal测试场$\phi(g)$为$(H_\phi+1)^{1/2}$有界，$\Pi(g)$同理；自由演化测试也保这些条件。本轮$D(H_0)$实际输入／输出均在所用测试场域。Weyl共轭把另一测试场平移一个标量；F具有有界一、二阶导数，使V的该共轭可在范数中微分，得到上述有界交换子。在有限粒子光滑芯先用有界Weyl算符正则化Duhamel，再对Weyl参数在0求导，闭算符域及统一交换子界延拓到所用输入。等价地，这是随自由演化测试函数变化的有界交换子Duhamel。得到实际retarded公式及其有界差：

$$
U(t)^\dagger\phi_\epsilon(x)U(t)-\Phi_0(t,x)
=i\int_0^tU(s)^\dagger[V,\Phi_0(t-s,x)]U(s)ds,
\quad
\left\|\frac{\phi_{\epsilon,H}(t,x)-\Phi_0(t,x)}A\right\|
\le b_{back}(t):=\frac{\lambda LC_*q_\epsilon t^2}{2A^2}.
\tag{24}
$$

每个测试场域在该比较中受有界差控制。特别$x\in\operatorname{supp}\eta$及完整输入$\Psi$上有

$$
\|(\phi_{\epsilon,H}(t,x)/A-t)\Psi\|
\le v_\epsilon/A+b_{back}(t).
\tag{25}
$$

在$t=0,1$，$F(t)=0$。以同一个标量$t$作谱演算，再用(25)；等时的Heisenberg材料和源因子仍对易。因此实际端点均交互能满足

$$
\boxed{|\langle V(0)\rangle|,|\langle V(1)\rangle|
\le\lambda C_*L\left[\frac3{\epsilon A}+\frac{\lambda LC_*(8\epsilon^{-3})}{2A^2}\right]
<8.385846\times10^{-7}.}
\tag{26}
$$

实际交互没有被强行时间截断；源涨落使窗外仍可有交互尾。(26)独立来自retarded场响应，不能用(18)代替。

真实力须以测试函数描述。对实$g\in C_c^\infty$，$\Pi(g)=\int g\pi$，共同芯及闭形式延拓上

$$
\frac d{dt}\Pi_H(g)=\phi_H(\Delta g)
-\frac\lambda A\int\eta(x)(\rho_\epsilon*g)(x)
\bigl[F'(\phi_\epsilon(x)/A)J_\epsilon(x)\bigr]_H\,dx.
\tag{27}
$$

交互力范数不超过$\lambda LC_*\|\rho_\epsilon*g\|_\infty/A$；$\phi_H(\Delta g)$是自由波动力，不能删掉。有限粒子芯、光滑测试和有界交互交换子给积分式的域解释，不把点动量密度视为有界算符。模型显式保留这一相互作用算符力，不将来源替换为冻结参数。本特定码的某些平均力可能因对称性抵消，本轮不要求每个测试g或每个未知输入都有非零平均推力；反作用的算符、关联、实际后态及能源始终保留。

## 8. 同H的完整空间动量责任式

令$P=P_D+P_\phi$为材料全部粒子／空穴和标量场的共同物理空间平移生成元。固定约定$T(a)=e^{-ia\cdot P}$；场密度$W(x)=F(\phi_\epsilon(x)/A)J_\epsilon(x)$满足$T(a)W(x)T(a)^\dagger=W(x+a)$。核平移协变，故

$$
T(a)VT(a)^\dagger=\lambda\int\eta(y-a)W(y)dy,
\quad
\boxed{\frac d{dt}P_{j,H}(t)=-\lambda\int(\partial_j\eta)(x)W_H(t,x)dx.}
\tag{28}
$$

由于$\eta$在$L^1$中可微，右侧有界，范数不超过$16\lambda\,8q_\epsilon\|\partial_j\eta\|_1$；这一论证只对总P作平移，未把分别的$[P_\phi,V]$或材料场导数无证宣称为有界。符号由$\partial_a(TVT^\dagger)=-i[P,V]=-\lambda\int(\partial\eta)W$固定。

正自由能源给$|P_{D,j}|\le H_F$、$|P_{\phi,j}|\le H_\phi$，因此本轮$D(H_0)$输入／输出具有P域；有界交换子积分式控制实际动量变化。(28)明确把非守恒总场动量归给静态部署函数的外力。若没有动态定位支撑的相反反冲账，就不能把场能源守恒加此责任式称为闭合四动量。本轮没有构造支撑应力、实际CO/SR钟、boost生成元或共同Lorentz变换。

## 9. 复算与停止线

[check.py](check.py)使用已有Python和NumPy，默认只读比较[results.json](results.json)，仅`--write-results`登记结果。四组检查分别为：精确有理预算；bump核、谱投影和有限跳跃诊断；含负谱真空关联的完整CAR码及复相位校正；有界有限源模型的联合演化、真实源后态和共同能源账。后两组只检代数与机制，不具有连续标量场的CCR、真空或UV证书。

本例通过(18)、(20)、(22)、(26)、(28)共同准入。拒绝只优化A、ε、误差常数再算一轮；后续须连接实际缺失的钟、定位支撑、准备／提取仪器或共同变换中的一个实质接口。详见[selection_audit.md](selection_audit.md)、[NEXT.md](NEXT.md)。
