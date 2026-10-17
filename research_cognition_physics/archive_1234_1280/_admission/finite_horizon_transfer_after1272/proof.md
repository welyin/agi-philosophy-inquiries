# 正常相干初值到有限视界段：同一最小标量的自由能流接续

日期：2026-10-17。非编号解析接续，沿[来源基线](source_baseline.json)记录的范围开展。不新增认知公理、实际设备、正式研究轮或科学计数；CT仍未采用。单位 $c=\hbar=k_B=1$，度规号差 $(-+++)$。

## 1. 公开材料选择与本次净增

本稿选择固定四维 de Sitter 背景上的**最小耦合**实标量：

$$
ds^2=-F(r)dt^2+F(r)^{-1}dr^2+r^2d\Omega_2^2,\quad
F(r)=1-H^2r^2,\quad 0\le r<H^{-1},
\tag{1}
$$

$$
I_{\rm m}=-\frac12\int\sqrt{-g}
\left[(\nabla\Phi)^2+m^2\Phi^2\right]d^4x,\qquad m^2=2H^2.
\tag{2}
$$

$H$是本次固定背景参数；$m$是事前选定的常数材料参数。对度规变分时固定 $m$，不能随变分后的曲率把 $m^2$更新成 $R/6$。因此本稿与[上一静态来源桥](../static_reference_source_bridge_after1272/proof.md)中的共形作用是**不同物理完成**：同一固定背景上自由方程及 BD 二点函数相同，但物理 Hilbert 应力不同。本稿不继承其改进应力或物质 Wald 接触项。

参考取这份最小标量的 Bunch–Davies 状态 $\omega_{\rm BD}$；$\zeta=\partial_t$在中心 $r=0$单位归一，$t$是该中心测地线固有时。复用[Fröb–Much–Papadopoulos，式(66)、(70)—(78)](https://arxiv.org/pdf/2310.12185)的完整静态楔结果：对楔内正常局部 Weyl 相干态，其相对熵等于 $\beta_0=2\pi/H$乘同一最小标量的静态来源荷。这里不是重新证明或原创发现该成熟定理。

旧结果直接复用：

- [1004，§4](../../../archive_990_1008/research_note_1004.md)：相干应力增量、输入与来源不能删除；
- [318](../../../archive_301_341/research_note_318.md)：有限子段不能任意替代完整模 Hamiltonian；
- [旧 de Sitter 有限字典，§3、6](../thermal_mechanism_admission_after1270/de_sitter_finite_dictionary_after1270.md)：有限线性面积端点及其范围；其类光数据未认证为正常相干态；
- [1250](../../1250/proof.md)：另一材料的独立有限几何/来源比较，不把其环面、钟或读口移到本模型；
- [上一静态来源桥](../static_reference_source_bridge_after1272/proof.md)：正常紧支初值与整片模荷的接入。

本次项目内净增是把**同一正常初态的实际自由传播**明确接到一份有限时间球面通量和一份有限仿射视界段。它仍是成熟定理与显式解析解的非编号接续，不是已完成的自洽量子几何任务。

## 2. 同一物理方程的无势径向解

只取球对称部门，定义

$$
x=\frac{\operatorname{atanh}(Hr)}H,\qquad
\frac{dr}{dx}=F(r),\qquad \psi=r\varphi .
\tag{3}
$$

由 $(\Box-m^2)\varphi=0$直接得

$$
\psi_{tt}-\psi_{xx}
+F(r)\left(m^2+\frac{F'(r)}r\right)\psi=0 .
\tag{4}
$$

这里 $F'(r)=-2H^2r$，所以所选质量使势项恒为零。这一消去只在当前固定背景上成立，不自动适用于其反作用后的度规。

取实、非零 $f\in C_c^\infty([a,b])$，$0<a<b<\infty$。给定

$$
\psi(0,x)=f(x),\qquad \psi_t(0,x)=-f'(x).
\tag{5}
$$

满足中心正则性 $\psi(t,0)=0$的全静态半线解是

$$
\boxed{\psi(t,x)=f(x-t)-f(-x-t).}
\tag{6}
$$

对 $t\ge0$，第二项恒零，故

$$
\varphi(t,r)=\frac{f(x-t)}r .
\tag{7}
$$

不能把(7)无条件延伸到所有负时刻，否则当 $-t$落在初始支撑中时会破坏中心正则性。式(6)在中心为 $x$的奇函数，除以同样为奇函数且一阶非零的 $r(x)$后给光滑球对称场。

## 3. 紧支初值、正常态及输入边界

在 $t=0$，初值支撑于

$$
r\in\left[H^{-1}\tanh(Ha),\,H^{-1}\tanh(Hb)\right]
\Subset(0,H^{-1}).
\tag{8}
$$

场及其时间导数在中心和视界附近均为零；把它们在完整全局 Cauchy 面的其余部分置零，得到光滑紧支数据。线性 KG 演化与有限传播给正常的全局场解，不需新的 PDE 假设。

取 $0<\delta<a/2$，令 $\chi(t)$光滑地从过去的0过渡到未来的1，导数支撑于 $(-\delta,\delta)$。在这份有限因果管内，支撑仍有 $x\le b+\delta<\infty$，因而严格在 patch 内。定义

$$
P=\Box-m^2,\qquad
g=P(\chi\varphi)=[P,\chi]\varphi .
\tag{9}
$$

于是 $g\in C_c^\infty(\mathrm{patch})$。以 $E=E_{\rm ret}-E_{\rm adv}$为约定，
$E_{\rm ret}g=\chi\varphi$、$E_{\rm adv}g=(\chi-1)\varphi$，故 $Eg=\varphi$。采用使均值为 $Eg$的 Weyl 符号约定，$W(g)\Omega_{\rm BD}$是 BD 表示中的正常局部向量态，且保 Hadamard 短程结构。

这份初态可作为公开的有限能源输入。其全部能源和与其它系统的关联须保留；不能把纯相干初态称为免费热浴。依[已有目标范围，§2](../goal_scope_and_joint_delivery_after1250.md)，这里**不要求为未调用的初态输入再构造一个准备工厂**。若日后实际调用外源去实施(9)，该外源及其交换另属于新过程，不在本稿自动获得来源资格。

对同一固定背景、同一重整化处方，零均值参考与相干态相减后，

$$
\Delta\langle T_{\mu\nu}\rangle
=\partial_\mu\varphi\,\partial_\nu\varphi
-\frac12g_{\mu\nu}\left[(\nabla\varphi)^2+m^2\varphi^2\right].
\tag{10}
$$

以下使用这份完整最小 Hilbert 来源，而非仅保传播方程、不保质量项。背景几何的参考相减项不构成新增激发源；这不等于已认证变化背景后的绝对真空来源。

## 4. 初始静态荷：所有质量与梯度项共同相消

定义正静态能源

$$
Q=\int_{\Sigma_0}n^\mu\zeta^\nu
\Delta\langle T_{\mu\nu}\rangle\,d\Sigma .
\tag{11}
$$

将(5)代入(10)，用 $dr=Fdx$，得

$$
Q=2\pi\int_0^\infty
\left[
f'^2+\left(f'-\frac Frf\right)^2+m^2Ff^2
\right]dx .
\tag{12}
$$

令 $q(x)=F/r$。积分交叉项 $-q(f^2)'$并保边界，有

$$
Q=4\pi\int f'^2dx
+2\pi\int F\left(m^2+\frac{F'}r\right)f^2dx
-2\pi\left[\frac Frf^2\right]_{0}^{\infty}.
\tag{13}
$$

势项由所选材料消去；边界项由紧支初值消去。因此

$$
\boxed{Q=4\pi\int_a^b f'^2dx>0,\qquad
S_{\mathcal A(\mathrm{patch})}(\omega_g\Vert\omega_{\rm BD})
=\beta_0Q,\quad \beta_0=\frac{2\pi}H.}
\tag{14}
$$

质量能没有被丢弃；它参与了(13)的精确相消。相对熵等式引用的是完整静态楔的成熟定理，不把有限球自行当作相同区域代数。

## 5. 有限球面的真实自由通量

取有限 $L>b$及 $R=H^{-1}\tanh(HL)$。外向正 Killing 功率为

$$
P_R(t)=-4\pi R^2\,\Delta T^r{}_t(t,R).
\tag{15}
$$

因为 $\partial_t\varphi=-f'/r$、$\partial_r\varphi=f'/(Fr)-f/r^2$，有

$$
\boxed{P_R(t)=4\pi\left[
f'(L-t)^2-\frac{F(R)}R f(L-t)f'(L-t)
\right].}
\tag{16}
$$

全脉冲通过此球面的时间恰为

$$
t\in[L-b,L-a]\Subset(0,\infty).
\tag{17}
$$

将 $y=L-t$代入，$f(a)=f(b)=0$，所以

$$
\boxed{\int_{L-b}^{L-a}P_R(t)\,dt
=4\pi\int_a^b f'^2dy
-2\pi\frac{F(R)}R[f^2]_a^b
=Q.}
\tag{18}
$$

瞬时的第二项不能删除；本等式只在完整脉冲两端为零时使其积分消失。任意截短窗口须另保该端点项。式(18)是同一完整物理来源的有限时间流出量，不是从全荷重新定义功率。

本稿没有安装球面接收器。有限 $L$使通过时刻有限；若要把它变成实际末端记录，仍须指定合格读口及其来源。因果上也不能把视界事件直接称为中心在同一有限窗口收到的事件。

## 6. 未来视界：正确仿射方向与有限支撑

令出射坐标 $v=t-x$。度规成为

$$
ds^2=-Fdv^2-2dv\,dr+r^2d\Omega_2^2 .
\tag{19}
$$

未来宇宙学视界 $r=H^{-1}$上有

$$
\Gamma^v{}_{vv}=-\frac{F'}2=H,\qquad
\lambda=\frac{e^{Hv}}H,\qquad
k=\partial_\lambda,\qquad
\zeta=\partial_v=H\lambda k .
\tag{20}
$$

这是未来向仿射参数。也可从嵌入直接核验：
$X^0=X^j=\lim_{x\to\infty}(N/H)\sinh(Ht)=e^{H(t-x)}/H=\lambda$；
该未来视界生成线是嵌入空间的直类光线。取 $-e^{-Hv}/H$不能满足这份出射图的仿射关系。

同一解在视界上的光滑迹为

$$
\varphi_{\mathcal H^+}(v)=Hf(-v),\qquad
v\in[-b,-a],
\tag{21}
$$

区间外为零。于是

$$
\frac{d\varphi}{d\lambda}
=-\frac{f'(-v)}{\lambda},\qquad
\Delta T_{\lambda\lambda}
=\frac{f'(-v)^2}{\lambda^2}.
\tag{22}
$$

这里的质量项因 $g_{\lambda\lambda}=0$不进入该收缩；它仍属于完整应力(10)。背景面积元为 $dA_0=H^{-2}d\Omega_2$，脉冲的仿射支撑为

$$
\lambda\in[\lambda_a,\lambda_b]
=\left[\frac{e^{-Hb}}H,\frac{e^{-Ha}}H\right]
\Subset(0,\infty).
\tag{23}
$$

直接改变积分变量，$d\lambda/\lambda=Hdv$，得到

$$
\int_{\lambda_a}^{\lambda_b}
\Delta T_{\mu\nu}\zeta^\mu k^\nu\,d\lambda\,dA_0
=H\int_{\lambda_a}^{\lambda_b}
\lambda\Delta T_{\lambda\lambda}\,d\lambda\,dA_0=Q,
\tag{24}
$$

以及同一参考/时间归一的身份

$$
\boxed{\beta_0Q
=2\pi\int_{\lambda_a}^{\lambda_b}
\lambda\,\Delta T_{\lambda\lambda}\,d\lambda\,dA_0
=S_{\mathcal A(\mathrm{patch})}(\omega_g\Vert\omega_{\rm BD}).}
\tag{25}
$$

这是同一正常态的真实自由传播，不是先指定一段独立类光波形再假定它由该态产生。有限仿射支撑不等于有限中心固有时收信：到达视界在静态坐标中是 $t\to\infty$；有限球(17)与视界段(23)是不同的任务边界。

相干位移提供的是具有确定相位的、纯功样的激发资源。KMS参考下的加权能流(25)不能仅凭名称改成已经认证的热库交换 $\delta Q$；本稿未实现热接触、热指针或独立熵计。

## 7. 一份完全明确的光滑波包与正下界

取 $H=1,a=1,b=2$，并置

$$
f(x)=
\begin{cases}
\exp\!\left(1-\dfrac1{1-(2x-3)^2}\right),&1<x<2,\\
0,&\text{其余}.
\end{cases}
\tag{26}
$$

它在两个端点各阶导数为零，$f(5/4)=e^{-1/3}$。Cauchy–Schwarz在 $[1,5/4]$给

$$
\int_1^{5/4}f'^2dx
\ge\frac{(f(5/4)-f(1))^2}{1/4}
=4e^{-2/3}.
\tag{27}
$$

因此可保守登记

$$
\boxed{Q\ge16\pi e^{-2/3}
\ge\pi e^{-2/3}>0.}
\tag{28}
$$

例如固定 $L=3$时，完整球面脉冲在 $t\in[1,2]$经过；未来视界脉冲位于 $\lambda\in[e^{-2},e^{-1}]$。这些均是精确有限区间。无需以未认证数值求积来证明非零性。

## 8. 几何接法只引用旧条件结论，不新增其证书

若**另行**采用同一源的线性 Einstein/Raychaudhuri响应，并给定合法的类光面/仿射参数识别、独立面积初值与端点，则旧有限端点工具给

$$
\delta A''(\lambda)
=-8\pi G\int\Delta T_{\lambda\lambda}\,dA_0,\qquad
2\pi\int\lambda\Delta T_{\lambda\lambda}\,d\lambda\,dA_0
=\frac{[\delta A-\lambda\delta A']_{\lambda_a}^{\lambda_b}}{4G}.
\tag{29}
$$

式(29)仅引用[既有 de Sitter 有限字典§3](../thermal_mechanism_admission_after1270/de_sitter_finite_dictionary_after1270.md)在本仿射取向下的条件工具。这里没有把右端定义成面积测量，也没有证明一份包含该正常量子态、参考真空变化、全部约束初值和实际几何读口的共同反作用解。固定 $m$在新度规上的传播势一般不再为零；(7)不能原样充当自洽回授后的解。

没有获得有限振幅的剪切/膨胀非线性、真空极化及几何变化总误差。未调用的这些量也不被升级为所有有限共同模型的通用门槛。

## 9. 准入与停止线

- 新闭合的是正常态到有限球面/视界能源流的**同对象解析连接**；不是新引力生成定理。
- 本稿公开采用新最小材料，不能把1270/1272的实际端口、共动热参考、源完成及误差预算移入；也不能把其热流反向认作本稿的视界流。
- 初态有限能源输入合法，不增加准备工厂要求；但实际额外准备、交互或读取若被调用，其来源必须加入同一过程。
- 将 $f$乘幅度 $\eta$时，荷、应力和相对熵均乘 $\eta^2$。$\eta$和$-\eta$不是317的双向非零一阶饱和方向；这不否定有限/二阶选择，也不声称本稿已完成它。
- 当前可作为非编号成熟接续保留。只有另外取得同一对象的独立有限几何/来源/实际任务结论，才可考虑正式新轮；不靠重复(29)预约轮次，不补造热装置。
