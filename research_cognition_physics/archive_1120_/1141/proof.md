# 1141证明：同Ho相对区域源能与真实量子尾

2026-10-12。沿1140后的方向决定，停止为该候选全网逐门配置控制器。本轮复用原1129/1133/1137的正常场、原参数与b=0任务部门，在比较动力学H_o上建立一项新的有限证书：**同一过去的区域源能减去自由源基线后，其真实反作用贡献有不随A增长的主界。**

新增认知公理0、科学校准0。不改A、ε、正常Y、原码、模式、脉冲或保持材料。该中间结果不是实际BT动力学H的迁移，不是正功下界，也不是已造出的能源仪器。总研究目标不缩小。

## 1. 原对象及本轮问题

沿用[1129](../1129/proof.md)、[1133](../1133/proof.md)、[1137](../1137/proof.md)：

$$
\epsilon=10^{-5},\quad A=10^{34},\quad
\delta_c=10^{-6},\quad t_c=5/4,\quad
\lambda=\pi/2<11/7,\quad\kappa=\pi/4<11/14.
$$

原正常源a的Cauchy均值为φ_cl(0)=0、π_cl(0)=Aβ_0，β_0在B_7为1、支撑于B_8。保原完整正常Dirac海编码J、未知Q与任意被动R；实际正常Y、完整旧D、CP、b=0的Z及正常方向O均保留。

本轮唯一演化生成元是1137式(2)的

$$
H_o=H_F+H_\phi+\sqrt{(M+H_D+1)^2+P_Y^2}
+m_C+H_Z+m_O+H_O+V_{1129}(Y,b=0)+V_{\rm rec}(Y,O).
\tag{1}
$$

它与实际BT生成元H不同。H_o保正常量子源，不是经典脉冲过程；所有接收与保持项从t=0就在式(1)中。短窗内经典S、F_2、F_3尚未启动，不等于真实量子尾为零。

输入用纯化表示为原准备的归一Ψ；源因子是固定a，其他正常部门可以含原Q/R相关。下面的界统一这些任务输入，不改成海真空校准态。

## 2. 不能把非零功当成交接必然要求

理想均匀平台交接的全空间比较生成元为

$$
H_F+\lambda F(t)Q_x,\qquad
Q_x=q(I_{\rm space,spin}\otimes\sigma_x^{AB}\otimes I_{\rm rail}).
$$

因A/B有同一h，Q_x与H_F对易。原A编码及B正常海没有AB相干；理想全空间味旋转保持局部σ_x双线性，故其均值为0。于是原逻辑交接可以发生，而理想材料全局能量变化和经典系数下的平均源力均为0。

这只判断理想比较结构。实际η边缘、正常海响应、量子源相关、动态Y及真实保持/接收尾未被排除。尤其〈σ_x〉=0不能推出实际相关量〈F'(φ_ε/A)J_ε〉=0。也不能把原巨大自由源流称为交互做功。

因此本轮选取**净反作用的上界**，不先索取不适合该简并任务的正功下界。

## 3. 独立源力及其全空间列范数

写原φ_ε(x)=φ(ρ_{ε,x})，K_ε是与ρ_ε卷积，||ρ_ε||_1=1。原J_ε为完整旋量、A/B及双rail的偶CAR密度，||J_ε(x)||≤8q_ε，q_ε=||ρ_ε||_2²≤8ε^{-3}。

从原相互作用对φ的实谱平移求导，定义

$$
\begin{aligned}
\mathcal B_x={}&\lambda\eta(x-Y)
 F'(\phi_\epsilon(x)/A)J_\epsilon(x)\\
&+\chi(x-Y)\left[
-S'(\phi_\epsilon(x)/A)Q_0(Y,O)
+\lambda F_2'(\phi_\epsilon(x)/A)N_C(Y,O)X_P
+\kappa F_3'(\phi_\epsilon(x)/A)Y_g(Y,O)
\right],\\
\mathcal R(z)={}&A^{-1}\int\rho_\epsilon(x-z)\mathcal B_x\,d^3x .
\end{aligned}
\tag{2}
$$

这里ℛ=δV/δφ；波动方程的力是−ℛ。材料与标量因子在同一时刻对易。上述导数及核由同一H_o独立决定，未把待证明的区域流反向定义为总H交换子。

原界为

$$
\|Q_0(Y,O)\|<42,\quad\|N_C\|=1,\quad\|Y_g\|=1,
\quad\|F'\|_\infty,\|F_j'\|_\infty\le512,\quad
\|S'\|_\infty\le8/\delta_c.
$$

它们对所有Y、O和完整场/寄存器统一成立。令

$$
a_*=\frac{11}{7}\,512\,64\epsilon^{-3},\qquad
b_*=42\frac8{\delta_c}+
\left(\frac{11}{7}+\frac{11}{14}\right)512.
$$

先固定Y=y（及需要时的O纤维），对任意该纤维向量v作L²_x三角估计：

$$
\|\mathcal B_\cdot(y)v\|_{L^2_x}
\le\bigl(a_*\|\eta\|_2+b_*\|\chi\|_2\bigr)\|v\|.
\tag{3}
$$

η≤1支撑于B_5，||η||_2<sqrt(524)<23；χ=ρ为原尺度1归一核，||χ||_2²≤8，所以||χ||_2<3。再用Hilbert空间值Young不等式及||ρ_ε||_1=1，最后对Y/O直积分，得到

$$
\|\mathcal R(\cdot)\Psi\|_{L^2_z(\mathscr H)}
\le\frac D A\|\Psi\|,
\quad
D=23a_*+3b_*
=\frac{8290304000007056025344}{7}
<1.18433\times10^{21}<1.2\times10^{21}.
\tag{4}
$$

**顺序是先逐Y纤维L²估计，再直积分；不能先对Y取点态上确界后对全空间积分。** 所有正常Y尾均保留。该列界是整个b=0任务部门上的算符界，不依赖未知Q或被动R。

设U(t)=e^{-itH_o}。由于ℛ_H(s)=(I_{L²}⊗U(s)†)ℛU(s)，式(4)原样控制真实Heisenberg源力列；没有沿参考轨道替换真实反作用。

## 4. 独立退迟方程与联合波动列界

原1129的光滑测试场/Weyl论证逐纤维适用于式(2)。H_o的其余自由项与标量场对易，因此

$$
(\partial_t^2-\Delta)\phi_H=-\mathcal R_H.
$$

对测试g，独立退迟式为

$$
\phi_H(g,t)=\phi_{\rm free}(g,t)
-\int_0^t\mathcal R_H(s)
 \left(\frac{\sin(|D|(t-s))}{|D|}g\right)ds .
\tag{5}
$$

π的公式对时间求导。没有将实际BT H代入这条方程。

本轮能源读取采用原分辨率K_ε。记

$$
Z(t)=\bigl(\pi_{K,H}(t),\nabla\phi_{K,H}(t)\bigr),\quad
Z_0(t)=\bigl(\pi_{K,\rm free}(t),\nabla\phi_{K,\rm free}(t)\bigr),
\quad\delta Z=Z-Z_0 .
$$

联合δZ由式(5)中的cos和∇sin/|D|核组成。对每一动量，二者平方和为1；K_ε在L²上收缩。式(4)与Bochner积分三角不等式给

$$
\|\delta Z(t)\Psi\|_{L^2(\mathbb R^3;\mathscr H^4)}
\le \frac{tD}{A}\|\Psi\|,\qquad t\ge0.
\tag{6}
$$

这里是一个联合列界，没有把π与三个梯度分别粗加。卷积力已有一个ρ_ε；能源测试又带K_ε，二者均在公式中，未漏掉第二次平滑。

## 5. 独立区域源能源、闭形式与共同过去

取原R为归一bump累积，固定实验室区域权

$$
w(x)=\eta(x)\,[1-R(x_1+1/2)] .
\tag{7}
$$

0≤w≤1、w∈C_c∞，支撑在B_5，体积V_w<524。它的过渡界面穿过原交互区。w是固定lab权，不随Y移动，因此不遗漏移动边界功；Gaussian Y仍按原H_o动态演化。

独立定义正常序的有限分辨率自由源能源

$$
E_{\epsilon,w}
=\frac12\int w(x):
\left[\pi_K(x)^2+|\nabla\phi_K(x)|^2\right]:d^3x .
\tag{8}
$$

这是四个平滑Segal场列的平方范数闭二次形式，减去有限真空常数。它不是裸点时Wick区域能源，也不是从ℛ或总H的交换子定义出来。

原输入属于H_o的共同自伴域。H_o是正自由生成元的有界扰动，其传播保持这一域，故实际输出的H_φ域及式(8)的形式域均受控。对平滑测试，Π、∂φ及所需额外空间导数都相对(H_φ+1)^{1/2}有界。直接以列的闭图/闭形式定义式(8)，不要求一般有界函数自动保持高阶域。

测试退迟式先在共同平滑芯成立；式(4)给强可测且一致有界的L²列，从而式(5)可用强Bochner积分延拓。实际场列与自由场列的差在全Hilbert空间有界。以下只取合法向量上的二次形式期望，不冒称差场列的逐点乘积是有界能源算符。

定义相同初态、相同过去的相对区域源能

$$
\Delta E_{\epsilon,w}(t)
=\langle U(t)\Psi,E_{\epsilon,w}U(t)\Psi\rangle
-\langle e^{-itH_\phi}\Psi,
 E_{\epsilon,w}e^{-itH_\phi}\Psi\rangle .
\tag{9}
$$

第二项只作用来源因子；其他自由项对源观测无影响。它是已知原相干源的自由基线，没有重新准备实际输出，也不删掉来源后态和相关性。

由式(8)逐项展开，**不交换实际与自由场算符的顺序**，得到精确二次形式身份

$$
\Delta E_{\epsilon,w}(t)=
{\rm Re}\langle\sqrt w Z_0(t)\Psi,\sqrt w\delta Z(t)\Psi\rangle
+\frac12\|\sqrt w\delta Z(t)\Psi\|^2.
\tag{10}
$$

正常序常数及巨大的共同自由A²项已在这一步精确相消。这一身份保留交叉关联，不用小态差乘以巨大能源估计来代替。

## 6. 固定相干源的列常数与有限预算

采用酉Fourier规范，令

$$
\nu_\epsilon=\int |k|\,|\widehat\rho_\epsilon(k)|^2d^3k .
$$

它是π_K和三个∂φ_K真空方差之和。原核满足||ρ||_2²≤8、|∇ρ|<16、|B_1|<88/21，故

$$
\nu_\epsilon
\le\|\rho_\epsilon\|_2\|\nabla\rho_\epsilon\|_2
<100\epsilon^{-4}.
\tag{11}
$$

原||β_0||_2²<2146<47²。自由经典源列为

$$
A\bigl(K_\epsilon\cos(|D|t)\beta_0,\,
K_\epsilon\nabla\sin(|D|t)|D|^{-1}\beta_0\bigr),
$$

其联合全空间范数不超过A||β_0||_2。自由相干源的均值与真空涨落正交，故

$$
\|\sqrt w Z_0(t)\Psi\|
\le\sqrt{A^2\|\beta_0\|_2^2+V_w\nu_\epsilon}
\le C_w:=47A+230\epsilon^{-2}.
\tag{12}
$$

这些数对原Q/R及正常Y/O/Z/CP输入统一；固定来源a是本合同的原准备，没有改用海校准。

将式(6)、(10)、(12)合并：

$$
\boxed{
|\Delta E_{\epsilon,w}(t)|
\le C_w\frac{tD}{A}+\frac12\left(\frac{tD}{A}\right)^2 .
}
\tag{13}
$$

主项47tD独立于A；剩余项随1/A及1/A²衰减。没有调A来获得本结论。

固定t≤5/4时，精确有理预算为

$$
\frac{tD}{A}\le1.480412\times10^{-13},\qquad
C_w=470000000000000000000002300000000000,
$$

$$
|\Delta E_{\epsilon,w}(t)|
<6.958\times10^{22}<7.1\times10^{22}<10^{23}.
\tag{14}
$$

选择统一证书10^{23}。这与自由源O(A²)约10^{68}的绝对量不同：改进来自共同过去相消及1/A的真实退迟力，不来自改模型参数。它是上界，未给正功、非零可读反差或相对小于某项实际交互功的证书。

## 7. 独立局部流、交换功率与弱形式

由平滑场本身定义

$$
j_\epsilon=-\frac12:\{\pi_K,\nabla\phi_K\}:,\qquad
\Phi_{\epsilon,w}=-\int\nabla w\cdot j_\epsilon\,d^3x ,
$$

$$
P_{\epsilon,w}
=-\frac12\int w(x)\{\pi_{K,H}(t,x),
  (K_\epsilon\mathcal R_H)(t,x)\}\,d^3x .
\tag{15}
$$

Φ是向外流的权重约定；P是来源自由能受到的交互功率。反对易积不能改成未经证明可交换的普通乘积。K_εℛ包含原保持、记录和继续的全部量子尾。

式(15)的功率期望严格按二次形式定义：

$$
\langle P_{\epsilon,w}(t)\rangle_\Psi
=-{\rm Re}\left\langle
\sqrt w\,\pi_{K,H}(t)\Psi,\,
\sqrt w\,K_\epsilon\mathcal R_H(t)\Psi
\right\rangle_{L^2_x(\mathscr H)}.
$$

反对易符号是这一实内积的记法，不主张$\pi_{K,H}K_\epsilon\mathcal R_H\Psi$等算符乘积都在原Ψ上有定义。流的期望也按对应平滑场列的实内积理解。

从∂tφ_K=π_K、∂tπ_K=Δφ_K−K_εℛ直接作乘积法则与空间分部积分，得

$$
\frac d{dt}\langle E_{\epsilon,w,H}(t)\rangle
=-\langle\Phi_{\epsilon,w,H}(t)\rangle+
\langle P_{\epsilon,w}(t)\rangle .
\tag{16}
$$

严格理解为时间弱形式：对任意紧支光滑τ，

$$
-\int\tau'(t)\langle E_{\epsilon,w,H}(t)\rangle dt
=\int\tau(t)\left[-\langle\Phi_{\epsilon,w,H}(t)\rangle+
\langle P_{\epsilon,w}(t)\rangle\right]dt .
$$

空间平滑给列所需的有限空间导数。具体地，ρ_ε及其各阶空间导数在L¹，故$|k|m_\epsilon(k)$等有限阶乘子有界；$\mathcal R_H(s)$由有界列共轭得到，在每个向量上强连续。退迟Bochner积分因此给δZ在紧时间窗内的局部强C¹列，求一次时间导数只需这些有界空间乘子及端点力列，不需对力列求时间导数。自由平滑场测试在原$D(H_\phi^{1/2})$输入上同样局部强C¹。于是局部列内积可以作乘积法则和空间分部积分，得到所列时间弱形式；无需宣称H_o保持任意高阶幂域。真空减项为时间常数，自由真空流的减项为0。

自由基线满足同式而P=0，所以同一过去的净账为

$$
\Delta E_{\epsilon,w}(t)
=-\int_0^t\Delta\langle\Phi_{\epsilon,w}(s)\rangle ds
+\int_0^t\langle P_{\epsilon,w}(s)\rangle ds .
\tag{17}
$$

同一列界还直接分别控制功和流差，不需新假设。记$q=D/A$、$T\le5/4$、$C=47A+230\epsilon^{-2}$。由于$\|\sqrt w\,\pi_{K,H}(t)\Psi\|\le C+qt$、$\|\sqrt w\,K_\epsilon\mathcal R_H(t)\Psi\|\le q$，

$$
\int_0^T|\langle P_{\epsilon,w}(t)\rangle|dt
\le CqT+\frac12q^2T^2<7.1\times10^{22}<10^{23}.
\tag{17a}
$$

这里控制的是功率期望绝对值的时间积分，因此也控制带符号净功；不假设它为正或非零。

原η的径向导数及R的导数都不超过8，所以式(7)的$G=\|\nabla w\|_\infty\le16$。在$B_5$上令$Z=(\pi_K,\nabla\phi_K)$，流形式的矩阵具有两个非对角块$\nabla w/2$，其范数不超过G/2。对实际列$Z_0+\delta Z$与自由列$Z_0$展开，式(6)及在整个$B_5$上同样成立的自由列界C给

$$
|\Delta\langle\Phi_{\epsilon,w}(t)\rangle|
\le G\left(Cqt+\frac12q^2t^2\right),
$$

$$
\int_0^T|\Delta\langle\Phi_{\epsilon,w}(t)\rangle|dt
\le16\left(\frac{CqT^2}{2}+\frac{q^2T^3}{6}\right)
<7.0\times10^{23}<10^{24}.
\tag{17b}
$$

这些是同H_o、同一过去、同原准备的分开预算，未将巨大自由通量算成交互功。流由独立波动能源乘积定义，力由原V的源导数定义，二者在式(16)相遇。界不证明正功、非零可读裕量或实际BT迁移。

## 8. 未解析源能与未完成的身份

全空间的分辨源能由Fourier直接给

$$
E_{\epsilon,\rm all}
=d\Gamma\bigl(|k|\,|m_\epsilon(k)|^2\bigr),\qquad
m_\epsilon=(2\pi)^{3/2}\widehat\rho_\epsilon .
$$

因为|m_ε|≤1，完整物理来源能源还含

$$
H_\phi=E_{\epsilon,\rm all}+E_{\rm rem},\qquad
E_{\rm rem}=d\Gamma\bigl(|k|(1-|m_\epsilon(k)|^2)\bigr)\ge0.
\tag{18}
$$

本轮保留E_rem责任，不删除它，也没有为它构造完整局部流。H_o的全局能源仍含完整H_φ、材料、支撑、O/Z、D/CP及V。式(17)因此是一份有限分辨率来源子账，不是整个交互区的完整局部能源闭合。

实际BT动力学$H=\sqrt{\mathsf M^2+K^2}$一般不直接满足式(5)。小的H–H_o态误差不能乘以无界E便得到本证书；1138外壳粗图预算也不自动解决该相对量。下一步必须使用同原源测试、同一过去和同一实际物理字典，独立控制BT差场列或相对功率。

式(8)定义有效能源量并证明合法域；没有制造终末能源探测器、免费授予其全部仪器、或证明读取后仍完成原任务。

## 9. 检查、净增量与停止点

[check.py](check.py)默认只读比对[results.json](results.json)，只有显式--write-results写本轮结果。四组为：

1. 原来源资产身份与范围冻结；
2. 精确有理Young常数、真空列常数和最终10^{23}预算；
3. 有限周期空间、量子支撑纤维及非交换矩阵力的Young/退迟列诊断；
4. 同一有限初态的非交换二次形式相消和Cauchy界诊断。

有限模型只检查代数和列界；没有连续Dirac海、原Gaussian支撑、CCR、17m模式或实际BT迁移的模拟。解析式(3)—(18)承担连续结论。未运行图像检验。

净增量是同H_o、同原准备、全部真实量子尾下的A⁰相对区域源能证书，并给出独立交换功与流差的分开预算。理想简并任务零功不是失败判据。下一步见[NEXT.md](NEXT.md)：直接审查该相对量能否接到实际H；不改回全网施工，也不以完成全部应力或UV极限为门槛。
