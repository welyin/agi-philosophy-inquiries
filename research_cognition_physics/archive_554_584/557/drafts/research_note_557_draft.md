# 第557轮研究：共同规范链路、实际读数与测量注能

日期：2026-09-30。接续[556](../../research_note_556.md)，本文件为终审前草稿；[代码](../joint_gauge_link_reference.py)、[结果](../joint_gauge_link_reference_results.json)、[预审](link_readout_review.txt)已保存。未冻结前不增加正式计数。

## 0. 问题、增量与范围

本轮在一个明确的两胞元候选中，证明同一正边势F能同时承担相邻物质作用、规范不变读数和可计算的物质／电场反作用；一个平均源能量E₁预算即控制该配置读数的二阶矩及非选择测量注能。

这是一处条件性共同实现，仍输入有限图、SU(2)、链路电动能、中性指针及理想测量耦合。没有证明这些输入来自认知公理，也没有产生空间维数、连续场论或引力。

## 1. 历史与成熟工具

- [360](../../../archive_342_369/research_note_360.md)、[362](../../../archive_342_369/research_note_362.md)、[365](../../../archive_342_369/research_note_365.md)的物理区域、Gauss部门和拼接边界直接继承，不重报“约束空间不是普通区域张量积”。
- [544](../../../archive_531_553/research_note_544.md)和[556](../../research_note_556.md)的同一正四次物质势继续使用，不重新拟合参数。本轮补一种相邻协变差，不等于[548](../../../archive_531_553/research_note_548.md)全部导数参考已实现。
- 群链路、局域作用及电动能属于成熟格点规范工具，参见[Zohar、Cirac、Reznik](https://arxiv.org/abs/1303.5040)第II、IV节。本轮采用其对象类型，不声称复现其冷原子装置或推出规范群。
- von Neumann指针耦合也是成熟工具，参见[Sokolovski](https://arxiv.org/abs/2009.12568)附录（第XIII节）。本轮增量是将它接到同一边势、局域Gauss部门和同一来源预算，并核清仪器失配反例。

认知动机是交互和记录须共用内部关系；具体群、图及连续配置空间是额外建模输入。Higgs四个实分量是内部表示，不能改称时空坐标。

## 2. 同一个源Hamiltonian与物理部门

两端i=1,2分别取Xᵢ∈ℝ⁴、sᵢ∈ℝ，rᵢ=|Xᵢ|，动态链路U∈SU(2)；R(U)为基本双重态的实四维正交表示。源空间为L²(ℝ¹⁰×SU(2))，Haar测度归一。令v、κ、b>0，a=ℏ²/(2v)，L严格正定，u=L⁻¹C，ℓ=λ_min(L)>0；继承u_h、u_s>0。以下正闭二次型定义Friedrichs算符：

$$
\begin{aligned}
H&=-a\sum_{i=1}^2(\Delta_{X_i}+\partial_{s_i}^2)
+b(-\Delta_G)+W_1+W_2+F,\\
W_i&=\frac v4
\begin{pmatrix}r_i^2-u_h\\s_i^2-u_s\end{pmatrix}^{T}
L\begin{pmatrix}r_i^2-u_h\\s_i^2-u_s\end{pmatrix},\qquad
F=\frac{\kappa v}{2}|X_1-R(U)X_2|^2,\qquad
-\Delta_G\chi_j=j(j+1)\chi_j .
\end{aligned}
\tag{1}
$$

这里只给定两个顶点一条边，无磁面项、无费米子；v并非已从几何导出的物理体积，b包含选定单位下的电动能系数。两端局域作用及物理投影是

$$
(X_1,X_2,U)\mapsto
(R(g_1)X_1,R(g_2)X_2,g_1Ug_2^{-1}),\qquad
P_{\rm phys}=\int dg_1dg_2\,\mathcal U(g_1,g_2),\qquad
[H,P_{\rm phys}]=[F,P_{\rm phys}]=0 .
\tag{2}
$$

s不变。正交性和双不变Laplacian给形式不变及算符强对易。物理态须支持在P_phys像中，不只是密度矩阵与群作用对易。

## 3. 平均能量直接控制边读口

令E₁=TrρH<∞，允许物质—链路相关、非Gaussian以及与任意旧参考纠缠。正二次型给〈F〉、〈W₁+W₂〉≤E₁。由z²≤2(z−u_h)²+2u_h²及平方和界，

$$
r_1^4+r_2^4\le4u_h^2+\frac{8(W_1+W_2)}{v\ell},\qquad
F^2\le2\kappa^2v^2(r_1^4+r_2^4),\qquad
\langle F^2\rangle\le B_F:=
8\kappa^2v^2u_h^2+\frac{16\kappa^2v}{\ell}E_1 .
\tag{3}
$$

这些是配置函数逐点界，没有使用错误的非对易推理“H≥F所以H²≥F²”。E₂不是这个纯配置读口的独立要求；556完整动量菜单的E₂条件不能一并取消。

## 4. 群元细读：正确概率不保证物理后态

考虑热核分辨仪器，t>0是分辨参数，不是H演化时间。共轭角χ∈[0,π]，类函数Haar积分为(2/π)∫sin²χdχ。Casimir约定给

$$
k_t(\chi)=\sum_{n=1}^{\infty}n e^{-t(n^2-1)/4}
\frac{\sin(n\chi)}{\sin\chi},\qquad
\int k_t\,dU=1,\qquad
\int k_t(U^{-1}\widehat U)R(\widehat U)d\widehat U
=e^{-3t/4}R(U)=:\eta_tR(U).
\tag{4}
$$

紧群热核严格正；归一和基本表示衰减可由角色正交性检验。两端位置再配独立中心Gaussian噪声，每分量方差t_q。全部效应是配置乘法；其乘积是环境空间上的正归一POVM。输出

$$
Y_{\rm fine}=\frac{\kappa v}{2}
\left(|\widehat X_1|^2+|\widehat X_2|^2-8t_q
-\frac{2}{\eta_t}\widehat X_1\cdot R(\widehat U)\widehat X_2\right),
\qquad \mathbb EY_{\rm fine}=\langle F\rangle .
\tag{5}
$$

两个范数各偏4t_q，交叉项由独立中心噪声及式(4)恢复原c=X₁·R(U)X₂。这是联合配置谱的条件积分，不要求Wigner正性，对任意相关源成立。输出整体规范不变，因此其粗粒化效应与群对易，压缩到物理部门仍为POVM。

单独考虑链路平方根仪器K_Uhat(U)=√k_t(U⁻¹Uhat)。取两端各自径向Schwartz物质态和链路常函数1，这是正规物理源且完整H的有限次矩均有限。结果后的链路函数不再常数；两端物质仍各自为singlet，所以联合Gauss投影只保常链路。每个结果及丢弃标签后的物理权重为

$$
p_{\rm phys}(t)=
\left(\int\sqrt{k_t(U)}\,dU\right)^2<1,\qquad 0<t<\infty .
\tag{6}
$$

严格性来自Cauchy–Schwarz和η_t>0保证的非恒定核。t=0.5时p_phys=0.44144936236508947。该比例仅对这份链路仪器和源成立，不是端点也被测后的通用泄漏率。密度仍可协变，物理支持却已丢失；这不排除更大的带荷探针实现或其它不变仪器。

该链路仪器的非选择源注能为

$$
D_t=\int|\nabla_G\sqrt{k_t}|^2dU
=\frac1{16}\frac2\pi\int_0^\pi
\sin^2\chi\,\frac{[k'_t(\chi)]^2}{k_t(\chi)}d\chi,\qquad
\Delta E_{\rm link\ instrument}=bD_t .
\tag{7}
$$

本约定Δ_G=(∂χ²+2cotχ∂χ)/4，连同平方根导数给1/16。乘法势和物质动能保持不变，群动能交叉项积分为零。D_0.5=0.5682053903670031，不包含端点仪器能量，也不使不保物理部门的操作合法。本轮不优化小t渐近。

## 5. 直接不变读口逐分支保Gauss

输入独立中性指针(Q,P)，[Q,P]=iℏ，实中心最小Gaussian，位置方差σ_q²。给定g>0，ν²=σ_q²/g²>0，理想测量酉及缩放输出y的条件算符为

$$
V_g=e^{-igF\otimes P/\hbar},\qquad
K_y=(2\pi\nu^2)^{-1/4}e^{-(y-F)^2/(4\nu^2)},\qquad
\int K_y^\dagger K_y\,dy=I,\qquad
[K_y,P_{\rm phys}]=0 .
\tag{8}
$$

F与P作用不同因子，指数由联合谱定义。K_y是F的函数，每个分支保物理支持，任意被动参考作恒等扩展。未知原态通常改变；f、f′间离对角因子为exp[−(f−f′)²/(8ν²)]，并非恒等操作。

$$
\mathbb EY=\langle F\rangle,\qquad
\mathbb EY^2=\langle F^2\rangle+\nu^2\le B_F+\nu^2,\qquad
\operatorname{Var}(Y)=\operatorname{Var}_{\rho}(F)+\nu^2 .
\tag{9}
$$

源涨落没有扣除，Var(Y)不能全叫仪器误差。明确的V_g不等于原源H已自主生成探针、选择或开关。

## 6. 同一边势控制物质与电场反作用

群梯度采用实反对称生成元J_A，−ΣJ_A²=3I/4。J_A R(U)X₂为三个互相正交、长度r₂/2的切向量。直接微分得

$$
\sum_{i=1}^2|\nabla_{X_i}F|^2=4\kappa v F,\qquad
|\nabla_GF|^2=\frac{(\kappa v)^2}{4}(r_1^2r_2^2-c^2),\qquad
G_F=a\sum_i|\nabla_{X_i}F|^2+b|\nabla_GF|^2 .
\tag{10}
$$

F无s依赖。对ρ′=∫K_yρK_y dy，非选择源能量精确满足

$$
\operatorname{Tr}(\rho'H)-E_1
=\frac{\langle G_F\rangle}{4\nu^2}
=\frac{a\kappa v}{\nu^2}\langle F\rangle
+\frac{b(\kappa v)^2}{16\nu^2}
\langle r_1^2r_2^2-c^2\rangle .
\tag{11}
$$

在形式核上，∇K_y=K_y(y−F)∇F/(2ν²)。所以∫K_y∇K_y dy=0、∫|∇K_y|²dy=|∇F|²/(4ν²)。展开∇(K_yψ)的平方，交叉项为零；势因对易及归一保持平均；两个动能分别给式(11)两项。再用r₁²r₂²≤(r₁⁴+r₂⁴)/2，

$$
\langle G_F\rangle\le B_G:=
4a\kappa vE_1+\frac{b\kappa^2v^2u_h^2}{2}
+\frac{b\kappa^2v}{\ell}E_1,\qquad
\operatorname{Tr}(\rho'H)\le E_1+\frac{B_G}{4\nu^2}.
\tag{12}
$$

有限E₁由式(3)保证梯度可积，按闭形式和正分解延拓至正常混态，不需要源在D(H)或E₂有限。非选择能量有限意味着条件能量几乎处处有限，但不保证所有选择结果有同一个能量上界。

对更一般的独立中心指针，同一理想V_g给ΔE=g²σ_p²〈G_F〉/ℏ²，因此

$$
\Delta E\,\nu^2=
\frac{\sigma_q^2\sigma_p^2}{\hbar^2}\langle G_F\rangle
\ge\frac{\langle G_F\rangle}{4}.
\tag{13}
$$

实最小Gaussian饱和。这是指定指针／脉冲族的源注能—读噪关系，不是所有仪器的普遍最优界、装置总耗散或Landauer成本。尚未构造有限脉冲期间源H与探针的自主联合演化。

## 7. 有限记录与适用边界

对预选δ∈(0,1)的一次末读，

$$
R=\sqrt{\frac{B_F+\nu^2}{\delta}},\qquad
\Pr\{|Y|>R\}\le\delta .
\tag{14}
$$

区间内可有限分箱加溢出标签，保持非选择能量与Gauss结论。δ是记录认证要求，不是宇宙演化阈值。矩界不保证足够坐标精度、图册满秩或带宽。未测H保持E₁，故任意等待后的单次末读可沿用预算；重复读取须更新能源账。

## 8. 复算

采用544同一L、C、u；ω=1.3，径向源乘s Gaussian，链路常函数。k=1表示一端相关非Gaussian径向激发。

|k；激发端；v；κ；b；ν²|E₁|读数二阶矩|物质注能|链路注能|测后源能量|
|---|---:|---:|---:|---:|---:|
|0；—；1；0.8；0.3；0.36|4.810852|2.632189|1.367521|0.059172|6.237545|
|1；1；0.7；1.1；0.4；0.09|8.893172|4.826183|10.858974|0.584714|20.336861|
|1；2；1.5；0.5；0.7；0.81|6.084636|5.303343|0.534188|0.107865|6.726689|

式(12)注能上界分别10.583752、122.590201、5.905655；实际值均在界内。界很保守，不以单独缩小常数自动增加后续轮次。

六组检查覆盖：局域群作用、Casimir与独立数值梯度；热核归一、衰减及无偏均值；常链路物理投影反例；同势E₁共同界；指针概率、离对角衰减与梯度积分；完整测后波函数导数的物质／电场能量。最后一项在指针动量表示中独立展开，最大注能残差1.78×10⁻¹⁴。

热核用两组阶数独立积分t=0.5、0.75、1、2、4，未裁剪负数修补核，未数值认证t→0渐近。有限样本不代替一般态和定义域的解析证明。

## 9. 条件账与下一项

本轮把C02局域物理部门、C14规范链路、C17相邻作用、C03／C21实际概率和后态、C19来源及C22源反作用接在同一有限候选。该配置边口共用E₁，减少了独立来源矩要求；给定图、群、探针身份、脉冲来源、尺度极限仍未消除。

后继先核同一链路／物质粗粒化是否共同保持规范约束、相邻读取和演化，哪些有效项必须进入同一作用。先复用成熟区域拼接及有效理论映射；只复述一般形式不增加轮次。不继续细扫单边精度，不从SU(2)内部维数签收空间三维。

原失败判据是效应合法但后态不守约束，或同能量预算不足以控制读数／注能。前者排除所列细群仪器，后者在直接不变读口由共同界解除；这不代表探针已自主实现。统一目标仍未完成。
