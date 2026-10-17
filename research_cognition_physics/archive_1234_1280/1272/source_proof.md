# 1272来源分项：同一固有钟、有限热载体与完整来源

2026-10-17。按[准入](admission.md)和[规格](candidate_spec.md)计1个作者核验组。九门的实际热任务由另一分项核验。新增认知公理0，CT未采用。

## 1. 对象、能隙和输入量词

旧对象 $O$ 包含1270原七对、转子、几何指针、旗标和输出环境，真实生成元为 $H_O(\theta)$，不冻结旧几何。
新增
$(A,B,S_A,S_B,Z_A,Z_B,K,L,M,E,J_{\rm out})$
的维数为 $(2,2,2,2,4,5,2,3,2,2,2)$，总维数15360。$J_{\rm out}$是新记录，不是钟共轭量 $J$。每载体正基线1，合计11。先独立指定

$$
h_D=|1\rangle_D\langle1|,\quad
H_{\rm th}=h_A+h_B+h_{S_A}+h_{S_B},\quad
h(\theta)=11I+H_{\rm th}+v(\theta).
\tag{1}
$$

另7个载体在能源字典中简并；其种子纯度、库存、保留日志和未来复位不免费。

原实际尺度保持

$$
\mathcal A=V_0s^2,\quad V_0=10^{-6},\quad
\varepsilon=\hbar/\mathcal A,\quad
s=1000m2^{2N},\quad m=\lceil\max(1,\hbar)\rceil,\quad\tau=s\theta .
$$

故新增总物理Hamiltonian及独立能隙为

$$
E_T=V_0s\varepsilon h=\frac{\hbar}{s}h,\qquad \Delta=\hbar/s>0.
\tag{2}
$$

与同一固有钟 $\tau=s\theta$匹配，不从目标人口反定义能量。
随后人口 $1/4,2/5$对应
$\beta_A^{\rm phys}=(s/\hbar)\log3$、
$\beta_B^{\rm phys}=(s/\hbar)\log(3/2)$（$k_B=1$）。
只称两个约化有限能源部门为Gibbs参考，不声称联合态Gibbs、加速温度或Jacobson温标。

初始 $A=B=0,S_A=S_B=1$，完整新侧支持在
$\Pi_2=\mathbf1_{\{2\}}(H_{\rm th})$。未知 $K/R$、种子纯化和允许的准备误差都留在新侧并保此支持。

**共同几何预测要求新热准备与原指定 $O$ 准备独立。**
完整通道和算符身份可以保任意参考；若预先跨新旧关联，再按新热结果条件化，则可能重估旧几何，不能直接沿用旧条件误差。

## 2. 内部辛项、约束与张量传播

取归一 $z\in\mathbb C^{15360}$的射影辛形式
$a_{\rm FS}=\tfrac i2(z^\dagger dz-dz^\dagger z)$，同流增量

$$
\Delta S=\varepsilon\int j^\mu(a_{\rm FS})_\mu
-\varepsilon\int\sqrt{-g}\,n\langle z,h(\theta)z\rangle .
\tag{3}
$$

原 $j^\mu J\partial_\mu\theta-\sqrt{-g}\,nJ$保留。变分给
$U^\mu\partial_\mu\theta=1$和
$iU^\mu\partial_\mu z=hz$（允许射影相位规范）。
齐次前因子 $\mathcal A$使新增辛项系数为
$\mathcal A\varepsilon=\hbar$；漏掉辛项的 $\varepsilon$会改变能量—时间字典。
这是一份明确的有限内部模型，不代表每个空间标签可独立访问同一集体算符。

真实有效完成选为

$$
H_{\rm tot}=H_O\otimes I+\varepsilon I\otimes h.
\tag{4}
$$

**新内热项不再乘旧 $\chi$。**
若改为 $H_O+\varepsilon\operatorname{Op}(\chi)\otimes h$，一般不能使用本件分离。
式(4)在原内核吻合同流作用，管外是明示有界完成。所有新旧因子及旧门交叉强对易，所以

$$
U_{\rm tot}=U_O\otimes U_T,\qquad i\partial_\theta U_T=h(\theta)U_T.
\tag{5}
$$

九门相互不必对易。即使与旧写入窗重叠也成立。
同一关系约束给

$$
J_{\rm new}=J_{\rm base}\otimes I-\varepsilon I\otimes h,\qquad
J_{\rm new}+\varepsilon h=J_{\rm base}\otimes I.
\tag{6}
$$

旧 $J_{\rm base}$仍可随旧写入、读取变化。

## 3. 一次提款与非对易门工作

受控SWAP交换等隙部门；人口控制记录只改简并载体。因而

$$
[G_j,H_{\rm th}]=0,\qquad U_T^\dagger H_{\rm th}U_T=H_{\rm th}.
\tag{7}
$$

初片 $v(0)=0$及实际 $\Pi_2$支持给一次固定提款

$$
J_{\rm new}(0)=J_{\rm base}(0)-13\varepsilon .
\tag{8}
$$

11支付载体基线，2支付库存激发；种子与空白纯度仍是明确初始资源，不提供免费热制备。提款不依赖未知控制、程序或结果。

不同门非对易时，正确Heisenberg身份为

$$
\frac d{d\theta}(U_T^\dagger hU_T)
=U_T^\dagger(\partial_\theta v)U_T,
$$

$$
\frac d{d\theta}(J_{\rm new}^H-J_{\rm base}^H)
=-\varepsilon U_T^\dagger(\partial_\theta v)U_T,
$$

$$
J_{\rm new}^H-J_{\rm base}^H
=-\varepsilon(11+H_{\rm th})-\varepsilon U_T^\dagger vU_T .
\tag{9}
$$

不能按裸矩阵积分 $v'$。门与门之间及全部门结束时 $v=0$，实际支持仍有 $-13\varepsilon$；门功归零不退还基线和库存。实际幅误差、起点偏移用同一个 $v$补偿，不以理想门代替实际门。

全部种子、库存、记录、撤销前后的关联和参考保留。末读 $J_{\rm out}$基效果与 $H_{\rm th}$对易，所有非零概率条件支仍在 $\Pi_2$；不免费擦除已公布记录。

## 4. 独立来源与同完成抵销

固定稠密流 $j^\mu$和内部变量作度规变分。新增质量不依赖杆导数或机械速度，式(3)在解约束前给

$$
\Delta T^{\mu\nu}=n\varepsilon hU^\mu U^\nu,\qquad
T_{c,\rm new}^{\mu\nu}
=n(J_{\rm new}+H_{0,\rm mat}+\varepsilon h)U^\mu U^\nu.
\tag{10}
$$

这是共动质量源，无新增杆压力；不是1269体积型反馈。先独立定义源再代入式(6)，完整质量恢复原 $J_{\rm base}+H_{0,\rm mat}$，不以 $G_{\mu\nu}/\kappa$定义物质。

对旧来源系数 $\mathcal B$，钟与新质量必须用同一系数及固定线性Weyl／完成映射 $\mathcal L$：

$$
\mathcal L[\mathcal B(J_{\rm new}+\varepsilon h)]
=\mathcal L[\mathcal B J_{\rm base}]\otimes I,\qquad
\mathcal L[\mathcal B\varepsilon h]
=\varepsilon\mathcal L[\mathcal B]\otimes h.
\tag{11}
$$

即使 $\mathcal B$与旧钟不对易，旧排序也原样保留。只对钟乘截止或换系数通常会破坏抵销。代码核非对易实例及删除检查；一般身份由线性和张量结构承担。

这保所有相干和条件支，不靠平均能源为零。各部门能源仍改变：$A,B$、库存和 $J_{\rm new}$分别变化。新增载体不是零源，而是有已记账钟补偿；本模型没有新净热致几何响应。

## 5. 峰值、真实有限资源与正域

九窗各长 $1/16$、相继执行，归一 $C^2$脉冲峰35，
$|\xi_j|\le10^{-9}$、$\|G_j\|\le\pi<22/7$。因此

$$
\|v\|\le35\pi(1+10^{-9})<111,\quad
\|h\|<126,\quad \|h\Pi_2\|<124.
\tag{12}
$$

本候选的撤销门再次使用同一个正CNOT生成元（理想CNOT自逆），不使用负生成元；$1+\xi_j>0$，所以每窗 $v\ge0$、全新空间 $h\ge11I$、实际支持上 $h\ge13I$。仍用式(12)上界控制钟提款；门功可有正负，不能与交互能的正性混同。
足够的导数界为 $\|g'_j\|_\infty\le6720$，
故 $\|\partial_\theta v\|<21121$。物理交互能峰 $<111\hbar/s$，工作率峰 $<21121\hbar/s^2$，初提款 $13\hbar/s$，每门固有时长 $s/16$；有限作用量没有消失。

直接读取[1249原资源整数](../1249/quantum_bounds_results.json)和[1267统一证书](../1267/geometry_results.json)，不重新挑小参数：

$$
M=2^{3059},\ E_H=3093,\ E_0=14008,\ E_F=95928+100M,
$$

$$
N=1+\max\left\{\left\lceil\frac{E_F+42}{3}\right\rceil,
\left\lceil\frac{E_H+E_F+42}{6}\right\rceil\right\},\quad
4N-7>100,\quad4N-15>100 .
\tag{13}
$$

于是 $\varepsilon\le2^{-4N}$、
$126\varepsilon<2^{-100}$、$21121\varepsilon<2^{-100}$。
固定提款也远小于保存的全参数正域阈值余量；代码用精确有理和整数核验，不仅说“足够小”。

在指定时刻，对任意正常态和阈值 $j$，实际 $\Pi_2$上有谱平移

$$
\Pr(J_{\rm base}\le j-124\varepsilon)
\le\Pr(J_{\rm new}\le j)
\le\Pr(J_{\rm base}\le j+124\varepsilon).
\tag{14}
$$

这是强对易因子的谱序身份；全新有限空间将124换126。门外在实际支持上恰为 $13\varepsilon$平移。
旧实际谱／坏域预算须在对应平移阈值调用；经典 $J>0.098$不是全Hilbert空间正谱或连续监测永不越界证明。原完整质量及几何不变，钟这一单项明确改变。

## 6. 同一真实后态和原几何承接

指定输入 $\rho_O\otimes\rho_{TKR}$始终演化为

$$
\rho_{\rm tot}(\theta)=\rho_O(\theta)\otimes\rho_{TKR}(\theta).
\tag{15}
$$

两侧实际输出仪器分别作用旧侧和新侧，每个联合未归一结果是两份未归一输出的张量积。任何非零概率新热结果不改变旧条件后态。因此原正常旗标的 $p_\pm\ge1/4$及其 $2^{-80}$来源／几何预算保留，无需按新热结果概率再次放大旧误差。未知 $KR$边缘不必不变；未参与参考 $R$单独边缘不变。

九门最晚在 $17/16+1/256$完成，早于旧几何采样 $[2-10^{-8},2]$；公共未来 $17/8$提供新热记录及原几何记录。具名有效边界权限由另一分项承担，不授窗内任意细端口，也不继承K1219运动、捕获或反冲合同。

式(11)、(15)精确承接[1270](../1270/continuation_proof.md)的完整有界源及独立比较几何。实际幅度和输入误差只在新侧，源身份仍成立；不用“热概率误差乘巨大源范数”证明来源小。
对原有限测试 $\sum_{ab}\int|w^{ab}|\le1$，

$$
\left|\sum_{ab}\int w^{ab}
(G_{ab}[g_\kappa]-3\langle\widehat T^c_{{\rm tot},ab}\rangle)\right|
\le\frac{111}{100}\Delta_\kappa+3\cdot2^{-80},
\quad\Delta_\kappa=\frac1{15000}+5\cdot10^{-16}.
\tag{16}
$$

这是原比较几何测试，不是量子度规算符方程。旧条件旗标的测试支持仍限原实际发布后的窗口；新热条件化不解除时间范围。额外外部读者、任意未知几何态或跨新旧关联准备不自动领取本误差。

## 7. 核验与停止线

[source_check.py](source_check.py)计1组：有理峰值／时窗／提款、旧真实整数余量、守恒激发支持、非对易门及同系数来源抵销的有限矩阵。默认只读比对[source_results.json](source_results.json)，首次 --write 排他创建。一般全态结论由解析张量和谱身份承担；不重跑旧ODE／量子传播或热参考分类。

本轮具体材料和门是明示输入，治理机制未唯一选择它们。与真实热任务共同验收后停止增加部门、库存和门。未选 $\kappa$，未证几何熵、面积系数、Jacobson平衡或新热致几何响应；返回这些剩余身份，不继续修辅助装置。
