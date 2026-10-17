# 1256后局部程序：双流辅助源拆分与Leray准入筛查

2026-10-16。作者：`/root/spatial_parent_next`。非编号准入说明；只处理原连续候选的一个数学接口，不增加认知公理，不改变CT状态，不修改导航或冻结文件。

**判定：双流源可以作一份等价的辅助拆分，使所需最高阶符号成为因果的准对角Leray形式；两种流速重合本身不是这个改写的障碍。** 原C²脉冲的导数阶也没有立即失配。但本文未完成该混合系统的有限正则能量估计、解的连续依赖及几何唯一性，故没有证明局部程序扰动存在至原写窗末端，也不把这些候选技术条件列为认知公理。停止继续修PDE。

## 1. 原对象、旧接口及本次范围

继承[1245](../../1245/proof.md)的CT条件选后Einstein分支，以及[1246](../../1246/proof.md)、[1247](../../1247/proof.md)同一双标签流、标量、杆、程序、实际钟与K记忆。[上一接口](../local_program_domain_interface_after1256.md)已经说明局部B准备的初始全源不变和条件性域依赖，不重复记为新结果。

本次新增仅为：不把两个沿流导数乘到同一波方程上，而按来源引入固定初数据的辅助波变量，检查等价性与最高导数。辅助变量不是新增物质，不分别分配物理应力，不拥有可调的实际准备参数。

原实际菜单只认证均匀B=1、2。局部B=1+b(x)仍是候选新增准备。本文不构造该准备器、量子父过程或新终末仪器。

## 2. 先以完整质量替换J，保留所有原内部量

用y表示原记忆坐标Y，避免与下文共同波变量混淆。第一流质量仍为

$$
M=\frac{71}{100}+J+\frac{q^2+r^2+K^2}{2}
+\alpha\epsilon w(\theta)qp+g f(\theta)Xy,
\qquad \alpha=B-1,\quad \epsilon=\frac14,\quad g=\frac1{10}.
\tag{1}
$$

由于J的系数恰为1，可以用M替代J作为演化变量，最后按(1)恢复J。这是可逆变量替换，没有删去钟余额。

记$T_c=U_c^\mu\partial_\mu$。原内部方程给

$$
\begin{aligned}
T_c\theta&=1,&T_cB&=0,&T_c\psi&=\epsilon wqp,\\
T_cq&=r,&T_cr&=-q-\alpha\epsilon wp,\\
T_cy&=K,&T_cK&=-g fX,\\
T_cJ&=-\alpha\epsilon w'qp-g f'Xy.
\end{aligned}
\tag{2}
$$

对(1)求沿流导数，振子、记忆及交互的内部项两两抵消，得到精确身份

$$
T_cM=\alpha\epsilon wq\,T_cp+g fy\,T_cX.
\tag{3}
$$

反过来，(3)与其余(2)恢复原J方程。因此以M为变量的第一阶内部运输方程只用w、f，不用w′、f′。第二流保持原I、$\vartheta$，$T_rI=0$、$T_r\vartheta=(1+p)/10$，$M_r=(1+(1+p)I)/10$。两股流的全部内部量都保留。

质量在位置变分中的偏导仍是原作用的$M_{s,p}$、$M_{s,X}$，不能在改变量后误设为0：

$$
M_{c,p}=\alpha\epsilon wq,\quad M_{c,X}=g fy,
\qquad M_{r,p}=I/10,\quad M_{r,X}=0.
\tag{4}
$$

## 3. 将完整约化方程写成按流分开的源

在同一谐和规范中，令

$$
\mathcal Y=(g_{\mu\nu},p,X,\phi^1,\phi^2,\phi^3),
\qquad P_g=g^{\mu\nu}\partial_\mu\partial_\nu.
$$

杆的固定绕数可留在已知仿射部分，用周期偏差作局部变量；这不改变以下导数阶。F=1+p/5>0使标量动能矩阵可逆。完整约化波方程具有

$$
P_g\mathcal Y=\mathcal B(\mathcal Y,\partial\mathcal Y)
+\nu_c\mathcal A_c(\mathcal Y,U_c,z_c)
+\nu_r\mathcal A_r(\mathcal Y,U_r,z_r).
\tag{5}
$$

$\nu_s$为原proper数密度；$z_s$为保留下来的内部量，第一流使用M而非J。$\mathcal B$含完整几何的一阶导数项、标量势、F的波映射项和必要的规范项。每份$\mathcal A_s$不含来源密度或场的导数。以$\kappa=3$的四维trace-reversed约定为例，其对应分量是

$$
(\mathcal A_s)_{g_{\mu\nu}}
=-2\kappa M_s(U_{s\mu}U_{s\nu}+g_{\mu\nu}/2),
\quad (\mathcal A_s)_p=M_{s,p},\quad
(\mathcal A_s)_X=M_{s,X},\quad
(\mathcal A_s)_{\phi^A}=0.
\tag{6}
$$

使用未约化的协变材料方程时，标量的契约Christoffel项并入$\mathcal B$，不能为凑主符号删除它。物质运动及数守恒仍是

$$
T_sU_s^\nu=-\Gamma^\nu_{\alpha\beta}U_s^\alpha U_s^\beta
-\frac{g^{\nu\mu}+U_s^\nu U_s^\mu}{M_s}
(M_{s,p}\partial_\mu p+M_{s,X}\partial_\mu X),
$$
$$
T_s\nu_s=-\nu_s\nabla_\mu U_s^\mu.
\tag{7}
$$

交互通过同一p、X和g发生；没有新增跨流密度梯度，也没有把总T倒定义成$G/\kappa$。

## 4. 辅助源拆分及等价初值

取

$$
\mathcal Y=\mathcal Y_0+\mathcal Y_c+\mathcal Y_r,
$$
$$
P_g\mathcal Y_0=\mathcal B,
\qquad
T_s(P_g\mathcal Y_s)=\nu_s
\bigl[T_s\mathcal A_s-\mathcal A_s\nabla_\mu U_s^\mu\bigr],
\quad s=c,r.
\tag{8}
$$

在右侧用原(2)—(4)、(7)代入$T_s\mathcal A_s$。它只含总$\mathcal Y$的一阶导数、该流速度的一阶导数和内部量；没有$T_c\nu_r$或$T_r\nu_c$。两份源在共享场中的相互影响仍全部保留。

在原初片固定

$$
\mathcal Y_s|_0=\partial_t\mathcal Y_s|_0=0,
\qquad
g^{00}\partial_t^2\mathcal Y_s|_0
=\nu_s\mathcal A_s|_0,
\tag{9}
$$

并将原物理$\mathcal Y$、$\partial_t\mathcal Y$全部赋给$\mathcal Y_0$。零辅助Cauchy数据使(9)中空间及混合导数为零；$\mathcal Y_0$的二阶初导数由其波方程确定。故总二阶初导数恰是(5)，不修改原几何或约束初值。

令

$$
E_s=P_g\mathcal Y_s-\nu_s\mathcal A_s.
$$

由(7)—(8)精确得到$T_sE_s=0$，且(9)给$E_s|_0=0$。只要所讨论的解具有足够正则性，运输唯一性给$E_s=0$，相加即恢复原(5)。反向给定一份原解，在该固定g上解$P_g\mathcal Y_s=\nu_s\mathcal A_s$、零位置/速度初数据，并定义$\mathcal Y_0$为差，即恢复(8)。这是固定辅助初值下的等价扩张；不为辅助分量分别要求Lorentz号差或正能源。

这份代数等价不自动证明约化Einstein方程的规范约束传播；该项仍需完整系统的存在性与原总源守恒一起验收。

## 5. 准对角最高阶符号：速度重合不要求排除

取未知量的Leray指数

$$
\ell(\mathcal Y_0,\mathcal Y_c,\mathcal Y_r)=3,
\quad \ell(U_s,z_s)=2,\quad \ell(\nu_s)=1.
$$

方程指标分别为

|方程|指标m|对角最高阶|其它出现的最高导数|
|---|---:|---|---|
|$\mathcal Y_0$波方程|1|$P_g$，二阶|总$\mathcal Y$一阶|
|$\mathcal Y_s$运输波方程|0|$T_sP_g$，三阶|其它波变量一阶、$U_s$一阶、内部及密度零阶|
|$U_s$运输|1|$T_s$，一阶|总$\mathcal Y$一阶、内部零阶|
|$z_s$运输|1|$T_s$，一阶|总$\mathcal Y$一阶、其余内部零阶|
|$\nu_s$连续|0|$T_s$，一阶|$U_s$一阶、总$\mathcal Y$一阶|

其它列的导数阶均严格小于$\ell(列)-m(行)$，因此主矩阵是对角的。$T_s$作用于$P_g$系数产生的$(\partial g)\partial^2\mathcal Y_s$也是低阶项，不能漏计于后续能量估计。

在一点的局部Lorentz标架、空间频率$\xi\ne0$下，归一化对角三阶多项式为

$$
p_s(\lambda,\xi)=(\lambda-v_s\cdot\xi)
(\lambda^2-|\xi|^2).
\tag{10}
$$

每流严格类时，即$|v_s|<1$，便使三个根分离。两流即使$v_c=v_r$，也是两个独立对角块；没有把$T_cT_rP_g$放到同一个块上，故不用为它们的相交引入速度差下界。这一符号判定只处理传播主部，不等于完成低正则非线性适定性。

## 6. 原C²脉冲：有限阶账本可闭合，但定理尚未签署

原$w=64\theta³(1-\theta)³$在[0,1]内，外部为0；f是同一函数平移至[1,2]。零延拓满足

$$
w,f\in C^{2,1}\cap W^{3,\infty},
\qquad w',f'\in C^{1,1}\cap W^{2,\infty},
$$

而不是C³或C∞。例如w'''(0+)=384、w'''(1−)=−384，外部三阶导数为0。不能为了引用光滑定理偷偷修改原脉冲。

三空间维的一份候选混合尺度是

$$
\mathcal Y_\bullet\in H^4,\qquad
U_s,z_s\in H^3,\qquad \nu_s\in H^2.
\tag{11}
$$

|估计对象|需控制的源阶|脉冲最高要求|当前检查|
|---|---|---|---|
|二阶$P_g\mathcal Y_0$|H³|无原脉冲直接源|$\mathcal B$的一阶场乘积阶数相容|
|三阶$T_sP_g\mathcal Y_s$|H²|w′、f′的两阶空间组成|只需w/f到三阶弱导数；$\partial U_s$亦为H²|
|$U_s$及$z_s$一阶运输|H³|w、f的三阶空间组成|质量替换后不含w′、f′的一阶源|
|$\nu_s$一阶运输|H²|无新增脉冲导数|$\nabla U_s$为H²|

在有界H³球上，弱链式法则给w(u)的H³有界性；H²差分可用w''的Lipschitz性控制。w′的对应阶数下降一阶。**这些是源的组成账本，不是完整Leray能量估计。** 特别不能把H³有界性称作H³局部Lipschitz性。

尚未完成：该混合尺度下三阶块的变系数能量估计、统一有限正则光滑化极限、低阶差分及高阶连续依赖，以及完整规范/几何唯一性。原C²脉冲在非均匀钟龄下不允许直接按固定坐标t=1、2逐段拼接；若使用$\theta$切片，还需核其空间性、共同规范和接续估计。

## 7. 与现成定理及公共时间窗的关系

[Choquet-Bruhat–Friedrich, Motion of Isolated Bodies](https://arxiv.org/abs/gr-qc/0605096)的定理2.1与证明给单尘流的T_uP_g改写和混合Sobolev指标，三维局部存在要求s>7/2；其几何唯一性推论2.2采用更强s>9/2。本文复用方法，不能把该单流结论直接当作已涵盖(8)的定理，尤其不能用s=4自动领取较强几何唯一性推论。

[Friedrich–Rendall, The Cauchy Problem for the Einstein Equations](https://arxiv.org/abs/gr-qc/0002074)说明尘的密度导数损失和Leray路径。[Disconzi–Shao, Strongly hyperbolic quasilinear systems revisited](https://arxiv.org/abs/2308.09851)提供明确的局部存在、连续依赖和延拓框架，但针对其已声明的一阶系统假设；本文未把(8)完整归约并核实那些假设，不能只引用题名完成准入。

若(11)的完整有限正则适定性/连续依赖日后得到证明，则原均匀[0,9/4]解的严格正域、有限范数与无焦散余量，可由紧时间区间的有限次延拓给一个非零初数据邻域。固定$b_0\in C_c^\infty$，足够小的$b=a b_0$即可获得共同有限发展；用原9/4与写窗2之间的余量控制实际$\theta$末端，而不是默认扰动后的坐标t仍等于$\theta$。这只是说明所需连接是什么，本文没有给出该邻域半径或公共窗证书。

较大平台准备可能只需局部发展包含对应内域diamond，而不需整环面统一末时片；那是另一个需核局部唯一性/延拓的任务，不能用小数据条件顺带覆盖。

## 8. 可复算内容、净收益和停止线

[check.py](check.py)只用Python标准库Fraction，默认只读比对[results.json](results.json)。它精确核：

1. (3)的完整内部能量抵消及源残差运输代数；
2. 原脉冲的端点二阶匹配与三阶跳变；
3. 上述Leray导数阶表的非对角严格降阶；
4. 一组有理同速/异速样本的三阶主块特征向量及严格根间隔。

它不求解PDE，不认证所有主符号点，不证明能量估计、公共时间窗、量子因果性或实际K误差。全符号量词由(10)的解析不等式承担，有限样本只防实现错误。

**净收益：**原双流与共同场并非只能停在“主锥看起来因果”；已得到保留同一来源和全部内部量的等价准对角改写，明确了一个有限正则接口。**停止线：**本次到此，不继续修此候选的PDE，不把未完成技术接口升级为认知公理或统一研究失败，也不宣布原局部程序菜单已准入。
