# 第747轮：原跳跃对异地实际读口的六阶贡献

2026-10-04。接[746](../../research_note_746.md)及[已执行入口](research_note_747_working.md)。[代码](../joint_remote_record_response.py)、[结果](../joint_remote_record_response_results.json)、[核验](../research_round_747_checks.json)、[完整条件账](../unified_physics_condition_ledger_747.md)。

**本轮严格结论限于一个分量：原费米跳跃对远端原Higgs读口产生非零的六阶响应系数。原标量边也能在同阶贡献，尚未求总和，故未完成异地实际记录。** 不将局部读出概率与占据转移概率相乘，不改原边、不添加中间测量。

## 1. 认知动机、建模输入与历史去重

统一候选需要不同主体通过同一量子过程留下记录。577已完成原相位输入的传播，746完成原未知sterile偶编码的本地T读口；本轮继续同一空／对输入的异地读取。空间382—386、425、522—523及604、649／699保持原范围。

采用原有限图、固定正外部几何和正常Gauss准备。分解一条边时取原两节点一边实例，全部物种、质量、标量边、链路和动能存在：

$$
H=T_A+T_B+T_{\rm link}+V_{\rm sc}+M_A+M_B+C_{AB},
\qquad T_v=-\frac{\hbar^2}{2w_v}\Delta_{K,v}.
\tag{1}
$$

M含原Dirac与Majorana，C含所有保持物种的跳跃。固定几何不是已解出的动力学引力；w与k遵守原共同几何。原Yukawa校准值与具体spin矩阵是建模输入，非认知原则的推论。

$$
\Psi_0=\Phi_A\Phi_B\Phi_{\rm link}|{\rm vac}\rangle,\quad
\Psi_1=\Phi_A\Phi_B\Phi_{\rm link}
\nu_{A,\uparrow}^{\dagger}\nu_{A,\downarrow}^{\dagger}|{\rm vac}\rangle,
\qquad E_{B,+}=\frac12+\frac14\sin T_B^{\rm obs}.
\tag{2}
$$

配对整体相位不影响这两种输入的概率比较。Phi为反射偶、规范不变准备，链路为正规化Haar常函数。旧菜单T^obs=|X|²/2区别于式(1)动能：

$$
x=\phi/\sqrt F,\quad a=|x_H|^2,\ b=x_5^2,\ u=1+(a+b)/6,\qquad
T^{\rm obs}=a/u,\quad \nabla_K T^{\rm obs}=(2x_H/u,0).
\tag{3}
$$

入口已在原64模、非平凡规范链接中核验

$$
\Delta\langle N_B^\nu(t)\rangle
=\frac{t^2}{\hbar^2}\operatorname{tr}(\tau^\dagger\tau)+O(t^3),\qquad
\tau=\begin{pmatrix}0.4&0.2+0.13i\\-0.17+0.11i&0.31\end{pmatrix},
\quad\operatorname{tr}\tau^\dagger\tau=\frac{177}{500}.
\tag{4}
$$

原B端Majorana共同背景保留并相消。式(4)不能替代E_B的响应，过程中没有插入N_B仪器。

## 2. 六阶跳跃分量的解析筛查

用lambda仅标记跳跃次数，实际物理模型始终lambda=1：

$$
H(\lambda)=H-C_{AB}+\lambda C_{AB},\qquad
d_6(\lambda)=\left.\partial_t^6
\left(p_{1,+}^B(t;\lambda)-p_{0,+}^B(t;\lambda)\right)\right|_{t=0}.
\tag{5}
$$

这不是独立改变原物理速度或共同几何的许可。各项可能干涉，一个lambda系数非零不代表d6(1)非零。

取A_B=E_B,+−1/2。右端首个非零对易必为T_B。在B端质量出现前，跳跃与B标量函数／微分系数对易。首个质量力为

$$
[M_B,[T_B,A_B]]
=\frac{\hbar^2}{w_B}\langle\nabla_K A_B,\nabla M_B\rangle_K
=\frac{\hbar^2}{w_B}\frac{2f'(T_B^{\rm obs})}{u_B}D_B,
\qquad f(z)=\frac14\sin z .
\tag{6}
$$

至少有一个B端Dirac因子。奇数次跨边跳跃在A留下奇CAR次数，不能与原空态或配对态对角闭合。两个跳跃加一个质量，仍有一次Dirac物种切换或一次异常配对，不能形成对角差。

六个算符中，两个跳跃的非零分量恰有两个动能和两个质量。若仅一个动能、三个质量，奇Majorana次数改变总粒子数，偶Majorana次数留下奇Dirac次数，都不闭合。标量势替换第二动能时，所有可存活的后续算符为乘法矩阵，标量势与它们对易，贡献零。因此

$$
[\lambda^2]d_6(\lambda)
=\frac{i^6}{\hbar^6}\Delta\!\left\langle
\sum_{\#T=2,\,\#M=2,\,\#C=2}
\operatorname{ad}_{O_6}\cdots\operatorname{ad}_{O_1}A_B
\right\rangle .
\tag{7}
$$

一质量在A、一质量在B的项对x_A整体反演为奇；原A准备和T_A为偶，其积分为零。没有B质量的项从右端即为零。所以存活项两质量均在B，两个动能也为T_B；T_A与链路动能与此存活子式对易。A的输入差迫使两次相连跳跃属于sterile块，规范表示平凡。

这只是**指定系数的缩约**，不将实际完整H改成34模模型。四跳加一质量、一动能同样不能对角闭合；奇跳全部为零。于是

$$
d_6(\lambda)=a_6+\lambda^2c_6,\qquad c_6=[\lambda^2]d_6(\lambda).
\tag{8}
$$

**a6包含原标量边传播，尚未计算，不能设为零。** 低阶完整远端响应也不靠式(8)单独签收。

## 3. 可复算的精确六阶密度

[精确程序](remote_hop_coefficient.py)保B的完整32模质量，附加A两sterile模，在五个真实玻色坐标上无截断地计算整数CAR及多项式。暂取hbar=w_B=1、单位spin跳跃。对Gaussian乘齐次多项式，原H5动能为

$$
\Phi^{-1}\left(-\frac12\Delta_K\right)(\Phi P_d)
=-\frac12\Delta_{\mathbb R^5}P_d
+\left(\frac52+\frac{2d}{3}-\frac{d^2}{12}
+\frac d6|x|^2-\frac1{12}|x|^4\right)P_d,
\quad\Phi=e^{-|x|^2/2}.
\tag{9}
$$

质量分母100，统一每步分母300。向量jet保(T,M,C)次数，最后取(2,2,2)。若H^k Psi=Phi P_j,k，则读口密度由

$$
Q_6(x)=\left[
i^6\sum_{k=0}^6(-1)^k\binom6k
\left(P_{1,6-k}^{\dagger}P_{1,k}
-P_{0,6-k}^{\dagger}P_{0,k}\right)
\right]_{(2,2,2)}
\tag{10}
$$

给出。B规范不变性允许最后把Higgs四向转到径向r；微分先在全部五坐标完成，未删横向导数。a=r²,b=s²时，得到

$$
\begin{aligned}
Q_6(a,b)={}&-\frac{7912}{625}-\frac{5987}{1250}b
-\frac{6773}{5625}b^2-\frac{521}{2250}b^3
-\frac{3201}{1250}a+\frac{507}{625}ab-\frac{17}{225}ab^2\\
&+\frac{11336}{5625}a^2+\frac{1223}{2250}a^2b+\frac{436}{1125}a^3 .
\end{aligned}
\tag{11}
$$

它不是746局部P4的常数倍，“先转移再本地读取”的概率乘法没有复现原连续演化的全部时间次序。第二次精确计算使用式(4)原复spin矩阵，逐个有理系数核实

$$
Q_{6,\tau}=\frac{177}{1000}Q_{6,I_2}.
\tag{12}
$$

本轮只使用这个已声明tau，不把两个矩阵的核验冒充任意新跳跃结构的分类定理。

## 4. 严格符号与物理核准备

复用746区间积分算法，输入新多项式。约去共同S3角体积：

$$
\begin{aligned}
Z&=\frac43\int_0^\infty\frac{R^4e^{-R^2}}{\sqrt{1+R^2/6}}\,dR>0,\\
J_6&=\int_0^\infty\frac{R^4e^{-R^2}}{\sqrt{1+R^2/6}}
\int_{-1}^{1}(1-z^2)Q_6(R^2(1-z^2),R^2z^2)
\frac{\sin\!\left(\frac{6R^2}{6+R^2}(1-z^2)\right)}4\,dz\,dR .
\end{aligned}
\tag{13}
$$

[证书](remote_hop_sign_certificate_results.json)给出含解析余项、sin尾项、R>10尾部和向外舍入的区间：

$$
-0.505363087653744092<J_6<-0.505363087653695236<0,\qquad
c_6=\frac{177}{1000}\frac{J_6}{Z\hbar^2w_B^2}<0.
\tag{14}
$$

三个误差上界分别小于2.440×10^-14、3.625×10^-17、5.245×10^-21。这是系数认证，没有模拟原全图玻色时间传播。

Gaussian是系数比较向量。按746取光滑偶径向紧支撑逼近，落在原Gauss核内，全部有限H矩存在。本系数阶数有限、微分系数至多多项式增长，带权L²收敛使

$$
c_{6,R}\longrightarrow c_6<0 .
\tag{15}
$$

这里只转移跳跃分量符号，总系数a6,R+c6,R仍待核。未给截止半径或可用等待窗。

## 5. 同一后态、资源及下一步

实际输出继续是

$$
\mathcal J_{B,r,t}(\rho)=
L_{B,r}U(t)V\rho V^\dagger U(t)^\dagger L_{B,r},\quad
L_{B,r}=\sqrt{E_{B,r}},\qquad
0\le\mathcal I_B^*(H)-H\le\frac{\hbar^2}{9w_B}.
\tag{16}
$$

未知输入未换成经典消息；完整后态及来源沿623—625、718、746核算。准备、末端仪器和几何仍是声明输入。自治终端、无限系统、连续极限及引力闭合未获证明。

两组核验分别复算入口完整CAR首步和本轮精确多项式／严格积分；入口此次纳入正式核验，没有额外计一轮。关掉Ys的探索诊断保存在草稿，属于改变模型的机制比较，没有据此替换原参数或宣布解析抵消。

**新增结果：** 原跳跃确实进入旧远端读口的非零六阶分量，且不能把两项旧概率简单相乘。**待解缺口：** 原标量边会加强还是抵消该分量。

下一项[748](../../748/drafts/STATUS.md)：保原共同几何和同一准备，求a6及d6(1)。先从原Duhamel展开连接746局部密度与原测地边力，保全部时间次序；不改弱边来避开总和。如果只得到数值证据，明确与严格符号证明区分。
