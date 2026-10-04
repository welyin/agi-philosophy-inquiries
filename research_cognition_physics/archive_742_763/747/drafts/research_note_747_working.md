# 第747轮工作报告：原占据转移与实际异地读口之间

2026-10-04。正式仍[746](../../research_note_746.md)／3430。本入口接[任务范围](STATUS.md)，[代码](remote_population_entry.py)、[结果](remote_population_entry_results.json)、[检查](entry_checks.json)；不计完成轮次。

## 1. 已回查的接口

577已证明原A端相位准备能通过测地边势改变B端实际读数；不能重复计算一个相位响应再声称解决本题。746则只证明同节点未知sterile空／对输入可经原T菜单读出。本轮缺口是**同一未知输入是否能在另一节点留下实际读口记录**。

保持746的正常严格Gauss核准备，输入在A，B初始CAR真空。令原sterile跳跃为

$$
C_{AB}^{\nu}=\nu_B^\dagger\tau_{BA}\nu_A+
\nu_A^\dagger\tau_{BA}^\dagger\nu_B,\qquad
N_B^\nu=\nu_{B,\uparrow}^\dagger\nu_{B,\uparrow}
+\nu_{B,\downarrow}^\dagger\nu_{B,\downarrow}.
\tag{1}
$$

τ来自原保持物种的边spin矩阵。ν是原规范singlet，所以此块不依赖规范链路；其他所有物种、原位质量、图势、链路动能与跳跃都保留。

## 2. 原全H首步已经给出远端占据差

两个输入在B均无粒子，故N_B^(1/2)Ψ_j=0，直接得到

$$
\langle N_B^\nu(t)\rangle_j
=\frac{t^2}{\hbar^2}
\|(N_B^\nu)^{1/2}H\Psi_j\|^2+O(t^3).
\tag{2}
$$

玻色部分保原B占据为0，A的质量也不改变B占据。B的Majorana在两个输入都创建同样的B对；A有对时原跳跃还可将其中一个送到B。这些像的A、B粒子数不同，彼此正交，故

$$
\langle N_B^\nu(t)\rangle_1-\langle N_B^\nu(t)\rangle_0
=\frac{t^2}{\hbar^2}\operatorname{tr}(\tau_{BA}^\dagger\tau_{BA})
+O(t^3).
\tag{3}
$$

多条同端点边必须先合并其实际τ再取范数，不能把干涉误作概率相加。其他邻居的初始空态不提供本轮差异。若τ=0，该机制没有信号，不排除标量边通道。

8组完整64模、非平凡规范链路及复spin矩阵的校准，原τ矩阵的Hilbert–Schmidt平方均为0.354；B的共同Majorana背景项随原x₅,B改变，差精确到校准精度。它是全H首步的CAR校准，不是全玻色时间传播仿真。

## 3. 尚不能签收异地实际记录

式(3)使用N_B这个可观测量来诊断转移。当前问题恰恰要求由原T_B读口实际读取，因此**不在中间插入一个尚未实现的N_B理想仪器**。

真实待算的是

$$
p_{j,+}^{B}(t)=
\langle U(t)\Psi_j,E_{T_B,+}U(t)\Psi_j\rangle,\qquad
E_{T_B,+}=\frac12+\frac14\sin(|X_B|^2/2).
\tag{4}
$$

746的局部四阶系数不能直接与式(3)相乘：未知态、相关和连续演化中的所有时间次序需要同算。原标量边传播和费米跳跃可能在相同阶次出现，不能只保其中一条正贡献而忽略另一条。

下一步首先做算符词的阶次／支持筛查，并在原完整图中确定实际非零系数。若必须增加条件控制，则把控制作为新操作输入记账，不将其冒称原自治传播。原同态来源和资源继续沿623—625、718、746处理；目标与旧空间结论保持。
