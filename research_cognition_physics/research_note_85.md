# 第八十五轮：共同取向结构怎样给出复接口

日期：2026-09-15。接续第 83 轮的受限整体路线与第 84 轮“保留并扩展能力”原则。本轮得到条件性对应，同时找出一条不能省略的额外要求。

## 1 明确新增结构

沿第 56 轮符号约定，令

\[
J=\begin{pmatrix}0&-1\\1&0\end{pmatrix}=-iY,\qquad
\mathcal J_d=J\otimes I_d,\qquad \mathcal J_d^2=-I.
\]

它是一个声明的共同参考结构。尚未证明认知组合或任务增益必然提供它。

考察与 $\mathcal J_d$ 对易的实对称状态。逐块解对易关系得到

\[
\omega=\begin{pmatrix}A&B\\-B&A\end{pmatrix},
\quad A^{\mathsf T}=A,\quad B^{\mathsf T}=-B.
\]

定义复矩阵 $\rho=2(A+iB)$，则 $\omega\succeq0,\operatorname{Tr}\omega=1$ 当且仅当 $\rho\succeq0,\operatorname{Tr}\rho=1$。反向编码为

\[
\mathcal E(\rho)=\tfrac12R(\rho),\quad
R(M)=I\otimes\operatorname{Re}M-J\otimes\operatorname{Im}M.
\]

因此这个受限实整体的齐次维数恰为 $d^2$，等价于复 d 维密度矩阵。编码中的额外两维是参考自由度，不是新增时空坐标；复纯态编码为实秩二态。

对任意复效应 F，

\[
\operatorname{Tr}_{\mathbb R}[\mathcal E(\rho)R(F)]
=\operatorname{Tr}_{\mathbb C}(\rho F).
\]

概率保持来自 $R(MN)=R(M)R(N)$ 及实迹为复实部迹的两倍。此处是表示与受限状态族，不把分化定义成删除整体状态。

## 2 哪些操作保持这个取向

任意与 $\mathcal J_d$ 对易的实矩阵都是 R(K)。再要求正交，便等价于 $K^\dagger K=I$：

\[
\{O\in O(2d):[O,\mathcal J_d]=0\}=R(U(d)).
\]

同理，若实 Kraus 算符每个都与共同结构对易，保迹条件化为 $\sum_kK_k^\dagger K_k=I$，并且

\[
\sum_kR(K_k)\mathcal E(\rho)R(K_k)^{\mathsf T}
=\mathcal E\!\left(\sum_kK_k\rho K_k^\dagger\right).
\]

这给出复态、效应和完全正通道的准确有效接口。但“每个 Kraus 算符保持共同取向”强于“输出状态还在同一集合”，不能自动替换成后者。

一般维数的对应是矩阵层面的证明，没有新增任意维数旧门的有限编译。d=2 的相关旧实控制与提升可复用第 50、56、69 轮。

## 3 保持状态集合还不够

取实正交矩阵

\[
Q=\operatorname{diag}(I_d,-I_d).
\]

它满足 $Q\mathcal J_dQ^{\mathsf T}=-\mathcal J_d$，但依然把每个允许编码态送回允许编码态：

\[
Q\mathcal E(\rho)Q^{\mathsf T}=\mathcal E(\rho^*).
\]

因此仅要求保持同一个状态集合，仍容许整体取向反转，对应逻辑复共轭。

这不是可在任意复纠缠背景中自由调用的局部通道。若将复共轭的线性扩展——转置——只作用于 Bell 态的一半，反对称投影的期望为 −1/2。不存在负概率问题的实际实 Q 操作，作用的是整个公共参考编码；把它误当作独立复子系统的局部通道，才错误改变了组合条件。

所以“同一结构 + 能力扩展”尚未自动给出复 CPTP 理论：共同结构的存在、取向保持及其跨子系统使用方式都需要依据。

## 4 来源与解释边界

该实表示与对易约束是已有方法。本轮读取 [Real-Vector-Space Quantum Theory with a Universal Quantum Bit，第 II 节](https://arxiv.org/pdf/1210.4535)，核对共同 J、块矩阵、状态的 1/2 因子及对易分解。我们的符号沿第 56 轮；全局 J 的正负选择只是约定。

没有复现该论文以环境动力学逼近对易约束的后续模型，也没有将其特定极限当作认知推导。这里的新工作是：把“保持状态集合”和“保持取向的可组合操作”分开，避免从相似结构直接跳到量子理论唯一性。

## 5 交付

[common_orientation_structure.py](common_orientation_structure.py)；[结果 JSON](common_orientation_structure_results.json)。8 项检查覆盖状态锥维数、正性和解码、效应概率、控制与 Kraus 提升、实取向翻转及 Bell 部分转置反例。

    python -X utf8 research_cognition_physics/common_orientation_structure.py --write-results

下一轮检查两个独立编码能否直接当成同一共同参考下的复合系统。
