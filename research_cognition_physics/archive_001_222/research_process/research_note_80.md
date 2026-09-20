# 第八十轮：组合后仍是同类认知结构，能否选择复数域

日期：2026-09-15。接续 [第 79 轮](research_note_79.md)。用户进一步提出：**整体采用复结构，会不会是因为组合以后必须复现原本的认知结构？**

本轮检验一个较弱而明确的版本：整体和子系统使用同一类状态、操作及组合规则，允许容量不同，改换组合顺序不改变预测。结果是：**完整实矩阵模型也满足这种递归一致性。**

## 1 同一类实结构怎样递归组合

全部实矩阵分成对称部分 S 与反对称部分 A：

\[
M_d(\mathbb R)=S_d\oplus A_d,\qquad
K(d)=\frac{d(d+1)}2,\quad L(d)=\frac{d(d-1)}2.
\]

K、L 分别为两个空间的维数。局部实密度矩阵属于 S；A 不是额外可准备的实状态。转置作用于每个因子，故

\[
S_{ab}=(S_a\otimes S_b)\oplus(A_a\otimes A_b),\qquad
A_{ab}=(S_a\otimes A_b)\oplus(A_a\otimes S_b).
\]

因此

\[
(K,L)\star(M,N)=(KM+LN,KN+LM).
\]

三个块的对称部分由 SSS、SAA、ASA、AAS 构成，与先合并哪两个无关。矩阵张量积的结合律也保证逐状态、逐操作的描述一致，不只是维数相等。

| 实块维数 d | 对称空间 K | 反对称空间 L |
|:---|:---|:---|
| 2 | 3 | 1 |
| 4 | 10 | 6 |
| 8 | 36 | 28 |
| 16 | 136 | 120 |

K 包含归一化坐标，L 描述组合时的关系方向；它们不是时空维数。上述组合式与任意划分下的实成对层析已有文献，见 [Hardy–Wootters，第 4.2—4.4 节](https://arxiv.org/pdf/1005.4870)。本轮实际重读式 (31)—(33) 及对称、反对称基的证明，将它用于用户的递归结构问题，不声称该规律原创。

## 2 任意大小的块都需要完整操作接口

实 Kraus 操作 $\Phi(M)=\sum_jK_jMK_j^{\mathsf T}$ 分别保持 S 和 A，完整接口为

\[
\widehat\Phi=
\begin{pmatrix}B_\Phi&0\\0&C_\Phi\end{pmatrix}.
\]

d=2 时 C 就是第 52 轮的标量 κ；d=4 时，B 为 $10\times10$，C 为 $6\times6$。它们受完全正性及相应保迹条件限制，不是任意自由参数。

两操作的 B、C 均相同，则在全部实矩阵单位上相同，与任意辅助组合后也相同。接口的复合和混合分别按矩阵乘法和凸组合执行。因此实模型可以一致地把多个个体看成新的整体；只是块自身的实状态统计未必足够校准其完整操作身份。

## 3 较大块也确实存在缺口

对任意 d≥2，构造两个保迹映射

\[
\Phi_\pm(M)=\operatorname{Tr}(M)\frac Id
\pm\frac{M-M^{\mathsf T}}{2d^2}.
\]

它们在 S 上相同，在 A 上分别为 $\pm A/d^2$，且都完全正。其未归一化 Choi 矩阵为

\[
J_\pm=I_{d^2}/d
\pm\frac{|\Omega\rangle\langle\Omega|-F}{2d^2},
\quad\Omega=\sum_i|ii\rangle.
\]

F 是交换矩阵。括号除以 2 的特征值在 Ω、其正交对称空间和反对称空间上分别为 $(d-1)/2,-1/2,+1/2$。所以

\[
J_\pm\ge\left(\frac1d-\frac{d-1}{2d^2}\right)I>0.
\]

实正 Choi 矩阵的谱分解给出实 Kraus，代码逐矩阵单位核对。

对任意非零反对称 A，只加一个实二能级辅助，定义

\[
\rho_A=\frac{I_{2d}+\frac12(A/\|A\|)\otimes J}{2d},
\quad J=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\]

它严格正且归一化。两映射输出之差非零，故能被整体实效应区分。更一般地，任何在 S 上相同、在某个 A 上不同的操作对，都由这个探针暴露区别。**一个额外 rebit 足以检测任意大小实块的某个缺失接口方向。** 没有给出任意维数的最少原门或样本成本。

## 4 哪种更强的复现要求会排除它

若要求组合后的全部预测状态由局部实验的联合记录确定，还需局部记录映射没有核。实模型缺口为

\[
K(ab)-K(a)K(b)=L(a)L(b)>0\qquad(a,b\ge2).
\]

复模型的齐次维数为 $d^2$，满足乘法；经典 d 状态概率模型的维数为 d，也满足乘法。所以这个更强条件也不能单独选择复量子理论。

“恢复全部状态预测”和“分开后执行所有原共同操作”又是不同要求。[第 83 轮](research_note_83.md) 将进一步审计共享参考在这里的作用。

## 5 交付

[partition_interface_closure.py](partition_interface_closure.py)；[结果 JSON](partition_interface_closure_results.json)。新增 8 项检查，覆盖整数基、张量结合、完整操作重建、较大块实 CP 反例、辅助探针和接口复合混合。矩阵与张量规则仍为候选输入，未从认知推出。

    python -X utf8 research_cognition_physics/partition_interface_closure.py --write-results
