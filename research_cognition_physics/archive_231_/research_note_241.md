# 第241轮：从有限量子操作接入历史框架——记录究竟改变什么

日期：2026-09-22。阶段：231起的物理生成研究。

## 1. 接续问题、已有工作与假设

[第240轮](research_note_240.md)停止调参常规CSG的整体二维取样分支，提出转向量子历史。这里首先解决接口：已有CP仪器能否直接给出归一化退相干泛函？“删除出生记录”究竟指忽略、重置、后选择还是测量后反馈？

这接续[第225轮](../archive_223_230/research_note_225.md)已明确留下的相干实现相位问题，不重做量子理论的重建。[第193轮](../archive_001_222/research_process/research_note_193.md)已经说明，概率开平方的坐标替换不等于量子化。

直接采用已有框架：

- Sorkin，[Quantum Mechanics as Quantum Measure Theory](https://arxiv.org/html/gr-qc/9401003)：二阶干涉的量子测度。
- Dowker、Johnston、Sorkin，[Hilbert Spaces from Path Integrals](https://arxiv.org/html/1002.0589)，第2节及第3.2节：强正退相干泛函、历史Hilbert空间及给定酉模型的构造。
- Dowker、Halliwell，[Quantum mechanics of history](https://doi.org/10.1103/PhysRevD.46.1580)：历史的退相干与概率求和规则。

**待检验假设：** CP操作可接入该框架，但仪器完备性本身不确定相干历史权重；记录的实际可访问性会改变干涉与恢复能力。

**认知动机：** “共同认可历史”需要明确哪些记录可被实际访问和操作。**额外输入：** 以下给定干涉线路、历史分解、记录耦合及反馈权限；它们尚未由认知合同选定。

## 2. 采用的数学接口及一个归一化反例

设有限历史的类算子为C_h，初态为密度矩阵ρ。采用第一变量共轭线性的约定：

$$
D(h,k)=\operatorname{Tr}(C_h^\dagger C_k\rho),\qquad
D(A,B)=\sum_{h\in A,\,k\in B}D(h,k),\qquad
\mu(A)=D(A,A).
\tag{1}
$$

C_h√ρ在Hilbert–Schmidt空间中的Gram矩阵就是D，因此强正；事件向量是相应列向量之和。若所有历史算子之和为酉算子U，则：

$$
\sum_h C_h=U
\quad\Longrightarrow\quad
D(\Omega,\Omega)=\operatorname{Tr}(U^\dagger U\rho)=1.
\tag{2}
$$

式(2)是充分条件，不是所有退相干泛函的必要表示。二次型结构使互不相交A、B、C满足：

$$
\begin{aligned}
I_3(A,B,C)
={}&\mu(A\cup B\cup C)-\mu(A\cup B)-\mu(A\cup C)\\
&-\mu(B\cup C)+\mu(A)+\mu(B)+\mu(C)=0.
\end{aligned}
\tag{3}
$$

这是已有量子测度规则，展开后逐项相消即可证明。本轮只核对它如何接入已获得的操作结构。

**CP仪器的完备性不是式(2)。** 取两个带不同经典输出的Kraus算子K₀＝K₁＝I／√2，则：

$$
\sum_a K_a^\dagger K_a=I,\qquad
D_{\rm bare}=\frac12
\begin{pmatrix}1&1\\1&1\end{pmatrix},\qquad
\operatorname{Tr}D_{\rm bare}=1,\quad
D_{\rm bare}(\Omega,\Omega)=2.
\tag{4}
$$

将这两个Kraus算子直接当作式(1)的历史类算子会失败。保留正交输出记录得到对角D、总权重1；或另行指定相干实现。简单把式(4)整体除以2会改变原仪器的分支权重，不能冒充原操作的无代价重写。

## 3. 四条历史与记录重叠

目标初态为|+〉，路径投影P₀、P₁，重组门H为Hadamard门，最终输出端口x∈{0,1}。定义：

$$
C_{a,x}=P_xHP_a,\qquad
\sum_{a,x}C_{a,x}=H,\qquad
|\psi\rangle=|+\rangle.
\tag{5}
$$

按(0,0)、(1,0)、(0,1)、(1,1)排列，直接得到：

$$
D_0=\frac14
\begin{pmatrix}
1&1&0&0\\
1&1&0&0\\
0&0&1&-1\\
0&0&-1&1
\end{pmatrix}.
\tag{6}
$$

端口0两条历史的并集测度为1，而对角之和只有1／2。还存在三历史子集的测度5／4。**量子测度不是对任意历史提问都能直接使用的普通概率分布。** 实际最终端口构成退相干的划分，可使用Born概率；任意非退相干子集的μ不能自动这样解释。

现在路径a写入归一化标记态|r_a〉。令γ＝〈r₀|r₁〉，则：

$$
D_\gamma((a,x),(b,y))
=D_0((a,x),(b,y))\langle r_a|r_b\rangle,\qquad
p(x=0)=\frac{1+\operatorname{Re}\gamma}{2}.
\tag{7}
$$

γ＝0给D＝I₄／4，端口各半；γ＝1恢复式(6)。Gram秩由2变为4，是因为把记录系统纳入了总系统，不是目标量子比特或物理空间突然增加维数。对γ＝0、0.4、1、0.5i，代码核验强正、归一化与式(3)。

## 4. 四种“抹除”必须分开

### 4.1 忽略或局部重置记录

任意目标A、参考R、记录M的联合态σ，以及只作用于M的迹保持通道Λ，都满足：

$$
\operatorname{Tr}_M
[(\operatorname{id}_{AR}\otimes\Lambda_M)(\sigma)]
=\operatorname{Tr}_M\sigma.
\tag{8}
$$

证明：对任意AR效果E，两边期望相同，因为Λ*把恒等效果映为恒等效果。因此只操作记录而不按结果处理目标，不能改变目标及其参考的边缘态。

所以“我不看它”和“把它重置成0”都不能恢复本例干涉。后者若通过环境实现，路径信息可以转移到环境；不能把局部重置等同于逆转此前的联合耦合。

### 4.2 测量记录并后选择

对正交标记，联合编码为V|ψ〉＝Σ_a P_a|ψ〉⊗|a〉。测量标记的X基，诱导目标算子：

$$
K_+=\frac{I}{\sqrt2},\qquad
K_-=\frac{Z}{\sqrt2},\qquad
p_+=p_-=\frac12.
\tag{9}
$$

对初态|+〉，选择“+”子样本可恢复亮端口，选择“−”子样本则得到暗端口。两者不加区分地合并仍各端口一半。这里没有丢掉后选择的成功概率。

### 4.3 测量后反馈可以确定性恢复

若“−”结果后允许在目标施加Z，“+”后不动，则：

$$
L_+=K_+=\frac{I}{\sqrt2},\qquad
L_-=ZK_-=\frac{I}{\sqrt2},\qquad
\sum_s L_s\rho L_s^\dagger=\rho.
\tag{10}
$$

这是完整恒等通道，故也保持目标与任意外部参考的未知纠缠。恢复并非必须后选择；它要求标记仍可相干访问，并实际对目标反馈。式(8)不禁止式(10)，因为这里不再只操作标记。

若另一个不可访问环境保留正交路径副本，本协议只能给Z退相干通道。对于Bell态输入，与原态的迹距离为1／2。这个计算是该协议在指定噪声后的失败见证；它不声称任何环境、任何编码均不可恢复。此例目标到可访问系统的通道已经丢失输入的非对角项，因而无法对所有未知输入作通用逆转。

## 5. 可复算结果

[代码](quantum_history_record_audit.py)、[结果](quantum_history_record_audit_results.json)、[证据检查](research_round_241_checks.json)、[公式检查](round241_math_checks.json)。

| 检查对象 | 结果 |
|---|---|
| 无标记／忽略正交标记／仅重置标记，端口0概率 | 1／1⁄2／1⁄2 |
| 选择“+”记录，端口0条件概率／成功率 | 1／1⁄2 |
| 两种结果均接受、目标反馈后，端口0概率 | 1 |
| 相干标记反馈，Bell参考恢复的迹距离 | 5.56×10⁻¹⁷以下 |
| 额外不可访问正交副本，同一恢复的迹距离 | 1⁄2 |
| 四历史所有256种互斥三事件分配，I₃最大残差 | 2.23×10⁻¹⁶以下 |

8项检查通过。有限数值负责检验矩阵实现与反例；一般结论由式(2)、(8)、(10)等证明，不由采样替代。

    python -B -X utf8 research_cognition_physics/archive_231_/quantum_history_record_audit.py
    python -B -X utf8 research_cognition_physics/archive_231_/verify_history_rounds.py 241

## 6. 结论与下一步

**接口成立，但条件必须写全：** 给定相干类算子与记录模型，既有量子操作可以产生规范化强正历史泛函；普通CP仪器完备性本身不能选定它。

认知解释上，“共同拥有记录”是可操作的物理条件，不能用主观知晓与否替代。以上没有选出宇宙的自然增长权重，也未导出时空或引力。

下一轮处理独立问题：同一无标签历史的不同自然编号，在什么条件下可以视为重复描述？通道层面的交换是否足够？接续225轮的相干相位缺口与235轮的可交换接口，避免把编号重复计数制造的权重变化当作真实干涉。
