# 第242轮：相干过程与自然编号——哪些历史可以合并

日期：2026-09-22。阶段：231起的物理生成研究。

## 1. 本轮问题与已有依据

[241轮](research_note_241.md)接上强正历史框架，但没有决定历史标签的物理意义。本轮检验：**两个编号给相同CP通道，是否足以把它们视为同一个相干历史？**

[225轮](../archive_223_230/research_note_225.md)已经指出：张量重结合的共轭映射相容，不自动固定实现子的相干控制相位。[235轮](research_note_235.md)要求真实联合接口与可交换性；[239轮](research_note_239.md)区分无标签取样与带出生编号的分布。本轮把这些已有缺口连到量子历史中，不将相位问题当作新发现。

直接使用的既有结果：

- Araújo等，[Quantum circuits cannot control unknown operations](https://arxiv.org/html/1309.7976)：未知操作的普通黑箱说明不足以任意添加相干控制。
- Araújo、Costa、Brukner，[Computational advantage from quantum-controlled ordering of gates](https://arxiv.org/html/1401.8127)：量子控制操作次序可以区分交换与反交换。
- Chiribella等，[Quantum computations without definite causal structure](https://arxiv.org/abs/0912.0195)：量子开关的高阶操作框架。

本地贡献是将这些已知结构用于本路线的编号审计，并给出一个有限偏序上的充分条件及完整复算例。**认知动机：** 更换描述编号应保持所有已获授权的测量关系。**额外输入：** 明示事件偏序、具体分支算子、相干控制实现；不从“可融入”直接假定任意高阶连接权限。

## 2. 相同普通通道不固定相干控制

U和带整体相位的U给出相同目标通道：

$$
\operatorname{Ad}_{e^{i\phi}U}
=\operatorname{Ad}_{U},\qquad
\operatorname{Ad}_U(\rho)=U\rho U^\dagger.
\tag{1}
$$

但在给定控制寄存器后：

$$
W_\phi
=|0\rangle\langle0|\otimes I
+|1\rangle\langle1|\otimes e^{i\phi}I,\qquad
p(+\mid |+\rangle)=\frac{1+\cos\phi}{2}.
\tag{2}
$$

φ＝0、π／2、π分别得到1、1／2、0。普通目标通道始终是恒等通道。式(2)的可观测相位属于扩大后的联合操作说明；不是同一完整物理操作给出了互相矛盾的预测。

第二阶段的有限门权限足以任意逼近这些已知联合矩阵；R_seed下可精确实现。它并不从一个未指定实现的黑箱U自动选出其中哪一个。这里正是225轮保留的说明层次差别。

### 2.1 同一组合通道，仍可能有不同次序干涉

采用已有量子开关的酉形式：

$$
S(U,V)
=|0\rangle\langle0|\otimes VU
+|1\rangle\langle1|\otimes UV.
\tag{3}
$$

取U＝X、V＝Z。因ZX＝−XZ，两种组合的普通目标通道相同；相干控制的|+〉却变成|−〉。作为对照，U＝V＝X给|+〉。因此“两个约化通道交换”不足以把两个可控制次序宣称为纯编号差异。

**这不把同一目标上的X、Z称为类空事件。** 它是检查过弱判据的反例。真实独立张量因子的算子严格交换，见下一节。

### 2.2 不能把受控U的不定性错误推广给量子开关

开关两分支各使用U、V一次，故分别改变两门整体相位，只给S乘共同相位：

$$
S(e^{i\alpha}U,e^{i\beta}V)
=e^{i(\alpha+\beta)}S(U,V).
\tag{4}
$$

其联合通道不变。更一般地，给A、B两通道的Kraus族，标准开关使用：

$$
S_{ij}
=|0\rangle\langle0|\otimes B_jA_i
+|1\rangle\langle1|\otimes A_iB_j.
\tag{5}
$$

两族分别作酉混合后，S族按张量积酉矩阵混合，因而表示同一通道。代码同时核对迹保持与Kraus基变换不变性。

这里使用的是已有且明示的高阶构造。没有宣称普通固定次序黑箱线路、每个未知门只调用一次，就已经实现了量子开关。

## 3. 有限自然编号独立性的充分条件

给定有限事件偏序P，给每个事件及其物理结果指定同一总Hilbert空间上的分支算子K_{e,a}。要求不可比事件具有严格的算子交换关系：

$$
e\parallel f
\quad\Longrightarrow\quad
K_{e,a}K_{f,b}=K_{f,b}K_{e,a}
\quad\text{对所有允许结果 }a,b.
\tag{6}
$$

类型改变或带多时隙记忆时，需要先给相容的共同接口；这里的命题限于式(6)已经有定义的有限分支。

**命题。** 对固定事件结果赋值h，任一线性扩展σ得到的历史乘积均相同：

$$
C_{h,\sigma}
=K_{\sigma(n),a_{\sigma(n)}}\cdots
K_{\sigma(1),a_{\sigma(1)}}
=C_h.
\tag{7}
$$

**证明。** 任意两个有限偏序线性扩展可以通过交换相邻的不可比元素连接。构造上，取目标排序的首元素，在当前排序中将其左移：它前面的元素既不可能是它的前驱，也不可能是它的后继，因此均不可比。首元素对齐后，对余下元素归纳。每次相邻交换由式(6)保持算子乘积，于是式(7)成立。

因此在物理结果历史空间上，241轮的D可以用任意一个线性扩展计算；无需按该历史具有多少种自然编号额外乘权重。此处是**充分条件**，不声称所有编号独立模型必须具有这种表示，也不声称由此得到广义相对论的协变性。

### 3.1 六种编号的非平凡正对照

采用两条独立事件链A₀≺A₁、B₀≺B₁，其余不可比。A链用H、Z；B链用实旋转R(0.37)、X。链内操作一般不交换，跨链因处于不同张量因子而严格交换。

六个线性扩展的整体算子全部一致。进一步将A₁换成投影仪器的结果P_a，将B₁换成比特翻转仪器结果√q_b X^b，q＝(0.7,0.3)。四种结果赋值各自在六个编号中保持相同分支算子，且整个结果族满足仪器完备性。

这是给定局域接口的编号无关性，不是预先给定偏序已经从认知原则中涌现。

## 4. 重复描述不应伪装成路径振幅

两个物理结果的正确概率均为1／2。若同一结果分别使用1个、2个标签重复书写，直接加标签概率后再归一化会变成(1／3,2／3)；直接加相同相位的标签振幅再平方则变成(1／5,4／5)：

$$
p^{\rm count}_h
=\frac{m_hp_h}{\sum_km_kp_k},\qquad
p^{\rm amplitude}_h
=\frac{m_h^2p_h}{\sum_km_k^2p_k}.
\tag{8}
$$

这是把重复描述当新增物理路径造成的偏置；不是所有路径积分的统一反例，也不是对某个既有因果集测度的替代公式。选择一次代表可避开本例重复；真正的带标签测度还需明确其商映射和权重合同。

若m种调度确实由物理寄存器相干控制，就必须包括归一化的寄存器态。满足式(7)时：

$$
\begin{aligned}
|\mathrm{u}_m\rangle&=\frac1{\sqrt m}\sum_{\sigma=1}^m|\sigma\rangle,\\
\left(\sum_\sigma|\sigma\rangle\langle\sigma|\otimes C_h\right)
(|\mathrm{u}_m\rangle\otimes|\psi\rangle)
&=|\mathrm{u}_m\rangle\otimes C_h|\psi\rangle.
\end{aligned}
\tag{9}
$$

没有额外m倍振幅。无选择地读取或丢弃该调度寄存器也不制造m²权重。六调度酉例在纠缠输入上直接核验式(9)，统一叠加态的回读概率为1。

## 5. 证据与范围

[代码](242/coherent_label_audit.py)、[结果](242/coherent_label_audit_results.json)、[证据检查](242/research_round_242_checks.json)、[公式检查](242/round242_math_checks.json)。代码复用241轮基本矩阵工具。

| 核验 | 结果 |
|---|---|
| 目标通道相同，φ＝0、π／2、π的控制“+”概率 | 1、1⁄2、0 |
| X、Z两次序的普通通道差／开关“+”概率 | 0／0 |
| 两通道的Kraus酉混合前后，开关Choi矩阵最大差 | 8.35×10⁻¹⁷以下 |
| 两链六种编号的酉乘积最大差 | 0 |
| 带选择结果的历史乘积最大差 | 2.78×10⁻¹⁷以下 |
| 归一化六调度寄存器的因子化误差 | 3.47×10⁻¹⁸以下 |

8项检查通过。本轮一般命题是式(6)⇒式(7)，其余量子控制机制均按既有工作使用和复算。未选定增长权重、空间维数或物理作用量。

    python -B -X utf8 research_cognition_physics/archive_231_/coherent_label_audit.py
    python -B -X utf8 research_cognition_physics/archive_231_/verify_history_rounds.py 242

## 6. 下一步：真正进入量子增长候选

241—242轮接口已足够用于下一步，不继续重复干涉玩具实验。下轮优先核对已有复顺序增长模型，在明确的有限历史树上计算强正、归一化和跨步相容：

$$
D_n(h,k)
=\sum_{h'\succ h,\,k'\succ k}D_{n+1}(h',k').
\tag{10}
$$

这里h′、k′各比父历史多一步。首先区分经典记录模型与复振幅模型，再按文献实际条件处理编号及可观测事件；不假定有限层式(10)已保证无限历史上的测度存在。

下轮指定入口为Dowker、Johnston、Surya的[On extending the Quantum Measure](https://arxiv.org/abs/1007.2725)，以及Surya、Zalel的[A Criterion for Covariance in Complex Sequential Growth Models](https://arxiv.org/abs/2003.11311)。本轮只核对其摘要所述的延拓问题及研究目标，尚未应用其具体判据；下一轮先读正文，随后复用公式。两轮接口的完成不意味着已有自然动力学或时空生成证明。
