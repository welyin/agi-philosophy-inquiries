# 第一百二十七轮：保留旧群体的确定性近似接入

日期：2026-09-16。沿第 126 轮继续，允许共同输出与目标有非零误差，所有分支都接受。本轮给出构造；第 128 轮检验全局最优，第 129 轮检验双方怎样分担误差。

**结果：旧群体的未知逻辑状态可以完全保留。新的两能级主体确定性接入后，其状态变为 $2\sigma/3+I/6$，最坏迹距离误差为 1/6。**

## 1 输入、目标和误差

输入仍是单份独立标准实编码

\[
\Omega=\mathcal E(\rho_A)\otimes\mathcal E(\sigma_B),\qquad
\mathcal E(\rho)=\operatorname{real\_lift}(\rho)/2.
\]

ρ_A 可以是旧群体内任意维数的未知状态，σ_B 是 d 维未知状态。目标为 $\mathcal E(\rho_A\otimes\sigma_B)$。不增加未知态副本，不读取其经典描述。

误差采用迹距离

\[
D(\omega,\tau)=\tfrac12\|\omega-\tau\|_1.
\]

它界定任一效应的概率差，并界定对同一输入执行相同后续物理协议时的完整记录总变差。这里比较指定共同接口的输出；额外环境记录的访问范围必须另行说明。

## 2 已有的近似转置通道

使用

\[
\widetilde T_d(X)=\frac{\operatorname{tr}(X)I+X^{\mathsf T}}{d+1}.
\]

这不是新的量子通道：其最优近似转置性质及 Kraus 实现见 [Buscemi 等，第 2 节，式 (15)—(20)](https://arxiv.org/pdf/quant-ph/0304175)；量子设计实现见 [Kalev 与 Bae，式 (1)—(2)](https://arxiv.org/pdf/1303.3096)。本轮只把这一已有工具接入本项目的独立编码任务，不能直接把单系统转置最优性当作整个接入过程的最优性。

可取全部为实的 Kraus 矩阵：

\[
A_{ii}=\sqrt{\frac2{d+1}}|i\rangle\langle i|,
\qquad
A_{ij}=\frac{|i\rangle\langle j|+|j\rangle\langle i|}{\sqrt{d+1}}
\quad(i<j).
\]

直接相加得到 $\sum A_{ij}^{\mathsf T}A_{ij}=I$ 及上述通道公式。共有 d(d+1)/2 个 Kraus 标签。代码对矩阵单位逐项核对线性映射，不只测试几份密度矩阵。

## 3 完整接入通道

第 125 轮两分支分别为

\[
\omega_+=\tfrac12\mathcal E(\rho_A\otimes\sigma_B),\qquad
\omega_-=\tfrac12\mathcal E(\rho_A\otimes\sigma_B^*).
\]

正分支原样接受，负分支在 B 上施加 $\widetilde T_d$。σ 为 Hermitian 矩阵，故 $(\sigma^*)^{\mathsf T}=\sigma$，从而负分支变为

\[
\widetilde\omega_-=
\tfrac12\mathcal E\!\left(\rho_A\otimes\frac{\sigma_B+I}{d+1}\right).
\]

完整实 Kraus 集合是 $K_+$ 和

\[
(I_{R A}\otimes A_{ij})K_-.
\]

其伴随积之和等于 I，故不需要成功筛选。取当前共同输出边缘：

\[
\boxed{\omega_{\rm out}=\mathcal E(\rho_A\otimes\Lambda_d(\sigma_B)),\quad
\Lambda_d(\sigma)=\frac{(d+2)\sigma+I}{2(d+1)}.}
\]

ρ_A 精确保留，包括其内部已纳入 A 的记忆与纠缠。这里不声称对任意未纳入输入模型的外部关联都获得相同保证。

## 4 状态和旧能力损失

实编码保持 Hermitian 差的迹范数，且 $\|\rho_A\otimes X\|_1=\|X\|_1$，所以

\[
D(\omega_{\rm out},\mathcal E(\rho_A\otimes\sigma))
=\frac{\|I-d\sigma\|_1}{4(d+1)}
\le\frac{d-1}{2(d+1)}.
\]

最后一步由谱分解和凸性得到，每个纯 σ 都达到上界。

| 新主体维数 d | 本构造最坏误差 |
|:---|:---|
| 2 | 1/6 |
| 3 | 1/4 |
| 4 | 3/10 |
| 8 | 7/18 |

这些是本构造的界；本轮尚未证明全部 d 的一般最优。

对于 d=2，

\[
\Lambda_2(\sigma)=\tfrac23\sigma+\tfrac16I.
\]

任意两份新主体状态的区别收缩为原来的 2/3：

\[
D(\Lambda_2(\sigma),\Lambda_2(\tau))=\tfrac23D(\sigma,\tau).
\]

因此原先相互正交、可完全判别的两种状态，在此共同接口上的最优等先验判别成功率变为 5/6。旧群体 A 的局部能力则保持；B 的任意固定后续仪器仍能执行，但其原预测不能全部精确复现。对同一份已接入输出的任意有限后续协议，记录误差由迹距离收缩性控制；这不是无限份新接入任务的累计误差界。

## 5 保留全部分支，不虚报条件精度

保留粗取向标记 h 时，可与具有同一个均匀 h 的理想目标比较：

\[
\omega_{\rm flagged}=
|0\rangle\langle0|\otimes\omega_+
+|1\rangle\langle1|\otimes\widetilde\omega_-.
\]

纯新主体的条件误差分别为 0 和 1/3，平均及带 h 的联合迹距离为 1/6。因此**不能给每条历史都贴上“误差至多 1/6”**。

更细的 Kraus 标签 ij 对输入有条件依赖。它们保留于环境中；本轮没有声称读取并条件化这些细标签后，所有条件态或任意新日志比较也满足同一 1/6 界。h 或细标签不是丢弃失败数据的理由。

## 6 “整体还在”的物理实现账

把完整 Kraus 列表写成等距嵌入

\[
V=\sum_r |r\rangle_E\otimes L_r,\qquad V^{\mathsf T}V=I.
\]

V 将输入内容保存在活跃输出和环境的整体中，代码核对了 $V^{\mathsf T}(V\Omega V^{\mathsf T})V=\Omega$。取环境边缘表示该接口暂不访问环境，不是物理删除。

两量子位例子的输入维数为 16，活跃共同输出维数为 8，四个 Kraus 标签用一个四维环境保存。V 是 32×16 实等距矩阵，可用**一个额外纯 rebit**加理想实正交延拓实现；这是充分容量，未证明全部实现的最小值。环境包含两位容量，其中也承载取向分支信息；复制经典 h 另计记录空间。

纯态初始化、保持环境相干、联合控制及分布式访问都有成本。这里没有把 V 编译为项目的有限噪声原门和原读取，也没有沿用第 124 轮的任务门数。

若反向访问整个环境，可以撤销接入并恢复旧输入。若要求在撤销之外仍输出更准确的指定共同编码，该解码又是一个输入到共同输出的通道，须受下一轮的一般下界约束。

## 7 复核与接续

[脚本](approximate_subject_joining.py)；[结果](approximate_subject_joining_results.json)。8 项检查覆盖转置通道的完整线性作用、全部 Kraus 的迹保持、未知混态公式、旧群体精确边缘、纯态达到误差、标记条件精度、三轴区别收缩和整体等距保存。

    python -X utf8 research_cognition_physics/approximate_subject_joining.py --write-results

接续：[第 128 轮](research_note_128.md)检验 1/6 是否为任意通道的必要界；[第 129 轮](research_note_129.md)检验能否让双方分担。

三轮合并后，目录全部 **1226 项检查通过**，耗时 143.341 秒；旧脚本与旧结果保留。
