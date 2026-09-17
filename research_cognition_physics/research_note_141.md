# 第一百四十一轮：经典历史恢复外部关联的锐界是 1/2

日期：2026-09-16。承接第 138—140 轮，将保护目标从主体自身的边缘态提高为它与外部系统已有的未知关联。

**结果：第一次固定平分接入之后，只开放三位共同接口和任意经典环境记录，最优最坏外部恢复误差恰为 1/2。第 138 轮方案已经达到这个界。换测量基底、增加经典消息、采用不同条件恢复，都不能突破。**

这是指定访问模型中的条件定理，不是认知原则或复数结构的推导。

## 1 任务与量词

原两主体标准独立编码所在空间 S 的维数为 d=16。第一次最优平分接入后，共同接口 A 的维数为 q=8；其余四位环境 E 全部仍在。

允许对 E 作任意测量，将任意多的经典结果交给恢复器；恢复器对 A 作任意依结果选择的 CPTP 操作，可以加入独立准备的辅助系统。恢复输出 S' 与原 S 同维。

不向恢复器开放源纯化 R，不提供额外环境到恢复器的量子传输，也不免费加入跨环境端与恢复端的预共享纠缠。环境测量后，未开放的量子余部不再反馈。这一类包含此前的有限含噪经典记录方案，也放宽了原实门和实际读数限制。

目标误差为

\[
\epsilon_{\rm ext}=\inf_{\text{允许恢复}}
\sup_{\omega_{SR}:\ \omega_S=E(\rho_A)\otimes E(\rho_B)}
D\bigl((\mathcal C\otimes\mathrm{id}_R)(\omega_{SR}),\omega_{SR}\bigr).
\]

其中 D(ρ,σ)=‖ρ−σ‖_1/2。经典记录保留在整体；上式比较恢复系统与原外部系统的联合边缘，没有要求物理删除记录。

## 2 所有经典测量与恢复都受同一个秩界约束

把任意环境测量细分，得到接入仪器的 Kraus 算子 K_a：16→8。给定报告 a，恢复器的 Kraus 算子记为 L_{aj}：8→16。完整恢复的每个算子

\[
A_{aj}=L_{aj}K_a,\qquad \operatorname{rank}A_{aj}\le8.
\]

辅助系统、额外局部随机性和粗记录都包含在 Kraus 索引中；独立辅助不会提高这个因子分解的中间维数。

对任意方阵 A，奇异值三角不等式与 Cauchy–Schwarz 给出

\[
|\operatorname{Tr}A|\le\|A\|_1
\le\sqrt{\operatorname{rank}A}\,\|A\|_2.
\]

令 \(|\Phi_d\rangle=d^{-1/2}\sum_i|i\rangle|i\rangle\)。恢复后的目标重叠满足

\[
F_\Phi=\frac1{d^2}\sum_{a,j}|\operatorname{Tr}A_{aj}|^2
\le\frac q{d^2}\sum_{a,j}\operatorname{Tr}A_{aj}^\dagger A_{aj}
=\frac qd=\frac12.
\]

这里 F 是与纯目标的重叠，未取平方根。测量目标投影便有 D≥1−F，故任何允许恢复都满足 ε_ext≥1/2。

见证输入的边缘 I_16/16=E(I_2/2)⊗E(I_2/2)，完全符合独立编码承诺；该纯化也可解释为两个主体分别与各自外部记忆的纯化，随后仅重排坐标。外部记忆从未交给解码器。

## 3 旧经典恢复恰好达到

第 138 轮完整经典记录的恢复通道为

\[
\Delta(\omega)=P_+\omega P_++P_-\omega P_-.
\]

定义 Q=P_+−P_-，有 Q²=I、Q 实正交，且

\[
\Delta=\tfrac12\mathrm{id}+\tfrac12\operatorname{Ad}_Q.
\]

对任意外部系统和任意联合态，

\[
D((\Delta\otimes\mathrm{id})(\omega),\omega)
=\tfrac12D((Q\otimes I)\omega(Q\otimes I),\omega)\le\tfrac12.
\]

最大纠缠见证达到等号。因此

\[
\boxed{\epsilon_{\rm ext}^{\rm classical}=1/2}.
\]

下界甚至允许复 CPTP 处理；达到方案只用既有实操作。更换数域本身不会移除当前输入和访问权限造成的限制。

## 4 对“整体还在”的意义

主体自身的承诺状态仍能完全恢复，因为它原本没有两取向块间的相干项。但外部关系可以包含这些项；恢复边缘无法检验它们是否仍可访问。

完整源、接入环境、测量装置与记录仍在整体中。定理限制的是仅有经典反馈的恢复器；重新开放量子历史或提供额外共享资源会改变任务，不能再直接套用该界。

量子环境测量与条件纠错的框架参见 [Gregoratti–Werner，第 II 节](https://arxiv.org/pdf/quant-ph/0209025)。秩限制与最大纠缠重叠的关系也见 [Terhal–Horodecki，Lemma 1](https://arxiv.org/pdf/quant-ph/9911117)。本轮使用标准工具完成当前模型的匹配上下界，没有把这些原理当作新发现。

## 5 复算与推进

- 脚本：[external_correlation_recovery_bound.py](external_correlation_recovery_bound.py)。
- 结果：[external_correlation_recovery_bound_results.json](external_correlation_recovery_bound_results.json)。
- 8 项检查覆盖非正规复矩阵秩界、实际恢复通道、其他实与复环境基底、合法纯化见证、完整联合态误差、实反射表示、独立辅助和维数条件。

随机基底检查不是全局优化证据；普遍下界来自解析秩证明。下一轮推广到多主体与相干记忆容量。
