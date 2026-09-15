# 第六十六轮：固定两源网络，写出完整统计与噪声分数

日期：2026-09-15。接续 [第 65 轮](research_note_65.md)。

## 1 本轮改变了什么

现在比较的是同一个两源网络的同一份记录，来源为 A–B₁ 与 B₂–C，不另加待测共同态。两个来源独立；Alice、Charlie 分别选本地设置 x、z，Bob 没有设置输入、输出四种标签之一 b。各方测量不接收其他方的消息，设置与准备独立，全部结果保留。最后把记录汇总计算分数是允许的。

本轮先条件性采用复矩阵状态、通常张量组合与迹概率规则。它们仍是候选模型的输入。

这种用两源纠缠交换区分实、复矩阵模型的框架已有文献：Renou 等的 [Quantum theory based on real numbers can be experimentally falsified](https://arxiv.org/pdf/2101.10873)，图 2、正文网络公式及三组 CHSH 分数。本项目在此使用已有网络，接上旧门与原噪声，不声称发现了这一分离现象。

## 2 来源与 Bob 标签

两份独立准备为

\[
\rho=\Phi^+_{AB_1}\otimes\Phi^+_{B_2C},\qquad
\Phi^+=|\Phi^+\rangle\langle\Phi^+|,
\quad |\Phi^+\rangle=(|00\rangle+|11\rangle)/\sqrt2.
\]

它们均为实态，可以用已有实准备和一次旧 YX(π/2) 流生成。Bob 用该门的逆作用后读取两个 Z 标签。因此理想 Bell 基就是旧门的四列；允许的相关符号如下。

| b | X 相关 q₁ | Y 相关 q₂ | Z 相关 q₃ |
|:---|---:|---:|---:|
| 0 | 1 | −1 | 1 |
| 1 | 1 | 1 | −1 |
| 2 | −1 | −1 | −1 |
| 3 | −1 | 1 | 1 |

每行均满足 q₁q₂q₃=−1。下文数学用 1、2、3 表示 X、Y、Z；代码索引用 0、1、2。

Alice 测 Pᵢ；Charlie 对每个 i<j 测 (Pᵢ+sPⱼ)/√2，s=±1，共六种设置。完整概率为

\[
p(a,b,c\mid x,(i,j,s))
=\frac1{16}\left[1+\frac{ac\,q_x^b}{\sqrt2}
(\delta_{xi}+s\delta_{xj})\right].
\]

3×6×2×4×2=288 项都已保存。每组设置归一化，p(b)=1/4；忽略 b 后两端完全均匀。因此差异存在于联合记录中，不能从单边摘要读取。

## 3 一个不丢结果的线性分数

令

\[
E^b_{x,z}=\sum_{a,c=\pm1}ac\,p(a,b,c\mid x,z),
\]

注意它没有除以 p(b)。定义

\[
T=\sum_b\sum_{i<j}
\left[q_i^b(E^b_{i,ij+}+E^b_{i,ij-})
+q_j^b(E^b_{j,ij+}-E^b_{j,ij-})\right].
\]

代入完整概率，得到 **T=6√2≈8.48528137424**。每个 b 的贡献已经带着该结果的实际发生概率；不是把四份条件分数直接相加。

## 4 原读取噪声怎样进入

设两端最终读取可见度分别为 β_A、β_C，有限私人 Y 参考偏置为 ν_A、ν_C。两端实际方向为

\[
(X,Y,Z)\longmapsto\beta_L(X,\nu_LY,Z),\qquad L=A,C.
\]

Charlie 的和差方向使用同一线性替换。Bob 两次 Z 标签读取可见度均为 γ_B；其效应为

\[
B_{uv}=U_{YX}(\pi/2)
\left[\frac{I+(-1)^u\gamma_BZ}{2}\otimes
\frac{I+(-1)^v\gamma_BZ}{2}\right]U_{YX}(\pi/2)^\dagger.
\]

按上表，q_X=(−1)^u，q_Y=−(−1)^v，q_Z=(−1)^{u+v}。所以 X、Y 相关衰减 γ_B，Z 相关衰减 γ_B²。精确分数为

\[
\boxed{T=2\sqrt2\,\beta_A\beta_C
\left[\gamma_B(1+\nu_A\nu_C)+\gamma_B^2\right].}
\]

这也解释了为什么不能把 Bob 的整个 Bell 测量统一乘一个可见度。若参考先理想化、四个目标都只读一次原仪器，β_A=β_C=γ_B=α=4sin(1/4)，得到 T≈8.19521382227。第 69 轮将把理想参考和读取都换成有限流程。

两份来源若独立退极化，可见度为 v_L、v_R，则上式再乘 v_Lv_R。这是该噪声模型的公式，并未覆盖任意控制偏差。

## 5 验证与接续

- [脚本](bell_network_statistics.py)：9 项检查通过，包括 288 项完整概率、旧 Bell 门、四输出噪声译码与非各向同性分数。
- [结果](bell_network_statistics_results.json)：保存全部理想记录和诊断值。
- [第 67 轮](research_note_67.md)：固定逐来源实参考方案，允许 Bob 使用任意测量，求其最优值。

本轮没有测量真实设备，也没有从认知要求推出来源独立、通信限制或新增 Y 能力。
