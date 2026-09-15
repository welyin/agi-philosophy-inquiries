# 第十一轮：公开执行模式后的条件预测与连续事件记录

日期：2026-09-15。接续第九、十轮，固定读取方向为 0、纵向参数为 1。

## 1 问题与接口

第十轮发现：只公开微弱二元结果时，记录区分力可以趋零，而状态仍发生有限改变。现在公开每一步究竟是 idle（保留输入并产生公平位）还是 active（实际读取并重新准备），以及该步二元结果。

这是同一个经典核的记录接口扩充。这里的“来源”只指执行原语；没有公开初始准备标签、潜变量 x 或输出窗口中的随机位置，也没有增加输入备份。idle 位与来源无关，可从有效记录中删去。

## 2 条件更新

沿用 h=1/4、α=sin(h)/h、κ=√(1−α²)。归一化预测摘要为 z=(1,u,v)，u²+v²≤α²。用 s=±1 表示 active 结果：

\[
A_s=\frac12\begin{pmatrix}
1&s\kappa&0\\
\alpha s\kappa&\alpha&0\\
0&0&\alpha^2
\end{pmatrix},\qquad p_s(z)=\frac{1+s\kappa u}{2}.
\]

已知 active 及结果后的预测是

\[
\Phi_s(z)=\left(1,\frac{\alpha(u+s\kappa)}{1+s\kappa u},
\frac{\alpha^2v}{1+s\kappa u}\right).
\]

已知 idle 后 z 不变。对于一般未知初态，以上是给定先验的前向预测规则；观察一个事件并不会告诉观察者精确的初态摘要。

每步累计强度增量 δτ≤κ，以 δτ/κ 的概率 active。因此把模式重新隐藏就恢复第九轮 B_s(δτ)，没有改变底层核。

## 3 连续事件极限

把累计无量纲强度 τ 分成 n 步。active 计数服从 Binomial(n,τ/(κn))，极限是 Poisson(τ/κ)。不相交段的模式选择独立，因此得到率 1/κ 的泊松事件过程；等待强度间隔服从相应指数分布。这里使用标准稀有事件极限，泊松定义、计数与等待分布可对照 [MIT 概率课程第 22 讲](https://ocw.mit.edu/courses/res-6-012-introduction-to-probability-spring-2018/d92a00a2bc9d20d84ad11d43be0e7ae0_MITRES_6_012S18_L22.pdf)。本轮针对这些内容核对讲义。

给定当前预测，带结果的强度为

\[
\lambda_s(z)=\frac{1+s\kappa u}{2\kappa},\qquad
dz=\sum_{s=\pm1}(\Phi_s(z)-z)\,dN_s.
\]

事件之间状态不变。N_s 是由状态相关结果标记的计数过程；两个标记过程不能额外假定为相互独立、常强度的泊松过程。

对 0<t₁<⋯<t_k<τ，指定有序标记的路径密度为

\[
e^{-\tau/\kappa}\kappa^{-k}
e_0^T A_{s_k}\cdots A_{s_1}z_0\,dt_1\cdots dt_k.
\]

对标记求和、在有序时间单纯形上积分，就得到泊松计数概率。固定方向时，给定计数的事件位置分布不依赖输入，全部输入信息在有序 active 标记中。

忽略记录，S=Σ_s A_s=diag(1,α,α²)，平均演化满足

\[
\bar z(\tau)=\exp[\tau(S-I)/\kappa]z_0
=\operatorname{diag}(1,e^{-(1-\alpha)\tau/\kappa},e^{-\kappa\tau})z_0.
\]

这与第十轮的均值极限相同。τ 是累计操作参数；这一步没有定义物理时间或时空对称性。

## 4 全部 active 标记能多保留多少信息

比较输入 μ₀ 与 μ_π。令 D_k 为恰有 k 次 active 时，全部 2^k 个有序标记串分布的总变差距离（TV）。因后续记录可被丢弃，D_{k+1}≥D_k。完整事件记录距离为

\[
D(\tau)=\sum_{k=0}^{\infty}e^{-m}\frac{m^k}{k!}D_k,
\qquad m=\tau/\kappa.
\]

脚本精确枚举到 k=18；这不是轨迹蒙特卡洛。截断以后只用 0≤D_k≤1。若 K+2>m，余项上界是

\[
\sum_{k>K}p_k\le\frac{p_{K+1}}{1-m/(K+2)}.
\]

| 累计强度 τ | 仅第一次 active 的 TV | 全部 active 标记的 TV |
|:---|---:|---:|
| 0.1 | 0.0713046972 | 约 0.0736098929 |
| 0.3 | 0.1246006715 | 约 0.1511520309 |
| 1 | 0.1421096033 | 0.2793107 至 0.2794325 |

最后一行的解析截断误差上界为 0.000121699。表中端点已向外取整；浮点累加并非区间算术认证。第十轮“隐藏模式时区分力趋零”的结果不受影响，因为两个结论使用不同的公开接口。

## 5 本轮结果与接续

完整执行记录确实可以提高关于初始来源的统计信息，并给出自洽的经典跳跃过程。它还没有证明能够恢复原状态。下一轮将检查：A_s 的代数可逆性，是否对应一个对未知输入普遍合法的恢复操作。

代码：[marked_event_limit.py](marked_event_limit.py)。结果：[marked_event_limit_results.json](marked_event_limit_results.json)。新增 10 项检查通过，包括记录归一化、条件更新、模式粗粒化、均值极限、路径密度积分和截断余项。

```text
python -X utf8 research_cognition_physics/marked_event_limit.py --write-results
```
