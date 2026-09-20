# 第十四轮：缩小准备宽度后，经典模型也允许有信息的二阶弱扰动

日期：2026-09-15。接续第十三轮，只改变准备宽度的标度；这是模型敏感性检验，不是原固定宽度问题的解。

## 1 先说明改变了什么

第十三轮证明固定 h>0 时，全方向非选择扰动不能是 O(η²)。系数 1−sin(h)/h 来自准备规则。现在考虑一族实验：第 n 个实验运行 n 步，各步都用同一个准备宽度 h_n，原语为第七轮 λ=1 的 fresh 更新，且

\[
h_n=\sqrt{\nu/n},\qquad \eta_n=\sigma/\sqrt n,
\qquad \alpha_n=\frac{\sin h_n}{h_n}.
\]

ν≥0、σ>0 固定；n 足够大，使 η_n≤1、h_n<π/2。一个实验内部不改变准备集合；n 不同表示模型不同。ν=0 需要显式允许点准备及其混合，不属于原来的正宽度窗口集合。

这是人为增加的分辨率资源：实验越长，准备可以越精细。现有认知原则没有推出这种标度，也没有提供其成本模型。

## 2 每个模型都有实际经典实现

每一步使用

\[
A_{s,n}=\frac12\begin{pmatrix}
1&s\eta_n&0\\
\alpha_n s\eta_n&\alpha_n&0\\
0&0&\alpha_n\sqrt{1-\eta_n^2}
\end{pmatrix}.
\]

其点态核以 f_s(x)=(1+sη_n cosx)/2 的分支权重，输出一个半宽 h_n 的窗口，中心为

\[
\phi_s(x)=\operatorname{atan2}
\left(\sqrt{1-\eta_n^2}\sin x,\cos x+s\eta_n\right).
\]

零权重点的中心可任意选。任意输入都输出原宽度窗口的非负混合，完整准备保持成立；窗口积分直接给出上述矩阵。所有执行都是 active，只公开二元结果，无需增加执行模式信息。

## 3 二阶单步扰动与非平凡均值极限

忽略结果，N_n=diag(1,α_n,α_n√(1−η_n²))。令 a=ν/6，则

\[
\alpha_n=1-\frac{a}{n}+O(n^{-2}),\qquad
\alpha_n\sqrt{1-\eta_n^2}
=1-\frac{a+\sigma^2/2}{n}+O(n^{-2}).
\]

因此 N_n−I=O(1/n)=O(η_n²)，且

\[
N_n^n\longrightarrow
\operatorname{diag}(1,e^{-a},e^{-a-\sigma^2/2}).
\]

ν=σ=1 时两个均值增益分别为 **0.8464817249、0.5134171190**，有限且非零。固定 h=1/4 的同一 fresh 家族则含 α^n→0；改变准备标度的影响不能忽略。

这里的二阶扰动专指忽略结果后的更新接近恒等；每个已知结果的条件变化一般仍是 O(η_n)。第十二轮关于精确恢复的结论没有被推翻。

## 4 用全部二元记录的解析下界证明信息没有消失

比较第 n 个模型的 μ₀ 与 μ_π，初态分别 (1,±α_n,0)。令 s_j=±1，记录统计量

\[
Y_n=\frac1{\sqrt n}\sum_{j=1}^n s_j.
\]

由均值矩阵与结果差矩阵 J_n=A_{+,n}−A_{-,n} 直接得到

\[
\mathbb E_\pm s_j=\pm\eta_n\alpha_n^j,
\quad
\mathbb E(s_i s_j)=\eta_n^2\alpha_n^{j-i}\quad(i<j).
\]

第二个等式对这两个输入相同，且无需独立结果假设。因此

\[
\mathbb E_\pm Y_n=\pm m_n,\qquad
m_n=\frac{\sigma\alpha_n}{n}
\sum_{j=0}^{n-1}\alpha_n^j
\longrightarrow \sigma\frac{1-e^{-a}}{a}.
\]

a=0 时将最后比值连续延拓为 1。二阶矩满足

\[
\mathbb E_\pm Y_n^2
=1+\frac{2\eta_n^2}{n}\sum_{d=1}^{n-1}(n-d)\alpha_n^d
\le1+\sigma^2\frac{n-1}{n}=:M_n.
\]

设 P_n、Q_n 是全部二元串的分布。对它们相对于共同控制测度的密度，用 Cauchy–Schwarz：

\[
|\mathbb E_PY-\mathbb E_QY|^2
\le (\mathbb E_PY^2+\mathbb E_QY^2)
\int\frac{(p-q)^2}{p+q}
\le4M\,\operatorname{TV}(P,Q).
\]

代入均值差 2m_n，得到有限 n 的界 TV(P_n,Q_n)≥m_n²/M_n，继而

\[
\boxed{\liminf_n\operatorname{TV}(P_n,Q_n)
\ge\frac{\sigma^2}{1+\sigma^2}
\left(\frac{1-e^{-a}}a\right)^2>0.}
\]

ν=σ=1 时下界为 **0.4242214943**。这是全部公开二元记录的信息下界，不是某次模拟的显著性，也不是精确 TV 的数值。它只使用直接算出的记录矩，因此不依赖下面尚未完成的扩散极限定理。

| 步数 n | 准备半宽 h_n | 纵向均值增益 | 横向均值增益 | 全记录 TV 下界 |
|:---|---:|---:|---:|---:|
| 100 | 0.1 | 0.8464346696 | 0.5120981105 | 0.4256199243 |
| 1000 | 0.0316227766 | 0.8464770219 | 0.5132858434 | 0.4243606825 |
| 10000 | 0.01 | 0.8464812546 | 0.5134039977 | 0.4242354066 |

## 5 条件过程的候选连续形式

在归一化状态 (1,u,v) 上，对一步条件变化直接做有限结果求和与展开，每单位步参数的漂移、协方差趋向

\[
b(u,v)=(-au,-(a+\sigma^2/2)v),\qquad
G(u,v)=\sigma(1-u^2,-uv),\quad GG^T.
\]

这与下列候选随机微分方程的局部矩一致：

\[
du=-au\,dt+\sigma(1-u^2)\,dW,\qquad
dv=-(a+\sigma^2/2)v\,dt-\sigma uv\,dW.
\]

W 是候选创新噪声，t 是缩放的操作步数。脚本核对了条件漂移与协方差的收敛。**本轮没有完成完整记录与状态过程的路径弱收敛证明**，所以这些式子暂作为候选扩散过程；不能只靠一步局部矩一致就宣称全过程等价。

## 6 结论与下一步

本轮给出一个完整经典反例：调整准备分辨率，就能让每步非选择扰动为二阶，同时累计二元记录保有非零输入信息。因此“有效弱读取”也不能独立作为量子性的判据。原固定宽度一般横向优化问题保持未解；这里更换了准备规则。

下一轮优先补上联合记录—状态的路径极限，检验候选过程的边界保持、创新噪声与可识别参数；随后比较不同 ν 的过程能否仅由现有认知原则选定。若仍有任意噪声参数，就明确指出欠确定性，进一步寻找可观测而非形式类比的约束。

代码：[shrinking_preparations.py](shrinking_preparations.py)。结果：[shrinking_preparations_results.json](shrinking_preparations_results.json)。新增 10 项检查通过，包括独立输出窗口积分、完整二元树记录矩、信息下界、均值极限及条件局部矩。

```text
python -X utf8 research_cognition_physics/shrinking_preparations.py --write-results
```
