# 第二十五轮：夹住旧显式核的可行范围，区分估计余量和真正失效

日期：2026-09-15。只分析第二十四轮的 Poisson 加性修正核，不据此限制其他实现。

## 1 结果

原先 η≤0.00174982 的充分范围可以扩大到 **η≤0.004048809588**。但同一个核在 η=0.0041 时已出现严格负值。因此旧界确实保守，同时这个具体构造也确实很快失效，不能只靠继续细化网格无限改善。

这没有排除更好的准备非情境更新核。

## 2 利用两个角度之间的几何关系

设输入角 x、输出角 y、c=cos(y−x)、t=1−√(1−η²)。旧核乘以 4π 后为

\[
Q=P_\alpha(y-x)+s\eta(\cos y+2\alpha\cos x)-2\alpha t\sin y\sin x.
\]

固定 c 时，线性项幅度及正弦积满足

\[
|\cos y+2\alpha\cos x|\le B(c):=\sqrt{1+4\alpha^2+4\alpha c},\qquad
\sin y\sin x\le(1+c)/2.
\]

因此 Q≥P(c)−ηB(c)−αt(1+c)。精确求线性项比值最小值：

\[
\min_{-1\le c\le1}\frac{P(c)}{B(c)}
=\ell:=\frac{1-\alpha^2}{(1+2\alpha^2)^{3/2}}.
\]

证明：令 D=1+α²−2αc，最大化 D²B²；其导数为 −12α²D(α+2c)，唯一内部最大点是 c=−α/2。

另有 (1+c)/B(c)≤2/(1+2α)，可直接求导验证。于是

\[
Q\ge B(c)\left[\ell-\eta-\frac{2\alpha}{1+2\alpha}t\right].
\]

用 t≤η²，令 d=2α/(1+2α)，得到充分条件 η+dη²≤ℓ，即

\[
\boxed{\eta\le\frac{2\ell}{1+\sqrt{1+4d\ell}}=0.004048809588\ldots}.
\]

这是对全部连续角度的证明，不是从有限角网格中挑出最小值。

## 3 真正负值的明确见证

取 s=+1、y=3π/2、x=5π/6，则

\[
Q=\frac{1-\alpha^2}{1+\alpha+\alpha^2}-\sqrt3\alpha\eta+\alpha(1-\sqrt{1-\eta^2}).
\]

在 η=41/10000，利用 α≥95/96、α≤95/96+1/30720、√3>433/250，以及 1−√(1−η²)≤η²，得到纯有理数上界

\[
Q<-0.00002985<0.
\]

代码以精确分数验证了上界为负。因此该核不能实现包含 η=0.0041 的强度族。

还有一个更接近线性最坏角度的点见证，其浮点零点约为 0.004064501546。该零点只属于一个指定角对，没有证明是全核的全局失效阈值。本文的严格结论是全局充分界加明确负值点，而非全局精确最优性。

## 4 接续

旧构造用一阶强度修正消耗 Poisson 核的正余量。下一轮改用随结果变换的角度核，再补偿归一化偏差，尝试将消耗降至二阶。状态编码、原读取矩阵和准备混合规则保持同样定义。

代码：[kernel_positivity_bounds.py](kernel_positivity_bounds.py)。结果：[kernel_positivity_bounds_results.json](kernel_positivity_bounds_results.json)。8 项检查通过，包含全局比值极值、几何界、改进正性证书和有理数负值见证。

```text
python -X utf8 research_cognition_physics/kernel_positivity_bounds.py --write-results
```
