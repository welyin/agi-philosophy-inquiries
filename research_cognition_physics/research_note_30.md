# 第三十轮：精确旋转把 0.144739 上界推广到任意准备非情境编码

日期：2026-09-15。接续第二十九轮。仍固定 h=1/4、α=sinc(h)，研究全部精确静默旋转和原 λ=1 fresh 读取。

## 1 主要结论

第二十九轮在一个固定编码下得到的条件

\[
\boxed{\alpha G(\eta)\le1,\qquad
G(\eta)=\sqrt{1-\eta^2}+\eta\arcsin\eta}
\]

实际上约束**所有**尊重真实准备混合、条件准备与全部精确旋转的准备非情境模型。隐藏变量可以连续、有任意额外标签，不预设它就是一个单位角度。

因此一般可行强度上限满足

\[
\eta_*\le\eta_B\approx0.144739232095,
\qquad \alpha G(\eta_B)=1.
\]

η_* 指同一模型允许全部 η≤η_max、所有方向及任意有限组合时，可行 η_max 的上确界。此结果把第二十七轮的 0.17 上界收紧；没有证明 η_B 能达到。

## 2 一般准备编码的仿射系数

用 t=(u,v)/α 表示归一化准备圆盘。沿用第二十二轮的仿射表示论证，准备非情境性加真实混合要求

\[
d\mu_t(\lambda)=[1+a(\lambda)\cdot t]\,d\mu_0(\lambda),
\qquad |a(\lambda)|\le1,\qquad \int a\,d\mu_0=0.
\]

a 是准备变化在隐藏概率中的两个系数，以下称为“准备系数向量”。这里不是额外假定密度光滑：对任意可测集合，准备概率关于 t 仿射；中心准备支配所有 μ_t，Radon–Nikodym 系数给出 a。对圆盘所有方向的准备要求正性，便得到 |a|≤1。

特别地，不能未经论证就令 |a|=1，或把不同半径的隐藏编码都视为第二十九轮的 F。

## 3 精确旋转迫使系数向量的分布旋转不变

令 R_θ 为平面旋转，T_θ 为任意实现该静默操作的随机核。它必须满足

\[
T_\theta\mu_t=\mu_{R_\theta t},\qquad T_\theta\mu_0=\mu_0.
\]

在联合分布 μ_0(dλ)T_θ(dλ′|λ) 下比较 t 的系数：

\[
E[a(\lambda)\mid\lambda']=R_\theta^T a(\lambda').
\]

两边输入、输出边缘都是 μ_0，而 R_θ 保持模长。因此

\[
\begin{aligned}
E|a(\lambda)-R_\theta^Ta(\lambda')|^2
&=E|a(\lambda)|^2-E|R_\theta^Ta(\lambda')|^2\\
&=0.
\end{aligned}
\]

于是 a(λ′)=R_θa(λ) 几乎处处。完整隐藏状态的更新仍可以随机，但它的准备系数向量只能精确旋转。

由于每个 θ 都是允许操作，a 在 μ_0 下的分布旋转不变。因此它是半径 r∈[0,1] 的均匀圆分布的混合，记半径分布 κ(dr)；允许 r=0 的原子。这里没有要求隐藏空间上的核构成群，也没有给隐藏标签预设旋转对称性。

## 4 选择性更新给出的 Jensen 不等式

以 (w,t_x,t_y) 作为齐次坐标时，读取方向为零的 fresh 矩阵为

\[
B_s=\frac12
\begin{pmatrix}
1&s\alpha\eta&0\\
s\eta&\alpha&0\\
0&0&\alpha q
\end{pmatrix}.
\]

操作后的条件准备仍在同一圆盘里，准备非情境性要求整个输出隐藏测度等于这一状态的编码，而不只是读取概率相同。于是均匀中心准备输入后的结果—隐藏状态边缘为

\[
f_s(a')\,d\mu_0(\lambda'),\qquad
f_s(a')=\frac{1+s\eta a'_x}{2}.
\]

比较准备系数，得到反向条件均值

\[
E[a_{\rm in}\mid s,\lambda']
=b_s(a')
=\frac{\alpha(a'_x+s\eta,\ q a'_y)}{1+s\eta a'_x}.
\]

η<1 时分母严格为正。对 x 投影的绝对值使用条件 Jensen 不等式，再把结果相加：

\[
\begin{aligned}
E_{\mu_0}|a_x|
&\ge\sum_s\int f_s(a)|b_{s,x}(a)|\,d\mu_0\\
&=\frac\alpha2 E[|a_x+\eta|+|a_x-\eta|]\\
&=\alpha E\max\{|a_x|,\eta\}.
\end{aligned}
\]

这条必要条件从一般模型直接得出，没有用单位圆编码。

## 5 消除任意半径分布

在半径 r 的均匀圆上，

\[
E|a_x|=\frac{2r}{\pi},\qquad
M(r,\eta):=E\max\{|a_x|,\eta\}
=\begin{cases}
\eta,&r\le\eta,\\
\dfrac2\pi\left[\sqrt{r^2-\eta^2}
+\eta\arcsin(\eta/r)\right],&r>\eta.
\end{cases}
\]

对 r>0，把积分中的 r 提出，M(r,η)/r 是 E max{|cos φ|,η/r}。它随 η/r 单调增加，故对 r≤1，

\[
M(r,\eta)\ge\frac{2r}{\pi}G(\eta).
\]

若 E_κr>0，代入上一节立即得到 αG(η)≤1。若 E_κr=0，则 a=0 几乎处处，原 Jensen 不等式变成 0≥αη，排除任何 η>0。两种情况覆盖全部编码。

所以 η>η_B 不可能存在这样的准备非情境模型。对于允许全部较低强度的上限问题，这也排除了把 η_max=1 当作例外。

## 6 严格的 0.145 排除证书

浮点求根给出 η_B 的位置，严格排除另用有理数验证。α 的交错 Taylor 下界取第二十七轮已有证书

\[
\alpha_L=1-\frac1{96}+\frac1{30720}-\frac1{4096\cdot5040}<\alpha.
\]

由 G′(η)=arcsin η 的正系数级数，

\[
G(\eta)=1+\frac{\eta^2}{2}+\frac{\eta^4}{24}
+\frac{\eta^6}{80}+\cdots.
\]

取 η=29/200=0.145，精确分数运算得到

\[
\alpha_L\left(1+\frac{\eta^2}{2}+\frac{\eta^4}{24}\right)-1
=\frac{29648629958351}{792723456000000000}
>3.74009\times10^{-5}>0.
\]

因此该点的排除不依赖浮点精度、角网格或有限停止策略的优劣。η_B 本身仍由解析等式定义，约数不能替代等式。

## 7 等号还要求什么

若有模型在 η_B 实现读取，上一节的各步必须取等号。对于 η_B>0，r<1 时

\[
\alpha M(r,\eta_B)-2r/\pi>0.
\]

因此所有质量必须落在 r=1：|a|=1 几乎处处，其角分布均匀。条件绝对值 Jensen 取等号又要求：给定每个 (s,λ′)，输入 a_x 不能同时有正、负的非零质量；若条件均值为零，则 a_x=0 几乎处处。

这给端点构造增加了具体限制，但尚不证明这些限制能够同时满足。它也不意味着在低于端点的所有模型中都能提前假设 |a|=1。一般隐藏变量纤维及其他凸函数的约束仍应谨慎处理。

## 8 复算与下一轮

代码：[rotation_noncontextual_bound.py](rotation_noncontextual_bound.py)。结果：[rotation_noncontextual_bound_results.json](rotation_noncontextual_bound_results.json)。8 项检查通过，包括一般半径积分、原分支系数对照、严格有理数余量、等号半径限制和旋转方差控制例。一般证明在本笔记中，有限网格控制例不替代证明。

~~~
python -X utf8 research_cognition_physics/rotation_noncontextual_bound.py --write-results
~~~

下一轮回到正面构造：尝试用小角度对称转移代替大量重新采样，利用第二十九轮的条件均值接口提高可行下界。准备非情境性及完整条件准备保持仍是明确的附加要求；这里没有排除保留不可访问准备差异的原经典模型，也没有推出量子理论。
