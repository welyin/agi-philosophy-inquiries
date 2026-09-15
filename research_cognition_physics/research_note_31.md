# 第三十一轮：小角度循环转移把完整实现下界提高到 0.11

日期：2026-09-15。接续第二十九轮的条件均值接口及第三十轮的一般上界。仍固定原准备宽度与 λ=1 fresh 仪器。

## 1 主要结果

构造连续隐藏角上的正 Markov 仪器，**完整实现全部 η≤0.11 的 fresh 读取、全部方向和任意有限自适应组合**。核包含每个观察结果的四个确定性输出原子，以及一个很小但不能删除的连续项。

全角度非负性由解析界和精确有理数证书保证。结合第三十轮，

\[
\boxed{0.11\le\eta_*\le\eta_B\approx0.144739232095},
\qquad
\alpha[\sqrt{1-\eta_B^2}+\eta_B\arcsin\eta_B]=1.
\]

上一阶段的区间约为 [0.08346824,0.17]。新下界是充分构造，上界是必要条件，中间仍未解决。

## 2 从反向条件分布设计核

沿用第二十九轮的 n_ψ、R_η、f_s、ψ_s、T_s。需要一个反向条件核 P(dx|ψ)，使

\[
\int n_x P(dx\mid\psi)=\alpha n_\psi,\qquad
\int R(\psi)P(dx\mid\psi)\frac{d\psi}{2\pi}=\frac{dx}{2\pi}.
\]

直观上，第一式限制每次转移保留多少方向信息，第二式要求把略不均匀的角质量搬成均匀分布。相比大幅重新采样，小角度搬运可以用较少的一阶矩损失调整局部质量。

取 N=12、d=π/N、c=cos d。令未归一化模式质量 A、U、C、V 都是 ψ 的 π 周期函数，定义

\[
P(dx\mid\psi)=\frac1{R(\psi)}
\left[
A(\psi)\delta_\psi(dx)
+U(\psi)\delta_{\psi+d}(dx)
+U(\psi)\delta_{\psi-d}(dx)
+C(\psi)\delta_{\psi+\pi}(dx)
+V(\psi)\frac{dx}{2\pi}
\right].
\]

只要

\[
\begin{aligned}
C&=\frac{(1-\alpha)R-V}{2}-(1-c)U,\\
A&=\frac{(1+\alpha)R-V}{2}-(1+c)U,
\end{aligned}
\]

便有 A+2U+C+V=R、A+2cU−C=αR；对称的两次偏移又使切向均值相消。剩下的任务是选择 U、V，使输入边缘均匀且全部质量非负。

## 3 十二项循环 Green 函数

定义离散差分

\[
(LU)(\psi)=2U(\psi)-U(\psi+d)-U(\psi-d),
\]

以及循环平均

\[
H(\psi)=\frac1N\sum_{j=0}^{N-1}R(\psi+jd).
\]

由于 R 是 π 周期函数，这里的有限循环用于解函数方程；ψ 仍连续，并没有把隐藏空间换成十二个标签。

取

\[
G_j=\frac{N^2-1}{12N}-\frac{j(N-j)}{2N},
\qquad
U_p(\psi)=\sum_{j=0}^{N-1}G_j[R(\psi+jd)-1].
\]

直接差分得到 ΣG_j=0，2G_j−G_{j-1}−G_{j+1}=δ_{j0}−1/N，下标按 N 取模。于是

\[
LU_p=R-H.
\]

还需保留 L 无法改变的谐波。第二十六轮已给出

\[
R(\psi)=1+2\sum_{k\ge1}c_k\cos(2k\psi),\quad
c_k=t^{2k}(1+2kq)>0,\quad
t=\frac{\eta}{1+q}.
\]

H 保留其中 k 为 N 倍数的项。记 a=t^{2N}、z=a e^{2iNψ}，则

\[
H-1=2\operatorname{Re}\left[\frac z{1-z}+
\frac{2Nqz}{(1-z)^2}\right],
\quad
V_0=2\left[\frac a{1-a}+\frac{2Nqa}{(1-a)^2}\right].
\]

令 V=V_0+H−1，由正 Fourier 系数有 0≤V≤2V_0，平均 V=V_0。该闭式对所有 η≤0.11 都有效，没有截断无限级数。

## 4 选取非负偏移并证明均匀边缘

对 k 不是 N 倍数的谐波，L 的特征值为 4sin²(kd)。U_p 的对应余弦振幅为 c_k/[2sin²(kd)]。因为 sin²(kd)≥sin²d，

\[
|U_p|\le U_0:=
\frac{R_{\max}-1}{4\sin^2d},\qquad
R_{\max}=\frac{1+\eta^2}{q}.
\]

于是取 U=U_0+U_p，得到 0≤U≤2U_0，且 LU=R−H。

反向核的输入角边缘密度，相对于 dx/(2π)，为

\[
\begin{aligned}
h(x)
&=A(x)+C(x-\pi)+U(x-d)+U(x+d)+\overline V\\
&=R(x)-V(x)-LU(x)+V_0\\
&=R-[V_0+H-1]-(R-H)+V_0=1.
\end{aligned}
\]

这同时证明边缘匹配和总核的逐输入归一化，不依赖角度网格精度。

## 5 全区间非负性及有理数证书

R_min=q³，且 0≤U≤2U_0、0≤V≤2V_0。因此足够要求

\[
\begin{aligned}
2C&\ge(1-\alpha)q^3-2V_0-4(1-c)U_0\ge0,\\
2A&\ge(1+\alpha)q^3-2V_0-4(1+c)U_0\ge0.
\end{aligned}
\]

在 η=0.11，直接浮点代入给出这两个全角度下界分别约为 0.0008981798、1.41715646。严格证书进一步用有理数包住每一项。

令 e=11/100，使用第二十七轮的 α_L<α<α_U，以及

\[
q_L=1-\frac{e^2}{2}-\frac{e^4}{2},\quad
\ell=\frac{6698}{100000}<\sin^2(\pi/12),\quad
c_L=\frac{9659}{10000}<c<\frac{966}{1000}=c_U.
\]

q_L≤√(1−η²) 对全部 η≤e 成立：将 1−x/2−x²/2 平方，与 1−x 比较即可验证 x∈[0,1] 的下界。三角界由 sin²(π/12)=(2−√3)/4、cos(π/12)=(√6+√2)/4 及有理数平方比较取得，代码逐项核对。

再取

\[
U_U=\frac{(1+e^2)/q_L-1}{4\ell},\quad
a_U=\left(\frac e{1+q_L}\right)^{2N},\quad
V_U=2\left[\frac{a_U}{1-a_U}
+\frac{2Na_U}{(1-a_U)^2}\right].
\]

这些界同时适用于所有 η≤e；V_U 使用 q≤1，所以不必依赖 Fourier 系数对 η 的单调性。精确分数运算得到

\[
\begin{aligned}
(1-\alpha_U)q_L^3-2V_U-4(1-c_L)U_U
&>0.00085985>0,\\
(1+\alpha_L)q_L^3-2V_U-4(1+c_U)U_U
&>1.41510>0.
\end{aligned}
\]

完整分子、分母写入结果 JSON。这是整个连续角域、整个强度区间的证书；有限采样只作交叉检查。

## 6 写出正向可执行仪器

令偏移 j 取 {0,d,−d,π}，对应质量 W_j 为 {A,U,U,C}。给定实际输入隐藏角 x，结果 s 的四个原子为

\[
\begin{aligned}
\psi&=x-j,\\
p_{s,j}(x)&=
\frac{q^3}{2(1-s\eta\cos\psi)^2}
\frac{W_j(\psi)}{R(\psi)},\\
y&=T_s(\psi).
\end{aligned}
\]

连续部分与输入 x 无关，密度为

\[
k_s^{\rm reset}(y)=
\frac{f_s(y)}{2\pi}\,
\frac{V(\psi_s(y))}{R(\psi_s(y))}.
\]

所以

\[
\boxed{K_s(dy\mid x)=
\sum_j p_{s,j}(x)\delta_{T_s(x-j)}(dy)
+k_s^{\rm reset}(y)dy.}
\]

连续项的合计质量为 V_0，原子总质量逐点为 1−V_0。给定 x 即可在八个原子及连续模式间抽样；若进入连续模式，使用上述联合结果—输出密度除以 V_0。η=0 时 V_0=0，无须进入该模式。模式属于隐藏实现，不额外公开。

η=0.11 时 V_0≈3.13914×10^−29。它虽非常小，却负责循环差分不能处理的谐波，**不能为了简化把它删掉后仍称精确实现**。代码直接计算 H−1 的闭式；没有用接近 1 的浮点平均再减去 1 来求它。

其他读取方向 θ：输入、输出都先相对 θ 表示，再将输出角加回 θ。静默旋转用 x↦x+θ。

## 7 完整条件分布，而非仅低阶矩

正向核来自第二十九轮的联合分布翻译。其反向条件均值已严格等于 αn_ψ，因此对任意合法 z，

\[
\begin{aligned}
(K_sF_z)(y)
&=\frac{f_s(y)}{2\pi}
\left[w+\frac{(u,v)}{\alpha}\cdot
E[n_x\mid s,y]\right]\\
&=\frac{f_s(y)}{2\pi}[w+(u,v)\cdot n_{\psi_s(y)}]
=F_{A_sz}(y).
\end{aligned}
\]

代码还从正向原子逐项换元验证整个输出密度，并直接积分正向测度，检查结果质量、一阶矩及应为零的二、三、六阶余弦矩。这避免只把目标矩阵再算一遍。

完整密度恒等式、真实混合的线性性和旋转恒等式一起，逐分支推出任意有限顺序、自适应、后选择协议保持同一准备编码。

η=0 时，U=V=0，A=(1+α)/2、C=(1−α)/2。因此非选择隐藏更新是留在原角或跳到对径点，其第一谐波保留率仍为 α。原 fresh 的零强度扰动没有被改成恒等操作。该隐藏实现与第二十八轮可以在高阶隐藏过程上不同，但允许准备与协议无法区分它们。

## 8 复算与接续入口

代码：[cyclic_transport_kernel.py](cyclic_transport_kernel.py)。结果：[cyclic_transport_kernel_results.json](cyclic_transport_kernel_results.json)。12 项检查通过，覆盖循环方程、极小连续质量的相对归一化、有理数三角界、全部模式正性、条件向量均值、两种边缘、整个条件密度、直接测度积分和混合一致性。

~~~
python -X utf8 research_cognition_physics/cyclic_transport_kernel.py --write-results
~~~

**下一项任务**：判断 [0.11,η_B] 中首先出现的障碍是什么。可先改进 U_p 的常数偏移与角步长估计，或构造其他保持条件均值的局部转移；同时用第三十轮端点必须满足的单位半径及条件符号限制寻找更强凸函数见证。任何有限扫描的通过都不能充当存在证明，固定编码的失败仍需区分一般编码。

本轮缩小的是固定模型加准备非情境条件后的可行范围。准备非情境性的认知来源仍未解决；第七轮 η=0.5 的一般横向最优保留区间 [0.8570324548,0.8660254038] 也保持为另一未解问题。结果未推出量子理论或引力。
