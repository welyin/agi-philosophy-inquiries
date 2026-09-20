# 第三十二轮：保留角度相关性，把严格可行下界推进到 0.1182

日期：2026-09-15。接续第三十一轮。仍固定原窗口半宽 h=1/4、α=sinc(h)，研究第七轮 λ=1 的 fresh 族；不改变原准备和零强度读取规则。

## 1 结果

通过更紧的 Fourier 常数偏移与联合角度估计，构造完整的准备非情境核，实现**全部 η≤591/5000=0.1182** 的读取、全部方向及任意有限自适应组合。

本轮可行性由全角度解析证明和有理数证书保证，不是从采样最小值推断。一般上界仍为第三十轮的 η_B≈0.144739232095，因此

\[
\boxed{0.1182\le\eta_*\le\eta_B}.
\]

另在 η=0.12 证明本轮具体公式出现负原子质量。这限制当前构造，并不排除其他更新核。

## 2 上一轮估计哪里偏松

第三十一轮分别用 R 的最小值和 U 的最大值估计非负性，把两个通常不会出现在同一角度的极端同时代入。数值探索提示：危险点靠近读取轴，此处 R 也较大；保留这一相关性能够改善证书。

另一个偏松处是：把全部 Fourier 模式都按最小差分特征值估计。实际从第二个模式起，主要系数对应较大的特征值，可以单独处理。

最初扫描不同循环阶数只用于选择候选参数，不用于证明。最终固定 N=20，d=π/20，c=cos d；N 是函数方程中的循环项数，隐藏角仍连续。

## 3 同一反向耦合，使用更小的偏移

沿用第三十一轮的 R、V、V_0 和循环 Green 函数：

\[
\begin{aligned}
R(\psi)&=1+2\sum_{k\ge1}c_k\cos(2k\psi),
&c_k&=r^k(1+2kq),\\
q&=\sqrt{1-\eta^2},
&r&=\frac{1-q}{1+q}
=\left(\frac{\eta}{1+q}\right)^2,\\
U_p(\psi)&=\sum_{j=0}^{N-1}G_j[R(\psi+jd)-1],
&G_j&=\frac{N^2-1}{12N}-\frac{j(N-j)}{2N}.
\end{aligned}
\]

对于 k 不是 N 倍数的模式，U_p 的余弦振幅为

\[
a_k=\frac{c_k}{2\sin^2(kd)}>0.
\]

所有非零振幅为正，所以它们的总和正好等于 U_p(0)。取

\[
\boxed{S=U_p(0),\qquad U(\psi)=S+U_p(\psi)\ge0.}
\]

这里的 S 通过有限 Green 求和精确表示，不需要截断无限 Fourier 级数。S 比第三十一轮使用的统一振幅上界小；它仍不一定是保证 U≥0 的最小常数，后者不在本轮宣称的最优范围内。

继续令

\[
V=V_0+H-1,\quad
C=\frac{(1-\alpha)R-V}{2}-(1-c)U,\quad
A=\frac{(1+\alpha)R-V}{2}-(1+c)U.
\]

H 是 R 的 N 项循环平均，V_0 和 H−1 使用第三十一轮的闭式。由于常数偏移不改变 LU=R−H，边缘均匀性和条件向量均值的证明保持成立：

\[
A+2U+C+V=R,\qquad
A+2cU-C=\alpha R,\qquad
\int R(\psi)P(dx\mid\psi)\frac{d\psi}{2\pi}=\frac{dx}{2\pi}.
\]

## 4 联合角度估计

记

\[
E_2=2\sum_{k\ge2}c_k,\quad E_3=2\sum_{k\ge3}c_k,\quad
B=\frac{c_2}{2\sin^2(2d)}+\frac{E_3}{4\sin^2d}.
\]

忽略已被循环平均移除的谐波只会增大这个上界，故

\[
|U_p-a_1\cos2\psi|\le B,\quad
S\le a_1+B,\quad
U\le a_1(1+\cos2\psi)+2B,
\]

同时 R≥1+2c_1cos2ψ−E_2，V≤2V_0。

令 δ=1−α、H_α=1+α。代入 2C 与 2A 后，cos2ψ 的系数分别为

\[
\left(2\delta-\frac1{1+c}\right)c_1,\qquad
\left(2H_\alpha-\frac1{1-c}\right)c_1.
\]

对 N=20 和本项目 α，两者均为负，因此这些线性下界在 cos2ψ=1 时最小。利用 sin²d=(1−c)(1+c)、sin²(2d)=4c²sin²d，得到全角度充分界

\[
\begin{aligned}
2C\ge\;&
\delta-\left(\frac2{1+c}-2\delta\right)c_1
-\delta E_2-\frac{c_2}{2(1+c)c^2}
-\frac{E_3}{1+c}-2V_0,\\
2A\ge\;&
H_\alpha-\left(\frac2{1-c}-2H_\alpha\right)c_1
-H_\alpha E_2-\frac{c_2}{2(1-c)c^2}
-\frac{E_3}{1-c}-2V_0.
\end{aligned}
\]

这一步没有假定实际 C 或 A 的最小值一定在读取轴上；只需所构造的线性下界在那里最小。

## 5 有理数证书覆盖全部较低强度

令 e=591/5000。取

\[
q_L=1-\frac{e^2}{2}-\frac{e^4}{8}-\frac{e^6}{8},\qquad
r_U=\frac{1-q_L}{1+q_L}.
\]

用精确分数直接检查 q_L>0、q_L²<1−e²，便有 q≥q_L、r≤r_U 对全部 η≤e 成立。无需对 η 作有限扫描。

第一系数满足

\[
c_1(q)=\frac{(1-q)(1+2q)}{1+q},\qquad
\frac{dc_1}{dq}=-\frac{2q(q+2)}{(1+q)^2}<0,
\]

所以 c_1≤r_U(1+2q_L)。其余项用 q≤1 作统一上界：

\[
\begin{aligned}
c_2&\le5r_U^2,\\
E_2&\le2r_U^2\left[\frac5{1-r_U}+\frac{2r_U}{(1-r_U)^2}\right],\\
E_3&\le2r_U^3\left[\frac7{1-r_U}+\frac{2r_U}{(1-r_U)^2}\right],\\
V_0&\le2\left[\frac{r_U^N}{1-r_U^N}
+\frac{2Nr_U^N}{(1-r_U^N)^2}\right].
\end{aligned}
\]

α 的上下界沿用第二十七轮。c 的有理数区间使用

\[
\cos(\pi/20)=
\sqrt{\frac{1+\frac14\sqrt{10+2\sqrt5}}2},
\]

对每层平方根作精确整数平方比较，得到 30 位小数精度的包围区间。整个证书只有有理数运算与整数平方根。

将这些界按单调方向代入第四节，有

\[
\boxed{2C>0.0000164782656,\qquad 2A>0.3223161594.}
\]

完整分子、分母保存在 JSON 中；这里的十进制只是展示。U、V 的非负性已由正 Fourier 系数保证。因此全部 η≤0.1182 都存在明确正核。

在 η=0.1182，连续补偿的总质量 V_0≈6.8443051×10^−48。代码仍计算并验证它的相对归一化，没有把它删除。

## 6 正向实现与严格失败点

对应的正向仪器仍为四个偏移 j∈{0,d,−d,π} 的原子加连续项。令 W_j={A,U,U,C}，ψ=x−j，则

\[
p_{s,j}(x)=\frac{q^3W_j(\psi)}
{2(1-s\eta\cos\psi)^2R(\psi)},\qquad
y=T_s(\psi),
\]

\[
k_s^{\rm reset}(y)=
\frac{f_s(y)V(\psi_s(y))}{2\pi R(\psi_s(y))}.
\]

任意方向先转入读取轴坐标，输出后转回。第二十九轮的联合测度论证或逐原子换元均给出

\[
\boxed{K_sF_z=F_{A_sz}}
\]

作为整个密度恒等式，进而包含完整记录和任意有限自适应协议。

另一方面，η=3/25=0.12 时，U_p(0)≥a_1，故

\[
2C(0)\le
(1-\alpha)\frac{1+\eta^2}{q}
-\frac{2c_1}{1+c}.
\]

使用 q、α、c 的有理数区间向上估计，严格得到

\[
2C(0)<-0.0002827505<0.
\]

这在正向核中对应一个与其他原子位置不同的负原子，不能由连续正质量抵消。因此失效是真正的概率测度失效，而非充分估计失效。它只否定这一个公式，不能把 0.12 当作一般上界。

## 7 复算与后续

代码：[refined_cyclic_transport.py](refined_cyclic_transport.py)。结果：[refined_cyclic_transport_results.json](refined_cyclic_transport_results.json)。10 项检查通过，包括精确平方根区间、全局有理数余量、Fourier 偏移、完整双边缘与条件均值、整个选择性输出密度、独立测度积分、高阶矩、真实混合和负原子证书。

~~~
python -X utf8 research_cognition_physics/refined_cyclic_transport.py --write-results
~~~

原第三十一轮代码与结果保持其当时的 0.11 证书，本轮单独保存改进实现。下一轮转向上界见证的表达能力及端点条件，避免只在同一个充分估计里反复微调。准备非情境性依然是新增条件，原认知原则没有因此导出量子理论。
