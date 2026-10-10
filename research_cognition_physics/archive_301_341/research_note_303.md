# 第303轮：标量作用量、完整应力与局部能量账

日期：2026-09-22。承接[302轮](research_note_302.md)的标量传播，在[当前主线](301/spacetime_mainline_priority.md)下建立物质能量接口。本轮仍使用给定背景；几何反作用由[304轮](research_note_304.md)继续。

## 1. 假设与已有工具

**认知动机：** 传播携带哪些资源，变化的背景中怎样交代其能量和压力，避免把未记录的项当作凭空产生或消失？

**待检验假设：** 302轮的方程可由标准协变标量作用量产生；完整变分应力在场方程成立时守恒。该守恒不要求唯一的曲率耦合ξ，也不自动给出背景的引力动力学。

复用Hilbert应力定义与微分同胚不变性的Noether恒等式，参见[Tong，广义相对论§4.5](https://www.damtp.cam.ac.uk/user/tong/gr/grhtml/S4.html)。非最小耦合中的几何应力及不同移项约定可核对[Energy conditions and classical scalar fields，§2.1](https://arxiv.org/abs/hep-th/0106168)。本轮选择保留全部曲率耦合项的Hilbert定义，不把另一种移项后的张量与它混用。

**额外输入：** 光滑Lorentz背景、实标量场、局部二阶作用量、常数m和ξ、标准变分原则。采用c＝1与(−＋＋＋)号差。它们不是仅从有限维量子操作或认知定义已经推得的物质定律。

## 2. 从作用量到传播和应力

令边界变分消失，取：

$$
S_\phi=-\frac12\int d^4x\sqrt{-g}\left[(\nabla\phi)^2+m^2\phi^2+\xi R\phi^2\right],
\qquad \frac{1}{\sqrt{-g}}\frac{\delta S_\phi}{\delta\phi}
=(\Box_g-m^2-\xi R)\phi=:E_\phi.
\tag{1}
$$

对场变分并分部积分就得到302轮方程。对逆度量变分，要同时变化体积、动能收缩和R，得到：

$$
T_{\mu\nu}:=-\frac2{\sqrt{-g}}\frac{\delta S_\phi}{\delta g^{\mu\nu}}
=\nabla_\mu\phi\nabla_\nu\phi
-\frac12g_{\mu\nu}\left[(\nabla\phi)^2+m^2\phi^2\right]
+\xi\left[G_{\mu\nu}\phi^2+(g_{\mu\nu}\Box_g-\nabla_\mu\nabla_\nu)\phi^2\right].
\tag{2}
$$

含ξ的最后一组不能在采用非最小耦合传播时任意丢弃。这里出现Einstein张量G是曲率变分的结果，不等于已经假设或推出Einstein场方程。

将$\delta\phi=X^\mu\nabla_\mu\phi$、$\delta g^{\mu\nu}=-2\nabla^{(\mu}X^{\nu)}$代入作用量变分，利用任意紧支撑X并分部积分，得到不要求场方程成立的恒等式：

$$
\nabla^\mu T_{\mu\nu}=E_\phi\nabla_\nu\phi.
\tag{3}
$$

因此E＝0时应力协变守恒。ξ＝0、1/6及其他常数都满足这个结论；守恒本身不选择ξ。

## 3. 对302轮模式建立可复算能量账

沿用$a=-1/(H\eta)$，记$\mathcal H=a'/a$。取真实场$\phi=\operatorname{Re}[q(\eta)e^{ikx}]$，对一个空间周期平均：

$$
f:=\langle\phi^2\rangle=\frac{\lvert q\rvert^2}{2},\qquad
K:=\langle(\phi')^2\rangle=\frac{\lvert q'\rvert^2}{2},\qquad
\langle(\boldsymbol\nabla\phi)^2\rangle=k^2 f.
\tag{4}
$$

定义共动观察者所见的ρ及空间应力迹的三分之一p。单个有方向的模式通常有各向异性应力；p只是平均压力，不把该模式宣称为各向同性完美流体。

将式(2)平均，空间全导数在周期边界消失，得到：

$$
\begin{aligned}
\rho&=\frac{K+k^2f}{2a^2}+\frac{m^2f}{2}
+\frac{\xi}{a^2}(3\mathcal H^2f+3\mathcal H f'),\\
p&=\frac{K-k^2f/3}{2a^2}-\frac{m^2f}{2}
+\frac{\xi}{a^2}\left[-(2\mathcal H'+\mathcal H^2)f-f''-\mathcal H f'\right].
\end{aligned}
\tag{5}
$$

直接求导，或平均式(3)，得到：

$$
\rho'+3\mathcal H(\rho+p)
=\frac{\operatorname{Re}[Dq\,q'^*]}{2a^2},\qquad
Dq=q''+2\mathcal H q'+[k^2+a^2(m^2+\xi R)]q.
\tag{6}
$$

代码先用不满足运动方程的复多项式验证两边的非零值，再用302轮解析解检查右侧消失。若只保留式(2)的最小耦合部分，却仍按非最小耦合方程传播，残余为：

$$
\rho_0'+3\mathcal H(\rho_0+p_0)=-\frac{\xi R}{2}f'
\quad(Dq=0).
\tag{7}
$$

这给出一个可定位的漏项反例，而非把非零数值误差解释为真实能量不守恒。

## 4. 局部守恒与共动能量的变化

对单位共动体积，物理体积为a³。式(6)在场方程成立时等价于：

$$
E_{\rm com}=a^3\rho,\qquad
\frac{dE_{\rm com}}{d\eta}=-p\frac{d(a^3)}{d\eta},\qquad
\Delta E_{\rm com}+\int p\,d(a^3)=0.
\tag{8}
$$

压力功使共动物质能量随膨胀改变。它不自动定义一个“物质＋引力”的全宇宙守恒总能量；在一般曲时空中，构造全局能量还需对称性和边界条件。本轮测试场的背景仍外部给定，也没有借此完成封闭宇宙的资源证明。

用302轮相同初值，H＝1/2、k＝2、m＝0、η从−4到−1：

| 量 | ξ＝0 | ξ＝1/6 |
|---|---:|---:|
| 初始ρ | 8.0625 | 8 |
| 末态ρ | 0.03531996107503528 | 0.03125 |
| 初始共动能量 | 1.0078125 | 1 |
| 末态共动能量 | 0.28255968860028224 | 0.25 |
| 积分压力功 | 0.7252528113997179 | 0.75 |
| 积分能量账误差 | 1.12×10⁻¹⁶以内 | 机器数值为0 |
| 最大局部守恒残差 | 1.78×10⁻¹⁵以内 | 1.34×10⁻¹⁵以内 |

共形无质量分支还有$-\rho+3p=0$、$a^4\rho=1/2$，分别核验至约1.78×10⁻¹⁵和2.22×10⁻¹⁶。但它的$a^3\rho$从1降至1/4，说明不能把“协变守恒”误写成“任何共动区域的物质能量不变”。若漏掉ξ应力，守恒残差最大为0.0625。

## 5. 独立变分检查

除守恒核验外，代码直接积分带lapse的作用量。固定空间尺度a，令$N=a+\epsilon v$，其中v及v′在区间端点均为0；另做独立场变分。对应的解析基准为：

$$
\left.\frac{dS}{d\epsilon}\right|_{\delta q=v}
=-\frac12\int a^2(Dq)v\,d\eta,\qquad
\left.\frac{dS}{d\epsilon}\right|_{\delta N=v}
=-\int a^3\rho v\,d\eta.
\tag{9}
$$

选择不在解上的实多项式q、m＝0.3和ξ＝0、1/6、1/2，数值变化整个作用量，分别核对波方程与完整ρ。中心差分步长10⁻⁵，最大绝对误差约1.11×10⁻¹¹；两种变分没有用应力公式反过来定义被比较的作用量。

9项检查覆盖两种变分、非解与解上的Noether恒等式、漏项反例、共形迹与红移、积分压力功、实场平均归一化和幅值平方的资源缩放。无图像检验。

## 6. 结论与接续

已把302轮传播接到同一作用量的完整应力及能量账；这为物质源项提供可复算接口。守恒在不同ξ下均成立，不能以“内部资源须自洽”单独选出耦合或引力方程。

若把$\xi G_{\mu\nu}\phi^2$移到几何侧，需同时改变物质／几何划分与后续熵的解释。本轮不把含此项的源不加区分地塞入常系数面积熵假设。下一轮先取ξ＝0的明示比较分支，对接Jacobson的成熟条件性推导，并运行物质和几何共同演化的检验。

[代码](303/scalar_energy_balance_audit.py)、[结果](303/scalar_energy_balance_audit_results.json)、[核验](303/research_round_303_checks.json)。

    python -B -X utf8 research_cognition_physics/archive_231_/scalar_energy_balance_audit.py
