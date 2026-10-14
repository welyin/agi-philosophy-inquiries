# 第1186轮证明：复合背景、完整局部接触与有限闭合量子对照

## 0. 命题范围与旧成果复用

本文件证明两件限定清楚的事：

1. 对指定复合标量候选的局部母作用，同一背景的一阶驻定固定完整零导数 Hessian，并固定它与树交换张量的代数关系；这一接触不是独立可调的局部系数。
2. 对同类四维共形标量的一个实际正量子模式，保留或删除该模式的二阶接触，会改变一个有限时间、有限来源幅度下的内部联合量子过程及末端记录；两者可以拥有相同的基点 Hamiltonian 和一阶来源。

它们没有证明动态几何生成、正的无质量公共张量极点、完整广义相对论，或该文献的全部调节与复合算符已有实际量子实现。

直接复用625的共同二阶来源接口、626及704对固定准备与随来源重准备的区分、632的旧费米背景与接触匹配、739—740的因果反馈范围。本轮不重新证明通用 CTP 或 Ward 恒等式。净增量是新复合标量母作用的背景—接触—树张量身份，以及一个完整有限闭合量子对照。第一部分的形式背景期望不被宣布为第二部分振子的 Lorentz 真空。

主要文献是 Batz、Erlich、Mrini，Background-Independent Composite Gravity，Class. Quantum Grav. 38, 095008 (2021)：
https://arxiv.org/html/2011.06541

其中式（39）已经保留相关局部二次项。本轮补的是它与完整背景、来源及实际有限响应的共同核验，不指称文献最终传播错误。

## 1. 固定测度下的同一局部母作用

令背景度规为非退化的 $G_{\mu\nu}$，$D>2$。考虑

$$
S_{\rm loc}[g,\phi]
=\int d^Dx\,\sqrt{|g|}
\left[\frac12g^{\mu\nu}A_{\mu\nu}-P_{\rm pot}\right],
$$

其中

$$
A_{\mu\nu}=\sum_a\partial_\mu\phi^a\partial_\nu\phi^a,
\qquad
P_{\rm pot}=V_0'+\frac{m^2}{2}\sum_a(\phi^a)^2.
$$

本节的来源微分保持材料变量、参考态及已选定的测度不变。复合算符的减除、测度 Jacobian、额外局部作用及来源依赖准备若存在，其导数须另外加入；不能通过本节的固定测度身份将它们删去。

假定在所选背景，已经给定同一处方下的局部期望

$$
\langle A_{\mu\nu}\rangle=aG_{\mu\nu},
\qquad
\langle T^{\rm total}_{\mu\nu}\rangle=0.
$$

这里的 $a$ 是同一背景的一次来源数据，不是引入的新动力学场。由

$$
T^{\rm total}_{\mu\nu}
=A_{\mu\nu}
-g_{\mu\nu}
\left(\frac12g^{\alpha\beta}A_{\alpha\beta}-P_{\rm pot}\right)
$$

立即得到

$$
\langle P_{\rm pot}\rangle
=\left(\frac D2-1\right)a.
\tag{1}
$$

这是一阶总驻定条件。它不意味着二阶来源等于零，也不固定来源改变以后的新状态。

这些局部期望可以是指定重整化处方的量。它们不是已经得到的正有限 Gram 期望构造；复合算符的签名与正量子实现不在此假定中偷渡完成。

## 2. 完整零导数二阶张量

令 $g=G+h$，并将一个指标用固定背景升起：

$$
H^\mu{}_\nu=G^{\mu\alpha}h_{\alpha\nu}.
$$

使用

$$
\sqrt{|\det(I+H)|}
=1+\frac12\operatorname{tr}H
+\frac18(\operatorname{tr}H)^2
-\frac14\operatorname{tr}H^2+O(H^3),
$$

以及

$$
\left\langle g^{\mu\nu}A_{\mu\nu}\right\rangle
=a\,\operatorname{tr}(I+H)^{-1}
=a\left[D-\operatorname{tr}H+\operatorname{tr}H^2\right]
+O(H^3),
$$

将式（1）代入同一作用，得到

$$
\langle S_{\rm loc}[G+h]\rangle
=\int d^Dx\,\sqrt{|G|}
\left[
a+\frac a4
\left(\operatorname{tr}H^2-\frac12(\operatorname{tr}H)^2\right)
\right]+O(h^3).
\tag{2}
$$

线性项消失，二次项一般不消失。

注意式（2）中的 $a/4$ 是作用展开的二次系数。若定义 Hessian $C$ 使

$$
S^{(2)}_{\rm loc}=\frac12\int\sqrt{|G|}\,h:C:h,
$$

则

$$
C=\frac a2Q,
\qquad
Q=I-\frac12J,
\qquad
J(h)=G\,\operatorname{tr}(G^{-1}h).
\tag{3}
$$

因此对两个变化 $h,k$，完整局部 Hessian 的双线性形式为

$$
\delta_h\delta_k\langle S_{\rm loc}\rangle
=\frac a2\int\sqrt{|G|}
\left[
\operatorname{tr}(G^{-1}hG^{-1}k)
-\frac12\operatorname{tr}(G^{-1}h)\operatorname{tr}(G^{-1}k)
\right].
\tag{4}
$$

局部 delta 函数与测度由这里的积分约定确定。式（4）只覆盖无导数部分，不包含另外给定或诱导的曲率项。

## 3. 与原文树张量的代数逆一致

原文式（31）的背景条件，在其声明的低导数匹配范围内给出

$$
a=\frac{2(V_0+\Lambda)}{D-2}.
\tag{5}
$$

定义原文式（38）的张量算子

$$
P=\frac{D-2}{2}I-\frac12J.
$$

由于 $J^2=DJ$，

$$
Q^{-1}=I-\frac{J}{D-2}
=\frac{2P}{D-2}.
$$

在 $a\ne0$、$D>2$ 时，式（3）—（5）给

$$
C^{-1}=\frac{2P}{V_0+\Lambda}.
\tag{6}
$$

在同一背景升降指标约定下，原文式（44）的树振幅为

$$
A_0=-\frac{P}{2(V_0+\Lambda)}
=-\frac14C^{-1}.
\tag{7}
$$

这个四分之一因子与两端 $h:T/2$ 耦合相符。等价地，从

$$
\mathcal L^{(2)}
=\frac12h:C:h-\frac12h:T
$$

代数消去 $h$，得到

$$
\mathcal L_{\rm int}
=-\frac18T:C^{-1}:T
=-\frac{1}{4(V_0+\Lambda)}T:P:T,
\tag{8}
$$

对应原文式（40）。式（7）的振幅与式（8）的作用系数不能混用。

**结论：** 在这一指定母作用和背景匹配下，零导数树接触、背景一次期望和树反馈系数不是三份独立输入。这不固定二导数 Einstein 系数，也不把代数逆自动解释为正量子传播子。 此处 $C$ 是所选局部母作用的单顶角接触 Hessian，不是加上连通涨落之后的总有效作用二阶核。后者的零动量常数可以发生抵消；不能由 $C\ne0$ 推断引力子质量非零。

## 4. 共形方向、正规序及保留的接触

在四维取光滑有限窗 $f(x)$，并令

$$
g_{\mu\nu}(s,x)=(1+sf(x))G_{\mu\nu}(x),
\qquad |s|\,\|f\|_\infty<1.
$$

因为最小标量作用不含度规导数，固定材料变量时

$$
S_{\rm free}(s)
=\int\sqrt{|G|}
\left[
\frac{1+sf}{2}A
-\frac{(1+sf)^2m^2}{2}B
\right],
$$

其中 $A=G^{\mu\nu}A_{\mu\nu}$，$B=\sum_a(\phi^a)^2$。故

$$
S_{\rm free}''(0)
=-m^2\int\sqrt{|G|}f^2B.
\tag{9}
$$

加入同一个常数势 $V_0'=V_0+\Delta V_{\rm ct}+\Lambda$ 后，

$$
S_{\rm total}''(0)
=-\int\sqrt{|G|}f^2(m^2B+2V_0').
\tag{10}
$$

原文式（35）的匹配为
$\Delta V_{\rm ct}=-m^2\langle B\rangle/2$。它可以使式（9）的期望与该反项的接触相消，但同一总作用仍留下

$$
\langle S_{\rm total}''(0)\rangle
=-2(V_0+\Lambda)\int\sqrt{|G|}f^2
=-2a\int\sqrt{|G|}f^2.
\tag{11}
$$

这正是式（4）在 $h=fG$ 下的结果。不能一方面沿用背景匹配，另一方面把全部二阶接触删除。

如果采用来源依赖的正规序或测度，必须对其减除项及 Jacobian 一起求导。式（9）单独非零不证明原文最终振幅错误；式（11）给出的是总作用中必须共同保留的局部项。这里也没有用一次驻定取代二阶因果响应。

## 5. 同类四维共形标量的一个正量子模式

取 $\hbar=m=1$、空间体积归一的零动量模式，令四维共形参数 $y>0$。相应 Lagrangian 为

$$
L_y=\frac y2\dot q^2-\frac{y^2}{2}q^2+\frac{y^2}{8}.
$$

作同一 Legendre 变换，$p=y\dot q$，得到实际 Hamiltonian

$$
H_y=\frac{p^2}{2y}+\frac{y^2q^2}{2}-\frac{y^2}{8},
\qquad [q,p]=i.
\tag{12}
$$

这保留完整 $L^2(\mathbb R)$ 振子，不截断成有限能级。它是一个给定共形来源的有限量子模型，不是完整协变调节器或几何场动力学。 本模式的来源导数固定 canonical 变量 $(q,p)$；前文协变作用的导数固定材料字段 $\phi$。两者由明确 Legendre 变换联系，但其二阶局部项不逐项相等，不能直接交换使用。

固定 $y=1$ 的基态 $|0\rangle$，其能量为 $E_1=3/8$。基点来源及二阶接触算符为

$$
O=H'_1=-\frac{p^2}{2}+q^2-\frac14,
\qquad
K_{\rm op}=H''_1=p^2+q^2-\frac14.
\tag{13}
$$

直接计算：

$$
\langle O\rangle=0,\qquad
K=\langle K_{\rm op}\rangle=\frac34.
\tag{14}
$$

常数项已经使一次来源驻定，但二阶接触非零。来源依赖的常数项是真实相位的一部分；在后面的量子控制过程里，它不能被当成所有分支共有的不可观测相位。

## 6. 固定初态的因果响应

以 $a_{\rm osc}$ 表示基点振子的降算符。式（13）等价于

$$
O=\frac34(a_{\rm osc}^2+a_{\rm osc}^{\dagger2})
+\frac12a_{\rm osc}^\dagger a_{\rm osc}.
$$

因此

$$
|\langle2|O|0\rangle|^2=\frac98,
$$

而跃迁频率为 $2$。对正号 Hamiltonian 扰动 $+\delta y\,O$，实际来源 $H'_y$ 的固定准备线性响应为

$$
\chi_{\rm ret}(t)
=\frac34\delta(t)-\frac94\theta(t)\sin(2t).
\tag{15}
$$

第一项来自来源算符本身的变化；第二项来自同一过程的对易子。这里没有重新准备一个随 $y$ 改变的基态。

静态绝热核验可用收敛因子积分：

$$
\chi_{\rm ret}(0)
=\frac34-\frac94\frac12=-\frac38.
$$

另一方面，式（12）的精确基态能量是

$$
E_0(y)=\frac{\sqrt y}{2}-\frac{y^2}{8},
\qquad E_0''(1)=-\frac38.
\tag{16}
$$

两式一致。这个检验不把一般静态重准备等同于有限时间动力过程；它只是本有隙模式的独立校核。

## 7. 只改变二阶项的对照动力学

在同一参数区间 $1\le y\le3/2$ 定义

$$
\widetilde H_y
=H_y-\frac{(y-1)^2}{2}K_{\rm op}.
\tag{17}
$$

两个族有相同的 $H_1$、$H'_1$，三阶及更高阶来源导数也相同；只改变二阶项：

$$
\widetilde H''_1=0.
$$

这不是只修改读取定义。后续整个联合过程使用式（17）的真实演化。

写成

$$
H=\frac a2p^2+\frac b2q^2+c,
$$

完整族的系数是

$$
a=y^{-1},\quad b=y^2,\quad c=-y^2/8;
$$

删除族的系数是

$$
\widetilde a=y^{-1}-(y-1)^2,\quad
\widetilde b=2y-1,\quad
\widetilde c=-(2y-1)/8.
\tag{18}
$$

在整个区间，

$$
a\ge\frac23,\quad b\ge1,\quad c\ge-\frac9{32},
$$

$$
\widetilde a\ge\frac5{12},\quad
\widetilde b\ge1,\quad
\widetilde c\ge-\frac14.
$$

两族在 Schwartz 核上本质自伴，共用谐振子算符域；谱下界为 $\sqrt{ab}/2+c$。上述界给出

$$
\inf\sigma(H_y)
\ge\frac12\sqrt{\frac23}-\frac9{32}
>\frac{19}{160},
$$

$$
\inf\sigma(\widetilde H_y)
\ge\frac12\sqrt{\frac5{12}}-\frac14
>\frac1{16}.
\tag{19}
$$

因此两族均有共同严格正下界，特别是都大于 $3/160$。这里的负常数不造成无下界 Hamiltonian。证明不依赖有限矩阵截断。

## 8. 有限闭合量子控制与未知资料保全

加入量子控制比特 $B$，$P=|1\rangle\langle1|$。固定联合规律为

$$
H_{\rm full}
=(I-P)\otimes H_1+P\otimes H_{3/2},
$$

$$
H_{\rm delete}
=(I-P)\otimes H_1+P\otimes\widetilde H_{3/2}.
\tag{20}
$$

两个联合 Hamiltonian 都自伴且下有界。$P$ 守恒，但控制比特的相干会受到与振子相互作用的反作用。

对任意未知的 $\rho_{BR}$，演化是 $I_R\otimes U_{Bq}$ 的正当酉过程；不复制、预先读取或经典化未知控制资料。后续的具体概率见证选择 $|+\rangle_B|0\rangle_q$，不限制式（20）对其它输入的定义。

分支常数、控制系统、振子及末端效应均在同一数学联合过程中交代。固定选择的系数、总 Hamiltonian 和读取效应仍是模型输入；这里没有认证任意物理控制硬件，也没有宣称 $y$ 已成为动态时空几何。

## 9. 有限时间实际效应

对一个分支 $H=ap^2/2+bq^2/2+c$，令

$$
\Omega=\sqrt{ab},\qquad r=\frac{a+b}{2\Omega}.
$$

基点真空的精确重叠是

$$
L(t)
=e^{iE_1t}\langle0|e^{-iHt}|0\rangle
=
\frac{e^{i(3/8-c)t}}
{\sqrt{\cos(\Omega t)+ir\sin(\Omega t)}}.
\tag{21}
$$

平方根沿 $t=0$ 的值 $1$ 连续延拓。本例 $t=1$ 时两个分支的余弦均为正，主平方根与该连续选择一致。

初态为 $|+\rangle_B|0\rangle_q$。末端只读取控制比特效应

$$
E_{Y+}=\frac{I+Y}{2}\otimes I_q.
$$

由控制比特约化态得到

$$
p(Y+)=\frac{1+\operatorname{Im}L(t)}2.
\tag{22}
$$

在 $y=3/2,t=1$：

|过程|$a$|$b$|$c$|$p(Y+)$|
|---|---:|---:|---:|---:|
|完整|$2/3$|$9/4$|$-9/32$|$0.5082573380439783$|
|删除二阶接触|$5/12$|$2$|$-1/4$|$0.5470527495625932$|

概率差约为

$$
p_{\rm delete}(Y+)-p_{\rm full}(Y+)
=0.0387954115186149.
\tag{23}
$$

这些小数是式（21）—（22）的近似值；严格有限精度下界由本轮独立区间核验给出，不能以浮点小数本身代替认证。

两套过程共享基点和一阶来源，却在有限时间、有限幅度下产生不同实际效应概率。记录差异来自不同的完整酉过程，不是给同一过程换一个来源名称。

## 10. 判定与停止线

已经成立：

- 指定复合母作用中，一次背景匹配固定完整零导数 Hessian，并固定 $A_0=-C^{-1}/4$。
- 正规序不能授权删除同一总作用的剩余接触。
- 同类正量子模式的二阶接触参与实际因果响应，删除它产生可区分的有限闭合联合记录。

尚未成立：

- 该文献复合度规的正有限量子定义、Lorentz 签名与实际参考共同实现；
- 无预置二导数项条件下的正无质量公共张量极点；
- 全来源普适匹配、完整非线性引力与暗部门的生成。

后续应检查新候选能否在同一量子定义中保持背景、完整来源及响应；不得把本轮有限振子当作该候选的全部协变实现，也不应重复优化这个概率见证以拖延上述机制问题。
