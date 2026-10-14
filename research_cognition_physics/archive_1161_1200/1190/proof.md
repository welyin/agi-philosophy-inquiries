# 1190：同一Cartan自旋零材料的双源、钟频率与完整量子响应

2026-10-14。状态：作者证明稿，待独立审阅。本文是限定材料与作用类别下的一轮机制审计，不是GR生成定理。原文入口与采用范围见[sources](sources.md)。

## 1. 本轮命题、输入及结论范围

研究的问题是：已经给定一个多项式Cartan材料作用后，它能否同时确定动态范数的材料反馈和该材料的钟频率；仅校准钟频率是否足以替代完整来源。

本轮输入如下：

1. 四维光滑有效域、$SO(1,4)$联络$A^{AB}$及内部向量$V^A$；内度规为$\eta_{AB}=\operatorname{diag}(-1,1,1,1,1)$。
2. 正范数分支$V^2=\phi^2>0$，并在局部选$V^A=\phi\delta^A_4$；四标架$e^I=\phi A^{I4}$非退化，$g_{\mu\nu}=\eta_{IJ}e^I_\mu e^J_\nu$。
3. 一个实、自旋零、无其他规范荷的五分量材料$\Phi^A$，及下节明确给定的作用。$m>0,\zeta>0,\mu>0$是常数输入。这个材料不是已从认知公理选择出的标准模型物种。
4. 对钟与量子读口另声明固定度规背景中的单模数学任务；使用$\hbar=1$。其有限区域、背景改变与读口均未认证为原共同量子过程中的完整内部仪器。

主要结论是：辅助分量可精确消去；同一材料的几何变分源和钟参数有确定联系。但钟能隙没有穷尽完整量子来源；允许的多项式质量函数还保留不同材料的相对钟率自由。本轮不选择引力作用或几何动力学，也不证明实际暗部门。目标未完成。

## 2. 原始作用、符号与归一化

令$E^A=DV^A$，$\epsilon_{01234}=+1$，并定义

$$
\Sigma=\frac{\mu}{4!}\epsilon_{ABCDE}V^E E^A\wedge E^B\wedge E^C\wedge E^D,
$$

$$
\Sigma_D=\frac{\mu}{3!}\epsilon_{ABCDE}V^E E^A\wedge E^B\wedge E^C.
$$

这是1209.5358式(22)—(23)的形式记号；这里的$\Sigma$名称已包含原文所写的星号，不额外实施Hodge运算。令

$$
\chi=V_A\Phi^A,\qquad s=V^2,\qquad
U(s,\chi)=m^2 f(s)\chi^2.
$$

本轮选实标量，定义作用为

$$
S_m=-\frac{\zeta}{\mu}\int
\left[
\Sigma_A\wedge
\left(\Phi^A D\Phi^B+D\Phi^B\,\Phi^A\right)V_B
+\frac12\Sigma\,U(s,\chi)
\right].
\tag{1}
$$

两项实标量动能相同，但均保留以对应原式(27)中场及其共轭项。等价地，在本取向下选择原文常数$\kappa_{HC}=-\zeta/\mu$。此定义使约化后的实标量动能具有通常的正能源符号；不是从规范对称性推出该符号。原文展示的Higgs势在这里被替换为正二次势，属于明确的材料选择。

先取$f=1$，随后检验合法多项式族

$$
f_c(s)=1-c+c\,s/v^2,\qquad 0\le c\le1,\quad v>0.
\tag{2}
$$

$v$是固定参考数值，与后文共形换元常数采用相同标记。所有系数均是常数，故(1)—(2)仍为原始场中的多项式。在$s>0$上$f_c(s)>0$，包括$c=1$；没有声称在$s=0$或变号分支上保持本约化。

## 3. 动态范数下的辅助消元：独立推导

定义

$$
q=\Phi^4,\qquad y^I=\Phi^I,\qquad \chi=\phi q,\qquad
X=g^{\mu\nu}\partial_\mu q\,\partial_\nu q.
$$

虽然$E^4=d\phi$可非零，epsilon中已有$V^4$，故

$$
\Sigma=\mu\phi\,\mathrm{vol}_g,\qquad
\Sigma_4=0,\qquad
\Sigma_I\wedge e^J=\mu\phi\,\delta_I^J\,\mathrm{vol}_g.
\tag{3}
$$

由$A^4{}_I=-e_I/\phi$，

$$
D\Phi^4=dq-\frac{e_Iy^I}{\phi},\qquad
V_B D\Phi^B=\phi\,dq-e_Iy^I.
\tag{4}
$$

这里没有设$d\phi=0$。若采用$\chi$作变量，(4)右边应写为$d\chi-\chi\,d\log\phi-e_Iy^I$，不能漏掉中间项。

把(3)—(4)直接代入(1)，记$\partial_Iq=e_I{}^\mu\partial_\mu q$，得到

$$
S_m=-\zeta\int\mathrm{vol}_g\,\phi
\left[
2\phi y^I\partial_Iq-2y_Iy^I
+\frac12m^2 f(\phi^2)\phi^2q^2
\right].
\tag{5}
$$

$y^I$没有导数。其方程是

$$
2\phi\,\partial_Iq-4y_I=0,\qquad
y^I=\frac{\phi}{2}g^{\mu\nu}e^I_\nu\partial_\mu q.
\tag{6}
$$

代数二次型在$\phi>0$时非退化，因而解唯一；这不是把四个独立传播模式截掉。代回(5)得到精确作用

$$
S_m[e,\phi,q]=-\frac12\int\mathrm{vol}_g\,Z(\phi)
\left[X+\mathfrak m^2(\phi)q^2\right],
\qquad
Z=\zeta\phi^3,\quad
\mathfrak m^2=m^2 f(\phi^2).
\tag{7}
$$

在$\chi$表示中同一作用为

$$
S_m=-\frac{\zeta}{2}\int\mathrm{vol}_g\,\phi
\left[
g^{\mu\nu}(\partial_\mu\chi-\chi\partial_\mu\log\phi)
(\partial_\nu\chi-\chi\partial_\nu\log\phi)
+\mathfrak m^2\chi^2
\right].
\tag{8}
$$

(7)—(8)是本轮推导，不能引用原论文的常$V^2$恢复结果代替。这个子类约化后不依赖$\omega^{IJ}=A^{IJ}$，故材料的独立Lorentz自旋源为零；这不等于断言完整Cartan引力的挠率为零。

### 3.1 消元与变分相容

把$y=y_*(e,\phi,dq)$代回原作用，任意变分满足

$$
\delta S_{\mathrm{red}}
=\left.\delta S_m\right|_{y_*}
+\left.\int E_{y^I}\,\delta y_*^I\right|_{y_*}.
$$

第二项因辅助方程恒为零。故完整材料变分可从(7)计算；这个结论不要求$q$在壳。不能对未消去的其他物质Euler项使用同样的删除。

## 4. 固定$A$、固定$e$及固定$h$的径向源

采用清楚的归一：

$$
\delta S_m=\int\mathrm{vol}_g
\left[
\frac12T^{\mu\nu}\delta g_{\mu\nu}
+j_e\,\delta\phi+\mathcal E_q\,\delta q
\right]
+\text{边界项}.
\tag{9}
$$

$j_e$固定$e$与$q$；$\mathcal E_q$是材料Euler导数。定义$T=g_{\mu\nu}T^{\mu\nu}$，则

$$
T_{\mu\nu}
=Z\left[
\partial_\mu q\,\partial_\nu q
-\frac12g_{\mu\nu}(X+\mathfrak m^2q^2)
\right],
\qquad
T=-Z(X+2\mathfrak m^2q^2),
\tag{10}
$$

$$
j_e=-\frac{Z}{2\phi}
\left[
3X+(3\mathfrak m^2+\phi\,\partial_\phi\mathfrak m^2)q^2
\right].
\tag{11}
$$

### 4.1 固定原联络

在同一单位规范中保持$A^{I4}$、$q$不变，径向变分给

$$
\delta e^I=\frac{e^I}{\phi}\delta\phi,\qquad
\delta g_{\mu\nu}=\frac{2g_{\mu\nu}}{\phi}\delta\phi.
$$

故固定$A$的径向源为

$$
j_A=j_e+\frac{T}{\phi}
=-\frac{Z}{2\phi}
\left[
5X+(7\mathfrak m^2+\phi\,\partial_\phi\mathfrak m^2)q^2
\right].
\tag{12}
$$

这是对原$\delta V$源的径向投影，而不是将全部五个$V$方程重新证明。其余规范方向还涉及被运输的物质及联络变分。

为说明联络源没有丢失，若写独立四标架变分为$\int\mathrm{vol}_g\,\tau_I{}^\mu\delta e^I_\mu$，则固定$\phi$时$\delta e^I_\mu=\phi\,\delta A^{I4}_\mu$，独立分量$A^{I4}$的源为$\phi\tau_I{}^\mu$。若对全部反对称指标重复求和，需同时保留其$1/2$约定；本文不混用两种计数。

### 4.2 换到Einstein表示

仅作变量换元，取

$$
h_{\mu\nu}=\frac{\phi}{v}g_{\mu\nu},\quad
\mathrm{vol}_g=\left(\frac v\phi\right)^2\mathrm{vol}_h,\quad
g^{\mu\nu}=\frac{\phi}{v}h^{\mu\nu}.
\tag{13}
$$

此处“Einstein表示”沿用已选引力子类的名字，不表示本文已引入或推导Einstein方程。材料作用成为

$$
S_m=-\frac12\int\mathrm{vol}_h
\left[
\zeta v\phi^2Y+\zeta v^2\phi\,\mathfrak m^2q^2
\right],
\qquad
Y=h^{\mu\nu}\partial_\mu q\,\partial_\nu q.
\tag{14}
$$

固定$h$和$q$时$\delta g_{\mu\nu}=-g_{\mu\nu}\delta\phi/\phi$。以$\mathrm{vol}_h$归一的源因此满足

$$
j_h=\left(\frac v\phi\right)^2
\left(j_e-\frac{T}{2\phi}\right)
=\left(\frac v\phi\right)^2
\left(j_A-\frac{3T}{2\phi}\right),
\tag{15}
$$

$$
j_h=-\zeta v\phi\,Y
-\frac{\zeta v^2}{2}
\left(\mathfrak m^2+\phi\,\partial_\phi\mathfrak m^2\right)q^2.
\tag{16}
$$

直接对(14)求导给出同一式(16)。这是同一作用的具体双源运输，不把一般Ward恒等式当成新结果。它也不单独证明完整引力—物质在壳解或总能守恒。

### 4.3 物质重定义时的Euler项

上述三个源始终固定$q=\Phi^4$。如果改为固定$\chi=\phi q$，则$\delta q=-q\,\delta\phi/\phi$，所以在任一已经固定的几何表示中，

$$
j_{\text{固定}\chi}
=j_{\text{固定}q}-\frac q\phi\,\mathcal E_q.
\tag{17}
$$

对$\phi$依赖的规范化材料变量也有相应项。只有声明材料在壳时才可使$\mathcal E_q=0$；在量子算符、源插入或任意未知状态比较中不能无条件丢弃它。即使作参数依赖的幺正换元，有限时演化也必须运输态、可观测量及该换元的连接项。

## 5. 钟频率与完整量子来源

### 5.1 声明的数学读口

取固定、静态的$h$，在声明的有限体积$V_h$中选一个空间常模，例如周期边界的零模。先把$\phi$当作空间及时间常数。由(14)得到

$$
L_{\mathrm{mode}}=\frac M2(\dot q^2-\omega^2q^2),\qquad
M=\zeta v\phi^2V_h,\qquad
\omega^2=\frac v\phi\,\mathfrak m^2(\phi).
\tag{18}
$$

这里的体积、背景、边界及读口是数学任务输入，不是免费内部仪器。空间非零模另有其动量项；本轮不把零模结论无条件扩成所有钟。

在$g$的常范数片上，静止模频率为$\mathfrak m$；在$h$上为$\omega=\mathfrak m\sqrt{v/\phi}$。二者由$d\tau_h=\sqrt{\phi/v}\,d\tau_g$保持同一相位。这是表示相容，不意味着动态范数下的任意材料都会维持同样频率比。

量子化使用共同Hilbert空间$L^2(\mathbb R,dq)$及固定的$q,p$，$[q,p]=i$：

$$
H(\phi)=\frac{p^2}{2M}+\frac12M\omega^2q^2
=T_{\mathrm{osc}}+W_{\mathrm{osc}}.
\tag{19}
$$

有限正$\phi$区间上，$H$是正的自伴振子族；二次型有共同形式域，Schwartz空间为共同核心。以下导数是该核心上的算符或二次型身份，不声称无界来源在全部Hilbert空间上有统一有限范数。

### 5.2 常质量子类

先取$f=1$。固定$q,p$求导，

$$
\partial_\phi H
=\frac{W_{\mathrm{osc}}-2T_{\mathrm{osc}}}{\phi}.
\tag{20}
$$

定义当前参数的$a_\phi$为

$$
a_\phi=\sqrt{\frac{M\omega}{2}}q+
\frac{i}{\sqrt{2M\omega}}p.
$$

利用$W_{\mathrm{osc}}-T_{\mathrm{osc}}
=\frac{\omega}{2}(a_\phi^2+a_\phi^{\dagger2})$，得到

$$
\boxed{
\partial_\phi H
=-\frac{H}{2\phi}
+\frac{3\omega}{4\phi}
(a_\phi^2+a_\phi^{\dagger2})
}.
\tag{21}
$$

$\partial_\phi H$固定的是原始$q,p$，不是把$a_\phi$错误地当成参数无关变量求导。Hamiltonian反馈源的符号是$-\partial_\phi H$；与(16)积分后的Lagrange源相容。

能级$E_n=\omega(n+\frac12)$仅给

$$
\partial_\phi E_n=-\frac{E_n}{2\phi}
=\langle n;\phi|\partial_\phi H|n;\phi\rangle.
\tag{22}
$$

但同一完整来源还有

$$
\langle 2;\phi|\partial_\phi H|0;\phi\rangle
=\frac{3\sqrt2\,\omega}{4\phi}\ne0.
\tag{23}
$$

故只根据钟能隙领取$-\!H/(2\phi)$会丢掉未知相干态的反馈。身份可张量延拓为$(\partial_\phi H)\otimes I_R$，不改变外部参考$R$；对有定义的能量域内输入成立。这不等于证明此前原量子父过程已经拥有本材料或读口。有限初码的来源压缩也不是全部输出的Fock截断。

## 6. 全Fock有限读数：同谱替代确实改变结果

这一节给式(21)中额外项的一个解析有限见证，不使用有限矩阵截断作为证明。

归一化$\zeta=v=m=V_h=1$，于是

$$
M(\phi)=\phi^2,\qquad \omega(\phi)=\phi^{-1/2},\qquad
H(1)=H_1=\frac12(p^2+q^2).
$$

准备原振子真空$|0;1\rangle$，在$t=0$把背景参数指定为

$$
\phi_*=\frac{121}{100},\qquad
t_*=\frac{11\pi}{20}.
$$

这是明示的数学突变方案，不是已计入能源的内部驱动；不声称瞬时突变是可实施的物理操作。比较的两套Hamiltonian是

$$
H_{\mathrm{Cartan}}=H(\phi_*),\qquad
H_{\mathrm{diag}}=\phi_*^{-1/2}H_1.
\tag{24}
$$

二者有相同的全部振子能谱，并在$\phi=1$一致。两者都在同一$q,p$与原读口下比较；没有把读口随候选更换。

令

$$
k=M(\phi_*)\omega(\phi_*)=\frac{1331}{1000},\qquad
\omega(\phi_*)t_*=\frac{\pi}{2}.
$$

Cartan振子的Heisenberg演化在此时给

$$
q(t_*)=\frac{p}{k},\qquad p(t_*)=-kq.
\tag{25}
$$

相对于原始$a_1=(q+ip)/\sqrt2$，Bogoliubov系数为

$$
\alpha=-\frac i2(k+k^{-1}),\qquad
\beta=-\frac i2(k-k^{-1}),\qquad
n=|\beta|^2
=\left[\frac{1331/1000-1000/1331}{2}\right]^2.
\tag{26}
$$

演化后的真空是归一化的纯挤压真空。解其线性湮灭方程得到相邻偶数系数递推，或展开该态的高斯波函数，均给

$$
P_{2j}
=\frac{(2j)!}{2^{2j}(j!)^2}
\frac{\left(n/(1+n)\right)^j}{\sqrt{1+n}},
\qquad P_{2j+1}=0.
$$

特别地，

$$
\boxed{
P_2^{\mathrm{Cartan}}
=\frac{n}{2(1+n)^{3/2}}> \frac3{100}
}.
\tag{27}
$$

该严格下界可完全用有理数核验：

$$
\frac2{25}<n=
\left(\frac{771561}{2662000}\right)^2<\frac1{10},
\qquad
\left(\frac{11}{10}\right)^3<\frac{16}{9}.
$$

所以$(1+n)^{3/2}<4/3$，代入(27)即得严格$P_2>3/100$。另一方面，$H_{\mathrm{diag}}$只给原真空一个相位，

$$
P_2^{\mathrm{diag}}=0.
\tag{28}
$$

这一差异在完整Fock空间精确成立；有限截断只可作为诊断。它表明把同一材料换成“频率相同的对角振子”会改变原材料读口。它没有证明完整几何反作用、内部参数驱动、实际时空仪器、引力生成或新的物种。

## 7. 合法多项式留下共同钟率自由

对(2)的质量函数，在任意正$\phi$区间材料都具有正动能与正质量平方。常片频率为

$$
\omega_c(\phi)
=m\sqrt{\frac v\phi}\,
\sqrt{1-c+c\frac{\phi^2}{v^2}}.
\tag{29}
$$

在$\phi=v$，所有$c$不仅有同一频率，而且在相同$\zeta,m,V_h$下有同一初始Hamiltonian及$q,p$读口归一：

$$
\omega_c(v)=m,\qquad
\alpha_c:=
\left.\partial_\phi\log\omega_c\right|_{\phi=v}
=\frac{c-\frac12}{v}.
\tag{30}
$$

两个材料$c_1,c_2$的相对钟率响应因此为

$$
\left.\partial_\phi\log\frac{\omega_{c_1}}{\omega_{c_2}}\right|_{\phi=v}
=\frac{c_1-c_2}{v}.
\tag{31}
$$

共同共形因子在比值中消去，故不能靠换时间单位抹掉这项差异。完整来源亦随(11)—(16)中的$\partial_\phi\mathfrak m^2$改变；不能一面更改材料质量函数，一面继续使用$c=0$来源。

结论限于明确的多项式材料类：规范对称性、该作用语法和单点钟校准不选择共同钟率对范数的响应。它不是对完整认知公理或任意可能的共同几何机制的反证。若另有原则限制所有材料的$c$，需要独立提出、审计并检验，不能从“共同钟应当成立”倒写选择结果。

## 8. 本轮新增与停止线

本轮净内容是同一明确材料作用的动态范数消元、双源与表示运输、完整量子来源，以及在相同能谱下不同原读数的全Fock有限见证。没有把一般来源守恒、标准Ward恒等式或任意振子换名重新计为GR生成。

下一步若只有：

- 再选择一个$c$以便得到期望钟率；
- 将$\phi$直接叫作实际时钟或宇宙尺度因子；
- 把已知Einstein表示看成已经生成Einstein动力学；
- 继续解决该候选的全部内部技术问题而不减少原目标的输入，

则应停止这一候选的进一步扩建。本轮保留的净问题是：是否存在独立、可删弱检验的认知材料机制，约束(31)中的自由，同时保留(21)的完整量子来源和全部实际资源。本文没有给出该机制，也没有把它加入现行公理。
