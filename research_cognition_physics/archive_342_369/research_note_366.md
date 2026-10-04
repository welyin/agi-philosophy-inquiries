# 第366轮：切片不变性与形变代数是否已经选择引力

日期：2026-09-23。完整独立轮次，科学基线冻结至364；接续[351轮](research_note_351.md)与[364轮](research_note_364.md)，不依赖365轮。已读根README、三份研究导航、相关最终稿与结果，并检索旧档确认没有已完成的参数化场论检验。本轮15项检查。

## 1. 问题与结论

364轮说明普通多时间可积性还不是引力的超曲面形变。本轮进一步问：**如果已经升级为真实的连续场、弯曲切片、非恒定空间度规结构函数和完整超曲面形变代数，是否就得到GR？**

答案仍是否定的，范围非常具体：已有的参数化场论可以全部满足这些描述一致性条件，物质有非零能量且正常传播，却没有独立的度规动力学。本轮复算这个已知接口，给出连续解析闭合证明、真正沿两条切片历史积分的代码，以及保持平直背景的非零应力见证。

这关闭的是“重切片不变性／HDA本身足以生成引力”的捷径。它既不是F＋U＋C＋P全部条件的反模型，也不证明认知到有效GR不可能。这里选定的背景、场、相空间与作用量均是额外输入。

| 层次 | 本轮内容 |
|---|---|
| 认知动机 | 同一物理历史能否被不同主体、不同切片一致地描述 |
| 文献接口 | Dirac–Kuchař参数化场论及后续严格构造 |
| 额外模型输入 | 固定Minkowski背景、连续标量场、光滑类空嵌入、正则Poisson括号 |
| 解析结论 | 两个Witt代数、完整HDA、嵌入方程冗余及物理自由度计数 |
| 数值验证 | 非约束曲面上的括号、实际泛函差分、六场Hamilton积分、非零应力和平直曲率 |
| 未推出 | 自然选择该模型、独立度规变量、物质对度规的动力学反作用、Einstein方程 |

## 2. 已有研究接口与归属

- [Kuchař，Phys. Rev. D 39, 1579 (1989)](https://doi.org/10.1103/PhysRevD.39.1579)是圆柱参数化标量场的原始工作之一。本轮核对出版信息；其全文接口未成功读取，不把后续文献的论证冒称为已逐页核读该文。
- 本轮实际读取的原始论文是[Laddha–Varadarajan，*Polymer Parametrised Field Theory*，§2](https://arxiv.org/pdf/0805.0208)与[其Hamiltonian约束论文，§2.1](https://arxiv.org/html/1011.2463)。这里只使用经典作用量、嵌入变量与约束接口，不借用后续聚合量子化结果。
- [Hájíček–Isham，J. Math. Phys. 37, 3505 (1996)，§II—IV](https://arxiv.org/pdf/gr-qc/9510028)在固定四维整体双曲背景和紧Cauchy面上研究嵌入、初值映射及第一类约束的几何结构。这明确说明参数化机制并非仅在1＋1维才存在。

以下重新推导所用符号、因子和括号，并实现独立数值检查；不把参数化场论及其HDA当作本项目原创定理。特别是本轮统一采用正的类空长度平方q以及由两个手征约束直接相减得到的物质系数1／2，不依赖单独展开式的转写。

## 3. 相空间：空间度规来自嵌入，背景度规保持固定

取Minkowski圆柱，空间周长设为2π，背景符号为(−,+)。参数空间坐标为(τ,x)，x属于圆周。每条切片用两个函数描述：

$$
\begin{gathered}
T(x+2\pi)=T(x),\qquad X(x)=x+\xi(x),\qquad
\xi(x+2\pi)=\xi(x),\\
X'>0,\qquad q=X'^2-T'^2>0,\qquad r=\sqrt q,\qquad
n^\mu=\frac{(X',T')}{r},\\
\{T(x),P_T(y)\}=\{\xi(x),P_X(y)\}
=\{\phi(x),\pi(y)\}=\delta(x-y).
\end{gathered}
\tag{1}
$$

其它基本括号为零，素数表示x导数，所有正则变分周期。用ξ而不是带绕数的X做FFT，可以避免把坐标接缝误算成尖峰。

给定参数化标量作用量的正则形式为

$$
\begin{aligned}
S&=\int d\tau\,dx\,[P_T\dot T+P_X\dot X+\pi\dot\phi
                    -N\mathcal H-v\mathcal D],\\
\mathcal D&=P_TT'+P_XX'+\pi\phi',\\
\mathcal H&=\frac{\mathcal A}{r},\qquad
\mathcal A=P_TX'+P_XT'+\frac{\pi^2+\phi'^2}{2}.
\end{aligned}
\tag{2}
$$

N是法向切片推进，v是切向移动。q会随切片改变，但此处q不是带有独立共轭动量的引力场。

定义手征组合可以检查所有归一化：

$$
\begin{aligned}
X^\pm&=T\pm X,&P_\pm&=(P_T\pm P_X)/2,&Y^\pm&=\pi\pm\phi',\\
\mathcal C_\pm&=P_\pm X^{\pm\prime}\pm\frac{(Y^\pm)^2}{4},&
\mathcal D&=\mathcal C_++\mathcal C_-,&
\mathcal H&=\frac{\mathcal C_+-\mathcal C_-}{r},\\
q&=-X^{+\prime}X^{-\prime}>0.
\end{aligned}
\tag{3}
$$

本轮始终保留φ、π这对完整变量，不把Y⁺、Y⁻误当作可逆正则变换而丢掉φ的零模。

## 4. 解析闭合：不是仅在约束曲面上碰巧为零

### 4.1 手征括号

对于任意光滑周期试探函数f、g，记C[f]为积分。手征嵌入变量相互独立，物质流满足

$$
\{Y^\pm(x),Y^\pm(y)\}=\pm2\partial_x\delta(x-y),
\qquad \{Y^+(x),Y^-(y)\}=0,
\qquad
\{\mathcal C_\pm[f],\mathcal C_\pm[g]\}
 =\mathcal C_\pm[fg'-gf'],\quad
\{\mathcal C_+[f],\mathcal C_-[g]\}=0.
\tag{4}
$$

证明可直接用泛函梯度。对符号s＝±1，C_s[f]关于X^s、P_s、φ、π的梯度依次为−(fP_s)'、fX^{s′}、−(fY^s/2)'、sfY^s/2。代入正则括号后，嵌入项给出P_sX^{s′}(fg'−gf')，物质项给出s(Y^s)²(fg'−gf')／4；周期分部积分不留边界项。相反手征的物质交叉项为全导数，嵌入交叉项为零。

### 4.2 法向括号的度规依赖必须参与变分

不能只把f换成N／r，然后把这个依赖相空间的试探函数当成常数。完整H[N]梯度为

$$
\begin{aligned}
\frac{\delta H[N]}{\delta(P_T,P_X,\pi)}
 &=N\left(\frac{X'}r,\frac{T'}r,\frac\pi r\right),\\
\frac{\delta H[N]}{\delta T}
 &=-\partial_x\left[N\left(\frac{P_X}{r}
                         +\frac{\mathcal A T'}{r^3}\right)\right],\\
\frac{\delta H[N]}{\delta \xi}
 &=-\partial_x\left[N\left(\frac{P_T}{r}
                         -\frac{\mathcal A X'}{r^3}\right)\right],\\
\frac{\delta H[N]}{\delta\phi}
 &=-\partial_x(N\phi'/r).
\end{aligned}
\tag{5}
$$

给出一个短的完整推导。记w＝(T',X',φ')、p＝(P_T,P_X,π)、h(w,p)＝A／r。凡密度只依赖这些变量，同一密度的两次涂抹满足

$$
\begin{aligned}
\{H[N],H[M]\}
&=\int dx\,(NM'-MN')\,h_p\cdot h_w,\\
h_p\cdot h_w
&=\frac{X'P_X+T'P_T+\pi\phi'}{q}
  +\frac{\mathcal A(X'T'-T'X')}{r^4}
=\frac{\mathcal D}{q}.
\end{aligned}
\tag{6}
$$

第一行由−(Nh_w)'和Nh_p配对得到；所有含h_w'的项精确相消。第二行显示法向括号中的全部度规变分如何保留并相消，结果不需要D＝H＝0。

混合括号也可直接验证：h对(p,w)联合缩放为一次齐次，故p·h_p＋w·h_w＝h，且h'＝h_p·p'＋h_w·w'。将D[v]梯度−(vp)'、vw代入、分部积分，得到H[vM']。因此完整代数为

$$
\begin{aligned}
\{D[v],D[w]\}&=D[vw'-wv'],\\
\{D[v],H[N]\}&=H[vN'],\\
\{H[N],H[M]\}&=D[q^{-1}(NM'-MN')].
\end{aligned}
\tag{7}
$$

这些恒等式对全部满足式(1)的光滑正则数据和任意光滑周期涂抹函数成立，是经典连续场的强等式。它包含真正的、非恒定的q⁻¹结构函数，并非364轮那种所有多时间方向平坦对易的要求。

## 5. 同一物理解沿不同切片推进

### 5.1 精确解与约束动量

取固定背景上的同一个无质量标量解：

$$
\begin{aligned}
\Psi(T,X)&=f(T+X)+g(T-X),\\
f(u)&=0.3\sin(2u)+0.1\cos(3u),\\
g(v)&=0.2\cos v-0.07\sin(4v),\qquad
(-\partial_T^2+\partial_X^2)\Psi=0.
\end{aligned}
\tag{8}
$$

拉回任意允许切片后，全部六个正则场有明确表达：

$$
\begin{aligned}
\phi&=f(X^+)+g(X^-),\\
\pi&=X^{+\prime}f'(X^+)-X^{-\prime}g'(X^-),\\
P_+&=-X^{+\prime}[f'(X^+)]^2,\qquad
P_-=X^{-\prime}[g'(X^-)]^2,\\
P_T&=P_++P_-,\qquad P_X=P_+-P_-.
\end{aligned}
\tag{9}
$$

由Y⁺＝2X^{+′}f'、Y⁻＝−2X^{−′}g'立刻得两个C_±恒为零。这里动量不是为了逐点拟合输出临时添加；它们就是同一参数化作用量的约束解。

### 5.2 两条真正不同的切片历史

两条路径共享τ＝0与τ＝1端点，期间切片不同：

$$
\begin{aligned}
T(\tau,x)&=0.4\tau+a(\tau)\sin x,\qquad
X(\tau,x)=x+b(\tau)\sin(2x),\\
a_A&=0.15\tau,&b_A&=0.08\tau,\\
a_B&=0.15\tau+0.025\sin(\pi\tau),&
b_B&=0.08\tau-0.02\sin(\pi\tau),\\
N&=\frac{X'\dot T-T'\dot X}{r},&
v&=\frac{-T'\dot T+X'\dot X}{q}.
\end{aligned}
\tag{10}
$$

两条路径的a、b在[0,1]上单调，范围分别为[0,0.15]与[0,0.08]。因此q至少为0.84²−0.15²＝0.6831；Jacobian行列式至少为0.84(0.25−0.025π)−0.15(0.08＋0.02π)＞0。这给出全部连续参数上的类空、未来定向和可逆性保证，数值采样只是额外检查。

代码固定这两组外加N、v，对全部六个正则场积分Hamilton方程，而不是在每一步用式(9)覆盖状态。动量方程包括

$$
\begin{aligned}
\dot P_T&=\partial_x\left[
N\left(P_X/r+\mathcal A T'/r^3\right)+vP_T\right],\\
\dot P_X&=\partial_x\left[
N\left(P_T/r-\mathcal A X'/r^3\right)+vP_X\right],\\
\dot\pi&=\partial_x(N\phi'/r+v\pi),\qquad
\dot\phi=N\pi/r+v\phi'.
\end{aligned}
\tag{11}
$$

在取泛函导数时，N、v被视为已给定的函数；不能再对其目标路径表达式求相空间导数。代码用式(9)的独立时间差分核对整个Hamilton向量场，并用四阶Runge–Kutta独立推进。

这两个端点初值映射相同，有解析理由：它们都来自式(8)同一个背景场解，且类空Cauchy数据确定唯一解。数值四阶收敛是该解析结论在所选数据上的复算。

## 6. 非零物质能量，并没有生成引力反作用

固定背景中的物质应力为

$$
\begin{aligned}
T_{TT}&=[f'(T+X)]^2+[g'(T-X)]^2,\qquad
T_{TX}=[f'(T+X)]^2-[g'(T-X)]^2,\\
E_{\partial_T}
&=\int dx\,[X'T_{TT}+T'T_{TX}]
=-\int dx\,P_T\\
&=\pi(0.6^2+0.3^2+0.2^2+0.28^2)
=1.7856812643004385.
\end{aligned}
\tag{12}
$$

由光锥坐标的单调绕数，两个积分各覆盖一整周期，故所有允许切片得到同一能量。零约束并不意味着物质没有能量；P_T记录的是嵌入变量的正则补偿，不能将它误称为产生背反应的引力能。

参数空间中的度规分量会变化，但它只是固定η的拉回。令J为可逆Jacobian：

$$
g_{ab}[X]=\eta_{\mu\nu}\partial_aX^\mu\partial_bX^\nu,
\qquad
\Gamma^a{}_{bc}=(J^{-1})^a{}_\mu\partial_b\partial_cX^\mu,
\qquad
R^a{}_{bcd}[g[X]]=0.
\tag{13}
$$

最后的零是坐标拉回保持曲率的恒等式。代码另用度规分量差分计算Christoffel符号，与Jacobian公式比较，再差分连接计算Riemann张量。非零连接系数与零曲率同时出现，防止把“分量变化”当成“引力曲率生成”。

## 7. 为什么变化嵌入不等于变化独立度规

这个区别不限于1＋1维。对于任意维数的局部可逆坐标映射，给定固定平直背景和标量作用量：

$$
S[X,\phi]=-\frac12\int d^{d+1}x\,
\sqrt{-g[X]}\,g^{ab}[X]\partial_a\phi\partial_b\phi,
\qquad g[X]=X^*\eta.
\tag{14}
$$

记ξ^a＝(J⁻¹)^a_μδX^μ。嵌入变化只能引起沿坐标重描述的度规变化，不是任意δg：

$$
\begin{aligned}
\delta_Xg_{ab}&=\mathcal L_\xi g_{ab}
=\nabla_a\xi_b+\nabla_b\xi_a,\\
\delta_XS
&=\frac12\int\sqrt{-g}\,T^{ab}\delta_Xg_{ab}
=-\int\sqrt{-g}\,(\nabla_aT^{ab})\xi_b,\\
\nabla_aT^{ab}&=(\Box_g\phi)\partial^b\phi.
\end{aligned}
\tag{15}
$$

边界项在周期空间及固定时间端点下消失。因此X的方程已由物质场方程保证，不产生Einstein方程。它与“把g作为独立场变分”不同；后者的引力作用量和动力学正是尚待推出的内容。

在正则约束保持独立的局部区域，自由度计数也吻合：

$$
\begin{aligned}
1+1:\quad &3\ \text{个正则对}
 -2\ \text{个独立第一类约束}
 =1\ \text{个物理正则对},\\
d+1:\quad &(d+1+1)-(d+1)=1
 \quad\text{个物理标量正则对}.
\end{aligned}
\tag{16}
$$

这里只作局部自由度计数；全局绕数、零模、可允许嵌入与边界条件仍需保留。去掉的嵌入自由度没有留下独立引力度规对。

**本轮不以“1＋1维Einstein引力没有局部传播自由度”作为3＋1维反证。** 非蕴涵的原因是式(15)的变量及变分方向限制；固定四维背景参数化场论已有独立原始构造支撑。即使在3＋1维重写同一背景标量场，坐标拉回仍不会使固定背景获得自主的物质反作用。

## 8. 数值结果与证据强度

128点周期网格上的非约束曲面数据给出：

| 括号 | 数值结果 | 与解析右侧之差 |
|---|---:|---:|
| D—D | 0.08991690171217027 | 1.39×10⁻¹⁷ |
| D—H | 0.039983980277373196 | 2.78×10⁻¹⁷ |
| H—H | 0.08458491244474886 | 1.39×10⁻¹⁷ |
| C₊—C₊ | 0.08881647838193868 | 1.39×10⁻¹⁷ |
| C₋—C₋ | 0.001100423330231587 | 1.02×10⁻¹⁷ |
| C₊—C₋ | −8.52×10⁻¹⁹ | 8.52×10⁻¹⁹ |

非零右侧表明没有靠先设置所有约束为零来掩盖符号或系数错误。另对实际非线性H[N]沿H[M]的Hamilton向量作中心差分，步长0.01、0.005、0.0025的误差依次为7.40×10⁻⁷、1.85×10⁻⁷、4.62×10⁻⁸，显示独立二阶收敛。删除q⁻¹的对照明确失败。

两条完整切片积分与精确终点的六场最大误差：

| 步数 | 路径A | 路径B |
|---:|---:|---:|
| 40 | 1.29×10⁻⁷ | 1.53×10⁻⁷ |
| 80 | 8.05×10⁻⁹ | 9.62×10⁻⁹ |
| 160 | 5.05×10⁻¹⁰ | 6.02×10⁻¹⁰ |

160步两路径端点之差为1.79×10⁻¹⁰，均保持能量到浮点误差。精确拉回数据上的最大约束残差4.42×10⁻¹⁴；采样最小q为0.7056、最小N为0.17146。一个非平凡坐标点的连接最大分量0.4070，连接两种算法之差9.50×10⁻¹³，Riemann最大分量1.96×10⁻¹³。

**验证限制：** 一般量词由第4节连续推导承担。FFT配置有限、乘积及1／r会有高频分量；不能把这些误差报告推广为任意有限网格数据上的精确约束代数。数值也没有量子化约束，更没有证明Virasoro中心项消失。量子表示、正则化与异常属于另外的问题。

## 9. 与351轮及认知主线的关系

351轮从已给定的独立纯度规正则变量、局域动能／曲率ansatz及离壳任意lapse闭合，得到受限Einstein选择。这里相空间是嵌入与物质，q是拉回的复合量，没有独立纯度规动量；两者前提不同，没有矛盾。

本轮把下一道门槛收紧为：**必须说明可操作几何变量为何超出坐标／嵌入描述自由，成为有自主响应的物理场，以及其响应为何落入351或358轮已核验的引力条件。** 只有HDA、同一历史或共同背景上的切片一致性，不能替代这一接口。

下一步可独立检验“只有度规可观测量的关系描述”是否保留足够物理自由度和非平凡应力响应；已有条件定理可复用，但必须先识别哪些变量是规范描述、哪些变量承载实际可测的几何变化。G、Λ的具体数值不构成本轮及最终方程族重建的额外完成门槛。

## 10. 文件与复算

- [代码](366/parametrized_field_constraint_audit.py)
- [结果](366/parametrized_field_constraint_audit_results.json)
- [本轮统一核验入口](366/research_round_366_checks.json)由整合步骤生成。

~~~powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 -m unittest parametrized_field_constraint_audit.Checks -v
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 parametrized_field_constraint_audit.py
~~~

在archive_231_目录执行。15项检查已通过，随后才用既有growing_stream_audit.main的--write-results保存结果。使用Python 3.12.14、NumPy 2.3.5；没有新增依赖、图像或改写旧轮次。
