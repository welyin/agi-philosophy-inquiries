# 1284 证明：Gauss初值的局部接续如何约束同一物质顶点

2026-10-18。[报告](../research_note_1284.md)、[复算](check.py)、[结果](results.json)、[独审](independent_review.md)、[后继](NEXT.md)。本轮不采用新角色；验收以正式独审及主线记录为准。

## 1. 新桥、旧结果与明示输入

[1016](../../archive_1009_1043/research_note_1016.md)已分类指定复标量顶点：对任意平滑场及紧支参数的局部作用不变，等价于其五项系数关系。本轮不重报这个分类，而检查一个较接近接续任务的前件：**同一原拉格朗日量产生的全部主约束—Gauss次约束，对任意允许局部初值和任意时间连接描述保持，不另收缩物质初值，不增加约束或固定乘子。**

继承1016的四维平直有效域、正则复标量 $\phi\in\mathbb C^n$、健康已给Yang–Mills矢量作用、实紧群Lie代数及归一 $g$。固定Hermitian物质生成元 $R_a$，本轮已输入

$$
[R_a,R_b]=ig f^c{}_{ab}R_c.\tag{1}
$$

它使下面的原Gauss生成元形成合法同一表示，**不再统计为新选出的关系**。群、表示、Lorentz主部、原生操作字典及本节有效类别不由本轮生成。

开放的物质拉格朗日量仍是1016的类别，$C_a,S_a,K_{ab}$ 为常Hermitian矩阵，$K_{ab}=K_{ba}$，$V$ 为光滑实无导数势：

$$
\begin{aligned}
\mathcal L_m={}&-\partial_\mu\phi^\dagger\partial^\mu\phi-V
-iA^{a\mu}[\phi^\dagger C_a\partial_\mu\phi-(\partial_\mu\phi)^\dagger C_a\phi]\\
&-A^a_\mu A^{b\mu}\phi^\dagger K_{ab}\phi
+(\partial_\mu A^{a\mu})\phi^\dagger S_a\phi,
\qquad \eta=(-,+,+,+).\tag{2}
\end{aligned}
$$

若原式同时有 $A\cdot\partial(\phi^\dagger T\phi)$，这里的 $S$ 是 $S_{\rm eff}=S_{\rm raw}-T$；纯全导数不被排除。固定合法边界变分问题，全部改写连动量、来源和端点一起运输。

## 2. 完整Legendre字典，不能用主约束直接删掉S

原式含 $-\dot A_0^a\phi^\dagger S_a\phi$，因此原主约束是 $\Pi^0_a+\phi^\dagger S_a\phi=0$。现在只作同一作用的合法分部积分，取边界生成函数

$$
B=\int d^3x\,A_0^a\phi^\dagger S_a\phi,
\qquad \mathcal L'_m=\mathcal L_m+\partial_t(A_0^a\phi^\dagger S_a\phi)
-\partial_i(A_i^a\phi^\dagger S_a\phi).\tag{3}
$$

令 $C_0=A_0^aC_a,S_0=A_0^aS_a,K_{00}=A_0^aA_0^bK_{ab}$。时间部分变为

$$
\mathcal L'_{\rm time}=\dot\phi^\dagger\dot\phi
+\dot\phi^\dagger(S_0-iC_0)\phi
+\phi^\dagger(S_0+iC_0)\dot\phi+\phi^\dagger K_{00}\phi.
$$

以行向量 $\pi$ 为 $\phi$ 的共轭动量，得到

$$
\pi=\dot\phi^\dagger+\phi^\dagger(S_0+iC_0),\qquad
\Pi^0_a=0,\qquad
\dot\phi=\pi^\dagger-(S_0-iC_0)\phi.\tag{4}
$$

物质速度Hessian可逆，时间连接保持原退化方向。新旧辛势之差正是 $dB$；$S$ 已移进物质动量和次约束，未被丢掉。

记 $Q_{C,a}=i(\pi C_a\phi-\phi^\dagger C_a\pi^\dagger)$、$K_{S,a}=\pi S_a\phi+\phi^\dagger S_a\pi^\dagger$，以及

$$
M(A_0)=(S_0+iC_0)(S_0-iC_0)-K_{00}.
$$

包含原矢量部分的正则Hamiltonian为

$$
H_c=H_0(C,K,S,V)+\int d^3x\,
\left[A_0^a(-D_iE^i_a+Q_{C,a}-K_{S,a})+\phi^\dagger M(A_0)\phi\right].\tag{5}
$$

$H_0$ 含全部 $\pi\pi^\dagger$、矢量电磁能及

$$
\mathcal H_{m,\rm sp}=|\partial_i\phi|^2+V
+iA_i^a[\phi^\dagger C_a\partial_i\phi-(\partial_i\phi)^\dagger C_a\phi]
+A_i^aA_i^b\phi^\dagger K_{ab}\phi+A_i^a\partial_i(\phi^\dagger S_a\phi).
$$

使用由原L实际得到的 **total Hamiltonian** $H_T=H_c+\int u^a\Pi^0_a$，不任意把所有次约束再配上新独立乘子。原预定Gauss约束及其空间生成作用为

$$
G_{R,a}=-D_iE^i_a+Q_{R,a}=0,
\quad Q_{R,a}=i(\pi R_a\phi-\phi^\dagger R_a\pi^\dagger),
\quad \delta_\epsilon\phi=i\epsilon^aR_a\phi,\quad
\delta_\epsilon A_i^a=(D_i\epsilon)^a.\tag{6}
$$

这是在式(3)共同正则代表中指定的纯内部旋转电荷，不是允许事后把任意 $S,C$ 都吸收进一个重新命名的Gauss电荷。

## 3. Gauss壳的局部投影：物理初值不限制这份系数缺陷

在开空间小盒内，固定任意光滑局部 $\phi,\pi,A_i$，任取 $E^2,E^3$ 与横截面上的有限 $E^1$ 初值，沿 $x^1$ 解

$$
D_1E^1=Q_R-D_2E^2-D_3E^3.\tag{7}
$$

这是一阶线性常微分方程，光滑系数给局部光滑解。因而Gauss壳投影到任意这些物质／连接有限Cauchy jet是满射。该事实直接接[1027 §4](../../archive_1009_1043/research_note_1027.md)的相容场资料方法；不需要把任意非物理量子态当作可准备态。

在固定紧子盒上，平滑有界初值、系数及有限端点数据使 $E$ 与全部保留阶局部源有限。例如紧群正内积下平行运输正交，沿长度 $L$ 有 $|E^1(x)|\le |E^1(0)|+\int_0^x|Q_R-D_2E^2-D_3E^3|ds$。这给该盒内有限电场能源；其它场和动量亦光滑有界。不声称任意资料能以同一固定全局边界荷完成：若预先把所有边界通量和全域物质都固定，式(7)可能有整体相容条件。

本轮的“任意物理初值”指这个**经典有效局部约束域**。它不是所有这些jet的量子制备／读取权限，也不是全局有限能解或量子态存在定理。$A_0$ 是共同描述的任意初始连接值，$u$ 保持原主约束乘子自由，不能为救坏系数而限制它们。

## 4. 主约束保持先选择S、C及接触项

在 $\Pi^0=G_R=0$ 上，主约束保持要求

$$
0=\Delta_a:=Q_{C,a}-Q_{R,a}-K_{S,a}
+\phi^\dagger\frac{\partial M(A_0)}{\partial A_0^a}\phi.\tag{8}
$$

整体负号来自 $\dot\Pi^0=-\delta H_c/\delta A_0$，不影响零条件。该缺陷不含 $E$。式(7)因此把“对全部合格Gauss初值弱零”提升为式(8)对任意局部 $\phi,\pi,A_0$ 的恒等式。

先取 $A_0=0$，二次项导数消失，得

$$
\Delta_a=2\operatorname{Re}\{\pi[i(C_a-R_a)-S_a]\phi\}=0.
$$

任意复 $\pi,\phi$ 迫使 $i(C_a-R_a)-S_a=0$。其中两部分分别为反Hermitian和Hermitian，故

$$
C_a=R_a,\qquad S_a=0.\tag{9}
$$

再让 $A_0$ 任意，式(8)化为

$$
2A_0^b\phi^\dagger\left[\tfrac12\{C_a,C_b\}-K_{ab}\right]\phi=0,
\qquad K_{ab}=\tfrac12\{C_a,C_b\}.\tag{10}
$$

这不是单凭 $\Pi^0=0$ 删掉 $S$，也不是预先把 $A_0$ 乘子形式写入正确答案。正则Hamiltonian由开放的式(2)导出，坏系数留下可核的主约束漂移。

## 5. 次约束保持、完整历史生成元及1016接续

式(9)—(10)后，$H_c=H_0+\int A_0^aG_{R,a}$，物质空间动能成为 $|D_i\phi|^2$，$D_i\phi=\partial_i\phi-iA_i^aR_a\phi$。给定式(1)后，除势外各项在原Gauss变换下不变。次约束在壳上保持因此要求

$$
\delta_RV=0.\tag{11}
$$

证明仍只用局部投影：剩余缺陷不含 $E$，任意 $\phi$ 可由式(7)补齐；再用任意紧支参数得到点式(11)。若容许新增次级约束，例如只准许某些 $\phi$ 取值，便不再是本轮原初值域的接续。

反过来，式(1)、(9)—(11)使完整原主—次链保持，不固定 $u$，无新增约束。完整时间依赖的历史变换由包含 $\Pi^0$ 和 $G_R$ 的组合生成，例如本符号下

$$
\mathcal G[\epsilon]=G_R[\epsilon]+\int d^3x\,(D_0\epsilon)^a\Pi^0_a,
\qquad \delta A_0=D_0\epsilon.\tag{12}
$$

不是把一个孤立次约束的固定时作用当成完整历史规范变换。Legendre可投影、total Hamiltonian与完整生成元的区别见[Pons §§2、4、8.2](https://arxiv.org/html/physics/0409076)；本轮不证明一般Dirac猜想或另加extended-H动力学。

此时原式(2)恰进入1016已经证明的局部作用类别：$S_{\rm eff}=0,C=R,K=\{C,C\}/2,\delta_CV=0$，表示关系继承式(1)。可直接复用其作用模合法边界的充分性和物质顶点结论。不再把相同五项分类称为新发现。整个桥仍是**明示正则局部有效类别内的接续条件定理**；E-RED及桥接角色的采用由主线另行决定。

## 6. 有限核验、删弱与来源范围

代码从式(2)的时间L直接数值微分求物质动量，分别核原／新辛势的 $dB$、Legendre变换、固定动量下 $A_0$ 导数及式(8)。错误 $C,S,K$ 分别在可补Gauss的资料上给非零主约束漂移；非不变正质量势给次约束漂移。正确族另核任意局部参数的一阶方向变分。非零常连接中的解析电场补齐也被核验，非仅把已知正确系数代回自身公式。所有有限数据是代数／接续证书，不是量子实验或新器件。

删弱边界：

- 只保一条零物质轨迹或少量物理读数，坏顶点可不可见；物理通道在Gauss子空间平凡不够。
- 允许额外约束或任意缩小物质初值域，可以避开上述极化；允许固定／改换 $A_0$ 描述也改变合同。
- 固定全局边界荷时，局部补齐不自动延成全局完成；不得调用该情况下未证的满射。
- 若删去同一Legendre及边界字典，可把 $S$ 的来源藏进动量或改名电荷；该操作不属于本轮固定比较。
- 放开原局域阶数，可加约束平方等弱零项而保同一物理过程；本轮不选所有未约化Hamiltonian代表。

凡把作用作分部积分或正则运输，其控制、记录与几何来源都应由同一全泛函运输，不能只保正则概率而省略来源差。旧[920](../../archive_894_922/research_note_920.md)、[1029](../../archive_1009_1043/research_note_1029.md)分别提供固定字典及固定作用的提升范围；它们不自动替本轮授予新量子量测、反常消除或实际共同设备。
