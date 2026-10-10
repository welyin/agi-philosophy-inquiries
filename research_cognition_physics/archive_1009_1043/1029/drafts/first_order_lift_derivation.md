# 1029推导资产：固定顶点后的一阶Noether提升等价

2026-10-08。仅承担本轮推导、科学代码与结果；正式报告／导航／交付收据另由根代理管理。新增科学校准1组，累计3806，新认知公理0，整体生成目标未完成。

## 1. 问题与严格结论

1028的同父对象桥采用了完整最小一阶变换及引力括号。本轮不分类所有第一阶相互作用，而固定同一份

$$
S_0=\sum_aS_{\rm PF}[\gamma^a]+S_{\rm SM,flat}+S_{\sigma},
\qquad
S_1=\sum_aS_{{\rm EH},3}^a+\frac12\sum_aq_a\int\gamma^a_{\mu\nu}T^{\mu\nu}.
\tag{1}
$$

物质部分、无导数Yukawa／声明阶Weinberg、可选B993门户以及1027领先源类别保持。物质相互作用不随形式引力参数 $\varepsilon$ 关闭；平直Lorentz结构、字段、非退化动能、$S_1$顶点均仍是输入。

设两个正规局域一阶提升分别满足

$$
\delta_0S_1+\delta_1S_0=0,
\qquad
\delta_0S_1+\delta^{\rm std}_1S_0=0
\quad(\bmod\ \text{边界}).\tag{2}
$$

其差对任意局部参数 $\xi$ 是 $S_0$ 的齐次规范对称。核验下述正规条件后，Noether第二定理的完备生成结论给

$$
\boxed{\;
R_1(\xi)=R^{\rm std}_1(\xi)+R_0(Z\xi)+M_\xi E_0.
\;}\tag{3}
$$

$Z$允许依赖字段及其有限jet，并含任意固定有限阶的参数导数；$M_\xi$为分部积分意义下的graded反自伴局域算符。结论是模原规范生成元与Euler平凡变换的等价，**不是离壳逐点唯一**。正文第4节给出作用于 $O=\Phi^\dagger\Phi$ 的显式非零平凡提升。

成熟依据是[Barnich—Brandt—Henneaux](https://arxiv.org/html/hep-th/0002245)第5.1.3节的独立／依赖Euler坐标、第6.6节定理6.9和推论6.3。第6.4.3节说明标准模型的正常局域结构；这里还须逐项核验所采用父作用，不能以“完备生成元按定义完备”代替证明。

## 2. 本父作用为什么满足正规条件

以下在完整、未删分量的局部字段jet中进行，不先取unitary规范、径向截面或把规范场置零。规范轨道在特殊字段值可有稳定子，不等于下述Euler主部秩降低。

### 2.1 规范、标量与一阶Weyl方程

把规范Euler式写成 $E_A^\mu$，实化Higgs及 $\sigma$ 写成 $E_h,E_\sigma$，Weyl和其独立共轭写成 $E_\psi,E_{\bar\psi}$。选择独立式：

- 全部标量Euler式及其导数；
- 全部Weyl／共轭Euler式及其导数；
- 规范 $E_A^i$ 的全部导数；
- $E_A^0$ 及其纯空间导数。

它们分别解出下列主jet及相应递推导数：

$$
\partial_0^2h,\ \partial_0^2\sigma;
\qquad \partial_0\psi,\ \partial_0\bar\psi;
\qquad\partial_0^2A_i;
\qquad\partial_1^2A_0.\tag{4}
$$

规范和标量的系数是非零常动能矩阵；Weyl时间主符号是 $i\sigma^0=iI_2$（共轭式符号相应改变）。它是可逆的一阶主部，**不是**二阶费米速度Hessian。若用Hamilton写法，Weyl代数动量约束属二类，并不带来新规范参数；graded jet写法直接求解奇时间jet，其可逆矩阵具有非零常数body。

规范流至多含一阶标量和零阶费米子；无导数Yukawa、Weinberg、Higgs／$\sigma$势及门户只改变上述递推的低阶项。它们不改变主pivot，无需Yukawa满秩、质量不简并或Higgs背景非零。对任何固定有限jet阶，先消去已确定的低阶／时间物质jet，再递推规范式，坐标变换保持三角和可逆。

完整内部Noether身份为

$$
D_\mu E_A^\mu+E_h\cdot t_Ah
+E_\psi\cdot T_A\psi+\text{共轭项}=0
\tag{5}
$$

（整体符号随既定Euler约定一致选择）。它以系数1解出 $\partial_0E_A^0$；连续微分后解出所有含时间导数的 $E_A^0$。依赖式与这些身份的变换三角、对角1，恰好提供剩余依赖，不留另一份任意局部参数的Euler身份。

超荷 $U(1)$ 保留，并且有带荷物质。本论证没有使用“纯半单群”假设，也没有把自由Abelian场的额外特征上同调类别误当成这里的新局部规范生成元。全局守恒流乘任意 $\xi(x)$ 一般会留下参数导数项，不能据此制造齐次局部对称。

### 2.2 自由Pauli–Fierz块

各独立PF块可用线性ADM变量检查。空间速度的DeWitt主形式在六个对称空间分量上可逆：归一为

$$
L_{\rm time}=\frac12\big(v_{ij}v_{ij}-(\operatorname{tr}v)^2\big)
\tag{6}
$$

时，六维Hessian行列式为 $-16$。未约化的迹方向负号不表示物理螺旋度出现新鬼；约束和规范仍须保留。

四个约束的主形式可取

$$
\mathcal C_0=\partial_i\partial_j\gamma_{ij}-\Delta\gamma_{ii},
\qquad\mathcal C_i=\partial_j\pi_{ij}.\tag{7}
$$

沿空间方向1，可取独立pivot $\gamma_{22,11}$、$\partial_1\pi_{11}$、$\partial_1\pi_{21}$、$\partial_1\pi_{31}$；对应系数矩阵为 $\operatorname{diag}(-1,1,1,1)$。所有空间延拓同样可解。线性Bianchi身份 $\partial_\mu E^{\mu\nu}_{\rm PF}=0$ 恰好把 $\partial_0E^{0\nu}_{\rm PF}$ 及其延拓变成依赖式。空间速度与动量间为可逆线性变换，故该说明也给Lagrangian jet的正规坐标。

PF与物质在 $S_0$ 是直和。合并各块的独立Euler坐标即可；跨块反对称Euler组合属于式(3)的平凡项，不增加非平凡身份。此处采用metric-PF与固定平直自旋标架，原生成元是PF及YM。局部Lorentz在对称vierbein写法中是自旋提升的补偿，不是额外自由传播规范内容；如改用未定规vierbein，须连同原有纯标架冗余处理，而不能只对 $\psi$ 任意局部转动并保持外部自旋矩阵不变。

### 2.3 正规性及结论的边界

限定平滑／局部形式、每阶有限jet和参数微分阶，不允许 $1/E$、逆微分算符、奇异参数投影或改变动能秩的附加菜单。完整作用的无导数相互作用没有破坏此限定。一般高导数主部、无穷非局域算符或退化动能须另审；本轮不将它们全部排除。

在上述graded正规坐标中，局部壳上为零的函数可按独立Euler坐标及其延拓展开；Noether第二定理把齐次对称与Euler身份对应。因此式(3)是固定父作用的推论，不是偷偷再加入“其它规范对称不存在”的新物理原则。有限矩阵校准仅检查列出的主系数；完备量词由本节拆分与成熟定理承担。

## 3. 必须运输完整生成元，不能只在反冲资料上丢 $E_0$

考虑一份任意正规完成

$$
S(\varepsilon)=S_0+\varepsilon S_1+\cdots,
\quad R(\varepsilon)=R_0+\varepsilon R_1+\cdots,
\quad E(S)=E_0+\varepsilon E_1+\cdots.\tag{8}
$$

第一步作整个参数空间的正规形式重定义

$$
\xi=P(\varepsilon)\eta,
\qquad P=I-\varepsilon Z+O(\varepsilon^2).\tag{9}
$$

单位主项使它按 $\varepsilon$ 可形式递归求逆；每个固定阶仍是有限jet。这里没有声称整个形式级数收敛或任意高能有效。变換整个生成元后，一阶的 $R_0Z$ 被移除；场依赖的 $P$ 也必须一并运输参数括号，包括变分作用在 $P$ 的附加项，不能只将旧括号左右乘矩阵。

第二步从已运输的生成元减去完整平凡对称

$$
\varepsilon M_\eta E(S),\tag{10}
$$

而不是只减 $\varepsilon M_\eta E_0$。graded反自伴性使式(10)对完整作用保持平凡；其一阶恰消去式(3)的 $M E_0$，其后继 $-\varepsilon^2M E_1$ 等项进入新的高阶生成元。参数括号与开放Euler项也按该等价变化共同运输。作用本身和物理等价类未变，所需一阶代表成为 $R_1^{\rm std}$。

1028的首反冲资料满足 $E(S)=O(\varepsilon^2)$，但其自由PF式只有 $E_{0,\rm PF}=O(\varepsilon)$。所以裸项 $\varepsilon M E_0$ 可在二阶非零；不能把“Euler平凡”误读成在此资料上立即为零。式(10)才是 $O(\varepsilon^3)$，它把原二阶贡献保留在被运输的 $R_2$ 中。这是规范代表等价，不是忽略真实反作用。

同一论证应用到内部参数：固定的 $S_1$ 在原YM行动下不变，因此任何一阶内部修正也是 $S_0$ 的齐次规范对称。可以在同一完整参数空间规范化到

$$
R_{0,\rm int}O=0,\qquad R_{1,\rm int}O=0.\tag{11}
$$

这里只需零阶和一阶，不再要求所有未知高阶内部行动都原样不变。

## 4. 真正的离壳非唯一例及反冲对照

实化Higgs为 $h\in\mathbb R^4$，$O=h^Th/2$，实超荷生成元 $t_Y$ 反对称且与内部作用相容。令

$$
b_\xi=\xi^\mu\partial_\mu O,
\qquad\Delta h=b_\xi t_Y E_h,
\qquad\Delta\text{其余字段}=0.\tag{12}
$$

由于 $h,E_h,b_\xi$ 都是偶量，

$$
E_h^T\Delta h=b_\xi E_h^Tt_YE_h=0,
\qquad\Delta O=b_\xi h^Tt_YE_h
\tag{13}
$$

后者一般非零。这对任意参数 $\xi(x)$ 恒成立，是完整 $S_0$ 的平凡规范对称，不是只在所选背景上的偶合。

代码采用 $h=(0,5/4,0,0)$、$\partial_yh=(0,2/7,0,0)$，选择真正的离壳加速度jet使

$$
E_h=(1/3,-2/5,1/7,3/2),\quad b_{\partial_y}=5/14.
\tag{14}
$$

具体令全部空间二阶导数为零、$A=\psi=\sigma=0$，并取 $\partial_t^2h=-E_h-\nabla U$，于是从真实Euler式 $E_h=\Box h-\nabla U$复得指定值。校准给 $E_h^T\Delta h=0$，而

$$
\Delta O=-75/224.\tag{15}
$$

没有声称这份离壳资料满足完整规范方程；它恰用于否定“裸提升唯一”。

另为反冲计阶，取一个真实跨Higgs／PF块的反对称平凡算符：

$$
\Delta h_i=b_\xi h_i\eta_{\mu\nu}E_{\rm PF}^{\mu\nu},
\quad
\Delta\gamma_{\mu\nu}=-b_\xi\eta_{\mu\nu}h_iE_{h_i}.
\tag{16}
$$

两块Euler收缩恒相消，$\Delta O=b_\xi h^2\operatorname{tr}E_{\rm PF}$。1028径向分支在校准点有

$$
\operatorname{tr}T=-(f')^2-4U_{\rm rad},\qquad
E_{0,\rm PF}=\varepsilon r,\quad r=-qT/2.
\tag{17}
$$

故 $\varepsilon M E_0$ 的 $O$ 变化在二阶确实非零。运输后的 $-\varepsilon^2M E_1$精确保留此阶，减去完整 $\varepsilon ME(S)$仅剩三阶。代码使用实际径向应力和首反冲系数；未知完整二阶Euler系数保留为符号 $c$，不假装已求得某个非线性完成。

## 5. 完整引力一阶括号为何不再独立选定

先完成第3节正规化，再比较一阶闭合。标准代表给纯PF参数的引力分量

$$
(B_{1,\rm std}^{\rm grav})^c
=A^c{}_{ab}[\xi^a,\zeta^b],
\qquad A^c{}_{ab}=\delta_{ab}\delta_{ac}\kappa_a.
\tag{18}
$$

任意另一个闭合括号与其差记为 $C^a(\xi,\zeta)$。减去两份闭合身份，保留开放项，给

$$
2\partial_{(\mu}C^a_{\nu)}\approx0,\tag{19}
$$

$\approx$表示对完整 $S_0$ 的无限延拓壳取商。内部YM参数不作用于metric-PF块，所以不能帮助式(19)。

不能把一个特殊Killing参数当成新自由：$C$须对**任意局部参数函数**成立，并是其有限阶局部微分表达。将 $\zeta$ 及其jets暂作系数，考察最高 $\xi$ 导数阶。PF参数主符号为

$$
v_\nu\longmapsto k_\mu v_\nu+k_\nu v_\mu.\tag{20}
$$

对任何非零协向量 $k$，选 $k_a\ne0$，对角式先给 $v_a=0$，混合式再给所有 $v_b=0$；光样 $k$ 也没有例外。于是最高参数jet系数在壳上为零。它及导数在延拓壳上均为零，再递降参数阶，得到

$$
C^a(\xi,\zeta)\approx0.\tag{21}
$$

这证明所需的是弱等价，而非一个未经取商的唯一括号。正规Euler坐标使弱零系数是 $E_0$及有限导数的局部组合；1028首反冲上它们为 $O(\varepsilon)$，因此其在二阶 $R_1O$ 项中的影响为三阶。不存在把任意局部参数投影到非零全局Killing场的有限局部算符；若靠边界积分或逆微分来做，已离开本合同。

完成此规范化后，1028二阶障碍只需式(11)、(18)—(21)，不再额外规定未知高阶内部行动。$B_{2,\rm int}$只经 $R_{0,\rm int}O=0$作用；$B_{1,\rm int}$经 $R_{1,\rm int}O=0$作用；$B_2^{\rm grav}$经 $R_0O=0$作用。Killing与首反冲引理继续排除剩余正规修补。

## 6. 恒等守恒改进与常量源：扩大必要性的自然推论

式(1)的 $S_1$还可加

$$
\frac12\sum_a\int\gamma^a_{\mu\nu}
\big(I_a^{\mu\nu}+C_a\eta^{\mu\nu}\big),\tag{22}
$$

其中 $I_a$只依固定物质jets而不含 $\gamma$，为正规局部、对称、内部规范不变的张量，且 $\partial_\mu I_a^{\mu\nu}\equiv0$；$C_a$是真正常数。它们不是仅在某个特殊解上守恒。

新增项的 $\delta_0$变分只是边界，原标准 $R_1$仍是一个提升。固定这份扩大的 $S_1$ 后，同样适用式(3)；首反冲源 $q_aT+I_a+C_a\eta$在真实物质解上守恒。由于 $I_a$不含引力字段，物质Euler项在 $\gamma=\varepsilon\gamma_1$上仍从二阶开始。因而1028动态权的必要条件保持：

$$
q_aq_b-\delta_{ab}\kappa_aq_a=0.\tag{23}
$$

这移除1028对**必要性**采用 $I=C=0$的限制。没有证明任意 $I$都有高阶完成，没有禁止其它引力理想通过非最小结构产生作用，也没有证明每项物理响应只依一份几何。$q=0$只表示这类最小一阶动态源为零。

同组校准取规范不变量 $O$与非零改进

$$
I^{\mu\nu}=(\partial^\mu\partial^\nu-\eta^{\mu\nu}\Box)O.\tag{24}
$$

令 $\Phi=(0,f)/\sqrt2$、$f=1+t+xy+z^2/2$，以精确多项式求全部分量和散度；无需任何方程便得到 $\partial_\mu I^{\mu\nu}=0$。这是一项真实非零改进的身份校准，不是一般改进分类。

## 7. 可复算交付与停止线

科学代码 `first_order_lift_equivalence.py` 默认只读重算对照 `first_order_lift_equivalence_results.json`；提供 `run()`、`compare()`，`--write`仅排他首次创建。历史源码只读导入1027的有理表示工具，未重跑旧科学扫描、未改历史文件。

本组包含：由动能形式得到的41维规范空间／标量速度主块；Weyl一阶符号；PF六维DeWitt及四约束pivot；含光样点的PF参数符号；式(12)真Euler平凡反例；式(16)—(17)完整Euler运输；非零恒等守恒改进。它们校准证明中的独立接口，不组成穷举全部局域变换的求解器。

实际收益是：**固定父作用与一阶顶点后，不需再将完整一阶提升和作用于 $O$ 的相关括号另当独立物理输入。** 它们可在同一等价类中规范化，1028必要条件的适用范围相应扩大。没有选定顶点本身、非最小参数、自由模式、共同Lorentz主部、维数或完整量子理论。

本轮结束后停止全集 $S_1$分类和更多等价样本；回到独立输入账选择真正的上游问题。科学1组，累计3806，认知公理新增0，整体目标保持未完成。
