# 第683轮：同时间片物理读口、完整来源极限与辅助反射边界

2026-10-02。接[682](../../research_note_682.md)及[已执行入口](time_local_source_entry.md)。[代码](../joint_time_local_source_interface.py)、[结果](../joint_time_local_source_interface_results.json)、[核验](../research_round_683_checks.json)、[联合账](../unified_physics_condition_ledger_683.md)。两组新核验、十六式；主代理审查，无独立代理审查。

## 1. 问题、旧成果与假设

上一目标轮发布682并复算683入口，属于有效进展。本轮核对README、research_direction、RESEARCH_STATE、最新笔记与结果；进程列表未发现旧Python研究。CIM命令行枚举被系统拒绝，随后使用Get-Process核查，没有据此重启旧运行。

682证明辅助逆块可以非零甚至非退化，但保持全部原物理来源并不自动给正性。本轮沿实际物理时间接口推进：678原观测含B=2I+aX，因此仍跨相邻时间片。我们补出严格同时间片的辅助读口，证明它规范协变，并在原完整平均下回到同一物理来源极限；同时实际核查补偿部门，排除一种看似直接的反射接法。

|层次|内容|
|---|---|
|认知动机|同一操作读口须同时遵守时间支撑、规范关系及原状态权重|
|继承输入|原有限盒、16通道Wilson核、M、Pφ、H_b、双Haar和S⁹平均|
|新增有限调节选择|观测端点B替换为2I，质量来源同步替换；动能K不改|
|解析增量|层数一致的全来源误差；完整平均极限；实际规范字典；指定补偿反射的精确障碍|
|数值|原非平坦完整通道、实际局部群变换；原自由盒的解析矩阵交叉复核|
|未完成|真正辅助正实现、原Q₀正性、指定归一、H_F身份、共同连续及量子GR|

空间接口按[679 §1](../../research_note_679.md)逐项继承：382上界、383自由反向对合、384已消去额外Lipschitz的下界、386真实邻域及有限认证、425替代性的Lie坐标路线、522绝对差坐标热下界、523仪器与下界共同实现。既有证明不重做，不把386和425叠加；当前h、s、CAR/Gauss与旧实际端点、参考及仪器的身份映射仍单独记账。

## 2. 严格同时间片的读口

记γ₅=Q₊−Q₋、H=γ₅X，ε为677有限层软符号。沿678原约束，端点解为

$$
B=2I+aX,\quad
F_0=B^{-1}(I+\varepsilon)/2,\quad
F_L=B^{-1}(I-\varepsilon)/2,\quad
\varepsilon^\dagger=\varepsilon,\quad\|\varepsilon\|\le1.
\tag{1}
$$

原端点观测O=(B,0,…,−B)。仅在观测中改为Ō=2(I,0,…,−I)，令T为682式(9)原物理端点系数，‖T‖=1/2。新观测为

$$
\bar Y=
\begin{pmatrix}
\tfrac12J_-^T\psi+J_-^T(z_0-z_L)\\
-\tfrac34J_+^TM\psi+J_+^T\chi+
\tfrac12J_+^TM(z_0-z_L)
\end{pmatrix}.
\tag{2}
$$

它只用同一物理格点上的ψ、χ、z₀、z_L和局部M。第五方向标签不是物理时间，故正时间来源现在确实落在辅助正半区代数。积分后有效Φ̄仍可非局部；不能把这个支撑结论误写为原裸ξ字典局部。

代回端点后，相对于原Φ的差为

$$
\bar\Phi-\Phi=
\bigl(-aTB^{-1}X\varepsilon,\ 0\bigr),\qquad
\|\bar\Phi-\Phi\|\le\delta_a:=
\frac{a\|X\|}{2(2-a)}.
\tag{3}
$$

证明只用Re X≥−I给‖B⁻¹‖≤1/(2−a)、‖ε‖≤1及‖T‖=1/2。原两传播方向m₀=1给‖X‖≤3，各方向的互补自旋投影和酉链路块范数不超过1。它是当前模型方向数的输入，不是维数生成。另有

$$
\left\|\begin{pmatrix}F_0\\F_L\end{pmatrix}\right\|
\le\frac1{2-a},\qquad
\left(\frac{I+\varepsilon}{2}\right)^2+
\left(\frac{I-\varepsilon}{2}\right)^2
=\frac{I+\varepsilon^2}{2}\le I.
\tag{4}
$$

因此误差一致于L和背景Wilson谱隙；a是辅助有理调节，不是物理时间格距。

## 3. 同步质量、全部来源及完整平均

原质量依赖同一Y。必须同步使用Φ̄，不能只换测试来源。记γ=(1+√5)/2为670旧观测上界：

$$
\bar N=N_0+\lambda\bar\Phi^TP_\phi\bar\Phi,\qquad
\|\bar N-N\|\le
|\lambda|\|P_\phi\|\delta_a(2\gamma+\delta_a).
\tag{5}
$$

有限a时这份调节与原调节不完全相同。对固定有限阶来源Z，增广矩阵与系数定义为

$$
\mathcal A_Z=\begin{pmatrix}N&\Phi^TZ\\-Z^T\Phi&0\end{pmatrix},
\qquad C_Z=\operatorname{Pf}\mathcal A_Z;
\qquad \bar{\mathcal A}_Z\text{ 使用 }\bar N,\bar\Phi.
\tag{6}
$$

其固定偶数尺寸为q。对a≤a_*<1，元素界R=C(1+‖Pφ‖)，差值界Δ=C′a(1+‖Pφ‖)；C、C′可依赖盒、来源和λ，不依赖L、谱隙。逐乘积伸缩得

$$
|\bar C_Z-C_Z|
\le(q-1)!!\frac q2\Delta R^{q/2-1}
\le C''a(1+\|P_\phi\|)^{q/2}.
\tag{7}
$$

先按678精确积分辅助与补偿部门，再用此固定尺寸界，避免把随L增加的辅助维数放进粗Pfaffian常数。669全阶热矩及673完整双Haar、S⁹可积性直接给

$$
\int|\bar C_{a,L,Z}-C_{a,L,Z}|\,d\varpi_\tau\le C_{Z,\lambda,\tau}a,
\qquad
\int|\bar C_{a,L,Z}-C_{{\rm original},Z}|\,d\varpi_\tau\longrightarrow0
\quad(a\to0,\ aL\to\infty).
\tag{8}
$$

第二式复用677原来源极限，不新增统一Wilson谱隙。仅为固定有限盒、τ>0、固定来源阶数的未归一结论；不涉及归一化、体积一致、背景导数或物理连续时空。682不可见族也直接保留：Ō=2B⁻¹O使OΣOᵀ=0蕴含ŌΣŌᵀ=0。

## 4. 实际局部规范合同

U为原全16通道的局部变换，V为两分量手征框表示。ψ、z₀、z_L按U变换，χ按U*变换；局部M′=U*MU†。定义Uξ=diag(U,U*)、UY=diag(V,V*)，则

$$
\bar Y'\operatorname{diag}(U,U^*,U,U)=U_Y\bar Y,\qquad
\bar\Phi'U_\xi=U_Y\bar\Phi,\qquad
U_\xi^T\bar N'U_\xi=\bar N.
\tag{9}
$$

第一式逐块用原J±和M合同；第二式还用K、端点解协变；第三式同步输送原Pφ和质量。这是真实局部颜色、弱和超荷变换，不是只核全局相位。实际非平坦盒的最大元素残差约1.14×10⁻¹⁵。因此来源变形可以进入原完整Gauss平均，不需要重选规范群或物质表示。

## 5. 反射只在极限恢复，不能提前签收

679已证原全部来源反射。记Zθ为原反线性反序来源；新调节在有限a下没有同一精确证明，但由三角不等式及式(7)—(8)有

$$
|\bar C_{b^\theta,Z^\theta}-\overline{\bar C_{b,Z}}|
\le |\bar C_{b^\theta,Z^\theta}-C_{b^\theta,Z^\theta}|
+|\bar C_{b,Z}-C_{b,Z}|,
\qquad\int|\text{反射差}|\,d\varpi_\tau=O(a).
\tag{10}
$$

适用同一同步反射的背景及有界协变来源剖面。极限继承原Hermitian物理型，有限调节的正性没有因此成立。实际两个来源在a=.2、.1的反射差约2.805×10⁻²⁰、6.990×10⁻²²；这类微小系数只作绝对诊断，不靠倒权重或相对误差宣判物理正负性。

## 6. 原补偿矩阵并非成熟Wilson体的自动实例

重新核读[Kikukawa—Usui，§II及§VI](https://arxiv.org/html/1005.3751v3#S6)：其自由overlap和非规范Yukawa结果要求实际半区反射；文中另指出特定动态规范Pauli–Villars玻色作用存在反射障碍。这里不把文献的Wilson五维核等同本项目约束K，而直接计算后者。

678补偿复玻色作用为−b†Gb，G=K†K>0。指定组件逐点时间反射b↦P_t b、再对函数复共轭，不混合辅助层和自旋。其作用反射要求

$$
\mathcal P^\dagger G^\theta\mathcal P=G,
\qquad\mathcal P=I_{\rm layer}\otimes P_t\otimes I_{\rm spin,int}.
\tag{11}
$$

仅检查此明确选择，不宣称所有可能反射都已排除。原L=1时K=[[B,B],[−(B−aH),B+aH]]；B†H=HB给

$$
G=\begin{pmatrix}
C-2aHB&a^2H^2\\a^2H^2&C+2aHB
\end{pmatrix},\qquad C=2B^\dagger B+a^2H^2.
\tag{12}
$$

现在取原2×2盒的单位链路，仍保原16通道、时间反周期。令Q_t=[[0,1],[−1,0]]作用在时间指标，W=I−空间交换，A=γ₅γ₀Q_t。各恒等式都在原实际核内：

$$
X=W+\gamma_0Q_t,\quad H=\gamma_5W+A,\quad
A^2=I,\quad[W,A]=0,\quad\operatorname{spec}W=\{0,2\},\quad
P_tAP_t=-A.
\tag{13}
$$

将其代入(12)，得到精确差和算符范数：

$$
\mathcal P^\dagger G\mathcal P-G
=\operatorname{diag}\bigl(8a(I+aW)A,-8a(I+aW)A\bigr),
\qquad\|\mathcal P^\dagger G\mathcal P-G\|=8a(1+2a)>0.
\tag{14}
$$

这不是随机背景上的近似反例。a=.2、.1分别给2.24、.96，直接原矩阵复核误差小于2×10⁻¹⁵；G最小特征值仍为约7.315、7.629。原非平坦背景亦给反射差范数2.374、.987。正定高斯精度与指定反射不变是不同条件。

**所排除的是把这个组件逐点反射直接套到原补偿作用的接法。** 不是完整Gauss/S⁹物理泛函反例，更不是原自由overlap正性的反例。全部辅助变量积分后，原确定体因子仍严格抵消：

$$
\frac{\det K\,\det K^\dagger}{\det(K^\dagger K)}=1.
\tag{15}
$$

未积分的辅助部门不良性质不能在积分后无条件移给物理子代数；反过来该数值抵消也不能证明完整辅助反射正性。

## 7. 增量、复算与下一项

本轮同时减少C01、C04、C14、C16、C17、C20、C22之间的独立选择：新同时间片读口、质量、规范输送及完整平均只能共同变换。真正关闭的是678读口跨相邻时间片的支撑问题，在明确调节极限意义下连接回原物理候选。

$$
\text{严格半区读口＋原来源极限＋规范协变}
\quad\not\Longrightarrow\quad
\text{辅助正性或原 }Q_0\ge0.
\tag{16}
$$

第一组复算入口四个a、L端点字典、质量界及来源；第二组核原局部群、反射差和补偿精确自由反例。未构造L=1600巨大体矩阵，未数值声称完成全群积分；原热矩与来源多项式承担解析平均结论。既有Python/NumPy，无图像检查。旧入口原样保留。

下一项[684](../../684/drafts/STATUS.md)检验适配K的辅助反射是否能同时满足作用协变、反线性对合及正半区合同，尤其区分“作用有某种对称矩阵”与真正物理反射。候选失败后检查替代表示或直接物理子代数；不按同一种错误反射重复扫描a。四个共同分支、原Q₀、H_F与记录身份、连续和量子GR继续开放，目标不改。
