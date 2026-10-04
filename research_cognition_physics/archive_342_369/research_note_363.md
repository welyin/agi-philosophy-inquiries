# 第363轮：边界能量能否驱动内部认知——Marolf 条件、低能压缩与区域映射

## 0. 本轮问题与结果

本轮以冻结至360轮的成果为基线，与361、362轮独立。已读取根 README、研究目录 README、research_direction、RESEARCH_STATE，回读[359轮](research_note_359.md)、[360轮](research_note_360.md)并检索旧轮中的边界 Hamiltonian／运动学非局域性。未发现此前已完成本轮接口的审计。先开展原文和代码，待357—360整合完成通知后才建立本笔记。

问题是：**如果完整能量在物理扇区中由边界可观测量生成，同时微观“内部操作”保持严格局域，内部还能有非平凡的同钟动力学吗？**

本轮得到四个层次明确的结果：

1. Marolf 原定理确实给出一个受限冻结结论；它不需要 Lorentz 对称或短程 Hamiltonian，但需要共同时间与相容的边界／内部映射。
2. 低能投影不是免费通过此障碍的通行证。两个原本对易的算符，压缩到子空间后可以不对易；必须追踪保持物理扇区的条件。
3. 构造真正由边界产生的非平凡编码动力学，但可变化的逻辑相位不能由严格内部算符重建，其最佳算符范数误差至少为1。
4. 给出能量代表误差、算符越出码空间及时间窗口的联合界；小误差不保证无限长时间冻结。

这些是**边界生成元与操作映射的有限矩阵结论**。本轮没有生成时空或证明某个矩阵模型就是引力，也没有证伪整个认知生成纲领。

| 层次 | 内容 |
|---|---|
| 认知动机 | 整体能量和内部变化必须在同一操作描述中记账 |
| 额外输入 | 预先指定的空间支集、共同时间、Hamiltonian、物理码空间及候选边界能量 |
| 文献接口 | Marolf 的有条件内部冻结定理 |
| 本轮解析结果 | 压缩对易子缺陷、一般误差界、一参数反例、非局部编码与重建下界 |
| 数值验证 | 15项检查，包含完整局部代数、一般随机矩阵及未知参考 |
| 未得到 | 引力边界能量律、其认知来源、自然编码、Einstein 方程或全息对偶 |

## 1. 原始文献的精确条件

直接读取[Marolf，2014预印本](https://arxiv.org/abs/1409.2509)及[2015年 PRL 114, 031104 的 accepted manuscript](https://link.aps.org/accepted/10.1103/PhysRevLett.114.031104)，以正式接受稿的定义 I、II 和定理为准。

这里的“能量普适耦合”采用一个强的操作定义：在使 Hamiltonian 存在的边界条件下，总能量是边界引力场及其导数的局域积分，通量须对边界条件允许的规范变换不变。Hamiltonian 可以显含时间：

$$
H_{\mathrm{phys}}(t)
=\int_{\partial\Sigma_t}\mathcal F_{\mathrm{grav}}(t),
\qquad
\mathcal F_{\mathrm{grav}}\ \text{为规范不变的玻色可观测量}.
\tag{1}
$$

这不是“某个守恒荷可在边界读出”，也不是只令物质能量满足边界关系：H 必须生成所讨论的完整时间演化，包括引力自身的能量贡献。

其余关键条件是：底层存在背景定义的等时划分；不同位置的规范不变局部算符在至少一个为玻色时等时对易；底层与有效理论使用同一个时间演化；有效边界通量映射后仍支撑在微观边界附近。定理结论是该极限中远离边界的局部观察不再随时间变化。

正式稿同时明确：不要求平移／Lorentz 对称，不排除瞬时远程相互作用；它不处理整个时间演化本身涌现的方案，也不因此限制严格线性的自旋2激发。故本轮不同于359轮的表示论与应力协变性障碍。

本项目使用式(1)时必须提供自己的来源，不能把“所有资源都在整体内”直接改写成“所有能量必为规范不变边界通量”。前者是资源闭合要求，后者是强得多的生成元结构。

## 2. 原定理在可观测代数上的核心

设同一物理描述中的边界代数与内部代数在等时对易，且边界生成元确实等于完整 H。对没有另外显式时间参数的内部观察 O，

$$
\begin{aligned}
[H_\partial(t),O(t)]&=0,\qquad H(t)=H_\partial(t),\\
\frac{dO(t)}{dt}&=i[H(t),O(t)]=0.
\end{aligned}
\tag{2}
$$

这里取 ℏ＝1。该结论无需假设一个严格张量分解，只需对应可观测代数的交换关系和正确的生成元。因而360轮展示的“区域不独立因子化”**本身还不是 Marolf 条件被破坏的证明**；必须检查具体物理算符。

相反，仅把 H 改成长程相互作用不会自动破坏等时运动学局域性：共同酉演化总保持

$$
[U^\dagger A U,U^\dagger B U]=U^\dagger[A,B]U.
\tag{3}
$$

程序用一个跨两因子的相互作用 H 验证：不同因子的等时观察仍对易，但不同时刻的观察可以不对易且各自发生变化。它说明“空间远程动力学”和“等时可观测代数非局域”是不同条件。

**有界矩阵只是本轮的精确接口。** 连续场论中的无界算符、无限空间边界、约束量子化和极限域问题，不能由有限矩阵单例宣称全部处理完毕。

## 3. 低能投影的陷阱：压缩不保持乘积

令 V 把有效空间等距嵌入运动学空间，P 为其投影，Q 为补投影，φ 为算符压缩。直接插入 I＝P＋Q 得到

$$
\begin{aligned}
V^\dagger V&=I,\quad P=VV^\dagger,\quad Q=I-P,
\quad\phi(A)=V^\dagger A V,\\
[\phi(A),\phi(B)]
&=\phi([A,B])
-V^\dagger A Q B V+V^\dagger B Q A V.
\end{aligned}
\tag{4}
$$

因此原空间中的对易子为零，不保证压缩后为零。对 Hermitian 算符还可精确计入越出码空间的部分：

$$
\begin{aligned}
\|\,[\phi(A),\phi(B)]-\phi([A,B])\,\|
&\le2\|QAV\|\,\|QBV\|,\\
\phi(B^2)-\phi(B)^2
&=V^\dagger BQB V\succeq0.
\end{aligned}
\tag{5}
$$

式(5)第二行不只是代数细节：原测量的第一矩可以压缩为 φ(B)，原测量的第二矩却是 φ(B²)，不能任意替换成 φ(B)²。直接把 φ(B) 的谱测量当成同一个微观测量，会改变仪器的统计。

一个足以恢复精确冻结的子空间条件是：

$$
\begin{gathered}
QHV=0,\qquad QH_\partial V=0,\qquad
V^\dagger(H-H_\partial)V=0,\qquad [H_\partial,O]=0,\\
h:=V^\dagger H V,\quad o:=V^\dagger O V
\quad\Longrightarrow\quad
[h,o]=0.
\end{gathered}
\tag{6}
$$

第一项令 h 真实生成封闭码空间的动力学，第二项令边界能量代表保持该空间，第三项匹配能量，第四项保留空间局域性。由式(4)可见第二项已经消去对易子压缩缺陷，不需要再额外假设 O 保持码空间。若边界与内部算符本来均是保持物理扇区的观察，这些要求自然更加明确。

需要区分两种“能量相同”：只要求 PH P＝PH_∂P，是矩阵元匹配；要求 H_∂本身为保持物理空间的边界可观测量，还需检查其越出空间的部分。后者不能由前者推出。

## 4. 近似接口：误差、码空间保持与时间窗口

以下为本轮独立的有限维命题，不将它归为原文已经给出的近似定理。取自主 Hermitian H、H_∂、O，固定 V，全部范数为算符范数，定义

$$
\begin{aligned}
h&=\phi(H),\quad o=\phi(O),\\
\epsilon&=\|\phi(H-H_\partial)\|,\\
\beta_\partial&=\|QH_\partial V\|,\qquad
\beta_O=\|QOV\|,\qquad
\lambda_H=\|QHV\|,\\
c&=\|\phi([H_\partial,O])\|.
\end{aligned}
\tag{7}
$$

用式(4)、三角不等式和对易子范数界，得到

$$
\begin{aligned}
\|[h,o]\|
&\le r:=c+2\beta_\partial\beta_O+2\epsilon\|o\|,\\
\|e^{iht}o e^{-iht}-o\|
&\le |t|r.
\end{aligned}
\tag{8}
$$

若还要把 h 的演化解释成完整 H 在码空间的实际近似，则必须计入 λ_H。由 Duhamel 积分及 HV−Vh＝QHV，

$$
\begin{aligned}
\|e^{-iHt}V-Ve^{-iht}\|
&\le |t|\lambda_H,\\
\|\phi(e^{iHt}Oe^{-iHt})-e^{iht}o e^{-iht}\|
&\le2|t|\,\|O\|\lambda_H,\\
\|\phi(e^{iHt}Oe^{-iHt})-o\|
&\le |t|\bigl(r+2\|O\|\lambda_H\bigr).
\end{aligned}
\tag{9}
$$

算符变化另有平凡上界2‖O‖。这些是对全部码空间初态统一的界；张量上任意不参与演化的参考后，等距映射差的算符范数不变。因此含相关参考的纯态误差及观察期望也受相同控制。它不允许把上界外推到无限 t，更不说明误差常数在连续极限中自动有界。

这里 c 为区域映射缺陷，β 为算符离开指定空间的幅度，ε 为压缩能量误差，λ_H 为实际演化的空间保持误差。它们不是相同概念，也不能仅检查其中一项就声称完成引力映射。

## 5. 一参数精确反例：能量矩阵元全部吻合仍不够

两比特中第一个标为边界、第二个标为内部。取0＜θ＜π/2，构造两个正交编码向量：

$$
\begin{aligned}
V_\theta|0\rangle
&=\frac{|00\rangle-|01\rangle}{\sqrt2},\\
V_\theta|1\rangle
&=\cos\theta\,\frac{|00\rangle+|01\rangle}{\sqrt2}
-\sin\theta\,|10\rangle,\\
A&=Z\otimes I,\qquad B=I\otimes Z,\qquad [A,B]=0.
\end{aligned}
\tag{10}
$$

直接压缩得

$$
a=
\begin{pmatrix}1&0\\0&\cos2\theta\end{pmatrix},
\qquad
b=
\begin{pmatrix}
0&\cos\theta\\
\cos\theta&\sin^2\theta
\end{pmatrix}.
\tag{11}
$$

所以

$$
\begin{aligned}
\|[a,b]\|&=2\sin^2\theta\cos\theta>0,\\
\|QAV_\theta\|&=|\sin2\theta|,\\
\|QBV_\theta\|&=\sin\theta\sqrt{1+\cos^2\theta}.
\end{aligned}
\tag{12}
$$

这不是随机找到的一次例外，而是整个参数区间的显式族。给该空间一个真正不变、隔离的低能扇区：

$$
\begin{aligned}
H_\theta&=PAP+4Q,\qquad
\operatorname{spec}H_\theta=\{\cos2\theta,1,4,4\},\\
QH_\theta V_\theta&=0,\qquad
\phi(H_\theta)=\phi(A)=a,\\
\phi(H_\theta-A)&=0,\qquad
\|(H_\theta-A)V_\theta\|=|\sin2\theta|.
\end{aligned}
\tag{13}
$$

低能的所有能量矩阵元与边界 A 精确吻合，内部压缩 b 却发生变化。但 A 不保持码空间，B 也不保持码空间；H_θ 在全空间中不是边界算符。于是**它不是 Marolf 定理的反例**，而是“仅检查压缩能量匹配”的反例。

误把压缩当同态，会产生一个明确的错误动力学预言：

$$
\begin{aligned}
\phi(e^{iAt}Be^{-iAt})&=b,\\
\|e^{iat}b e^{-iat}-b\|
&=2\cos\theta\,|\sin(t\sin^2\theta)|.
\end{aligned}
\tag{14}
$$

第二行确实是 H_θ 在不变码空间中的演化，却不是由原边界 A 演化后再压缩的结果。取 cosθ＝1/√3，对易子范数为4/(3√3)≈0.769800；t＝0.9时错误预测的算符差约0.651993。矩阵元匹配误差仅约6×10⁻¹⁶，不能用它解释这个有限偏差。

测量乘积也给出独立操作见证：

$$
\begin{aligned}
\phi(B^2)&=I,\qquad
\operatorname{spec}b=\{1,-\cos^2\theta\},\\
\|\phi(B^2)-b^2\|&=1-\cos^4\theta.
\end{aligned}
\tag{15}
$$

在 b 的负本征态上，原 B 测量第二矩为1，直接测量 b 的谱所得第二矩为 cos⁴θ；上述角度时分别为1与1/9。若改用保持空间的完整算符 PBP＋QBQ，它已不再只有严格内部支集：程序检验其与边界 X 的对易子非零。真实可测映射必须保留这些差别。

## 6. 正构造：边界能驱动逻辑相位，但该相位不能继续叫严格内部局部观察

取 n≥2 个比特组成指定链，首位为边界，其余为内部。选一个重复编码及能量惩罚：

$$
\begin{aligned}
V|0\rangle&=|0\rangle^{\otimes n},\qquad
V|1\rangle=|1\rangle^{\otimes n},\\
H_\partial&=\frac12Z_0,\\
H_{\mathrm{pen}}&=\frac32\sum_{i=0}^{n-2}(I-Z_iZ_{i+1}),\\
H&=H_\partial+H_{\mathrm{pen}},\qquad
HV=H_\partial V=V\frac Z2.
\end{aligned}
\tag{16}
$$

两个码态能量为±1/2，离开该扇区的最小能量为5/2，故相对上端能量的间隔为2。边界生成元真实保持物理码空间，不只匹配其期望值。

等距编码保持任意未知逻辑态及其外部参考。逻辑相位可以随时间变化：

$$
\begin{aligned}
X_L&=V^\dagger\left(\prod_{i=0}^{n-1}X_i\right)V=X,\\
e^{-iHt}V&=V e^{-itZ/2},\\
X_L(t)&=\cos t\,X-\sin t\,Y,\qquad
\langle X_L(t)\rangle_{|+\rangle}=\cos t.
\end{aligned}
\tag{17}
$$

这里动态的 X_L 需要跨过边界的完整串算符。任意严格内部算符 I_0⊗M 的压缩均为对角矩阵：其两个码态之间的矩阵元含边界重叠〈0|1〉＝0。因此

$$
\begin{aligned}
\phi(I_0\otimes M)&\in\operatorname{span}_{\mathbb C}\{I,Z\},\\
\inf_M\|X-\phi(I_0\otimes M)\|&=1,\\
[Z/2,\phi(I_0\otimes M)]&=0.
\end{aligned}
\tag{18}
$$

下界1只需考察矩阵作用于|0〉的非对角分量即可证明，M＝0取到。故再复杂的严格内部测量，也不能把重建误差压到零。这是区域映射必须改变的具体代价，而不只是把逻辑比特更名为“内部”。

与360轮相比，本例不重新论证拼接或共有通量，而是增加真实边界 Hamiltonian、孤立低能扇区、非平凡逻辑演化及重建精度下界。全部严格内部压缩观察确实冻结，非平凡动态观察跨边界，恰与冻结命题相容。

链、编码、惩罚和边界 H 全是输入。这个例子没有产生引力的体积、度规或 Gauss 律，也不证明任意认知系统自然选择这类编码。边界能产生某种完整逻辑动力学，尚远不足以成为全息引力模型。

## 7. 近似冻结不能交换时间极限

即使没有子空间问题，近似边界生成元也只给受控窗口。令

$$
\begin{aligned}
H_\epsilon&=\frac12 Z\otimes I+\frac\epsilon2 I\otimes X,\qquad
O=I\otimes Z,\\
\|H_\epsilon-H_\partial\|&=\frac\epsilon2,\\
\|O(t)-O\|&=2\left|\sin\frac{\epsilon t}{2}\right|
\le\epsilon|t|,\\
t=\epsilon^{-1}
&\quad\Longrightarrow\quad
\|O(t)-O\|=2\sin(1/2)\approx0.958851.
\end{aligned}
\tag{19}
$$

固定时间下 ε→0 的冻结成立；观察时间同步增长到1/ε时，变化仍为有限量。重标定慢时间可以保留动力学，但此时必须重新核对所采用的同一时间生成元和边界能量极限。本轮不借此否定原定理，也不重启微观时间均匀性问题。

## 8. 数值检查、结果及可复算范围

代码：[emergent_boundary_dynamics_audit.py](363/emergent_boundary_dynamics_audit.py)；结果：[emergent_boundary_dynamics_audit_results.json](363/emergent_boundary_dynamics_audit_results.json)；交付检查：[research_round_363_checks.json](363/research_round_363_checks.json)。

15项 Checks 覆盖完整内部 Pauli 代数冻结、远程 H 下等时与非等时对易的区别、一般非对易矩阵的压缩恒等式、一参数弱能量匹配反例、真实与错误演化对比、一般误差界、码空间泄漏与参考、精确编码动力学、全部内部算符压缩秩、重建下界、低能隙、非一致时间极限和测量第二矩。

| 对象 | 数值结果 |
|---|---:|
| 两比特原边界／内部对易子 | 0 |
| 一参数压缩恒等式残差 | 不超过1.1×10⁻¹⁶ |
| cosθ＝1/√3时压缩对易子范数 | 0.769800 |
| 同点压缩能量误差 | 5.98×10⁻¹⁶ |
| 同点边界代表越出码空间幅度 | 0.942809 |
| 同点 t＝0.9 时错误演化预测差 | 0.651993 |
| 同点测量乘积缺陷 | 8/9 |
| n＝2、3、5 时全部严格内部压缩代数维数 | 2 |
| 同上完整逻辑代数维数 | 4 |
| 边界生成逻辑演化的等距误差 | 数值为0 |
| t＝0.8时动态逻辑 X 的算符变化 | 0.778837 |
| 严格内部重建逻辑 X 的最佳误差 | 1 |
| ε＝0.2、0.04、0.008，t＝1/ε时变化 | 都为0.958851 |

一般界另用确定种子的6—9维随机矩阵核对；一般量词由式(4)—(9)证明，随机样本不是证明本身。整个模型为有限矩阵，没有隐藏空间无限大、连续场论或自由度数极限。

运行方式，在本目录执行：

~~~powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 -m unittest emergent_boundary_dynamics_audit.Checks -v
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 emergent_boundary_dynamics_audit.py --write-results
~~~

最终15项检查均通过后写结果。补入测量乘积检查前的14项临时结果由保存保护拒绝直接替换，已完整保存在[本轮临时结果](363/drafts/emergent_boundary_dynamics_audit_results_pre_square_check.json)，再正常生成最终结果；这是未冻结本轮的完善过程，没有覆盖旧轮。复用 Python 3.12.14、NumPy 2.3.5，无新增依赖或图像检查。

## 9. 对认知生成路线的约束与下一步

本轮排除一条具体捷径：**“把某个低能 Hamiltonian 的矩阵元写成边界矩阵元，就可直接沿用原微观局部观察，并宣称已有非平凡内部引力。”** 需要同时满足物理空间保持、正确的算符乘积和同钟动力学；只核能谱或第一矩不足。

同时得到一个正向设计接口：若采用编码关系产生物理区域，动态的内部引力观察可能必须在微观描述中带到边界的支集；这种区域映射需要被实际构造，而不是从微观对象的普通复量子组合中自动读取。360轮的区域代数审计与本轮的能量生成元审计在这里相接，但二者仍不是 Einstein 约束的推导。

后续应对一个具体候选编码同时给出：物理低能空间、有效区域的可观测代数、完整能量的可测边界代表、这些算符的空间保持误差，以及随尺度与观察时间的控制。若缺少引力能量律，只得到有趣的编码动力学，应按其实际范围记录；若需要改变时间概念，应另立可检验模型，不能仅以“时间也涌现”一句话取消现有约束。

**当前总目标既未完成证明，也未被这一定理整体证伪。** 认知原则是否会选择所需的非局部物理映射、边界能量结构及有效时空，仍是开放问题。
