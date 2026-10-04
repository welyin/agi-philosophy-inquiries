# 第372轮：量子位正锥、Lorentz因果锥与过滤操作的代价

日期：2026-09-23。完整独立轮次，科学基线冻结至367，不依赖370—371轮。18项检查通过。已读研究导航，并全档检索twistor、spinor、Lorentz cone、Pauli与行列式相关内容。

## 1. 成熟接口及本轮结论

用户提出“两个回波时间加端口”的测量结构。一个可以直接对接的数学事实是：**复量子位的未归一化正算符锥，恰与3＋1维Minkowski未来锥线性同构；其秩一射线构成一个二维方向球，SL(2,C)合同变换给出Lorentz群的双覆盖。**

这比单凭参数数量猜测三维更具体，但须区分：

1. 数学上，两个谱权重和一个量子位效应方向如何组织成Lorentz锥？
2. 物理上，为什么这些量必须代表事件位移和光线，而不是状态与测量的内部参数？

本轮完成第一个接口及其操作审计。一般boost在固定量子归一化规则下不是确定性信道；作为实际过滤实施，必须缩放、保留成功／失败结果，并计入依赖未知输入的成功率。这个实施限制不禁止把同一公式作为被动坐标变换；二者的物理含义不同。

| 层次 | 本轮处理 |
|---|---|
| 认知动机 | 少量可测数据、方向效应与主体间描述变换，能否给因果几何 |
| 已有数学 | Herm₂(C)与Lorentz锥、SL(2,C)与SO⁺(1,3) |
| 本地推导／复算 | 正性、行列式、双覆盖、谱权重、过滤概率、任意参考及非仿射反例 |
| 额外输入 | 指定二能级对象、Pauli标架、过滤强度、空白结果寄存器和读出 |
| 未提供的桥梁 | 将事件位移识别为正算符、将效应方向识别为传播方向、物理钟尺与平移 |
| 未得到 | 物理三维的必然性、自然Lorentz动力学、场论、Einstein方程 |

## 2. 旧结果和文献归属

[183轮](../archive_171_189/research_note_183.md)已研究任意维spin factor的Lorentz型锥、齐次性及单体结构选择，且明确没有由此推出物理时空。本轮不重复Jordan重建；新增的是**特定复量子位的旋量对应，以及锥自同构与合法量子仪器的区别**。

[早期阶段稿](../archive_217_222/_shared/article8_hypothesis.md)已指出正定Fisher度量的实拉回不能给Lorentz号差，并要求归一化过滤补全测量。本轮用的是行列式二次型，不是将Fisher度量改成不定度量；旧障碍仍保留。

实际核读的原始研究：

- [Arrighi–Patricot，*A Note on the correspondence between Qubit Quantum Operations and Special Relativity*，J. Phys. A 36 (2003)，Lemma 1及Propositions 1、3](https://arxiv.org/html/quant-ph/0212135)：给出正锥、旋量映射、过滤分支与缩放Lorentz变换的对应，并讨论归一化和测量解释。
- [Verstraete–Dehaene–De Moor，*Local filtering operations on two qubits*，式(1)及Theorem 2](https://arxiv.org/pdf/quant-ph/0011111)：局部过滤须满足A†A≤I，归一化分母不可删；Pauli系数与Lorentz作用相连。

这是一套已知接口。本轮将其适用条件、参考量词和资源账落实为可复算证据，不宣称发现了新的旋量表示。

## 3. 行列式准确给出1＋3号差

用σ₀＝I及三个Pauli矩阵作实基，写

$$
X=x^0I+x^i\sigma_i
=\begin{pmatrix}x^0+x^3&x^1-ix^2\\x^1+ix^2&x^0-x^3\end{pmatrix},
\qquad
x^\mu=\frac12\operatorname{Tr}(\sigma_\mu X),
\qquad
\det X=(x^0)^2-|\boldsymbol x|^2.
\tag{1}
$$

由(x·σ)²＝|x|²I可得

$$
\lambda_\pm(X)=x^0\pm|\boldsymbol x|,
\qquad
X\ge0\ \Longleftrightarrow\ x^0\ge|\boldsymbol x|,
\qquad
\rho=\frac{I+\boldsymbol r\cdot\boldsymbol\sigma}{2}
\Longrightarrow x^0=\frac12,\quad
\det\rho=\frac{1-|\boldsymbol r|^2}{4}.
\tag{2}
$$

“未来”首先是正锥分支的数学命名，还没有被识别为物理事件的未来。归一化态只取固定迹的三维球截面，不等于全部四维时空。

非零秩一正算符可用复二分量旋量或效应表示：

$$
\begin{aligned}
X&=\zeta\zeta^\dagger,\qquad
x^0=\|\zeta\|^2/2,\qquad
x^i=\zeta^\dagger\sigma_i\zeta/2,\qquad
(x^0)^2-|\boldsymbol x|^2=0,\\
E&=wP_{\boldsymbol n},\qquad
P_{\boldsymbol n}=\frac{I+\boldsymbol n\cdot\boldsymbol\sigma}{2},
\quad|\boldsymbol n|=1,\quad0<w\le1,\qquad
\mathbb{CP}^1\simeq S^2.
\end{aligned}
\tag{3}
$$

旋量整体相位不改变X，去掉正尺度后得到S²。效应的w≤1来自E≤I，射线尺度不能当作不受限单次概率。S²是内部测量方向球；解释为实际光线方向仍需操作映射。

## 4. SL(2,C)给出Lorentz群，但先是代数作用

令detA＝1。合同作用及对应实矩阵为

$$
X'=AXA^\dagger,\qquad
\Lambda(A)^\mu{}_\nu
=\frac12\operatorname{Tr}(\sigma_\mu A\sigma_\nu A^\dagger),
\qquad
x'=\Lambda(A)x,\qquad
\Lambda(AB)=\Lambda(A)\Lambda(B).
\tag{4}
$$

正性和行列式由合同作用保持，对式(1)极化即得度规保持。极分解A＝HU中，H正定且detH＝1，U属于SU(2)，二者都可连续接到I。因此

$$
\Lambda(A)^T\eta\Lambda(A)=\eta,\qquad
\eta=\operatorname{diag}(1,-1,-1,-1),\qquad
\det\Lambda(A)=1,\qquad
\Lambda(A)^0{}_0=\tfrac12\operatorname{Tr}(AA^\dagger)\ge1.
\tag{5}
$$

两个显式子类固定符号：

$$
\begin{aligned}
U&=e^{-i\theta\boldsymbol n\cdot\boldsymbol\sigma/2}
&&\Longrightarrow\quad
\Lambda(U)=\operatorname{diag}(1,R_{\boldsymbol n}(\theta)),\\
B&=e^{\chi\boldsymbol n\cdot\boldsymbol\sigma/2}
&&\Longrightarrow\quad
\Lambda(B)=
\begin{pmatrix}
\cosh\chi&\sinh\chi\,\boldsymbol n^T\\
\sinh\chi\,\boldsymbol n&
I_3+(\cosh\chi-1)\boldsymbol n\boldsymbol n^T
\end{pmatrix}.
\end{aligned}
\tag{6}
$$

每个proper orthochronous Lorentz矩阵把时间单位向量送到某个未来单位类时向量；一个上述boost可做到同样的事。消去boost后，剩余变换固定时间轴，必为SO(3)旋转，而轴角形式由SU(2)覆盖，故映射满射。

若Λ(A)恒等，则AIA†＝I使A酉，继而A与全部Hermitian矩阵对易，只能为标量；detA＝1再给A＝±I。这是二对一覆盖，不是一对一相等。

## 5. 两个非负数与方向：用户猜想的候选结构

给定s₊≥s₋≥0和一个效应方向n，有

$$
X=s_+P_{\boldsymbol n}+s_-P_{-\boldsymbol n}
=tI+r\,\boldsymbol n\cdot\boldsymbol\sigma,\qquad
t=\frac{s_++s_-}{2},\quad r=\frac{s_+-s_-}{2},\quad
\det X=s_+s_-=t^2-r^2.
\tag{7}
$$

两个数选定谱，方向选定本征投影；s₋＝0落在类光边界。这与“两标量＋方向”的形式吻合。

若s₊、s₋要取为回波钟读数或其差值，还须说明时间原点、单位、传播校准，以及为什么这些操作量等于式(7)的本征值。本轮没有从回波记录自动得到S²；它由指定量子位效应提供。一个无内部方向结构的端口标签，也不能直接替代n的全部测量内容。

## 6. 确定性量子合同变换只留下旋转

固定归一化效应I后，X↦AXA†是完全正映射，但保迹要求

$$
\operatorname{Tr}(A\rho A^\dagger)=1\quad
\text{对全部 }\rho
\quad\Longleftrightarrow\quad A^\dagger A=I.
\tag{8}
$$

对A∈SL(2,C)，这只留下SU(2)旋转。一般boost的某些输入输出迹超过1，因此不是固定量子位上的确定性CPTP信道。

合法CPTP信道也通常不是Lorentz变换。例如γ＝0.4的振幅阻尼将I／2送到Bloch向量(0,0,0.4)，行列式从0.25变为0.21。CP／TP与行列式保持不能互相替代。

### 6.1 实际过滤须计成功和失败

对任意可逆A，取最大允许尺度并补足失败分支：

$$
K=cA,\qquad c=\frac1{\|A\|},\qquad
L=\sqrt{I-K^\dagger K},\qquad
K^\dagger K+L^\dagger L=I,\qquad
p_{\rm s}(\rho)=\operatorname{Tr}(K\rho K^\dagger).
\tag{9}
$$

更小正c也合法。成功分支的未归一化四向量是c²Λ(A)x，行列式乘c⁴；保留实际概率时不能删掉缩放。

对χ≥0的boost，负χ可通过翻转n处理：

$$
c=e^{-\chi/2},\qquad
p_{\rm s}(\boldsymbol r)
=e^{-\chi}\big(\cosh\chi+\sinh\chi\,\boldsymbol n\cdot\boldsymbol r\big),
\qquad e^{-2\chi}\le p_{\rm s}\le1.
\tag{10}
$$

每个有限χ下所有输入的成功概率非零，但最坏值随χ增大而下降。未知输入失败后不能假定可免费重新准备同一个未知态；这里只计一个完整试次，不假设复制或无限免费重试。

取n＝z、χ＝ln2：

$$
A=\operatorname{diag}(\sqrt2,1/\sqrt2),\qquad
K=\operatorname{diag}(1,1/2),\qquad
L=\operatorname{diag}(0,\sqrt3/2),\qquad
p_{\rm s}=\frac58+\frac38 r_z.
\tag{11}
$$

北极、南极、中心态的成功率分别为1、1／4、5／8。这是实际过滤的输入依赖，不能当成观察者之间纯粹重新标记坐标。

## 7. 任意未知参考：完整分支一致，不是无扰动

对任意有限R和任意ρ_AR，带记录instrument为

$$
\begin{aligned}
\mathcal I(\rho_{AR})
&=|{\rm s}\rangle\langle{\rm s}|\otimes
(K\otimes I_R)\rho_{AR}(K^\dagger\otimes I_R)\\
&\quad+|{\rm f}\rangle\langle{\rm f}|\otimes
(L\otimes I_R)\rho_{AR}(L^\dagger\otimes I_R),\\
\operatorname{Tr}\mathcal I(\rho_{AR})&=1,\qquad
\operatorname{Tr}_A\!\left[\sum_{b={\rm s,f}}\rho_{AR}^{\,b}\right]=\rho_R .
\end{aligned}
\tag{12}
$$

完全正性来自Kraus形式，任意R量词不由有限样本代替。最后一式使用完备性；**单独成功分支的条件参考一般改变**。

Bell输入和式(11)给

$$
|\Phi\rangle=\frac{|00\rangle+|11\rangle}{\sqrt2},\qquad
p_{\rm s}=\frac58,\qquad
|\Phi_{\rm s}\rangle=\frac{2|00\rangle+|11\rangle}{\sqrt5},\qquad
\rho_{R|{\rm s}}=\operatorname{diag}(4/5,1/5),\qquad
\frac58\rho_{R|{\rm s}}+\frac38|1\rangle\langle1|=\frac I2 .
\tag{13}
$$

成功参考相对I／2改变迹距离0.3，平均参考不变。删掉失败样本又不报告成功标签，会误报一个不存在的确定性过程。

### 7.1 完整有限实现与资源

一个空白两能级结果寄存器足以承载两分支。代码同时核验合同算符的标准酉扩张：

$$
\mathcal U_K=
\begin{pmatrix}
K&-\sqrt{I-KK^\dagger}\\
\sqrt{I-K^\dagger K}&K^\dagger
\end{pmatrix},\qquad
\mathcal U_K^\dagger\mathcal U_K=I,\qquad
\mathcal U_K(|0\rangle\otimes|\psi\rangle)
=|0\rangle K|\psi\rangle+|1\rangle L|\psi\rangle .
\tag{14}
$$

酉性用Kf(K†K)＝f(KK†)K。保留整个相干结果寄存器时是等距过程；读成经典分支则给式(12)，读出与环境记录不再可被省略。

这是指定有限操作的实现；控制耦合、χ、寄存器准备和结果存储均需资源。没有由此选出自然Hamiltonian或实际钟尺。

## 8. 归一化不能修补确定性boost

将AρA†除以其迹给成功条件态，所得映射一般不仿射。式(11)中两极纯态分别保持，但中心态变成diag(4／5,1／5)：

$$
\widehat{\mathcal B}_A(\rho)
=\frac{A\rho A^\dagger}{\operatorname{Tr}(A\rho A^\dagger)},\qquad
D\!\left(
\widehat{\mathcal B}_A(I/2),
\frac{\widehat{\mathcal B}_A(|0\rangle\langle0|)
+\widehat{\mathcal B}_A(|1\rangle\langle1|)}2
\right)=0.3.
\tag{15}
$$

这不是未知混态上的线性确定性信道。作为带条件的测量更新没有问题，不能删掉条件事件及其概率。

若保留所有试次再忽略结果，确定性通道为

$$
\rho=\begin{pmatrix}a&b\\b^*&d\end{pmatrix}
\longmapsto K\rho K^\dagger+L\rho L^\dagger
=\begin{pmatrix}a&b/2\\b^*/2&d\end{pmatrix}.
\tag{16}
$$

它是相位退相干，既不是未缩放boost，也不是成功条件映射。

## 9. 被动改变表示可以保持概率，但归一化效应也要改变

有序线性空间的表示变换可以连同效应和单位效应一起变：

$$
X'=AXA^\dagger,\qquad
E'=A^{-\dagger}EA^{-1},\qquad
u'=A^{-\dagger}IA^{-1},\qquad
\operatorname{Tr}(E'X')=\operatorname{Tr}(EX),\qquad
\sum_bE'_b=u'.
\tag{17}
$$

归一化由Tr(u'X')＝1表达。一般u'不等于I，因此这不是在固定归一化仪器上主动实施无失败boost。数值实例的概率配对误差2.78×10⁻¹⁷，但u'与I算子距离1.014。

式(17)给“不同认知描述等价”的具体候选；是否就是物理惯性观察者变换，仍要连同事件、钟尺及全部测量验证。物理Lorentz协变并不要求每个有限量子位单独承载确定性boost信道，故过滤障碍不反驳相对论。

## 10. 不能由锥同构直接跨越的两道选择问题

### 10.1 归一化态之间没有非平凡正序

若直接把归一化状态看成事件，并用正算符次序定义未来，就会退化：

$$
\operatorname{Tr}\rho=\operatorname{Tr}\sigma=1,\qquad
\sigma-\rho\ge0
\Longrightarrow
\operatorname{Tr}(\sigma-\rho)=0
\Longrightarrow \sigma=\rho .
\tag{18}
$$

要把Loewner序用作因果序，必须改为未归一化增量，或以Herm₂为模型空间的事件仿射空间，并解释操作意义。若预先定义“Y在X未来当且仅当Y−X≥0”，式(1)立刻给Minkowski因果锥；但事件集合、平移和位移识别是新增前提。

### 10.2 为什么选择量子位

复矩阵结构允许不同大小。这里选择指定二能级对象或可操作二能级部门；不能反推所有基本事件都必须如此表示，也不预设原合同保证每一种整数维对象实际存在。

更高矩阵正锥一般不再是普通Lorentz锥。以n≥3为例：

$$
P_1=|1\rangle\langle1|,\qquad
P_2=|2\rangle\langle2|,\qquad
\operatorname{rank}(P_1+P_2)=2<n .
\tag{19}
$$

两个非比例极射线之和仍在正锥边界。Lorentz锥中两个非比例未来空射线之和却在内部；线性锥同构保持极射线和边界，因此整个Herm_n⁺不能一概识别成某个更高维Lorentz锥。

已有量子位部门提供局部数学接口；为何它代表基本因果方向，复合后怎样拼接成同一个时空，仍需证明。

## 11. 数值结果与范围

| 检查 | 结果 |
|---|---:|
| 行列式二次型误差 | ≤9.11×10⁻¹⁶ |
| 本征值公式误差 | ≤1.11×10⁻¹⁵ |
| 合同作用与四向量作用之差 | ≤4.45×10⁻¹⁶ |
| 旋转＋boost的Lorentz度规残差 | 7.78×10⁻¹⁶ |
| Λ(A)与Λ(−A)之差 | 0 |
| χ＝ln2：北／南／中心成功率 | 1、1／4、5／8 |
| Bell条件参考改变迹距 | 0.3 |
| 全分支平均参考改变迹距 | 1.67×10⁻¹⁶ |
| 完整酉扩张误差 | 1.12×10⁻¹⁶ |
| 归一化boost非仿射见证 | 0.3 |

一般锥、群、CP和参考量词由解析推导承担，抽样不替代证明。一般酉扩张交叉检验使用略低于最大尺度的合同矩阵，避免饱和零本征值开方的舍入放大；解析式(14)适用于全部合同矩阵，式(11)的饱和对角实例另有验证。

## 12. 下一步及文件

值得沿用的候选是将可辨认的二元谱数据与秩一效应方向连接起来。但必须提出并检验“事件位移／可达增量为何是这种正算符”的桥梁；不能把迹权重命名为时间、Pauli方向命名为空间便宣称完成生成。

若能给出该桥梁，还需检验跨对象拼接、平移、钟尺、局域动力学及物质反作用。单个Lorentz锥和变换群不承担这些结论；G、Λ具体数值也不被另设为引力方程族重建门槛。

- [代码](372/qubit_lorentz_cone_audit.py)
- [结果](372/qubit_lorentz_cone_audit_results.json)
- [统一核验入口](372/research_round_372_checks.json)，由整合步骤生成。

~~~powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 -m unittest qubit_lorentz_cone_audit.Checks -v
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 qubit_lorentz_cone_audit.py
~~~

在archive_231_目录执行。18项检查通过后才保存结果；Python 3.12.14、NumPy 2.3.5；没有安装依赖、生成图像、修改旧轮次或导航。

