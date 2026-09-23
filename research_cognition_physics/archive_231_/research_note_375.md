# 第375轮：共同方向输运能确定多少连接，完整传播又能辨认多少

日期：2026-09-23。完整独立轮次，科学基线冻结至373；不使用374轮的结果。16项检查通过。已读根与研究README、research_direction、RESEARCH_STATE，并复核349、358、359、362、372、373及相关结果；先完成历史去重与原始文献核对，再进行本轮分类。

## 1. 新问题与结论

[373轮](research_note_373.md)在指定常系数Weyl演化中，把内部概率角与实际传播方向相连。现在增加一项明确条件：

> 不同位置的内部方向做平行输运时，必须与同一给定空间度量的向量输运相容。

这比“概率守恒”强，但仍需指定向量如何输运。本轮把它作为额外建模合同，不称为已从认知原则证明的要求。

结论分三层：

1. **固定Levi-Civita连接后，兼容的Hermitian旋量连接只剩一个中心相位一形式。** 若另外限为SU(2)连接，这个中心也被去掉。这里的无挠、给定空间度量和不可约二分量表示是前提。
2. **连接与完整传播算子不能重复计数。** 在三维非退化符号下，中心一形式经Dirac收缩是单射；但是再任意加入Hermitian零阶项时，其中三个分量已能重新分配给中心连接，不能当成另外三个独立物理量。
3. **不固定无挠条件时，完整传播甚至不能识别全部连接。** 本轮给出五维有挠连接族：连接及其曲率可以改变，完整Weyl算子却严格不变。因此，没有另行定义平行输运的独立读出时，不能把这些连接称为不同的可观测传播几何。

此外，一个有非零中心曲率的给定例子会改变Landau能级，说明“存在不可辨识连接”不等于“所有剩余连接都只是记号”。关键在于完整算子与实际操作，而非系数个数。

## 2. 去重、额外输入与文献边界

[349轮](research_note_349.md)已证明概率守恒只固定半散度项，保留任意Hermitian势；把这一式机械推广为三维，不足以构成新轮。本轮新增的是**兼容连接的仿射空间、连接与零阶项的重参数化、以及有挠连接到完整算子的五维不可辨识核**。

[358轮](research_note_358.md)的Levi-Civita条件来自已给定的无挠Palatini作用量及其变分，不能用它反推本轮认知动机已选择无挠。[359轮](research_note_359.md)讨论规范代表与物理曲率，[362轮](research_note_362.md)讨论Z₂带衣操作；二者没有完成本轮的微分连接分类。[360轮](research_note_360.md)的Wilson读出和[364轮](research_note_364.md)的平坦描述连接不重复做实验。

| 层次 | 本轮内容 |
|---|---|
| 认知动机 | 不同位置的内部方向与向量方向采用相容输运 |
| 额外输入 | 定向三维光滑空间、正定度量及体积、非退化不可约Clifford符号、指定向量连接；主定理再选Levi-Civita |
| 解析证明 | 中心差分类、收缩单射与重参数化、五维有挠核、完整参考扩张不变、Landau封闭块 |
| 数值验证 | 非平直兼容方程独立求解、曲率收缩、核维数、同算符演化、Landau谱及跃迁 |
| 未证明 | 这些结构的认知来源、三维选择、物理体积来源、无挠必然性、Maxwell或Einstein方程 |

直接使用的原始文献：

- [Nicolaescu，*Geometric connections and geometric Dirac operators on contact manifolds*，§1.1，Propositions 1.3、1.5、1.6](https://arxiv.org/pdf/math/0101155)：固定度量连接下的兼容连接空间，以及允许挠率后的Dirac等价。其§1.1的一般局部结果不要求本轮空间具有接触结构。本轮用Pauli矩阵独立给出所需三维证明和明确见证。
- [Fang–Vassiliev，*Analysis as a source of geometry: a non-geometric representation of the Dirac equation*，§5—6、§9](https://arxiv.org/html/1401.3160)：提供一阶算子与几何／电磁势的表示接口。原文采用四流形、正密度和非退化符号；其几何表示涉及Levi-Civita连接。**Remark 5.2中某个由主符号构造的修正项的唯一性仍是猜想，本文不采用该猜想。** 下文固定连接的分类由兼容方程直接证明。
- [Li–Roy–Das Sarma，*Weyl fermions with arbitrary monopoles in magnetic fields: Landau levels, longitudinal magnetotransport, and density-wave ordering*，§III的单Weyl子类](https://arxiv.org/pdf/1608.06632)：作为最后谱接口的已有原型。本文只复算线性、各向同性单节点的封闭块，不处理材料、输运系数或异常。

这些是对已有数学与物理接口的复算，不称为本项目新发现的一般定理。

## 3. 固定几何合同与兼容方程

在一个可平凡化的定向三维空间局部片上，给定实可逆三标架系数E与正定度量。使用Hermitian Pauli约定：

$$
A^j=E_{aj}\sigma_a,\qquad
g^{jk}=(E^{\mathsf T}E)_{jk},\qquad
A^jA^k+A^kA^j=2g^{jk}I,\qquad
d\mathrm{vol}_g=\mu\,d^3x,\quad
\mu=\sqrt{\det(g_{jk})}.
\tag{1}
$$

E的列编号为空间上指标j，行编号为内部方向a。这里“Clifford符号”就是式(1)这组矩阵及其乘法关系，并不是从矩阵乘法推导出了空间的存在。

先指定一个向量连接，约定∇_i v^j=∂_i v^j+Γ^j_(ik)v^k。旋量连接写为∇^S_i=∂_i+Ω_i。要求Hermitian内积与Clifford映射兼容：

$$
\Omega_i^\dagger=-\Omega_i,\qquad
\partial_iA^j+[\Omega_i,A^j]
+\Gamma^j{}_{ik}A^k=0.
\tag{2}
$$

它表示：先把向量／余向量映成内部矩阵再输运，与先用指定向量连接输运再映射相同。**指定Γ是合同的一部分。** 只知道g时，选择其Levi-Civita连接还用到了无挠条件；一般度量兼容连接可以带挠率。

局部正交标架的旋转连接可以通过so(3)与su(2)的二对一表示提升，给出式(2)的一组解。这里只讨论已经提供局部Clifford模的情形；不据此证明任意全球拓扑上已选定spin／spinᶜ结构。

## 4. 固定向量连接后，剩余恰为中心一形式

设Ω与Ω̃满足同一式(2)，令C_i=Ω̃_i−Ω_i。相减即得：

$$
[C_i,A^j]=0\quad\forall j.
\qquad
E\text{可逆}\ \Longrightarrow\
[C_i,\sigma_a]=0\quad\forall a.
\qquad
C_i^\dagger=-C_i
\ \Longrightarrow\
C_i=i a_iI,\quad a_i\in\mathbb R.
\tag{3}
$$

证明最后一步不需要数值抽样：任何2×2复矩阵都能写成c₀I+c·σ；与三个Pauli都对易迫使c=0，反Hermitian性使c₀纯虚。反过来，加上i a_i I不改变式(2)。因此：

$$
\Omega_i=\Omega_i^{\rm LC}+i a_i I
\quad\text{当固定Levi-Civita连接时};
\qquad
\operatorname{tr}\Omega_i=0
\ \Longrightarrow\ a_i=0
\quad\text{在选定SU(2)旋量标架中}.
\tag{4}
$$

每个空间方向有一个实中心自由度，不是整个问题只剩一个实常数。SU(2)限制是对行列式丛连接的额外限制，并非单靠相位重选总能实现；非零中心曲率无法这样消除。它不能由“存在共同方向”自动得到。

局部换标架必须同时改变符号、旋量、连接与观测表示。对静态U(x)∈SU(2)：

$$
\Psi'=U\Psi,\qquad
A'^j=UA^jU^\dagger,\qquad
\Omega'_i=U\Omega_iU^\dagger-(\partial_iU)U^\dagger.
\tag{5}
$$

于是式(2)保持成立。仅更换Ω而不改变对应数据，通常不是同一实验的换记号。对中心相位U=e^(iχ)I，同一约定给a'=a−dχ。中心曲率与完整曲率为：

$$
f_{ij}=\partial_i a_j-\partial_j a_i,\qquad
F^S_{ij}=\partial_i\Omega_j-\partial_j\Omega_i+[\Omega_i,\Omega_j]
=F^{S,\rm LC}_{ij}+i f_{ij}I.
\tag{6}
$$

局部相位换标架不能改变f。反过来，f=0只给局部纯规范结论；全球拓扑与整体相位输运须另行交代。本轮不做新的Wilson环实验。

## 5. 收缩后的算子：哪些参数重复，哪些仍可测

定义空间Dirac／Weyl算子：

$$
D_a=-iA^j(\partial_j+\Omega_j^{\rm LC}+i a_jI)
=D_{\rm LC}+A^j a_j,\qquad
\|A^j\delta a_j\|_{\rm op}^2
=\delta a_jg^{jk}\delta a_k.
\tag{7}
$$

正定性说明，在固定LC、符号与表示后，a到算子差的映射是单射。

对无边界局部支撑的试验态，式(2)与LC体积恒等式给：

$$
\mu^{-1}\partial_j(\mu A^j)+[\Omega_j^{\rm LC},A^j]=0,
\qquad
D_{\rm LC}^\dagger=D_{\rm LC}
\quad\text{作为形式微分算子}.
\tag{8}
$$

这里使用带μ的内积并通过分部积分；全局自伴域及边界条件是另一个问题。实中心a给Hermitian乘法项，保持形式自伴。**一般有挠度量连接不自动满足这个结论**；第7节单独检查其条件。

若又允许任意Hermitian零阶项V，则它与中心连接不能全部视作独立参数：

$$
V=\phi I+\mathbf b\cdot\boldsymbol\sigma,\qquad
\phi=\tfrac12\operatorname{tr}V,\qquad
\boldsymbol\alpha=E^{-1}\mathbf b,\qquad
D_a+V=D_{a+\alpha}+\phi I.
\tag{9}
$$

这是完整算子恒等式，不是只比较主符号或能谱。给定一个参考分解D_LC后，每个Hermitian零阶差可唯一写成中心一形式加标量；但若同时自由指定a和V，分解本身有三函数冗余。标量项不是相对论质量已被推导；在这个单Weyl部门中应按其实际作用解释。

算子相等还比量子通道相等强。对空间常数φ₀：

$$
e^{-it(D+\phi_0I)}=e^{-it\phi_0}e^{-itD},
\qquad
\mathcal U_{D+\phi_0I,t}\otimes\operatorname{id}_R
=\mathcal U_{D,t}\otimes\operatorname{id}_R
\quad\text{对任意参考 }R.
\tag{10}
$$

所以孤立传播通道也不能确定一个绝对常能量零点。位置依赖的φ一般不能以同样的整体相位消去；若引入可相干控制的演化分支或额外能量标定，其物理资源与生成元也须说明。

## 6. 非平直正例：独立恢复LC旋量连接

为避免只在常矩阵上验证交换关系，本轮选取一个明确的非平直空间度量：

$$
g_{ij}=e^{2f(\mathbf x)}\delta_{ij},\quad
A^j=e^{-f}\sigma_j,\quad
f=\tfrac12\mathbf x^{\mathsf T}Q\mathbf x,\quad
Q=
\begin{pmatrix}
0.2&0.07&0\\
0.07&0.4&0\\
0&0&0.6
\end{pmatrix}.
\tag{11}
$$

它是输入的背景，不由物质反作用求出。由度量导数得到Γ，再代入兼容方程，可解得：

$$
\begin{aligned}
\Gamma^j{}_{ik}
&=\delta^j_i f_k+\delta^j_k f_i-\delta_{ik}f_j,\\
\Omega_i^{\rm LC}
&=\frac i2(\mathbf e_i\times\nabla f)\cdot\boldsymbol\sigma,\\
D_{\rm LC}
&=-i e^{-f}\bigl(\boldsymbol\sigma\cdot\nabla
+\boldsymbol\sigma\cdot\nabla f\bigr),\\
R^{\rm LC}
&=e^{-2f}\bigl(-4\Delta f-2|\nabla f|^2\bigr).
\end{aligned}
\tag{12}
$$

程序不用第二行作为线性求解器输入：它从∂A及Γ，分别在iI、iσ₁、iσ₂、iσ₃基中解式(2)，得到每个方向秩3、中心核维数1；再与解析Ω比较。另独立计算向量曲率和旋量曲率：

$$
(R_{ij})^k{}_\ell
=\partial_i\Gamma^k{}_{j\ell}
-\partial_j\Gamma^k{}_{i\ell}
+[\Gamma_i,\Gamma_j]^k{}_\ell,\qquad
[F^S_{ij},A^k]+(R_{ij})^k{}_\ell A^\ell=0.
\tag{13}
$$

在x=(0.2,−0.3,0.4)处，直接缩并R约−4.32944854516；连接恢复误差约5.55×10⁻¹⁷，兼容残差约8.33×10⁻¹⁷。非零曲率说明这项检验没有依赖纯换基的平直背景。

代码还检查带μ=e^(3f)的形式伴随关系。忽略体积权重会得到另一套错误的连接条件；不能以349中的平坦概率测度公式直接替换这里的式(8)。

## 7. 放宽无挠：五维连接变化对完整传播不可见

现在明确改变合同：仍取平直g_ij=δ_ij和A^j=σ_j，但向量连接不再固定为LC。令b为任意实3×3矩阵，定义：

$$
\Omega_i=i\sum_a b_{ia}\sigma_a,\qquad
\Gamma^j{}_{ik}=2\sum_a b_{ia}\epsilon_{ajk},\qquad
T^k{}_{ij}=\Gamma^k{}_{ij}-\Gamma^k{}_{ji}.
\tag{14}
$$

Γ_i在向量指标j、k上反对称，所以它保持平直度量；直接计算[Ω_i,σ_j]验证式(2)。通常T不为零。因此这里是一对相容的有挠向量／旋量连接，**不是**同一LC连接下式(4)的反例。

其Dirac收缩可完整算出：

$$
D_b=D_0+C_b,\qquad
D_0=-i\boldsymbol\sigma\cdot\nabla,\qquad
C_b=-i\sum_i\sigma_i\Omega_i
=(\operatorname{tr}b)I
+i\sum_{i,a,k}b_{ia}\epsilon_{iak}\sigma_k.
\tag{15}
$$

对实b，最后一项反Hermitian；其中三个系数恰好检测b的反对称部分。因此：

$$
D_b\text{形式自伴}
\ \Longleftrightarrow\ b=b^{\mathsf T},\qquad
D_b=D_0
\ \Longleftrightarrow\
b=b^{\mathsf T},\quad\operatorname{tr}b=0.
\tag{16}
$$

这里的充要条件限于式(14)这个平坦常系数子类。对称矩阵有6维，去掉迹后有5维；从全部实b到C_b的实线性映射秩为4。代码直接构造这张矩阵求秩，并逐项核对五个正交的对称无迹基，不以几个随机b替代核空间证明。

尤其取一个明确非零例：

$$
b=\operatorname{diag}(1,2,-3),\qquad
\Omega_x=i\sigma_x,\quad
\Omega_y=2i\sigma_y,\quad
\Omega_z=-3i\sigma_z,\qquad
C_b=0.
\tag{17}
$$

这些系数为常数，旋量曲率来自对易子：

$$
F^S_{xy}=-4i\sigma_z,\qquad
F^S_{yz}=12i\sigma_x,\qquad
F^S_{zx}=6i\sigma_y,\qquad
R^{\rm LC}_{ijkl}=0.
\tag{18}
$$

Γ自身的有挠曲率也非零，并与式(18)满足式(13)。它不是平直度量的Levi-Civita曲率；把两者都叫“时空曲率”会制造假矛盾。

与此同时，D_b与D₀是逐系数相同的微分算子。在指定平直三维环面和周期spin结构上取共同自伴域后：

$$
D_b=D_0
\ \Longrightarrow\
e^{-itD_b}=e^{-itD_0}
\ \Longrightarrow\
(\mathcal U_{b,t}\otimes\operatorname{id}_R)(\rho)
=(\mathcal U_{0,t}\otimes\operatorname{id}_R)(\rho)
\quad\forall t,\ R,\ \rho.
\tag{19}
$$

因此不只是能谱相同：所有传播通道、任意参考关联和固定读出均相同。若中间插入同样的仪器序列，每段传播仍相同，完整实验记录也相同。

数值选四个环面Fourier动量，形成精确封闭的8维旋量子空间；它是同一常系数算子的约化子空间，不是假定有限网格具有严格正则对易关系。直接构造两个生成元、分别指数化，再作用于最大纠缠参考，差值均为0。任意输入与参考的量词由式(19)证明。

本例中连接给出的独立平行输运当然可以不同；但“能够单独读取该连接的平行输运”必须是额外的操作合同。若唯一给出的自然实验就是D生成的传播，就没有依据把这些连接解释为不同的可观测几何。这是本轮与“同一主符号但势不同”的旧式反例最主要的差别。

## 8. 中心曲率并非全部不可见：一个实际谱与演化接口

回到固定平直LC连接，取同一Weyl主符号并指定中心一形式：

$$
\mathbf a=(0,-Bx,0),\qquad B>0,\qquad
f_{xy}=-B,\qquad
\boldsymbol\Pi=-i\nabla+\mathbf a,\qquad
[\Pi_x,\Pi_y]=iB,\qquad
D_B=\boldsymbol\sigma\cdot\boldsymbol\Pi.
\tag{20}
$$

此处B是输入的中心曲率大小。称它为电磁形式接口，需要另行给出荷、场源及仪器实现；本轮没有推导Maxwell方程或引力源。

固定纵向动量k_z，使用横向振子表示，约定：

$$
c=\frac{\Pi_x+i\Pi_y}{\sqrt{2B}},\qquad
[c,c^\dagger]=1,\qquad
D_B=
\begin{pmatrix}
k_z&\sqrt{2B}\,c^\dagger\\
\sqrt{2B}\,c&-k_z
\end{pmatrix}.
\tag{21}
$$

对n≥1，子空间由上分量|n;B⟩与下分量|n−1;B⟩张成，保持不变：

$$
D_B\big|_n=
\begin{pmatrix}
k_z&\sqrt{2Bn}\\
\sqrt{2Bn}&-k_z
\end{pmatrix},\qquad
E_{n,\pm}=\pm\sqrt{k_z^2+2Bn},\qquad
D_B\binom{|0;B\rangle}{0}
=k_z\binom{|0;B\rangle}{0}.
\tag{22}
$$

这给出一个真正有可测谱差的中心自由度。不同B具有不同f与能级间隔，不能通过局部相位规范变换消除。它与第7节改变连接但完整D不变的情形不同。

若在对应B下制备上分量|n;B⟩，演化t后测内部下分量，则：

$$
P_{\downarrow}(t\mid n,B)
=\frac{2Bn}{k_z^2+2Bn}
\sin^2\!\left(t\sqrt{k_z^2+2Bn}\right).
\tag{23}
$$

代码从完整旋量—振子矩阵独立本征分解推进，再投影到下分量；没有用式(23)当作模拟输出。

**两种B下的磁长度与振子基随B变化。** 相同n标签不是同一个坐标空间波函数；这里是各自匹配场背景的制备。内部自旋读出可以相同，场源、谱标定与制备条件均为输入。这些概率不能冒称对同一未知空间输入比较得到的通道距离。主要可辨识证据是规范不变曲率及能级间隔；所有换标架比较仍须同步变换状态与读出。

数值截取N个振子态时，也没有把有限矩阵的正则对易关系冒充精确连续关系：

$$
[c_N,c_N^\dagger]
=I_N-N|N-1\rangle\langle N-1|.
\tag{24}
$$

只有n<N的完整二态块与上分量n=0模在本次表示中逐项认证；最高下分量出现的额外孤立态是截断伪影，明确剔除。把N从8增加到11时，认证块的谱保持一致。未计算无限横向简并、边缘态或完整异常。

## 9. 可复算结果

代码：[clifford_connection_compatibility_audit.py](clifford_connection_compatibility_audit.py)；结果：[clifford_connection_compatibility_audit_results.json](clifford_connection_compatibility_audit_results.json)；主代理冻结入口：[research_round_375_checks.json](research_round_375_checks.json)。

使用既有Python 3.12.14与NumPy 2.3.5，先运行Checks，16项全部通过后执行程序的--write-results；未安装依赖、未修改历史科学文件，未生成图像。

~~~powershell
& 'C:/Users/admin/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B -X utf8 'research_cognition_physics/archive_231_/clifford_connection_compatibility_audit.py'
~~~

复算保存追加--write-results；程序通过检查后保存，拒绝覆盖不同结果。

| 项目 | 结果 |
|---|---:|
| 固定向量连接：每方向4维反Hermitian未知量的方程秩 | 3 |
| 每方向中心核维数 | 1 |
| 非平直LC连接独立恢复最大误差 | 5.55×10⁻¹⁷ |
| 曲例的直接标量曲率 | −4.32944854516 |
| 一般实b收缩的秩／核维数 | 4／5 |
| b=diag(1,2,−3)的最大挠率分量 | 6 |
| 三个旋量曲率块算符范数，xy／yz／zx | 4／12／6 |
| 同一有挠例的LC曲率、完整生成元差、传播差 | 0／0／0 |
| 含最大纠缠参考的输出差 | 0 |
| a与V重分配后的算子差 | 0 |

Landau例取k_z=0.4、n=2、t=0.6、N=9：

| B | 正能级E₊ | 对应制备的下分量概率 |
|---:|---:|---:|
| 0.4 | 1.32664991614 | 0.46417363104 |
| 0.7 | 1.72046505341 | 0.69713456652 |

两者主符号与空间度量完全相同；跃迁计算与式(23)相差约10⁻¹⁶。这里的零与小误差只指明确矩阵计算，所需一般量词由式(3)、(7)、(9)、(15)—(19)、(21)—(24)证明。

16项检查分别覆盖：Clifford恒等式；中心核；曲连接恢复；无挠／度量相容／曲率；带体积的形式伴随；局部SU(2)标架变换；中心收缩的正定范数；零阶重参数化；常数能量的参考通道；有挠相容；五维核；非零有挠曲率；同算符含参考演化；Landau封闭块和跃迁；截断边界；非零中心曲率的谱差。

## 10. 本轮关闭与剩余接口

本轮关闭的是一个明确分类问题：**在已经给定几何与方向输运兼容性的条件下，哪些连接被迫相同，哪些差异只是同一传播算子的不同分解，哪些会改变实验。**

它没有从“互相认同”证明向量连接必为LC，也没有把本轮中心U(1)等同于标准模型的U(1)_Y。几何和可测作用量的来源、物质普适耦合、连接或度量的反作用仍需另外推导。

后续若要说“测量恢复了连接”，必须先明确测量恢复的是哪一个等价类：

- 仅由完整Weyl传播恢复算子，不能识别第7节的有挠核。
- 加入独立平行输运读出，可以提出更强的连接识别任务，但这种能力及资源必须来自明确模型。
- 已经选择LC、体积及最小连接耦合后，中心曲率可以通过指定传播谱辨识；这仍不提供该曲率或度量的场方程。

因此下一步有价值的是检验**新增输运能力或几何反作用如何由同一内部操作资源实现**，而不是继续换一个V展示同一主符号不唯一。若没有这一新的操作输入，应保留本轮条件分类，不以额外编号重复“不足以得到GR”的结论。
