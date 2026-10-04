# 第674轮：互补手征配对、完整规范边界与标量泛函的实性

日期：2026-10-02。接[673](research_note_673.md)、[实际入口](674/drafts/unnormalized_source_entry.md)。[代码](674/joint_gauss_reflection_reality.py)、[结果](674/joint_gauss_reflection_reality_results.json)、[核验](674/research_round_674_checks.json)、[条件账](674/unified_physics_condition_ledger_674.md)。两组新核验、十六式；主代理审查，无独立代理审查。

## 1. 本轮关闭哪一项

上一轮完成673并实际执行674入口，属于进展。已核README、direction、state、报告／结果和发布后核验；没有运行中的Python研究。回查660、668、670、672—673；382—386、425、522—523仍按667 §1.1继承，不恢复已经消去的空间假设。

673撤回了一般Hermitian／实性表述，因为有限样本不足以承担证明。本轮采用互补Pfaffian恒等式，给**原标量权重与局部质量来源**的解析反射共轭关系，再接到673的完整双边界Haar积分。任意外加Grassmann来源的反射、全多项式RP、正归一及原H_F身份仍分别开放。

成熟代数工具见[Knuth《Overlapping Pfaffians》§2，式(2.6)—(2.8)](https://arxiv.org/pdf/math/9503234)。下文写出本项目使用的互补子式形式和方向约定；不是把它当作新代数定理。只读取原文文字，无图像检验。

|层次|本轮地位|
|---|---|
|认知动机|同一历史及其反向组合必须使用同一共轭关系|
|继承输入|原有限Wilson／overlap、原16通道、S⁹配对、原质量及673矩形字典|
|解析增量|原权重跨手征秩的反射身份；奇异质量多项式延伸；完整标量Gauss核Hermitian及总权重实性|
|数值核验|互补子式、原强规范场、原质量双重关系及奇异Higgs质量|
|未完成|正性、非零配分函数、任意费米来源的反射合同、实际时间过程与连续／量子GR|

## 2. 反射会交换两个手征谱空间

R为原时间反射的格点置换，L=Rγ₀，T=γ₅L；L†=L、T†=−T、det T=1。对反射后的完整规范链路，

$$
D^\theta=LD^\dagger L,\qquad
H^\theta=-THT^\dagger,\qquad
(P_u^\theta,P_v^\theta)=(TP_vT^\dagger,TP_uT^\dagger).
\tag{1}
$$

因此可取uθ=Tv、vθ=Tu。令V=[u,v]、d=det V、维数为(rᵤ,rᵥ)，rᵤ+rᵥ=n=2r。偶数rᵤ、rᵥ时

$$
(r_u^\theta,r_v^\theta)=(r_v,r_u),\qquad
\det[u^\theta,v^\theta]=d,\qquad |d|=1.
\tag{2}
$$

原r为32的倍数。若rᵤ、rᵥ为奇数，673已证明全部纯物理来源前有奇数辅助积分，标量权重为零；反射后也为零，下文只需处理偶数情况。空块Pfaffian取1。

原辅助M(E)在固定原手征方向约定下满足

$$
M^{\mathsf T}=-M,\quad M^\dagger M=I,\quad
\operatorname{Pf}M=1,\quad T^{\mathsf T}M^\theta T=M,\quad
\det[J_+,J_-]=\det[J_-^*,J_+^*]=1.
\tag{3}
$$

单位性和手征配对继承660／670。Pf M的方向由原E=e₀固定，原Clifford配对的行列式恒为1、S⁹连通，故其符号不随E改变。格点反射只重排偶数尺寸块，不添加隐藏的方向符号。

## 3. 用互补子式代替逆传播子

记Aᵤ=Pf(uᵀMu)、Aᵥ=Pf(vᵀMv)。由可逆反对称矩阵的互补子式公式及M⁻¹=−M*，

$$
A_u=d\,\overline{A_v}.
\tag{4}
$$

该式不要求uᵀMu或vᵀMv可逆，包含它们的零模。

把实质量尺度吸收到原物理配对中：P=diag(p,q)，其中

$$
p^{\mathsf T}=-p,\qquad q=-p^*,\qquad
N_W=
\begin{pmatrix}
w^{\mathsf T}pw&-K_\ell^{\mathsf T}\\
K_\ell&q
\end{pmatrix},\quad
w=J_-^\dagger v,\quad K_\ell=J_+^\dagger v.
\tag{5}
$$

这里p包括原Dirac和Majorana矩阵，允许各格点不同。先在p可逆处定义

$$
\mathscr B=J_-^*pJ_-^\dagger+J_+^*q^{-1}J_+^\dagger,\qquad
a=\operatorname{Pf}q,\quad
Q_v=\operatorname{Pf}(v^{\mathsf T}\mathscr Bv),\qquad
F=\frac{A_u\,a\,Q_v}{d}.
\tag{6}
$$

这是673的矩形N_W对固定q块作Schur积分后的同一权重。临时使用的是质量q的逆，既没有用Kₗ⁻¹、D⁻¹，也没有除以F。

原质量P=diag(p,−p*)和T的手征交换给

$$
T^{\mathsf T}\mathscr B^\theta T=-(\mathscr B^*)^{-1},\qquad
\operatorname{Pf}\mathscr B=\frac{\operatorname{Pf}p}{\operatorname{Pf}q},
\qquad
\frac{a}{\overline{\operatorname{Pf}\mathscr B}}=\bar a.
\tag{7}
$$

r/2为偶数，Pf(q)=Pf(−p*)=overline(Pf p)，没有遗失(−1)^(r/2)。反射仅重排局部质量块，aθ=a。

对任意可逆反对称B和上述酉V，互补子式公式写成

$$
\operatorname{Pf}\!\left(u^{\mathsf T}[-(B^*)^{-1}]u\right)
=
\frac{\overline{\operatorname{Pf}(v^{\mathsf T}Bv)}}
{\bar d\,\overline{\operatorname{Pf}B}}.
\tag{8}
$$

可以先在可逆压缩块上通过Pfaffian Schur分解验证，再乘去可逆整矩阵的Pfaffian，以多项式恒等式延伸到奇异压缩块。因而式(8)无需对近零压缩块取逆，方向符号由前u后v的排序固定。

## 4. 原权重的反射共轭关系

反射后辅助Pf变成Aᵥ，质量压缩由式(7)—(8)变成Qᵤ的互补式。代入(6)得

$$
F^\theta
=\frac{A_v\,a}{d}\,
  \frac{\bar Q_v}{\bar d\,\overline{\operatorname{Pf}\mathscr B}}
=A_v\bar a\,\bar Q_v
=\overline{\frac{A_u aQ_v}{d}}
=\bar F.
\tag{9}
$$

这是保持原相位的身份，没有取绝对值，也没有对不同配置独立归一。

固定Wilson谱投影时，原F=Pf Nλ是p实部、虚部的有限多项式；反射关系同样如此。可逆反对称p在允许的反对称矩阵空间中稠密。因此

$$
F^\theta(p,-p^*)=\overline{F(p,-p^*)}
\quad\text{对包括奇异质量、零质量在内的全部 }p.
\tag{10}
$$

这项延伸只用于证明，最终定义仍是673不含逆质量的固定尺寸矩阵。H有零模的背景依673 §3属于两时间边界中的Haar零测度，不影响积分；未由此领取背景可微性或连续局域性。

式(9)—(10)覆盖局部质量来源：允许把p在有限个格点作任意保持q=−p*的实参数变形，反射时同步输送该变形；逐多项式系数给相应质量插入身份。它不自动证明任意单费米源j的同一反射作用。

## 5. 接入完整双边界Gauss平均

使用673的全部原H_b热核和群𝒢=G^V：

$$
\mathcal K(x_i,x_j)=\int da\,db\
k_\tau(q_i,a q_j)k_\tau(q_j,b q_i)
F(x_i,x_j;a,b).
\tag{11}
$$

反射交换xᵢ、xⱼ，各自把两条时间接口变为逆元；原反周期符号保持：

$$
F(x_j,x_i;a^{-1},b^{-1})
=\overline{F(x_i,x_j;a,b)}.
\tag{12}
$$

由于H_b为原实标量自伴算符，kτ实、对称且规范协变；Haar在取逆下不变。改变积分变量(a,b)→(a⁻¹,b⁻¹)，由673的绝对可积性可合法换序，得到

$$
\mathcal K(x_j,x_i)=\overline{\mathcal K(x_i,x_j)}.
\tag{13}
$$

再以原配置与辅助测度对两个半区积分，定义Z；交换半区给

$$
Z=\bar Z\in\mathbb R,\qquad
\partial_\eta Z\in\mathbb R
\quad\text{对反射配对的实局部质量来源 }\eta.
\tag{14}
$$

这关闭673此前留下的标量Hermitian／实性缺口。没有恢复任意Grassmann来源的一般断言，也没有把实数Z当成正数或非零数。

## 6. 计算规则与数值核验

式(12)还允许为同一Haar积分成对计算反射配置：

$$
\frac12\left[
F(x_i,x_j;a,b)+
\overline{F(x_j,x_i;a^{-1},b^{-1})}
\right]=F(x_i,x_j;a,b).
\tag{15}
$$

这可以核对被积式或构造成对估计器；它不允许删除负特征值、对整个矩阵人为取正部，或用若干样本代替完整群积分。

第一组在原64维配对上检查rᵤ=0、30、31、32、34、64，覆盖空块、不同秩和奇数辅助块。它们是代数验证样本，不冒充实际Wilson拓扑配置。非零权重的反射相对误差最大约1.0×10⁻¹⁴，Schur权重误差约2.5×10⁻¹⁴，互补质量子式残差约6.6×10⁻¹⁴。原Higgs取零、只保s时物理配对秩为4，直接无逆质量计算仍满足反射身份，相对误差约5.3×10⁻¹⁵。

第二组回到四份原全群、不同φ及E的实际2×2时间空间配置；不是只用低维子群。Wilson反射残差为机器零，overlap残差<1.6×10⁻¹⁵，质量反射双重关系误差<10⁻¹⁵，质量Pf前因子误差<2.8×10⁻¹⁵。

[首次代码](674/drafts/first_reflection_attempt.py)曾对近零Pfaffian继续使用相对互补误差，检查失败；[诊断](674/drafts/first_reflection_diagnostic.json)保留。与673相同的近奇异样本中，某些未缩放质量Pfaffian因大尺度因子达到10³³—10³⁶，其相对舍入误差仍无意义。最终保持全部样本，使用原矩阵双重身份和规模归一后的绝对多项式余量；相对权重／子式误差只在可分辨样本报告。归一后的很小余量不是数值精度证明，解析结论由式(4)—(10)承担。

[首次通过版本](674/drafts/first_pass_reflection.py)和[结果](674/drafts/first_pass_reflection_results.json)保留；最终版增加空谱块约定与端点检查，并分开记录规范幅度和Pfaffian平衡尺度。没有改动历史稿件，无图像检查。

## 7. 下一项仍是正性和同一过程

$$
\mathcal K=\mathcal K^\dagger,\quad Z\in\mathbb R
\quad\not\Longrightarrow\quad
\sum_{ij}\bar c_i\mathcal K(x_i,x_j)c_j\ge0
\quad\text{或}\quad Z>0.
\tag{16}
$$

673的可积性、674的实性合在一起给一个明确的Hermitian候选；672的固定边界负方向仍不能未经全部群平均就升级为物理反例。尚缺完整正半区代数、严格归一、与643原有序CAR影响／624仪器的身份及共同时间参数。

接[675入口](675/drafts/STATUS.md)：对完整群平均寻找真实正性见证或反例，同时核原H_b及质量／时间参数。可研究原H_b基态主导时对边界核的受控数学约化，但必须声明参数如何保持，不能把候选参数族直接称作原H_F的低温演化。研究目标、四分支边界和旧空间接口保持不变。
