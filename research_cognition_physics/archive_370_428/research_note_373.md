# 第373轮：量子方向怎样与实际传播方向校准

## 0. 本轮问题、基线与结论

**科学基线：已冻结第372轮。** 本轮只处理一个新接口：在明确给定的三维 Weyl 传播模型中，能否通过位置与时间读数，把二能级方向的 Born 概率角校准为传播角？不再重复检查三方向的代数来源，也不从字段数、Bloch 球或矩阵锥单独推断物理维数。

**结论：可以有条件地闭合这个接口。** 六个已校准的内部方向态，用同一包络制备并测量初始质心速度，可以恢复传播系数与漂移。它们确定的空间度量使概率角等于相对传播速度的角。若要进一步把这种初速称为持续传播的射线，需要另行制备上谱支窄带波包；本轮给出质心与完整波函数的误差界。任意固定内部纯态没有这种保证。

已有工作与新增内容的分界：

| 已有结果 | 本轮沿用而不重做 | 本轮新接口 |
|---|---|---|
| [348：量子行走与给定 Dirac 极限](../archive_342_369/research_note_348.md) | 局域酉步与低动量传播的条件接口 | 给定三维生成元上的方向—质心校准 |
| [349：给定变系数传播与有效几何](../archive_342_369/research_note_349.md) | 概率连接项、背景输入与引力动力学的区别 | 常系数三维中，实际流与内部方向的可辨识关系 |
| [370：雷达字段与测量秩](research_note_370.md) | 字段名称不决定维数，端口需有可辨识结构 | 内部端口如何经传播读数得到几何意义 |
| [371：方向比特及条件角度恢复](research_note_371.md) | 内部概率角、额外方向任务假设 | 校准概率角与物理位置读数之间的映射 |
| [372：二能级正锥与 Lorentz 锥](research_note_372.md) | 正矩阵的锥同构本身不是事件字典 | 从指定演化推出守恒概率流，再识别该流的锥 |

第183、210—212轮的高维反例与三生成元问题，以及第282、289—290、343轮的图度数、时钟、环路和维数非唯一性不作为本轮新增实验。特别地，本轮**没有从认知原则选出三维位置空间或这个自然哈密顿量**。

## 1. 模型输入与原始文献接口

额外输入包括：位置流形为三维欧氏坐标空间、可操作的时间与位置测量、二分量波函数、固定实可逆矩阵 E、固定漂移向量 β，以及下述生成元。取 ħ=1；E 与 β 的单位均为速度。首先在足够光滑、衰减充分的态上计算，再由紧支撑动量及适当位置矩条件扩展。

$$
H(\mathbf k)=\boldsymbol\beta\cdot\mathbf k\,I
+(E\mathbf k)\cdot\boldsymbol\sigma,
\qquad
V_i=\frac{\partial H}{\partial k_i}
=\beta_i I+\sum_a E_{ai}\sigma_a,
\qquad
i\partial_t\Psi=-i\sum_iV_i\partial_i\Psi .
\tag{1}
$$

线性二能级 Weyl 哈密顿量及其倾斜、各向异性系数是已有结构，不是本轮发现。[Soluyanov 等，2015，式(1)—(2)](https://arxiv.org/html/1507.01603) 给出这一类谱；[Nissinen–Volovik，2017，式(1)—(7)](https://arxiv.org/html/1702.04624) 讨论其 tetrad 与有效度量接口。本轮只使用后者的普通线性哈密顿量子类，不引入含频率自能的额外分支。

量子行走对接也有明确前提：[D’Ariano–Perinotti，2014，§II、§VII](https://arxiv.org/html/1306.1934) 在局域性、均匀性、各向同性等输入下讨论 Weyl/Dirac 极限；其波矢、内部极化方向、群速度是不同对象。本轮补足项目中的校准步骤，不声称重新推导该分类。

记 G 为正定空间逆度量矩阵，g 为其逆：

$$
G=E^{\mathsf T}E,\qquad g=G^{-1},\qquad
\omega_\pm(\mathbf k)=\boldsymbol\beta\cdot\mathbf k\pm|E\mathbf k|,
\qquad
\det(\omega I-H)=(\omega-\boldsymbol\beta\cdot\mathbf k)^2
-\mathbf k^{\mathsf T}G\mathbf k .
\tag{2}
$$

一般只能称“上谱支”，不能无条件称“正能量支”：若倾斜足够大，上支也可能在部分波矢处为负。数值例另外满足 βᵀgβ<1，因此其上支在非零波矢处确为正。谱公式只给出指定模型的特征结构，不给出 Einstein 场方程。

## 2. 从演化推出概率流，而不是把内部态直接命名为位置

对式(1)与其伴随逐项相乘，由各 V_i 为常 Hermitian 矩阵得到：

$$
\rho=\Psi^\dagger\Psi,\qquad
s_a=\Psi^\dagger\sigma_a\Psi,\qquad
\mathbf j=\boldsymbol\beta\rho+E^{\mathsf T}\mathbf s,\qquad
\partial_t\rho+\boldsymbol\nabla\cdot\mathbf j=0 .
\tag{3}
$$

这里 j 是给定演化的实际概率流。对混合态的局部二分量密度矩阵 D，同一关系逐项线性成立。其正性给出：

$$
D=\frac12(\rho I+\mathbf s\cdot\boldsymbol\sigma)\succeq0
\ \Longrightarrow\ 
|\mathbf s|\le\rho,\qquad
(\mathbf j-\boldsymbol\beta\rho)^{\mathsf T}g
(\mathbf j-\boldsymbol\beta\rho)=|\mathbf s|^2\le\rho^2,
\qquad
\rho^2-|\mathbf s|^2=4\det D .
\tag{4}
$$

因此，第372轮的矩阵锥在本模型中有了**概率流**这一操作含义。纯局部旋量对应锥面，局部混合态对应锥内。这个证明不把一个内部密度矩阵直接等同于一个时空事件；也不意味着波包质心必在锥面上。

式(2)—(4)对应一个有效线元代表，采用负、正、正、正号约定：

$$
ds^2=-dt^2+
(d\mathbf x-\boldsymbol\beta\,dt)^{\mathsf T}
g(d\mathbf x-\boldsymbol\beta\,dt).
\tag{5}
$$

常系数情形是平直有效几何。其钟尺单位、整体尺度和坐标解释仍依赖已经给定的操作。这里既没有几何反作用，也没有证明其他物质种类共享这个度量。

## 3. 六个内部方向怎样校准传播系数

准备同一归一化空间包络 f 与已知纯内部态 χ_n 的乘积。内部方向 n 是单位 Bloch 向量，分别取三个内部坐标轴的正、负方向。在全空间且位置一阶矩存在时，由式(3)积分或位置算符的定义得到：

$$
\Psi_{\mathbf n}(\mathbf x,0)=f(\mathbf x)\chi_{\mathbf n},
\qquad
\chi_{\mathbf n}^\dagger\boldsymbol\sigma\chi_{\mathbf n}=\mathbf n,
\qquad
\left.\frac{d\langle\mathbf x\rangle_{\mathbf n}}{dt}\right|_{0}
=\boldsymbol\beta+E^{\mathsf T}\mathbf n .
\tag{6}
$$

令这六个初速为 v_(a,+)、v_(a,−)，则无需预先指定 E 的数值，就能从测量值恢复：

$$
\boldsymbol\beta=\frac16\sum_{a=1}^{3}
(\mathbf v_{a,+}+\mathbf v_{a,-}),
\qquad
E_{ai}=\frac{v_{a,+,i}-v_{a,-,i}}2,
\qquad
g=(E^{\mathsf T}E)^{-1}.
\tag{7}
$$

这里“测量”指独立重复制备的总体均值，不是对一个未知量子态进行六种无扰动读取。需要已知内部方向设置、共同包络制备、位置测量与时间读数。这些资源是本轮的额外输入。六个初态一般含上下两支，其初速可以用于校准，却不能因此说每个初态以后都沿此方向直行。

### 3.1 只用正向时间读数的误差

实际不需要负时间演化。用 0、h、2h 三个时刻的质心构成二阶前向差：

$$
\widehat{\mathbf v}_{h}
=\frac{-3\mathbf X(0)+4\mathbf X(h)-\mathbf X(2h)}{2h},
\qquad h>0.
\tag{8}
$$

对每个分量，误差是 X_i 的三阶导数与下述 Peano 核的积分。核恒非正：

$$
K_h(u)=
\begin{cases}
-u+\dfrac{3u^2}{4h},&0\le u\le h,\\[4pt]
-\dfrac{(2h-u)^2}{4h},&h\le u\le2h,
\end{cases}
\qquad
\int_0^{2h}|K_h(u)|\,du=\frac{h^2}{3}.
\tag{9}
$$

若动量包络支撑在 |k−k₀|≤ε，记 M=|Ek₀|+‖E‖₂ε。支撑不随这个常系数演化变化；标量漂移与 V_i 对易。由速度的二阶导数及双对易子估计：

$$
X'''_i(t)=-\langle[H,[H,V_i]]\rangle_t,\qquad
|X'''_i(t)|\le4M^2\|E_{\cdot i}\|_2,
\qquad
\|\widehat{\mathbf v}_h-\dot{\mathbf X}(0)\|_2
\le b_h:=\frac43h^2M^2\|E\|_{\mathrm F}.
\tag{10}
$$

紧动量支撑与有限位置矩保证这里所用导数存在。把每个速度的误差代入式(7)可得：

$$
\|\widehat{\boldsymbol\beta}-\boldsymbol\beta\|_2\le b_h,\qquad
\|\widehat E-E\|_2\le\sqrt3\,b_h,\qquad
\|\widehat G-G\|_2
\le2\sqrt3\,\|E\|_2b_h+3b_h^2 .
\tag{11}
$$

式(10)—(11)只控制有限演化时间误差，不包含有限样本、位置计量、积分与数值微分误差。若每次质心读数还有范数不超过 η 的误差，式(8)会额外放大至不超过 4η/h；所以把 h 取到零并非无成本的实验操作。数值程序中的正负动量差分仅用于积分验证，不代表物理时间倒流。

## 4. 概率角与传播角：度量必须经过校准

对单位方向 n、m，令 P_n、P_m 为纯态投影，w_n 为减去漂移后的初速。式(7)恢复的 g 给出：

$$
P_{\mathbf n}=\frac12(I+\mathbf n\cdot\boldsymbol\sigma),\qquad
\mathbf w_{\mathbf n}=E^{\mathsf T}\mathbf n,\qquad
\|\mathbf w_{\mathbf n}\|_g=1,\qquad
\cos\theta_g
=\mathbf w_{\mathbf n}^{\mathsf T}g\mathbf w_{\mathbf m}
=\mathbf n\cdot\mathbf m
=2\operatorname{tr}(P_{\mathbf n}P_{\mathbf m})-1 .
\tag{12}
$$

这回答了本轮的核心问题：在指定传播规律成立且校准可执行时，内部概率角与相对传播角有一条可检验的桥梁。普通坐标欧氏角并不自动满足这个等式。对所有方向都保持坐标角的充要条件是 E Eᵀ 为单位阵的正倍数：必要性可先用任意正交向量对得到 E Eᵀ 的非对角项为零，再用两个坐标轴的和、差得到所有特征值相等。

一个不能靠坐标轴掩盖的反例是：

$$
E=\operatorname{diag}(1,2,3),\qquad
\mathbf n=\frac{(1,1,0)}{\sqrt2},\qquad
\mathbf m=\frac{(1,-1,0)}{\sqrt2},\qquad
\cos\theta_{\mathrm{Born}}=0,\quad
\cos\theta_{\mathrm{coordinate}}=-\frac35,\quad
\cos\theta_g=0.
\tag{13}
$$

即使 E 各向同性，也必须先减去 β，不能把含共同漂移的原始速度角直接称为内部概率角。

内部参考轴的共同旋转不改变这个度量。若内部设置和 E 同步变换：

$$
R\in SO(3),\qquad E'=RE,\quad \mathbf n'=R\mathbf n,
\qquad
E'^{\mathsf T}\mathbf n'=E^{\mathsf T}\mathbf n,\qquad
E'^{\mathsf T}E'=E^{\mathsf T}E .
\tag{14}
$$

因此内部绝对轴标签不影响同一实验室的校准结果。式(14)没有解决不同主体之间的位置、时钟、姿态或平移如何对齐；那些仍须通过真实传播数据检验。

## 5. 从初速到持续传播：必须额外选择谱支

固定内部旋量与同一包络的乘积态，一般在各波矢处同时占据上下两支。自旋会随时间转动，流也会变化。本轮数值例中，初始 +X 方向的总流在 t=0.7 时改变约 1.40335；因此不能用式(6)把所有纯内部态解释成永远直行的射线。

另行制备上谱支包络：

$$
\widehat\Psi_0(\mathbf k)=f(\mathbf k)u_+(\mathbf k),
\qquad
(E\mathbf k)\cdot\boldsymbol\sigma\,u_+
=|E\mathbf k|u_+,\qquad
\widehat\Psi_t(\mathbf k)=e^{-it\omega_+(\mathbf k)}
\widehat\Psi_0(\mathbf k).
\tag{15}
$$

这是额外的带选择制备条件，不是无成本地把任意未知输入投影后宣布其能力保留。取包络支撑于半径 ε 的动量球，且远离交叉点：

$$
m_0=|E\mathbf k_0|,\qquad
\mu=m_0-\|E\|_2\epsilon>0,\qquad
\mathbf n(\mathbf k)=\frac{E\mathbf k}{|E\mathbf k|},\qquad
L=\frac{\|E\|_2^2}{\mu}.
\tag{16}
$$

支撑凸包内的谱隙至少为 2μ。上支群速度与 Hessian 为：

$$
\mathbf v(\mathbf k)=\boldsymbol\nabla_{\!k}\omega_+
=\boldsymbol\beta+E^{\mathsf T}\mathbf n(\mathbf k),\qquad
\nabla_{\!k}^2\omega_+
=\frac{E^{\mathsf T}(I-\mathbf n\mathbf n^{\mathsf T})E}
{|E\mathbf k|},\qquad
\|\nabla_{\!k}^2\omega_+\|_2\le L .
\tag{17}
$$

令 v₀=v(k₀)。对该支撑上的所有归一化上谱支态，均有：

$$
\sup_{\mathrm{supp}f}\|\mathbf v(\mathbf k)-\mathbf v_0\|_2
\le L\epsilon,\qquad
\|\overline{\mathbf v}-\mathbf v_0\|_2\le L\epsilon,\qquad
\overline{\mathbf v}
=\int |f(\mathbf k)|^2\mathbf v(\mathbf k)\,d^3k.
\tag{18}
$$

总概率流恒等于该平均速度；只要位置一阶矩存在，质心满足精确关系：

$$
\langle\mathbf x\rangle_t-\langle\mathbf x\rangle_0
=t\,\overline{\mathbf v},\qquad
\|\langle\mathbf x\rangle_t-\langle\mathbf x\rangle_0-t\mathbf v_0\|_2
\le |t|L\epsilon .
\tag{19}
$$

这仍不是波包不形变。若以**完整初态**平移作比较，利用式(17)的相位 Taylor 余项和两个单位相位的差估计，可得到更强的受控陈述：

$$
\left\|\Psi_t(\mathbf x)-
e^{-it(\omega_+(\mathbf k_0)-\mathbf v_0\cdot\mathbf k_0)}
\Psi_0(\mathbf x-t\mathbf v_0)\right\|_{L^2}
\le\min\left\{2,\frac{|t|L\epsilon^2}{2}\right\}.
\tag{20}
$$

比较项保留 u_+(k) 的全部动量依赖，不把它替换成中心波矢处的常旋量；后者需要另加极化误差。本界对指定谱支与支撑的态族一致，不是对整个 Hilbert 空间的无条件算符范数收敛。

有限宽度包络的内部约化态一般是混合态。其平均相对流可严格在式(4)的锥内；局部纯旋量流在锥面，与平均质心不在锥面并不矛盾。紧动量支撑还意味着不能把这些波包误称为实空间紧支撑的点事件。

## 6. 数值实现与资源代价

### 6.1 全空间位置积分，避免周期边界伪影

数值采用下列模型参数：

$$
E=
\begin{pmatrix}
1&0.25&0\\
0&1.6&0.1\\
0.15&0&0.8
\end{pmatrix},
\qquad
\boldsymbol\beta=(0.15,-0.08,0.05)^{\mathsf T},
\qquad
\mathbf k_0=(1.2,-0.7,0.9)^{\mathsf T}.
\tag{21}
$$

这里 βᵀgβ≈0.0296566。采用紧支撑、归一化的实动量包络：

$$
q_i=\frac{k_i-k_{0i}}{\delta},\qquad
f_\delta(\mathbf k)=
\left(\frac{315}{256}\right)^{3/2}\delta^{-3/2}
\prod_{i=1}^3(1-q_i^2)^2\,\mathbf1_{\{|q_i|<1\ \forall i\}},
\qquad
\epsilon=\sqrt3\,\delta .
\tag{22}
$$

归一化来自一维积分 ∫_(−1)^1(1−q²)^4dq=256/315。包络在支撑边界为 C¹，具有所需位置矩；所有算例的本征旋量都处在同一光滑局部图中。程序对实际矩阵指数演化后的完整旋量进行 Gauss–Legendre 张量求积，而不通过预先写入群速度来制造质心位移：

$$
\langle x_i\rangle
=\operatorname{Re}\int\widehat\Psi^\dagger
i\partial_{k_i}\widehat\Psi\,d^3k,\qquad
\langle x_i^2\rangle=\int|\partial_{k_i}\widehat\Psi|^2\,d^3k .
\tag{23}
$$

式(23)是全空间 Fourier 变换下的位置算符，不是有限周期盒中锯齿坐标的对易子。积分没有周期绕回问题。程序用独立动量差分计算导数，并以直接积分总流作对照；还检查求积阶数与动量差分步长的细化。数值矩阵指数在另一检查中与独立本征分解比较。

### 6.2 收窄动量的空间代价

对常旋量，式(22)的每个坐标位置方差精确为 3/δ²。对归一化且光滑的 u_+(k)，实包络导数与旋量导数的实交叉项为零；剩余旋量项减去 Berry 均值的平方非负，由 Cauchy–Schwarz 即可证明：

$$
\sum_i\operatorname{Var}(x_i)\ge\frac{9}{\delta^2}.
\tag{24}
$$

因此，把群速度方向逼近得更准确需要更宽的空间包络。完整资源账还包括谱支选择、方向与包络准备、重复样本、位置与时间计量；本轮不把这些控制资源说成已经从认知公理免费产生。

### 6.3 主要结果

六方向实验使用 δ=0.18。所有时刻均为 0、h、2h：

| h | 恢复 E 的算符范数误差 | 式(11)的 E 上界 | 恢复 g 的算符范数误差 |
|---:|---:|---:|---:|
| 0.004 | 9.66354×10⁻⁵ | 3.77152×10⁻⁴ | 1.98207×10⁻⁴ |
| 0.001 | 6.03306×10⁻⁶ | 2.35720×10⁻⁵ | 1.23759×10⁻⁵ |

h=0.002 也保存在完整结果中；两次减半均符合本例的二阶时间误差。β 的误差约 8.2×10⁻¹⁰，较小是本例对称正负内部态的额外抵消，不代替一般上界。

上谱支波包在 t=0.7 的结果如下。其中心群速度约为 (0.828661,−0.894247,0.410978)：

| δ | 平均流相对中心群速度误差 | 完整态刚性平移误差 | 初始位置方差之和 | 理论位置方差下界 |
|---:|---:|---:|---:|---:|
| 0.24 | 4.03147×10⁻³ | 4.73292×10⁻³ | 156.51096 | 156.25 |
| 0.12 | 1.01137×10⁻³ | 1.18406×10⁻³ | 625.26040 | 625 |
| 0.06 | 2.53061×10⁻⁴ | 2.96067×10⁻⁴ | 2500.26004 | 2500 |

这些相位误差均在式(20)界内。本例对称包络的平均速度误差观察到约二阶缩小；普遍证明采用式(18)的保守一阶界，不能由三行数据把后者升级为所有制备的二阶定理。

独立位置积分给出的位移除以时间，与直接总流之差分别约 3.36×10⁻⁹、1.34×10⁻⁸、5.36×10⁻⁸。固定动量差分步长对更窄包络的误差会变大；减半该步长能改善误差，不能把这点掩盖成严格零。δ=0.12 时求积阶数由18升至24，质心速度差约 1.16×10⁻¹²。概率守恒达到双精度舍入范围。

## 7. 可复算代码、检查与限制

本轮文件：

- [程序：propagation_direction_calibration_audit.py](373/propagation_direction_calibration_audit.py)
- [结果：propagation_direction_calibration_audit_results.json](373/propagation_direction_calibration_audit_results.json)
- [单轮核验：research_round_373_checks.json](373/research_round_373_checks.json)

复算命令，使用项目既有 Python 与 NumPy：

~~~powershell
& 'C:/Users/admin/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe' -B -X utf8 'research_cognition_physics/archive_231_/propagation_direction_calibration_audit.py'
~~~

先执行检查，再以同一命令追加 --write-results 保存结果。实际执行 **15项检查全部通过**，Python 3.12.14、NumPy 2.3.5。本轮共 **24个连续编号公式**。

| 检查 | 所验证内容 |
|---|---|
| 1 | Hermitian 谱、特征行列式与独立本征分解的实际传播子 |
| 2 | 三个平面波叠加的独立时空差分连续性方程 |
| 3 | 纯、混合局部密度的概率流锥与行列式关系 |
| 4—5 | 六方向真实质心校准、前向时间二阶收敛和解析界 |
| 6—8 | 内部参考旋转不变性、概率角与 g 角、非轴向各向异性反例 |
| 9—10 | 谱梯度、带保持、范数与上支总流守恒 |
| 11 | 任意固定纯内部方向不是持续射线的反例 |
| 12—13 | 窄带速度、完整态平移误差与独立位置积分 |
| 14—15 | 积分与差分细化、窄带的空间方差代价 |

**已解析证明：** 对指定常系数模型和相应定义域，式(3)—(24)中的守恒、校准、角关系、有限时间误差及窄带误差。数值只验证明确有限参数，不代替这些证明。

**没有证明：** FUCP 必须选择式(1)；物理维数必须为3；位置、平移与时钟从内部态自动出现；所有物质共享同一 E；常系数必须推广成哪个变系数方程；几何必须服从 Einstein 动力学；内部概率锥已经是现实宇宙全部事件的因果锥。

## 8. 物理解释与真正下一步

本轮把“内部量子方向”与“传播几何”之间的语言类比，推进为一项有前提、有校准操作、有误差及资源账的数学对应。给定传播规律可以提供这一字典；认知原则为什么选择该传播规律，仍是不同问题。

最直接的后继判据是：**两个不同的可操作传播通道，分别经过本轮校准，是否必须得到同一个漂移和空间度量？** 同一个 Bloch 球、同样 Born 规则或两个通道能够通信，本身都没有在本轮被证明足以强制这个一致性。应先回顾现有多物质共同度规研究，只有找到真正未覆盖的操作条件才开新轮，避免把任意两个 E 不相同的旧式反例再编号。

即使以后建立共同度量，仍须另行证明尺度极限、动态几何及其场方程。本轮既未完成总目标，也未证明总目标不可能；它完成的是“概率角如何接受实际传播读数检验”这一有限而必要的接口。
