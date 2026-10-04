# 第611轮：原物种质量、旋量投影与共同来源字典

日期：2026-10-01。接[610](../../research_note_610.md)及[611入口](STATUS.md)。[代码](../joint_original_mass_spinor_bridge.py)、[结果](../joint_original_mass_spinor_bridge_results.json)、[核验](../research_round_611_checks.json)、[条件账](../unified_physics_condition_ledger_611.md)。三组复算、十二式；主代理审查，无新增独立代理审查。

## 1. 从任意质量诊断返回原32模式

609给出了正常／配对项的共同补全公式，但其四模式质量诊断不是598原物种。此次直接导入598原32模式、原Higgs／单态标量依赖、全部原复Y系数，以及606各模块的实际电荷投影，不重新拟合参数。

**新增连接：** 在平坦零动量背景和明确旋量字典下，原Dirac与Majorana矩阵、完整64维BdG质量谱及规范协变均可同时保持。直接把两分量自旋恒等矩阵换成四分量恒等矩阵，会把所有Dirac质量错误删掉。采用指定β映射可修复静止接口，但规范场改变后，同一投影一般改变质量块和相应来源；保谱需要额外改变场依赖质量，不能称为免费换记号。

这里“原质量”指598已保存的候选耦合，不是实验拟合的标准模型参数。四维Euclidean纤维尚不等于实时间CAR态空间，本轮同时给出不能逐点识别其共轭场的明确限制。

|层次|地位|
|---|---|
|继承|598原群／表示／32模式及复Y、604原BdG约定、606带电投影、609补全|
|附加输入|四分量旋量、β=γ₄、左右手基及指定时间方向；平坦零动量片|
|解析结果|原质量与规范共同字典；错误嵌入排除；变化投影下质量和来源联合变化|
|数值验证|原32模式全部矩阵与64维BdG；全部模块真实Wilson零动量投影；原标量来源|
|未得到|物种、旋量或时间方向必然性、物理极点质量、完整实时间手征理论、测度或GR|

去重复：604的物种复制与609的任意配对泄漏不重做。本轮补的是原真实数据到当前候选的接口。早期完整过程要求继续约束后续，静止质量相等不是过程等价。

## 2. 原模块和明确旋量字典

沿598物理右手场约定，模块Q、u、d、L、e、ν的内部重数分别6、3、3、2、1、1；每个有两个spin分量。电荷分别1、4、−2、−3、−6、0。定义

$$
\mathfrak f_{\rm old}=\bigoplus_A R_A\otimes\mathbb C^2,\quad
\dim\mathfrak f_{\rm old}=32,\qquad
\mathfrak f_{\rm amb}=\bigoplus_A R_A\otimes\mathbb C^4,\quad
\dim\mathfrak f_{\rm amb}=64 .
\tag{1}
$$

后者是嵌入环境，不是新增一倍物理Weyl物种。允许子空间仍32维；全Fock维数也不能作物种计数。

使用606的γ表示，γ₅=diag(I₂,−I₂)。取

$$
V_L=\begin{pmatrix}0\\I_2\end{pmatrix},\qquad
\beta=\gamma_4=\begin{pmatrix}0&I_2\\I_2&0\end{pmatrix},\qquad
V_R=\beta V_L,\qquad
V_L^\dagger\beta V_R=I_2,\qquad
V_L^\dagger V_R=0 .
\tag{2}
$$

按Q、L取V_L，其余取V_R，构成全等距J₀。这个选择明确了两套spin指标如何对应。内部群不作用spin，所以环境表示R_amb=⊕R_A⊗I₄满足R_amb J₀=J₀R_old。β是选定旋量／时间结构的一部分，不是从FUCP推出。

## 3. 质量、配对、规范与来源一起映射

原Dirac块每个内部矩阵元具有y_ab(φ)I₂，其中y包含原Higgs、复Y以及F^(−1/2)。定义环境正常质量为

$$
(b_\beta)_{ab}=y_{ab}(\phi)\beta
\quad(a\text{左手},\,b\text{右手}),
\qquad (b_\beta)_{ba}=(b_\beta)_{ab}^\dagger .
\tag{3}
$$

这由原y与一份β确定，不是为每个物种另造新参数。对原中性ν_R Majorana项，用ε=iσ₂，环境自旋配对取

$$
\varepsilon_4=\operatorname{diag}(\varepsilon,\varepsilon),\qquad
(\Delta_{\rm amb})_{\nu\nu}
=\frac{Y_s s}{\sqrt F}\varepsilon_4,\qquad
\Delta_{\rm lift}=p_0\Delta_{\rm amb}p_0^T,\quad p_0=J_0J_0^\dagger .
\tag{4}
$$

最后一步正是609所要求删除补模式配对；保原允许项而改变环境补项，不能遗漏。由(2)逐块得到

$$
J_0^\dagger b_\beta J_0=h_{\rm old},\qquad
J_0^\dagger\Delta_{\rm lift}J_0^*=\Delta_{\rm old},\qquad
\mathcal J_0^\dagger
\begin{pmatrix}b_\beta&\Delta_{\rm lift}\\
\Delta_{\rm lift}^\dagger&-b_\beta^T\end{pmatrix}
\mathcal J_0
=
\begin{pmatrix}h_{\rm old}&\Delta_{\rm old}\\
\Delta_{\rm old}^\dagger&-h_{\rm old}^T\end{pmatrix},
\quad\mathcal J_0=\operatorname{diag}(J_0,J_0^*) .
\tag{5}
$$

因此允许CAR质量算符、原BdG矩阵与其谱全部保持；外代数提升使用相同系数，不需构造2⁶⁴全矩阵来“数值证明”一般CAR恒等式。原规范协变同时继承；ε只用于规范中性的ν，Dirac块的Higgs荷仍按598补偿。对φ或几何参数的偏导，因J₀固定也同步映射。

错误的直接I₄嵌入则给

$$
(b_I)_{ab}=y_{ab}I_4
\quad\Longrightarrow\quad
J_0^\dagger b_IJ_0=0,\qquad h_{\rm old}\ne0 .
\tag{6}
$$

这不是原Yukawa不允许，而是对象字典错误。诊断φ=(.4,.7,−.2,.1,.5)时原正常质量范数为.76258669，被错误映射全部删去；β字典的正常、配对及BdG误差均为零。原完整群的随机变换核验最大残差1.69×10⁻¹⁵。

## 4. 变化投影不再保原质量块

沿606实际平坦单方向holonomy族，对模块电荷q_A，

$$
V_A(\theta)=U_A(\theta)V_{\chi_A},\qquad
U_A=e^{iq_A\theta\gamma_\mu/2},\qquad
J(\theta)=\bigoplus_A I_{R_A}\otimes V_A(\theta),\qquad
h_{\rm eff}=J^\dagger b_\beta J .
\tag{7}
$$

这里h_eff是质量部分的压缩，尚非含完整空间动能的物理色散或极点质量。标量背景固定，原y未被重拟合。空间方向μ=1、2、3时β反对易γμ；指定Euclidean时间方向μ=4时二者相同。因此

$$
(h_{\rm eff})_{ab}=f_{ab}(\theta)(h_{\rm old})_{ab},\qquad
f_{ab}(\theta)=
\begin{cases}
\cos((q_a+q_b)\theta/2),&\mu=1,2,3,\\
\cos((q_b-q_a)\theta/2),&\mu=4.
\end{cases}
\tag{8}
$$

推导直接把U_a†βU_b夹在V_L†、V_R之间；含一个额外γμ的项因手征性为零。空间方向四类Dirac块u、d、ν、e分别出现5、−1、−3、−9的电荷组合；时间方向差荷绝对值均为3。原中性ν_R配对的q=0，故在这个族上不变。

代码对所有原模块，从完整64维Wilson核抽取零动量投影，核验(7)，不是另选任意旋转曲线。空间／时间方向的差异依赖当前指定β与投影片；没有由此宣布Lorentz破缺的实验预测或原标准模型质量测量失配。它表明静止接口不能直接推广成完整理论字典。

## 5. 原质量的来源也随投影一起改变

θ是链路相位，不是钟或度规。虽然环境bβ在该诊断中不依赖θ，压缩质量却有

$$
\partial_\theta h_{\rm eff}
=(\partial_\theta J^\dagger)b_\beta J+
J^\dagger b_\beta(\partial_\theta J)
\ne J^\dagger(\partial_\theta b_\beta)J=0 .
\tag{9}
$$

本平坦标架J†∂θJ=0，所以它也是该质量块在此标架中的协变导数；相应力取负号。一般标架须包括连接，不能把普通导数当作规范不变量。θ=.17时导数范数.78607585，直接差分误差3.88×10⁻¹¹。

对原标量φ，J不依赖它，所以

$$
\partial_{\phi_a}h_{\rm eff}
=J^\dagger(\partial_{\phi_a}b_\beta)J,\qquad
(\partial_{\phi_a}h_{\rm eff})_{LR}
=f_{LR}(\theta)(\partial_{\phi_a}h_{\rm old})_{LR}.
\tag{10}
$$

原五个标量方向、相同F与相同Y的差分核验通过，最大残差4.56×10⁻¹¹。不能一边用投影后的质量，一边继续用未经同样映射的旧标量力或链路来源。这把C15质量与C22反作用变成同一个验收条件。

## 6. 强行保谱的代价是什么

在此局部片中可以显式选择随场共转的质量。令U=⊕I⊗U_A，取

$$
b_{\rm co}(\theta)=U(\theta)b_\beta U(\theta)^\dagger,\qquad
\Delta_{\rm co}(\theta)=U(\theta)\Delta_{\rm amb}U(\theta)^T,\qquad
J(\theta)^\dagger b_{\rm co}(\theta)J(\theta)=h_{\rm old}.
\tag{11}
$$

配对限制后同样恢复原允许Δ。数值恢复误差1.69×10⁻¹⁶，但环境正常质量已改变.32166566。这是一个场依赖新质量接法；在其余结构固定的比较中，不是原bβ的免费换基。若主张完整酉等价，须同步转换动能、链路动量、来源、仪器及边界。

更重要的是：本片的U来自全方向平坦holonomy及特定基，尚未构造在一般规范配置上的全局光滑、规范协变、局部字典。不能仅因(11)能形式写出就消除测度、拓扑、局域性或相互作用条件。

## 7. Euclidean共轭场仍不是已构造的Hilbert伴随

[Lüscher原文§2，式(6)—(10)](https://arxiv.org/pdf/hep-lat/9909150)对格点Weyl场使用场依赖的hat P_−，而bar ψ使用固定P_+；两个Grassmann变量在Euclidean积分中分别定义。本项目p_L与hat P_−的符号约定一致，但不能因此直接令bar ψ=ψ†β并获得完整态空间。

在本平坦族，朴素逐点识别要求βp_L(θ)β=P_R(0)。实际对单个电荷q有

$$
\|\beta p_L(\theta)\beta-P_R(0)\|
=\left|\sin\frac{q\theta}{2}\right|,
\tag{12}
$$

由两个旋转秩2投影的主角直接得到，非零θ一般不满足。θ=0时它恢复(2)，正好说明本轮静止字典的边界。这个等式是解析检查，不单独增加数值组数。

它只排除这条朴素逐点识别，不排除带时间／场重定义的重建，也不否认成熟Euclidean手征理论。需要反射正性、转移结构或明确实时间正则化后，才能谈完整过程是否是同一个。598本身已有正Hamiltonian，不代表它自动重建指定的四维GW关联函数。

## 8. 验收、保留条件与下一步

三组检查均通过：

1. 原32模式静止质量、Majorana、64维BdG与完整内部群字典；错误I₄映射及补配对项均明确检出。
2. 全部原模块的真实Wilson平坦投影、空间／时间两类系数和质量块谱变化；未改Y、未改标量背景。
3. 链路／标量来源与场依赖保谱替代；恢复旧质量时确实改变环境质量项。

同一候选现在有了一个**原物种质量的明确静止连接**，以及离开该片必须共同调整的质量、来源和旋量字典条件。没有完成一般背景、相互作用态空间、手征量子测度、物理极点质量或GR。

[612入口](../../612/drafts/STATUS.md)：对接成熟Euclidean到实时间重建，先核自由GW核的反射正性／传播函数及可能所需Hilbert扩展，区分静止质量相等与完整过程相同。复用现有有限谱／过程反例，不重做抽象嵌入论证；若成熟定理只能覆盖自由或向量样分支，范围明确保留。条件整合优先，认知装置设计后置。
