# 同一弱衰变的总角能流与渐近引力记忆

2026-10-09。**成熟GR接口采用，计0，不占新轮。** 承接[准入审计](../_admission/after1058_memory_source/selection.md)和[1056](../research_note_1056.md)，将其同一电子记录—总角能流关系接到一个明确的渐近位移记忆系数。本项没有新的拟合、实验数据或认知选择定理；尚待独立审阅。

## 1. 采用哪一条共同接口

物质部分严格保留1056的对象：静止μ⁺，质量m；领先标准V−A作用；条件于 $e^+\nu_e\bar\nu_\mu$ 衰变；忽略三个末态质量及辐射修正；同一相位空间效应 $F(d\alpha)$ 与完整电子标签。没有重新拟合弱流或源系数。

另采用：3+1渐近平直GR、给定Bondi框架、无额外入射辐射记忆的迟滞分支，以及领先G的远区响应。成熟约束要求**总物质能流和初末质量aspect同时进入记忆**。本项选择的零质量末态支中，后一个量的l=2部分为零；这使1056已有总源矩可以直接运输。该边界并不适用于任意散射。

这个采用把M4的一项来源报告接到渐近几何响应；不把1056的完整daughter相干态、引力探测器或有限距离误差列为已完成。

**源条件化合同：** 下文先给每个已分辨硬末态动量配置α赋予领先迟滞记忆响应系数，再用1056的同一算符值测度加权。采用的是这种分支source—memory字典，不是仅以无条件 $\langle T_{\mu\nu}\rangle$ 驱动的一条确定性半经典几何。后者在本例只有各向同性平均源，不能单凭该平均场推出电子标签—记忆相关。1056概率效应本身没有构造这些分支与量子引力场的实际联合仪器，本项也不补造该仪器。

## 2. Bondi约定与完整守恒式

单位 $c=\hbar=1$。本节用度规 $(-+++)$，与1056粒子振幅的 $(+---)$ 书写习惯不同；正能量与三维方向不变，不重新改变其μ⁺极化号。Bondi展开定义

$$
ds^2=-du^2-2\,du\,dr+r^2q_{AB}d\theta^Ad\theta^B
+\frac{2m_B}{r}du^2+rC_{AB}d\theta^Ad\theta^B+\cdots,
\quad q^{AB}C_{AB}=0,\quad N_{AB}=\partial_uC_{AB}.
\tag{1}
$$

q是单位球度规，D是其协变导数。$m_B$有长度量纲；静止质量m的领先质量aspect是Gm。$C_{AB}$也有长度量纲。物理应变为 $C_{AB}/r$，没有用正则场 $h_{\mu\nu}=(g_{\mu\nu}-\eta_{\mu\nu})/\sqrt{32\pi G}$ 代替它。

在这一约定下，约束为

$$
\partial_um_B
=\frac14D_AD_BN^{AB}
-\frac18N_{AB}N^{AB}
-4\pi G\,\mathcal T,\qquad
\mathcal T=\lim_{r\to\infty}r^2T^{\rm matter}_{uu}.
\tag{2}
$$

定义

$$
\mathcal F_{\rm tot}(\Omega)
=\int du\left(\mathcal T+
\frac{N_{AB}N^{AB}}{32\pi G}\right),
$$

则

$$
D_AD_B\Delta C^{AB}
=4\Delta m_B+16\pi G\mathcal F_{\rm tot}.
\tag{3}
$$

这是同一渐近守恒式的两种写法，不允许将通量和末端质量变化分别任意拟合。[Strominger—Zhiboedov，式(2.2)、(3.6)](https://arxiv.org/pdf/1411.5745)

这里保留引力波能源项是为了说明完整账。下文只取领先G：由所选物质过程激发的C、N为O(G)，其能源对式(3)的贡献从O(G²)起。此为采用的微扰阶数；没有给某个有限G、有限窗的统一余项界。独立入射经典引力波不能靠这套计数免费忽略。

## 3. 平滑角任务，避免把单方向源等同完整度规

对任意光滑实f定义无迹Hessian

$$
H_{AB}[f]=D_AD_Bf-\tfrac12q_{AB}D^2f,\qquad
\mathfrak M_f=\int d\Omega\,H_{AB}[f]\Delta C^{AB}.
\tag{4}
$$

在无边界球面两次分部积分，式(3)给

$$
\boxed{\mathfrak M_f
=4\int f\,\Delta m_B\,d\Omega
+16\pi G\int f\,\mathcal F_{\rm tot}\,d\Omega.}
\tag{5}
$$

l=0、1的f有 $H_{AB}[f]=0$，对应总能量／动量平衡，不是记忆中可额外保留的四个角模式。电型l≥2才贡献本任务；(4)不测磁型记忆。

取

$$
f(\mathbf n)=n_xn_z,\qquad D^2f=-6f,\qquad
\int f^2d\Omega=\frac{4\pi}{15},\qquad
\int H_{AB}[f]H^{AB}[f]d\Omega=\frac{16\pi}{5}.
\tag{6}
$$

最后一式由球面Hessian恒等式
$\int|H[f]|^2=\frac12\ell(\ell+1)[\ell(\ell+1)-2]\int f^2$
在l=2得到。相同l≥2结构也见[Bieri—Garfinkle，式(64)—(73)](https://arxiv.org/pdf/1312.6871)，但其由Weyl积分定义的记忆张量与本C有不同符号／归一，不能直接混用。

$\mathfrak M_f$是一个平滑全角加权的渐近剪切变化，长度量纲；它不是单点应变，也不等于一台有限角覆盖探测器的原始输出。

## 4. 为什么1056的这个支没有额外l=2 Coulomb项

对每个1056硬末态事件 $\alpha=(e,k,l)$，同一相位空间逐事件给

$$
\sum_a E_a=m,\qquad \sum_aE_a\mathbf n_a=0,\qquad
\mathcal F_{\rm mat}(\Omega;\alpha)=
\sum_{a=e,k,l}E_a\delta^2(\Omega-\mathbf n_a).
\tag{7}
$$

三者都计入；不能只保未探测中微子或只保带电粒子。

现**另声明**该硬事件的领先迟滞引力边界：初态为静止母体的Coulomb场，末态没有有质量残余或外加支撑源，全部能源进入三个无质量出射体。因此

$$
m_B^{\rm early}=Gm,\quad m_B^{\rm late}=0
\quad\Longrightarrow\quad
\int n_xn_z\,\Delta m_B\,d\Omega=0.
\tag{8}
$$

母体自旋影响角动量部门，不使这一静止领先质量aspect产生l=2；这里没有计算亚领先软项或自旋记忆。逐事件总动量为0也使(7)的l=1部分消失，与(3)的零模约束一致。

式(8)不是1056概率数据单独证明的完整量子几何初态；它是本次明示采用的迟滞边界。没有独立自由引力记忆、不同前后Bondi参考或另一有质量残余。若改变这些条件，须返回式(5)，不能继续把质量aspect项设为零。

特别是，真实电子和中微子非零质量会把末态送到类时无穷远；极小质量与无限远／长时极限不能无条件交换。恢复它们时，必须共同更新衰变核和有质量体的普通记忆，不能仅把null通量公式换几个质量数值。本项不认证这一误差。[Tolish等，§II—III](https://arxiv.org/pdf/1405.6396)

在所声明无质量、领先G支中，由(5)—(8)得

$$
\boxed{\mathfrak M_{xz}(\alpha)=16\pi GQ_{xz}(\alpha),\qquad
Q_{xz}=\sum_aE_an_{a,x}n_{a,z}.}
\tag{9}
$$

这是所采用引力响应的领先系数关系。

### 独立归一核对

物理剪切的迟滞null射线核可写为

$$
\Delta C_{ij}^{\rm TT}(\mathbf N)
=4G\sum_aE_a
\frac{\big[n_{a,i}n_{a,j}\big]^{\rm TT}_{\mathbf N}}
{1-\mathbf N\cdot\mathbf n_a}.
\tag{10}
$$

初态静止体的空间TT项为0。以一条沿z的射线作线性核诊断，
$C_{\theta\theta}=2GE(1+\cos\theta)$，
$C_{\phi\phi}^{\rm orth}=-C_{\theta\theta}$。
对 $f_z=\cos^2\theta-1/3$ 积分，

$$
\int H[f_z]:\Delta C\,d\Omega
=8\pi GE\int_{-1}^1(1+u)(1-u^2)\,du
=\frac{32\pi GE}{3}
=16\pi GE f_z(\hat z).
\tag{11}
$$

这是格林核系数检查，单条射线并未被称作一份完整守恒衰变。式(10)可由[Strominger—Zhiboedov，§5及附录B](https://arxiv.org/pdf/1411.5745)的软／迟滞结果转换；其附B位移张量含 $2G/r$，物理度规应变含 $4G/r$，两者相差1/2。

实际惯性测试体在适用远区极限有
$\Delta L/L=(2r)^{-1}e^Ae^B\Delta C_{AB}$。
这里没有将此渐近式升级成有限r、任意长观测时的误差定理；也不领取固定角BMS观察者与惯性观察者的时钟记忆相同，原文§4.2明确区分二者。

## 5. 同一电子标签与记忆系数，保完整补集

复用[1056 proof B、E—H](../1056/proof.md)：

$$
b=\{2/5\le2E_e/m\le1/2,\ n_{e,z}>0\},\qquad
E_b=\frac{851}{20000}I-\frac{113}{120000}\sigma_z,\quad
E_{\bar b}=I-E_b,
$$


$$
W_b=\int_bQ_{xz}F(d\alpha)=C_b\sigma_x,\qquad
C_b<-m/4500,\qquad C_b/m\simeq-0.001138097578044878.
\tag{12}
$$

采用(9)后，同一硬记录加权的记忆响应系数为

$$
\boxed{\mathfrak W_b:=\int_b\mathfrak M_{xz}(\alpha)F(d\alpha)
=16\pi GC_b\sigma_x.}
\tag{13}
$$

对 $\rho_\pm=(I\pm\sigma_x)/2$，两份完整电子二标签分布相同，而

$$
\mathbb E_\pm[\mathbf1_b\mathfrak M_{xz}]
=\pm16\pi GC_b.
\tag{14}
$$

这是未归一分支矩，不抹去稀少事件。条件于b时才除以原 $p_b=851/20000$。如果只报告C在(6)这个张量方向的投影系数
$a_{xz}=\mathfrak M_{xz}/(16\pi/5)$，相应加权效应是 $5GC_b\sigma_x$；该5与(13)的16π不同，二者不能混用。

这里的b首先是同一渐近硬末态的电子动量标签。若现实装置在有限距离吸收电子再写下记录，电子与吸收器后续的总能源／动量轨迹会改变；不得仍把原三粒子自由逃逸到无穷远当作该装置的完整源。本项采用分辨硬通道的加权响应系数，不宣称已实现一个既写记录又不改变全部未来来源的探测过程。

这里没有新增恢复禁阻：仅电子标签不能重构(13)，只是1056源矩不足结果经同一线性字典运输。不能再次计一组科学。

### 无条件抵消与参考量词

同一1056单粒子谱给

$$
\frac{d\langle\mathcal E_e\rangle}{d\Omega}
=\frac{d\langle\mathcal E_{\bar\nu_\mu}\rangle}{d\Omega}
=\frac m{80\pi}(7+3\mathbf s\cdot\mathbf n),\quad
\frac{d\langle\mathcal E_{\nu_e}\rangle}{d\Omega}
=\frac m{80\pi}(6-6\mathbf s\cdot\mathbf n).
\tag{15}
$$

总和是 $m/(4\pi)$，所以无条件平均的本l=2记忆矩为0，且
$\mathfrak W_{\bar b}=-\mathfrak W_b$。差异只在保留真实标签的相关中，不能只看“成功”而删补集。

对任意母自旋／被动R输入，硬记录加权的R算符为

$$
16\pi G\,\operatorname{Tr}_\mu[(W_b\otimes I_R)\rho_{\mu R}].
\tag{16}
$$

这是线性第一矩字典，沿用1056的经典记录／R边缘。它不是正概率态，也不确定R与量子引力场的完整联合后态。更没有从振幅模方补造未知daughter的非对角相干。

## 6. 采用的层级与不能领取的结论

|已采用并接上的内容|仍未认证|
|---|---|
|同一弱衰变、完整三体硬源与电子标签|完整正常入射包、制备装置、实际量子出射相干态|
|同一Bondi守恒式中的物质通量／末端项／记忆|任意几何边界或独立自由引力初态下仍有(9)|
|无质量支的平滑l=2渐近响应系数|非零末态质量修正的统一误差、有限G及高阶辐射预算|
|领先迟滞响应和原源矩的同一线性运输|实际引力读出POVM、噪声、探测器后态及总装置来源|
|经典标签与被动R加权矩|完整量子策略范数、任意后续反馈或未来历史|
|成熟记忆／软关系|精确不稳定μ渐近S矩阵、IR完成或新的认知定律|

软定理中的无穷远、零频和波区极限不能在固定r上任意交换；原文§5脚注明确其波区条件。μ本来是不稳定粒子，本项使用已声明衰变分支的硬运动学及迟滞源采用，没有认证一个把μ当稳定渐近入态的精确S矩阵。也没有把角概率分布直接同义替换成相干局域应力。

## 7. 复算与交付

[轻量检查](gravitational_memory_source_adoption/check.py)默认只读，[结果](gravitational_memory_source_adoption/results.json)保1056三份输入SHA。实际核：

- 球面 $f^2$、$H[f]^2$、横向及无迹身份；
- 独立null射线TT核的16πG归一；
- 原三体总能源与自旋能流精确抵消；
- 从已冻结1056系数作同单位线性运输，没有重跑衰变或引力传播。

诊断取 $GE=1$ 只核无量纲系数；没有猜实际实验距离、信号强度或可探测性。默认检查通过；作者文件哈希另交独审签收。

采用的净覆盖是**同一个已算物质对象的渐近来源—记忆接口**。完整M4／M5及ROADMAP仍未完成；不因此启动更精细仪器或把全部UV／IR工程增加为当前目标门槛。

