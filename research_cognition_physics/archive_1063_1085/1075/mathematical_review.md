# 1075 独立数学审阅

## 0. 审阅对象、结论与版本状态

**结论：最终证明、来源范围及以下解析命题通过，无必须修改项。** 本审阅已实际全文核读 proof.md 与 source_scope.md，并核验完整量子比特通道分类、纯输入谱界、逐对可达范围，以及“六个一致输入＋两个互换输入”的有限误差判据。数学核验为独立解析推导，并非仅转述作者的数值结果。

最终绑定的 proof.md SHA256：

`47982c8485a1e50b1b38633fed32868adbd527b873a5a3905067cb79b47d9168`

作者 check.py 已实际只读复跑：1164 次断言全部通过，最大残差 `7.771561172376096e-16`；stdout JSON 与已保存 results.json 完全一致。程序只作有限数值校验，不替代下面的解析证明。最初草稿在作者文件落盘前写成；本版已经取代该未绑定状态。

本轮可记为一项限定来源合同的兼容性校准：检验“双方对等、保留未知纯态一致、固定自主处理、确定性单量子比特输出”是否相容于始终得到纯态会合结果。成熟线性性、完全正性、对称子空间结构及既有密度矩阵平均通道均不另计原创成果。本报告不主张相关分类在全部文献中首次出现。

## 1. 合同必须明确的对象与量词

输入为两个量子比特，输出为同一个指定量子比特。固定通道为

$$
\Phi:\mathcal B(\mathbb C^2\otimes\mathbb C^2)\longrightarrow
\mathcal B(\mathbb C^2),
$$

并为完全正、迹保持（CPTP）。交换算符为 $S$，对称与反对称投影为

$$
P_+=(I+S)/2,\qquad P_-=(I-S)/2.
$$

完整分类采用以下两个条件：

1. **全通道交换对称：**
   

$$
   \Phi(S\omega S)=\Phi(\omega)\quad\text{对全部输入 }\omega.
   
$$

   这不是仅在若干抽样输入上检验输出相等，也不是另一个具有两个输出的交换协变条件。

2. **未知纯态一致保持：**
   

$$
   \Phi(P\otimes P)=P
   \quad\text{对全部秩一量子比特态 }P.
   
$$

也可将第二项的“全部”换为纯态球面 $\mathbb CP^1\simeq S^2$ 中一个非空开片上的全部纯态，然后用多项式恒等延拓。这里的开片是二维纯态流形的开集；仅一条实大圆上的区间不足以直接调用该延拓论证。

固定辅助态、固定环境、内部测量及自适应处理均可吸收进上述固定 CPTP 通道。辅助系统不能额外取得未知输入的经典描述或额外副本；最终输出不能只在后选成功分支上定义。若保留外部记录，结论仍约束迹掉该记录后的指定输出量子比特；它不自动分类记录与输出构成的整个实际任务。输入相关程序、额外活动记录、保留多个输出或允许混合输出属于另一份合同。

## 2. 完整分类：唯一自由量是反对称扇区输出态

**命题成立。** 恰好所有符合第 1 节合同的通道为

$$
\boxed{\displaystyle
\Phi_\tau(\omega)=
\operatorname{Tr}_2(P_+\omega P_+)
+\operatorname{Tr}(P_-\omega)\,\tau,
\qquad \tau\in\mathcal D(\mathbb C^2).}
\tag{2.1}
$$

### 2.1 交换对称消去交叉块

对称—反对称交叉算符在 $S(\cdot)S$ 下变号。因此交换对称与线性性给出

$$
\Phi(P_+\omega P_-)=\Phi(P_-\omega P_+)=0.
$$

这一步采用全通道交换对称；不能从有限样点直接推出。

### 2.2 纯态一致输入张满对称扇区的算符空间

在对称空间基

$$
|00\rangle,\quad (|01\rangle+|10\rangle)/\sqrt2,\quad |11\rangle
$$

中，未归一化向量 $(|0\rangle+z|1\rangle)^{\otimes2}$ 的坐标为

$$
(1,\sqrt2 z,z^2).
$$

其投影矩阵的各项给出 $z^a\bar z^b$，$0\le a,b\le2$。这些多项式线性独立，故 $P^{\otimes2}$ 的实线性张成空间是全部
$\operatorname{Herm}(\mathrm{Sym}^2\mathbb C^2)$。

更直接地，对任意正交于所有这些投影的厄米算符，其上述多项式恒为零，只能是零算符。若一致条件仅在非空纯态开片成立，多项式在对应的实二维开集为零，仍推出恒等式。

由于部分迹通道也满足

$$
\operatorname{Tr}_2(P^{\otimes2})=P,
$$

$\Phi$ 与部分迹在整个对称扇区上一致。

### 2.3 反对称扇区及逆向验证

两量子比特反对称扇区是一维，且 $P_-$ 本身是归一化纯态。因此

$$
\tau=\Phi(P_-)
$$

是任意可能的量子比特密度矩阵；这个扇区上的作用只能为
$\omega\mapsto\operatorname{Tr}(P_-\omega)\tau$。

反过来，式 (2.1) 是两个完全正映射之和，且两项的迹相加为输入迹。它满足交换对称与纯态一致保持。因此分类同时具有必要性与充分性，而非只给出若干实现例子。

一个等价表达式是

$$
\Phi_\tau(\omega)=
\frac{\operatorname{Tr}_1\omega+\operatorname{Tr}_2\omega}{2}
+\operatorname{Tr}(P_-\omega)\left(\tau-\frac I2\right).
\tag{2.2}
$$

这里两份部分迹的输出空间按固定量子比特接口作同一识别。

## 3. 任意纯输入对的谱界及逐对可达性

取

$$
\rho=|u\rangle\langle u|,\qquad
\sigma=|v\rangle\langle v|,\qquad
c=|\langle u|v\rangle|\in[0,1].
$$

定义半迹距离 $D(A,B)=\tfrac12\|A-B\|_1$。

对称部分的输出为

$$
R=\operatorname{Tr}_2\!\left(P_+(\rho\otimes\sigma)P_+\right)
=\frac14(\rho+\sigma+\rho\sigma+\sigma\rho),
$$

其两个特征值为

$$
\lambda_\pm(R)=\frac{(1\pm c)^2}{4}.
$$

反对称部分权重为

$$
a=\operatorname{Tr}\!\left(P_-(\rho\otimes\sigma)\right)
=\frac{1-c^2}{2}.
$$

于是

$$
\Phi_\tau(\rho\otimes\sigma)=R+a\tau.
$$

因 $a\tau\ge0$，输出最小特征值至少为 $(1-c)^2/4$。对任意归一化量子比特态 $\eta$ 与纯态 $\pi$，

$$
D(\eta,\pi)\ge1-\operatorname{Tr}(\pi\eta)
\ge1-\lambda_{\max}(\eta)=\lambda_{\min}(\eta).
$$

因此

$$
\boxed{\displaystyle
D\!\left(\Phi_\tau(\rho\otimes\sigma),\pi\right)
\ge\frac{(1-c)^2}{4}
\quad\text{对每个纯目标 }\pi.}
\tag{3.1}
$$

只要两输入纯态不同，即 $c<1$，输出就不可能是任何纯态。该结论无需指定“应当会合到哪一个纯态”。

**可达性的量词：** 对每一对预先选定的输入，可以选择一个固定 $\tau$，令它等于 $R$ 的最大特征值本征纯态，并以该纯态为目标，从而达到式 (3.1) 的下界。这里是“逐对选择通道参数后的可达性”；不意味着一个固定 $\tau$ 对所有输入对同时达到下界，也不授权执行时根据未知输入描述调整 $\tau$。

## 4. 共同 SU(2) 协变的额外结论

如果另要求

$$
\Phi((U\otimes U)\omega(U^\dagger\otimes U^\dagger))
=U\Phi(\omega)U^\dagger
$$

对全部 $U\in SU(2)$ 成立，则反对称纯态 $P_-$ 的不变性给出

$$
U\tau U^\dagger=\tau\quad\text{对全部 }U,
$$

故 $\tau=I/2$。

此时式 (2.2) 化为成熟的边缘平均通道：

$$
\Phi(\omega)=\frac{\operatorname{Tr}_1\omega+\operatorname{Tr}_2\omega}{2},
\qquad
\Phi(\rho\otimes\sigma)=\frac{\rho+\sigma}{2}.
$$

对纯输入对，其特征值为 $(1\pm c)/2$，故

$$
\boxed{\displaystyle
\min_{\pi\text{ pure}}
D\!\left(\Phi(\rho\otimes\sigma),\pi\right)=\frac{1-c}{2}.}
\tag{4.1}
$$

这是到**最近**纯态的距离；不能将等号写成对每个纯目标均成立。

共同 SU(2) 协变是额外输入，不能仅由交换对称或纯態一致保持推出。平均通道本身也不是本轮新构造。

## 5. 六个一致输入与两个互换输入的有限判据

本节独立于第 2 节的全域分类，**不要求全通道交换对称，也不要求全部纯态一致保持**。

设固定 CPTP 通道 $\Phi$ 满足六个实际输入的误差界

$$
D\!\left(\Phi(P_j^{\otimes2}),P_j\right)\le\delta,
\qquad j\in\{+x,-x,+y,-y,+z,-z\},
\tag{5.1}
$$

以及两个实际互换输入的输出差异界

$$
D\!\left(\Phi(|01\rangle\langle01|),
         \Phi(|10\rangle\langle10|)\right)\le\kappa.
\tag{5.2}
$$

这里 $P_j$ 为六个 Pauli 本征纯态；$\delta,\kappa$ 是对应实际输出半迹距离的上界，不是有限样本频数自动给出的置信结论。

### 5.1 六态分解

令 $|\psi^+\rangle=(|01\rangle+|10\rangle)/\sqrt2$，则精确恒等式为

$$
A=\frac12|\psi^+\rangle\langle\psi^+|
=\frac14\left(
\sum_{j=\pm x,\pm y}P_j^{\otimes2}
-\sum_{j=\pm z}P_j^{\otimes2}\right).
\tag{5.3}
$$

例如两边均等于

$$
\frac18(I\otimes I+X\otimes X+Y\otimes Y-Z\otimes Z).
$$

由式 (5.1)、线性性及三角不等式，

$$
\Phi(A)=\frac I4+\Delta,\qquad
\operatorname{Tr}\Delta=0,\qquad
\frac12\|\Delta\|_1\le\frac32\delta.
$$

由于 $\Delta$ 是二阶无迹厄米矩阵，

$$
\|\Delta\|_{\mathrm{op}}=\frac12\|\Delta\|_1\le\frac32\delta.
\tag{5.4}
$$

$\Delta$ 不是密度矩阵，因此此处应写半迹范数，避免将其未经说明称为两量子态之间的距离。

### 5.2 局部交换误差足够

记

$$
\eta_{01}=\Phi(|01\rangle\langle01|),\quad
\eta_{10}=\Phi(|10\rangle\langle10|),\quad
\tau=\Phi(P_-).
$$

不使用任何全域交换条件，就有

$$
\bar\eta=\frac{\eta_{01}+\eta_{10}}2
=\Phi(A)+\frac{\tau}{2}
=\frac I4+\Delta+\frac{\tau}{2}.
$$

因此

$$
\lambda_{\max}(\bar\eta)\le\frac34+\frac32\delta.
$$

式 (5.2) 又给出

$$
D(\eta_{01},\bar\eta)
=D(\eta_{10},\bar\eta)
\le\frac{\kappa}{2}.
$$

所以对两种输入方向、对每个纯态目标 $\pi$，均有

$$
\boxed{\displaystyle
D(\eta_{01},\pi),\ D(\eta_{10},\pi)
\ge\max\left\{0,\frac14-\frac32\delta-\frac{\kappa}{2}\right\}.}
\tag{5.5}
$$

若确有全通道交换对称，则 $\kappa=0$，但有限判据只需要这两个实际输出之间的约束。取 $\delta=0.01,\kappa=0$ 得到至少 $0.235$ 的距离下界。这是有限精度、有限输入的可检验预测差，不依赖对一个连续开片实行精确验证。本报告不声称 $\delta,\kappa$ 系数是全局最优。

## 6. 项目既有结果与文献合同差异

### 6.1 项目去重

旧 [429](D:/workspace/AGI的哲学思考/research_cognition_physics/archive_429_466/research_note_429.md) 的主要合同是封闭两体幺正演化、全部混态的相同输入、两个输出边缘均保持原态，由此约束为部分交换结构。它没有直接给出本轮的“任意二入一出 CPTP、仅纯态一致保持、单输出交换对称”的 $\tau$ 分类，也没有式 (5.5) 的有限见证。

此前“全部密度矩阵一致保持＋交换对称”经极化得到算术平均的论证应继续按成熟结果处理。仅要求纯态一致保持时，反对称扇区的自由态 $\tau$ 确实保留，不能无条件删去。

### 6.2 已核原始文献

- [Alvarez-Rodriguez 等，The Forbidden Quantum Adder](https://arxiv.org/pdf/1411.4534)：已有可实现的密度矩阵平均通道；其主要未知态禁阻针对指定的态矢量相加。该结果应引用并复用，不能把平均通道或一般“量子加法不可能”重新计为本轮成果。
- [Oszmaniec 等，Creation of superposition of unknown quantum states](https://arxiv.org/pdf/1505.04955)，定理 1：禁止对任意未知纯输入产生固定非零系数的相干叠加，甚至考虑概率完全正实现及相位代表的选择。它指定了输出叠加规则；本轮式 (3.1) 对任意纯输出目标成立，但另加纯态一致保持与交换对称。二者合同不同。
- [Ticozzi，Symmetrizing quantum dynamics beyond gossip-type algorithms](https://arxiv.org/pdf/1509.01621)：量子网络共识与所有输入上的纯输出/守恒限制提供成熟背景；其量词不等于本轮纯乘积输入上的二入一出条件。
- 协作文献复核还比对了 [Mazzarella–Sarlette–Ticozzi 的量子共识](https://arxiv.org/pdf/1303.4077) 及 [Czartowski 等的量子态同步](https://arxiv.org/pdf/2103.02031)。相同边缘、网络同步及平均保持均不能直接替代本轮未知纯态一致保持合同。

本次检索未找到以相同量词陈述的完整 $\tau$ 分类或式 (5.5) 的有限判据。这只支持明确区别来源合同，**不足以据此宣称文献首创**。本报告是项目内独立数学审核，不是外部同行评审。

## 7. 来源解释与未证范围

本轮可以排除下述具体合取：两份未知纯量子比特状态作为全部可用输入；双方对等且保留已经形成的未知纯态一致；以一个固定、确定性自主过程给出单个无条件量子比特结果；并要求所有不同纯输入也总是得到纯态会合结果。式 (5.5) 进一步给出这类来源设想的有限误差检验。

该合取本身是待检验的认知来源假说，不能从一般“会合”“共同承诺”或“可接续”直接推出。

本结果没有证明这些量子态就是实际位置，没有生成实际接触等价类、实际位置更新、完整不变壳、位置协变或三维空间。它也不反证已验收的 [1074](D:/workspace/AGI的哲学思考/research_cognition_physics/archive_1063_1085/1074/proof.md)：1074 的实际会合定义在同一实际位置对象上，采用连续性、实际作用协变与局部唯一性等合同，并未要求通过上述固定二入一出量子通道将所有未知纯输入变成纯输出。

允许混合结果、额外活动历史/记录、具有输入相关资源的过程或更丰富的关系输出，均需另列实际任务条件。本轮不对这些未覆盖机制下结论，也不将修复这一特定量子接口合同提升为所有空间来源必须经过的门槛。

## 8. 签收范围

- 完整 $\tau$ 分类：通过。
- 纯输入谱界及逐对可达限定：通过。
- 共同 SU(2) 协变后的最近纯态距离：通过。
- 六个一致输入＋两个互换输入的 $(\delta,\kappa)$ 判据：通过。
- 认知来源与实际空间的范围隔离：按上述限定通过。
- 最终 proof.md 与 source_scope.md：实际全文核读，通过，证明字节绑定如第 0 节。
- 作者 check.py/results.json：实际核读、只读复跑及结构化结果逐项比对通过；不将数值抽样视为全域定理。

作者代码与结果的 SHA256：

- check.py：`db528c78b7c22a6315b817da99090f705c0d84b4736da5cdcf3246a89064fd8d`
- results.json：`52e73c09330ce64118aa4038a6f967a52d90a3e31d78bce6821b4893bb2a6c33`
- source_scope.md：`942a79a4be09a9e1abcb702076824ed1f5f4e84f281d4aa4612c6f4ff8e56ffd`

本次仅写入本 mathematical_review.md，没有修改主证明、作者程序、结果或任何冻结轮次；没有以数值脚本替代上述解析证明。
