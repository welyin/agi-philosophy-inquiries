# 1073独立审阅：稳定复核、完整仪器与反向分离

2026-10-10。独立构造；不导入本轮作者代码。本文只审1073已准入合同，不新增实际位置、反向操作或空间重定向权限。

## 1. 审阅结论与量词

**通过，须保留以下范围。** 同一固定二结果仪器在完整qubit任务上闭合，每一分支非零；每次都由原仪器实际重读，而不是重新显示旧记录。对所有输入及其每个非零分支，条件重读错误不超过 $\varepsilon<1/2$。该合同支持：

- 效果谱窗口与近固定迹；
- 到同一高本征投影所确定Lüders仪器的完整过程误差；
- 任意被动参考和有限自适应接续的误差预算；
- 实际完整壳上严格反向效果间隔大于 $\varepsilon$ 时的无迹分离。

这仍是新增认知候选到旧物理接口的**条件桥**。它不是认知必然性证明，也不是维数结论的全部前提。后态会被测量改变；本合同不承诺保留初次测量前的未知态。

## 2. 独立解析核对

记分支为 $\Phi_+,\Phi_-$，效果 $E=\Phi_+^*(I)$，$\Phi_-^*(I)=I-E$。全输入条件重读等价于

$$
\Phi_+^*(I-E)\le\varepsilon E,\qquad
\Phi_-^*(E)\le\varepsilon(I-E).
$$

这是算符序，不能由几个固定输入的统计代替。

### 2.1 谱、投影与漏出

设 $\alpha\le\beta$ 是 $E$ 的两本征值。每个分支至少有一份可达后态，分别给

$$
\alpha\le\varepsilon,\quad\beta\ge1-\varepsilon,\quad
|\operatorname{Tr}E-1|\le\varepsilon.
$$

由于 $\varepsilon<1/2$，高本征投影 $P$ 唯一、秩一；$P^\perp=I-P$。关键不等式为

$$
I-E\ge(1-\varepsilon)P^\perp,\qquad
E\ge(1-\varepsilon)P.
$$

因此，以 $q=\varepsilon/(1-\varepsilon)$ 记，

$$
D=\Phi_+^*(P^\perp)+\Phi_-^*(P)\le qI.
$$

这保证保留结果标记的任意参考扩展输出，在“$+$对应$P$、$-$对应$P^\perp$”子空间外的总概率不超过 $q$。这里只需对输入边缘应用算符不等式，没有假定输入与参考不纠缠。

### 2.2 压缩及完整过程误差

用 $\Lambda$ 表示保留标记和qubit后态的通道，用 $\Pi$ 投影到上述正确子空间，令 $\Gamma(\rho)=\Pi\Lambda(\rho)\Pi$。$\Gamma$ 是CP但通常不保迹。温和投影估计逐任意参考输入给

$$
\|\Lambda-\Gamma\|_\diamond\le2\sqrt q.
$$

设按原效果测量后重备为

$$
M_E(\rho)=|+\rangle\langle+|\otimes\operatorname{Tr}(E\rho)P
+|-\rangle\langle-|\otimes\operatorname{Tr}((I-E)\rho)P^\perp.
$$

因正确输出子空间每一块秩一，$M_E-\Gamma$ 是CP；其对偶作用于单位的结果正是 $D$，所以

$$
\|M_E-\Gamma\|_\diamond=\|D\|\le q.
$$

又 $\|E-P\|\le\varepsilon$，而同一标记、同一准备态的二值测量图满足

$$
\|M_E-\mathcal L_P\|_\diamond=2\|E-P\|\le2\varepsilon.
$$

故

$$
\frac12\|\Lambda-\mathcal L_P\|_\diamond
\le b(\varepsilon)
=\min\left\{1,\sqrt{\frac{\varepsilon}{1-\varepsilon}}
+\frac{\varepsilon}{2(1-\varepsilon)}+\varepsilon\right\}.
$$

这些估计不要求分支本身是重备通道，不追求最优常数。$\varepsilon=0$时，正确后态支撑和效果均为秩一投影，因此完整qubit仪器确为 $\mathcal L_P$。这不禁止全硬件中无后续反馈的环境，也不允许删去将会反馈的私有记忆。

若每次替换都满足相同完整输入／输出合同，后续控制、反馈和其他系统在两流程中相同，则逐槽替换、通道收缩性及三角不等式给最多 $N$ 次接续的半diamond误差不超过 $\min(1,Nb)$。无限重复、错误随历史无界增长、或未纳入任务的反馈记忆不在此结论内。

### 2.3 平方根障碍

取精确锐测 $P$，随后对两分支后态统一施

$$
V_\theta=\begin{pmatrix}\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta\end{pmatrix},\qquad
\varepsilon=\sin^2\theta.
$$

效果仍精确为 $P$，重读条件错误恰为 $\varepsilon$；但对**原效果投影所确定**的 $\mathcal L_P$，

$$
\tfrac12\|\Lambda_\theta-\mathcal L_P\|_\diamond
=\sin\theta=\sqrt\varepsilon.
$$

标记使两块正交；每块参考态乘上相同的纯态差，迹范数可直接相加。因此等式对任意输入参考态成立，无须数值优化。这里没有证明在所有其他投影中重新优化目标后的最优常数；仅此固定目标例已经排除一般 $O(\varepsilon)$ 完整过程保证。一次结果率甚至可以完全不变，而后态过程相差 $\sqrt\varepsilon$。

### 2.4 严格反向阈值

若 $E,F$ 都满足谱窗口，则

$$
\left|\frac{\operatorname{Tr}(E-F)}2\right|\le\varepsilon,\qquad
\|\operatorname{tf}(E-F)\|\ge\|E-F\|-\varepsilon.
$$

因此 $\|E-F\|>\varepsilon$ 足以推出无迹分离。锐阈值为

$$
E=(1-\varepsilon)P,\quad
F=P+\varepsilon P^\perp,\quad F-E=\varepsilon I.
$$

两者都有合法的按本征分支重备仪器满足原复核合同，但无迹部分完全相同。故不能把严格大于改为大于等于。

对原完整实际壳上每对实际反向点使用此结论，直接提供383和1067所需的无迹分离。壳、自由连续反向及同源动作仍须独立建立；有限样点不能认证全壳不等式。若效果图连续，统一谱隙使 $P_x$ 连续，但本轮无需把该辅助投影当作一个实际免费新仪器。

## 3. 独立数值设计与结果

[独立脚本](independent_check.py)默认只向标准输出打印JSON；没有作者代码导入、图像检查、额外依赖或默认写文件。原字节输出保存于[独立结果](independent_results.json)。

首轮复算：**1692项全部通过**，最大等式残差约 $1.12\times10^{-15}$，校验容差 $2\times10^{-11}$。

覆盖：

1. 七个 $\varepsilon$（含0与0.49）、三个非共轴基，共21组仪器。
2. 分支含真正保留输入相干的Kraus项：
   

$$
   \Phi_s(\rho)=(1-t)\operatorname{Tr}(E_s\rho)P_s
   +tU_s\sqrt{E_s}\rho\sqrt{E_s}U_s^\dagger.
   
$$

   $E_s$ 的两本征值为 $\varepsilon/4$ 与 $1-\varepsilon/4$，
   $t=0.9(\varepsilon-\varepsilon/4)/(1-\varepsilon/2)$。
   非零 $\varepsilon$ 下直接核对支配本征基非对角输入的响应，排除只用测后重备特例验证一般命题。
3. 全输入重读算符序、漏出、压缩图与重备图之差的CP性、效果投影界。
4. 独立的完整范数上证书：对差图的未归一Choi矩阵 $J$ 作正负分解，使用
   

$$
   \|\Delta\|_\diamond\le
   \|\operatorname{Tr}_{\rm out}|J(\Delta)|\|.
   
$$

   数值证书均不超过解析界；这不是仅以参考态抽样估计diamond范数。
5. 十份参考输入（含最大纠缠、混合态、固定种子的正矩阵）另核温和压缩及过程距离，抽样只作实现检查。
6. 六个后态旋转例的 $\sqrt\varepsilon$ 等式、五个偏置阈值例及一般谱窗口方向样本。
7. 最多4次保留完整结果历史、按旧结果选择下次仪器和反馈的接续。
8. 同一6维仪器中的三态记忆反例：首用对qubit的效果为 $0.3I$，保留记忆后可以零错误重放结果；如果把记忆抹掉后仍称其为完整qubit仪器，则下一次调用的合同已改变。反例说明真实任务闭包不是措辞要求。

数值浮点结果不替代解析的全输入、全参考或全历史量词。

## 4. 去重与选择力

- [383](../../archive_370_428/research_note_383.md)：自由实际反向、完整壳和固定迹／等价无迹分离的上界直接复用，不重证拓扑定理。
- [384](../../archive_370_428/research_note_384.md)：旧极小性合同已推出固定迹；不把它重新列成旧路线缺口。
- [1067](../1067/proof.md)第7节：当前删极小性的路线允许迹随位置变化，四维反例正是本轮新合同要区别的旧边界；上界接口不能偷领固定迹。
- [1069](../1069/proof.md)第3节：sharp读出和精确可重复性已给方向任务闭包。它没有证明每个实际定位效果都由该内部sharp菜单产生，也没有本轮有限错误到完整定位过程预算的连接。

故成熟谱条件与精确Lüders结论计复用；本轮可计的研究增量是**新增认知复核候选对当前实际壳仪器的有限过程约束及可判反向阈值**，不是声称可重复测量理论原创。

特别提醒：有限 $\varepsilon$ 不会让允许仪器族的数学维数自动坍缩。额外过程参数可仍然存在，只是对所列任务的影响受误差界限制。不能因此说五维通道族已变成真实三维空间。本轮三维上界的有效途径是同一完整效果图的严格反向阈值；空间下界、实际端点生成与固定基点稳定子来源仍按1067、1071及后1072审计保留。

