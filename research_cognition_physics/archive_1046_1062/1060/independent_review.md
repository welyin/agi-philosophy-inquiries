# 1060独立终审：Pati–Salam全局映射与完整背景拉回

2026-10-09。审阅者：`/root/next_structural_selection`。**结论：通过，无待修科学问题。** 本审阅独立于作者使用的偶Fock代码；不增加科学校准组、认知公理或经验组。主线验收和整体ROADMAP完成情况另判。

已逐节读[正式报告](../research_note_1060.md)、[证明](proof.md)、[准入](selection.md)、[来源](sources.md)、[依赖账](dependency_update.md)及科学代码和保存结果。作者已按审阅建议，在证明§5明确区分所选Higgs真空中的满秩质量与所有拓扑背景上的异常消去。最终快照见§6。

## 1. 独立的全局群映射

采用Clifford乘积分解对应的真实子群

$$
P=\frac{SU(4)\times SU(2)_L\times SU(2)_R}
{\langle(-I_4,-I_2,-I_2)\rangle}\subset\mathrm{Spin}(10),
\qquad16=(4,2,1)\oplus(\bar4,1,2).
$$

记 $a=1+6k$。在覆盖群上直接定义

$$
\widetilde\phi_k(C,W,e^{i\theta})=
\left[
\operatorname{diag}(e^{ia\theta}C,e^{-3ia\theta}),\;
W,\;\operatorname{diag}(e^{-3i\theta},e^{3i\theta})
\right].
\tag{R1}
$$

所有权都是整数，故有真实 $2\pi$ 周期；三个块是群同态。若其在 $P$ 中为单位，则存在 $\epsilon\in\{0,1\}$，使三个因子同时等于 $(-1)^\epsilon$。第三因子给 $\theta=j\pi/3$，且 $j\equiv\epsilon\pmod2$。第一和第二因子遂强制

$$
C=e^{i(3-a)\theta}I_3=e^{2i\theta}I_3,
\qquad W=(-1)^jI_2=e^{-3i\theta}I_2.
$$

反过来，这六个元素在式(R1)中确实是 $((-1)^jI_4,(-1)^jI_2,(-1)^jI_2)$。所以核**恰为**原固定 $Z_6$，得到 $\phi_k:G_6\hookrightarrow\mathrm{Spin}(10)$。此证明既不只检查荷的模6同余，也不把向量表示中的单位误作Spin提升的单位。

## 2. 六类荷、中心和Higgs

在 $SU(4)$ 中取 $b=\operatorname{diag}(1,1,1,-3)$，则式(R1)的荷生成元为

$$
q_k=(1+6k)b+6T_R^3.
$$

$(4,2,1)$ 给 $Q:1+6k$、$L:-3-18k$；$(\bar4,1,2)$ 的两个 $6T_R^3=\mp3$ 权给

$$
u^c:-4-6k,\quad d^c:2-6k,\quad N:18k,\quad e^c:6+18k.
$$

这是1035完整16分量的同一表示。$10=(6,1,1)\oplus(1,2,2)$ 中的正荷弱双重态具有 $q_H=3$，不随 $k$ 改变。

取

$$
c=[iI_4,I_2,-I_2].
$$

它在 $(4,2,1)$ 上作用 $i$，在 $(\bar4,1,2)$ 上作用 $(-i)(-1)=i$；在向量10两块上均为 $-1$。其平方 $[-I_4,I_2,I_2]$ 在16上为 $-1$、在10上为 $+1$，正是Spin中心的 $-1$。这与作者Fock中心 $\exp(i\pi X_0/2)$ 的作用一致。

不需新增实际Spin(10)规范玻色子，亦不需补齐完整10标量或统一Yukawa。上式只提供原费米表示及Higgs子表示的认证容器。

## 3. 全部声明背景的运输

有与切向 $SO(d)$ 投影相容的群同态

$$
([s,r^j],g)\longmapsto[s,c^j\phi_k(g)]
\quad\text{从}\quad
\mathrm{Spin}^{Z_4}(d)\times G_6
\quad\text{到}\quad
(\mathrm{Spin}(d)\times\mathrm{Spin}(10))/Z_2.
\tag{R2}
$$

$(-1,r^2)$ 被送到被商掉的 $(-1,-1)$，所以(R2)在真实商群上成立。任意原结构主丛均可延伸结构群，费米关联束、联络、Dirac算符和异常沿同一个映射拉回。这里没有先要求原 $G_6$ 丛提升到直积覆盖群，也没有补选独立spin结构。因此作者的非提升规范丛及相容非spin背景量词正确。

已重读[Wang–Wen–Witten，arXiv:1810.00844v4 §5.1.3](https://arxiv.org/html/1810.00844v4#S5.SS1.SSS3)。原文处理的正是上述扭合Spin–Spin(10)背景：它既检查微扰异常，也检查可存在的 $w_2w_3$ 五维异常，普通16在该项上的系数为零。母背景分类并非因而变成零群；平凡的是当前费米表示的异常。这个成熟结论的自然拉回逐代为零，故三代亦为零。

此处签收的是四维局部异常及闭五维异常相位，不把它等同于带边界完整仪器构造、任意高形式背景、Pin/时间反演或非微扰UV存在定理。也不声称每个非spin流形都允许原 $\mathrm{Spin}^{Z_4}\times G_6$ 结构。

## 4. 独立复算与作者入口

[独立程序](independent_check.py)不导入作者或614函数，直接构造Pati–Salam的16与10矩阵。它核了：

- 六类荷的常数项与 $k$ 系数；
- 18个有理中心相位控制及其全 $k$ 整数系数身份；
- 12组非Abelian群乘法、周期、Higgs子表示和幺正性诊断；
- $c$、$c^2$、Pati–Salam商掉的对角中心及切向中心抵消。

[独立结果](independent_checks.json)全部通过，最大矩阵残差 $4.45\times10^{-14}$。有限采样不是全整数或全背景证明；全整数核由§1解析证明，全背景异常由§3承担。

作者 `global_anomaly_pullback.py` 默认只读运行通过，与已有结果一致。其42组矩阵样本的最大残差约 $2.84\times10^{-13}$。本审阅没有运行其保存选项、没有改变作者资产或历史。

## 5. 净增量、限定自由与停止线

- 614只给原 $k=0$ 容器；627仅采用普通spin的旧SU(2)检测；628完成原标准荷、普通spin的全局规范丛接口。这些旧结论不被重新计作缺口。
- 本轮真正增加的是1035**全部整数族**到更广Spin–$Z_4$联合背景的同态，因而排除了“该类尚未完成的混合费米异常或许只留下 $k=0$”这一选择路线。
- 所选满秩Yukawa／质量分支不表示每个非平凡背景都有处处非零Higgs或同一有隙真空。证明§5已明确这点。
- 仍采用规范群、16物种、四维Weyl形式及额外 $r$。$k\ne0$ 的带荷 $N$ 不是现实标准模型，亦未认证完整认知A合同；不能把此结果写成全部认知理论或现实物理的不唯一性。
- 原中性实 $SNN$ 条件和经验电荷约束继续有效。对这条整数族应按准入停止，不再靠更多群元、拓扑实例或GUT修补增加轮次。

## 6. 签收作者快照

七份作者资产的SHA256如下；独立结果也保存同一快照。后续主线收据应单独登记，不由此表提前认证尚未写成的文件。

|资产|SHA256|
|---|---|
|research_note_1060.md|`c64282eb868bbc2f34fac94ea20b950d880616921ec6b3233e875f775394b307`|
|proof.md|`67c446380245f91b96a0814591e923823ddc976d74be9a044723b7f77a3a1c6a`|
|selection.md|`664e60a0f6da674fafac109746f7f7efea38694aa486c7d46d4133c31b537882`|
|sources.md|`b6c48765c0ab89ca8f52a454ae354a8dc6f23bca26fd8be28eb2df60dfeb621a`|
|dependency_update.md|`1654b065c67d47ac8cf3c550033d487b91075b52493da5c10faa99d145621d01`|
|global_anomaly_pullback.py|`b0862c6af60eeace84e5e529fd8324ad93c3f7ea9b7f60b69c88b22bb616b15d`|
|results.json|`165fc7284ed44b596c07be224745b6203f7bc2443c5e6a9016672f82037f33cd`|
