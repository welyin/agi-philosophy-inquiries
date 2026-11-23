# 独审　中性弱流的零三次矩使希格斯阈值有更强的下界

在主文规定的单双重态非超对称费米子匹配里，重质量的 Takagi 基为 $m=\operatorname{diag}(M_i)>0$，Dirac 质量端口为 $d=uY$，取 $F=d\,m^{-1}$。忽略更高阶电弱混合、并在同一标尺和端口比较，则

$$
W=FF^\dagger,\quad m_\nu=-F m F^T,\quad T=\operatorname{Tr}(F m^4 F^\dagger)=8\pi^2u^2\Delta m_H^2 .
$$

最后一个等号沿用第1558轮的逐极有限匹配，取相同正权。$W$ 是六维弱流算符乘上电弱真空因子的权重。在第1689轮，规范基中质量的 $NN$ 块为零，且只有 $N$ 直接读到活性端口，故 $F m^3 F^T=0$。性能相同的任意重态数都可应用以下结论。[原图及动量矩](../1689/proof.md) · [希格斯匹配定义](../../archive_1545_1585/1558/proof.md) · [弱流低能约定](../1689/independent_review.md)

## 精确常数

对任意归一活性方向取行向量 $f$，记 $w=\sum_i|f_i|^2,\ t=\sum_i|f_i|^2M_i^4,\ a=|\sum_i f_i^2 M_i|$。$w=0$ 时该行平凡。对每个 $x\ge0,\ s>0$，有

$$
\left|x-\frac{x^3}{s^2}\right|\le\frac{x^4+s^4}{2\sqrt2\,s^3}.
$$

把右侧与左侧括号的两个符号分别相加减，乘 $2\sqrt2 s^3$ 后恰为 $(x^2\pm\sqrt2 s x-s^2)^2$。利用 $\sum_i f_i^2M_i^3=0$，得 $a\le(t+s^4w)/(2\sqrt2s^3)$。取 $s^4=3t/w$，给出 $t w^3\ge27a^4/4$。在轻质量 Takagi 基逐行相加，使用 $\sum_j a_j^4/w_j^3\ge(\sum_j a_j)^4/(\sum_j w_j)^3$，最终得到

$$
T(\operatorname{Tr}W)^3\ \ge\ \frac{27}{4}\left(\sum_j m_{\nu j}\right)^4.
$$

常数可达到。取每个有质量轻方向各接一个实重块 $\left(\begin{smallmatrix}0&M\\M&\sqrt2 M\end{smallmatrix}\right)$，共用 $M>0$，各行端口强度的平方取正比于目标轻质量。逐行和迹的两个不等式都取等。没有零三次矩限制时，单用加权 Hölder 只能统一保证常数1；任意多个等质量的直接中介和正交的实端口即饱和。这两个最佳常数覆盖相应矩阵类别，取等构造无须具备原候选指定的每个味顶角。原图具体归一比为318.694，另有32个随机复图在独立重场基旋转后通过数值检验。[复算](independent_check.py) · [结果](independent_results.json)

## 同一物理端口的可辨量

在 $N,H$ 基，重矩阵为 $\left(\begin{smallmatrix}0&C\\C^T&H\end{smallmatrix}\right)$，活性只接秩二的 $D$。令 $P=C^*C^T$，则 $T=\operatorname{Tr}(DPD^\dagger)$，相对于同 $D$ 同 Schur 质量的直接模型，$\delta W=DP^{-1}D^\dagger$。利用满列秩得到

$$
 T=\operatorname{Tr}\!\big[(DD^\dagger)(\delta W)^+(DD^\dagger)\big].
$$

由此，若直接端口与弱流增量均可独立取得，可反推出正希格斯谱矩；单凭轻质量矩阵不足以做到这一点。用第1678轮质量见证、$u=174$ GeV、零独立希格斯边界并让该部门供应 $125.11^2$ GeV$^2$ 的全部有限匹配，尖锐下界给 $\operatorname{Tr}W\ge1.2892\times10^{-17}$。运输权放宽至 $[1/2,2]$，用最小权得到 $1.0232\times10^{-17}$。数值是给定输入后的必要条件。独立边界、软超对称破缺阈值及更高阶匹配会改变从希格斯极点到正谱矩的联系；零三次矩和正矩不等式可由相应中性矩阵本身单独检验。
