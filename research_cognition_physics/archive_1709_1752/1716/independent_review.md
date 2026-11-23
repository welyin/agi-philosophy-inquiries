# 第1716轮独审　两项振荡频率对补充味方向的必要条件

在第1712轮同一带电成功点，取两源 Dirac 图实际运行至 $10$ PeV 后的左零线 $v$。按第1702轮的记号 $N$ 存储 Dirac 质量列的复共轭，所以 $\|N^\dagger v\|=1.2\times10^{-17}$；电荷本征基中此线的电子分率为 $0.67034$。不使用任何轻子混合角约束，只用两个振荡频率，仍能否定把这条线交给第1711轮最小 Schur 补的候选。下面的界覆盖三条上型行相位，无需相位扫描。[原 Schur 实现](../1711/independent_review.md) · [带电共同点](../1712/results.json) · [尺度合同](../1715/inputs.json)

## 从两个非零质量比去掉全部自由行相位

令 $D$ 为同一味基的下型 Yukawa，$u_i$ 为上型正奇异值，$C=D^T\operatorname{diag}(e^{-i\alpha_i}/u_i)D$。取 $W$ 的两列为 $v^\perp$ 的 Hermitian 正交基，$X=D^{-1}\operatorname{diag}(\sqrt{u_i})$，$x_i$ 为其列。若 $v^TCv\ne0$，令

$$
C_v=C-\frac{Cv\,v^TC}{v^TCv},\qquad
K=(W^TC_vW)^{-1}=W^\dagger C^{-1}W^*
=\sum_i e^{i\alpha_i}f_if_i^T,\qquad f_i=W^\dagger x_i.
$$

中间恒等式由在酉基 $(W,v)$ 中取 Schur 补的块逆公式得到。$K$ 的大、小奇异值分别是 $C_v$ 两个非零奇异值的倒数，因而所求质量比 $R=m_3/m_2=\|K\|_2^2/|\det K|$。

对每个 $j$ 令 $L_j=x_jx_j^\dagger-\sum_{i\ne j}x_ix_i^\dagger$。任一单位二维列 $a$ 给出

$$
\|K\|_2\ge |a^\dagger K a^*|
\ge |a^\dagger f_j|^2-\sum_{i\ne j}|a^\dagger f_i|^2 .
$$

故 $\ell=\max_j\lambda_{\max}(W^\dagger L_jW)$ 是统一下界，实际由粲行达到，$\ell=7243.53219$。另一方面，Cauchy–Binet 展开 $\det K$，三项振幅模为 $|\det[f_i,f_j]|^2=(668529.75,1118707.90,4931521.96)$。因此对于全部相位

$$
3144284.31\le|\det K|\le6718759.62,\qquad R\ge7.8092924.
$$

下界还排除了 $v^TCv=0$ 的例外，因为未先取逆即可对 $W^\dagger C^{-1}W^*$ 使用上述行列式展开，其模等于 $|\det C^{-1}|\,|v^TCv|$。本证书与给裸质量选择的公共归一无关。

正常序、严格一条零质量线时，两项 JUNO 一宽边际输入均满足的必要条件为 $R\le\sqrt{0.002536/0.0000731}=5.890009$。若有效质量从高至低按 $G^TC_vG$ 输运，奇异值乘积不等式给 $R_{\rm low}\ge R/\operatorname{cond}(G)^2$。为仅凭输运挽救此候选，至少需要 $\operatorname{cond}(G)\ge1.15146$。第1715轮登记的有活跃 Dirac 端、无独立轻子数破缺阈值的一圈合同比为 $1.001141$，给低端 $R\ge7.79150$。另取从高端直接运行 Weinberg 算符、无活跃重 Yukawa 的竞争类别，已有统一界给 $R_{\rm low}\ge7.80768$。有限味阈值、独立附加质量源、带电输入盒外的图或改变排序均须重新核验。[频率口径](../../archive_1632_1674/1668/inputs.json) · [无活跃重态的一圈界](../1715/independent_review.md) · [本轮数值](independent_results.json)

## 改哪条线、至少改多少

另外给任意候选方向 $v'$ 的连续稳定界。记两线的射影角为 $\theta=\arccos|v^\dagger v'|$，令 $L$ 取上节胜出的粲行。将 $v^\perp$ 中取得 $\ell$ 的单位向量投影回 $v'^\perp$ 并归一，其位移角至多 $\theta$，Hermitian 期望值改变至多 $\operatorname{spread}(L)\sin\theta$，其中 spread 为最大减最小本征值。又令

$$
Q=\sum_{i<j}\overline{x_i\times x_j}\,(\overline{x_i\times x_j})^\dagger,\qquad d=v^\dagger Qv .
$$

外积余子式恒等式给 $\sum_{i<j}|\det[f_i,f_j]|^2=d$。改变零线时，$v'^\dagger Qv'$ 至多增 $\operatorname{spread}(Q)\sin\theta$。于是

$$
R_{\rm low}(v')\ge
\frac{[\ell-\operatorname{spread}(L)\sin\theta]_+^2}
{\big(d+\operatorname{spread}(Q)\sin\theta\big)\operatorname{cond}(G)^2}.
$$

这里 $\operatorname{spread}(L)=30718.37$、$\operatorname{spread}(Q)=35137574.63$。只用已登记的有源一圈输运时，要进入频率盒，新的具名方向至少要转 **$1.1584^\circ$**；若另给满足 $\operatorname{cond}(G)\le1.01$ 的任意合同，上述必要角度仍为 **$1.0863^\circ$**。这允许改变方向以后重新选择全部相位。若方向改变同时改变了有源 RG，应核验新的条件数预算。上述角度为必要下界，尚未给出能达到窗口的替代方向。

若在低能端直接加任意质量扰动，按 Weyl 不等式，若原低端谱比下界为 $R_0$，新的两个最大奇异值欲达比值 $R_t$，则 $\|\delta M\|_2/m_3\ge(1-R_t/R_0)/(1+R_t)$。对有源输运为 **3.542%**。保持 $m_1=0$ 时才可从两项频率直接采用上述 $R_t$。尤其，四单态图只给新增中介添一个同向自身质量 $\mu$，其有效矩阵成为 $C_\mu=C+(Cv)(Cv)^T/(\mu-v^TCv)$，满足 $\det C_\mu/\det C=\mu/(\mu-v^TCv)$。在重块可积出的范围内，如继续要求精确秩二，只能取 $\mu=0$。允许第三个非零质量是另一个可检验的模型类别。[复算程序](independent_check.py)
