# 真实 Majorana 标量能给出多少 Higgs 阈值修正

## 比较对象与前例

第536—542轮给出带权动能与四次矩的中性补全、轻质量限制与特定保护对的谱运行。第543轮已经用一个实 Majorana 标量、一只 top 和两个等强重极检验门户修正及共同负二次项，第544轮还给出在两个尺度间运行可打破二次相等的实例。本轮的新问题是，允许三代任意复矩阵与实际九份带电质量后，多少修正能由**同一个给右手中微子质量的场**产生；其物理质量能否落在第1805轮匹配所需的尺度。这也检验从正的共同物质迹建造不同物理系数的联动，而不改变认知原则。[早期中性补全](../../archive_531_553/research_note_536.md) · [早期保护对](../../archive_531_553/research_note_543.md) · [跨尺度竞争例](../../archive_531_553/research_note_544.md) · [第1805轮实际输入](../1805/proof.md)

取一个轻 Higgs 双重态，其规范化实径向场 $h^2=2H^\dagger H$，另给一个实规范单态 $s$，Majorana 质量来自 $sR$，$R=R^T$。可取至多三个实际右手态，没有另加固定 Majorana 质量。$s$ 的存在是额外物质输入；原来的普通有限代数内涨落无法单独产生这个分量。以下在平直背景、至维四的共同谱作用及单一树匹配面上检验，另行说明能标和误差。Chamseddine 与 Connes 已有该标量降低 Higgs 四次项的机制，其权、谱、运行和匹配条件与本轮有别。[原始模型](https://arxiv.org/abs/1208.1030) · [产生该场的扩展](https://arxiv.org/abs/1304.0415)

## 从同一有限算符算出三个耦合

写 $A=Y_\nu^\dagger Y_\nu$、$G=R^\dagger R$。两者都是作用在同一右手标架的正半定矩阵。内部对称中性质量块在左、右基上取

$$
{\cal M}(h,s)=\begin{pmatrix}
0 & hY_\nu/\sqrt{2}\\ hY_\nu^T/\sqrt{2} & sR
\end{pmatrix}.
$$

连同粒子和共轭，有限算符的偶次迹为 $2\operatorname{tr}({\cal M}^\dagger{\cal M})^j$。加入全部带电谱与颜色多重数，取第1805轮的 $\omega=w_\ell/w_q$、$S=2\pi^2/f_0$，规范化两场动能要求

$$
A_{\rm ch}+\omega\operatorname{tr}A=S,\qquad
\omega\operatorname{tr}G=S .
$$

记 $D=S-A_{\rm ch}$。规范化后势为

$$
V=-\frac12(\mu_H^2 h^2+\mu_s^2s^2)
+\frac14(\lambda_Hh^4+\lambda_p h^2s^2+\lambda_s s^4).
$$

同一热核的四阶矩及 $R$ 的对称性给出

$$
\lambda_H=\frac{B_{\rm ch}+\omega\operatorname{tr}A^2}{S},\qquad
\lambda_p=\frac{4\omega}{S}\operatorname{tr}(AG),\qquad
\lambda_s=\frac{2\omega}{S}\operatorname{tr}G^2 .
$$

$\operatorname{tr}A^2$ 表示 $\operatorname{tr}(A^2)$。在第543轮的记号中 $p=\lambda_p/2$。随机复、非对易的 $A,G$ 从显式对称质量块重算二次、四次多项式；五组背景最大绝对残差 $1.8\times10^{-15}$。同一热核的二次矩正是两个动能分别使用的平方迹，故若连相关算符也无独立反项，就有 $\mu_H^2=\mu_s^2=\mu^2>0$。这个要求和只让边际四次项共用谱矩可以分别检验。[复算](check.py) · [来源约定](../../archive_531_553/research_note_533.md)

## 门户修正的极限及其真实重谱

先在电弱对称背景 $h=0$ 展开。$s$ 有非零驻点 $w^2=\mu_s^2/\lambda_s$。树级径向质量与各右手质量为 $m_s^2=2\mu_s^2$、$M_i=w\sqrt{g_i}$，$g_i$ 是 $G$ 的本征值。消去有质量径向模的近零动量响应，得到

$$
\lambda_{\rm eff}=\lambda_H-\frac{\lambda_p^2}{4\lambda_s},\qquad
\mu_{\rm SM}^2=\mu_H^2-\frac{\lambda_p}{2\lambda_s}\mu_s^2.
$$

低能写作 $-\mu_{\rm SM}^2 H^\dagger H$。该四次匹配是 $h$ 小于改变 $s$ 相的场幅、外部动量远低于重极的树级系数；有限圈阈、外动量依赖、物理极与上述树质量的位移另需计算。

由 Hilbert–Schmidt 不等式，门户向下修正不超过 $2\omega\operatorname{tr}A^2/S$，因而

$$
\lambda_{\rm eff}\ge\frac{B_{\rm ch}-D^2/\omega}{S}.
$$

忽略额外轻质量与重极位置约束时，$A$ 趋近秩一且 $G\propto A$ 给出此界；若坚持三个重质量非零，等号仅作极限。中央全谱给 $\lambda_{\rm eff}\ge0.011611$，同尺实测量的 SM 一圈上行给 $-0.014942$，还差 $0.026553$。门户最大的单独修正 $0.310282$ 不能与第1805轮最小的 $\lambda_H$ 独立搭配，因为二者读取同一 $A$。这个中央界没有外推为任意方向的全数据盒结果。[中央数值](results.json)

对于三个相等的右手质量，$G=S I/(3\omega)$，向下修正精确等于 $2D^2/(3\omega S)=0.103427$。全部中性方向与奇异值都允许时，中央 $\lambda_{\rm eff}$ 仍只在 $[0.11504,0.21847]$。若只知 $M_{\max}/M_{\min}\le K$，则 $\operatorname{tr}G^2\ge g_{\max}^2(1+2/K^4)$，同时 $\operatorname{tr}(AG)\le g_{\max}\operatorname{tr}A$、$\operatorname{tr}A^2\ge(\operatorname{tr}A)^2/3$。结合上一段全域界，得到

$$
\lambda_{\rm eff}\ge
\frac{B_{\rm ch}}S-\frac{D^2}{\omega S}
\min\left(1,\frac{2}{1+2/K^4}-\frac13\right).
$$

第1805轮完整一宽盒里，等重三极时还差至少 $0.09688$，$K\le1.25$ 时还差至少 $0.02090$。更大分裂仅凭此外包可能失去分辨力。独审另以独立轻中微子数据检验强相消所需的实际重谱分裂。[全盒数值](results.json) · [独审](independent_review.md)

标量本身还有独立质量读数。消去 $w$ 与规范化后，

$$
m_s^2=4\,\frac{\sum_i M_i^4}{\sum_i M_i^2}.
$$

给定 $n$ 只右端、记 $x_i=M_i^2/M_{\max}^2$ 并置 $x_1=1$，有 $\sum x_i/\sum x_i^2\le(1+\sqrt n)/2$。可由对 $\sum(x_i-u)^2$ 配方，或对固定 $\sum_{i>1}x_i$ 将其余项均分后单变量微分得到。最大值由 $x_{i>1}=1/(1+\sqrt n)$ 达到。因此对三只右端有 $1.711199\le m_s/M_{\max}\le2$。全等时 $m_s=2M$。可以选择 $M_i=\mu_*$、$m_s=2\mu_*$，在共同重整化尺度 $\mu_*$ 积去，两组物理质量各自保留；有限匹配会含质量比的对数。高度分裂时须逐阈运行，不能用一段纯 SM 流替代多个活跃区。[显式达界谱](check.py)

## 共同二次项的独立弱尺度读数

一旦二次项也按共同谱作用识别，让右手尺度很重而 Higgs 保持轻，须使 $r=\lambda_p/(2\lambda_s)$ 很接近 $1$。对所有复相位、混合角与三个重谱，本模型却有

$$
r=\frac{\operatorname{tr}(AG)}{\operatorname{tr}G^2}
\le\frac{D}{S}\,\frac{g_{\max}\operatorname{tr}G}{\operatorname{tr}G^2}
\le \frac{D}{S}\frac{1+\sqrt3}{2}.
$$

中央信息给 $r\le0.626284$。在第1805轮全盒，用同一数据得 $D/S=1-A_{\rm ch}/S\le0.480595$，故 $r<0.65651$。这个比例不依赖弱尺度运行四次项的符号，也没有对哪个中微子态负责主要 Yukawa 权重作额外指定。忽略轻谱约束时，令 $A$ 指向最大重极、其余重极平方取上一段的最优比例，确实可达。

将这一界与另一个可测量量比较。纯 SM 下，$\mu_{\rm SM}^2$ 从顶质量标尺的运行是乘性的，

$$
16\pi^2\,\frac{d\ln|\mu_{\rm SM}^2|}{d\ln\mu}
=12\lambda+2T-\frac92g_2^2-\frac32g_Y^2 .
$$

原始多圈匹配的中央 $m(M_t)=131.55$ GeV，按其势约定 $|\mu_{\rm SM}^2|=m^2/2$。本次仅取极宽的 $|\mu_{\rm SM}^2(M_t)|\le(1\ {\rm TeV})^2$。用第1369轮 $T<3.001$ 及第1805轮全程 $\lambda<0.25$，在 $\mu_*$ 的允许区内最多放大 $2.823$ 倍。共同二次项则要求 $|\mu_{\rm SM}^2|=(1-r)m_s^2/2$，从而给 $m_s<4.055$ TeV，而注册的共同尺度至少 $1.1487\times10^{10}$ GeV。该冲突排的是将真实标量放在此高尺度并直接匹配的方案；它没有宣告一个 $4$ TeV 物理极存在。若果真将它放低，标量留在中间理论，必须另送整个谱到规范交点。[实测 Higgs 二次匹配式(56)](https://arxiv.org/abs/1307.3536) · [原规范输运](../../archive_1364_1409/1369/results.json) · [同一 Higgs 演化](../1805/check.py)

若 $\mu_H^2$ 和 $\mu_s^2$ 可分别给定，则须 $|\mu_H^2/\mu_s^2-r|\le4.28\times10^{-14}$ 才能达到上述极宽弱尺度条件（$m_s\ge\mu_{*,\min}$）。相对于原来两个相等的边界，必须先有至少约 $34.34\%$ 的分裂，再控制这份残差。若只增加相同正模块权的右手线、仍由同一 $s$ 和同一迹产生全部质量，则必要的计数条件为 $n\ge(2S/D-1)^2$。中央至少十二线，全盒至少十线；这一计数只是令 $r$ 达到一的必要条件，新的规范谱、Yukawa 来源和实际轻质量并未因此成立。具有任意权的额外不可见 singlet 或独立裸质量也不在此计数中。

所有四次系数非负，$h^2,s^2\ge0$ 上的树势因此有界。两凝聚的驻点还要求 Hessian 正定；其条件之一是 $\lambda_H\lambda_s-\lambda_p^2/4>0$。SM 从低能运行到高标尺的负 $\lambda$ 不是完整量子真空寿命计算，本轮只用同尺度匹配系数与独立的二次质量。若谱边界与标量脱耦相距很远、增加其它双重态/有色伙伴、有限大阈值或独立二次算符，第544轮已有运行分裂的理由，本轮同面界须改为实际多阶段输运。本轮不为这条受限候选继续引入自由量，新增认知公理零。

