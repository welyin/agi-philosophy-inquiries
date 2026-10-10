# 绝对弱归一的共同匹配：μ 子寿命到 W 质量（成熟采用，计 0）

## 1. 本批补什么

在 [P981](../../archive_956_989/research_note_981.md) 的**可重整化标准模型分支**中，低能四费米强度与规范/Higgs 部门使用同一组参数。树阶有 $G_\mu=g^2/(4\sqrt2 M_W^2)=1/(\sqrt2 v^2)$；精密比较必须共同匹配辐射修正，不能分别调整两端的归一。

[1052](../research_note_1052.md) 检验弱流比值，[1061](../research_note_1061.md) 把 $G_\mu$ 用作 CKM 抽取输入，[μ 子读口采用](muon_decay_readout_adoption.md) 处理衰变记录；它们未单独完成这条绝对归一匹配。本批接入成熟关系及一个直接质量读口，**新增科学/经验组均为 0**。未建立新的独立参数后验或显著性，亦不宣称全部 P981 高维算符分支共享最小 SM 结果。

## 2. 两层修正属于同一匹配链，不能重复计入

采用自然单位。MuLan 的寿命合同为

$$
 \tau_\mu^{-1}=\frac{G_\mu^2m_\mu^5}{192\pi^3}(1+\Delta q).
 \tag{1}
$$

$\Delta q$ 包含相空间及低能辐射修正。采用该文 v2 的已发表提取：
$\tau_\mu=2196980.3(2.2)\ {\rm ps}$，
$G_\mu=1.1663788(7)\times10^{-5}\ {\rm GeV}^{-2}$。
寿命直接定义的是 μ 衰变专属耦合；跨过程普适性仍是待检的物理合同，不由符号 $G_F$ 赋予。[MuLan，式 (1)、(4)、(5)](https://arxiv.org/pdf/1010.0991v2)

在此低能归一已经定义之后，采用 Awramik 等 v3 的 SM 匹配：

$$
 M_W^2\left(1-\frac{M_W^2}{M_Z^2}\right)
 =\frac{\pi\alpha(0)}{\sqrt2G_\mu}(1+\Delta r),
 \tag{2}
$$


$$
 \Delta r^{(\alpha)}
 =\Delta\alpha-\frac{c_W^2}{s_W^2}\Delta\rho
   +\Delta r_{\rm rem},\qquad
 \Delta\alpha=\Delta\alpha_{\rm lept}+\Delta\alpha_{\rm had}^{(5)}.
 \tag{3}
$$

这里的 $\Delta r$ 按原文**展开**，不能直接改成 $1/(1-\Delta r)$；$\Delta q$ 也不能再加一次。式 (3) 只说明结构，计算不只保留领先 top 项。v3 式 (4) 采用完整电弱两圈及所列高阶 QCD/领先三圈项，并非全部阶数。

式 (2) 的根号正支对应已采用的 $M_W>M_Z/\sqrt2$ 物理支；这项分支选择不是三个数自动导出的唯一性。因 $\Delta r$ 本身依赖 $M_W$，原计算迭代求解。[Awramik 等 v3，式 (1)—(5)](https://arxiv.org/html/hep-ph/0311148v3)

## 3. 质量定义、可复算近似与版本

running-width 分母 $s-M^2+i s\Gamma/M$ 的复极点为

$$
 s_{\rm pole}=\frac{M^2}{1+i\Gamma/M}
 =\bar M^2-i\bar M\bar\Gamma,\quad
 (\bar M,\bar\Gamma)=
 \frac{(M,\Gamma)}{\sqrt{1+(\Gamma/M)^2}}.
 \tag{4}
$$

这是代数转换，不是两个物理质量测量。Awramik 内部使用复极点定义，但最终 $M_W$ 参数化已转换回 running-width；CMS 也使用 running-width。二者不能再转换一次，原文表 1 的内部 $\Delta r$ 数值也不能直接当成任意外部质量约定下的修正。

复算固定使用 v3 式 (6)、(7)、(9)：

$$
\begin{aligned}
M_W={}&M_W^0-c_1H-c_2H^2+c_3H^4+c_4(h-1)
-c_5a+c_6t-c_7t^2\\
&-c_8Ht+c_9ht-c_{10}s+c_{11}z,\\
H={}&\log(M_H/100),\quad h=(M_H/100)^2,\quad
t=(m_t/174.3)^2-1,\\
z={}&M_Z/91.1875-1,\quad
a=\Delta\alpha/0.05907-1,\quad
s=\alpha_s(M_Z)/0.119-1 .
\end{aligned}\tag{5}
$$

质量输入以 GeV 计。全部系数在 [检查程序](electroweak_absolute_normalization/check.py) 明列；本批只核原参考点和 $M_H=125\ {\rm GeV}$、其余保持原参考的诊断点及导数。**125 这一设定不构成现实全输入拟合；不将任一诊断输出减去 CMS 值计算偏差。**

版本与误差边界：

- v3 日期为 2021-11-08，式 (5) 已用 $G_\mu=1.166379\times10^{-5}\ {\rm GeV}^{-2}$；不是更早的 $1.16637$。参数化没有任意 $G_\mu$ 调节项。
- 式 (9) 的 $M_W^0=80.3779\ {\rm GeV}$。表 2 后保留的 $80.3799$ 叙述不作为本批绝对基准。
- 原文 $0.25\ {\rm MeV}$ 是其 $100\le M_H/{\rm GeV}\le1000$、其他输入原联合 $2\sigma$ 域内对所保留计算的近似精度；约 $4\ {\rm MeV}$ 是轻 Higgs 区遗漏高阶估计，既不是严格界，也不是自动独立 Gaussian 标准差。二者不可混为总预测误差。[原文式 (9)、(10) 及相邻讨论](https://arxiv.org/html/hep-ph/0311148v3)

## 4. 一个直接读口：采用 CMS v2 的条件性结果

固定 [CMS arXiv:2412.13872v2](https://arxiv.org/html/2412.13872v2)，采用其 2016 年数据的发表结果：

$$
 M_W^{\rm CMS}=80360.2\pm9.9\ {\rm MeV}.
 \tag{6}
$$

它来自带模拟/校准/理论 nuisance 参数的 $(p_T^\mu,\eta^\mu,q^\mu)$ 模板似然。质量参数在拟合中不受约束；模板中心使用已有实验值不等于强制 SM 预测。然而 J/ψ 标尺、Z 闭合/分辨率、PDF 和产生模型仍参与提取。尤其 LEP $M_Z$ 的闭合误差贡献为 $1.7\ {\rm MeV}$；宽度随质量按 SM 的 $\Gamma_W\propto M_W^3$ 变化并有理论先验，其质量误差贡献小于 $0.2\ {\rm MeV}$。名义结果还采用两种 W 电荷质量相等。

原文称其与所引 EW 拟合相容。本批只采用这一发表判断，**未重建该拟合，未把共享 $M_Z$、宽度或模型输入当作独立，未平方相加分项 impacts 生成新误差，也不拼 CDF/ATLAS 平均**。[CMS §0.1、§0.5—0.5.3、§.7.5、§.7.12—13](https://arxiv.org/html/2412.13872v2)

## 5. 现在签收的共同义务与停止线

同一 SM 分支不能任意给低能 $G_\mu$ 与高能质量各配一套不相容参数。若将来作新的回顾检验，须固定不以被检 $M_W$ 校准的 $\alpha(0),G_\mu,M_Z,M_H,m_t,\Delta\alpha_{\rm had}^{(5)},\alpha_s$ 输入、top 质量方案及共享误差，并运输质量/宽度约定；直接重建质量的 MC 参数不能无误差当作匹配中的 on-shell top 质量。只有该共同合同下预测与读口不能同时满足既定覆盖/误差要求，才能否决该合同。没有共同协方差时不能臆造独立显著性。

本批的小检查只验证式 (1)、(4)、(5) 的单位、代数及局部敏感性；[结果](electroweak_absolute_normalization/results.json) 不含新的现实参数预测或 $p$ 值。原参考 $G_\mu$ 与 MuLan 中心的微小差只以“固定 $\Delta r$”导数诊断，不代替全匹配运输。

这是 M3B/M5/M6 的**共同参数映射**补入；不要求全部装置同体，也不授予完整未知输入/参考后态、总应力、认知独立选择或 ROADMAP 完成。[共同范围](joint_scope_after1062.md) 维持。来源、版本与检查入口见 [sources](electroweak_absolute_normalization/sources.md)。
