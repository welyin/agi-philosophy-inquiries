# 强角的背景依赖与可测的 Higgs 味响应

## 共同输入

用同一径向背景 $r=v+h$ 比较一次质量测量与一次 Higgs 顶角，场在匹配尺度动能正则。设 $x=(r/v)^2$，每个带色部门有非奇异复矩阵 $Y_f$ 和无量纲 $A_f$，

$$
M_f(r)=rY_f(I+xA_f)/\sqrt2 .
$$

这可来自 $H^\dagger H\,\bar Q H f$ 及上型的共轭双态。首先将右式作为给定的树级质量函数；只保留一般 SMEFT 的 $1/\Lambda^2$ 精度时，涉及多个 $A$ 的乘积没有独立精度。以下不预定认知原则会给这个函数。匹配的常强角及 $\arg\det(Y_uY_d)$ 之和取零，其他带色态、轴子、$H^\dagger H\,G\widetilde G$ 及阈值修正单列。$r=0$ 的夸克无质量，只用 $r\to0^+$ 定多项式常项，未将无质量点的强角当成可测量。

令 $U_L^\dagger M U_R=D>0$，$G=U_L^\dagger M' U_R$，$\kappa_i+i\eta_i=vG_{ii}/m_i$。定义 $P(x)=\prod_{f=u,d}\det(I+xA_f)$。在行列式非零且连续选相位的分支，

$$
\bar\theta(v)=\arg P(1),\qquad
v\partial_v\bar\theta\big|_{\text{固定作用系数}}
=\operatorname{Im}\sum_f\operatorname{tr}(M_f^{-1}v M_f')
=\sum_{\rm quarks}\eta_i .
$$

这里的导数是沿 $r$ 改变背景再取 $r=v$，将 $A_f$ 中的参考匹配尺度固定。恒等式允许左右味基不同。$\tau$ 没有色荷，它的质量相位不能补偿色反常。单个非零复本征值也不能直接辨认为顶味，下面的单味结论另加实际质量方向对齐条件。

## 单个带色味、一个真空

只给对角顶味 $z$。顶质量的相对相位为 $\arg(1+z)$，重相后 $\kappa_t+i\eta_t=3-2/(1+z)$，故 $\bar\theta=\operatorname{atan2}(\eta_t,3-\kappa_t)$，相位按 $2\pi$ 计。小强角支有 $3-\kappa_t>0$。按第1526轮的 PSI 上限与两种格点斜率，各自将两项格点误差线加取两宽，给出 $|\eta_t|/(3-\kappa_t)\le\tan\theta_{\max}$。在 $\kappa_t=1$ 的例子中，上界分别为 $6.21\times10^{-10}$ 和 $1.39\times10^{-10}$。这使用近零强角处的线性中子响应，并要求其余中子 EDM、常强角和阈值对相消的贡献可忽略。若其余中子 EDM 有幅度预算 $d_*$，则用 $(d_n^{\max}+d_*)/|d'_n(0)|$ 代替原角上限。独立常强角可以任意抵消此相位。若顶与 τ 共享同一复 $z$，这项检验也约束 τ 赝标量，无须输入电子 Yukawa。[中子测量](https://arxiv.org/abs/2001.11966) · [格点斜率一](https://arxiv.org/abs/2301.04331) · [格点斜率二](https://arxiv.org/abs/2608.22587)

## 持续保护与有限条件

$P(0)=1$。其虚部在任一非空开区间恒为零，当且仅当 $P$ 所有系数都为实数。由于 $\det(I+xA)=\prod_j(1+xz_j)$，这等价于全部非实本征值按代数重数成共轭对。此结论对非正规矩阵及 Jordan 块成立。回到原点的无零连通区间上实 $P$ 保持正号；其它区间须另查符号。相位随背景不变于是强迫 $\sum_i\eta_i=0$。在一般的一阶有效场论，只能在该阶约束 $\operatorname{Im}\sum_f\operatorname{tr} A_f$；要在精度内辨认共轭配对，须有完整质量函数或对更高维算符另设独立预算。

若带色变形只在已辨认的顶、底两个对角味中，强保护且有非零赝标量迫使 $z_b=\bar z_t$。这个判据甚至有有限版本。记 $P=1+(z_t+z_b)x+z_tz_bx^2$，取 $u=\operatorname{Im}P(1)$、$w=\operatorname{Im}P'(1)$。两个虚系数恰是 $2u-w,w-u$。因此规范匹配常项确定后，一个真空处的零相位和零背景相位导数已经足以核明整个二次多项式。有限误差也可传递。令 $\theta=\arg P(1)$，$R+iS/2=P'(1)/P(1)$，则 $u=|P(1)|\sin\theta$，$w=|P(1)|(R\sin\theta+S\cos\theta/2)$。已知 $p=\operatorname{Im}z_t\ne0$ 时，$\operatorname{Im}(z_b-\bar z_t)=2u-w$，$\operatorname{Re}(z_b-\bar z_t)=[w-u-\operatorname{Re}(z_t)(2u-w)]/p$。单次中子测量没有给 $S$；探测在壳 Higgs 转为这项零动量导数还要计算各味传播的形式因子。

只使一个真空强角为零的反例取 $z_t=it,\ z_b=k(1-it)-1$、$k>0$。$P(1)=k(1+t^2)>0$，而径向角导数为 $2t(1-1/k)/(1+t^2)$。在单顶质量函数另容许 $1+i a x(x-1)$ 时，同样能在 $x=1$ 有零质量相位和非零响应。这登记了高阶输入如何改变结论。

## 独立的底、τ 衰变比

在上述两带色对角味的保护支再加第1532轮检验的共同**完整复变形** $z_\tau=z_t$。于同一匹配面得到 $\kappa_b=\kappa_\tau,\ \eta_b=-\eta_\tau$。这是额外的跨部门来源条件。可让 $Y_d$ 承载原有 CKM 酉框架并使 $A_d$ 在右质量基对齐，因而弱混合可保留为独立输入。树级归一化的衰变宽比为

$$
{\Gamma_{bb}/\Gamma_{\tau\tau}\over(\Gamma_{bb}/\Gamma_{\tau\tau})_{\rm SM}}
={\kappa_\tau^2+\eta_\tau^2/\beta_b^2\over
  \kappa_\tau^2+\eta_\tau^2/\beta_\tau^2},\qquad
\beta_f^2=1-4m_f^2/m_h^2 .
$$

对非零信号，这是 1 与 $\beta_\tau^2/\beta_b^2$ 之间的加权平均。用第1435轮的 ATLAS 表7 同一浮动产生率拟合，$ZZ^*$ 和总宽抵消。脚本逐项重算中央双比和边际盒，未把两个衰变当成独立似然。$m_b$ 用 $[2.5,5]$ GeV 作树级相空间的宽工作区；QCD、运行、其他顶角和分析模板对**双比**的总转移系数另容许 0.90—1.10。这是条件预算，观测精度与预算都大于有限质量效应。当前读数容许该预测。第1532轮的电子偶极约束仍须在同一实物模型上共同满足。[ATLAS 表7](https://arxiv.org/pdf/1909.02845#page=27) · [标量、赝标量夸克宽度的 QCD 修正](https://arxiv.org/abs/hep-ph/9406431)
