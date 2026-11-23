# 第1415轮独审：125 GeV 标量能避开哪些变味读数？

作者定理通过。逐因子比较四元乘积、秩一投影的算子范数及三中性态正交求和均成立。正文按无序味对 $i<j$ 定义 $S_f$。下文针对单粒子另用有序和 $F_{hf}$，因为左右手有两个耦合。弱尺度的三代 CKM 与 Yukawa 取自第1386轮同一匹配资料；边际宽度盒只作为输入扰动范围。三个非简并质量和实际弱流给味标签以物理意义。[弱尺度资料](https://arxiv.org/html/2510.01312v2) · [作者定义](proof.md) · [独立结果](independent_results.json)

## 1. 由带电流约束真空相位与味因子的乘积

记 $\phi=2\theta$，$p_i=(n_f)_i^2$，$\chi_f=1-\sum_i p_i^2$。对四元乘积的四个因子 $z_j$（含相应复共轭）取 $W^2=\sum_j\prod_{k\ne j}|z_k|^2$。比较 $\prod z_j$ 与其复共轭并依次替换因子，由 Cauchy–Schwarz 得 $|J|\le W\|\operatorname{Im} V\|_F$。这对每个合法重相的 $V$ 都成立。逐行最优重相给

$$
\min_{\alpha}\|\operatorname{Im}(D_\alpha V)\|_F^2
=\frac12\sum_i\left(1-\left|\sum_j V_{ij}^2\right|\right)
=\frac12\sum_i\left(1-\sqrt{1-4p_i(1-p_i)\sin^2\phi}\right)
\le 2\chi_u\sin^2\phi .
$$

第二步用 $VV^T=I+(e^{2i\phi}-1)n_u n_u^T$。逐列应用于 $V^T$ 给 $\chi_d$ 同一界。由此两部门**各自**满足

$$
\chi_f\sin^2\phi\ge L=\frac{J^2}{2W^2}.
$$

中央从九个四元组中选 $W$ 最小者得 $L=5.76187\times10^{-6}$。两边际宽度的整盒用固定的 $(u,c;d,b)$，取 $|V_{ud}|\le1,\ |V_{cb}|\le s_{23}^{+},\ |V_{ub}|\le s_{13}^{+},\ |V_{cd}|\le s_{12}^{+}+s_{23}^{+}s_{13}^{+}$，得 $L\ge4.16367\times10^{-6}$。本界保留 $\phi$，因此能接入标量谱。

## 2. 指定原实标量势后，轻标量的总归一变味强度也有下界

复用原模型的树级实势、软破缺、两重态非零真空和 $h=125$ GeV 为最轻中性本征态。将标量质量平方矩阵记为 $A$，取原文约定 $R^T A R=\operatorname{diag}(m_h^2,m_H^2,m_A^2)$，$\Delta=\max(m_H^2,m_A^2)-m_h^2$。式(67)给

$$
A_{13}=-v^2\lambda_5\sin(2\beta)\sin\phi,\qquad
A_{33}=2v^2\lambda_5\sin^2\theta\ge m_h^2.
$$

$A_{13}$ 的谱分解中减去 $m_h^2 I$，由 Cauchy–Schwarz 得 $|A_{13}|\le\Delta\sqrt{1-R_{1h}^2}$。因 $2v^2\lambda_5\ge m_h^2$，有序和满足

$$
F_{hf}=\sum_{i\ne j}|v g^h_{ij}/m_j|^2
=\frac{4(1-R_{1h}^2)\chi_f}{\sin^2(2\beta)}
\ge\left(\frac{m_h^2}{\Delta}\right)^2\chi_f\sin^2\phi
\ge\left(\frac{m_h^2}{\Delta}\right)^2 L .
$$

严格简并 $\Delta=0$ 会给 $A_{13}=0$，与 $J\ne0$ 及这些真空条件不相容。两宽盒内，若独立确定最重中性态至多 600 GeV，得到每部门 $F_{hf}\ge8.5714\times10^{-9}$。若只采用原文树级微扰幺正条件，式(117)在 $\Lambda=16\pi$ 及正 $\lambda_1,\lambda_2,\lambda_5$ 下给 $\max(\lambda_1,\lambda_2)<\Lambda/6$ 与 $\lambda_5<\Lambda/2$。式(70)的迹给 $\Delta<4\Lambda v^2/3-3m_h^2$，相应 $F_{hf}>6.3021\times10^{-11}$。两数的适用范围取决于所给标量势、质量次序及匹配精度。改变势或放宽能标范围时必须重算。[原文式(67)、(70)及附录B](https://arxiv.org/pdf/1808.00493#page=14)

## 3. 实际125态可同时在两个顶味变通道中缺席

用完整实测中央 CKM 做更换标签后的 PDG 分解。先将行排作 $(t,u,c)$、列作 $(b,d,s)$，从模和有符号 $J$ 重建角度。除去允许的右重相，矩阵可写作 $R_{23}\operatorname{diag}(1,1,e^{i\phi})R_{13}R_{12}$，其中 $\phi=2.73201425$。恢复物理次序后上型方向为 $(0.22528983,0.97429179,0)$。这样 $t\leftrightarrow u$ 和 $t\leftrightarrow c$ 在全部三个中性态上严格为零，$u\leftrightarrow c$ 承担上型变味。原 CKM 九项模的最大复算偏差为 $1.2\times10^{-16}$ 以内，$J=+3.0875323\times10^{-5}$。选择左右实种子配给实际夸克奇异谱不改变这项构造。

同相位可用原文的势取得。取 $\tan\beta=1$，初始 $\lambda_1=\lambda_2=0.8,\ \lambda_{345}=-0.55,\ \lambda_5=0.9,\ \lambda_4=0$，其中 $\lambda_{345}=\lambda_3+\lambda_4-\lambda_5$。用共同正系数缩放全部四次耦合以使最轻质量恰为125 GeV。对 $X_i=\Phi_i^\dagger\Phi_i,\ Z=\Phi_1^\dagger\Phi_2$，选择实二次项系数使势相对其最小值为

$$
(X-X_0)^T
\begin{pmatrix}\lambda_1&\lambda_{345}\\ \lambda_{345}&\lambda_2\end{pmatrix}
(X-X_0)
+4\lambda_5(\operatorname{Re} Z-x_0)^2
+2(\lambda_5-\lambda_4)(X_1X_2-|Z|^2).
$$

$X_{i0}=v_i^2/2,\ x_0=v_1v_2\cos\theta/2$。三项非负，两中性真空的相对相位 $\pm\theta$ 使其全取零。树级中性质量为 $(125,297.65,344.84)$ GeV，带电态为 $243.03$ GeV，$R_{1h}^2=0.94513$，最大绝对散射本征值 $6.712<16\pi$。这样对标量混合的算例同步满足势、本征谱、真空和带电流。

ATLAS 对 $t\to hu,hc$ 的组合观测95%上限分别为 $2.6\times10^{-4}$ 和 $3.4\times10^{-4}$，各自采用另一耦合关闭的假设。上面的构造说明这两种搜索即使同时无信号，也容许此层检验中的 CKM 和非零归一总强度。该点的 $u\leftrightarrow c$ 对 $D^0$ 混合有贡献，全局可行性仍需检查介子约束、不同标量交换中的相消、Higgs率和直接搜索。后续原作者也专门考察了顶方向为零时的中性介子问题。即使物理质量有限，可探测率仍依赖产生和衰变、相空间、宽度及实验灵敏度。[ATLAS 检索](https://arxiv.org/abs/2404.02123) · [后续论文附录C](https://arxiv.org/pdf/2105.14054#page=17)

复算使用同目录的[独立程序](independent_check.py)，默认模式核对保存的结果。独立程序只读第1386轮的味输入和第1404轮的 $G_F$，未导入作者程序。
