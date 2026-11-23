# 第1449轮　有符号质量、真空对齐与一个强子留出

## 1　同一紫外基底的两个 CP 问题

两个三阶 Hermitian 纹理的重本征值取正，全部质量非零。对每部门 $\operatorname{sgn}\det M=\operatorname{sgn}(m_1/m_2)$。八组中央解令 $\sigma_u=\operatorname{sgn}(m_u/m_c)$、$\sigma_d=\operatorname{sgn}(m_d/m_s)$，便有

$$
\bar\theta=\theta_{\rm UV}+\arg\det(M_uM_d)\pmod{2\pi},\qquad
\exp[i(\bar\theta-\theta_{\rm UV})]=\sigma_u\sigma_d .
$$

在指定手性基底取 $\theta_{\rm UV}=0$、纯矢量 QCD 输运、无另增复阈值时，四组在零支、四组在 $\pi$ 支。八组的弱 Jarlskog 不变量都非零。复共轭翻转 $J$，保留两个行列式。把重味变为正实再积掉时，其原有相位累加进低能拓扑项。右手改号若同时计入测度反常便保持 $\bar\theta$；只改有符号质量而固定原裸角会改变物理理论。任意许可的裸角平移 $\pi$ 足以对调两组。Hermitian 基底的存在本身不给裸角。[第1445轮](../1445/proof.md) · [原纹理](https://arxiv.org/abs/2501.18508) · [第1364轮](../../archive_1364_1409/1364/proof.md) · [第1389轮](../../archive_1364_1409/1389/proof.md)

$\bar\theta=\pi$ 在作用量层面具有 CP 对称。真空既可能保持也可能自发破坏这一对称。零附近的中子电偶极矩斜率不能外推至此。以下检验绕开在 $\theta=0$ 下用强子谱拟合出来的轻夸克质量。[Aoki 与 Creutz 的有符号质量相图](https://arxiv.org/abs/1402.1837)

## 2　守 CP 的真空由三条质量不等式筛选

采用三味八重态的最低阶手征作用量、无电荷凝聚的对角真空，单独登记 $\eta'$ 消去、强手征高阶和电磁高阶修正。取 $h_i=B_0|m_i|>0$，真空极大化 $\sum h_i\cos\phi_i$，约束 $\sum\phi_i=\bar\theta\bmod 2\pi$。设 $X_i=h_i\cos\phi_i$。驻定性要求各 $h_i\sin\phi_i=L$ 共用。守 CP 的 $\pi$ 支满足 $L=0$ 且奇数个 $X_i$ 为负。[三味手征展开](https://doi.org/10.1016/0550-3213(85)90392-5)

令 $A=X_u+X_d,C=X_u+X_s,N_K=X_d+X_s$。最低阶共同的电磁位移 $\Delta$ 加在带电 $\pi,K$ 的质量平方上。中性质量平方矩阵按 $\pi_3,\eta_8$ 排列为

$$
N=\begin{pmatrix}A&(X_u-X_d)/\sqrt3\\(X_u-X_d)/\sqrt3&(A+4X_s)/3\end{pmatrix}.
$$

其小本征值 $p$ 满足 $\zeta=A-p\ge0$。保留中性混合，对实测质量平方构造

$$
I_u=2m_{\pi^0}^2-m_{\pi^+}^2+m_{K^+}^2-m_{K^0}^2,\quad
I_d=m_{\pi^+}^2+m_{K^0}^2-m_{K^+}^2,\quad
I_s=m_{K^+}^2+m_{K^0}^2-m_{\pi^+}^2 .
$$

最低阶分别为 $2X_u-2\zeta,2X_d,2X_s$，任何 $\Delta$ 都消去。实测宽工作盒按系数逐端传播，下界分别为 $0.01233719,0.02290191,0.47130628\ {\rm GeV}^2$。若每个组合相对该最低阶字典的总偏移绝对值小于其对应下界，则 $X_u,X_d,X_s>0$。程序用共同充分预算 $0.010\ {\rm GeV}^2$，未给未计入的电磁低能常数做拟合。质量输入为强子极点及实验运动学，带电粒子/反粒子同质量按 CPT，$K_S,K_L$ 小劈裂在宽盒内。[带电π原子光谱](https://arxiv.org/abs/1605.03300) · [π质量差](https://doi.org/10.1103/PhysRevLett.56.1043) · [带电K光谱](https://doi.org/10.1103/PhysRevLett.60.186) · [KLOE 的中性 K 与 η](https://arxiv.org/abs/0710.5892)

## 3　同样的谱可以来自真正的破缺真空

若 $L\ne0$ 且 $\bar\theta=\pi$，$\operatorname{Im}\prod_i(X_i+iL)=0$ 给出

$$
L^2=X_uX_d+X_dX_s+X_sX_u=\tfrac34\det N.
$$

给任意三个正 $X_i$，取 $L$ 为此正根，令 $h_i=\sqrt{X_i^2+L^2}$。三相角 $\arg(X_i+iL)$ 的和为 $\pi$，因为乘积的实部为 $-(X_u+X_d)(X_u+X_s)(X_d+X_s)$。这些新正质量在 $\pi$ 处与零角的 $h_i=X_i$ 给出完全相同的带电质量及中性质量矩阵。

该构造取到了所述真空类别的全域极值。固定 $c=\cos\phi_u$ 后，将其余两项的线性正弦、余弦在相对相位上极大化，得到 $W(c)=h_uc+\sqrt{h_d^2+h_s^2-2h_dh_sc}$。所构造的相角在内部取该极大（固定相位和时两项的二阶导为 $-X_d-X_s<0$），且 $W'(c)=0$。严格凹性保证全域最大。相反号相位给另一个相反取向真空。例取 $X=(.0065,.012,.232)\ {\rm GeV}^2$，得到 $\pi$ 支 $h=(.066425,.067186,.241234)\ {\rm GeV}^2$。两套最低阶中性质量都是 $(.13589,.56172)$ GeV。例子验证可辨性，未拟合实际介子谱全套数据。

## 4　独立的双π留出

取八重态矩阵指数 $\exp(i\sqrt2\Phi/F)$，动能正则归一。质量作用量三次项的绝对系数为 $\sqrt2|L|/(3F)$ 乘 $\operatorname{tr}\Phi^3$。迹中 $\eta_8\pi^+\pi^-$ 的系数为 $\sqrt6$，$\pi_3\pi^+\pi^-$ 的系数为零。首阶电磁项 $\operatorname{tr}(QUQU^\dagger)$ 对对角真空相位不敏感，且对 $\Phi\mapsto-\Phi$ 为偶。若重中性本征态在八重态方向的投影为 $R_{\eta8}$，得到

$$
|\mathcal A_{\rm LO}|=\frac{2|L|\,|R_{\eta8}|}{\sqrt3 F}
=\frac{m_{\pi^0,\rm LO}m_{\eta,\rm LO}|R_{\eta8}|}{F}.
$$

相反取向的均匀真空有同一衰变率，普通宏观畴的非相干平均保留该率。理想有限体积、精确 CP 猫态的选择定则需独立处理。实际双π实验使用可作为局域准粒子的窄 $\eta$、通常的辐射标记与接受度。质量工作盒和该过程是不同读数。

独立弱流可约束 $F$。左流正比于 $(\partial U)U^\dagger$，在线性阶仅将带电场乘以单位相位。对角对齐不改纯左流的模。采用 $f_{\pi^+}=\sqrt2 F$ 归一，在标准 $G_\mu$ 标尺、$|V_{ud}|\ge.90$、净带电流幅保留至少 .80 时，以 $\Gamma_{\mu\nu}\le\hbar/\tau_{\pi^+}$ 得 $F\le.127$ GeV。程序用独立的宽条件 $F\le.20$ GeV。采用单调性以低端 $m_{\pi^+}$ 和 $\tau_{\pi^+}\ge25.9$ ns 放宽界。[直接π寿命](https://doi.org/10.1103/PhysRevD.52.4855) · [弱流标尺的适用范围](../../archive_1364_1409/1400/proof.md)

为允许有限精度，登记 $m_{\pi^0,\rm LO}\ge.120$ GeV、$m_{\eta,\rm LO}\ge.45$ GeV、$|R_{\eta8}|\ge.50$，得到 $|\mathcal A_{\rm LO}|\ge.135$ GeV。这里的质量与投影指共同最低阶模型；辛格态混合、实际 $\eta$ 的外腿归一、终态相移及所有干涉可进入总幅的匹配。对总幅写出需要的量，不用自然性估计替代误差控制。

KLOE 一份 $1.61\ \mathrm{fb}^{-1}$ 样本给 $\mathcal B_{\pi^+\pi^-}<4.9\times10^{-6}$（90%）。保守采用 $\Gamma_{\gamma\gamma}<600$ eV、$\mathcal B_{\gamma\gamma}>.30$，来自 KLOE $520\pm20\pm13$ eV、CLEO $(38.45\pm.40\pm.36)\%$，从而 $\Gamma_\eta<2.0$ keV。这些广义边界是工作范围，未合成置信水平。CLEO 测量按五个主要通道归一；若有占优势的未标记新道须重新标定。两体运动学给

$$
|\mathcal A_{\rm phys}|^2
=\frac{16\pi m_\eta\Gamma_\eta\mathcal B_{\pi^+\pi^-}}{\sqrt{1-4m_{\pi^+}^2/m_\eta^2}}.
$$

实测上界对应 $1.7723\times10^{-5}$ GeV。若 $|\mathcal A_{\rm phys}|/|\mathcal A_{\rm LO}|\ge.001$，破缺 $\pi$ 支预言 $\mathcal B_{\pi^+\pi^-}\ge2.8432\times10^{-4}$。要保住该支必须至少放松一项条件，例如容许总幅保留比例低于 $1.313\times10^{-4}$。并无独立定理给这个比例的下界。上述两体式只在实测窄共振的工作范围内反推幅，未把首阶大耦合外推成超过一的分支比。[KLOE 双π检验](https://arxiv.org/abs/2006.14710) · [KLOE 双光子部分宽度](https://arxiv.org/abs/1211.1845) · [CLEO 分支比](https://arxiv.org/abs/0707.1601)

## 5　结论的逻辑范围

在列明的匹配范围，守 CP 与自发破缺的 $\pi$ 真空分别与两组留出冲突，八个局部纹理解中筛走四组。1445保留的中间正、$m_d/m_s<0$ 两个弱支中须选 $m_u/m_c<0$。弱共轭分歧仍由弱测量判。更一般的裸角、额外 UV 相位、轴子动力学、不同真空有效理论或可抵消的奇幅可解除排除。强作用本身未确定裸角，也未证明 Hermitian 纹理的量子稳定性；小变形满足 $\delta\arg\det M=\operatorname{ImTr}(M^{-1}\delta M)$。三个组合的误差和 $\eta$ 的幅传递若要成为对精确 QCD 的无条件结论，须有在 $\bar\theta=\pi$ 适用的独立控制。
