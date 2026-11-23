# 第1492轮独立复核　味绝对和、热初态与质量分离的必要界

取 $u=v/\sqrt2$、$v=246.22$ GeV，势中二次项为 $-M_H^2 H^\dagger H/2$。登记每份一圈阈值 $M_i^2(Y^\dagger Y)_{ii}/(8\pi^2)$，乘正运输系数 $[1/2,2]$ 后由两份阈值供给整个弱标量质量参数。改变势中二次项的二倍约定时须连同阈值一起转换。两次实际轻质量平方差来自第1484轮；$n=3$、列 $(0,1,1)$ 与 $(1,3,1)e^{i\eta/2}$ 及 $\cos\eta=-1/2$ 是候选的输入。令 $\tilde m_i=u^2(Y^\dagger Y)_{ii}/M_i$，两列的归一内积模是 $\rho=4/\sqrt{22}$。对于给定轻谱，$\tilde m_a\tilde m_b=11m_2m_3/3$，$(\tilde m_a+\tilde m_b)^2=m_2^2+m_3^2+10m_2m_3$。故两个根对调仍有相同的总强度。即使不考虑任何热史，正立方匹配曲线 $(\sum_i F_i\tilde m_i M_i^3)/(8\pi^2u^2)=M_H^2$ 在固定 $F_i$ 后通常仍留一维，足以给许多重质量，无法由弱尺度和振荡唯一反推。[Littlest Seesaw](https://arxiv.org/abs/1512.07531) · [同类阈值与热史的原研究](https://arxiv.org/abs/1905.12642)

### 允许任意味取向的锐界

以下使用真空、微扰、互相分开的重质量，在SM一圈衰变的约定中取 $x=M_j^2/M_i^2$。记 $h=Y^\dagger Y$、$y_i=Y_{:i}$，令

$$
f(x)=\sqrt{x}\,[1/(1-x)+1-(1+x)\log(1+1/x)],\qquad
g(x)=1/(1-x),\qquad z=f h_{ij}+g h_{ji}.
$$

各味不对称 $\epsilon_{i\alpha}$，在一致的整体符号下，是矩阵 $\mathcal A/(8\pi h_{ii})$ 的对角元，其中 $\mathcal A=[z y_jy_i^\dagger-(z y_jy_i^\dagger)^\dagger]/(2i)$。任何酉基下，对角绝对和都至多为迹范数，选其本征基则达到。二维支撑上 $\det\mathcal A=-|z|^2(h_{ii}h_{jj}-|h_{ij}|^2)/4$，故单个源可达的最大绝对和为

$$
\sup_U\sum_\alpha|\epsilon_{i\alpha}|=
\frac{\sqrt{\{\mathrm{Im}(z h_{ij})\}^2+|z|^2[h_{ii}h_{jj}-|h_{ij}|^2]}}{8\pi h_{ii}}
\le\frac{\rho h_{jj}}{8\pi}(|f|+|g|).
$$

第二界亦直接来自 $\sum_\alpha|Y_{\alpha i}^*Y_{\alpha j}|\le\sqrt{h_{ii}h_{jj}}$。前式允许复杂带电旋转，一个基底未必同时使两个源都达到最大。$g$ 项在迹上抵消、对全味绝对和可有贡献，把非分味 Davidson–Ibarra 界直接用于绝对和会漏掉该项。脚本用不同质量根的实际三行复Yukawa算全部本征值，分别检验等式及原味读数。[分味衰变公式及有限环函数，式120–122](https://scoap3-prod-backend.s3.cern.ch/media/files/19429/10.1016/j.nuclphysb.2017.03.015_a.pdf) · [分味界与求和界之别](https://arxiv.org/abs/hep-ph/0601083)

### 连续分离范围

记 $r=M_\mathrm{heavy}/M_\mathrm{light}>1$，$J(x)=(1+x)\log(1+1/x)-1=\int_0^1(1-y)/(x+y)\,dy$。因而 $xJ(x)\le1/2$，以及 $\sqrt x J(x)\le\int_0^1(1-y)/(2\sqrt y)\,dy=2/3$。自能部分保留完整有理函数，对同一物质的轻、重内线分别给

$$
r(|f(r^2)|+|g(r^2)|)\le\frac{r}{r-1}+\frac12,\qquad
|f(r^{-2})|+|g(r^{-2})|\le\frac{r}{r-1}+\frac23 .
$$

两者共同压到 $M_\mathrm{light}\sum\tilde m_i$。进一步放开每列的实际强度至从轻谱重建上界的两倍，保持列方向；令 $t_{\rm hi}$ 为其和的全盒上界。正阈值和运输下端给出

$$
M_\mathrm{light}\sum_i\tilde m_i
 \le (16\pi^2 u^2 M_H^2)^{1/3}t_{\rm hi}^{2/3}.
$$

保守取 $M_H\le126$ GeV，得到右端 $0.00132013\ {\rm GeV}^2$。每个重Majorana源初始数与熵之比至多 $135\zeta(3)/(4\pi^4g_*)$，$g_*=106.75$，每味净效率绝对值至多一；再乘 $28/79$ 与 $s_0/n_{\gamma0}=7.04$ 得公共转换 $0.00973498$。合计两个源的绝对贡献后，在 $r\ge3$ 全段有 $|\eta_B|\le3.1170\times10^{-11}$。反过来取保守比较目标 $5.8\times10^{-10}$，其余合同不变，则使用这些非共振环函数的解须有 $r\le1.025874$。这是偏宽的上界；在 $r=3$ 的两个中央根，逐源谱的锐界给出的上界之和各约 $3.05$、$2.94$ 乘 $10^{-12}$。这些数均是最多可产生量，尚未解出实际输运。重态极近时需要有限宽度/相干输运计算；2.59%的必要近邻性远不能证明实际宽度已重合。非热初态、其它破荷源、额外熵流、改变列方向或显著有限温度效应也会改变检验。[宇宙学原始测量](https://arxiv.org/abs/1807.06209) · [独立复算](independent_results.json)
