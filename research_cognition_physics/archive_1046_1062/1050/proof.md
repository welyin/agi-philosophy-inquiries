# 1050解析证明：非谱Coulomb材料的自然后态电荷仪器

2026-10-08。这里的物理采用、有限任务与证明范围均为本轮结论的一部分。配套[代码](coulomb_natural_instrument.py)只核轨道积分、矩阵身份及有理常数，不用一个六维代理Hamiltonian代替连续Coulomb演化。

## 1. 父对象与本轮要连接的量

采用三维非相对论电子、固定两个单位正核，原子单位取$\hbar=m_e=e^2=1$。核在$R_\pm=\pm R e_x/2$，电子空间为含自旋的两粒子反对称子空间。**全文H只指电子Hamiltonian**：

$$
 H=\sum_{a=1}^2\left[-\frac12\Delta_a
 -\frac1{|r_a-R_+|}-\frac1{|r_a-R_-|}\right]
 +\frac1{|r_1-r_2|},\qquad X=x_1+x_2.
 \tag{1}
$$

固定核的完整Born–Oppenheimer能源另含$1/R$；它在本轮演化中只加共同相位，但对R变化或能源来源不能省略。固定核支持、电子质量、电荷、空间维数、受控均匀偶极作用、辅助准备和末读均为采用。本轮不从认知原则推出Coulomb作用，也不声称这就是完整P981父对象。

初态空间是两个正常实1s Slater轨道的Löwdin正交化及其全部两电子态：

$$
 \phi_L(r)=\pi^{-1/2}e^{-|r+Re_x/2|},\quad
 \phi_R(r)=\pi^{-1/2}e^{-|r-Re_x/2|},\quad
 S=e^{-R}(1+R+R^2/3),
$$


$$
 G=\begin{pmatrix}1&S\\S&1\end{pmatrix},\quad
 W_1=(\phi_L,\phi_R)G^{-1/2},\quad
 W=\bigwedge^2(W_1\otimes I_{\rm spin}):\mathbb C^6\longrightarrow\mathcal H.
 \tag{2}
$$

这些轨道正常、在$H^2$域内，有Coulomb核尖点，未称其在核处无限光滑。$P=WW^\dagger$**不是**H的谱投影，也不要求H保P不变。

记$Q=(n_R-n_L)/2$、$D=Q^2$。Q谱为$-1,0,1$，D保留全部双占据态而非把六维先裁为四维自旋空间。反演对称性与(2)给

$$
 W^\dagger XW=x=\ell Q,\qquad \ell=R/\sqrt{1-S^2}.
 \tag{3}
$$

任务的理想最终算符是$W_T(I-D)$、$W_TD$，其中$W_T=e^{-iHT}W$。因此任务读取**初始轨道电荷奇偶，物理后态由同一个H自然运送**；不是T时刻未经运输的Q测量，不要求材料冻结在原P。

## 2. 真实Coulomb域、能源与电流界

三维Hardy界$\|\psi/|r-r_0|\|\leq2\|\nabla\psi\|$及动量对Laplacian的无穷小相对界，使有限个Coulomb项相对于自由动能相对界为0。相对坐标同样适用。故(1)在反对称$H^2$域上自伴；这里不是在任意奇异势上任选边界延拓。

对任意单位态，每个电子的两个核吸引满足

$$
 \langle r_{a+}^{-1}+r_{a-}^{-1}\rangle
 \leq4\|\nabla_a\psi\|=4\sqrt{2\langle T_a\rangle},\qquad
 T_a-4\sqrt{2T_a}\geq T_a/2-16.
$$

电子间斥力非负，所以按二次型有

$$
 H+32\geq T_{\rm kin}/2,\qquad
 j=i[H,X]=p_{1x}+p_{2x},\qquad
 j^2\leq4T_{\rm kin}\leq8(H+32).
 \tag{4}
$$

这里j是实际微观电流，所有局部Coulomb势与X对易。非谱P下不能改用$i[W^\dagger HW,W^\dagger XW]$；两者差含P外交叉项。

### 2.1 对全部六维输入的一致初值界

令$r_*=\sqrt{2/(1-S)}$。一轨道有$\|\nabla\phi\|=1$、$\|(x-x_{\rm center})\phi\|=1$，而$\|G^{-1/2}\|=(1-S)^{-1/2}$。两列算子的粗范数界给

$$
 b_X=\|(I-P)XW\|\leq 2r_*.
 \tag{5}
$$

具体而言，将$(\phi_L,\phi_R)$各列中心乘法分出，其像在原轨道空间内；剩余两列各范数1，再正交化。一粒子P外范数至多$r_*$，两粒子和至多$2r_*$。外积限制不增大范数。代码还校准了更紧的全六维泄漏Gram，但后续严格窗只用(5)。

利用每列在其本核下的氢原子恒等式，单电子双核算符的每列范数至多$1/2+2=5/2$，因此两个电子之和范数至多$5r_*$。对电子间项以$r=r_1-r_2$、$\nabla_r=(\nabla_1-\nabla_2)/2$应用Hardy，得到

$$
 \|r_{12}^{-1}\Psi\|\leq\|\nabla_1\Psi\|+\|\nabla_2\Psi\|\leq2r_*,
 \qquad \|HW\|\leq7r_*.
 \tag{6}
$$

这里H仍是电子部分；总BO算符另加$1/R$。令

$$
 E_* =\sup_{\|\psi\|=1}\langle W\psi|(H+32)|W\psi\rangle
 \leq32+7r_*,\qquad J_0=\sqrt{8E_*}.
 \tag{7}
$$

由H的保能演化，$\|j e^{-iHt}W\|\leq J_0$对所有t成立，不要求初态为基态或稳定分子，也不引入删除轨道能隙。

## 3. 实际有限脉冲及传播存在性

辅助qubit准备为$|+\rangle$，其所列裸Hamiltonian为0。采用同一有限脉冲

$$
 H_F(t)=H\otimes I+\beta(t)X\otimes P_1,
 \qquad P_1=|1\rangle\langle1|,
$$


$$
 \beta(t)=\frac{2\pi}{\ell T}\sin^2\frac{\pi t}{T}\quad(0\leq t\leq T),
 \qquad \beta(t)=0\quad\hbox{其余时刻}.
 \tag{8}
$$

这是一个已采用的受控偶极Hamiltonian，不是已构造的有限Maxwell量子场或自治控制器。均匀长度规线性势没有全空间下界，不能把它称作封闭正能源装置；未驱动H的下界(4)保持有效。

传播不是靠将X假定为有界得到。对辅助1分支令$A(t)=\int_0^t\beta(s)ds$，以$\psi=e^{-iA(t)X}\chi$变换，得到

$$
 i\dot\chi=K(t)\chi,\qquad K(t)=H-A(t)j+A(t)^2.
 \tag{9}
$$

最后常数系数为两电子总$\sum q_a^2/(2m_a)=1$。j相对H界为0，$A(t)$有界且光滑，故K(t)有共同$H^2$域，图范数在有限窗一致等价，$K(t)(K(s)-i)^{-1}$范数连续可微。$-iK(t)$及反向时间均生成酉群，满足[Kato 1953原文§1，条件C1—C4与定理3—4](https://www.jstage.jst.go.jp/article/jmath1948/5/2/5_2_208/_pdf/-char/en)的共同域演化条件。由此得到唯一有限时间酉传播，再以(9)运输回长度规。

Slater准备另有有限位置矩。先用有界位置截断在共同域上建立电流／位置积分式，再用(4)和下文驱动能源界作一致控制，移去截断。因此下文对X的Duhamel式具有强向量意义，不是在一般态上随意使用无界算符交换。

[Avron—Herbst 1977原文§II](https://phsites.technion.ac.il/avron/wp-content/uploads/sites/3/2013/05/Commun_math_Phys_52_239%E2%80%94254_1977.pdf)提供常Stark场的规范／传播背景。本轮不将其常场散射结论用于认证这个有限脉冲或整台装置；实际时变共同域核查由(9)承担。

## 4. 关键桥：保全部自然演化，而非强制原轨道不动

在H相互作用图景中$X(t)=e^{iHt}Xe^{-iHt}$，有

$$
 X(t)W-XW=\int_0^t e^{iHs}j e^{-iHs}W\,ds,
$$


$$
 \|X(t)W-Wx\|\leq b_X+tJ_0.
 \tag{10}
$$

选有限逻辑比较传播$V(t)=\exp[-i x A(t)\otimes P_1]$。这是比较工具，未宣称它是实际Coulomb材料Hamiltonian。由Duhamel，作用在任意六维与辅助输入上的差满足

$$
 \boxed{\|U_F(T)(W\otimes I)-(W_T\otimes I)V(T)\|
 \leq b_X L_B+J_0 M_B},
$$


$$
 L_B=\int_0^T|\beta(t)|dt,\qquad
 M_B=\int_0^Tt|\beta(t)|dt.
 \tag{11}
$$

证明只需把(10)作用于$V(t)$旋转后的有限输入，再用全输入算子范数界；真实传播的酉性保持差范数。因此没有遗漏逻辑输入因V而改变的误差项。对(8)的非负对称包络，**恰有**

$$
 L_B=\pi/\ell,\qquad M_B=TL_B/2,
 \qquad \epsilon=\frac\pi\ell(b_X+TJ_0/2).
 \tag{12}
$$

一般同面积脉冲不可免费使用1/2。

因为$e^{-i\pi Q}=I-2D$，对辅助末端X基PVM，理想两个材料算符为

$$
 K_+=W_T(I-D),\qquad K_-=W_TD.
 \tag{13}
$$

这是完整Kraus后态；包括轨道相干、电子关联、真实自然演化和离开P的成分。对任意六维混态及任意被动参考，把初态纯化后应用等距范数界；纯态半迹距离不大于向量差，末读是共同CPTP操作，故保材料、参考、经典结果的联合输出半迹距离至多ε。此为所有输入证明；代码的纠缠参考校准不是用几个样本替代该量词。单个低概率分支归一化的相对误差另需概率下界。

后续若双方继续同一个H自由演化，通道距离不增加；但此时的逻辑嵌入是$W_t$，旧固定Q读口、旧956长期h、旧时变资源均不自动随之运输。

## 5. 同一实际微观来源的独立矩界

半迹距离本身不能控制无界j。对实际驱动态（含辅助／参考），令$Y(t)=\langle H+32\rangle$，(4)与(8)给

$$
 |\dot Y|=|\beta\langle j\otimes P_1\rangle|
 \leq |\beta|\sqrt{8Y},\qquad
 \sqrt{Y(t)}\leq\sqrt{E_*}+\sqrt2 L_B(t),
$$


$$
 \|j\Psi_F(t)\|\leq J_0+4L_B(t).
 \tag{14}
$$

理想自然输出$W_TV(T)$的电流二阶矩范数至多$J_0$，因为V只把输入留在同一个全六维准备空间。对两输出向量用Cauchy，得到

$$
 |\Delta\langle j\rangle|\leq\epsilon\,[2J_0+4L_B].
 \tag{15}
$$

辅助末读和j作用于不同因子，故该式同样控制未归一分支电流及绝对值不大于1的结果标签加权电流，不给小概率条件化的免费精度。它只认证一个实际微观来源的一阶矩；不是完整应力、全部频率响应或动态引力误差。

位置方面，$Z_0=\|XW\|\leq\sqrt{\ell^2+b_X^2}$。实际驱动态有

$$
 \|X\Psi_F(t)\|\leq Z_0+J_0t+4\int_0^tL_B(s)ds.
 \tag{16}
$$

因$\operatorname{TV}(\beta)=4\pi/(\ell T)$，外控工作可由(16)的最大值乘该全变差控制。末端$\beta(T)=0$，所列辅助裸H为0，所以所列模型的辅助PVM与材料H对易，不改变非选择材料能源期望；这不认证测量装置、记录擦除或控制源的总能源闭合。核间能源$1/R$始终保留于BO源账。

### 5.1 代码中实际j的直接字典

原Slater交叉元满足

$$
 \langle L|p_x|R\rangle=iS'(R)
 =-i e^{-R}R(R+1)/3.
$$

Löwdin后为$\kappa\sigma_y/\sqrt{1-S^2}$，再按CAR升为六维$d\Gamma(p_x)$。代码直接积分这个微观p，不用$i[h,x]$替代。

若$x_1=W_1^\dagger x W_1$、$x_2=W_1^\dagger x^2W_1$，则完整两电子泄漏Gram为$d\Gamma(x_2-x_1^2)$；电流二阶压缩为

$$
 W^\dagger j^2W=d\Gamma(p)^2+d\Gamma(p_2-p^2).
 \tag{17}
$$

保留第二项是物理矩与仅压缩矩阵平方的区别。

## 6. 固定有限窗的有理证书

取$R=1000$。S在R>0单调下降，$e^{10}$的正Taylor下界足以证明$S(10)<1/100$，故本轮S更小。精确比较

$$
 r_*<10/7,\quad b_X<20/7,\quad \|HW\|<10,
 \quad E_*<42,\quad J_0<37/2,\quad \ell\geq1000.
$$

同一脉冲取$T=1000/137$，并用$\pi<22/7$：

$$
 \boxed{\epsilon<\frac{74239}{335650}<\frac9{40}}.
 \tag{18}
$$

实际双占据和单占据在同一末端读口的概率差至少

$$
 1-2\epsilon>\frac{93586}{167825}>\frac{11}{20}.
 \tag{19}
$$

因此不是用$T\to0$把任务消成恒等。对双占据与单占据等权相干，理想非选择仪器产生的半迹扰动为1/2；完整实际输出仍由(18)约束。

源界(15)在这一粗证书下约为8.19原子动量，具有单位；它是有限但偏松的一阶源误差，不能与无量纲状态误差0.22118混称为同一精度，也不能宣称精密电流测量已经实现。

在另采用$c\geq137$的原子单位字典时，所选T不短于R/c。**这只排除所选保守时间尺度检查与仪器精度之间的直接矛盾**。并未证明T≥R/c是所有预置共同控制的普遍必要条件，也未把T≈R/c当作严格长波条件；真实Maxwell包络、光锥、量子辅助控制和固定核支持仍未构造。正常Slater态虽有无穷空间尾，却有有限位置／能源矩；本轮没有称R为精确紧支撑直径。

## 7. 去重、剩余输入与反向失败判据

- [956](../../archive_956_989/research_note_956.md)、[960](../../archive_956_989/research_note_960.md)提供六维charge任务；本轮复用$Q^2=D$，不重新设计有限仪器。
- [965](../../archive_956_989/research_note_965.md)、[981](../../archive_956_989/research_note_981.md)的原生模式与来源任务有额外时变／量子场条件，本轮没有自动接入。
- [1048](../research_note_1048.md)是谱P及实际h/x/能隙匹配桥。本轮用具体真实Coulomb准备、自然$W_T$和微观j矩替换谱P条件，取得不同而明确的全输入仪器任务。不能把这写成1048全部假设已消除或旧长期h已匹配。
- [准入](../_admission/material_embedding/selection.md)中的K>0、CAS非谱与交叉来源障碍仍在。自然后态保留H全部作用，没有靠改U/v消掉真实pair transfer。
- 若坚持最终码仍为固定W，(11)必须加入自然演化相对该码的匹配误差，本轮没有界。若换任意包络，必须使用真实$M_B$。若删源矩前提，状态误差不能推出(15)。若要求Maxwell／自治装置完成，(8)不再足够。

本轮关闭的是一项**已采用真实NR Coulomb父作用中的有限完整后态桥**，不是所有物理输入生成、全P981识别或ROADMAP完成。完成后不继续轨道、脉冲和器件精度优化。
