# 1146：同一正常运动材料的钟速、标量操作与量子来源

2026-10-12。依[准入](../_shared/round1146_admission.md)。主线完整证明；波包估计详见[两谱包界](packet_bound.md)，正常CAR与来源的独立资格见[来源审查](source_car_review.md)。正式验收以本轮独立审阅及 acceptance.json 为准。

## 1. 要补的共同接口及明确材料输入

941—943已给内部质量、局部Dirac材料和标量来源投影；953已给世界线钟/源共同性；1129已给正常量子来源和全CAR交接。本轮不重证这些一般工具。净问题是：在同一正常运动材料、有限空间模式、同一真实量子来源上，能否同时看到内部钟速与标量操作权重，并保完整未知输入、全部负谱余部及来源反作用？

取自然单位，固定

$$
T=1,\quad t_c=1/4,\quad m=10^{10},\quad
\delta=10^{-3},\quad p_0=3m/4,\quad
v=3/5,\quad r=4/5,\quad E_0=5m/4,
\quad\epsilon=10^{-16},\quad A=10^{100},\quad\lambda=\pi/2.
$$

自由材料是一份两味Dirac场，其Cauchy空间为 $L^2(\mathbb R^3,\mathbb C^4)\otimes\mathbb C^2$，

$$
h_D=\boldsymbol\alpha\cdot p+\beta(mI+\delta P_1),
\quad P_1=|1\rangle\langle1|,
\quad S=1_{h_D<0}.
$$

使用此新异质量自由海的正常粒子/空穴Fock表示及正能源 $H_D$；不声称它与旧 $m=1$ 海有全空间酉等价。完整材料没有删去负频或只保正能一粒子。标准Dirac量子化、标量与矢量耦合的类别区分直接复用[942](../../archive_935_955/research_note_942.md)、[943](../../archive_935_955/research_note_943.md)及[Tong §5.2、5.6](https://davidtong.org/teaching/quantum-field-theory/qfthtml/S5)。这些材料、质量、耦合和准备是新增物理合同，并未从认知公理唯一选出。

来源为1129相同的独立无质量KG正常Fock场。核 $\rho_\epsilon$、卷积 $K_\epsilon$、平台 $\eta=1$ 于 $B_4$、支撑 $B_5$，及源准备平台 $\beta_0=1$ 于 $B_7$、支撑 $B_8$，逐字复用[1129 §1](../1129/proof.md)，只改变明确的有限分辨率和源尺度。脉冲函数 $F$ 支撑 $[1/4,3/4]$，非负、积分1、$\|F\|_\infty\le16$、$\|F'\|_\infty\le512$。

## 2. 同一正常来源下的两份真实H

令 $D_s=\beta$、$D_v=I_4$，定义有限分辨率双点标量核及对照

$$
J_{a,\epsilon}(x)=\sum_{s,t=1}^4(D_a)_{st}
\big[c_0^\dagger(\rho_{\epsilon,x}e_s)c_1(\rho_{\epsilon,x}e_t)
+c_1^\dagger(\rho_{\epsilon,x}e_s)c_0(\rho_{\epsilon,x}e_t)\big],
\qquad a=s,v.
$$

可取标准Hermitian对角 $\beta$；写成一般Hermitian矩阵时第二项按共轭指标配对。两味真空交叉收缩为0，故没有隐藏的跨味正常序常数。CAR及有限矩阵对角化给保守统一界 $\|J_{a,\epsilon}(x)\|\le8q_\epsilon$，$q_\epsilon\le8\epsilon^{-3}$。实际自治过程为

$$
H_a=H_D+H_\phi+V_a,\qquad
V_a=\lambda\int d^3x\,\eta(x)
F(\phi_\epsilon(x)/A)J_{a,\epsilon}(x).
\tag{1}
$$

$V_a$是有界自伴强积分，

$$
C_*=8\cdot8\cdot524\epsilon^{-3}=3.3536\times10^{52},
\qquad \|V_a\|\le16\lambda C_*.
$$

因此两份 $H_a$ 均在 $D(H_D+H_\phi)$ 自伴、半有界。此处 $a=s$ 使用真实 $\beta\sigma_x$ 标量旋量结构，$a=v$ 使用时间分量矢量结构；有限空间涂抹和固定部署已明列，不将其改称精确点态Yukawa或精确Lorentz不变的完整理论。

源准备是自由Cauchy均值 $\phi_{cl}(0)=0$、$\pi_{cl}(0)=A\beta_0$ 的正常相干向量 $|a\rangle$。在声明空间/时间管内其自由平滑均值为 $At$，方差平方根 $v_\epsilon<3/\epsilon$。源场能

$$
E_\phi(0)\le1073 A^2=1.073\times10^{203}
$$

巨大而有限，不能省略。来源仍为量子态并随同一 $H_a$ 变化；均值仅用于比较估计。该源在本实验室部署的均值是 $At$，并非1144世界管源的完整boost，不能据此领取1145全部来源配方运输。

## 3. 实际紧支运动模式与两谱包界

取1129归一实包络 $b(x)=\prod_jq(x_j)/\|q\|_2^3$，$q(s)=e^{-1/(1-s^2)}1_{|s|<1}$。其支撑位于 $B_2$，$\kappa_1=\|\nabla b\|<20$、$\kappa_2=\|\Delta b\|<1000$，导数解析界见[包证明](packet_bound.md)。取 $u_0$ 为 $\alpha_zp_0+\beta m$ 的单位正本征旋量。两个实际输入模式及移动目标模式为

$$
f_j(x)=e^{ip_0z}b(x)u_0\otimes e_j,
\qquad g_j(t,x)=e^{i(p_0z-E_0t)}b(x-vt e_z)u_0\otimes e_j.
\tag{2}
$$

所有模式归一、紧支且正交；不先作正能投影。源均值比较过程的一粒子生成元为

$$
h_{a,\epsilon}(t)=h_D+\lambda F(t)K_\epsilon\eta K_\epsilon\otimes D_a\sigma_x.
$$

两个有效内部酉都取 $R_a(0)=I$，由同一Dirac过程产生

$$
i\dot R_s=r(\delta P_1+\lambda F(t)\sigma_x)R_s,
\qquad
i\dot R_v=(r\delta P_1+\lambda F(t)\sigma_x)R_v.
\tag{3}
$$

钟质量分裂与标量操作共用 $r$，矢量操作则有另一权重；并未事后改钟刻度或把交互系数拟合为 $r$。

包证明在完整正负谱上逐动量用有限能隙估计，再以 $b\in H^2$ 控制全部动量尾，给出

$$
\left\|U_{a,\epsilon}(t)f_j-\sum_k g_k(t)(R_a(t))_{kj}\right\|
\le d_a(t).
\tag{4}
$$

统一有限窗可采用以下保守上界。置 $L_t=\int_0^t\|W_a(s)\|ds$，其中 $W_s=\beta\otimes(\delta P_1+\lambda F\sigma_x)$、$W_v=\beta\otimes\delta P_1+I_4\otimes\lambda F\sigma_x$，并令

$$
D_{pkt}(t)=\frac1m\left[
\frac{\|W_a(0)\|+\|W_a(t)\|+\int_0^t\|\dot W_a\|+3\int_0^t\|W_a\|^2}{2}
+(1+L_t)\kappa_1+\frac{t\kappa_2}{2}\right].
$$

对矢量支路这里只保更粗的同一上界；其正谱压缩严格为 $r(p)\delta P_1+\lambda F\sigma_x$，不把矢量项错误乘 $r(p)$。有限传播、平台和卷积的Duhamel估计另给

$$
d_a(t)\le D_{pkt}(t)+
2\lambda\epsilon(p_0+\kappa_1)\int_0^t F(s)ds.
\tag{5}
$$

真实齐次Dirac模式支撑于 $B_{2+t}$；平台覆盖它与卷积邻域，因此该核比较不依赖未证的实际轨道。这里 $\epsilon$ 没有取极限。本脉冲先增后减，故 $\int|F\prime|\le32$；代码以该更紧的总变差界代入上式，用有理上界 $\lambda<11/7$ 验 $d_a(1)<2.46\times10^{-6}$。在 $t_c=1/4$ 前 $F=0$，核项严格为0，$d_a(t_c)<1.451\times10^{-8}$；不以末端粗界覆盖较小钟信号。

## 4. 正常空码、完整CAR与量子源后态

令 $s_j=\langle f_j,Sf_j\rangle$。两味海没有跨味收缩，而

$$
s_j\le\frac{\kappa_1^2+\delta_j^2}{(m+\delta_j+E_0)^2}<8\times10^{-19},
\quad\delta_0=0,\ \delta_1=\delta,
\quad Z=(1-s_0)(1-s_1)>0.
$$

该界由 $P_j^-(h_j-E_0)=-(\omega_j+E_0)P_j^-$ 及Clifford关系得到，没有删除小负频。取

$$
e_\emptyset=(1-n_{f_0})(1-n_{f_1}),\quad
\Omega_0=e_\emptyset\Omega/\sqrt Z,\qquad
\mathcal J|j\rangle=c^\dagger(f_j)\Omega_0.
\tag{6}
$$

这是正常等距编码，全部输入与任意被动参考 $R$ 允许。准备权限明示，$Z$并不是已经制造了准备装置；尤其不误用1129旧 $Z=1/4$。

均值比较的完整CAR传播 $\widehat U_{a,\epsilon}$ 由同一有界材料生成元确定。根据(4)，先将第一目标列对齐，再在其正交补内对齐第二列，得到至多四个模式上的偶CAR酉 $\mathcal C_a(t)$，

$$
\|\mathcal C_a-I\|\le3d_a(t),\qquad
\mathcal C_a\widehat U_{a,\epsilon}\mathcal J|j\rangle
=\sum_k (R_a)_{kj}c^\dagger(g_k)\Xi_a(t),
\tag{7}
$$

其中 $\Xi_a=\mathcal C_a\widehat U_{a,\epsilon}\Omega_0$ 满足 $c(g_k)\Xi_a=0$。第一校正不超过 $d$，第二列残差不超过 $2d$，完整CAR酉范数相应相加。校正只是比较工具，不是隐藏的物理操作。实际余部 $\Xi_a$ 对未知输入独立，但随过程改变；不同 $a$ 的余部未被强制认作同一态。

源Lipschitz和有界扰动Duhamel逐字复用1129，但用新常数：对全部材料向量及任意参考一致有

$$
\left\|\left(e^{-itH_a}-\widehat U_{a,\epsilon}(t)\otimes e^{-itH_\phi}\right)
(I\otimes|a\rangle)\right\|
\le e_s(t):=\lambda C_*512\frac{3t}{\epsilon A}
<8.1\times10^{-29}\,t.
\tag{8}
$$

于是实际联合输出与(7)的理想内部操作、实际固定余部及自由比较源的等距距离不超过

$$
\varepsilon_a(t)=3d_a(t)+e_s(t),\qquad
\varepsilon_a(1)<7.38\times10^{-6},\quad
\varepsilon_a(t_c)<4.353\times10^{-8}.
\tag{9}
$$

该界包含完整CAR、源相位和后态，不只控制一粒子概率。任意 $Q/R$ 的半迹距离及其后同一CP仪器/继续都受相同界控制；并未声称比较源就是实际不受力的源。

## 5. 两种实际区域读数

终末模式均在有限球内；其单占据偶矩阵单位 $c^\dagger(g_i)c(g_j)$ 构成实际有限CAR逻辑代数。允许该代数的有界投影仪器和失败支路，是本轮明示终末权限；不据效果存在声称探测器供能/自主生成已经完成。全部仪器后态可由(9)经同一个CP映射控制，不后选择丢弃失败或场余部。

**早期实际钟。** 在 $t_c$ 的单占据码内测 $E_{-Y}=(I-\sigma_y)/2$，码外另给失败标签。输入 $|+\rangle$ 时，正确目标概率为

$$
p_c=\frac{1+\sin(r\delta t_c)}2.
$$

错误实验室钟字典给 $p_{lab}=[1+\sin(\delta t_c)]/2$。利用 $\cos x\ge1-x^2/2$，

$$
p_{lab}-p_{actual}\ge
\frac{(1-r)\delta t_c}{2}\left(1-\frac{(\delta t_c)^2}{2}\right)
-\varepsilon_s(t_c)>2.4\times10^{-5}.
\tag{10}
$$

这是同一移动材料的真实局部模式读数。没有把实际场态先投成平面波，也不靠预先赋给钟的重标度证明。

**末端源型。** 在 $t=1$ 取输入 $|0\rangle$ 和同一裸 $|1\rangle$ 码投影。忽略小 $\delta$ 的两个比较操作分别转过 $r\lambda$ 与 $\lambda$；恢复 $\delta$ 的每支内部酉差不超过 $r\delta$。故两份真实自主过程满足

$$
p_v-p_s\ge
\cos^2(2\pi/5)-2\delta-\varepsilon_s(1)-\varepsilon_v(1)
>0.0934>0.09.
\tag{11}
$$

代码以 $\sqrt5<223607/10^5$、$\cos^2(2\pi/5)=(3-\sqrt5)/8$ 和精确有理预算认证。两份过程的质量、正常准备、源、空间核及读口完全相同，只改变实际物理源型。若都被称为同一标量材料的共同预测，各1%的误差要求被拒收；矢量材料本身没有因此被排除。

(10)与(11)是同一H/部署家族的两个可选终止读窗，不宣称先作破坏性钟测量后又免费恢复末端未测过程。

## 6. 同H能源、反作用与有限因果范围

[来源资格](source_car_review.md)独立核新准备满足

$$
\sup_{\|z\|=1}\|H_D\mathcal Jz\|<9\times10^{10}.
$$

光滑有限CAR模式及其生成元交换子保域；源准备在 $D(H_\phi)$。有界 $V_a$ 使真实联合输出始终在 $D(H_D+H_\phi)$，

$$
\langle H_D+H_\phi\rangle_t
<9\times10^{10}+1073A^2+2\|V_a\|,
\quad\Delta\langle H_D\rangle+\Delta\langle H_\phi\rangle+\Delta\langle V_a\rangle=0.
$$

不能从小态距替代这个无界能源证明。1129的CCR与有界 $F'$ 推导给同一真实源力

$$
\frac d{dt}\Pi_H(g)=\phi_H(\Delta g)
-\frac\lambda A\int\eta(x)(\rho_\epsilon*g)(x)
[F'(\phi_\epsilon(x)/A)J_{a,\epsilon}(x)]_Hdx.
\tag{12}
$$

实际测试源相对自由源的归一退迟差不超过 $\lambda512C_*q_\epsilon t^2/(2A^2)<1.1\times10^{-96}$。结合真空方差，$t=0,1$ 的实际交互均能绝对值小于 $8.1\times10^{-29}$。这些界保实际算符力与完整相关性，未将反作用设为0。

同一材料和源共同平移的生成元 $P=P_D+P_\phi$ 仍满足

$$
\frac d{dt}P_{j,H}=-\lambda\int(\partial_j\eta)(x)
[F(\phi_\epsilon(x)/A)J_{a,\epsilon}(x)]_Hdx.
\tag{13}
$$

固定支撑仍承担外力责任，未纳入动态支撑就不能称全装置四P闭合。

有限因果菜单复用1129 §5：异质量自由Dirac主部速度仍为1，有限核一次最多跳 $2\epsilon$；$D_a\sigma_x$范数1不改变Dyson计数。在原声明的材料偶干预及单模式效果菜单、$\operatorname{dist}(X,Y)>t+20\epsilon$ 下，两支实际概率差上界

$$
4\frac{\lambda^{10}}{10!}\frac1{1-\lambda/11}+2e_s(1)<0.00012.
$$

不授予任意来源干预或全部局部代数diamond菜单，不将双点核运输称为点态壳流。完整局部交叉Ward仪器证书仍开放。

## 7. 准确增量与下一步

本轮在一个明示有限材料中把实际移动钟、标量操作权重、区域读口、完整正常来源/场后态及反作用连接起来。没有只把既有 $m/E$ 公式重新并列，也没有仅重跑1127固定矩阵反差。认知动机仍是同一操作、记录与资源责任须共同成立；新场类型、参数和末读权限全部作为输入登记。

未完成的共同接口包括原CO位置/288维钟身份、1144独立两站记录、全部源世界管的共同boost、动态支撑与终末探测器的四P，以及原点态局部Ward/应力/壳流。1015—1019与1027—1029的场论/引力定理只有其原前提共同满足后才可领取。有限分辨率是当前合同，既不预设微观连续，也不推定物理最小尺度。

本轮准入从待验的点态联合读窗转为可完成的有限分辨率共同材料证书，原因是它直接补实际钟/源权的同一身份，并保来源后态；没有据此将剩余点态流或完整阶段目标删除。停止继续缩小核、增大质量/来源或优化误差常数；下一项须补一个实际未连接接口。新认知公理0、科学校准0，阶段目标active。
