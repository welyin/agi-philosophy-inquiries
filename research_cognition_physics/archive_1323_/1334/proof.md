# 1334：固定矩阵重链的共同质量、弱流与门户约束

2026-10-19。[事前准入](../_admission/round1334_admission.md)。新公理0；本轮是固定候选的条件性联合约束，不是自主认知形成或现实标准模型完成。

## 1. 对象与输入

动机是层级复用同一步骤，两部门共享接续规则。具体语法仍是输入：$n=3$、同向开链$S=|1\rangle\langle0|+|2\rangle\langle1|$、共同整数$m\ge1$和共同连续比值$r>0$，

$$
D_f=r(I+q_fS),\qquad
K_f=M_f\begin{pmatrix}I&-D_f&&\\&I&-D_f&\\&&\ddots&\ddots\end{pmatrix},
\qquad \mathcal M_f=\begin{pmatrix}h g_fP_m\\K_f\end{pmatrix}. \tag{1}
$$

$P_m$取末端$R_m$，门户严格为$g_fI$，$M_f>0$、$g_f\ne0$、$h>0$。无任意入口矩阵、逐站权、换端点或自由动能。重左右站取各部门物理右夸克的规范表示，故其左站是弱单态；两部门共享原$Q_L$弱流。原规范群、三代、场表示和保护角色仍输入。此处$r$不是旧离散保护生成元，后者若需指称记$r_{\rm disc}$。

以下“完整质量”均指(1)在树级常Higgs背景下的全部奇异值，不是非微扰重整化后的精确极点。来源由同一个母核取导数；并未构造自治控制、动态认知参考或全量子仪器。

复用[1316](../../archive_1313_1322/1316/proof.md)的规范化零模、重响应、Schur极点和重隙方法；[1317](../../archive_1313_1322/1317/proof.md)的入口及完整来源限制继续保留。净增是取消任意$Y_0$后，两部门共享少标量语法的实际约束，非首次发现clockwork。

## 2. 六矩决定正则零模

因$S^3=0$，$D^k=r^k[I+kqS+\binom{k}{2}q^2S^2]$。定义

$$
(A,B,C,D,E,F)=\sum_{k=0}^m r^{2k}
\left(1,k,\binom{k}{2},k^2,k\binom{k}{2},\binom{k}{2}^2\right).
$$

本节$D$为标量矩，$D_f$始终表示一步矩阵。令$x=|q|$，直接相乘得

$$
Z=\sum_{k=0}^m(D_f^k)^\dagger D_f^k
=\begin{pmatrix}
A+Dx^2+Fx^4&q^*(B+Ex^2)&C(q^*)^2\\
q(B+Ex^2)&A+Dx^2&Bq^*\\
Cq^2&Bq&A
\end{pmatrix}\succ0. \tag{2}
$$

$Z\succeq I$。零模和领先耦合为

$$
V_0=\begin{pmatrix}D_f^m\\D_f^{m-1}\\\vdots\\I\end{pmatrix}Z^{-1/2},
\quad K_fV_0=0,\quad V_0^\dagger V_0=I,
\qquad Y_f=g_f Z_f^{-1/2}. \tag{3}
$$

$Z$由全部原正则动能固定，不能自由选来修谱；$D_f^k$也不是量子时间演化。对$q=xe^{i\phi}$，

$$
Q_\phi=\operatorname{diag}(1,e^{i\phi},e^{2i\phi}),\qquad
Z(q)=Q_\phi Z(x)Q_\phi^\dagger. \tag{4}
$$

因此单部门领先谱只依赖$m,r,x,|g|$，两部门相对相位须通过同一左字典保留。

## 3. 共同$r$是真约束：$m=2$的局部六维结果

定义形状不变量$I_2=\operatorname{tr}Z^2/(\operatorname{tr}Z)^2$、$I_3=\det Z/(\operatorname{tr}Z)^3$。它们可从$Y Y^\dagger$的逆矩阵归一恢复，不需先知道$g$。取$\tau=r^2,p=|q|^2$，在$m=2,\tau=p=1$，

$$
Z_* =\begin{pmatrix}9&5&1\\5&8&3\\1&3&3\end{pmatrix},\quad
\partial_\tau Z_* =\begin{pmatrix}14&9&2\\9&12&5\\2&5&3\end{pmatrix},\quad
\partial_p Z_* =\begin{pmatrix}7&9/2&1\\9/2&5&3/2\\1&3/2&0\end{pmatrix}.
$$

其迹20、平方迹224、行列式82给出

$$
\frac{\partial(I_2,I_3)}{\partial(\tau,p)}=
\begin{pmatrix}151/1000&123/1000\\-1017/80000&-363/40000\end{pmatrix},
\quad\det=\frac{3093}{16000000}>0. \tag{5}
$$

故每部门两质量形状在此邻域局部决定$r,|q|$。两部门共用$r$时，四个形状量的像局部为三维；不是仅减少参数计数。可取$(I_{2,u},I_{3,u},I_{2,d})$，其对$(\tau,p_u,p_d)$的Jacobian行列式为$(123/1000)(3093/16000000)\ne0$。

再加两部门整体质量尺度恢复$|g_u|,|g_d|$。剩余相位$\delta=\phi_d-\phi_u$在同一点有可观测不变量

$$
J_{\rm mix}=\frac{\operatorname{tr}(Z_uZ_d)}{\operatorname{tr}Z_u\operatorname{tr}Z_d}
=\frac{152+68\cos\delta+4\cos^2\delta}{400},\qquad
\left.\partial_\delta J_{\rm mix}\right|_{\pi/2}=-\frac{17}{100}. \tag{6}
$$

六个输入$|g_u|,|g_d|,r,x_u,x_d,\delta$既给上界，又有独立不变量给下界。因此在这一明确$m=2$邻域，普通**领先**Dirac质量／混合恰六个连续自由，较不受限三代两矩阵的十个量减少四个。结论不扩张到全部$m$的满秩分类、有限重混合或强CP。$M_u,M_d$不影响领先$Y$，却影响完整谱和来源，仍是独立输入；强CP及其它固定源相位不能借换基删除。

## 4. 非空CP见证及有限窗口

取$m=2,r=1,q_u=1,q_d=i$。$Z_u=Z_*$，

$$
Z_d=\begin{pmatrix}9&-5i&-1\\5i&8&-3i\\-1&3i&3\end{pmatrix},\quad
\det Z_u=\det Z_d=82,\quad
\operatorname{tr}[Z_u,Z_d]^3=25104i. \tag{7}
$$

对$H_f=Y_fY_f^\dagger=|g_f|^2Z_f^{-1}$，三阶无迹交换子及
$[A^{-1},B^{-1}]=A^{-1}B^{-1}[A,B]B^{-1}A^{-1}$给

$$
\operatorname{tr}[H_u,H_d]^3
=|g_ug_d|^6\frac{25104i}{82^4}\ne0. \tag{8}
$$

这是普通弱CP的充分证书，非强CP结论。两部门各自谱非简并，重隙为正；因此该非空领先点在足够小有限$h/M_f$下有连续的非简并轻谱及非零弱流相位窗口。

本轮事前选温和完整父例：$m=2,r=1,q_u=1,q_d=i,g_u=0.2,g_d=0.3,M_u=10,M_d=13,h=1$，所有量用同一有限模型单位。代码直接验完整九维质量核、原弱流及非零轻流四元积虚部，未将领先$Y$代替真实轻态。此点不是现实质量拟合，也不是自主形成过程。

## 5. 完整弱流：一条共同相位曲线及轻重补偿

$q$的相位在所有父站按$Q_\phi$运输。令$\mathcal Q_\phi=I_{m+1}\otimes Q_\phi$；在$g$已取正的字典中，

$$
\mathcal M(q)=\mathcal Q_\phi\mathcal M(x)\mathcal Q_\phi^\dagger. \tag{9}
$$

$g$的整体相位可另由右站／重左站换基运输；不据此删除其对强CP等外部源的影响。实$q$的完整左奇异向量取实正交矩阵$U_f^R$，按质量由小到大排序。令$P_Q$取原$Q_L$三行，$A_f=P_QU_{f,\rm light}^R$。它通常不是正交矩阵。实际轻带电弱流及完整流为

$$
V_{\ell\ell}=A_u^TQ_\delta A_d,
\qquad V_{\rm all}=(U_u^R)^TP_Q^TQ_\delta P_QU_d^R. \tag{10}
$$

每一轻／重矩阵元都是$\sum_{k=0}^2t_ke^{ik\delta}$，其中$t_k$由同一完整父矩阵的原弱投影给出、为实数。故$c=\cos\delta\in[-1,1]$对全部读口共同满足

$$
|V_{ij}|^2=(t_0-t_2)^2+t_1^2
 +2t_1(t_0+t_2)c+4t_0t_2c^2. \tag{11}
$$

取二次式在端点及区间内驻点的最小值，即得单项全相位下界。全部矩阵元共用**同一个**$c$。对非简并质量态，单项模方与质量本征态的重定相无关。若质量简并，式(11)在所选的实基及其相位运输字典中成立。仅按质量分辨时，可比较由简并谱投影定义的转移强度，即对应初末简并块内模方之和（等权平均初态再除以初态块维数）。若仪器能分别读出简并块内的状态，还须给出这套分辨字典。

令$G_f=(U_f^R)^TP_Q^TP_QU_f^R$，则

$$
V_{\rm all}V_{\rm all}^\dagger=G_u,\quad
V_{\rm all}^\dagger V_{\rm all}=G_d,\quad G_f^2=G_f.
$$

完整轻中性$T_3$投影$N_f=A_f^TA_f$满足

$$
N_u=V_{\ell\ell}V_{\ell\ell}^\dagger+V_{\ell H}V_{\ell H}^\dagger,
\quad N_d=V_{\ell\ell}^\dagger V_{\ell\ell}+V_{H\ell}^\dagger V_{H\ell}. \tag{12}
$$

轻中性流、带电流损失和轻重补偿不能分别校准。这里是带电弱流与$T_3$部分；颜色／电磁的同表示恒等发生器不受相同投影，不能笼统称全部规范流都投影。

## 6. 对全部$M_f$成立的门户必要界

$K_f$有三维零模，$\mathcal M_f^\dagger\mathcal M_f=K_f^\dagger K_f+h^2|g_f|^2P_m^\dagger P_m$。将Rayleigh商限制在$V_0$上，压缩为$h^2|g_f|^2Z_f^{-1}$。升序min–max于是给

$$
m_{f,3}^2\le \frac{h^2|g_f|^2}{\lambda_{\min}(Z_f)},
\qquad
\boxed{|g_f|\ge\frac{m_{f,3}\sqrt{\lambda_{\min}(Z_f)}}h}. \tag{13}
$$

这里$m_{f,3}$是完整树级母核第三小质量；将其认作最高已知轻质量是公开身份条件。界对所有$M_f>0$精确成立，增加重尺度不能降低它。不假设领先匹配已经准确。

由$Z\succeq(D_f^m)^\dagger D_f^m$及
$(I+qS)^{-m}=I-mqS+\binom{m+1}{2}q^2S^2$，得到无特征值计算的保守证书

$$
\sqrt{\lambda_{\min}Z}\ge
\frac{r^m}{\|(I+qS)^{-m}\|_F}
=\frac{r^m}{\sqrt{3+2m^2x^2+\binom{m+1}{2}^2x^4}}. \tag{14}
$$

此处Frobenius范数给出谱范数的上界。三角界$1+mx+\binom{m+1}{2}x^2$同样有效。下列两个校准点的Frobenius界较紧，两种估计的优劣随参数而变。例如$x=0$时三角界给出确切的谱范数1。

公开探索后固定样例为$m=5,r=8,x_u=10,x_d=8/5$，以$h=174$、历史同尺度最高轻质量171.7和2.89 GeV作标尺，质量来源仅承接[1333历史示例](../1333/proof.md)。不是盲测或最新拟合。完整堆叠零模矩阵的SVD给$\sqrt{\lambda_{\min}Z_u}\simeq227.2823$、$\sqrt{\lambda_{\min}Z_d}\simeq898.0456$，从而必要$|g_u|\gtrsim224.278$、$|g_d|\gtrsim14.9158$。

更保守的解析(14)也已使两个门户超过$4\pi$。$g^2/(16\pi^2)$为普通环展开的诊断；因此这点不能当作受控的普通微扰SM解释。它不是所有非微扰理论不存在的证明，也不是整类$(m,r,q)$排除。领先混合若另给全相位下界，仅用于说明同点的联合约束，不冒充失去微扰控制后的现实精确预测。

## 7. 同父重响应、有限误差与来源

令$G=KK^\dagger\succ0$、$\gamma^2=\lambda_{\min}G$。复用1316的分解，原$Q_L$字典中的重响应为

$$
C(z)=|g|^2P_mK^\dagger G^{-1}(G-zI)^{-1}KP_m^\dagger,
\qquad C_6=C(0)\succeq0. \tag{15}
$$

实际低极点的原弱分量$v$满足精确Schur方程

$$
\bigl[h^2|g|^2Z^{-1}-z(I+h^2C(z))\bigr]v=0. \tag{16}
$$

若$a=h^2|g|^2/\lambda_{\min}Z<\gamma^2$，三个最低平方质量不超过$a$，其余不小于$\gamma^2$。在$0\le z\le a$有$\|C(z)\|\le |g|^2/(\gamma^2-a)$。记$\kappa=h^2|g|^2/(\gamma^2-a)$，领先平方质量$a_i$与完整轻极点满足$a_i/(1+\kappa)\le z_i\le a_i$；按升序比较，不能换成逐矩阵元比例。

$K=M K_0$还给$C_6=|g|^2M^{-2}P_mK_0^\dagger(K_0K_0^\dagger)^{-2}K_0P_m^\dagger$。因此$M$保持领先$Y$却改变重响应；这个独立尺度必须保留。固定全部参数后，$C_6$、轻中性投影和弱流来自同一个核，不可另输入三个矩阵。

来源在原字典先取$\partial\mathcal M/\partial M,\partial r,\partial\Re q,\partial\Im q,\partial\Re g,\partial\Im g,\partial h$，连同共轭组成Hermitian母作用的源。共同$r$源同时作用两部门；不同$M_f,g_f$源不混为一份。代码核完整源及$C_6$链式导数，未将零模基变化忽略后称独立耦合。实际准备、控制和记录资源仍未构造；静态母核的来源恒等式不赋予免费仪器。

## 8. 层级极限的顺序

固定$m\ge2,r>0$，令$\Delta_2=DF-E^2$、

$$
\Delta_3=\det\begin{pmatrix}F&E&C\\E&D&B\\C&B&A\end{pmatrix}.
$$

它们是序列$\binom{k}{2},k,1$在正权$r^{2k}$下的Gram子式；$k=0,1,2$已保证严格正。令$x\to\infty$，分层Schur或特征多项式给$Z$降序谱

$$
\lambda_1\sim Fx^4,\quad\lambda_2\sim(\Delta_2/F)x^2,\quad
\lambda_3\to\Delta_3/\Delta_2.
$$

领先质量从轻到重记$m_1,m_2,m_3$，有

$$
\frac{m_1m_3}{m_2^2}\longrightarrow
\mathcal R_m(r)=\frac{\Delta_2^{3/2}}{F^{3/2}\sqrt{\Delta_3}}. \tag{17}
$$

共享$m,r$的两个深层级部门趋向同一个值，左方向同时趋向共同角色基，因而层级和近对齐可以共同存在。结论只在固定$r$极限成立，未给全$r$统一误差。$m=2$时$\mathcal R_2(r)=1$；$m=3,r=1$时为$\sqrt{6859/20000}$。

Cauchy–Binet展开显示$r\to\infty$时$F$由$k=m$主导、$\Delta_2$由$(m-1,m)$主导、$\Delta_3$由$(m-2,m-1,m)$主导，故$\mathcal R_m(r)\to1$。反过来先$r\to\infty$使最后一个$D_f^m$主导，再$x\to\infty$，倒数裸幂谱给$((m+1)/(m-1))^3$。两个顺序不同；不能套用1333后的裸幂排除。$m=1$没有这个三阶正Gram，保在一般公式内但不套用(17)。

## 9. 验证及结论范围

[check.py](check.py)直接构造全部有限母矩阵，核六矩、局部Jacobian和CP、有限弱流相位曲线／补偿／中性投影、全部$M$界的有限实例、同核重响应和源导数。温和点非空；探索点先因门户必要界拒绝普通微扰解释。所有数值是实现诊断，全类量词由以上证明承担；没有无截断场论或自主形成证书。

**新增选择力：** 两部门质量形状的真实局部约束、有限父弱流的共同相位曲线、轻重／中性同源身份及对所有$M$的门户下界。**剩余：** $n=3,m$、场角色、$r,q_f,g_f,M_f,h$及真实形成过程。不得将温和非空点称现实拟合，也不得将探索点失败称全类反证。

停止修已失败的强门户样例；后继先判断哪项独立形成规则或共同任务能减少剩余资料。当前未准入1335，不开启无目标的$m,r,q$扫描或器件扩建。
