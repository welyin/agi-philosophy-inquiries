# 1249：同一齐次父作用的内部钟规约与完整来源字典

2026-10-15。接[准入](admission.md)。本文只签正则与来源接口，有限量子误差由本轮另一份证明承担。CT仍为待验证猜想；先前材料、有效几何及条件Einstein分支不改。

## 1. 同一物理体积族，不能用坐标改名制造小量子参数

继承[1245](../1245/proof.md)、[1246](../1246/proof.md)、[1247](../1247/proof.md)的同一作用。固定$\kappa=3$、$\Lambda=0$、$k^2=1/10$、$n_c=1$、$n_r=1/10$及全部局部耦合。令空间周期为$10n$，$n$是有限正整数。每份compact杆的**靶周期保持$10k$不变**，空间绕数取$n$，于是

$$
\phi^A=kx^A,\qquad \partial_i\phi^A=k\delta_i^A,
\qquad V_n=1000n^3,\quad N_c=V_n,\quad N_r=V_n/10.
\tag{1}
$$

$N_c,N_r$是原连续有效标签测度的总量，不是声称有限点粒子精确形成均匀应力。此变化增加实际总体积、杆绕数、材料与准备资源；原杆归一及$a(0)=1$保持，故不是坐标重标。所有局部能源密度与Newton尺度不变。下文将$V_n$简记为$V$。

三杆的零模若先写成$\phi^A=kx^A+\chi^A(t)$，其动量贡献为$P_{\chi A}^2/(2Va^3F)$。原世界线机械动量为零、度规各向同性，空间约束给$kP_{\chi A}=0$；原$P_{\chi A}=0$及常$\chi^A$部门是一致不变部门。只在这个部门固定杆零模，**杆的绕数梯度能源仍保留**。非零杆动量或相对流动不在本次约化内。

## 2. 总动量、符号及未规约约束

在$ds^2=-N^2dt^2+a^2d\boldsymbol x^2$、shift为零的同质子族，选后引力Lagrangian为

$$
L_g=-\frac{Va\dot a^2}{N},\qquad
P_a=-\frac{2Va\dot a}{N},\qquad
C_g=-\frac{P_a^2}{4Va}.
\tag{2}
$$

负号不是物质的负质量；它是原引力约束的同一动能符号。不能为了量子稳定性在此改成正号。

场动量及世界线内部**总**动量为

$$
P_p=Va^3\dot p/N,\quad P_X=Va^3\dot X/N,
$$

$$
P_\theta=N_cJ,\quad P_q=N_cr,\quad P_\psi=N_cB,
\quad P_Y=N_cK,\quad P_\vartheta=N_r I.
\tag{3}
$$

这是把原标签辛项积分后的标准canonical配对；尤其$P_\vartheta$不能误写为$I$。在原共动且$M_c,M_r>0$的部门，世界线能量$\sqrt{M^2}=M$，同一作用成为

$$
S=\int dt\left[P_a\dot a+P_p\dot p+P_X\dot X
+P_q\dot q+P_\psi\dot\psi+P_Y\dot Y
+P_\vartheta\dot\vartheta+P_\theta\dot\theta-NC\right],
\qquad C=P_\theta+H_{\rm red}(\theta),
\tag{4}
$$

其中

$$
\begin{aligned}
H_{\rm red}={}&-\frac{P_a^2}{4Va}
+\frac{P_p^2+P_X^2}{2Va^3}+\frac{Va^3}{2}(p^2+X^2)
+\frac{3V}{20}aF(p)\\
&+N_c\left[\frac{71}{100}+\frac{q^2}{2}
+\frac{P_q^2+P_Y^2}{2N_c^2}
+\left(\frac{P_\psi}{N_c}-1\right)\frac{w(\theta)qp}{4}
+\frac{f(\theta)XY}{10}\right]\\
&+\frac{N_r}{10}+\frac{1+p}{10}P_\vartheta,
\qquad F=1+p/5.
\end{aligned}
\tag{5}
$$

$w,f$仍为1247的两个$C^2$紧时窗脉冲，未换成外部任意控制。$71/100=7/10+1/100$含钟基线与记忆静质量；最后两项含原记录器静质量及$I$能量。整个$H_{\rm red}$不必正，原量子算符的自伴性不能仅凭某个动能正性宣称。

## 3. 内部钟规约与七对强度变量

由$\{\theta,C\}=1$，有限单调钟片$\theta=t$给$N=1$。解$P_\theta=-H_{\rm red}$后，规约作用是其余七对变量的$\int(P\dot q-H_{\rm red})dt$。时变$H_{\rm red}(\theta)$来自原自治钟，不是在未规约作用外另加脉冲；不过这一步尚不证明量子钟实际达到指定时间分辨率。

取$Q=\log a$，$P_Q=aP_a$，把其余总动量也除以$V$。七对强度变量为

$$
(Q,\pi_Q),(p,\pi_p),(X,\pi_X),(q,\pi_q),
(\psi,\pi_\psi),(Y,\pi_Y),(\vartheta,\pi_\vartheta),
\tag{6}
$$

$$
\pi_Q=P_Q/V,\quad \pi_q=r,\quad\pi_\psi=B,
\quad\pi_Y=K,\quad\pi_\vartheta=I/10.
$$

令$h=H_{\rm red}/V$；为免混淆，原程序耦合始终写$1/4$，不用半经典参数符号表示它。则

$$
\begin{aligned}
h={}&e^{-3Q}\left[-\frac{\pi_Q^2}{4}
+\frac{\pi_p^2+\pi_X^2}{2}\right]
+\frac{e^{3Q}}2(p^2+X^2)+\frac3{20}e^QF\\
&+\frac{71}{100}+\frac{q^2+\pi_q^2+\pi_Y^2}{2}
+(\pi_\psi-1)\frac{w(\theta)qp}{4}
+\frac{f(\theta)XY}{10}
+\frac1{100}+\frac{1+p}{10}\pi_\vartheta.
\end{aligned}
\tag{7}
$$

规约作用为$V\int(\pi\dot q-h)dt$。因此强度经典方程与$V$无关；标准强度Hamilton方程使用$h$，对应原总作用的PB则带$1/V$。量子化时

$$
[\widehat q^i,\widehat\pi_j]=i\varepsilon\delta^i_j,
\qquad \varepsilon=\hbar/V,
\qquad i\varepsilon\partial_\theta\Psi=\operatorname{Op}_\varepsilon(h)\Psi.
\tag{8}
$$

式(8)是待给域与延拓的共同线性量子候选；几何、材料、控制、程序及记忆均在一个生成元内，不以状态依赖的期望值方程代替。

## 4. 完整旧方程与初始能源的还原

定义

$$
u=e^{-Q},\quad H=-\tfrac12e^{-3Q}\pi_Q,
\quad v=e^{-3Q}\pi_p,\quad V_X=e^{-3Q}\pi_X,
\quad J=-h,\quad z=\vartheta-\theta/10.
\tag{9}
$$

式(7)给$\dot Q=H$及原$p,X$方程、

$$
\dot q=\pi_q,\quad\dot\pi_q=-q-(\pi_\psi-1)wp/4,
\quad\dot Y=\pi_Y,\quad\dot\pi_Y=-fX/10,
$$

$$
\dot\psi=wqp/4,\quad\dot\pi_\psi=0,
\quad\dot\vartheta=(1+p)/10,\quad\dot\pi_\vartheta=0,
\quad\dot J=-\partial_\theta h
=-(B-1)w'qp/4-f'XY/10.
\tag{10}
$$

$\dot\psi$不乘$B-1$。$J$被重建而未从资源账消失。直接微分(9)还给

$$
\dot H=-\frac32(H^2+P);
\tag{11}
$$

配下节同一约束$H^2=\rho$，正是1247的完整几何反馈。所有旧变量均恢复；$L=\int e^{-Q}d\theta$保留为原null可达积分，未新增无源canonical记录器。

初值为$Q=p=\pi_p=\pi_X=\pi_q=Y=\pi_Y=\psi=\vartheta=0$，$X=1/10,q=1/2,\pi_\vartheta=1/10,\pi_\psi=1$或$2$，$\pi_Q=-2\sqrt{111/100}$。每支$h(0)=-1/10$，故$J(0)=1/10$。在原$V=1000$时，

$$
H_{\rm red}(0)=-1110+5+150+835+20=-100,
\qquad P_\theta+H_{\rm red}=100-100=0.
\tag{12}
$$

835是除去$J$后的钟/振子/记忆总能，20是原记录器；各项不可移掉。扩容仅把(12)全部项乘$n^3$。两支的局部密度、初$H$、全程强度历史及原严格$\Delta K$不变。生成元常数的量子全局相位自由不准许删除物理能源：若改$h$常数，$J=-h$及原约束/源字典也改变，已不再是同一已声明初值。

## 5. 独立材料来源、已解约束及漏源

在代入约束**之前**，来源由原材料变分给出，保留独立的$J$记账变量：

$$
\begin{aligned}
M_c={}&\frac{71}{100}+J+\frac{q^2+\pi_q^2+\pi_Y^2}{2}
+(\pi_\psi-1)wqp/4+fXY/10,\\
M_r={}&\frac1{10}+(1+p)\pi_\vartheta,\\
\rho_{\rm mat}={}&\frac12\left[e^{-6Q}(\pi_p^2+\pi_X^2)+p^2+X^2\right]
+\frac3{20}Fe^{-2Q}
+e^{-3Q}\left[M_c+\frac1{100}+\frac{1+p}{10}\pi_\vartheta\right],\\
P_{\rm mat}={}&\frac12\left[e^{-6Q}(\pi_p^2+\pi_X^2)-p^2-X^2\right]
-\frac1{20}Fe^{-2Q}.
\end{aligned}
\tag{13}
$$

这里$M_r$已经含$I=10\pi_\vartheta$，再乘数密度$1/10$进入$\rho$。所有控制、记忆和交互只计一次。式(10)给

$$
\dot M_c=(B-1)wqv/4+fYV_X/10,
\qquad \dot\rho=-3H(\rho+P).
\tag{14}
$$

代入$J=-h$后，(13)恒等于$\rho_{\rm mat}=H^2$。这是**已经求解的原Hamiltonian约束**，不能当作对同一方程的独立实验验证；(13)的材料分项与来源权限仍须保留，不能用$H^2$替换实际材料菜单。

1248的同历史漏源保留原$h,J$及全部轨道，只在账中删除

$$
m_K=\frac1{100}+\frac{\pi_Y^2}{2},\qquad
\rho_{\rm omitted}=\rho_{\rm mat}-m_Ke^{-3Q},\quad P_{\rm omitted}=P_{\rm mat}.
\tag{15}
$$

若先从生成元删$m_K$再重新解$J$，新的钟能会补偿删除；那是换了一份约束/历史，不是原漏源检验。

### 5.1 排序不能用已解约束掩盖

本轮共同量子桥须说明使用(13)的规约符号及同一Weyl字典，或给等价的未规约来源运输误差，不能另设$\widehat\rho=\widehat H^2$冒充独立定义。即使规约符号$\rho=H^2$恒等，通常仍有

$$
\operatorname{Op}_\varepsilon(H^2)\ne
[\operatorname{Op}_\varepsilon(H)]^2.
\tag{16}
$$

对未截断的局部符号$H=-e^{-3Q}\pi_Q/2$，Weyl乘积直接给

$$
[\operatorname{Op}_\varepsilon(H)]^2
=\operatorname{Op}_\varepsilon\!\left(H^2+\frac9{16}\varepsilon^2e^{-6Q}\right)
\tag{17}
$$

（在共同测试域上）。紧支延拓后的符号还含该延拓的导数项，不能把(17)当全局公式。单个算符平方的正性，也不自动赋予所有Weyl能源符号全局正谱。有限正域/能源尾预算由量子桥单独处理。式(15)的记忆删项只含彼此对易的$Q,\pi_Y$，其同历史定义不需重解$J$。

## 6. 可复算及适用范围

[canonical_check.py](canonical_check.py)、[canonical_results.json](canonical_results.json)合计一个科学组：有理Dual逐项核Hamilton向量场与来源抵消，核初片(12)、体积缩放及局部Weyl平方修正；仅读1247已保存的$t=1,2$端点作实现检查，不重积分。有限检查不替代本文的解析变分。默认只读；首次保存排他创建。

七对是**有限相空间维数**，不等同有限Hilbert空间。齐次约化在量子化之前进行，不证明全场非零模可忽略、全量子约束无异常或全CC/TR角色实现。局部action-angle图及有限单调$\theta$片不构成无限精确理想钟；$J,M_c,M_r,F$的经典正管道不会自动限制量子态的全部支撑。$K=\pi_Y$作为零模量子变量是集体强度，不自动等同每个局部标签的量子读数。

对$\psi,\vartheta$及关系钟的相位，本轮只选覆盖既定有限历史的展开坐标和明确正则延拓；这不是证明action-angle全局拓扑等于实直线，也不授予跨周期计圈。相应Weyl量子模型是这份有效延拓的模型，不能无证宣布与全部原compact相位量子化完全相同。

固定物理$\hbar$时，减小$\varepsilon=\hbar/V_n$靠有限但更多资源；集体均匀准备与共同读取是额外有效输入。未计算有限量子误差前，原固定$V_0=1000$不能被宣布满足旧记忆分辨率。本文不宣称生成完整量子引力、不新增认知公理，也不把量子桥的适用域问题变成所有引力路线必须先解决的门槛。
