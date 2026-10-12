# 1154证明：完整裸记录任务的累计广义支撑力

2026-10-12。按[准入](admission.md)，给同一1153父Hamiltonian族的一项新增观测合同。新实例的耦合与期限明确改变，不领取1153原固定参数的冲量精度。

## 1. 同一对象、参数与要证明的算符

复用[1153](../1153/proof.md)的全正常CAR、海、KG、两独立MIT站、双味/双spin、径向核、$\delta=1/16$及裸码等距$W_0$。在材料随站点的内部坐标表示中，取$y_A=0,y_B=d e_z$；$d$只出现在B场核$f_{B,x}(z)=\kappa_\delta(z-d e_z-x)/L$。变化$d$表示整份B材料、壁及其固定核一起改变位置，不改变站内材料谱和裸码坐标。

$d=3$为实际检验值；为取导数可在附近非重叠开区间定义同一族。各站自项不随$d$变：单站平移全部场变量后，自由KG及真空不变。只需实参数一阶光滑；没有假设紧支bump对复$d$解析。

定义实际广义力和码内积分：

$$
F_B(d)=-\partial_dH(\lambda,d)=-\lambda V_d,\qquad
V_d=\partial_d V_B,
$$

$$
\mathcal I_B(T)=W_0^\dagger\int_0^T
e^{itH}F_B e^{-itH}\,dt\,W_0,\qquad
T=\frac{\pi}{C(3)\lambda^2}.
\tag{1}
$$

取$\partial_d$时始终固定$\lambda$，不对任务期限求导；实际积分只在$d=3$、事先确定的$T$下进行，不能把每个$d$重新调时后的$TC(d)$恒定误读为力为零。

积分取共同能源域的弱形式。后面的统一有限界使其在16维码上成为有界自伴矩阵。对任意码态、与任意被动R的纠缠态，期望由同一矩阵决定。证明不要求实际测得每条微观轨迹，也没有给该力的测量仪器。

$F_B$是固定支撑位置的共轭力，不是已构造自主墙的动量算符。若另引入支撑坐标$X_B$和动量$p_B$、以同一$H(\lambda,X_B)$耦合，则其Heisenberg运动方程会含$-\partial_{X_B}H$；本轮不靠这一形式延伸领取自主墙或全四P。

## 2. 外区系数：直接复用1149，不重复发明静态力

1153实际码场剖面$q$非负、径向、光滑，支持半径$9/16$。由[1149 §3—§6](../1149/proof.md)已证明的精确球外Yukawa角积分，对$d>9/8$，

$$
C(d)=A_q\frac{e^{-d}}d,\qquad
A_q=\frac{Q_q^2}{4\pi}>0,\qquad
Q_q=\int q(x)\frac{\sinh|x|}{|x|}\,d^3x.
$$

因此

$$
\frac{C'(d)}{C(d)}=-1-\frac1d,\qquad
\frac{\pi C'(3)}{C(3)}=-\frac{4\pi}{3}.
\tag{2}
$$

二阶有效能源的交叉项是$-\lambda^2 C(d)P_{1,A}P_{1,B}$，所以$+\lambda^2 C'(d)P_{1,A}P_{1,B}$是领先B力，朝向A。小$\lambda$降低瞬时力，但任务时间$T$增加，(2)提示累计力不能直接忽略。

**本轮的新命题是(1)的真实全正常裸准备积分接近(2)，不是再次求静态势导数。** 仅有Hellmann–Feynman本征态公式或1153的态距上界均不足以证明它。

## 3. 力的图界和解析压缩

在每个净味荷部门$ab$，令$E_0=e_a+e_b$、$A=H_0-E_0\ge0$，沿用$0<g\le1$。定义固定有限常数

$$
\eta=\max\{1,\|\partial_z\kappa_\delta\|_2/L\},\quad
D=\frac{10}{\sqrt g}+\frac5g,\quad
R_0=\frac1{4D},\quad D_1=\eta D .
$$

固定平滑核使$\eta$有限。标准场一次图界给

$$
\|V_d\psi\|\le\eta\left(5\|A^{1/2}\psi\|+\frac52\|\psi\|\right),\qquad
\|V_d(A-z)^{-1}\|\le D_1\quad (|z|=g/2).
\tag{3}
$$

这里$H_\chi\le A$；仅B核取导，使用两站共同上界是保守放大。

令$P_\lambda$为真实rank-4最低簇，$P_0$为裸簇，$U_\lambda$为1153的解析直接旋转。对于复$|\lambda|\le R_0$，Riesz公式、$\|\lambda V(A-z)^{-1}\|\le1/4$给

$$
\|V_dP_\lambda\|
\le\frac g2\frac{D_1}{1-|\lambda|D}
\le\frac{2gD_1}{3}.
$$

实际力的固定空间压缩

$$
B_\lambda=P_0U_\lambda^{-1}F_B P_\lambda U_\lambda P_0
\tag{4}
$$

是解析矩阵；实$\lambda$时正是实际低簇的酉表示。由$\|U_\lambda^{\pm1}\|\le\sqrt2$，在圆盘上

$$
\|B_\lambda\|\le\frac{4g}{3}|\lambda|D_1
\le M_F:=2gR_0\eta D .
$$

KG奇偶将$F_B(\lambda)$变成$F_B(-\lambda)$，在裸簇恒等，因此$B_\lambda$为偶解析，且$B_0=0$。直接展开到二阶，记$R=Q_0A^{-1}Q_0$，

$$
B^{(2)}=P_0(V_dRV+VRV_d)P_0
=-\partial_dK^{(2)}=C'(d)ab\,I_4 .
\tag{5}
$$

$P_0,A$在这份内部坐标表示中不随$d$变化。自项的平移不变和1153跨项身份给最后等号。也可由有效矩阵导数得到(5)，但不能在高阶随意删掉直接旋转连接；(4)直接压缩真实力避开这个问题。

Cauchy系数界遂给

$$
\|B_\lambda-\lambda^2C'(d)abI_4\|
\le M_F\frac{(\lambda/R_0)^4}{1-(\lambda/R_0)^2}
=:r_F .
\tag{6}
$$

保全部高阶spin作用；没有把退化矩阵化成四条独立标量能源。谱工具可对照[Bravyi、DiVincenzo、Loss](https://arxiv.org/abs/1105.0675)，当前无界力的资格由(3)—(6)另证。

## 4. 原裸准备的全部尾部

以下力通式适用于实$0<\lambda\le R_0/4$；第7节的共同记录新合同会另取更小参数。设$P=P_\lambda,Q=1-P$，以

$$
p=\frac{\lambda D}{1-\lambda D}\le\frac43\lambda D
$$

作统一上界。除了$\|P-P_0\|\le p$，Riesz差公式还给

$$
\|A^{1/2}(P-P_0)\|\le\sqrt g\,p .
\tag{7}
$$

具体在同一圆周上，$\|A^{1/2}(A-z)^{-1}\|\le2/\sqrt g$，Riesz积分半径$g/2$，Neumann差因子至多$p$；而$A^{1/2}P_0=0$。

对原裸单位码向量$w$，尾$\chi_0=Qw$满足

$$
\|\chi_0\|\le p,\qquad \|A^{1/2}\chi_0\|\le\sqrt g\,p .
$$

这不是将初态换成衣着态。所有后续$\chi_t=e^{-itH}\chi_0$都保留。由1153能源下界（在本部门减去$E_0$）

$$
H-E_0\ge\tfrac12 A-K_\lambda,\qquad
K_\lambda=\tfrac{25}{2}\lambda^2+\tfrac52\lambda,
$$

及真实$H$期望守恒，

$$
\|A^{1/2}\chi_t\|^2
\le2p^2\left[g+\lambda(5\sqrt g+5)+\tfrac{25}{2}\lambda^2\right]
\le4gp^2 .
\tag{8}
$$

最后用$R_0\le g/20$、$g\le1$；括号乘2不超过$(49/16)g<4g$。这项图资源证明与小态距独立。

## 5. P/Q交叉项不随T线性积累

在真实P簇上，偶解析K给

$$
\|(H-E_0)P\|\le
2g\frac{(\lambda/R_0)^2}{1-(\lambda/R_0)^2}<g/4 .
$$

Riesz圆排除Q谱在$(-g/2,g/2)$，而$H-E_0\ge-K_\lambda>-g/2$排除圆左侧谱，所以$H_Q-E_0\ge g/2$。P/Q谱间隔至少$g/4$。

由(3)、(7)，$\|F_BP\|\le5\eta\lambda$。Sylvester方程

$$
H_QX-XH_P=QF_BP
$$

有有界解；将两谱间选一分割点后用正时间指数积分，可得

$$
\|X\|\le\frac4g\|QF_BP\|\le\frac{20\eta\lambda}{g}.
$$

时间导数把交叉积分化为两个端点，

$$
\int_0^T e^{itH_Q}QF_BP\,e^{-itH_P}\,dt
=\frac{e^{iTH_Q}Xe^{-iTH_P}-X}{i}.
$$

压回裸输入并加共轭交叉项，完整算符范数至多

$$
4p\|X\|
\le\frac{320}{3}\frac{\eta D}{g}\lambda^2
<128\frac{\eta D}{g}\lambda^2 .
\tag{9}
$$

不能用$T$乘瞬时交叉误差代替这一估计；后者会丢失隙带来的振荡控制。

## 6. P/P主体与Q/Q尾的长期预算

对P部分，在真实簇表示中用(6)。领先力在每个$ab$部门为标量，与实际高阶K精确对易。原裸输入的P权重满足

$$
W_0^\dagger P W_0=I-W_0^\dagger QW_0,\qquad
\|W_0^\dagger QW_0\|\le p^2 .
$$

所以领先权重损失至多$(4\pi/3)p^2<8D^2\lambda^2$。剩余P/P积分至多

$$
Tr_F\le \frac{8M_F}{R_0^4}\frac{\lambda^2}{C},\qquad C=C(3).
\tag{10}
$$

此处用$\pi<4$及$1/(1-(\lambda/R_0)^2)<2$。

对Q/Q部分，(3)、(8)给任意单位裸输入

$$
|\langle\chi_t,F_B\chi_t\rangle|
\le\lambda\eta(10\sqrt g+\tfrac52)p^2 .
$$

Q/Q压回矩阵自伴，对所有单位向量的二次型界就是算符范数界。沿全任务积分，

$$
\left\|W_0^\dagger Q\int_0^T
e^{itH}F_Be^{-itH}\,dt\,QW_0\right\|
\le8\eta(10\sqrt g+\tfrac52)D^2\frac{\lambda}{C}.
\tag{11}
$$

**保留$\lambda/C$项。** 1153仅保证$\lambda^2/C$小，不能用奇偶或小态距把长期连续高谱尾凭空升级成四阶。无需散射完备性，也不要求尾回收。

四部门正交直和、力保净味荷，不放大范数。由(9)—(11)得全16码定理：

$$
\left\|\mathcal I_B(T)+\frac{4\pi}{3}P_{1,A}P_{1,B}\right\|
\le K_2\lambda^2+K_4\frac{\lambda^2}{C}
+K_t\frac{\lambda}{C},
\tag{12}
$$

$$
K_2=8D^2+\frac{128\eta D}{g},\qquad
K_4=\frac{8M_F}{R_0^4},\qquad
K_t=8\eta(10\sqrt g+\tfrac52)D^2 .
$$

所有常数由固定自由材料、核和谱隙确定，与$\lambda,T$无关。张量任意R保持(12)，覆盖未知味、自旋及其相关性。这里证明的是积分**期望**矩阵；没有控制未声明的力分布、方差或自主壁末态。

## 7. 一次性有限新合同及继承的强记录

设$\lambda_{1153}$是旧配方数值，令$\varepsilon=10^{-3}$，

$$
\lambda_*=\min\left\{
\lambda_{1153},
\sqrt{\frac{\varepsilon}{3K_2}},
\sqrt{\frac{\varepsilon C}{3K_4}},
\frac{\varepsilon C}{3K_t}\right\}>0,\qquad
T_*=\frac{\pi}{C\lambda_*^2}<\infty .
\tag{13}
$$

(12)三项各不超过$\varepsilon/3$。这项新实例中，裸11准备的累计B广义力在$-4\pi/3$的$\varepsilon$邻域，绝对值大于4.18779；对原记录准备A=1、B=+，期望在$-2\pi/3$的$\varepsilon$邻域。这里沿**无中插的记录实验**计算；中插历史是同父模型另一次允许实验，不继承这个力值。

由于$\lambda_*\le\lambda_{1153}$，1153所有均匀配方条件仍成立。在这个同一新实例，完整记录误差仍小于0.000107，实际钟误差小于0.000114，一次非对易继续小于0.000207；读口相位和长时限按新$\lambda_*,T_*$计算。整站短窗因果及独立能源域也保留。没有把旧固定实例的性能偷换为新实例；两者都属于同一明确族，新合同只额外要求累计力精度。

各站空间参数都开放时可同法得到相反的A领先广义力；对全字段总平移有相应真实动量交换弱身份。但本轮签收的是(12)的B位置力，不因此声称全装置四P已闭合。

## 8. 方向含义和边界

认知动机是：同一个支持系统既允许记录形成，也承担维持其关系位置的交换责任。数学结果显示，在本候选和本任务内，“瞬时交互弱，所以支撑反作用可免费省略”不是有效近似；等待时间与力强度不能分别忽略。

保持站点不动，允许其它装置补偿相反力；大累计冲量不推出必须有正净功、耗散下界、最低质量或不可实现性。本轮没有新增普遍认知公理，也未证明该结构由现行公理唯一推出。反证范围仅是本合同下把这份广义力替换为零的匹配，其对裸11积分期望至少有$4\pi/3-\varepsilon$误差。

下一步优先判断共同运动描述、真实支撑及旧局部应力前提能否在有限有效范围内连接，不为减小上述保守常数继续扫描参数。完整QFT、点态UV、MIT墙工程和全材料微观设计仍不是阶段前置门槛。
