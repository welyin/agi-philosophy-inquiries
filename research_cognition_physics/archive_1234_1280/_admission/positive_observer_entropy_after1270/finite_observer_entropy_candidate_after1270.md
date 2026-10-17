# 正能观察者的有限熵余项：待审准入计算

2026-10-16。只检验一种明确的有限材料代数；未预约1271，未采用角色，不修改1270或历史论文。它是CT选后热学描述的候选接口计算，不是引力生成、实际温度计或新的量子引力模型。

## 1. 问题与对象

Jensen–Sorce–Speranza（JSS）的[§5.3](https://arxiv.org/html/2306.01837v3#S5.SS3)在慢变观察者波包条件下将代数熵写成材料相对熵、观察者能源与观察者分布熵之和。其几何解释另外使用局部引力约束。这里不继承其III型局域场或约束资格，而问：**有限材料、正常正能观察者和任意非对易工作态中，这个热学部分能否有不借小迹距的显式熵误差？**

固定有限材料Hamiltonian $H=\sum_j E_jP_j$、$E_j\ge0$，$\beta>0$，$\sigma=e^{-\beta H}/Z$。$j$标记不同能量；简并空间完整保留。任意材料态$\rho$允许与被动R纯化，不删关联。观察者纯波包初始与MR独立；不把任意预存材料—观察者关联态纳入本式。观察者能源$q\ge0$有密度

$$
p(q)=\frac{q^{k-1}e^{-q/\vartheta}}{\Gamma(k)\vartheta^k},\qquad k>2,\quad\vartheta>0,
$$

取$f=\sqrt p$并在负半轴置零。它正常归一，$\langle q\rangle=k\vartheta$有限；$\|f'\|_2^2=[4(k-2)\vartheta^2]^{-1}$。先在整条实轴表示平移，再限制正能corner；不声称半轴上存在自伴正则时间算符。

考虑内作用dressing的直积分表示。纤维坐标$y$、材料矩阵和物理观察者能源$q=y-H$分别保留。令

$$
P(y)=\boldsymbol1_{[0,\infty)}(y-H),\qquad
D(y)=f(y-H)\rho f(y-H)^*,
$$

$$
\tau(A)=\frac{\beta}{Z}\int_{\mathbb R}e^{-\beta y}\operatorname{Tr}A(y)\,dy,
\qquad d(y)=\frac{Z}{\beta}e^{\beta y}D(y).
$$

直接有$PDP=D$、$\tau(P)=1$、$\tau(d)=1$。这是带中心的有限材料直积分corner；不是III型局域场交叉积的II因子。$\tau$是事前选择的共同参考迹，中心存在使其选择不能由抽象代数唯一性偷得。物理能源$H,q$先指定，不由态对数倒定义。

这份迹有明确的参考态实现：取材料参考$\sigma$与另一个正常观察者波包$f_\beta(q)=\sqrt\beta e^{-\beta q/2}\boldsymbol1_{q\ge0}$，其纤维密度恰为$D_\beta(y)=(\beta/Z)e^{-\beta y}P(y)$。因此该共同参考在区域代数上的态就是$\tau$，密度为$P$；并非从当前待解释工作态倒拟合迹。选择此参考仍是本模型输入，未证明现行认知公理唯一选中它，也未给实际制备设备。归一有限迹下$S_\tau\le0$且参考值为0，负值不与普通有限矩阵熵非负相矛盾。

## 2. 精确熵余项

设$\pi_j=\operatorname{Tr}P_j\rho$，$J$为此能量标签。组合变量$Y=q+E_J$有密度$p_Y(y)=\operatorname{Tr}D(y)$；J是证明所用的cq标签，不执行一次免费能量测量，Y也没有取得实际读口资格。对任一固定纯化$|\Phi\rangle_{MR}$，令

$$
\tau_R^j=\operatorname{Tr}_M[(P_j\otimes I)|\Phi\rangle\langle\Phi|(P_j\otimes I)]/\pi_j.
$$

零权重项省略。R取完整纯化；若原来已有伙伴，补充纯化后共同保留，不能只取原伙伴边缘冒充完整R。条件R态为$\sum_j\pi_jp(y-E_j)\tau_R^j/p_Y(y)$；相应MR条件向量纯，因此材料与R的条件熵相等。由经典—量子条件熵身份，

$$
-\int\operatorname{Tr}D(y)\log D(y)dy
=h(Y)+S(\rho)-I(Y:R).
$$

又有$h(Y)-h(q)=I(Y:J)$及$\int y p_Y(y)dy=\langle H\rangle_\rho+\langle q\rangle$。这里$h(q)$是能源分布的微分熵，不是纯观察者波包本身的von Neumann熵；固定能源单位后才有一致数值。故定义$S_\tau(d)=-\tau(d\log d)$时，**精确**得到

$$
S_\tau(d)=-D(\rho\Vert\sigma)-\beta\langle q\rangle+h(q)+\log\beta+\epsilon,
\qquad \epsilon=I(Y:J)-I(Y:R).
$$

条件链$Y\leftarrow J\rightarrow R$给

$$
\epsilon=I(Y:J\mid R),\qquad 0\le\epsilon\le I(Y:J).
$$

若$[\rho,H]=0$，R可完美区分各能量标签，$\epsilon=0$。若$\rho$纯，R为平凡系统，$\epsilon=I(Y:J)$。一般相干混态介于二者；没有先把它退相干。

## 3. 正能有限预算

令$E_*=\min_jE_j$，$a_j=E_j-E_*\ge0$。把输出参考密度固定为$p_0(y)=p(y-E_*)$。相对熵链式分解给

$$
I(Y:J)\le\sum_j\pi_jD\big(p(\cdot-E_j)\Vert p_0\big).
$$

对$a\ge0$，Gamma密度移位的相对熵为

$$
D_a=\frac a\vartheta-(k-1)\mathbb E\log(1+a/q).
$$

使用$\log(1+x)\ge x-x^2/2$、$\mathbb E(q^{-1})=[(k-1)\vartheta]^{-1}$与$\mathbb E(q^{-2})=[(k-1)(k-2)\vartheta^2]^{-1}$，得到

$$
0\le\epsilon\le
\frac{\operatorname{Tr}\rho(H-E_*)^2}{2(k-2)\vartheta^2}.
$$

正能、有限期望能源与有限熵误差因此有非空共同区间。若要$\epsilon\le\delta$，此构造的充分资源条件为

$$
\langle q\rangle=k\vartheta\ge
k\sqrt{\frac{\operatorname{Tr}\rho(H-E_*)^2}{2(k-2)\delta}}.
$$

这是所选波包的充分预算，不是所有观察者的最优下界。有限矩阵谱有界不需要任意细UV延拓；但此波包不是有硬上界能源的准备。

### 3.1 硬能窗版本

保持原正能corner及原迹，**只改变准备**为$Q_c=Q\mid 0<Q\le B$，$P_B=\Pr(Q\le B)>0$，密度$p_c=p/P_B$（窗外为零）。精确熵恒等式仍成立，把$\langle q\rangle,h(q)$换为截断分布的实际值。

仍以未截断$p_0(y)=p(y-E_*)$作共同比较。对$a\ge0$，

$$
D(p_c(\cdot-a)\Vert p)
=-\log P_B+\frac a\vartheta-(k-1)\mathbb E_c\log(1+a/Q).
$$

截去上尾后$\mathbb E_cQ^{-1}\ge\mathbb EQ^{-1}$，故线性余项非正；又有$\mathbb E_cQ^{-2}\le\mathbb EQ^{-2}/P_B$，从而

$$
0\le\epsilon_c\le-\log P_B+
\frac{\operatorname{Tr}\rho(H-E_*)^2}{2(k-2)\vartheta^2P_B}.
$$

此时能源逐态满足$0\le q\le B$。用较大的有限$B/\vartheta$使尾项小，再选有限$\vartheta$即可取得有限目标熵误差。没有把代数改成双边能窗corner后重新归一；若真的改成$0\le y-H\le B$的代数，其迹为$1-e^{-\beta B}$，必须另记常数。硬切波包仅要求$L^2$正常性，一般不在$H^1$；本结论不授有限共轭动量方差、平滑准备权限或实际制备成本。若这些也属任务要求，需另给平滑截止预算。

## 4. 尚不能领取的物理身份

本计算没有面积、实际几何读口、完整世界线来源或实际制备装置。dressing是表示／区域代数构造，不是免费实施一次改动能源的物理门。它不能把原1270的齐次量子块换名为局部场，不能从任意状态的小概率误差推出熵误差，也不能作为新TR或CT采用。

可复用的净增候选是精确余项及正能预算；真正后继问题是同一CT选后局部模型能否给上述共同迹、材料能源、观察者及约束的物理字典。通过此有限内部计算本身不意味着JSS原III型问题已被控制，更不意味着几何熵已获热身份。

## 5. 有限诊断与证据边界

主程序检查三个能层、对易态、纯相干态及非对易混态，保完整纯化，另从R边直接核条件熵；30例的主恒等式残差约$4.5\times10^{-14}$。硬窗程序另外直接检查截断状态的熵、能源与余项。两阶求积比较只作数值稳定性诊断，不能冒充严格尾积分或舍入误差证书；所有普遍不等式由上面的解析证明承担。
