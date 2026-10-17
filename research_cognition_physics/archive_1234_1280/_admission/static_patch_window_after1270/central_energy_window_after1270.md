# 总能中央窗口：固定迹的精确条件化与未知态范围

2026-10-16。非编号方向核查；不预约1271，不改变目标、CT状态或冻结结果。只处理有限内作用模型，不认定它已成为dS/SdS的完整物理来源。

## 1. 去重与固定对象

复用[正能观察者引理](../positive_observer_entropy_after1270/finite_observer_entropy_candidate_after1270.md)及§3.1硬观察者能窗。旧§3.1只改变准备，已经提醒$q\le B$角与迹归一的区别；本次核的是**总能源$Y=q+H$处于$[0,N]$**。已回查[1189后内作用范围](../../../archive_1161_1200/_admission/after1189_reference_scope/README.md)§§3—4、[1245指定dressing](../marked_cocycle_test_after1245.md)§§2—3及1201—1218的有限参考/正能钟；没有找到下列同一固定迹的中央窗口、参考条件关联与全未知输入支持判据。不重算Gamma余项，也不把一般角归一视为新物理定理。

设$d<\infty$，$H=\sum_jE_jP_j$，$E_j\ge0$，$d_j=\operatorname{Tr}P_j$，$\beta>0$，$Z=\operatorname{Tr}e^{-\beta H}$，$\sigma=Z^{-1}e^{-\beta H}$。保原正能角和固定迹：

$$
P(y)=\boldsymbol1_{y\ge H},\qquad
\tau(A)=\frac{\beta}{Z}\int e^{-\beta y}\operatorname{Tr}A(y)\,dy,\qquad \tau(P)=1.
$$

工作态的Lebesgue密度与$\tau$密度分别是

$$
D(y)=f(y-H)\rho f(y-H)^*,\qquad
d_\rho(y)=\frac Z\beta e^{\beta y}D(y).
$$

观察者最初与材料及保留伙伴独立；材料可与任意伙伴关联，所需纯化共同保留。

## 2. 中央窗、迹与一般工作态

在原角代数内定义$C_N(y)=\boldsymbol1_{y\le N}P(y)$。它是中央投影：$y$为中心变量，$P(y)$是允许纤维的单位。下端$Y\ge0$由正能保证。对$N>E_{\min}$，

$$
t_N:=\tau(C_N)
=\frac1Z\sum_{E_j<N}d_j\left(e^{-\beta E_j}-e^{-\beta N}\right)\in(0,1),
\qquad
\tau_N(A)=\frac{\tau(A)}{t_N},\quad A=C_NAC_N.
\tag{1}
$$

当$N\ge E_{\max}$，有$t_N=1-d e^{-\beta N}/Z$。等于$N$的能层仅贡献零测度端点。这个常数由原$H,\beta,N$确定，未另选任意中心权重。

令$p_N=\tau(C_Nd_\rho)>0$，则条件工作态在$\tau_N$下的密度为

$$
d_{\rho,N}=\frac{t_N}{p_N}C_Nd_\rho.
\tag{2}
$$

在熵积分有定义时，利用角内对数直接得到

$$
S_{\tau_N}(d_{\rho,N})
=-\frac{\tau(C_Nd_\rho\log d_\rho)}{p_N}
+\log p_N-\log t_N.
\tag{3}
$$

若$0<p_N<1$，将外角也归一、记其条件熵为$S_{\rm out}$，各项有限时更有

$$
S_\tau(d_\rho)
=p_NS_{\tau_N}(d_{\rho,N})+(1-p_N)S_{\rm out}
-D_{\rm Ber}(p_N\Vert t_N).
\tag{4}
$$

仅有小失败率不保证连续观察者的熵差小，外角熵/能源尾仍需控制；这里不套有限维Fannes界。$p_N<1$时条件化通常改变原材料、伙伴和观察者，不能继续将旧$\rho$和旧独立$q$无更改地代入原引理。

## 3. 全未知工作态不被改变的精确条件

事先选择正常观察者工作准备$\operatorname{supp}f\subset[0,B]$并满足

$$
B+E_{\max}\le N.
\tag{5}
$$

在dressing前的完整空间，$C_N$对应$\boldsymbol1_{q+H\le N}$。式(5)使它在$\mathcal H_M\otimes\operatorname{span}\{f\}\otimes\mathcal H_R$上恒等。因此对任意未知$\rho_{MR}$均有$p_N=1$，新增窗口不改变完整态及其关联，也没有后选。这个结论仅指**新增窗口**；不把表示变换/区域限制自身说成无扰动的实际门。

于是

$$
d_{\rho,N}=t_Nd_\rho,\qquad
S_{\tau_N}=S_\tau-\log t_N
=-D(\rho\Vert\sigma)-\beta\langle q\rangle+h(q)+\log\beta
+\epsilon_f-\log t_N.
\tag{6}
$$

$\epsilon_f=I(Y:J\mid R)$仍按实际有界支持准备计算。右侧$\sigma$是保留原代数表达的旧参考，不意味着新角参考的材料边缘仍等于$\sigma$。

设$m=N-E_{\max}>0$，由$Z\ge d e^{-\beta E_{\max}}$，

$$
0<1-t_N\le e^{-\beta m},\qquad
0<-\log t_N\le-\log(1-e^{-\beta m})
\le\frac{e^{-\beta m}}{1-e^{-\beta m}}.
\tag{7}
$$

同一固定$H,\beta,N$下比较两份满足式(5)的工作态，这个常数精确相消。若参数或窗口改变，不能继续把它当作无关常数。

## 4. 原参考条件化后必须保留的关联

原参考$\sigma$与纯观察者波包$f_\beta(q)=\sqrt\beta e^{-\beta q/2}\boldsymbol1_{q\ge0}$共同实现$\tau$。中央条件化后的完整材料—观察者状态为

$$
\Omega_N=\frac1{Zt_N}\sum_j e^{-\beta E_j}P_j
\otimes|f_{\beta,j}\rangle\langle f_{\beta,j}|,\qquad
f_{\beta,j}=\boldsymbol1_{0\le q\le N-E_j}f_\beta,
\tag{8}
$$

其中$E_j\ge N$项为零。它在角代数上恰实现$\tau_N$，其密度为$C_N$。材料边缘是

$$
\sigma_N=\frac1{Zt_N}\sum_{E_j<N}
\left(e^{-\beta E_j}-e^{-\beta N}\right)P_j.
\tag{9}
$$

当$N>E_{\max}$，$\sigma_N=(e^{-\beta H}-e^{-\beta N}I)/(Zt_N)$；通常不是原$\sigma$。条件观察者波包的上端$N-E_j$随能级而变，参考不再独立。若原$\sigma$含保留纯化伙伴，必须在完整纯化上作用同一投影，不能删掉伙伴人口/关联的变化。

有限总能支持、原未截断指数观察者参考、仍独立的原$\sigma\otimes f_\beta$三者不能同时精确保留，因为指数参考有严格正的高能尾。这只否定三项同时坚持；允许式(8)的已计账参考则没有矛盾，更不改变式(5)的工作态。

记$u_N=1-t_N$。角参考作为原代数上的态延拓为零后，与原$\tau$的半范数距离恰为$u_N$。在完整参考纯化上，归一投影前后的纯态半迹距为$\sqrt{u_N}$。材料边缘另有

$$
\tfrac12\|\sigma_N-\sigma\|_1\le u_N,
\tag{10}
$$

因为$\sigma=t_N\sigma_N+u_N\sigma_{\rm out}$。这些是不同对象上的误差，不互换。

## 5. $q\le B$不是同一中央角

$$
Q_B(y)=\boldsymbol1_{0\le y-H\le B},\qquad
\tau(Q_B)=1-e^{-\beta B}.
\tag{11}
$$

$Q_B$一般非中央：某些$y$纤维仅保留部分材料能层，不与所有材料非对角算符对易。把它当作角会改变可访问代数；旧引理只选择工作波包，并未自动这样做。

式(5)下$Q_B\le C_N$，但一般不相等。若另选$Q_B\mathcal A Q_B$并归一，熵常数是$-\log(1-e^{-\beta B})$；其最大参考可由原$\sigma$与独立截断指数观察者实现。这与式(8)的中央总能条件参考不同，不能混用效果、常数或参考。

## 6. 有限代数核验及不能任意领取的近似

[代码](central_energy_window_after1270.py)和[结果](central_energy_window_after1270.json)只核本次有限窗口身份，不重跑Gamma：

1. $H=\operatorname{diag}(0,1)$、$\beta=\log2$、$N=3$给$Z=3/2$、$t_N=5/6$。材料参考$(2/3,1/3)$变为$(7/10,3/10)$。两能层条件下事件$q>2$的概率分别为$1/7,0$，故关联非零。
2. 观察者支持$[0,2]$时，全部未知材料及伙伴满足式(5)。等权纯相干材料加均匀波包可分段精确积分，给$S_\tau=\log(2\log2/3)$、$S_{\tau_N}=\log(4\log2/5)$，差恰为$-\log(5/6)$。数值重放只是解析式的诊断。
3. 同一准备若改$N=2$，成功率$p_N=3/4$，纯相干材料条件态变成$\left(\begin{smallmatrix}2/3&1/3\\1/3&1/3\end{smallmatrix}\right)$；原Bell伙伴人口也从$(1/2,1/2)$变为$(2/3,1/3)$。一般中央后选不能声称保未知态。
4. 若两能层间隙$\Delta>B$、工作态为等权纯相干叠加、$q\in[0,B]$，则$Y$的两条件支撑不相交，$\epsilon_f=\log2$。总能窗口与未知输入保持仍可相容，但该工作态的慢变熵近似不能任意精确；迹归一不消除此余项。此例不是所有有界支持波包或所有热几何模型的反证。

## 7. 几何候选的接口与停止线

若具体几何候选独立给出允许总能上界$N$，本页提供式(5)的工作支持、式(8)的参考关联和式(6)—(7)的熵校正。能否选择所需的有限$N$并符合几何源、尺度和实际观察者要求，仍须由该候选检验；不能反用本页定义它的几何。

本页没有证明$Y$就是某份dS/SdS的完整质量，没有计入新观察者的制造/支撑源，也不把参考条件化当免费实际制备。若实际后选，须保失败记录和资源。硬切正常波包仍不自动有有限共轭动量方差。

**结论：中央窗口不必破坏未知工作态；必须修改并保留关联的是参考，熵归一修正可精确控制。** 固定小$N$不能同时随意领取任意宽慢变观察者，这属于有限预算问题，不是全研究的无限能门槛。到此停止，不加编号、器件或几何假说。

