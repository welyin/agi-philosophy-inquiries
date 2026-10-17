# 1245：同一杆、钟与记录材料的全三维约束接续

2026-10-15。接[准入](admission.md)。本轮是CT猜想下的**同一有效模型连接**：把1242的完整相互作用材料放入1188待选择的正则类别，证明选后全三维约束保持及协变动力学，再继承1242的非空有限记录。CT未采用，Jacobson所需独立热—几何身份也未由这里取得。

## 1. 问题、复用及新增范围

直接复用[351](../../archive_342_369/research_note_351.md)的纯度规括号与逆Legendre，[352](../../archive_342_369/research_note_352.md)的常矩阵多标量及局域混合项消去方法，[1188](../../archive_1161_1200/1188/proof.md)的CT五片系数选择，[1236](../1236/proof.md)的移动内部钟正则身份。不能将这些重新计为发现。

[1242](../1242/proof.md)实际给了杆、标量、正作用量校时钟和变质量记录器的协变作用、全部来源及一个有限FLRW反馈记录，但直接输入了Einstein分支。本轮净增是这些**相互耦合的同一材料**在一般空间依赖下的联合Hamiltonian、全部交叉项、选后总约束闭合，以及它们与1242完整初值/预测的同一性。不是先作均匀约化，再由一个常微分约束冒领全三维保持。

### 1.1 明示的有效类

取紧致无边界三环面，或全部变分有紧支撑的局部域。已有平滑四维Lorentz有效几何及共同ADM字典

$$
ds^2=-N^2dt^2+h_{ij}(dx^i+\beta^idt)(dx^j+\beta^jdt),
\qquad K_{ij}=\frac{\dot h_{ij}-{\cal L}_\beta h_{ij}}{2N}.
$$

几何正则变量仅为$(h_{ij},\pi^{ij})$；公开的待选引力密度为

$$
\mathcal H_g=\frac a{\sqrt h}(\pi^{ij}\pi_{ij}-\lambda\pi^2)
-b\sqrt h\,{}^{(3)}R+2\Lambda_H\sqrt h,
\qquad a,b>0.
\tag{1}
$$

不先取Einstein系数。空间协变动量生成元取$D_g[\xi]=\int\pi^{ij}{\cal L}_\xi h_{ij}$。把$(1)$作为独立动态度规的局域低阶有效类别、采用该正则结构以及允许lapse/shift乘子变分，都是物理建模输入；不是普通报告重排的推论。本文不排除类外额外场、非最小混合或高阶项。

材料完全保持1242：

$$
\mathcal L_f=-\frac12(\nabla p)^2-\frac12p^2
-\frac12F(p)\sum_{A=1}^3(\nabla\phi^A)^2,
\quad F(p)=1+\frac p5>0,
\tag{2}
$$

$$
M_c(J,K)=\frac7{10}+J+K,
\qquad M_r(p,I)=\frac1{10}+\frac I{10}+\frac{Ip}{10}>0.
\tag{3}
$$

钟内部辛项为$Jd\vartheta+Kd\chi$，记录器为$Id\theta$；$I,J,K\ge0$。为使整个非负$I$部门都有下界，本文进一步选$p>-1$的有效域，即$\omega+\lambda p=(1+p)/10>0$；1242的有限管道已满足它。在正作用量轨道使用action-angle局部坐标；$K=0$存能支改用正平面振子笛卡尔坐标，不把原点的角当物理读数。

两种worldline流体以固定有限正标签测度$d\mu_s(z)=n_s(z)d^3z$描述；$X_s(t,z)$在声明时间窗内光滑、无交叉焦散，空间Jacobian有正下界，单粒子质量及动量光滑有界。标签不是外部控制或新增能量，代表原有效连续材料的守恒个体数。该条件使下文delta推前得到光滑有限来源；**有限点粒子的奇异自引力不在本文范围内**。不同流体可相互穿越，无须两个种类的映射合起来一一。全部正则配对用同一个$d\mu_s$：$\{X^i(z),P_j(z')\}=\delta^i_j\delta_\mu(z,z')$，内部作用量/角亦然；$\int d\mu(z')\delta_\mu(z,z')f(z')=f(z)$。若改用$d^3z$，则PB含$1/n_s(z)$。代码含$n_s\ne1$的核验，不能Hamiltonian带密度却在PB漏掉它。

只研究以上已声明有效材料。实际准备控制器、信号接收器、全部未知量子输入及读出装置未加入；将来加入时必须进入同一Hamiltonian/来源，不能领取本轮已有的完整性。

## 2. 材料Legendre变换：存在相互作用而不含几何动量

记$P_p,\Pi_A$为场的密度动量，则

$$
P_p=\frac{\sqrt h}{N}(\dot p-\beta^i\partial_ip),\qquad
\Pi_A=\frac{\sqrt hF}{N}(\dot\phi^A-\beta^i\partial_i\phi^A).
\tag{4}
$$

$F>0,N>0$保证这份物质Legendre映射可逆，得到

$$
\begin{aligned}
\mathcal H_f={}&\frac{P_p^2}{2\sqrt h}
+\sum_A\frac{\Pi_A^2}{2\sqrt hF}
+\frac{\sqrt h}{2}\left[h^{ij}\partial_ip\partial_jp+p^2
+F\sum_Ah^{ij}\partial_i\phi^A\partial_j\phi^A\right],\\
\mathcal D_{f,i}={}&P_p\partial_ip+\sum_A\Pi_A\partial_i\phi^A.
\end{aligned}
\tag{5}
$$

对每个标签，沿用1236的移动Legendre恒等式，令$P_{s,i}$为世界线机械动量，

$$
E_s=\sqrt{h^{ij}(X_s)P_{s,i}P_{s,j}+M_s^2},
\quad H_s[N]=\int d\mu_s\,N(X_s)E_s,
\quad D_s[\xi]=-\int d\mu_s\,\xi^i(X_s)P_{s,i}.
\tag{6}
$$

同一正则作用是

$$
S=\int dt\left[\int\pi^{ij}\dot h_{ij}+P_p\dot p+\sum_A\Pi_A\dot\phi^A
+\sum_s\int d\mu_s(P_{s,i}\dot X_s^i+\Theta_s)
-H_{\rm tot}[N]-D_{\rm tot}[\beta]\right],
\tag{7}
$$

其中$\Theta_c=J\dot\vartheta+K\dot\chi$，$\Theta_r=I\dot\theta$，总$H=H_g+H_f+\sum_sH_s$，总$D=D_g+D_f+\sum_sD_s$。式(6)的负号与$dx+\beta dt$一致：$\dot X^i=Nh^{ij}P_j/E-\beta^i$；标量场的shift项却是$+\beta^iP_p\partial_ip$。不能把两者强制同号。

关键身份均来自(7)，不是另行指定：

$$
\frac{d\theta}{dt}=N\frac{M_r}{E_r}\frac{1+p(X_r)}{10},
\quad \frac{dI}{dt}=0,
\quad \frac{d\vartheta}{dt}=N\frac{M_c}{E_c},
\quad\frac{dJ}{dt}=\frac{dK}{dt}=0.
\tag{8}
$$

由于$d\tau_s/dt=NM_s/E_s$，这正是1242的proper-time相位律。记录器对$p$的力项为

$$
\left.\dot P_p(x)\right|_r
=-\int d\mu_r\,N(X_r)\frac{M_rM_{r,p}}{E_r}\delta^{(3)}(x-X_r),
\qquad M_{r,p}=I/10.
\tag{9}
$$

全部$\mathcal H_m$不含$\pi^{ij}$、$\dot h$或$\partial h$；其$h$依赖代数且与lapse无导数。这核实了1188所需最小物质资格，范围比单自由钟更完整，也不适用于1238的度规速度混合材料。

## 3. 同一法向来源、交互与协变守恒

对$N,\beta,h$作变分，(5)—(6)给同一个(2)—(3)协变材料作用的Hilbert应力：

$$
\mathcal H_m=\sqrt h\,T_{nn},\qquad
\mathcal D_{m,i}=-\sqrt h\,j_i,\qquad
j_i=-T(n,\partial_i).
\tag{10}
$$

例如一流体的密度和空间应力为

$$
\rho_s=\frac1{\sqrt h}\int d\mu_s E_s\delta^{(3)}(x-X_s),\quad
S_s^{ij}=\frac1{\sqrt h}\int d\mu_s\frac{P_s^iP_s^j}{E_s}\delta^{(3)}(x-X_s).
\tag{11}
$$

材料场的$\mathcal D_f$正号与(10)一致，因为$T(n,\partial_i)=P_p\partial_ip/\sqrt h+\sum_A\Pi_A\partial_i\phi^A/\sqrt h$。内部辛项无度规依赖；质量中的$Ip/10$已经计入来源，不能再把它加一次或删去。

令$\nu_s$为由标签守恒流得到的proper数密度。材料方程给

$$
\Box p-p-\frac{F_p}{2}\sum_A(\nabla\phi^A)^2-\nu_rM_{r,p}=0,
\quad \nabla_\mu(F\nabla^\mu\phi^A)=0.
\tag{12}
$$

记录器的变质量运动方程是$M_ra_r^\nu=-M_{r,p}(g^{\nu\alpha}+U_r^\nu U_r^\alpha)\nabla_\alpha p$。世界线位置变分固定该标签的内部$I$；若不同标签的$I$不同，不能用Eulerian总$\nabla M_r$额外加入$M_{r,I}\nabla I$力。钟质量虽可在标签间不同，$J,K$逐线守恒时也不由其空间分布产生额外自由落体力。因此

$$
\begin{aligned}
\nabla_\mu T_p^{\mu\nu}&=
\left[\tfrac12F_p\sum_A(\nabla\phi^A)^2+\nu_rM_{r,p}\right]\nabla^\nu p,\\
\nabla_\mu T_\phi^{\mu\nu}&=-\tfrac12F_p\sum_A(\nabla\phi^A)^2\nabla^\nu p,\\
\nabla_\mu T_r^{\mu\nu}&=-\nu_rM_{r,p}\nabla^\nu p,\qquad
\nabla_\mu T_c^{\mu\nu}=0.
\end{aligned}
\tag{13}
$$

总$T$协变守恒，不要求相互交换部门各自守恒。这些全三维身份并未由$T=G/\kappa$倒定义。

## 4. 全三维总括号，不先限制为均匀片

取任意光滑测试$N,M$，$v^i=h^{ij}(N\partial_jM-M\partial_jN)$，采用351的Poisson符号。

### 4.1 非常数$F$的场部门

把$q^a=(p,\phi^1,\phi^2,\phi^3)$，令$B_{ab}(q)=\operatorname{diag}(1,F,F,F)$，$A=B^{-1}$。场Hamiltonian正是

$$
\mathcal H_f=\frac{P_aA^{ab}P_b}{2\sqrt h}
+\frac{\sqrt h}{2}h^{ij}\partial_iq^aB_{ab}\partial_jq^b+\sqrt hV(q).
$$

泛函变分中，$\partial_aA,\partial_aB,V_a$项都正比$N$。它们与另一$\delta H/\delta P$相乘后是$NM$，在反对称括号中抵消；散度中的$\partial_iB$及$\partial h$也同样抵消。余下的仅是lapse梯度项，恰因$AB=I$得到

$$
\{H_f[N],H_f[M]\}=D_f[v].
\tag{14}
$$

这是352方法对该具体$F(p)$的应用，不是假定$F_p=0$。

### 4.2 世界线与场—世界线交叉项

固定背景场$h,p$及该标签的$I,J,K$，令$\partial_iE_s$包含$h,p$在$X_s$的完整空间变化，包括$M_{r,p}\partial_ip$，不对其他标签的内部动作分布求导。每标签的机械括号为

$$
\partial_{X^i}(NE)\,\partial_{P_i}(ME)
-\partial_{P_i}(NE)\,\partial_{X^i}(ME)
=P^i(M\partial_iN-N\partial_iM).
\tag{15}
$$

所有$NM\partial_iE$严格消去。质量虽依赖$I,p$，Hamiltonian不含内部角，故内部正则对不产生额外项。不同种类的机械/内部正则对独立。因此$\sum_s\{H_s[N],H_s[M]\}=\sum_sD_s[v]$。

场—记录器交叉项不能预先忽略。唯一非零的正则场配对来自$p,P_p$，并且

$$
\{H_f[N],H_r[M]\}
=-\int d\mu_r\,\frac{N(X_r)M(X_r)P_p(X_r)}{\sqrt{h(X_r)}}
\frac{M_rM_{r,p}}{E_r}.
\tag{16}
$$

交换$N,M$后的另一个交叉项是相反数；所以二者相加为零。它们逐项非零，不能用删除记录器反作用来解释消去。综合得

$$
\boxed{\{H_m[N],H_m[M]\}=D_m[v].}
\tag{17}
$$

引力—物质交叉括号中，$\delta H_g[N]/\delta\pi$正比$N$，而$\delta H_m[M]/\delta h$正比$M$且不含其导数，故反对称后也精确抵消。351的纯度规结果因此成为该联合对象的精确式

$$
\{H_{\rm tot}[N],H_{\rm tot}[M]\}
=abD_g[v]+D_m[v]
+2ab(1-2\lambda)\int v^i\nabla_i\pi.
\tag{18}
$$

这里$\pi$是密度，$\nabla_i\pi=\sqrt h\partial_i(\pi/\sqrt h)$。数值局部jet只核公式，(14)—(18)的全光滑场量词由解析变分及边界条件承担；不宣称任何有限格点模型具有精确连续HDA。

## 5. CT条件选后：全约束保持与Einstein方程

**条件选择：**若该共同对象满足1188的CT五切片、同一源与几何速度字典，则直接复用其结果

$$
\lambda=\frac12,\qquad ab=1.
\tag{19}
$$

CT仍是用户选择保留的猜想。本文证明其在此材料类别内有非空选后共同实例，不证明该猜想来自现行公理或任意内部酉。

这里使用CT的精确有效身份，而不是把1235的固定误差系数区间当精确等式。若只有有限切片预算，(18)一般仍有小的额外项；其对具体记录的稳定性/截断界须另行给出，本文不从“系数接近”自动领取1242的精确轨道或1243的场余项预算。

式(18)此时成为$\{H[N],H[M]\}=D[v]$。空间协变性给$\{D[\xi],H[N]\}=H[{\cal L}_\xi N]$及$\{D[\xi],D[\eta]\}=D[[\xi,\eta]]$。于是对总生成元$H[N]+D[\beta]$，任意初始满足的全部Hamiltonian与动量约束在平滑解存在的声明窗口内保持：每个约束的时间导数是同一约束的线性组合。这里没有把“系数满足(19)”偷换成“任意基片自动满足约束”；基片需单独验收，见第6节。

选后的动量映射及逆Legendre为

$$
\pi^{ij}=\frac{\sqrt h}{a}(K^{ij}-Kh^{ij}),\qquad
S_g=\frac1a\int N\sqrt h
\left[K_{ij}K^{ij}-K^2+{}^{(3)}R-2a\Lambda_H\right]dtd^3x.
\tag{20}
$$

给定上述时空重建和相应时间边界条件，这与$\frac1a\int\sqrt{-g}({}^{(4)}R-2\Lambda_{\rm phys})$相差ADM边界项，$\Lambda_{\rm phys}=a\Lambda_H$。同一材料的逆Legendre回到(2)—(3)，所以完整变分得到

$$
\boxed{G_{\mu\nu}+\Lambda_{\rm phys}g_{\mu\nu}
=\frac a2 T_{\mu\nu},\qquad \kappa=\frac a2=8\pi G.}
\tag{21}
$$

这次(21)是**一般父类别＋CT条件＋既有几何/正则输入**的选后结果，未将Einstein系数先写入(1)。不过其认知来源、低阶类别的唯一性、独立几何热熵、量子父过程和现实读口仍未得到。正$a$来自类别正性输入；$a$及$\Lambda_H$未由CT计算。

## 6. 与1242的非空基片和有限预测完全同一

取$a=6,b=1/6,\Lambda_H=0$，则$\kappa=3$、$8\pi G/3=1$、Einstein作用系数$1/6$。保持1242原三环面、compact杆绕数、两流体标签密度和初值，不增加新的钟或无源记录器。$\phi^A=kx^A$，$k^2=1/10$；初始$h_{ij}=\delta_{ij}$、$p=\dot p=0$，场动量全零，世界线共动且机械动量零。相位/内部作用量、质量和密度同原文，因此

$$
\rho_+(0)=\frac3{20}+\frac{51}{50}=\frac{117}{100},\qquad
\rho_-(0)=\frac3{20}+\frac{41}{50}=\frac{97}{100}.
$$

取$K_{ij}=H_\pm h_{ij}$，$H_\pm^2=\rho_\pm$，则选后法向约束为

$$
\frac16(-6H_\pm^2)+\rho_\pm=0.
\tag{22}
$$

所有动量约束也为零：几何动量空间恒定、全部场动量和机械动量零。故是同一全三维理论的合格基片；场的非零绕数梯度和正物质能量没有被删除。均匀各向同性子族满足全部空间演化方程，原协变材料方程及(21)逐项约化回1242的$u,H,p,v,z$系统。既有有限窗解析存在/区间证书可以直接继承，不需重新拟合或重算不同参数的轨道。

冻结[1242结果](../1242/results.json)于独立钟$T=1/5$给

$$
\Delta z\in
[8.4764435848022704888447232043,\;
8.4764435848325279284209288119]\times10^{-8}.
\tag{23}
$$

扣除每支$10^{-8}$的**完整组合读口合同**后，下界为$6.4764435848022704888447232043\times10^{-8}>0$。单杆梯度范数的两支分离大于$0.0020536142067384$。这些是原结果的只读继承，不新计数值科学组，也不把合同提升为实际仪器实现。原冻结几何删除支给$\Delta z=0$，仍不是另一份满足Einstein初始约束的解。

该非空正例满足选后协变动力学，因而在其所有合法切片保持同一约束/来源身份；这说明CT条件在该有效类内不是空条件。它不证明所有认知系统必须采用CT，不排除类外理论，也不提供同一实际内部装置在两支之间的准备转换。

## 7. 删除检验与剩余自由

1. **删除CT：**原(1)中的$\lambda,a b$仍可独立取值；351/1188的原异常和五片见证直接复用。本文不把一般$\lambda$类别都宣布不一致；新增CMC条件、限制lapse或类外模式属于不同范围。
2. **只删记录来源而保原读数：**若从记录器Hamiltonian把$Ip/10$删除，却仍要求$d\theta/d\tau=(1+p)/10$，在静止非零$p$处就有$\partial_I M=1/10$与要求差$p/10$。不存在本文同一Hamiltonian身份；完整force/source交换同样失配。若连相位律也同步改掉，则是另一个自洽材料，已经失去1242原记录菜单，不能叫同一对象。1242已给连续方程和约束残差，直接继承而不新计。
3. **自由保留：**$a>0$和宇宙学项、材料及其参数、态与初边界仍独立登记；同一切片代数不计算这些数值。CT未选额外模式/高阶项，也未证明独立热参考的几何面积熵。目标三项联合验收及阶段论文尚未完成。

## 8. 可复算范围

[check.py](check.py)有两组新的作者检查：一是带$F_p,M_p$、非零机械动量及非对角$h$的正则/局部jet括号与单向删除诊断；二是选后全张量Legendre恒等式、法向来源归一及1242两支全基片约束。全部证书算术使用有理数；局部jet不承担全场定理。另调用1242既有入口默认只读复算，核旧结果及哈希；该继承不新增科学组。保存仅限首次`--write`排他创建，之后为只读比对。

停止线：得到这一共同模型连接后，不继续扩大经典杆/转子装置或扫描Hamiltonian系数。下一项仍应处理CT机制来源或同一参考/区域的独立热—几何响应，不能用选后Einstein作用反定义面积熵再回证自己。
