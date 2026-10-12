# 1121证明：局部CAR混味操作的完整平移反冲

2026-10-11。接[1120下一步](../1120/NEXT.md)和[来源去重](../_shared/1120_1156/source_bridge_audit.md)。本轮只审一项权限连接：879／881的局部偶CAR操作若接入D-P的完整四动量账，能否把“同质量、同自由味谱”当成无动量来源的理由。答案在以下声明类别内是否定的，并有与区域半径相关的算符下界。

新增认知公理0、科学校准0。四原则、CO、SR及D-P均保原范围。本轮新增的是明确物理分支和必要反冲证书，不声称已经构造实际局部控制器、原CO区域字典、全部场论或引力。

## 1. 父材料、局域性与新增物理输入

[853§2—4](../../archive_819_853/research_note_853.md)已有正质量、内部规范单态、两味简并的Dirac接收材料及光滑Cauchy模式；[879§3—6](../../archive_873_893/research_note_879.md)给原CAR区域呈现和局部偶代数；[881§4—7](../../archive_873_893/research_note_881.md)保有限局部模式、原关联及来源菜单，但未认证其实际自主控制或物理能源统一界。

本轮**另行选定**这类材料的3+1平直、常质量、自由Dirac分支。它不是把879的既定曲背景改写后仍声称同一背景，也不把853原形式交互当作已收敛的有限强度仪器。共同自由平移、正能源真空及场的材料意义是明确物理输入。

取$\hbar=c=1$，单味Cauchy空间$\mathcal K_0=L^2(\mathbb R^3,\mathbb C^4)$，两味空间$\mathcal K=\mathcal K_0\otimes\mathbb C^2$。记

$$
h_D=-i\alpha^j\partial_j+\beta m,\qquad m>0,\qquad p_j=-i\partial_j,
\qquad S=1_{h_D<0}\otimes I_2.
\tag{1}
$$

$\alpha^j,\beta$是标准Hermitian Dirac矩阵，味指标由$I_2$作用。$h_D$有正负两支；**不把满$L^2$上的$d\Gamma(h_D)$称为正能物理Hamiltonian。** 采用投影$S$的通常准自由Dirac真空表示。其粒子／反粒子Fock实现的物理Hamiltonian为正谱部分的第二量子化加负谱部分空穴的正能源，因而$\mathbb H\ge0$；物理空间动量记$\mathbb P_j$，真空动量为零。

更具体，令$\mathcal K_\pm=1_{\pm h_D>0}\mathcal K$，Fock一粒子空间为$\mathcal K_+\oplus\overline{\mathcal K_-}$，则

$$
\mathbb H=d\Gamma(h_D|_{\mathcal K_+})+
d\Gamma(-\overline{h_D}|_{\overline{\mathcal K_-}}),\qquad
\mathbb P_j=d\Gamma(p_j|_{\mathcal K_+})+
d\Gamma(-\overline{p_j}|_{\overline{\mathcal K_-}}).
\tag{2}
$$

Cauchy场$c(f)$在此表示中由粒子湮灭及反粒子产生两部分组成，满足$\{c(f),c(g)^\dagger\}=\langle f,g\rangle I$，按$f$反线性约定。紧支$f$不是纯正能源粒子波函数，以下证明不使用这种不存在的假设。平移由真实物理$\mathbb P_j$实现：

$$
\alpha_a^{(j)}(c(f))=c(e^{ia p_j}f)
=e^{ia\mathbb P_j}c(f)e^{-ia\mathbb P_j}.
\tag{3}
$$

选归一$f\in C_c^\infty(B_R,\mathbb C^4)$，两味为$f\otimes e_1,f\otimes e_2$。它们在同一Cauchy球及其因果域内的CAR中。局部偶元与类空区的全部局部偶元对易；这里的局域性是场代数支持，不是已证明一个空间延展模式门能在任意短时自主完成。

## 2. 明确的局部混味操作

对有限秩一粒子算子$A$，用一个包含其像与伴随像的有限正交模式基定义有界CAR元

$$
q(A)=\sum_{ab}A_{ab}c_a^\dagger c_b.
\tag{4}
$$

该定义与基选择无关，$q$在有限秩算子上满足$[q(A),q(B)]=q([A,B])$，无中心项。取$P_f=|f\rangle\langle f|$、$X=\begin{pmatrix}0&1\\1&0\end{pmatrix}$及实$\theta$，定义

$$
A=P_f\otimes X,\qquad U_\theta=e^{-i\theta q(A)},\qquad
E=e^{-i\theta X},\qquad
V=I+P_f\otimes(E-I).
\tag{5}
$$

$U_\theta$是有界局部偶幺正，保全局Dirac电荷及宇称；在所选两味模式上作同一混味，在其CAR正交补上取恒等。它是可检验的**目标操作／新增操作合同**，不是853旧$p\bar R\sigma_yR$顶点已经在有限强度实现它的证明。$\theta\in2\pi\mathbb Z$时它是恒等，其余角度在完整CAR上非平凡。

其中$c_i^\dagger(f)c_k(f)$是$B_R$内的双点核，并不等于已给出一个点状Hamiltonian密度或零时长可执行门。其代数支撑与有限持续时间的实施权限必须分别验收；把$f$形式上替换成$f(x-Q)$也尚未制造完整定位控制器。

## 3. 从实际平移导数定义动量缺陷

因为$f\in D(p_j)$，$a\mapsto\alpha_a^{(j)}(U_\theta)$在算符范数下C¹：只需有限个模式向量在$L^2$下可微，再用CAR的$\|c(f)\|=\|f\|$及有界指数的导数。故$U_\theta$保持$D(\mathbb P_j)$，且以下差在该域上有唯一有界延拓：

$$
\Delta\mathbb P_j:=U_\theta^\dagger\mathbb P_jU_\theta-\mathbb P_j
=-iU_\theta^\dagger
\left.\frac{d}{da}\alpha_a^{(j)}(U_\theta)\right|_{a=0}
=q(B_j),\qquad B_j=V^\dagger p_jV-p_j.
\tag{6}
$$

最后一步可由$[\mathbb P_j,q(A)]=q([p_j,A])$和有限CAR对易恒等式积分得到；$[p_j,A]$有有限秩及有界延拓，故没有对无界算子作无依据的级数运算。也可先在共同光滑核上微分式(3)，再用有界延拓。这里$B_j$是Cauchy空间上的有限秩算子，$q(B_j)$是实际正能场表示内的有界元；二者不是物理全Hamiltonian的负能Fock替代。

令

$$
\mu_j=\langle f,p_jf\rangle\in\mathbb R,\qquad
r_j=(p_j-\mu_j)f,\qquad \sigma_j=\|r_j\|.
\tag{7}
$$

若$\sigma_j=0$，式(6)为零。否则在$(f\otimes\mathbb C^2)\oplus((r_j/\sigma_j)\otimes\mathbb C^2)$上，

$$
B_j=\sigma_j
\begin{pmatrix}0&E^\dagger-I\\E-I&0\end{pmatrix};
\tag{8}
$$

其余正交补上为零。$p_jr_j$或更远模式不进入此差，不能因$p_j$不保持四模子空间而否定(8)：只有含$P_f$的交叉项留下，$p_jf=\mu_j f+r_j$已完整给出它们。直接展开$V^\dagger p_jV-p_j$中$f$—$f$块为零，正是味简并下$E^\dagger E=I$的结果。

## 4. 完整CAR范数与正常序

置$s=|\sin(\theta/2)|$、$a_j=2\sigma_js$。因为$(E-I)^\dagger(E-I)=4s^2I_2$，$B_j$的非零四维谱为

$$
\operatorname{spec}B_j=(a_j,a_j,-a_j,-a_j),\qquad
\|B_j\|=a_j.
\tag{9}
$$

四个CAR模式生成$M_{16}$。在其标准矩阵表示中，$q(B_j)$的本征值是任意占据集合的上述本征值之和，因此

$$
\operatorname{spec}q(B_j)=(-2a_j,-a_j,0,a_j,2a_j),
\quad\text{重数为 }(1,4,6,4,1),
\quad
\boxed{\ \|\Delta\mathbb P_j\|=4\sigma_j|\sin(\theta/2)|\ }.
\tag{10}
$$

端点$\pm2a_j$在双占据、偶宇称部门取得，不是只在奇部门的形式值。$M_{16}$的任何含幺$*$表示均为忠实表示，因而式(10)在原正能场表示中保持。**一粒子矩阵范数(9)比完整CAR结果少一半，不能混用。** 张量任意被动参考不改该算符范数与身份。

若采用Dirac真空正常序，则

$$
q(B_j)=:q(B_j):_S+\operatorname{Tr}(SB_j)I.
\tag{11}
$$

迹有定义，因为$B_j$有限秩。它一般不为零。若用正常序的$q(A)$定义门，门只相差一个全局相位，(6)不变；但把(11)的第二项从**动量差**删掉，会改变实际平移导数。物理$\mathbb P_j$本身加共同常数也不会改变(6)。本轮不宣称真空能源或物理动量基准可任意重标。

例：仅作四模代数校准，取空间压缩$p=\begin{pmatrix}0&1\\1&2\end{pmatrix}$、其负谱投影$S_0$，再乘味单位阵，$\theta=\pi/2$。此时$\operatorname{Tr}((S_0\otimes I_2)B)=\sqrt2\ne0$。这是正常序记账的有限表示检查，不把此任意二维$p$认成真实连续Dirac谱或真空的数值逼近。

## 5. 有限区域给明确的非零下界

对所有归一$f\in H_0^1(B_R,\mathbb C^4)$，令$g(x)=e^{-i\boldsymbol\mu\cdot x}f(x)$，其中$\mu_j$沿(7)。相位不改支撑或范数，故

$$
\sum_{j=1}^3\sigma_j^2=\int_{B_R}|\nabla g|^2d^3x
\ge\frac{\pi^2}{R^2}\int_{B_R}|g|^2d^3x
=\frac{\pi^2}{R^2}.
\tag{12}
$$

为自足说明常数：先取光滑紧支标量分量，用球坐标丢掉非负角向梯度；令$w(r,\omega)=rg(r,\omega)$，径向分部积分给$\int_0^Rr^2|\partial_rg|^2dr=\int_0^R|\partial_rw|^2dr$。$w$在0及R为零，一维Dirichlet正弦展开给下界$\pi^2/R^2$。逐角、逐旋量分量积分，再用$H_0^1$稠密性得到(12)。

结合(10)，

$$
\left(\sum_j\|\Delta\mathbb P_j\|^2\right)^{1/2}
\ge\frac{4\pi}{R}s,
\qquad
\max_j\|\Delta\mathbb P_j\|
\ge\frac{4\pi}{\sqrt3R}s.
\tag{13}
$$

左边是三个分别定义的算符范数的Euclidean组合，不是宣称三个非对易差能在同一个态上同时达到极值。固定$R$及非平凡$\theta$时，至少一方向不守恒；两味质量简并没有消除局部操作的平移来源。一般动能、能源消耗、热量、最优制备成本或GR应力没有被(13)确定。

定理对$C_c^\infty$原材料模式成立，(6)—(13)只需相关$H^1$条件即可连续延拓。数值采用

$$
f_R(x)=C_R(1-|x|^2/R^2)_+^2\chi,\quad
C_R^2=\frac{3465}{512\pi R^3},\quad \|\chi\|=1,
\quad \mu_j=0,
\quad \sigma_j^2=\frac{11}{3R^2}.
\tag{14}
$$

此函数属$H_0^1\cap H^2$，**不是$C_c^\infty$**。径向归一积分为$4\pi R^3(128/3465)$，梯度积分在未归一时为$4\pi R(128/315)$，比为$11/R^2$；各方向相等。它可在$H^1$中被球内光滑紧支函数逼近，归一后$\mu_j,\sigma_j$及有限秩CAR动量差连续。数值例仅核公式，不把有限阶光滑样例当成全部旧Hadamard准备条件已经认证。

## 6. 与D-P相接时，控制器究竟需要承担什么

以下是**必要条件**，不是实施存在性证明。令$W$是包含原材料和全部控制／支撑／参考的完整过程，完整动量$\mathbb P_j^{\rm tot}=\mathbb P_j+\mathbb P_j^C$满足守恒。则在共同形式域上

$$
W^\dagger\mathbb P_j^CW-\mathbb P_j^C
=-\bigl(W^\dagger\mathbb P_jW-\mathbb P_j\bigr).
\tag{15}
$$

若某个输入等距$J$和物理实现合同还保证，在其全部声明输入上，材料动量变化的压缩恰为(6)，则控制器变化的同一压缩必须为$-\Delta\mathbb P_j$。该压缩是已知有界算符时，它的范数必须等于(10)。所以不能同时要求本全CAR混味任务、完整守恒及控制器在该输入合同上没有相应动量变化。

量词必须保留：

- 完整CAR算符检验覆盖四模代数，不等于879／881旧的固定准备菜单已能制备所有极值态。即使物理表示忠实，达到范数的准备成本与统一能源预算仍需另核。
- 若只操作一个受限码或固定输入，压缩$J^\dagger\Delta\mathbb P_jJ$可能变小或为零；不能把全代数下界无条件移给它。完成原未知输入及参考的任务身份需要另给$J$。
- 精确实现一个目标系统幺正、只在某些状态上复现其动量变化、带资源关联的完整仪器，是不同合同；本轮不声称任意有限控制器都能精确实现前者。参考依赖或记录纠缠不能被当成免费辅助复原。
- 近似通道的普通diamond误差不控制无界$\mathbb P_j$。若需要近似版(15)，须直接给动量差的有界压缩误差或加权矩证书，不能只用末端概率接近。

通过这些限定，本轮连接的是旧局部操作权限与D-P之间的一份**必要反冲账**。它没有把场的Cauchy球自动识别为原CO实际端点，也没有证明局部模式门的触发、处理时间、来源族或应力响应。

## 7. 固定辅助、全输入精确酉合同的更强边界

这是成熟的守恒和无信息泄漏论证在当前对象上的应用，不作为新的认知反证。设初始辅助为与未知场输入独立的固定正常态$\eta_C$，总过程幺正$W$完整保$\mathbb P_j+\mathbb P_j^C$，并对**全部场输入态**精确实现酉通道$\operatorname{Ad}_{U_\theta}$。此处允许辅助最终改变，不要求辅助复原。为合法使用守恒期望，还要求初始辅助、以下共同辅助终态及被检验的稠密场输入域均有有限绝对一阶动量矩；$U_\theta$保持该场域。

先取纯辅助$|a\rangle$。纯输入$|\psi\rangle$的场边缘输出是纯态，故

$$
W(|\psi\rangle\otimes|a\rangle)
=U_\theta|\psi\rangle\otimes|b_\psi\rangle.
\tag{16}
$$

对两个线性独立输入及其叠加使用线性性，得同一个$|b\rangle$，相位可吸收到此共同向量；因此辅助终态与输入无关。一般正常混合态先纯化，令$W$对纯化参考取恒等；同一论证给整个环境终态固定，特别是实际辅助终态固定。纯化参考不参与守恒动量账。

在上述有限矩输入上，完整守恒于是要求

$$
\langle\psi|\Delta\mathbb P_j|\psi\rangle
=\operatorname{Tr}(\eta_C\mathbb P_j^C)
-\operatorname{Tr}(\eta'_C\mathbb P_j^C)=k_j,
\tag{17}
$$

其中$k_j$与$\psi$无关。$\Delta\mathbb P_j$有界，稠密域上的此恒等式经极化及连续性给$\Delta\mathbb P_j=k_jI$。式(10)—(13)却保证非平凡$\theta$与有限球至少一方向有非标量正负谱，故这组**全输入、精确、固定输入独立辅助、完整守恒、有限初末矩**条件不能同时成立。

这里排除的不是有限有效阶段。有限输入码、带显式矩控制的近似门、输入相关资源或未满足有限矩条件的理想参考均不在此结论内；它与[1110](../../archive_1106_1119/research_note_1110.md)的有限资源近似完全相容。下一步仍需一份实际有限区域过程和定位反作用合同，而非宣布所有局部处理不可能。

## 8. 验证与停止线

[代码](check.py)在四模16维CAR中独立构造产生湮灭矩阵、模式混味酉、平移差及正常序检查；另用独立Gauss积分核(14)。同一算符张量被动参考保持，默认复算只读比较[结果](results.json)。解析证明承担全部光滑模式、球半径、角度及正能场表示的量词，有限矩阵不证明实际材料可制造。

后续优先检查一份具有明确有限输入与参考的局部控制器／定位材料如何满足(15)，或哪些受限码能实际消去本缺陷而保原任务。停止同一四模的角度／半径扫描，也不以继续优化这个自由分支代替整个阶段的局域过程和共同来源目标。见[下一步](NEXT.md)及[去重审计](selection_audit.md)。
