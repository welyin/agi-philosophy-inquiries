# 1131证明：两份正常源与原裸材料读口的梯度补全筛选

2026-10-11。接[1130](../1130/proof.md)的旧实际钟—连续交接接口。本轮固定1129的同一自治Hamiltonian、1130同一正窗口及完整旧钟，只增加两份事前列明的正常源准备。它们在同一个原B偶占据读口给出可分离概率区间；这对“用同源纯梯度取代部署法向”的候选提出有限误差下界。

新增认知公理0，科学校准0，阶段active。纯梯度的场重定义是成熟工具，本轮新量是**当前正常码、原裸端口、有限源准备菜单上的匹配判据**。没有构造正常量子梯度QFT，没有反证原来只含一份源的M1130合同，也没有判定全部协变共源方案失败。

## 1. 原对象及新增有限准备菜单

冻结[1129](../1129/proof.md)的完整正能Dirac Fock、质量零标量Fock、CAR源码J、全部正负谱、核、门函数、部署函数及

$$
H=H_F+H_\phi+V,\qquad
V=\lambda\int\eta(x)F(\phi_\epsilon(x)/A)J_\epsilon(x)\,d^3x,
\quad \lambda=\pi/2,\quad\epsilon=10^{-5},\quad A=10^{34}.
\tag{1}
$$

固定$T_*=5/4$。源均值的空间平台、材料因果管分别满足$5+T_*+\epsilon<7$、$2+T_*+\epsilon<4$。$F\ge0$支撑于$[1/4,3/4]$、$\int F=1$、$\|F'\|\le L=512$；$q_\epsilon\le8\epsilon^{-3}$、$C_*=8q_{\rm upper}\,524=3.3536\times10^{19}$。A是同一个固定模型参数，不随源准备重标。

**新增有效材料匹配条件M1131**：在原预部署菜单上允许两份正常相干源$a_s$，$s\in\{1,2\}$，其初始经典Cauchy均值是

$$
\phi_{cl,s}(0,x)=0,\qquad \pi_{cl,s}(0,x)=sA\beta_0(x).
\tag{2}
$$

其余准备、码、原端点及读口不变。s是可选的实际源准备标签，不是boost、坐标重参数、未知Q的函数，亦不开放任意源干预或中途控制。它是新增材料许可，不是新认知公理。相干位移仍为光滑紧支数据对应的正常Fock向量，$a_s\in D(H_\phi)$；原材料输入属于$D(H_F)$，实际源及完整余部均保留。

源能量明确改变：

$$
E_{\phi,s}(0)=\tfrac12s^2A^2\|\beta_0\|_2^2\le1073s^2A^2.
\tag{3}
$$

两准备不是同能量来源。本轮不优化巨大但有限的资源，也不把制造两相干态的仪器当作已完成。

## 2. 同自治H的全输出源比较

在$\operatorname{supp}\eta$及$0\le t\le T_*$，自由波动的有限传播给$\langle\phi_\epsilon(x)\rangle_{a_{s,t}}=sAt$，其中$a_{s,t}=e^{-itH_\phi}a_s$。相干位移不改变真空方差，故

$$
\|[\phi_\epsilon(x)/A-st]a_{s,t}\|\le3/(\epsilon A).
\tag{4}
$$

比较材料传播$\widehat U_s(t)$由完整有界CAR生成元$H_F+\lambda F(st)B_\epsilon$定义，$B_\epsilon=\int\eta J_\epsilon$。沿用1129的谱Lipschitz和Duhamel论证，得到

$$
\left\|\left[e^{-itH}-\widehat U_s(t)\otimes e^{-itH_\phi}\right]
(I\otimes|a_s\rangle)\right\|
\le \delta_{\rm src}(t):=\lambda C_*L\frac{3t}{\epsilon A}.
\tag{5}
$$

量词为全部材料输入以及任意被动参考，亦可把完整旧D放入参考；源初始独立且取所选同一$a_s$。这是实际源相位、后态、材料海和关联都保留的联合等距差，不是只比较一粒子脉冲或将源换成外场。误差不因s增大而增加；这来自(4)，不来自假设反作用为零。

1129的retarded场差(24)、动量差及能源式仍由同一H成立，统一界与s无关。例如

$$
\left\|[\phi_{\epsilon,H}(t,x)-\Phi_0(t,x)]/A\right\|
\le\frac{\lambda LC_*q_\epsilon t^2}{2A^2}.
\tag{6}
$$

初始源能量以(3)替换旧界，$\langle H_0\rangle_t\le268+1073s^2A^2+2\|V\|$，并保持$\Delta\langle H_F+H_\phi+V\rangle=0$。本轮没有把被比较自由源当成真实无反作用来源。

## 3. 连续模式误差及原B占据效果

写$u_t=e^{-iht}f$，$F_{Ai}=f\otimes e_A\otimes e_i$、$F_{Bi}=f\otimes e_B\otimes e_i$，固定终端实际模式

$$
G_i=u_{T_*}\otimes e_B\otimes e_i,\qquad N_B=c^\dagger(G_1)c(G_1).
\tag{7}
$$

$N_B$是原有限区域B的偶、有界占据effect，$0\le N_B\le I$。本轮不将它重新旋转成源依赖的新端口。理想比较的一粒子传播是

$$
U_s^\circ(t)=e^{-iht}e^{-i\alpha_s(t)\sigma_x},\qquad
\alpha_s(t)=\lambda\int_0^tF(s\tau)d\tau,\qquad
\alpha_s(T_*)=\lambda/s.
\tag{8}
$$

因两源都已穿过完整F支撑，末式精确成立。沿$u_t$的材料管内$\eta=1$，两次K平滑的误差不超过$2\epsilon\|\nabla f\|$。1129的模式Duhamel给任意内部AB方向

$$
\|(U_s(T_*)-U_s^\circ(T_*))(f\otimes v)\|
\le2\lambda\epsilon\|\nabla f\|/s=:d_m(s)\le40\lambda\epsilon/s.
\tag{9}
$$

特别取$v$使$U_s^\circ(T_*)(f\otimes v)=G_1$，再左乘$U_s^\dagger$，得到前像模式差也不超过$d_m(s)$。CAR的$\|n_v-n_w\|\le2\|v-w\|$因此给**全部正常材料输入**的读口误差$2d_m(s)$。只在两码模上的完整输出近似不是本步所需；没有把其余真空输出强行置为自由真空。

## 4. 正常海占据不能删：全Q/R的实际效果公式

源码仍为1129的

$$
\Omega_0=2(1-n_{A1})(1-n_{A2})\Omega,\qquad
J_i=c^\dagger(F_{Ai})\Omega_0,\qquad J^\dagger J=I_2.
\tag{10}
$$

极化$S=1_{h<0}$在AB及rail上分块，且$\langle f,Sf\rangle=1/2$。四个模式$F_{A1},F_{A2},F_{B1},F_{B2}$上的真空限制是协方差$\frac12I_4$的规范不变准自由态，即该16维CAR因子的最大混态。A两模空占据条件投影不改变B的$\frac14I_4$限制；源码只在A两模中加入一个未知双轨逻辑态。于是对所有逻辑矩阵元，

$$
J^\dagger n_{A1}J=|1\rangle\langle1|,\quad
J^\dagger n_{B1}J=\tfrac12I_2,\quad
J^\dagger c^\dagger(F_{A1})c(F_{B1})J=0.
\tag{11}
$$

这可以由有限模式Wick关系得到；代码另以完整16维CAR和B限制的纯化复核。它不是把连续真空换成空Fock真空，也不是声称整份无限海因子化。

理想前像$U_s^{\circ\dagger}G_1=\cos\alpha_s F_{B1}+i\sin\alpha_s F_{A1}$。因此参考效果在完整未知码上的精确压缩为

$$
E_s^\circ:=J^\dagger U_s^{\circ*}(N_B)J
=\tfrac12\cos^2\alpha_s\,I_2+\sin^2\alpha_s\,|1\rangle\langle1|.
\tag{12}
$$

此处$U_s^{\circ*}$仅表示理想CAR自同构对该有限effect的作用；无需授予全空间理想旋转为实际控制。对任意$\rho_{QR}$和$0\le E_R\le I$，参考联合概率即$\operatorname{tr}[(E_s^\circ\otimes E_R)\rho_{QR}]$。结合(5)、(9)，实际过程的同一个联合效果满足

$$
\left|p_{H,s}(N_B\otimes E_R)-\operatorname{tr}[(E_s^\circ\otimes E_R)\rho_{QR}]\right|
\le e_s:=80\lambda\epsilon/s+\delta_{\rm src}(T_*).
\tag{13}
$$

这是全Q/R的有界观测证书；没有把单个effect误差称作整个输出通道diamond误差。源比较(5)另有完整等距强界。

取旧合法逻辑rail1输入$J_1$，而不先测量B，精确参考概率为

$$
p_s^\circ=\sin^2(\lambda/s)+\tfrac12\cos^2(\lambda/s),\qquad
p_1^\circ=1,\quad p_2^\circ=3/4.
\tag{14}
$$

因B初始海占据为1/2，s=2的本读口概率是3/4，**不是**把B当空模所给的1/2。

用$\pi<22/7$取有理保守界，代码以Fraction独立复算：

|量|严格保守数值界|
|---|---:|
|$\delta_m=80\lambda\epsilon$|$<0.001257142857143$|
|$\delta_{\rm src}(T_*)$|$<0.000001011829029$|
|$e_1$|$<0.001258154686172$|
|$e_2$|$<0.000629583257601$|
|$p_{H,2}$|$\in[0.749370416742399,\ 0.750629583257601]$|
|$p_{H,1}-p_{H,2}$|$>0.248112262056227$|

数值展示向外取整；精确分数保存在results。它们来自同一个实际H的两次合法准备，不是两台外部脉冲器的比较。

## 5. 保1130旧钟和正窗口，不扩大其校准结论

完整288维旧D、$H_D$、原CO接收区、初始合法捕获后相位$\sigma_0=\pi/2$及1130的正窗口保持。令$K=H+H_D$，$\Delta_{\rm cal}=1/4$、

$$
W=\exp\{+i[\phi_\epsilon(0)/A]\otimes H_D\},\qquad
V_{\rm cal}=W e^{-i\Delta_{\rm cal}K}.
\tag{15}
$$

总过程仍$W e^{-iT_*K}$，所有D、源和材料输出保留。$[W,N_B]=0$严格成立，且$[H_D,N_B]=0$，故末端$N_B\otimes E_R$概率恰等于第4节的概率，**无须增加**旧$\varepsilon_{\rm cal}$。式(13)仍允许Q/D/R任意相关、源独立；实际许可准备仍是旧捕获后的D，并不因数学强界获得任意D态制备权。

同W在s=2时不再将旧钟校回$\sigma_0$。从(6)及谱有界性只能得到它接近$e^{+isT_*H_D}$，与等待合成后的旧内部相位为$\sigma_0-(s-1)T_*$；本轮不在这个新准备上复制1130的Mr=1/8结论，也不另开旧Mr未来仪器权限。

完整接触工作身份仍为

$$
W^\dagger K W-K=\Pi(\rho_\epsilon)H_D/A+q_\epsilon H_D^2/(2A^2).
\tag{16}
$$

初始源独立，且自由均值$\langle\Pi(\rho_\epsilon)\rangle=sA$，所以1130同一retarded界给$|\Delta E-s\langle H_D\rangle-q_\epsilon\langle H_D^2\rangle/(2A^2)|\le B\lambda LC_*q_\epsilon T_*/A^2$，$B=27$。工作随s改变，外部有效接触驱动的责任仍开放；本轮没有借读口对易删除它，也没有把正窗口变成自治微观触发。

## 6. 纯梯度候选的精确范围和端点字典

现在筛选一种候选，不将其宣布为已经存在的量子模型。取等质量两味Dirac字段，rail为旁观因子，$T_x=\sigma_x$，$\theta=\phi/A$，$G'=F$、$G(0)=0$。采用$+---$号差以及与正Hamiltonian脉冲一致的局部经典代表

$$
\mathcal L=\bar\psi(i\not\!\partial-m)\psi+\tfrac12(\partial\phi)^2
-\lambda\,\partial_\mu G(\theta)\,j_x^\mu,
\qquad j_x^\mu=\bar\psi\gamma^\mu T_x\psi.
\tag{17}
$$

若把耦合整体号换掉，下述端点旋转的号相反，本轮占据概率不变。令$g(\phi)=\lambda G(\phi/A)$，作**整份**非手征矢量味旋转

$$
\psi=e^{-ig(\phi)T_x}\chi.
\tag{18}
$$

动力学项产生$+\partial_\mu g\,j_x^\mu$，抵消(17)，质量项不变，因为$mI$与$T_x$对易。$j_x$在该旋转下不变，故经典方程上$\partial_\mu j_x^\mu=0$。相互作用可在壳写为全散度；(18)则直接给精确作用变换，不是只用低阶方程删除一个顶点。这一论证不用于不等质量、轴流、手征反常或非对易的多个味连接。

指定平滑经典源时，裸端口传播为$U_g=R_f U_0R_i^\dagger$，$R=e^{-igT_x}$。在平台$\theta=st$上，$G(sT_*)=1$、$G(0)=0$，因此两个s的端点角都是$\lambda$。**共同变换输入和读口**时可得到自由描述；**固定原A/B裸端口**时仍有真实可读的端点相对旋转。不能因可场重定义便宣称A到B转移不存在。

源字典也须一起变换。(17)的canonical源动量和形式Hamiltonian密度分别是

$$
\pi_\phi=\dot\phi-g'(\phi)j_x^0,\qquad
\mathcal H=\psi^\dagger h\psi+\tfrac12[\pi_\phi+g'j_x^0]^2
+\tfrac12|\nabla\phi|^2+g'\nabla\phi\cdot\mathbf j_x.
\tag{19}
$$

源Euler变分中$g''\partial\phi\cdot j$项相消，仅剩与$\partial\cdot j$成正比的项；完整局部等质量经典系统的物理$\phi$因而满足自由波方程。它不等于在原1129 H里宣布源反作用为零。保同canonical $\phi,\pi$而忘记(19)，一般不保同物理Cauchy数据或量子准备；平方中的交叉／二次流项亦不可删除。式(19)是经典或有完整正则化字典时的形式身份，逐点流平方未在1129正常Fock中定义。

成熟文献明确区分守恒矢量流和轴流：Di Luzio等的[原论文Appendix A、式(51)附近](https://arxiv.org/html/2305.11958v2)删除守恒对角矢量导数项；本例等质量$T_x$可由常味旋转对角化，具体(18)仍由本节直接核。Criado与Pérez-Victoria的[§2–3](https://arxiv.org/html/1811.09413)说明场变换应连同来源项、正则化及高阶项处理，单靠低阶EOM删除不能代替完整变换。这些文献不是本轮正常梯度QFT的存在定理。

## 7. 有限候选匹配定理：不冒充候选的正常实现

声明候选类$\mathfrak C_{\rm grad}(\kappa_{\rm grad})$：候选若主张在M1131两个准备、同原码、同裸B效果及固定终点上实现最小纯梯度补全，须另给合法来源／场字典、完整实际输出，并证明至少在$s=2,J_1,N_B$上

$$
|p_{{\cal C},2}-1|\le\kappa_{\rm grad}.
\tag{20}
$$

这里1是第6节的裸端点参考预测。候选可用直接有限effect界证明(20)；若另有保持源、海、相位和参考的全输出半迹距或等距界，也可推出它。**本轮不提供、设为零或从经典作用推出$\kappa_{\rm grad}$**。无此受控端点性质的方案不被定理排除，但也不能无说明继承纯梯度平台预测。

由实际原H的(13)–(14)，对于每个满足(20)的候选，

$$
\boxed{|p_{{\cal C},2}-p_{H,2}|\ge
\tfrac14-e_2-\kappa_{\rm grad}
>0.249370416742399-\kappa_{\rm grad}.}
\tag{21}
$$

下界为负时理解为平凡界。若候选同时要求在此effect与原H匹配至$\varepsilon_{\rm match}$，则必要条件为

$$
\varepsilon_{\rm match}+\kappa_{\rm grad}\ge\tfrac14-e_2.
\tag{22}
$$

例如两个预算各为0.01时不可能通过。它只需一个旧合法逻辑输入，但原H上(5)、(13)已保全未知Q/R量词。匹配全任务必匹配这个有界边缘；完整源或余部不可隐藏这个差异。反过来，通过本effect不足以证明全来源或全任务匹配。

这个结果不否定s=1：其原脉冲角与纯梯度端点角都为$\lambda$。它也不要求梯度候选复制原固定n模型所有可能源准备；只对**额外接受M1131这个有限菜单且声称保原操作合同**的完成施加约束。把A改成sA、重新归一G、按s换读口或重设终点都改变了此合同。

## 8. 部署边界、平滑及归一梯度不属同一排除类

1129含固定$\eta$及非局部平滑$K_\epsilon$。在局部形式层面，$a=\lambda\eta F(\theta)d\theta$一般有

$$
da=\lambda F(\theta)d\eta\wedge d\theta;
\qquad d[\lambda\eta G(\theta)]=\lambda\eta F(\theta)d\theta+\lambda G(\theta)d\eta.
\tag{23}
$$

所以带部署边界的第一式并非全局纯梯度；补成第二式必须保留额外空间边缘项。平滑CAR流也不是未经证明就有精确局部守恒的逐点$j_x$，场变换与K一般不交换。不能只把(1)中的n换成梯度而照抄(18)，更不能据此宣布共同相互作用四P闭合。

本轮筛选定理只要求候选对原码支撑管的端点性质给出(20)。边界、正则化、源码海后态和量子涨落的影响若能控制，可并入$\kappa_{\rm grad}$；若产生大效应就须明报，不能引用纯梯度工具代替估算。不要求先解决全部UV问题才能陈述这个有限必要条件，也不把没有该实现当作阶段失败。

另一个候选是归一梯度。用$X=\partial_\mu\theta\partial^\mu\theta>0$写

$$
\widetilde a=\lambda F(\theta)\frac{d\theta}{\sqrt X},\qquad
d\widetilde a=-\frac{\lambda F(\theta)}{2X^{3/2}}dX\wedge d\theta.
\tag{24}
$$

一般curl非零，不能套纯梯度移除。在均匀正速率平台$\theta=st$上，$d\theta/\sqrt X=dt$，其比较角正是$\lambda\int F(st)dt=\lambda/s$，**两份准备都与原角相同，逃过本轮见证**。这并非证明它完整正常、稳定、局域或与原来源后态匹配；X接近零、量子定义、实际反作用及总四P仍须另证。但本轮不为排除它另造第三份例子。

## 9. 去重、数值责任及停止线

[564](../../archive_554_584/research_note_564.md)已说明固定物质端点后的连接改变不能当作纯规范；[606](../../archive_585_628/research_note_606.md)已要求场依赖过程映射保留导数及完整来源；[1028](../../archive_1009_1043/research_note_1028.md)已区分转动Higgs相位和连同连接的纯规范配置。1015–1016的场重定义／边界商、1017混合坐标不等于不可分部门、1018完整反作用、1019实际有限effect、1027真实交换与1029高阶变换范围均直接继承。

新结果是(13)–(14)这个实际正常源有限菜单证书及(21)的候选匹配下界。代码四组检查：冻结旧来源和轮次身份；有理预算／平台／资源；保海的完整四模CAR、B余部纯化和逻辑R；实际F的两类脉冲面积以及局部curl／边界差。后两组是独立机制诊断，不把有限CAR或经典场导数冒充连续源证明。

共同四P、动态定位支撑、梯度量子存在性、旧Q_r接管及接触自治生成均未在本轮签收。只对M1131和$\mathfrak C_{\rm grad}(\kappa_{\rm grad})$作有选择力的必要检验；不要求修好这个单一候选才允许阶段推进。
