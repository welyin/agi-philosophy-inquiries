# 第1293轮证明：原 Dirac 质量、分离的左手记录与共同二次来源

2026-10-18。依据[准入](../_shared/round1293_admission.md)。本文证明一个明确二次有效对象的联合证书；完整相互作用 SM$_n$、实际末读装置及其误差未获认证。参数事前固定，不扫描。

## 1. 同字段字典、作用与量词

取标准荷 $k=0$、含原右手中微子共轭 $N$ 的 Dirac 分支，$n=3$ 或 $4$。原字段是 $L_a,e^c_a,N_a,H$ 及同一模型的其余字段。取 $\hbar=c=1$、$g_{\mu\nu}=\operatorname{diag}(1,-1,-1,-1)$。原 Higgs 代表及参数为

$$
H_0=(0,v/\sqrt2),\quad v=\sqrt2,\quad A_\mu=0,\qquad
Y_\nu|_{1,2}=I_2,\quad Y_e|_{1,2}=2I_2 .
\tag{1}
$$

两个指定代际不与其余代际混合；其余块为公开固定输入。带电轻子块也同值是为使完整相关二次作用在这两代的共同 $U(2)$ 下不变，不是推导现实简并谱。规范场、Higgs 及度规仍作为同一作用的背景变量保留。

在 $k=0$，$\widetilde H^\dagger L_a$ 是原全规范群的中性局部复合插值量，$N_a^\dagger$ 是其右手伙伴。字典为

$$
\chi_{L,a}=\frac{\sqrt2}{v}\widetilde H^\dagger L_a
=\nu_{L,a}+\text{含 Higgs 涨落的高阶项},\qquad
\chi_{R,a}=N_a^\dagger .
\tag{2}
$$

**只有完整复合量具有原理论的规范不变身份。下述精确 CAR、自由传播和有界效果属于声明的领先二次模型，不属于完整复合量的全阶算符身份。** 这一区别见[2011.02301，§II 式(3)—(4)](https://arxiv.org/html/2011.02301)；该文的 vectorial 格点例不是本 SM$_n$ 任务的实验或误差证书。这里不借它证明有限读口。

先保留原背景费米作用
$S_F^{(2)}[g,A,H]$ 的协变动能与
$-\overline L Y_\nu\widetilde H\nu_R-\overline L Y_eH e_R+\mathrm{h.c.}$，
**先变分再取(1)**。在(1)上，原两中微子字段的精确二次 Hamiltonian 为

$$
h_D=
\begin{pmatrix}
-\boldsymbol\sigma\cdot p&mI_2\\
mI_2&\boldsymbol\sigma\cdot p
\end{pmatrix},\quad m=1,\quad
h_D^2=(-\Delta+m^2)I_4,\qquad U_t=e^{-ith_D}.
\tag{3}
$$

旧离散 $r$ 的原字段字典为 $L\mapsto iL$、$H\mapsto-H$、$N\mapsto iN$，故(2)的 $\chi_L$ 与 $\chi_R$ 均取 $-i$；单占据两味码具有同一该离散荷及同一宇称。这里不重新分类全局群。

同一作用还含原带电及其余费米块；共同参考和各自来源保留，不能只取中微子结果然后把其他部门的来源设零。本文给固定背景上的精确二次场过程及可耦合来源，不宣称该平直背景已解这些非零来源的自引力方程。$n=3,4$ 各有本实例；不要求两模型的全部来源或预测相同，不选择三代。

旧工具直接继承：[730 §3—6](../../archive_702_741/research_note_730.md)的原场 Cauchy/Hadamard 与来源；[798 §3—5](../../archive_797_818/research_note_798.md)的原 sterile 偶记录；[813 §3—7](../../archive_797_818/research_note_813.md)的同模式局部化及来源支持；[818 §3、8](../../archive_797_818/research_note_818.md)的共同未知输入及来源；[1287 §3—6](../../archive_1281_1287/1287/proof.md)只复用 bump/严格时样核估计方法，不搬其 KG 信号或 Dirac tag。

## 2. 紧支模式与不改变初始远端边缘的准备

在全 Cauchy Dirac CAR 空间中工作，**不作正频投影后冒称紧支**。定义

$$
a=\frac1{32},\quad L=\frac14,\quad I=[7/16,1/2],\qquad
b_A(x)=
\begin{cases}
\exp[-1/(1-|x|^2/a^2)],&|x|<a,\\
0,&|x|\ge a ,
\end{cases}
\quad b_B(x)=b_A(x-Le_1),\quad g_{A,B}=\frac{b_{A,B}}{\|b_{A,B}\|_2}.
\tag{4}
$$

取相同单位二旋量 $\zeta$；$f_j=(0,g_A\zeta)\otimes e_j$ 为两味右手模式。记等距列 $V=(f_1,f_2)$、$F=VV^\dagger$、$Q=1-F$。令 $P$ 为完整二次模型的真空 Cauchy 协方差，中微子约定 $P=1_{(-\infty,0)}(h_D)$、$\omega(a^\dagger(k)a(h))=\langle h,Pk\rangle$。共同参考保前两代 $U(2)$。取它在 $\operatorname{CAR}(Q\mathcal K)$ 上的完整限制 $\sigma_Q$。

有限模式 CAR 有分级分解。以 $W|0\rangle=a^\dagger(f_1)|{\rm vac}\rangle_F$、$W|1\rangle=a^\dagger(f_2)|{\rm vac}\rangle_F$ 的单占据子空间编码，定义

$$
\mathcal E(\rho_{QR})
=(W\otimes I_R)\rho_{QR}(W^\dagger\otimes I_R)
\ \widehat\otimes\ \sigma_Q .
\tag{5}
$$

此处下标 $Q$ 在 $\rho_{QR}$ 表示逻辑输入，与投影 $Q=1-F$ 由语境区分；$R$ 是不受作用的任意伙伴。所有码态在同一奇宇称扇区，其密度矩阵是偶的，故(5)与共同偶参考的分级乘积有定义。它是有限模式的 CPTP 替换/编码合同，保全部未知逻辑资料；不假定原真空在 $F$ 与其余场之间因子化。原有限模式被弃去的关联不能声称仍在，而其余代数的**整个限制**保持。

更确切地，对任意初始接收偶代数
$\mathfrak B_0\subset\operatorname{CAR}(Q\mathcal K)$ 及任意 $R$ 效果 $B_R$，

$$
\omega_{\mathcal E(\rho)}(B\,B_R)
=\omega_P(B)\operatorname{Tr}(\rho_R B_R).
\tag{6}
$$

包含 $B$ 球的任意初始远端区域都满足此前件；整个初左手代数也与 $F$ 正交。故初接收端既不携带输入，也不与其伙伴相关。该事实由保持完整限制的替换映射证明，不由两点相同猜测全态相同。对支撑在 $A$ 的偶局部操作，因果外偶代数不变；(5)不承诺该替换已经由资源闭合的自治装置实现。尤其不能把裸奇 Kraus 写成无需额外资源的偶自治准备。固定单占据及共同参考本身是本轮登记的初态资源。

空码参考协方差为 $C_0=QPQ$。编码态通常**非准自由**；只对二次量可使用

$$
C_\rho=C_0+V\rho_QV^\dagger,\qquad
\Delta C_\rho=-FP-PF+FPF+V\rho_QV^\dagger .
\tag{7}
$$

完整未知伙伴及高阶记录按(5)处理，不能把 $C_\rho$ 代入任意阶 Wick 公式。例如两味混合单占据态满足 $\langle n_1n_2\rangle=0$，而其准自由替代会给 $\det\rho_Q$。

有限 CAR 替换在原表示中给正常态；(7)含有限个光滑 Cauchy 模及其 $P$ 像，传播后为平滑双解差，故保 Hadamard 短距及固定涂抹 Wick 域。旧730/818提供这一工具；这里不再证明一般 Hadamard 存在。

## 3. 原质量给严格分离的左手信号

定义左手接收模式 $h_j=(g_B\zeta,0)\otimes e_j$，局部偶投影
$E_j(t)=a_t^\dagger(h_j)a_t(h_j)$。它是所选二次模型的有界局部效果，不是完整 SM 的已实现末读装置。两输入 $|0\rangle,|1\rangle$ 的响应差为

$$
p_0(t)-p_1(t)
=\left|\langle h_1,U_tf_1\rangle\right|^2=|\alpha(t)|^2.
\tag{8}
$$

共同参考贡献相消；原演化和准备均使用同一场。更一般地，接收偶双线性 $a_t^\dagger(h_i)a_t(h_j)$ 的输入依赖为 $|\alpha(t)|^2(\rho_Q)_{ji}$，因此输入相干也进入同一接收代数；这不等于整 qubit 无噪恢复或已实现末读装置。由(3)，

$$
P_LU_tP_R=-im\,\frac{\sin(t\sqrt{-\Delta+m^2})}{\sqrt{-\Delta+m^2}}I_2.
\tag{9}
$$

三维 KG 初速度核在 $t>r$ 的尾为 $-mJ_1(ms)/(4\pi s)$，
$s=\sqrt{t^2-r^2}$；锥面项为 $\delta(t-r)/(4\pi r)$。所以(9)的严格时样尾为
$im^2J_1(ms)/(4\pi s)$。本窗全部源—收点满足

$$
3/16\le r\le5/16<t,\qquad
3/32\le s^2\le55/256<1/4 .
\tag{10}
$$

锥面项无贡献。利用 $J_1(s)/s=\sum_{j\ge0}(-1)^js^{2j}/[2^{2j+1}j!(j+1)!]$ 的交错递减项，

$$
\frac{J_1(s)}s\ge\frac12-\frac{s^2}{16}
>\frac{31}{64},\qquad
\left|\frac{iJ_1(s)}{4\pi s}\right|>\frac{31}{256\pi}.
\tag{11}
$$

级数与余项可复核[NIST DLMF 10.2.2](https://dlmf.nist.gov/10.2.E2)。核在整个积分域具有同一相位，正 bump 不发生抵消。

令 $Z=\int b_A$。内半球 $b_A>1/4$，且 $\sup b_A<1/2$，因此

$$
Z>\frac{\pi a^3}{24},\qquad
\int b_A^2<\frac12 Z,\qquad
\left(\int g_A\right)\left(\int g_B\right)>2Z>
\frac{\pi a^3}{12}.
\tag{12}
$$

于是全窗严格得到

$$
|\alpha(t)|>\frac{31a^3}{3072}
=\frac{31}{3\cdot2^{25}}>\frac{31}{2^{27}},
\qquad
p_0(t)-p_1(t)>\frac{961}{2^{54}}.
\tag{13}
$$

本结论无需数值传播、积分扫描或渐近小时间。任意事前非负、归一且支撑在 $I^\circ$ 的光滑时间权平均 $E_j^\chi=\int\chi(t)E_j(t)dt$ 仍是该有限世界管的偶有界效果并保(13)。这不自动授予实现该平均的经典仪器。$t<3/16$ 时局部因果传播给零输入差。

## 4. 能源及全部背景来源属于同一二次作用

在固定物理耦合和同一局部协变、保 Ward 的重整化处方下，定义

$$
J_\zeta=\delta_\zeta S_F^{(2)},\qquad
\delta S_F^{(2)}
=\int\sqrt{-g}\left[-\frac12 T^{\mu\nu}\delta g_{\mu\nu}
+j^\mu_a\delta A^a_\mu+\mathcal F_I\delta H^I\right].
\tag{14}
$$

$H^I$ 是四个实 Higgs 分量；规范及 spin 连接变分全部保留。真空基线和固定局部反项仍在，不因态差光滑就设为零。有限来源陈述是对每份固定 $C_c^\infty$ 测试的期望，未给全部应力涨落或无限菜单统一界。

前两代共同 $U(2)$ 在物理 $L,\nu_R=N^\dagger,e_R=(e^c)^\dagger$ 上取 $U$，在全左手记号 $N,e^c$ 上取 $U^*$；(1)的二次作用及参考不变；每个(14)的背景方向在代际中为单位。因此(5)的码源泛函为

$$
\mathcal E_t^*(J_\zeta)=j_\zeta(t)I_2,
\qquad
\omega_{\mathcal E_t(\rho_{QR})}(J_\zeta B_R)
=j_\zeta(t)\operatorname{Tr}(\rho_RB_R).
\tag{15}
$$

二次形式域内理解此式；不把无界来源称作有界通道效果。直接证明亦可：输入依赖仅为
$U_tV\rho_QV^\dagger U_t^\dagger$，所有中性—中性来源块在味中为单位；$W^\pm$/带电 Higgs 变分的中性—带电混合腿，在共同无交叉参考中压缩为零。剩余带电及其它物种来源为输入无关基线。先删去这些顶点再只核 $T$，不是此证明。中性 Higgs 相位（Goldstone）及径向来源仍保留；码上标量不表示来源为零。

任意独立 Yukawa 矩阵 spurion 的变分不在菜单内；它一般不与共同 $U(2)$ 对易，不能领取(15)。条件化末读或加入味敏感装置后，也不能继续无条件套(15)。

守恒使用背景 Ward，而非宣布费米部门的每一份荷独立守恒。对紧支规范变化 $\delta_\epsilon A=-D\epsilon$、
$\delta_\epsilon H=\epsilon^a t_aH$ 及紧支微分同胚 $X$，在壳有

$$
J_{\delta_\epsilon(g,A,H)}=0,\qquad
J_{\mathcal L_X(g,A,H)}=0 .
\tag{16}
$$

按(14)的来源符号，第一式为
$D_\mu j^\mu_a+\mathcal F_I(t_aH)^I=0$；
第二式保所有 $j\,\mathcal L_XA$ 和 $\mathcal F\,\mathcal L_XH$ 项。
在本固定平直、$A=0$、常量 $H$ 背景，第二式才化成
$\partial_\mu\langle T^{\mu\nu}\rangle=0$。
不可由此删去一般背景的 Higgs/规范交换项。全原字段的反常一致性按已验收[1060](../../archive_1046_1062/1060/proof.md)继承，未用两个中性领先分量替代全反常账。

对平直真空正常序能源，取两选定中微子块
$P=P_-$；径向实 $g_A$ 的纯右手模满足
$V^\dagger h_DV=0$。因 $h_DP_-=(h_D-|h_D|)/2$，

$$
E_{\rm rel}=\operatorname{Tr}(h_D\Delta C_\rho)
=\operatorname{Tr}(F|h_D|)
=2\langle g_A,\sqrt{-\Delta+1}\,g_A\rangle ,
\tag{17}
$$

与未知码及 $R$ 无关、严格正、时间守恒。这里不是正能紧支单粒子；正常序正能源由参考的粒子/空穴及码的共同真实态承担。

$|\nabla b_A|<2/a$、$\int b_A^2>\pi a^3/96$ 给
$\|\nabla g_A\|_2^2<512/a^2$，从而

$$
0<E_{\rm rel}
\le2\sqrt{1+\|\nabla g_A\|_2^2}
<2\sqrt{524289}<1450 .
\tag{18}
$$

所有固定能量矩也因光滑紧支初始模而有限；相应平滑来源测试可直接收缩(7)。这个正常序**相对**能源界不是准备器成本、局部正能断言、引力弱场界或绝对真空来源归零。

## 5. 删除、误差和本轮净增量

删除 Dirac Yukawa 即 $m=0$，并使用删除模型自身的正常参考及同样初始盲性合同，(9)严格为零：指定右→左任务没有输入信号。不是在原真空中免费突然 quench 质量；其它质量机制或其它任务未被排除。

在此二次模型内(5)—(18)精确成立。若以后给出实际两输入输出相对这些理想限制态的误差 $\varepsilon_0,\varepsilon_1$（以所有允许效果的差为界），则同一接收效果的差至少为
$961/2^{54}-\varepsilon_0-\varepsilon_1$。
**本轮没有提供完整 SM 相互作用、Higgs 复合插值高阶或实际读口的这些预算，不能将它们赋零。** 源的匹配误差还需相应来源域预算，不能只由有界效果误差推得。

新证书是：原物种中同一非零 Dirac 质量，在一份正常、有限资源的未知味输入合同上，给全分离时间窗的左手局部记录及同一完整二次背景来源。旧一般 CAR/来源/Schur 工具只作复用。数值检查是一组有限代数与严格常数校准，不替代连续证明，不证明整个 SM$_3$/SM$_4$ 的 O 资格，不排除其它代数或代际。

[复算代码](check.py)默认及 `--check` 只读，`--write` 仅首次排他创建[结果](results.json)。主线验收前本文不自授最终签收。

