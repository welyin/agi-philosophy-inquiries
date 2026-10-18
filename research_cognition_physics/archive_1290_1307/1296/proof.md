# 第1296轮证明：同一有限未知记录合同跨越Dirac与Majorana分叉

2026-10-18。依据[事前准入](../_shared/round1296_admission.md)，仅证明原字段的明确背景二次有效合同。一般自对偶CAR、Hadamard与传播工具复用[730 §3—5](../../archive_702_741/research_note_730.md)、[813 §3—7](../../archive_797_818/research_note_813.md)、[818 §3、8](../../archive_797_818/research_note_818.md)；连续尾核与非退化味资料复用[1295](../1295/proof.md)。本轮不采用新公理，不将实际有限记录等同于全相互作用SM、无噪恢复或已制造的读出装置。

## 1. 同一作用族及明确的保护量

取1295的 $k=0,n=3$、签名 $(+---)$、自然单位，原 $L,e^c,N,H$ 及其余物种全部保留。固定

$$
\alpha_*={31\over2^{27}},\quad \epsilon={\alpha_*\over128},\quad
Y_\nu=I_3+\epsilon(D_3+P_3),\quad Y_e=2I_3+\epsilon D_3,
\quad H_0=(0,1)^T,\quad A_\mu=0,
\tag{1}
$$

其中 $D_3=\operatorname{diag}(1,2,3)$，$P_3$ 为路径邻接矩阵。1295已证两谱非零简单且共同左荷对易核仅为标量；这里不重计该分类。加入[1292 §1](../1292/proof.md)原有的一份复 $S$ 及规范不变项 $-\tfrac12(SN^Ty_NN+\mathrm{h.c.})$。取同一实 $S_0=w/\sqrt2=1$、$w=\sqrt2$，正则标量动能及固定非负径向势，两个实 $S$ 方向均保留。令 $y_N=\delta I_3$，其中

$$
0\le\delta\le\epsilon .
\tag{2}
$$

按所选Weyl/共轭字典，背景上的物理配对矩阵固定为 $\Delta_{RR}=\delta i\sigma_2\otimes I_3$；相位约定不改变以下范数与源身份。不同 $\delta$ 是同父作用族的不同常耦合，**不是**在一份真空上突然改质量的quench。同一非零 $S_0$ 已选择原 $r:S\mapsto-S$ 的破缺背景；$\delta=0$ 只称费米Dirac点，不称未破缺 $r$ 相。原复 $S$ 相位模及来源不删除。

保护量取原含 $N$ 的费米 $B-L$。在物理中微子粒子/空穴列中，它的表示生成元可写为 $\mathcal B=\operatorname{diag}(-I,+I)$；整体符号不影响是否对易。它是**固定背景二次费米**的生成元。$\delta=0$ 时与演化对易，$\delta>0$ 时不对易。母作用可令 $S$ 承担 $B-L=-2$，保含 $S$ 的总荷；不能把费米数交换叫作总荷、规范Ward或能量守恒失败。群与反常仍按[1060](../../archive_1046_1062/1060/proof.md)的原背景范围继承，不重算或扩大。

## 2. 自对偶演化与每模型自身参考

以原复Cauchy空间 $\mathcal K$ 表示Dirac点，把 $\mathscr K=\mathcal K\oplus\overline{\mathcal K}$ 作为Nambu空间。反线性等距 $\Gamma(u,\bar v)=(v,\bar u)$ 满足 $\Gamma^2=1$。适当固定列约定下，

$$
\mathbb H_\delta=\begin{pmatrix}h_Y&\Delta_\delta\\
-\overline{\Delta_\delta}&-\overline{h_Y}\end{pmatrix},\qquad
\Delta_\delta^T=-\Delta_\delta,\quad
\mathbb H_\delta^*=\mathbb H_\delta,\quad
\Gamma\mathbb H_\delta\Gamma=-\mathbb H_\delta .
\tag{3}
$$

共轭微分项含动量反号；不能在一个固定 $p$ 符号上漏掉它。$\Delta_\delta$ 仅在原右手中微子块，$\|\Delta_\delta\|=\delta$，所以 $\|\mathbb H_\delta-\mathbb H_0\|=\delta$。$\mathbb U_\delta(t)=e^{-it\mathbb H_\delta}$ 与 $\Gamma$ 对易，保自对偶CAR与宇称，不保普通粒子数。直接有 $\|[\mathbb H_\delta,\mathcal B]\|=2\delta$；这是指定保护量两支的差别。

每支取自身常系数真空 $P_\delta=\mathbf1_{(-\infty,0)}(\mathbb H_\delta)$，满足 $\Gamma P_\delta\Gamma=1-P_\delta$。完整中性质量矩阵的Takagi分解给有限多个自由Majorana质量；由 $Y_\nu$ 的正下界与小 $\delta$，质量不穿零。因此该参考是本方程的Hadamard态；同方程正常有限CAR操作仍在其表示中。这里不要求不同质量真空在一个预定全局Fock表示中等价，也不复用别支的 $P$。730的一般背景构造仅作成熟适用性依据，非本轮重证。

## 3. 完整未知码、环境反常关联与初远端盲性

取1295原 $a=1/32,L=1/4,I=[7/16,1/2]$，同实径向归一紧支 $g_A,g_B$、同单位二旋量 $\zeta$。真实两模由右手 $f_i=g_A\zeta\otimes e_i$ 定义，左手接收模式为 $h_i=g_B\zeta\otimes e_i$，$i=1,2$。令 $V=(f_1,f_2)$ 嵌入Nambu粒子列，$F=VV^\dagger$，

$$
\widehat F=F+\Gamma F\Gamma,\qquad Q_\perp=1-\widehat F .
\tag{4}
$$

$\widehat F$ 秩4且 $\Gamma$ 不变，对应**两个**真实CAR模，不是四个物种。其代数为两模CAR；以 $W|0\rangle=|10\rangle,W|1\rangle=|01\rangle$ 编码任意 $\rho_{QR}$。补代数取自身参考的整个限制 $\sigma_{\delta,\perp}$，包括全部正常与反常关联。准备为

$$
\mathcal E_\delta(\rho_{QR})=(W\otimes I_R)\rho_{QR}(W^\dagger\otimes I_R)
\ \widehat\otimes\ \sigma_{\delta,\perp} .
\tag{5}
$$

固定单占据码的密度矩阵为偶，分级乘积与任意逻辑伙伴扩张合法。这是一份公开正常有限模式CP初态合同，不假定原真空因子化，不声称准备器已经自治制造。初远端完整偶接收代数在补空间内，故对其任意 $B$ 和伙伴效果 $B_R$，准确有 $\omega(BB_R)=\omega_\delta(B)\operatorname{Tr}(\rho_RB_R)$。有限传播给因果接触前无输入差；完整等待演化及全部后态保留伙伴，不只保某一边缘。

为二次量计算，码协方差写为

$$
P_{\rm code}(\rho)=V\rho V^\dagger+\Gamma V(1-\rho)V^\dagger\Gamma,
\quad
S_{\delta,\rho}=Q_\perp P_\delta Q_\perp+P_{\rm code}(\rho).
\tag{6}
$$

式中 $\Gamma$ 反线性，第二块等价于共轭基中的 $1-\rho^T$。完整混合单占据态一般**不是准自由态**；例如双占据为零，而错误Wick替代给 $\det\rho$。未知伙伴、全部高阶统计必须用(5)，不能用(6)虚构Gaussian态。$S_{\delta,\rho}-P_\delta$ 为光滑有限秩腿的组合；同方程短距和固定紧支来源域沿730/813/818保持。

## 4. 全窗正常/反常两腿及有限下界

把接收湮灭算符在初片分解为码上的正常腿 $A_\delta a_F$、反常腿 $B_\delta a_F^\dagger$ 和补空间项。码内异常二点为零；补空间的全部正常/异常贡献及 $B_\delta B_\delta^\dagger$ 是共同基线。因此对Hermitian输入差 $X$，接收二点依赖为

$$
\mathcal L_\delta(X)=A_\delta XA_\delta^\dagger-B_\delta X^TB_\delta^\dagger .
\tag{7}
$$

这是输入差的实线性映射，**不是**把负号项单独当成CP通道；完整CAR演化与共同环境才给合法输出。

相对于同一实际 $Y$、$\delta=0$ 的块对角演化，Duhamel给 $\|\mathbb U_\delta-\mathbb U_0\|\le t\delta$，故 $\|B_\delta\|\le t\delta$。再加1295相对 $Y_\nu=I$ 基点的 $3\epsilon$ 界和原严格尾核 $|\alpha_0(t)|>\alpha_*$，全 $I$ 有

$$
\|A_\delta-\alpha_0 I_2\|\le e:={7\alpha_*\over256},\qquad
\|B_\delta\|\le b:={\alpha_*\over256},\qquad
s_{\min}(A_\delta)>{249\alpha_*\over256}.
\tag{8}
$$

同一局部偶效果 $E_1(t)=a_t^\dagger(h_1)a_t(h_1)$ 比较两基输入，由(7)

$$
\begin{aligned}
\Delta p_\delta(t)&=|A_{11}|^2-|A_{12}|^2-|B_{11}|^2+|B_{12}|^2\\
&\ge(|\alpha_0|-e)^2-e^2-b^2
> {61951\over65536}\alpha_*^2>{15\over16}\alpha_*^2 .
\end{aligned}
\tag{9}
$$

$e<\alpha_*$ 保右式对 $|\alpha_0|$ 单调。对Hilbert—Schmidt范数及非零Hermitian $X$，

$$
\|\mathcal L_\delta(X)\|_2
\ge\big(s_{\min}(A_\delta)^2-\|B_\delta\|^2\big)\|X\|_2
>{3875\over4096}\alpha_*^2\|X\|_2 .
\tag{10}
$$

$X=0$ 时只能用非严格等式/不等式。因而接收二点资料保未知二维信息的单射，且延拓至伙伴；不意味着无噪恢复或层析仪器。事前非负归一时间权的概率平均保(9)，不需无限精确的单时刻。小矩阵只核代数和预算，不能替代1293原连续尾核及本解析全窗界。

## 5. g/A/H/S全部来源及真正守恒账

固定每模型所有物理耦合，先变分原同一局部作用，再代经典背景。来源包含度规/标架、全部规范方向、四实Higgs和两实 $S$ 方向：

$$
\delta S^{(2)}_F=\int\sqrt{-g}\left[-\tfrac12T^{\mu\nu}\delta g_{\mu\nu}
+j_a^\mu\delta A_\mu^a+\mathcal F_{H,I}\delta H^I
+\mathcal F_{S,J}\delta S^J\right].
\tag{11}
$$

同一局部协变重整化规则、旁观字段、参考绝对源、标量/规范背景自身与固定反项均保留。不能先在背景上删掉相位方向，不能逐模型重调真空常数隐藏变化。任意Yukawa spurion变分不属于该来源菜单。

记紧支测试 $\zeta$ 的运输二次核为 $\mathfrak j_{\delta,\zeta}$，满足 $\Gamma\mathfrak j\Gamma=-\mathfrak j$。物理期望用Nambu半权 $\tfrac12\operatorname{Tr}(\mathfrak j S)$ 加同一参考常数。令 $F_c=\Gamma F\Gamma$，由(6)准确得到

$$
\mathcal E_\delta^*(J_\zeta)=c_{\delta,\zeta}I_2+V^\dagger\mathfrak j_{\delta,\zeta}V,
\quad
c_{\delta,\zeta}=\omega_\delta(J_\zeta)+\tfrac12\operatorname{Tr}
\mathfrak j_{\delta,\zeta}(-\widehat FP_\delta-P_\delta\widehat F
+\widehat FP_\delta\widehat F+F_c).
\tag{12}
$$

迹只含光滑有限秩腿，不是无限海的裸迹。对伙伴效果仍按 $c\operatorname{Tr}\rho_R B_R+\operatorname{Tr}[\rho_{QR}(V^\dagger\mathfrak j V\otimes B_R)]$ 运输。一般来源为逻辑矩阵并与R关联，不要求标量。正常/异常源插入均由同一式计算。

规范Ward含Higgs和S的背景变换；在 $k=0$ 时S规范中性，但其全源仍存在。对费米 $B-L$，右手配对块有非零对易项，恰由 $\delta S=-2i\vartheta S$ 的源项配平：质量符号层满足 $i[\mathbb H,\mathcal B]+\delta_{B-L}\mathbb H=0$。这是含S母作用的联合Ward；不是声称固定S背景下费米B−L单独守恒。微分同胚Ward同样保背景交换；仅在所取常背景中完整二次费米应力才有普通平移守恒。

每份固定紧支源测试的有限性来自有限光滑腿和实际导数域，复用[805](../../archive_797_818/research_note_805.md)。有界传播差不能控制所有无界源、高阶噪声或跨模型应力误差。本轮给来源身份与有限性，不给尚未证明的全相互作用余项或自引力解。

## 6. Nambu半权下的正有限能源

原实径向包的动量平均为零；Dirac质量将右手送往正交左手；同一 $\zeta$ 的配对压缩满足 $\zeta^T i\sigma_2\zeta=0$。对包括共轭的四列同时检查即得

$$
\widehat F\mathbb H_\delta\widehat F=0.
\tag{13}
$$

这一步不能只查普通F或漏空穴。利用 $\mathbb H P_-=(\mathbb H-|\mathbb H|)/2$，(6)给任意未知输入的相对总能

$$
E_{\rm rel}=\tfrac12\operatorname{Tr}\mathbb H_\delta(S_{\delta,\rho}-P_\delta)
=-\operatorname{Tr}(\widehat F\mathbb H_\delta P_\delta)
=\tfrac12\operatorname{Tr}(\widehat F|\mathbb H_\delta|)>0.
\tag{14}
$$

配对质量是Lorentz标量零阶项；完整中性质量范数至多 $\|Y_\nu\|+\delta<2$。主部反对易或Takagi分解给质量平方界，且原 $\|\nabla g_A\|^2<524288$。秩4连同半权于是给

$$
E_{\rm rel}<2\sqrt{524288+4}<1450 .
\tag{15}
$$

等待演化保此总相对能，每个固定能量矩有限。此为准备态相对自身参考的能源，不是制造该准备的成本、局部正能或所有来源恒定。码上能源与输入无关不把(12)其余来源变成标量；绝对基线仍按(11)保留。

## 7. 新增资格与停止线

同一具名 $n=3$ 非退化原字段、同两球、同全窗、同未知码/伙伴及同 $g,A,H,S$ 来源字典，在Dirac点和非零Majorana区间都满足(9)—(15)，而指定二次费米 $B-L$ 的对易保护不同。**这份有限记录合同不强迫精确Dirac保护。** 它不否定更强独立保护任务，不选择现实质量机制，也不是全部认知公理的反模型。

730已覆盖一般Majorana参考和记录，813/818已覆盖自对偶局部化与完整未知准备；它们未给本两球/全窗和1295非退化父的上述共同有限下界。本轮净增只在这份跨质量分叉的联合身份。不同固定耦合与各自参考不被伪装为同初态同Hamiltonian的两个后继。常背景二次模型、初态CP合同和有限数学效果的权限不扩大成动态量子Higgs/S、全SM或实际末读设备。

[代码](check.py)首次 `--write` 排他创建[结果](results.json)，默认和 `--check` 只读复算。有理数核严格常数；有限BdG/Fock/来源诊断核共轭、Nambu半权、非Gaussian伙伴、异常腿及背景Ward，不替代连续证明。停止扫描最佳δ、更多代际和仪器精度；主线独审验收另存。
