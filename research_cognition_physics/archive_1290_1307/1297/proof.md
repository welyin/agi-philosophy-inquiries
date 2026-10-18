# 第1297轮证明：非零整数荷的同父中性区域记录

2026-10-18。依据[事前准入](../_shared/round1297_admission.md)。本轮直接复用[1296](../1296/proof.md)的连续传播、正常／反常两腿、未知伙伴及有限源域；新增的是将该具名合同接到1035／1292的非零整数荷分支。一般Wilson运输、CAR协变性和复合中性不另算新发现，不采用新公理。[1055后子系统审计§4—6](../../archive_1046_1062/_admission/physical_subsystems_after1055/selection.md)已明确固定局部荷双轨可以承载复量子资料，同时禁止把中性材料证书无证移到带电Gauss场。本轮只补当前背景二次层的具名两球／全窗／来源连接，不领取该页§6更强的动态Gauss与完整物理仪器合同。

## 1. 固定父作用、荷与背景

取1296的签名 $(+---)$、自然单位、$n=3$，固定

$$
\alpha_*={31\over2^{27}},\quad \epsilon={\alpha_*\over128},\quad
\delta={\epsilon\over2},\quad
Y_\nu=I_3+\epsilon(D_3+P_3),\quad Y_e=2I_3+\epsilon D_3,
\quad H_0=(0,1)^T,\quad S_0=1,\quad A_\mu=0 .
\tag{1}
$$

$D_3=\operatorname{diag}(1,2,3)$，$P_3$为路径邻接矩阵。原全部物种、同一规范／标量动能、非负径向势、固定物理耦合和重整化规则保留。取[1292式(2)—(3)](../1292/proof.md)的 $k=0,1$ 两支；它们是同父作用类型的不同表示，不是在一份模型内切换电荷。全字段满足

$$
Y_k=Y_0+3k(B-L),\qquad Y_S=-6k,\qquad
-\tfrac12(SN^T\delta I_3N+\mathrm{h.c.}) .
\tag{2}
$$

$S$始终是复字段，两个实方向及来源均保留。含$N$的费米反常按[1060](../../archive_1046_1062/1060/proof.md)原域继承，不重复计算。$S_0\ne0$已选择$r:S\mapsto-S$的破缺背景；并非固定背景仍有全部连续超荷对称。规范变换必须连背景一起运输。

在$H\ne0$的平凡局部补丁定义

$$
\chi_L={\widetilde H^\dagger L\over |H|},\qquad
\chi_R=N^\dagger,\qquad Y(\chi_L)=Y(\chi_R)=q=-3k .
\tag{3}
$$

这是同一经典背景的精确二次字典，不是动态量子Higgs复合场的精确CAR声明。$U_H=(\widetilde H,H)/|H|$连同带电分量完整保留；普通非均匀背景中的$U_H^\dagger D U_H$不能删去。本轮只使用(1)及其规范轨道，参见[已登记背景字典](../_admission/background_dressing_screen_after1293.md)。整数归一$6Y(\chi)=-18k$使该弱／色单态表示下降到原$G_6$中心商；没有另造带荷物种。原$r$对两手均给$-i$，双轨仍同宇称。

在(1)代表上保留费米二次生成元与1296完全相同：规范势为零，$Y,H,S,\delta$相同。相同的是这份费米生成元和其参考；不是完整作用、规范涨落谱、全部绝对来源基线或更广泛物理预测相同。

## 2. 球内树运输与真正中性菜单

保原两球半径$a=1/32$、中心距离$L=1/4$及$I=[7/16,1/2]$。每球$C=A,B$固定内部基点$x_C$，从基点到$x$的直径向路径始终在本球。记$a_\mu=g'B_\mu$，取$D=\partial-iqa$，并定义

$$
W_{q,C}(x)=\exp\left[-iq\int_{x_C}^{x}a\right],\qquad
\widetilde\chi_C(x)=W_{q,C}(x)\chi(x).
\tag{4}
$$

对$a\mapsto a+d\Lambda$、$\chi\mapsto e^{iq\Lambda}\chi$，有

$$
W'_{q,C}(x)=e^{-iq(\Lambda(x)-\Lambda(x_C))}W_{q,C}(x),\qquad
\widetilde\chi'_C=e^{iq\Lambda(x_C)}\widetilde\chi_C .
\tag{5}
$$

径向线积分是球内光滑函数，不引入中心坐标奇点。经典相位点乘酉不扩大支持、保持等时CAR。用1296原实紧支归一包和旋量，将两味右手输入模记为$a_1,a_2$，左手接收模记为$b_1,b_2$。同球两个模只承受相同基点相位。因此单占据矩阵单位$|10\rangle\langle01|$及其对角元、$b_i^\dagger b_j$均中性；$E_1=b_1^\dagger b_1$是$0\le E_1\le I$的球内偶效果。

这是**区域**菜单：双线性核含两条在球内闭合到同基点的路径，不是点局域的单费米插值场。它不需要A、B之间新增带荷相位参考。对动态量子连接不能把(4)当无成本实体Wilson装置；本轮仅处理给定经典背景及其字典。

## 3. 完整初态、伙伴与规范轨道

保1296的真实两模码$W|0\rangle=|10\rangle,W|1\rangle=|01\rangle$。Nambu反线性$\Gamma$及秩4投影$\widehat F=F+\Gamma F\Gamma$均保留；这仍是两个物理模。补代数取原自身Hadamard参考$P_\delta$的**全部**限制$\sigma_\perp$，包括反常关联。任意输入和伙伴按

$$
\mathcal E(\rho_{QR})=(W\otimes I_R)\rho_{QR}(W^\dagger\otimes I_R)
\ \widehat\otimes\ \sigma_\perp
\tag{6}
$$

定义完整正常态。它一般非Gaussian，不能仅由二点函数作Wick重建。远端球的初片偶代数位于补空间内，所以其与$R$的完整限制为原参考限制与$\rho_R$的乘积；不是假设原真空因子化。全部后态按完整演化保留。

码内两基向量具有同一模式荷；其密度中性，**不等于整场替换态是固定总荷本征态**，也不表示以后始终保单占据码。准备的各Kraus若重置原参考的零／双占据部门，可以转移荷。式(6)是公开初态CP合同；若调用其物理实现须保荷库存／背景交换与能源，不能把中性密度叫作免费制备。本轮不要求另造工厂，亦未签其自治实现。

以下$G(t)$是完整背景规范变换的字段表示，不是式(4)的球内Wilson读口字典。对时间依赖规范运输，令$G(t)$为完整原字段的点乘表示，$\mathbb G=\operatorname{diag}(G,\bar G)$。它与$\Gamma$相容。以主动字段变换约定，

$$
\mathbb h'=\mathbb G\mathbb h\mathbb G^\dagger+i\dot{\mathbb G}\mathbb G^\dagger,
\quad \mathbb U'(t,s)=\mathbb G(t)\mathbb U(t,s)\mathbb G(s)^\dagger,
\quad P'(t)=\mathbb G(t)P_\delta\mathbb G(t)^\dagger .
\tag{7}
$$

$S$同时取$S'=e^{-6ik\Lambda}S$，Higgs、模式、码、补态及来源均运输。这里$P'(t)$是同一参考的运输，不是重新选$\mathbf1_{(-\infty,0)}(\mathbb h'(t))$；后者通常不同。$\Gamma$不变码的CAR同构及伙伴恒等扩张保持(6)、初远端盲性与Hadamard／正常性，不要求在一个预先固定的全局Fock表示中另实现规范控制器。固定$S_0$仅有其残余对称；不能把(7)误读为固定背景全$U(1)$不变。

## 4. 正常／反常端点相位与全窗记录

1296的接收输入差为实线性映射

$$
\mathcal L_\delta(X)=A_\delta XA_\delta^\dagger-B_\delta X^TB_\delta^\dagger .
\tag{8}
$$

它不是单独CP通道，合法完整输出还含环境正常／反常基线。记$\lambda_A=\Lambda(0,x_A)$、$\lambda_B=\Lambda(t,x_B)$，则(5)给

$$
A'_\delta=e^{iq(\lambda_B-\lambda_A)}A_\delta,\qquad
B'_\delta=e^{iq(\lambda_B+\lambda_A)}B_\delta .
\tag{9}
$$

因此(8)、全部同球偶二点资料及其伙伴扩张不变。正常腿与反常腿的端点相位不同，不能只运输$A$而冻结$B$或补态。完整偶效果和态同构也保持所有高阶统计，证明不依赖仅比较协方差。

在代表(1)上$W=1$，$\delta=\epsilon/2$在1296已证区间内。直接继承其保守常数，不另优化：

$$
\Delta p_k(t)>{61951\over65536}\alpha_*^2>{15\over16}\alpha_*^2,
\qquad
\|\mathcal L_\delta(X)\|_2>{3875\over4096}\alpha_*^2\|X\|_2
\quad(X=X^\dagger\ne0),\quad t\in I .
\tag{10}
$$

$X=0$只用非严格等式。两式对$k=0,1$及各自完整规范轨道成立。任意事前非负归一时间权平均保概率下界；同码量子相干亦到达接收二点。此为结构可辨与信息单射，不是无噪恢复、已制造层析或实际采样误差达标。不把任意真实非规范等价背景或全SM高阶余项设成零。

## 5. 独立全源身份与背景Ward

对每支先在**原**同一局部作用中固定物理耦合，变分全部$g,A,H,S$，再代背景：

$$
\delta S_F^{(2)}=\int\sqrt{-g}\left[-\tfrac12T^{\mu\nu}\delta g_{\mu\nu}
+j_a^\mu\delta A_\mu^a+\mathcal F_{H,I}\delta H^I+\mathcal F_{S,J}\delta S^J\right].
\tag{11}
$$

四实Higgs与两实S方向、原旁观字段、标量／规范背景来源、固定反项及参考绝对基线都保留。任意Yukawa spurion变分不属于此菜单。相同费米传播不使$j_Y$相同：在完整原字段表示上$Y_k=Y_0+3k(B-L)$；第一超荷插入有对应荷权变化，二次规范接触响应也须重算。这里签的是**来源插入及响应字典不同**，不宣称每个当前码态、每个源测试的均值都不同，更不称其全部动力学等价。

例如原物理中微子左／右超荷为$-1/2-3k,-3k$；归一化左手字典的$H$依赖不能先删。原背景质量矩阵满足

$$
i[\mathbb h,\mathcal Q_Y]+\delta_Y\mathbb h=0,\qquad
\delta_YH={i\over2}H,\quad \delta_YS=-6ikS .
\tag{12}
$$

$k=1$的Majorana腿需要$S$相位源，Dirac腿需要Higgs相位源。只保费米电流而删标量交换就破坏Ward；这不是原总规范荷守恒失败。微分同胚Ward同样保背景交换；在本常代表上等待阶段可用原平移守恒。

令紧支测试$\zeta$的原物理源核为$\mathfrak j^{(k)}_\zeta$，$\Gamma\mathfrak j\Gamma=-\mathfrak j$。沿1296同一实际导数域，Nambu半权给

$$
\mathcal E^*(J^{(k)}_\zeta)=c^{(k)}_\zeta I_2+V^\dagger\mathfrak j^{(k)}_\zeta V,
\quad c^{(k)}_\zeta=\omega_k(J^{(k)}_\zeta)
+\tfrac12\operatorname{Tr}\mathfrak j^{(k)}_\zeta
(-\widehat FP_\delta-P_\delta\widehat F+\widehat FP_\delta\widehat F+F_c).
\tag{13}
$$

$F_c=\Gamma F\Gamma$。迹只涉及光滑有限秩腿；一般来源是逻辑矩阵并与$R$关联。规范轨道上固定物理源的拉回与态同时运输，(13)仍成立。有限性复用[805](../../archive_797_818/research_note_805.md)；有界传播误差不控制无界全源或噪声，没有全整数$k$统一来源预算。

由于字典依赖背景，直接在新字段固定时求导还含Euler、边界和响应接触项。尤其

$$
\delta_a W_{q,C}=-iq W_{q,C}\int_{x_C}^{x}\delta a
\tag{14}
$$

进入**准备／模式／读口**的显式背景响应，不能冒充(11)的物理电流。源先变分再拉回；不把(14)丢掉后签总响应相同。时间连接、Higgs取向及相位、二阶变分接触均按[背景运输筛查§4](../_admission/background_dressing_screen_after1293.md)保留。本轮没有额外声称路径积分测度变化恒为零。

## 6. 有限资源、删除与净结论

代表上的$\widehat F\mathbb h\widehat F=0$及原径向包不变，故1296的独立能源计算直接给

$$
E_{\rm rel}={1\over2}\operatorname{Tr}(\widehat F|\mathbb h|)
<2\sqrt{524288+4}<1450 .
\tag{15}
$$

这是每份准备相对其费米参考的有限能源，等待守恒，非制备总成本或局部正能。规范轨道若要比较物理能源，须运输原时间／能流任务；不能把时间依赖规范代表中裸$\mathbb h'$的瞬时期望当同一能源。其它绝对背景来源按(11)另保，不能因(15)输入无关便把全部来源称标量。

|删除或改变|实际影响|
|---|---|
|要求单N点局域不变线性插值|1292§6的非零k障碍仍成立；本区域偶菜单不属于该任务|
|双轨不再同荷|未知轨间相干及矩阵单位一般不再中性，需另给补偿关系|
|冻结S／参考却任意变换连接|不再是同一规范轨道，(7)—(10)不可照搬|
|保概率但省略S相位源或模式背景导数|分别破坏原Ward或完整响应字典|
|把背景场也量子化并要求动态Gauss实现|超出本有效合同；本轮未签其准备和仪器|

**限定结论。** 原SNN父作用在$k=0$和$k=1$都具备同一具名有限背景记录合同及各自完整源身份。单场中性不是这份共同记录的必要条件；该任务本身不选$k=0$。这是跨原结构分支的同对象资格连接，不证明所有k等价、全O完成或全部认知公理不足。正常态、补空间、伙伴、来源和旧连续窗口共同绑定；没有新设备、群或质量扫描。

[代码](check.py)的有理数核继承常数，有限矩阵核荷／相位、Nambu运输、完整CAR与伙伴和背景来源。有限诊断不代替连续传播证明，亦不代表制造末读器。首次 `--write` 排他生成[结果](results.json)，默认及 `--check` 只读复算。
