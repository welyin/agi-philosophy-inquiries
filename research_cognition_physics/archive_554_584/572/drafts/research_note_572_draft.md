# 第572轮研究：共同规范物质与全局Einstein初始约束

日期：2026-10-01。接续571。八组复算、十六式；正式完成以research_round_572_checks.json为准。代码joint_gauss_einstein_initial_data.py，结果joint_gauss_einstein_initial_data_results.json。

## 0. 本轮真正接通的部分

571把给定光滑物质源映到精确Gauss图初值，但Gauss成立不保证同一来源能产生合法引力初值。本轮得到一对必须同时保留的结论：

1. **原571来源的密度动量不变时，不能放入平坦共形周期种子、Einstein框架常平均曲率（CMC）的引力分支。** 障碍是非零总动量，与调节势常数、共形因子或CMC大小无关。
2. **用原Higgs径向和singlet动量作有限反流，保留Gauss及原规范配置后，可以同时满足全部经典Einstein初始约束。** 在给定正F、指定有效势及足够大的CMC下，证明全环面存在光滑正共形因子；固定来源、CMC及无迹资料后，该因子唯一。

这完成的是同一限定经典模型的共同实现。维数、引力作用、匹配参数、边界和改变后的来源都是明确输入，没有从认知原则生成它们。反流不是571的趋零保真修正。

## 1. 成熟接口与历史去重

[Isenberg与Maxwell，2106.15027v2](https://arxiv.org/html/2106.15027v2)的定理4.4、§5.2和§6.3，以物质相空间资料组织共形方法，并列出紧致CMC的共形Killing兼容条件。本文继承这个方法，在具体非Abelian/Higgs源上独立核算密度、符号与积分。

[Choquet-Bruhat、Isenberg与Pollack，gr-qc/0610045v2](https://arxiv.org/pdf/gr-qc/0610045)的§2、§4.1及§5.2定理3提供约束分解和椭圆上下解方法。这里有多标量与Yang–Mills项，须重新检查共形权与上下解，不能直接搬其单标量分类。

549已给无规范场、零动量流的局部解析几何与参考构造；本轮保留571周期源中的非零规范应力，处理全环面的两类引力约束。550—553的裸曲率、真空匹配及运行限制不被撤销。共形方法、Fourier逆和最大值原理是成熟工具；本轮增量是具体来源的失败证书及改变来源后同模型的完整初值连接。

|层次|本轮地位|
|---|---|
|认知动机|将同一物质、内部约束和几何反作用共同交代|
|附加输入|四维二导数引力、三维平坦环面种子、全局平凡化、CMC及Einstein相空间匹配|
|继承对象|569的同尺度规范系数和完整正势、571的平滑规范及Higgs配置|
|解析结果|原来源受限不相容、保Gauss反流、全局正解及固定数据下的唯一性|
|数值验证|独立应力积分、约束残差、两种初猜和正性检查|
|物理解释边界|不是实际制备、时空维数选择、引力生成、量子引力或全宇宙解|

## 2. 固定的是密度动量，不能混同框架和速度

取五实分量$\phi^A=(\operatorname{Re}X_0,\operatorname{Re}X_1,\operatorname{Im}X_0,\operatorname{Im}X_1,s)$。给定Jordan二导数标量动能及规范不变的F>0，在Einstein框架有

$$
g_E=F g_J,\qquad
\mathcal K_{AB}=\frac{\delta_{AB}}F+\frac{3F_AF_B}{2F^2},\qquad
U=\frac V{F^2},\qquad
p_A=\sqrt\gamma\,\mathcal K_{AB}v^B,\qquad
F_Ak^A=0,\quad \mathcal K_{AB}k^B=\frac{k_A}F .
\tag{1}
$$

k是Higgs规范生成向量，v是Einstein法向速度。指定把571的数值动量识别为这里的canonical坐标密度p及规范电密度E；不是固定v，也不是声称Jordan与Einstein的完整标量动量天然不变。框架变换一般混合标量与引力迹动量；规范荷不受该径向混合影响，不能据此推出径向动量不变。

Einstein框架中的sigma动能、规范曲率平方和U没有度规导数，才可对接上述相空间方法。Jordan的F R不能直接作为“无度规导数的物质”套入其前提。规范Gauss仍是p·k加协变电散度，解几何不会自动改变这份密度约束。

沿571的Hermitian约定，$t_a=\sigma_a/2$，弱连接a吸收耦合；圆群用Q=6Y及$g_0=g_Y/6$，圆下标0不是时间。动量源密度为

$$
\begin{aligned}
D_iX&=\partial_iX+i a_i^at_aX+3ia_{0i}X,&
\mathcal F_{ij}&=\partial_i a_j-\partial_j a_i-a_i\times a_j,\\
f_{0ij}&=\partial_i a_{0j}-\partial_j a_{0i},&
M_i&=\operatorname{Re}\langle p_X,D_iX\rangle+p_s\partial_i s+
\sum_j E^j\!\cdot\mathcal F_{ij}+\sum_j E_0^j f_{0ij}.
\end{aligned}
\tag{2}
$$

在Gauss成立且周期无边界时，分部积分中的连接荷项与物质荷精确相消，得到

$$
P_i:=\int_\Omega M_i\,d^3x
=\int_\Omega\left[
\operatorname{Re}\langle p_X,\partial_iX\rangle+p_s\partial_i s+
\sum_jE^j\!\cdot\partial_i a_j+\sum_jE_0^j\partial_i a_{0j}
\right]d^3x .
\tag{3}
$$

这是总积分恒等式，不声称两种被积函数逐点相同。

## 3. 原来源为什么不能直接接入

沿549归一化，Hamiltonian约束为R−|K|²+τ²=2ρ，动量约束为div K−dτ=−M/√γ。取平坦周期种子及常τ：

$$
\gamma_{ij}=\psi^4\delta_{ij},\qquad
K_{ij}=\psi^{-2}\widetilde A_{ij}+\frac{\tau}{3}\psi^4\delta_{ij},
\quad \operatorname{tr}_{\delta}\widetilde A=0,\qquad
\partial^j\widetilde A_{ij}=-M_i
\ \Longrightarrow\ P_i=0 .
\tag{4}
$$

这里K是外曲率，与标量度量$\mathcal K$不同。任意无迹无散（TT）资料不能改变积分条件。

对571保存的$\Omega=(\mathbb R/2\pi\mathbb Z)^3$来源，$h=h_*(1+.05\sin(x+y))$、$s=s_*(1+.06\cos y)$、径向动量$.04\cos(x+y)$。径向项在x、y各给平均$.001h_*$；圆项$E_0^y=.07\cos x$、$a_{0y}=.08\sin x$另给x平均.0028。其余canonical积分为零。因此

$$
\frac P{|\Omega|}=(.001h_*+.0028,\ .001h_*,\ 0),\qquad
h_*=.6654429655855137,\qquad
P=(.859603867286,\ .165063269647,\ 0)\ne0 .
\tag{5}
$$

这严格排除“原密度来源＋此共形类＋Einstein CMC”的组合，不排除一般Einstein初值、其它种子或非CMC。固定配置时，总动量是动量的连续线性泛函；趋零L²动量修正不能消掉这份固定非零P，因此571的二阶修正也不能自动解决此障碍。

## 4. 同一物质的有限反流与动量约束

只调整原h、s的实径向动量：

$$
\delta p_h=-\lambda^i\partial_i h,\qquad
\delta p_s=-\lambda^i\partial_i s,\qquad
Q_{ij}=\int_\Omega(\partial_i h\partial_j h+\partial_i s\partial_j s)\,d^3x,
\qquad \lambda=Q^+P,\quad P\in\operatorname{range}Q .
\tag{6}
$$

径向实方向与所有规范生成方向实正交，故Gauss不变；总动量变为P−Qλ=0。一般源须核range条件。本例Q秩为2、P_z=0，xy块正定：

$$
\frac{Q_{xy}}{|\Omega|}
=\begin{pmatrix}.000553517926&.000553517926\\
.000553517926&.001089858828\end{pmatrix},\quad
\lambda=(11.481320394,-5.220560253,0),\quad
\|\delta p\|_{L^2}=\sqrt{P^TQ^+P}=3.001277173 .
\tag{7}
$$

拉格朗日乘子或Hilbert空间投影表明，这是仅允许这两类径向动量变化时，平直L²范数最小的反流。它不是实际几何动能最优解，更不是最小制备功；修正有限，原571径向状态已改变。

现在M光滑且均值为零。取零均值W，由每个非零Fourier模得到

$$
\widehat W(k)=\frac1{|k|^2}
\left(I-\frac{kk^T}{4|k|^2}\right)\widehat M(k),\qquad
\widetilde A_{ij}=\partial_iW_j+\partial_jW_i-\frac23\delta_{ij}\partial_\ell W_\ell
+\widetilde A^{\rm TT}_{ij},\qquad k\ne0 .
\tag{8}
$$

直接求散度得式(4)，解W的常量不影响$\widetilde A$。TT资料仍自由；数值例选择零TT，不能称所有引力资料唯一。Fourier零模只丢弃舍入量，解析零总量由式(6)承担。

## 5. Hamiltonian约束：共形权不能省略

固定上述配置、密度、$\mathcal K$及$\widetilde A$，用平坦种子缩并定义

$$
\begin{aligned}
A&=|\widetilde A|^2+p_A\mathcal K^{AB}p_B\ge0,&
B&=\sum_i\mathcal K_{AB}(D_i\phi)^A(D_i\phi)^B\ge0,\\
Y&=\frac{g_w^2}{2}\sum_i|E^i|^2+
\frac{g_Y^2}{72}\sum_i(E_0^i)^2+
\frac1{2g_w^2}\sum_{i<j}|\mathcal F_{ij}|^2+
\frac{18}{g_Y^2}\sum_{i<j}f_{0ij}^2\ge0,&
C&=\frac23\tau^2-2U .
\end{aligned}
\tag{9}
$$

物理标量法向能为$\psi^{-12}p\mathcal K^{-1}p/2$，梯度能为$\psi^{-4}B/2$，规范电磁能为$\psi^{-8}Y$。结合R=−8ψ⁻⁵Δψ及无迹外曲率范数ψ⁻¹²|$\widetilde A$|²，完整Hamiltonian约束等价于

$$
-8\Delta\psi+C(x)\psi^5-B(x)\psi-A(x)\psi^{-7}-2Y(x)\psi^{-3}=0 .
\tag{10}
$$

这些权来自密度与同一物理应力，不是把571的平直动能照搬到变几何上。方程中的2Y与式(9)本身的能量系数各有来源。

## 6. 一个明确的正F匹配与全局存在证明

数值见证另声明一个有效匹配系数$M_0^2=2$，不冒充549旧裸数值：

$$
F=M_0^2-\frac{|\phi|^2}{6},\qquad
\mathcal K^{-1}=F\left(I-\frac{\phi\phi^T}{6M_0^2}\right),\qquad
V=\frac14 d^T Ld,\quad d=(h^2-u_h,\ s^2-u_s)^T .
\tag{11}
$$

L、u及规范耦合全部取569—571同尺度参数；保留完整平方势，包括展开后的常数，不作额外真空能平移。F和该有效势的匹配地位与550—553裸作用限制分别记账。

原弱曲率$\mathcal F_{xy}^3=-A_0B_0(1+u\cos x)(1+.2\sin y)$处处非零，$A_0=.41,B_0=.29,u=.38$。解析保守界为

$$
F\ge1.86283354655>0,\qquad
0\le U\le U_+,\quad U_+\simeq9.21969440615\,10^{-5},\qquad
Y\ge\frac{[A_0B_0(1-u)\,.8]^2}{2g_w^2}
>.00413679447>0 .
\tag{12}
$$

U上界来自L的最大本征值及h²、s²的振幅界；不依赖网格求最大值。一般光滑正F固定源在紧域U有界。选择

$$
\tau^2>3\sup_\Omega U\quad\Longrightarrow\quad C_->0;
\qquad\text{本例取 }\tau^2=3U_++3\simeq3.000276590832 .
\tag{13}
$$

这是额外初值选择，不是宇宙学膨胀率或低曲率预测。以下是连续存在证明，与数值收敛分开：记各系数上下确界为下标±，反应项$P_x(t)=Ct^5-Bt-At^{-7}-2Yt^{-3}$。可取常数

$$
\ell=\min\left\{\frac12,\left(\frac{Y_-}{2C_+}\right)^{1/8}\right\},
\qquad
u=\max\left\{1,\left(1+\frac{A_++B_++2Y_+}{C_-}\right)^{1/4}\right\},
\qquad P_x(\ell)<0<P_x(u) .
\tag{14}
$$

下解符号乘以ℓ³即可核验；上解利用u≥1，将所有负项上界为u(A_++B_++2Y_+)。因此−8Δ+P在正区间[ℓ,u]有有序常数上下解。负幂在此区间无奇点；加足够正的线性项使迭代单调，周期椭圆上下解定理给正解，正则性提升给光滑解。A包含已解出的光滑纵向张量及所选光滑TT，紧域上其确界有限。

代码的上屏障使用格点极值，**只验证离散计算**；解析证明使用连续确界。不能把格点上界写成已认证的全域数值界。

固定全部这些数据，若两个正解ψ₁、ψ₂的比值在最大点取t>1，则

$$
\begin{aligned}
8(\Delta\psi_1-t\Delta\psi_2)
&=C(t^5-t)\psi_2^5+
A(t-t^{-7})\psi_2^{-7}+2Y(t-t^{-3})\psi_2^{-3}>0,\\
\Delta\psi_1-t\Delta\psi_2
&=\psi_2\Delta(\psi_1/\psi_2)\le0 .
\end{aligned}
\tag{15}
$$

−B项精确抵消，矛盾。交换两解即得相等。因此唯一的是**固定反流来源、τ及TT选择后的正共形因子**。

## 7. 一个必须保留的失败分支

若坚持同一平坦共形周期类并取最大切片τ=0，因本势U≥0、Y>0，对任何正ψ都有

$$
\int_\Omega\left[-2U\psi^5-B\psi-A\psi^{-7}-2Y\psi^{-3}\right]d^3x<0,\qquad
\int_\Omega\Delta\psi\,d^3x=0 .
\tag{16}
$$

故式(10)无解，非恒定ψ也不能例外。本结论仅限制该种子、势与切片组合；不排除一般GR。足够大的非零CMC给正面解，而非所有CMC均可或均不可。

## 8. 可复算证据

使用项目既有Python与NumPy运行：python -B -X utf8 joint_gauss_einstein_initial_data.py。

默认完整复算并与已保存JSON比较，不覆盖结果。代码还继承并核对入口探针及冻结571的依赖哈希。

|检查|保存结果与含义|
|---|---|
|原来源|总P与式(5)一致；独立差分协变／canonical积分差6.08×10⁻¹²|
|有限反流|Gauss变化0；剩余总动量各分量小于1.1×10⁻¹⁴；Q秩2，z方向反例保留|
|目标度量与正性|$\mathcal K\mathcal K^{-1}$误差6.42×10⁻¹⁸，正F、正Y及C>2检查通过|
|动量约束|16³谱散度残差1.39×10⁻¹⁷，无迹残差3.91×10⁻¹⁸|
|完整Hamiltonian|12³、16³、24³分别直接回代原约束，最大残差4.57、4.78、4.96×10⁻¹²|
|正解|24³上ψ在[.91127108,1.01841531]；两初猜在12³得到的解差6.67×10⁻¹⁶|
|最大切片对照|三个常数探针积分均负；覆盖任意正ψ的结论由式(16)承担|

独立二阶实空间差分对纵向张量的残差约.00290，属于该粗网格导数截断诊断，不冒充谱残差。数值采用Fourier配点及保正Newton求解；残差小不等于连续存在或误差定理，后者范围由式(14)—(15)确定。本例空间标量曲率约在[−1.22024,1.21381]，没有证明相对于截止尺度的曲率展开受控。

## 9. 合并现状、剩余输入与下一步

本轮把**同一Higgs＋singlet＋规范源、Gauss、完整经典Einstein初始约束**放进同一明确候选。571初值采样定理也可用于反流后的新光滑Gauss源；但不能因此声称每个离散格点已满足Einstein约束，或图的完整时间演化收敛到该几何。

尚未合并掉：引力作用及维数来源、有效参数匹配、量子连续极限、手征费米子、实际来源制备、自治装置、同一新解上的可用关系参考和低曲率范围。物质配置不变不代表参考坐标的时间导数不变，故549的满秩参考不能直接移植。

按用户提出的顺序，先采用可复用的共同实现和成熟定理，再集中于跨部门接口。下一入口优先核：**本轮共同初值的受控局部发展与同一物质参考能否共存**，先回查548—549及已有参考菜单；成熟适定性若直接覆盖只作映射，不另计轮次。需要新研究的是同一改变后来源的实际参考秩、共同应力和有效范围，而非重复扫描网格或任意另添参考场。

只在同一对象、同一参数及适用尺度确实相容时合并；“共同实现”不等于“条件等价”或“从认知原则必然推出”。整体统一目标保持，所有历史冻结和被排除分支保留。
