# 1185证明：同父结合能、互易来源与原物质钟的有限连接

2026-10-13。作者计算完成，待独立核验。条件结论在下述有限材料、准备域和实际实施合同内成立；新增采用公理为0。

## 1. 四维中介生成结合谱，精确保留未知输入

令材料 $t$ 与重中介 $b$ 均为qubit，

$$
A=(I+X_t)/2,\quad n=(I-Z_b)/2,\quad
h_\mu=I+I\otimes n+\sqrt2 A\otimes X_b.
\tag{1}
$$

定义

$$
r=\begin{pmatrix}\sqrt{2/3}&\sqrt{1/3}\\-\sqrt{1/3}&\sqrt{2/3}\end{pmatrix},\quad
W=(I-A)\otimes I+A\otimes r,\quad
E=W(I\otimes|0\rangle).
\tag{2}
$$

在 $A=0$ 块，$h_\mu=I+n$；在 $A=1$ 块，$I+n+\sqrt2X_b$ 的谱为 $0,3$。因此

$$
W^\dagger W=I,\quad E^\dagger E=I,\quad
W^\dagger h_\mu W=(I-A)\otimes I+(I+2A)\otimes n,
\quad h_\mu E=E h_c,\quad h_c=I-A=(I-X_t)/2.
\tag{3}
$$

整个微观谱为 $0,1,2,3$，半有界性不依赖制备码。负结合项 $-A$ 来自同一实际站内相互作用的本征谱；不是将结合质量作为另一自由数加入。$E\otimes I_R$ 对任意未知材料及参考严格等距。

对任一原 $t$ 算符 $O$，令

$$
\widetilde O=W(O\otimes I)W^\dagger,
\qquad\widetilde O E=EO.
\tag{4}
$$

这尤其适用于原空间耦合的 $Z_t$。若错误保留裸 $Z_t\otimes I$，则其压缩为 $\sqrt{2/3}Z_t$，码外泄漏算符范数为 $1/\sqrt3$。本轮两规则都采用式(4)，完整保留全部微观输出，不通过事后投影消除泄漏。定义 $N_b=W(I\otimes n)W^\dagger$；下述全部已运输端口、$h_\mu$ 及两B都与各站各复制的 $N_b$ 对易，故保留它的所有扇区，制备码是严格 $N_b=0$ 扇区。裸 $n$ 本来不守恒，不能混称它的真空不变。

## 2. 完整有限父过程及全母能源

保[1176](../1176/proof.md)的 $65^3$ 周期格点、每站源CAR4模、探针CAR12模和全部母扇区。$V=65^3,\ell=16/65$，物质质量与空间二次量子化来源仍为

$$
F_{si}=\epsilon m_i+\nu\sum_bt_{bi},\quad F_{ti}=g(t_{xi}-t_{yi}),
\qquad(\epsilon,\nu,g,\Omega)=(.15,.35,.2,.5).
\tag{5}
$$

每复制 $a$ 的自由项改为 $h_a^\mu=H_{S,0,a}-\Omega\sum_iX_{si,a}+\sum_i h_{\mu,ia}$。原物质跨复制作用为

$$
W^0_{ab}=\sum_i\left[
Z_{si,a}F^S_{si,b}+F^S_{si,a}Z_{si,b}
+\widetilde Z_{ti,a}F^S_{ti,b}+F^S_{ti,a}\widetilde Z_{ti,b}\right].
\tag{6}
$$

固定两规则 $r\in\{m,w\}$：

$$
B^\mu_m=\epsilon h_\mu,\qquad B^\mu_w=\epsilon I,
\qquad
H^r_{BG}=\sum_ah_a^\mu+\frac1K\sum_{a<b}\left[W^0_{ab}+
\sum_i(Z_{si,a}B^\mu_{r,ib}+B^\mu_{r,ia}Z_{si,b})\right],
\tag{7}
$$

$$
H^r_\mu=C_{\mathrm{new}}I+H^r_{BG}+H_{P,0}
+\sum_i(F^P_{si}\overline Z_{si}+F^P_{ti}\overline{\widetilde Z}_{ti}).
\tag{8}
$$

复制标号不同于空间位置；不含 $a=b$。全部CAR作用为偶算符，空间支持保持原邻边。各项有限、自伴、时间无关，完整总能源守恒。准备与末读是另有资源和误差的有限实施合同，式(8)不自动制造这些仪器。

令 $C_K$ 为1176式(4)的全CAR母空间正谱常数。$\|h_\mu\|=3$、$\|\widetilde Z_t\|=1$，两规则均有 $\|B^\mu_r\|\le.45$。故统一选择

$$
C_{\mathrm{new}}=C_K+(3-\Omega)VK+.45V(K-1)
\tag{9}
$$

足以沿原逐项全范数估计使式(8)正。原每复制 $t$ 自能预算 $\Omega$ 换为3；新增互易范数不超过 $.45V(K-1)$。这是重做三角范数账，不是声称新旧自能的差算符范数至多 $3-\Omega$。这不是只在码内作正性判断，也没有新增振子Fock截断。两实际规则、编码比较器及后述背景比较均保同一 $C_{\mathrm{new}}$ 相位。未假设这一外加全局正谱基线本身普适参与几何来源。

令 $\mathcal E$ 为全部站点/复制的 $E$ 张量积，CAR、$G_s$、探针及参考上作用恒等。由每一局部式(3)—(4)，在**全部原CAR母扇区**有

$$
H^r_\mu\mathcal E=\mathcal E H^r_{\mathrm{eff}},\qquad
e^{-itH^r_\mu}\mathcal E=\mathcal E e^{-itH^r_{\mathrm{eff}}},
\tag{10}
$$

其中有效 $t$ 自能为 $h_c$，有效互易来源为

$$
B_m=\epsilon h_c=\epsilon(I-X_t)/2,\qquad B_w=\epsilon I.
\tag{11}
$$

因为父H自伴，$\operatorname{Ran}\mathcal E$ 是约化子空间。任意原端口、记录和声明的后继都须以同一 $W$ 运输；这种数学交织不另授未经声明的真实控制权限。式(10)强于平均能量或单态匹配，带任意旁观参考也成立。有效自由项中的 $+.5VK$ 常数及其相位原样保留；匹配来源中的 $+.075I$ 也原样保留，二者不能用“中央项不改通道”一并删除。

## 3. 来源由同一H独立给出，中央能源不删除

以下采用原已知源CAR真空，它精确不变，式(6)在其上为0。微观准备码同样精确不变，故可使用式(10)的完整有效量子父过程，不能先删中介再重新猜来源。

例如原 $G_s$ 的Heisenberg关系在真空扇区为

$$
\dot Y_{si,a}=2\Omega Z_{si,a}
+2X_{si,a}\left[\frac1K\sum_{b\ne a}B_{r,ib}
+\frac1K F^P_{si}\right].
\tag{12}
$$

实际来源含完整复合来源及探针反作用；正/删弱差异直接来自式(7)，没有先把总静止能定义为来源。所有源、关系自能、交互与后态在式(8)中保留。

每个已知背景cluster准备 $G_s=+X$、$\mathbf g_t=(q,0,\sqrt3/2)$，$q=\pm1/2$，再经同一 $\mathcal E$ 装入。Hartree背景只作比较器。其 $x_t=q$ 守恒，$G_s$ 的有效生成元是 $-\Omega X_s+bZ_s$，

$$
b_{m,q}=\epsilon(1-q)/2,
\qquad b_{w,q}=\epsilon.
\tag{13}
$$

匹配支的 $t$ 生成元为 $[1+\epsilon z_s(t)]h_c$；删弱支为 $h_c+\epsilon z_s(t)I$。二者自由复合谱完全相同，公共响应不同。已知 $q=+1/2,-1/2$ 本来有不同复合能源，不宣称两准备同总能谱；固定q跨两规则也不宣称完整相互作用Hamiltonian同谱。

每站Hartree总能为 $-\Omega x_s+(1-q)/2+z_s b$，由同一变分进动守恒。该恒等式只是完整自治守恒的比较器版本。中央 $1/2$ 给出真实 $\epsilon/2$ 来源，不能套用1178的 $E_0=0$ 数值，也不把它解释为已推得宇宙学常数。

## 4. 四个固定原钟任务

取 $K=10^{14}$、$T=\pi/2$。两规则各作 $q=\pm1/2$，不挑选时间和读口。由常场进动，令 $r_b=\sqrt{\Omega^2+b^2}$，

$$
z_s(t)=-\frac{\Omega b}{r_b^2}[1-\cos(2r_bt)],\qquad
\tau_b(T)=T-\frac{\epsilon\Omega b}{r_b^2}
\left[T-\frac{\sin(2r_bT)}{2r_b}\right].
\tag{14}
$$

完整探针比较器仍用同一原生成元

$$
h_{P,p}=N\beta M_P+\sum_jv_j\alpha_j k_j,\quad N=1+\epsilon z_s,
\quad v_{x,y}=1+\nu z_s\pm gz_t,\quad v_z=1+\nu z_s.
\tag{15}
$$

数值见证选原48维码中的 $p=0$、$\beta=+1$、两质量 $.3/.5$ 的等幅态。原 $P_{-Y}$ 给

$$
p_D(b)=\frac{1+\sin(.2\tau_b)}2.
\tag{16}
$$

均匀 $p=0$ 波包的空间站点细记录完整保留。此见证没有宣称位置边缘分离、单条共位世界线或新曲率重建；全部未知48维输入的保护由下一节给出。正/删弱中的空间端口均保持式(4)，不能因见证是 $p=0$ 就删除空间耦合。

固定事件差之差

$$
\mathcal D=\left(p_{D,m,+}-p_{D,m,-}\right)
-\left(p_{D,w,+}-p_{D,w,-}\right).
\tag{17}
$$

比较器数值中心为：匹配两支概率 $.653898079035356,.652686687621766$；删弱两支均约 $.652090769074997$。因此 $\mathcal D_{\mathrm{comp}}\simeq .0012113914135893$。这是固定解析任务的值，最终严格判据用[代码](check.py)和[结果](results.json)的有理包络，不以显示小数承担证明。

## 5. 有限K、完整未知输入与仪器预算

必须先证明式(10)，再使用 $\|B_m\|,\|B_w\|\le.15$ 的有效范数。两规则均可取

$$
w_{\mathrm{site}}=.3,\quad
a_s=4.340625,\quad a_t=1.625,\quad
\mathcal B^2=7a_s^2+5a_t^2.
\tag{18}
$$

$.3$ 不是完整未压缩中介的互易范数；全母正性已单独采用 $.45$。源真空下不同空间站的精确背景独立，站内的s/t及各复制通常有关联。复用1178的计数证明，围绕精确均值 $\mu_c$ 分解 $\overline Z_{ci}-z_c=(\overline Z_{ci}-\mu_c)+(\mu_c-z_c)$；跨站零均值涨落正交，同站s/t用Cauchy，均匀偏置用完整一粒子 $\|\sum_iF^P_{ci}\|\le a_c$，不作体积倍增。

以下明确约定 $U^r_{BG,K}(t)=\exp[-it(C_{\mathrm{new}}I+H^r_{BG,\mathrm{eff}})]$，其中 $H^r_{BG,\mathrm{eff}}$ 是式(7)严格编码后的背景H；$U^r_{P,q}$ 只由式(15)生成。这样比较器合计带实际父过程的同一全局常数相位，且背景自由项内的 $+.5VK$ 没有删除。

对全部原48维输入和任意 $R$，Duhamel给准备等距界

$$
\left\|U^r_\mu(T)\mathcal E(\Phi_{0,q}\otimes J_P)
-\mathcal E\left[(U^r_{BG,K}(T)\Phi_{0,q})\otimes U^r_{P,q}(T)J_P\right]\right\|
\le\eta_{P,K}
\le\frac{4Te^{6w_{\mathrm{site}}T}}{\sqrt K}
\left(\sqrt2\mathcal B+a_s+a_t\right)<.000245.
\tag{19}
$$

其中 $J_P$ 是原48码准备等距；两个比较分量合计保同一全局常数相位。右侧背景是**精确量子整体** $U^r_{BG,K}\Phi_{0,q}$，没有用全背景Hartree产品替代它。真实传播保探针反作用、跨质量/动量相干、参考和全部微观资源。

统一每个任务的实际装入、端口和完整终端记录误差 $\eta_{\mathrm{inst}}\le.00002$ 是新增有限实施合同。它包括新增E准备、W运输后的实际耦合/指针，以及声明有限任务中使用的许可未来，不能视为旧权限自动覆盖。原物质钟/站点仪器作用不删任何背景；新增微观控制必须在同一预算中交代。后接同一完整CP记录操作保持误差，未认证的新中途干预不能直接继承式(19)。未知资料只出现在一个实际探针及其参考，已知q背景复制不构成未知态克隆。

每个概率总误差可取 $\delta=.000265$。利用四个任务，

$$
\mathcal D_{\mathrm{actual}}
\ge\mathcal D_{\mathrm{comp}}-4\delta
>.00121139-.00106
>.00015.
\tag{20}
$$

作者代码已给三角函数和指数的有理包络，净下界的显示值为 $.00015139141358929263$；严格断言使用结果中的精确有理字段，待独立审阅。数值舍入不额外挤占仪器合同，经验抽样若开展另付统计误差。全部 $p_D\pm\delta$ 在1178事前 $[1/2,(1+\sin(.38))/2]$ 窗内；这避免以相位别名选择错误钟分支，虽式(17)本身不需要反解时刻。小规模 $K=2$ 的 $128\to32$ 准备交织只作全父矩阵诊断；实际 $K=10^{14}$、完整48码及R的资格来自式(10)、(19)和有理包络，不谎称直接模拟了巨大父态。

## 6. 新结果、剩余自由及停止线

本轮净结果是：由一份有限微观相互作用实际生成的结合项，可通过独立原公共关系来源、自主后态和同一原物质钟给出严格有限差之差。精确编码只承担完整接续；超误差原钟任务才完成本页此前缺少的来源读回。

仍输入了微观材料及谱、两类响应规则、原空间父模型、已知背景与实际微观权限。匹配规则未从四原则/CO/SR/CC中唯一推得，中央能源和其他响应自由未全部消除。没有新选择 $H_G$、非线性约束、普适引力自耦或Einstein方程。三项GR硬交付仍未全部完成。

不继续调 $K$、钟频或中介谱，不改旧轮及历史数值。保留额外模式和暗部门，后续研究按[NEXT](NEXT.md)决定。
