# 1053解析证明：同一真实Coulomb仪器的电子Newton来源字典

2026-10-08。物理采用和任务范围属于定理前提。本轮沿用[1050证明](../1050/proof.md)的原参数和真实传播，不新增控制器或几何探针。[复算代码](coulomb_newton_source.py)校准源积分、矩阵与有理证书，不数值求解连续Coulomb传播。

## 1. 同一父过程与严格范围

原子单位下，两个固定单位正核在 $R_L=-Re_x/2,R_R=Re_x/2$，电子Hamiltonian为

$$
H=\sum_{i=1}^2\left(-\tfrac12\Delta_i
-|r_i-R_L|^{-1}-|r_i-R_R|^{-1}\right)+|r_1-r_2|^{-1}.
\tag{1}
$$

完整固定核Born–Oppenheimer能源另含 $1/R$。两正常Slater轨道、Löwdin正交化、全部两电子反对称初态空间为

$$
\phi_{L,R}=\pi^{-1/2}e^{-|r-R_{L,R}|},\quad
S=e^{-R}(1+R+R^2/3),\quad G=I+S\sigma_x,
$$


$$
W_1=(\phi_L,\phi_R)G^{-1/2},\quad
W=\bigwedge^2(W_1\otimes I_{\rm spin}):\mathbb C^6\to\mathcal H.
\tag{2}
$$

采用1050原辅助qubit、初态 $|+\rangle$、脉冲

$$
H_F(t)=H\otimes I+\beta(t)X\otimes P_1,\quad
X=x_1+x_2,\quad
\beta(t)={2\pi\over\ell T}\sin^2{\pi t\over T},
$$


$$
\ell=R/\sqrt{1-S^2},\quad R=1000,\quad T=1000/137.
\tag{3}
$$

这里 $U_F(t)$ 是(3)真实连续Hamiltonian的酉传播；W并非H谱子空间。记 $Q=(n_R-n_L)/2,D=Q^2,V=W\otimes|+\rangle,J=U_F(T)V$。1050已证明共同域、全部未知六维输入与任意被动参考上的

$$
\|J-J_{\rm ideal}\|<\epsilon_0={74239\over335650},
\quad K_{+,\rm ideal}=e^{-iHT}W(I-D),
\quad K_{-,\rm ideal}=e^{-iHT}WD.
\tag{4}
$$

所列裸辅助H为0；核支持、外部受控偶极和末读采用仍在。**本轮另采用领先非相对论静质量及线性Newton约束；只证明电子贡献的来源测试系数，不领总装置来源、实际引力探测器或物理有限G反作用余项。** 既有三维、电子参数、Coulomb作用和Newton约束未由认知原则推出。

## 2. 平滑Newton约束对偶测试的明确归一化

固定测试尺度 $a=R/2=500$，定义

$$
k_a(r)=(|r|^2+a^2)^{-1/2},\quad
f(r)=k_a(r-R_R)-k_a(r-R_L),
$$


$$
A={a\over2}\,[f(r_1)+f(r_2)],\quad
d=1-{1\over\sqrt5}.
\tag{5}
$$

由 $k_a\in[0,1/a]$ 知 $\|A\|\le1$，且中心值为 $f(R_R)=d/a,f(R_L)=-d/a$。a是选择的有限测试尺度，不是截断或物理最小长度。

正密度

$$
w_a(r)={3a^2\over4\pi(|r|^2+a^2)^{5/2}}
$$

满足 $-\Delta k_a=4\pi w_a$，径向累积积分为
$r^3/(r^2+a^2)^{3/2}$，总积分为1。故在平直固定核参考坐标和无穷远衰减边界下，$k_a$为Newton Green核与该Plummer测试密度的卷积。将电子静质量密度的领先Newton势分别用以两个核位置为中心的 $w_a$ 测试，得到

$$
\boxed{\overline\Phi_{e,R}-\overline\Phi_{e,L}
=-{2Gm_e\over a}\langle A\rangle.}
\tag{6}
$$

等质量固定核对这份差值的领先静质量贡献严格相消，前提是使用同一核锚和同一测试权重。这里的势差不是未选参考的GR关系可观测量；无穷远边界和Newton规范已采用。控制、支持、指针等来源没有因此相消。电子动能、结合能的后Newton修正也未由(6)认证。

式(6)可读作算符值线性约束的对偶系数，或对给定电子态的领先平均源预测；本轮不指定一个全量子几何模型。因而它不宣称实际测得几何的完整Kraus后态。带原指针标签的未归一源矩按同一个线性式变换即可。Plummer权重有无穷空间尾；a有限不等于紧支撑或T内实际可访问，本轮未实现全空间积分检测。

## 3. 初态全六维字典：Löwdin中心项必须保留

令 $Z=\operatorname{diag}(-1,1)$。核导数直接给

$$
\|\nabla f\|_\infty\le L_f={4\over3\sqrt3a^2},\qquad
\|\Delta f\|_\infty\le{6\over a^3}.
\tag{7}
$$

单核核函数的梯度最大值在 $|r|=a/\sqrt2$，等于 $2/(3\sqrt3a^2)$。每列Slater态有 $\|(r-R_{L,R})\phi_{L,R}\|=\sqrt3$。写 $V_1=(\phi_L,\phi_R)$ 并将乘法算符f分成中心值和残差，得

$$
fW_1-W_1(d/a)Z
=[fV_1-V_1(d/a)Z]G^{-1/2}
+V_1[(d/a)Z,G^{-1/2}].
\tag{8}
$$

两列残差的算子范数至多 $\sqrt6L_f$，$\|G^{-1/2}\|=(1-S)^{-1/2}$。同时 $\|V_1\|=\sqrt{1+S}$，

$$
\|[Z,G^{-1/2}]\|=(1-S)^{-1/2}-(1+S)^{-1/2}
\le{S\over(1-S)^{3/2}}.
$$

所以一粒子误差至多
$r_*\sqrt3L_f+(d/a)S\sqrt{1+S}/(1-S)^{3/2}$，其中 $r_*=\sqrt{2/(1-S)}$。升为两个含自旋的反对称粒子后，二项乘法和的范数误差至多两倍；结合A的a/2因子，

$$
\boxed{\|AW-WdQ\|\le\delta_0
={4r_*\over3a}+d\,{S\sqrt{1+S}\over(1-S)^{3/2}}.}
\tag{9}
$$

特别地，(9)控制离开轨道码的成分，不只是 $W^\dagger AW$ 压缩矩阵的误差。

## 4. 同一实际驱动后的来源字典

沿用1050的真实能源界

$$
H+32\ge T_{\rm kin}/2,\qquad
E_*:=\|W^\dagger(H+32)W\|\le32+7r_*<42,
$$


$$
Y(t)=\langle H+32\rangle_t,
\quad \sqrt{Y(t)}\le\sqrt{E_*}+\sqrt2L_B(t),
\quad L_B(t)=\int_0^t\beta(s)ds.
\tag{10}
$$

这些界对所有初态及被动参考一致。A和 $X\otimes P_1$ 均为位置乘法，所以
$[H_F(t),A]=[H,A]$。对每个电子

$$
[-\Delta_i/2,f(r_i)]\Psi
=-\nabla f(r_i)\cdot\nabla_i\Psi-\tfrac12\Delta f(r_i)\Psi.
$$

又 $\sum_i\|\nabla_i\Psi\|\le2\sqrt{2Y}$，故

$$
\|[H,A]\Psi_F(t)\|
\le aL_f\sqrt{2Y(t)}+{3\over a^2}
={4\sqrt{2Y(t)}\over3\sqrt3a}+{3\over a^2}.
\tag{11}
$$

一、二阶导数有界的f保持共同Coulomb域；以1050有限能源控制，(11)在真实传播上可积分。对(8)原同一对称非负脉冲，

$$
L_B=\pi/\ell,\qquad \int_0^T L_B(t)dt=TL_B/2.
$$

对Heisenberg差积分并加(9)，得到本轮正向主结论

$$
\boxed{\|AJ-JdQ\|\le\delta}
\tag{12}
$$


$$
\delta\le\delta_0+T\left[
{4\sqrt{2E_*}\over3\sqrt3a}+{3\over a^2}
+{4L_B\over3\sqrt3a}\right].
\tag{13}
$$

这里 $J=U_F(T)(W\otimes|+\rangle)$，是**真实驱动后的全输入来源字典**。没有用旧ε去控制无界应力，也没有令实际传播停在W。去掉(13)末项即得原自然后态 $e^{-iHT}W$ 的界；这只说明两种实际传播均可测试，主结论始终用真实 $U_F$。

张量任意被动参考不改算子范数。对原指针投影 $P_\pm$ 令 $J_y=P_yJ$；由于A与指针对易，(12)给完整堆叠字典

$$
\left\|\bigoplus_y(AJ_y-J_ydQ)\right\|\le\delta.
\tag{14}
$$

更具体地，令真实逻辑效果 $E_y=J_y^\dagger J_y$、源效果 $F_y=J_y^\dagger AJ_y$。任意实标签 $|c_y|\le1$ 都满足

$$
\left\|\sum_y c_yF_y-
{d\over2}\left\{Q,\sum_y c_yE_y\right\}\right\|\le\delta.
\tag{15}
$$

由(14)左右乘等距再取Hermitian部分即得；不把源插入与效果随意交换。对全部未知混态或与参考纠缠的输入，(14)—(15)保持联合未归一来源测试误差；如需低概率条件化，必须另除真实分支概率，不领免费相对精度。式(14)不是新增几何测量的Kraus后态实现定理。

## 5. 严格有理有限窗

S在R>0单调下降。代码以分数验证

$$
\sum_{j=0}^{49}{20^j\over j!}>10^6(1+20+20^2/3),
$$

故 $0<S(1000)<S(20)<10^{-6}$，不把浮点 $e^{-1000}$ 下溢当S=0。结合 $r_*<10/7$ 和
$\sqrt{1+S}/(1-S)^{3/2}<2$，

$$
\delta_0<{40021\over10500000}.
$$

再用 $\sqrt{84}<55/6$、$1/\sqrt3<26/45$、$L_B<22/7000$ 得

$$
\boxed{\delta<{4156812479\over38839500000}
\simeq0.107025386<27/250.}
\tag{16}
$$

只需此界，不调参追求更小误差。原ε界和全部平方根比较由代码的有理断言核验。

## 6. 同次真实指针—来源见证及受限恢复后果

取原全部六维空间内的 $|LL\rangle,|RR\rangle$ 双占据态，其Q为 $-1,+1$，D均为1。**理想**指针都是minus；实际两分布不被假定相等。

由(4)的等距向量差，理想plus分量为0，所以两个实际输入分别满足

$$
p_{+,i}=\|P_+J|i\rangle\|^2<\epsilon_0^2,
\qquad \operatorname{TV}(p_{LL},p_{RR})<\epsilon_0^2.
\tag{17}
$$

后一个界是两个bad概率处于同一区间，而非粗取两次ε。令 $B=P_-\otimes A$；这是原真实指针末读后的有界来源测试，$\|B\|\le1$。由(12)，每态A均值距相应 $\mp d$ 至多δ。删除plus分支改变A均值至多 $p_+$，因此

$$
\boxed{\langle B\rangle_{RR}-\langle B\rangle_{LL}
\ge2(d-\delta-\epsilon_0^2)>791/1000.}
\tag{18}
$$

可将 $F_\pm=(I\pm A)/2$ 与原指针联合，作为一份**数学任务效果**，得到联合分布q。F的存在由 $\|A\|\le1$ 保证，未称实际引力探测器已经实现。测试值 $1_{y=-}z$ 的绝对值至多1，故

$$
\operatorname{TV}(q_{RR},q_{LL})
\ge d-\delta-\epsilon_0^2>791/2000.
\tag{19}
$$

任意只收到**单次原二值指针**的共同合法CPTP后处理／随机恢复字典，其两个恢复输出的距离至多 $\epsilon_0^2$。三角不等式给至少一端的恢复误差

$$
\boxed{\max_i \operatorname{TV}(q_i,\mathcal R(p_i))
\ge{d-\delta-2\epsilon_0^2\over2}>0.1733.}
\tag{20}
$$

若只预报来源矩 $\langle B\rangle$ 且预报范围保留[-1,1]，同理至少一端的绝对误差大于0.3466。无需额外实现F也可将这理解为有界源矩字典的必要条件。

以上保守数值使用 $d>221/400,\delta<27/250$ 和原确切ε；若统一改用 $\epsilon_0^2<49/1000$，仍给(19)的791/2000和(20)的693/4000。恢复后果仅限定单次该记录、同一规则及合法有界任务；不排增加空间记录、保留原材料、使用准备身份、重复采样或无界统计权重。不据此判全部摘要失败，也不把该具体仪器提升为认知纲领的必要实现。

## 7. 数值校准与解析证明的分工

程序直接积分正常右Slater轨道下的f和f²，采用径向 $t=2r$ 的Gauss–Laguerre及轴向Gauss–Legendre规则。R1000计算保留 $0<S<10^{-6}$ 的独立解析界，压缩系数写为

$$
d_{\rm orbital}={a\langle\phi_R|f|\phi_R\rangle\over\sqrt{1-S^2}}.
\tag{21}
$$

反演对称使原交叉f矩为0；$G^{-1/2}ZG^{-1/2}=Z/\sqrt{1-S^2}$，故完整六维 $W^\dagger AW=d_{\rm orbital}Q$。这项真实源积分约为 $d_{\rm orbital}=0.5527805$，不同于把轨道直接当中心点；准确保存值见[results](results.json)。两阶积分比较仅是数值稳定性检查，不是统一严格积分余项；正式有限误差(16)不依赖该近似积分。

另在不过浮点范围的R20检查同一公式的Löwdin和全部六维CAR升维，只作矩阵身份诊断，不改任务参数。若 $f_1=W_1^\dagger fW_1,f_2=W_1^\dagger f^2W_1$，真实平方压缩为

$$
W^\dagger A^2W=(a/2)^2\{d\Gamma(f_1)^2+d\Gamma(f_2-f_1^2)\}.
\tag{22}
$$

程序保留第二项，校验 $\|(AW-WdQ)\|^2$ 的Gram矩阵，不能以压缩矩阵平方代替P外泄漏。R20不是新增仪器样例或寻找更好参数。

## 8. 本轮完成与未完成

完成的是同一1050真实父过程上的**全输入／参考来源字典及带真实指针的未归一来源效果**；有限差见证是后果。已有899的一般对偶测试、1050的仪器与能源界均复用，不重计为发现。

物理解释只到电子静质量领先Newton源系数。没有计算完整装置来源、有限G的 $O(G^2)$ 回反、真实引力仪器、有限速传播或完整GR约束与演化；无需为本有限合同继续修造这些部分。它不是预置Newton结构的生成证明，不结项M4/M5或整个ROADMAP。
