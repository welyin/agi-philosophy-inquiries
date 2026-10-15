# 1227：原 Gauss 热材料的实体字段记录、正能钟与共同能源熵预算

2026-10-15。作者稿；[代码](check.py)、[结果](results.json)。接[1226下一步](../1226/NEXT.md)。复用1128/1203的正能钟 Duhamel 引理及723—725的原材料域工具；新增的是原字段门交换子、同一热材料二阶矩、无限秩材料下的实际钟熵界，以及它们在同一未截断过程保住非零注能的联合证书。

## 1. 同一材料和新增硬件

原材料完全沿1226：$5\times5\times5$ 周期图，$N=125$；原 H⁵、商规范群、全部 CAR、质量、配对、电磁/梯度交叉和跳跃均保留于 $H$。几何在该轮正紧域内，$w\in[1/2,2]$，$\hbar=1$。参考为同一完整 Gauss Gibbs 态 $\sigma_\beta$，$\beta=17/10$。在两个事前节点 $B=(0,0,0)$、$D=(1,0,0)$，各执行一次598原字段 instrument：

$$
L_\pm(s)=\sqrt{\tfrac12\pm\tfrac14\sin s}.
\tag{1}
$$

附加两个指针 $P_1,P_2$ 和两个保留的环境 $E_1,E_2$，各为正常二级系统、空白 $|0\rangle$、正退化裸能源 $I$。简并不授予免费准备/复位。在B、D各放一份独立准备的正常钟 $C_j=L^2(\mathbb R)$，$[X_j,P_j^{C}]=i$、$H_{C_j}=(P_j^{C})^2/(2M)+1$。下文在单钟公式中把 $P_j^{C}$ 简记为 $P$，不与实体指针混同。程序坐标 $X_j$ 和 $M$ 是内部有效参数，未取得原三维坐标、SR质量或加速固有时资格。

在每对 $P_jE_j$ 的 $\{|00\rangle,|11\rangle\}$ 子空间定义 Pauli $Y_j$，其余两维上为零。取

$$
\theta(s)=\tfrac12\arccos(\sin s/2),\quad
G_j=\theta(s_j)Y_j,\quad \|G_j\|\le\pi/3.
\tag{2}
$$

则 $e^{-iG_j}|\psi\rangle|00\rangle=L_+|\psi\rangle|00\rangle+L_-|\psi\rangle|11\rangle$。这是同一 instrument 的实际有界平滑 Stinespring 耦合。全酉态保留相干和环境；经典分布来自查询 $P_1,P_2$。仅在迹掉环境的约化态中，才得到原按结果分块的 instrument 后态。

新增输入为四份空白、两正常钟的独立且已知时序对齐准备、固定位置门及各站本地字段耦合；没有由认知公理推导其制造过程。每门仅耦合该站原字段、该站钟和记录，不需要一份共享钟跨站访问；原Gauss图本身仍未获完整Lorentz因果认证。末读权限明示。几何标定沿1225同参数字典，但不要求同一态兼任标定包和热参考。

## 2. 同一材料的多项式二阶矩

记 $\mathcal L=H+1072N\ge0$、$A=\mathcal L+1\ge1$。钟能源不计入 $A$。1226给 $Z_{\mathcal L}(s)\le e^{Na_*(s)}$、
$Z_{\mathcal L}(\beta)\ge e^{-\beta e_*N}$，$e_*=5446$。因此

$$
\Pr_\beta(\mathcal L\ge x)\le
\min\{1,\exp[N(a_*(\beta/2)+\beta e_*)-\beta x/2]\}.
\tag{3}
$$

阈值 $x_0=N[2e_*+2a_*(\beta/2)/\beta]=Ne_\beta$。积分尾部并用1226的 $\langle\mathcal L\rangle\le x_0$：

$$
\langle\mathcal L^2\rangle_\beta\le x_0^2+4x_0/\beta+8/\beta^2,\qquad
\langle A^2\rangle_\beta\le x_0^2+(4/\beta+2)x_0+8/\beta^2+1.
\tag{4}
$$

这避免把指数配分函数比直接用作二阶矩。对 $s\in[.85,1.7]$，1226热迹有保守解析护栏 $a_*(s)<500$：$sc_T>1/5000$；标量热核 $<10^{10}$，势积分 $<10^5$；规范比较迹 $<10^{42}$。后者由
$\overline S_4(t/3)<3\cdot10^{14}$、
$\overline S_2(t/4)<4\cdot10^9$、
$2\overline S_0(t)-1<1000$ 得到；放大标签求和仅用于上界，物理群不变。
标量势积分核心、尾分别 $<5000,60750$，沿1226的 $e^{bR}=81$、$\Omega_4<27,b>1,\alpha\ge1/5,a=1/200$。
所以 $\log2<1,\log10<3$ 给 $a_*<32+45+3\cdot126=455<500$。

据此可以统一取

$$
e_\beta<11481,\quad x_0\le x=1435125,\quad
\sqrt{\langle A\rangle_\beta}<1200,\quad
r_0:=\sqrt{\langle A^2\rangle_\beta}<1200^2.
\tag{5}
$$

代码以有理数检查这些粗护栏，不计算完整热谱或实际热矩。

## 3. 原字段门交换子和材料矩运输

原曲目标直接给
$\Gamma_s=|\nabla s|^2=F(1-s^2/12)\le2$、
$\Delta s=-Fs/6$、$|\Delta s|<4/3$。而

$$
\theta'=-\frac{\cos s}{4(1-\sin^2s/4)^{1/2}},\qquad
\theta''=\frac{3\sin s}{16(1-\sin^2s/4)^{3/2}},\qquad
|\theta'|,|\theta''|\le1/3.
$$

原势、质量和CAR均与字段乘法对易。只由原节点动能得到

$$
[H,G_j]=-\frac1{2w}\big((\Delta\theta)Y_j+
2Y_j\nabla\theta\cdot\nabla\big),\qquad
\boxed{\|[H,G_j]\Psi\|\le C\|A^{1/2}\Psi\|,\ C=64.}
\tag{6}
$$

因为 $|\Delta\theta|\le10/9<2$、$|\nabla\theta|\le\sqrt2/3$、
$c_TL_v\le A$，一阶系数最多 $(2/3)\sqrt{8000}<60$。这适用于其它材料/钟/参考上的Hilbert值波函数。平滑有界门保形式域；共同核上(6)及闭包给 $G_jD(H)\subset D(H)$。没有删除无界相互作用。

对生成元 $H+\sum f_jG_j$，$f_j$ 可以是钟上有界算符，令
$r(t)=\|A\Psi_t\|$、$b(t)=\sum\|f_j\|$，则

$$
r'(t)\le Cb(t)\sqrt{r(t)},\qquad
\sqrt{r(t)}\le\sqrt{r_0}+\tfrac C2\int_0^t b(s)\,ds.
\tag{7}
$$

由(6)和 $\|A^{1/2}\Psi\|\le\sqrt{\|A\Psi\|}$ 得到；先作能源规整，再以闭性和一致界延拓。钟自由项与 $A$ 对易。第一矩类似给
$\sqrt{\langle A\rangle_t}\le\sqrt{\langle A\rangle_0}+(C/2)\int b$。可对无限秩 $\sigma_\beta$ 纯化后应用，不借总钟能源或有限秩假设。

## 4. 同一自治规则与两种比较器

取 $T=4$、门中心 $t_1=1,t_2=3$、$0<\delta<1/4$，

$$
b(u)=\tfrac{15}{16}(1-u^2)^2\mathbf1_{|u|\le1},\qquad
g_j(x)=\delta^{-1}b((x-t_j)/\delta).
$$

门为 $C^1$，面积1，范数 $<1/\delta$，导数范数 $<3/(2\delta^2)$。实际不显含时间的规则是

$$
H_{\rm all}=H+\sum_{j=1}^2(P_j^{C})^2/(2M)+6+V,\qquad
V=g_1(X_1)G_1+g_2(X_2)G_2,\qquad
\|\partial_{X_j}V\|\le L_j:=2/\delta^2,\quad L_1+L_2=L:=4/\delta^2.
\tag{8}
$$

6含两钟各1及四个记录系统各1。$V$有界，故在自由总算符域上自伴、半有界，保Gauss/总宇称；未截断任何原自由度。两钟各使用同规格初态

$$
\varphi(x)=(2\pi\sigma^2)^{-1/4}e^{-x^2/(4\sigma^2)+iMx}
\tag{9}
$$

两钟相互独立并独立于材料和空白，均速1，没有到达/出射后选。第一站程序门在其坐标1，第二站在其坐标3；共同初态相位/位置对齐是所声明准备的一部分。所有作用保留到共同末时 $T$，不硬删 Gaussian 尾。

比较器一是自由钟及经典有限脉冲 $H+\sum g_j(t)G_j$；比较器二是原瞬时instrument和真实等待：

$$
U_{\rm id}=e^{-iH}e^{-iG_2}e^{-2iH}e^{-iG_1}e^{-iH}.
\tag{10}
$$

所有自由相位保留。热态平稳仅在相应边缘中使用，不删除纯化全态中的相位。

## 5. 分开概率与无界能源的两次比较

**两实际钟到有限脉冲。** 直接复用1128/1203。各自由钟宽
$\sigma_t=\sqrt{\sigma^2+t^2/(4M^2\sigma^2)}$，$b_T=\sup\sigma_t$；完整纯化振幅误差

$$
\epsilon_C\le(L_1+L_2)Tb_T=LTb_T.
\tag{11}
$$

实际传播位于Duhamel左侧，保钟反作用。实际算符系数仅取粗界
$\int b\le2T/\delta$，因此
$r_{\rm a}\le(1200+CT/\delta)^2$；不能把它擅自换成面积2。
经典脉冲两门面积为2，理想两门亦然，给
$r_{\rm f},r_{\rm id}\le R:=(1200+C)^2$。
两纯化的能源差满足
$|\langle A\rangle_\Psi-\langle A\rangle_\Phi|
\le\|\Psi-\Phi\|(\|A\Psi\|+\|A\Phi\|)$，
故此步误差至多 $(r_{\rm a}+R)\epsilon_C$。

**有限脉冲到原瞬时门。** 在原 $H$ 相互作用绘景中，
$G_j(t)-G_j(t_j)=\int_{t_j}^t i e^{iHs}[H,G_j]e^{-iHs}ds$。
窗口内 $|t-t_j|\le\delta$，两门第一矩运输保证相应混合路径
$\|A^{1/2}\Psi\|\le1200+2C$。用(6)、每门面积1及逐窗口Duhamel，

$$
\epsilon_{\rm p}\le D\delta,\quad D=2C(1200+2C)=169984,\qquad
\boxed{\epsilon_E\le(r_{\rm a}+R)\epsilon_C+2RD\delta.}
\tag{12}
$$

第二步两端矩都是 $R$，不把它的 $\delta$ 误差乘实际钟粗矩 $\delta^{-2}$。原 $H$ 在有限门内持续运行，不能免费删去。

全纯化半迹距、两指针联合分布总变差和任意保留后态询问均受
$\epsilon_C+\epsilon_{\rm p}$ 控制。若报告 classical-quantum instrument 后态，两比较态都先约化环境/查询指针；不把完整酉相干态直接比作已去相干态。无界材料能源另由(12)控制。

上述具体数值预算以材料边缘固定为原 $\sigma_\beta$ 为前提，可保留它的任意参考扩张；不声称对任意初始 Gauss 态有同一热矩和同一误差预算。

## 6. 原正注能的可计算严格留量

598非选择记录的完整能源变化为

$$
\Delta E_j=\langle J_j\rangle,\qquad
J_j=\Gamma_s\theta'(s)^2/(2w)\ge0,\qquad J_j\le2/9<1/4.
\tag{13}
$$

其它完整H项对易但未删除。两门总理想注能 $\le1/2$。第一门之前是平稳Gibbs，第二门对其真实理想后态仍非负，故第一门下界足够。

取每点半径 $a=1/4$ 的正常包 $(1-r^2/a^2)^4_+$，Haar常链路、CAR真空，构成同一Gauss试态 $f_a$。按1226：动能 $240/a^2=3840$，势2、梯度36、磁4096、质量/跳跃期望0，加1072移位，总
$\langle\mathcal L\rangle_{f_a}\le9046N$。仅要求形式域，不称其为H不变空间。

其支撑上 $s^2\le|\phi|^2\le1/8$，
$F\ge95/48$、$\Gamma_s\ge(95/48)(95/96)>1$、
$\cos s\ge15/16$；因此 $J_1\ge225/16384>1/128$。
以该支撑投影再支配归一秩一投影，有
$J_1\ge|f_a\rangle\langle f_a|/128$。
Jensen及 $a_*(\beta)<500$ 给

$$
\boxed{\Delta E_{\rm id}\ge\operatorname{Tr}\sigma_\beta J_1
\ge j_0:=2^{-7}e^{-1984775}>0,}
\tag{14}
$$

因为 $125[(17/10)9046+500]=1984775$。这是极保守的显式对数下界，不是实用测量精度。

## 7. 一份有限参数证书和实际钟损能

固定一份充分实例，不作参数优化：

$$
\begin{gathered}
R=1597696,\ D=169984,\ T=4,\quad
\delta=j_0/(32RD),\quad L=4/\delta^2,\quad
r_{\rm a}=(1200+256/\delta)^2,\\
\sigma=j_0/[128LT(r_{\rm a}+R+1)],\qquad
M=T/\sigma^2+16L^2T^2/j_0+1.
\end{gathered}
\tag{15}
$$

全部严格正且有限。$M\ge T/\sigma^2$ 给 $b_T<2\sigma$，故

$$
\epsilon_C<\frac{j_0}{64(r_{\rm a}+R+1)},\quad
\epsilon_{\rm p}\le\frac{j_0}{32R},\quad
\epsilon_E<5j_0/64,\quad
\epsilon_C+\epsilon_{\rm p}<j_0/(16R).
\tag{16}
$$

每钟有界实际力给 $\|P_j^{C}(t)-P_j^{C}(0)\|\le L_jt$、
$\|(X_j(t)-t)\Psi_{\rm in}\|\le\sigma_t+L_jt^2/(2M)$。
初末中心均在相应门外，且 $\sum_jL_j^2\le L^2$，所以

$$
\eta:=|\langle V\rangle_T-\langle V\rangle_0|
\le L(3\sigma+LT^2/(2M))<j_0/16.
\tag{17}
$$

初始 $Y$ 期望恰为零，但仍用此界覆盖两端，不删末端尾。初态总能源有限、自伴总H守恒：

$$
\Delta E_H+\Delta E_{C_1C_2}+\Delta E_{P,E}+\Delta\langle V\rangle=0,\qquad
\Delta E_{P,E}=0.
$$

于是

$$
\boxed{\Delta E_H^{\rm actual}>59j_0/64>0,\qquad
\Delta E_{C_1C_2}<-55j_0/64<0.}
\tag{18}
$$

另有 $\Delta E_H<1/2+5j_0/64$、
$\Delta E_{C_1C_2}>-1/2-9j_0/64$。这是两钟合计的损能方向，不分别宣称每钟都满足同一下界。没有以巨大钟初末能源的浮点相减声称算出微小功；用的是同一动力学、独立材料误差和完整守恒。钟损能和材料注能是同过程预测，未制造分辨此极小均值的额外能源仪表。

每钟初始平均能源为 $M/2+1/(8M\sigma^2)+1$，两钟总准备能源为其两倍；四指针/环境基线4，材料能源/熵由原热态给出；$\|V\|\le2\pi/(3\delta)$、力和瞬时材料矩均登记。资源量只以对数展示，避免机器下溢/溢出。此证书有限但极昂贵，不是工程规格或免费操作。

## 8. 无限秩热材料下的实际钟熵

当前热材料无限秩，不能复制1203的有限rank论证。对每个实际末钟使用反自由传播参考模（下式省略钟标签j）

$$
a_T=(X-TP/M)/(2\sigma)+i\sigma(P-M),\qquad [a_T,a_T^\dagger]=1.
$$

这是计算熵的正规参考模，不从能源账删除大基线。初始钟被
$a_0=X/(2\sigma)+i\sigma(P-M)$ 湮灭。实际Heisenberg运动中

$$
\|P(T)-P(0)\|\le LT,\qquad
\|X(T)-TP(T)/M-X(0)\|\le LT^2/(2M).
$$

第二项等于 $-M^{-1}\int_0^T s\,\dot P(s)ds$ 的范数界。对热材料任意纯化仍有初始 $a_0\Psi=0$，从而

$$
n_C:=\operatorname{Tr}\rho_C(T)a_T^\dagger a_T
\le(\sigma LT+LT^2/(4M\sigma))^2
\le j_0^2/4096=:\bar n<1.
\tag{19}
$$

由(15)第二项 $\le\sigma LT/4$，两项和 $\le5j_0/[512(r_{\rm a}+R+1)]$，足够推出该宽界。固定单模占有下的最大熵变分不要求后态Gaussian：

$$
\boxed{S(C_T)\le g(\bar n)
=(1+\bar n)\log(1+\bar n)-\bar n\log\bar n
\le\bar n(2-\log\bar n)<\infty.}
\tag{20}
$$

每钟可用更松的总力上界 $L\ge L_j$，故上式分别成立。末态两钟可能经材料产生关联，不假设其独立；熵次可加给

$$
S(C_1C_2)_T\le S(C_1)_T+S(C_2)_T
\le2g(\bar n)\le2\bar n(2-\log\bar n).
$$

这覆盖无限秩材料导致的钟混合，无有限维Fannes或rank假设。无界域由bounded-force的Heisenberg差和正常初态二阶矩延拓。指针/环境熵 $\le4\log2$；材料末态因有限 $\langle A\rangle$ 与原热迹而熵有限。全部记录仍保留，未授予免费复位。

作为旧1203/Gibbs热力学身份在本实际过程的复用，令 $\mathcal R=P_1E_1P_2E_2C_1C_2$。初始 $\mathcal R$ 纯且与热材料独立，完整幺正演化给

$$
\Delta S_{\mathcal R}+\beta\Delta E_H
=I(\mathcal R:\mathrm{material})_{\rm out}
+D(\rho_H'\Vert\sigma_\beta)\ge0.
\tag{21}
$$

由本轮有限 $S(C_1C_2)$、$S(PE)\le4\log2$ 和材料有限能源/热迹，式中没有 $\infty-\infty$。这不是新一组发现；它说明旧四律接口现在可用于这份同对象过程。$S_{\mathcal R}$ 是保留全部装置的实际量子熵，不能用末端指针的 Shannon 熵替代；初末交互能源仍按(17)计入能源总账。

## 9. 复算与停止线

四组检查：原字段门/交换子和热矩有理护栏；Stinespring原概率及注能函数；固定参数、分开的概率/能源/反作用留量；无限秩钟熵及资源对数预算。小指针矩阵只核读口身份，不替代完整Gauss或连续钟，也未运行截断散射模拟。

精确不等式、积分、有理护栏承担证明；浮点仅展示对数和恒等式诊断，不称区间外向或实际系统误差认证。删除 $\theta$ 的字段依赖同时删除原读取和正注能；跳过有限门期间H演化漏掉(12)；仅用迹距无(6)—(7)不能签材料能源；沿有限rank论证不能签(20)。

本轮为原热材料、内部触发、两实体记录及能源/钟熵提供一份共同有限充分证书。原图、正紧几何、热参考和全部准备/耦合/末读权限仍是登记输入。没有推导boost/Unruh、真实空间运动或引力反作用，没有证明全部硬件由认知公理产生，不签阶段完成。本实例之后停止钟参数优化。

