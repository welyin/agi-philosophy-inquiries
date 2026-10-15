# 1226：原完整 Gauss 热参考的规模一致边界面积上界

2026-10-15。作者证明；[代码](check.py)、[结果](results.json)。复用[1221](../1221/proof.md)的匹配比较态与双相对熵身份、[1225](../1225/proof.md)的独立面积字典。本轮新增的是**原无界材料的局部热资源由同一完整 Gibbs 态推出，且常数不随图的节点数增长**。不另设每边能源预算，不计算或替换完整 Gibbs 谱。

## 1. 定理域和继承输入

取原周期立方图 $V_L=(\mathbb Z/L\mathbb Z)^3$，$L\ge5$、$N=L^3$；每点三条正向边与三张正向面。保留589的全部 H⁵ 标量、同源非对角电项、测地梯度、完整磁星及598的全部 CAR 质量、配对和规范协变跳跃；群仍为

$$
\mathsf G=(SU(3)\times SU(2)\times U(1)_Q)/\mathbb Z_6.
$$

固定格距、$\hbar=1$、原 $M=2$；均匀几何在正紧集内，$\det G=1$。有限代数和耦合固定且平移不变。显式常数例采用

$$
\tfrac12\le w\le2,\quad \tfrac14I\le G\le4I,\quad f=32,
\quad \|t_e\|\le .17.
\tag{1}
$$

这里每边 $t_e$ 是原单粒子跳跃矩阵，允许保持全部原模块；.17 沿644的中性跳跃幅度并作为其余原允许跳跃的上限。若其它跳跃超过该范围，须用实际范数重算常数，不得静默删除。质量系数采用598已登记的一代复 $Y$；束缚势用574/598继承的同一 $L,u$，磁/电系数用569。数据与来源在代码中绑定哈希。任意另一固定有限代/耦合紧集，只要保同样正束缚和动能下界，以下符号证明成立，但不自动授予本例数值。

参考为**完整有限体积 Gauss 空间**上的 $\sigma_{\beta,N}=e^{-\beta H_N}/Z_N$，$\beta>0$，不额外挑选排除下述 Gauss 试态的全局扇区。它不是小矩阵压缩 Gibbs 态。均匀背景、原图/材料、固定温度选择都是旧候选的范围输入，不是新认知公理或热化装置。认知动机是同一热参考、来源、区域和独立面积的相容连接。

## 2. 每节点质量界与保留完整电动能

598给全目标界 $F^{-2}\le AU+B_0$。原表重建

$$
\lambda_0=\det L/\operatorname{tr}L>.09,\quad
D=12-u_h-u_s>11.25,\quad
A=\frac{32}{\lambda_0D^2}<3,\quad B_0=\frac{144}{D^2}<2.
\tag{2}
$$

这些护栏可按冻结十进制参数作有理数复算；不使用全图的 $C+cW^{1/4}$ 常数。原复系数满足
$|Y_u|<5/4,|Y_d|<3/5,|Y_e|<1/4,|Y_\nu|<1/2,|Y_s|<1/3$。逐个 CAR 单项范数至多1，双重态两分量的绝对值和不超过 $\sqrt{24}$，$|s|<\sqrt{12}$，故单点完整质量满足

$$
\|B_v\|\le M_Y F_v^{-1/2},\quad
M_Y=2\sqrt{24}\{6(|Y_u|+|Y_d|)+2(|Y_e|+|Y_\nu|)\}
       +2\sqrt{12}|Y_s|<128.
$$

色/旋量重数分别为6和2；最后一项含 Majorana 的产生与湮灭，未漏 $1/2$ 反对称求和。$3^{1/4}<4/3,2^{1/4}<6/5$ 给 $\|B_v\|\le171U_v^{1/4}+154$。Young 的逐点最大值为

$$
\sup_{x\ge0}(a x^{1/4}-\delta x)
 =\frac{3a^{4/3}}{4^{4/3}\delta^{1/3}}.
$$

取 $\delta=1/4$，$171^{1/3}<6$，得

$$
\boxed{\|B_v\|\le U_v/4+1024.}
\tag{3}
$$

每条原跳跃的范数 $\le2f\|t_e\|<16$，共有 $3N$ 条边。这里以单粒子奇异值分解和 CAR 范数估计，未把 $R(U)$ 的每个矩阵元素另乘一个 $f$。

原电项在均匀几何下为 $w^{-1/3}\sum_{v,c,A}b_c G_{ij}P_{vi,cA}P_{vj,cA}$。569原 $b_c\in(3/2000,7/10)$；由 $2^{1/3}<4/3$ 和(1)，其参考 Casimir 下界大于 $9/32000>1/4000$。标量动能系数 $1/(2w)\ge1/4$。因此，作为闭形式，

$$
\boxed{\mathcal L_N:=H_N+1072N\ \ge\ K_N
 =c_T\left(\sum_v L_v+\sum_e L_e\right)+u\sum_vU_v\ge0,\quad
c_T=1/4000,\quad u=1/4.}
\tag{4}
$$

$L_v=-\Delta_{H^5}$，$L_e$ 为637原规范群参考 Casimir。$1072=1024+3\cdot16$。非对角电项由正系数矩阵整体控制；梯度和磁势保持非负并仍在 $H_N$ 中。移位只是谱估计，未从物理能源/来源扣掉 $1072N$。

对形式域的控制复用598/1221。以下所有无界期望先在共同核上估计，再由闭形式与非负谱截断延拓；不对非对易算符使用指数单调性。

## 3. 显式正常 Gauss 试态与每点下界

每点取同一径向包 $f(r)\propto(1-r^2)^4_+$，$r=d_{H^5}(\phi,0)$；链路为归一 Haar 常函数，全 CAR 取空真空。此向量正常、Gauss 不变、总偶且在形式域内；包为 $C^3$，本轮不称它为 $C^\infty$ 或 $H^3$ 矩证书。它无需在完整 $H_N$ 下不变。

半径1内双曲径向 Jacobian 比 Euclidean 的
$J(r)=[\sinh(r/\sqrt6)/(r/\sqrt6)]^4<e^2<8$。分母用 $J\ge1$，分子用 $J<8$，Beta 积分给

$$
\frac{\int|f'|^2d\mu}{\int|f|^2d\mu}
\le64\cdot8\frac{B(7/2,7)}{B(5/2,9)}=240.
$$

于是标量动能每点至多240。原 $\|L\|<1,u_h+u_s<1$；$r\le1$ 时 $|\phi|^2\le2,F\ge5/3$，从而 $U<81/100<1$，$wU\le2$。所有标量目标间距离不超过2，原梯度星至多36。电项在 Haar 常函数上为零；正常排序的质量、Majorana 和跳跃在 CAR 真空中期望为零，虽其作用一般离开真空。

为保留全部磁交叉，记 $S_{\alpha d}=\sum_R\alpha_Rd_R$。由569原系数

$$
S_{\alpha d}\le\delta(12w_q+4w_l)+8r_c+3r_w+r_Y<16.
$$

右侧允许多计零作用的平凡表示。$B=\Lambda^2G^{-1}$ 谱上界 $g_-^{-2}$，三面的 $\|O_{R,f}\|_{HS}^2\le d_R$、$\|E_{R,f}\|_{HS}^2\le4d_R$，故**完整三面磁星**

$$
0\le W_{{\rm mag},v}\le
\frac{15S_{\alpha d}}{2w_-^{1/3}g_-^2}<2560<4096.
\tag{5}
$$

没有对各交叉项分别假设正性；使用整个正二次型的谱界。因此

$$
\langle H_N\rangle_{\rm trial}\le4374N,
\qquad \boxed{\langle\mathcal L_N\rangle_{\rm trial}\le e_*N,\quad e_*=5446.}
\tag{6}
$$

## 4. 可计算局部热迹上界

603/637已经给原无限维热迹存在。本轮把其因子按点数显式组织。由(4)的 min–max 比较，Gauss 子空间迹不超过全空间乘积迹：

$$
Z_{\mathcal L_N}(s)\le
\big[2^{32}z_\phi(s)z_{\mathsf G}(sc_T)^3\big]^N.
\tag{7}
$$

局部标量迹有 $z_\phi(s)\le k_{H^5}(sc_T)I(su)$，其中

$$
k_{H^5}(t)=(4\pi t)^{-5/2}e^{-2t/3}(1+t/9).
$$

603原尾界可保守取 $a=1/200,b=4/\sqrt6,R_0=\sqrt6\log3$，给

$$
I(\alpha)\le\operatorname{Vol}(B_{R_0})+
\frac{9\Omega_4}{4b\alpha a}\exp[-\alpha a e^{bR_0}],\quad
\Omega_4=8\pi^2/3.
\tag{8}
$$

这里保留 $b=4/\sqrt6$ 与双曲体积增长的匹配，不以较小指数代入旧尾积分公式。原 $a_{603}=\lambda_0D^2/(512M^2)>1/200$；$u_h+u_s<1$ 与 $\tanh^2(\log3)=16/25>13/24$ 保证所选 $R_0$ 大于603的尾部起点，不靠浮点反双曲函数来认证。此时 $e^{bR_0}=81$。球体积的闭式，$x=R_0/\sqrt6$，为
$\Omega_4(36\sqrt6/8)[\sinh(4x)/4-2\sinh(2x)+3x]$。

规范热迹仍是原商群的 $z_{\mathsf G}(t)=\sum_{\lambda\ {\mathrm{allowed}}}d_\lambda^2e^{-tC_\lambda}$。只为给显式**上界**可放大到覆盖群标签和：若

$$
S_k(\alpha):=\sum_{n\ge0}(n+1)^k e^{-\alpha n^2},\quad
\overline S_k(\alpha)=1+\sum_{j=0}^k {k\choose j}2^{k-j}
\frac{\Gamma((j+1)/2)}{2\alpha^{(j+1)/2}},
$$

则 $S_k\le\overline S_k$。证明在每个 $[n-1,n]$，$(x+2)^ke^{-\alpha x^2}\ge(n+1)^ke^{-\alpha n^2}$。SU(3) 维数平方 $\le(a+1)^4(b+1)^4$，其 Casimir $\ge(a^2+b^2)/3$；SU(2) 的维数平方 $(\ell+1)^2$，Casimir $\ge\ell^2/4$。故

$$
z_{\mathsf G}(t)\le\overline z_{\mathsf G}(t):=
\overline S_4(t/3)^2\overline S_2(t/4)
\big(2\overline S_0(t)-1\big).
\tag{9}
$$

用于物理态和标签的仍是 $\mathbb Z_6$ 商；(9)不是改变规范群。所有上界可直接标量求值，不需要全图谱。设

$$
a_*(s)=\max\{0,32\log2+\log[k_{H^5}(sc_T)I(su)]
                  +3\log\overline z_{\mathsf G}(sc_T)\}.
\tag{10}
$$

## 5. 同一完整热态给出规模一致的单点资源

谱测度 Jensen 和(6)给
$\log Z_{\mathcal L_N}(\beta)\ge-\beta e_*N$。令 $F_N=\log Z_{\mathcal L_N}$；603确保可微、凸。凸性在 $\beta/2,\beta$ 给

$$
\begin{aligned}
\langle K_N\rangle_\beta&\le\langle\mathcal L_N\rangle_\beta=-F_N'(\beta)\\
&\le\frac2\beta\big[F_N(\beta/2)-F_N(\beta)\big]
\le N e_\beta,\qquad
\boxed{e_\beta=2e_*+\frac{2a_*(\beta/2)}\beta.}
\end{aligned}
\tag{11}
$$

原平移由场/链路置换和 CAR 第二量子化置换实现，保持完整 $H_N$ 与 Gauss 空间。因此有限体积完整 Gibbs 态严格平移不变。把 $K_N=\sum_vk_v$，每点归其三条出边，则

$$
\boxed{\langle U_v\rangle_\beta\le e_\beta/u,
\qquad \sum_{i=1}^3\langle L_{(v,i)}\rangle_\beta\le e_\beta/c_T.}
\tag{12}
$$

不要求三个方向实际均值相等，也不选破坏平移的单一热力学极限态。对紧正 $\beta$ 区间取上式常数的最大值即可；常数不依赖 $L,N$，但可依赖格距、几何域、温度和原耦合。任意非平衡有限总能态不自动满足(12)。

## 6. 同一比较态、原无界边项与簇计数

取任意节点/格胞区域 $A\subset V_L$，切边数 $b_A$。完全沿1221的 $J$ 和 $\tau=\tau[\sigma_{\beta,N}]$。$\tau$保区域边缘；切边 Casimir 是同一表示标签的函数，故其期望保持；每点 $U_v$ 也保持。(12)因此对 $\tau$逐点成立，尽管 $\tau$本身一般不平移不变。

其有限完整形式能源先按1221已证固定图反向比较取得，因而以下差值均合法；不靠本轮面积结论循环假设定义域。原节点项留在所属区域。**完整电星留在源节点侧**：617切边运输把 $P_{vi}$ 作用到源半边，同源 $P_{vi}P_{vj}$ 仍是该侧区域偶算符。因此其期望差为零；保留全部非对角电项，未将两侧授予独立完整 Hamiltonian。这里不再同时使用另一套“电跨簇”计数。

完整正梯度星若支持四点跨界，就整星列入 $H_\partial$，未跨界的留在区域内。原群固定目标原点，故 $d(\phi_v,R(U)\phi_w)^2\le2(r_v^2+r_w^2)$。由 $F=2\operatorname{sech}^2(r/\sqrt6)$ 和 $x^2\le\sinh^2x\le\cosh^4x$，

$$
r_v^2\le24F_v^{-2}\le24(3U_v+2)\le72(1+U_v).
\tag{13}
$$

原完整梯度星（含交叉）于是有

$$
0\le W_{{\rm grad},v}\le\frac{w_+^{1/3}}{g_-}
          \left(3r_v^2+\sum_i r_{v+e_i}^2\right),\quad
\langle W_{{\rm grad},v}\rangle_{\sigma,\tau}
\le2304(1+4e_\beta).
\tag{14}
$$

计数使用**完整簇支持**，不是逐个抹去交叉：

|跨界项|连通支撑|每条切边可归属的簇数|单簇估计|
|---|---|---:|---|
|完整梯度出星|$v,v+e_1,v+e_2,v+e_3$，3条边|1|每态均值见(14)|
|完整三面磁星|7个点、9条边，保全部面间 $O$ 交叉和 $E$ 余项|3|范数 $\le4096$|
|原 CAR 跳跃|一条边两个端点|1|范数 $\le16$|

跨星的连通支撑必含切边。每条梯度边只有其源基点；固定正向 $i$ 边属于磁星基点 $v$、$v-e_j$（两个 $j\ne i$），所以磁星数 $\le3b_A$。周期 $L\ge5$ 避免本地顶点重合。该关联证明对所有 $N,A$ 成立；有限图程序检查只是实现核对。

两态期望差按两个绝对值保守相加，得

$$
\boxed{0\le\Delta E_\partial:=\operatorname{Tr}(\tau-\sigma)H_\partial
\le C_\partial(\beta)b_A,\quad
C_\partial=18432e_\beta+29216.}
\tag{15}
$$

其中 $C_\partial=2\cdot2304(1+4e_\beta)+6\cdot4096+2\cdot16$；非负性来自1221双相对熵身份，非来自假设每项差非负。此为同一完整热参考导出的资源结论，不再输入每边预算。

## 7. 独立面积与双相对熵上界

1225定义并给原字段有限菜单的条件性稳定字典：$h=w^{2/3}G$、$P=\operatorname{cof}(h)$，坐标法向面面积 $a_i=\sqrt{P_{ii}}$。固定该格胞解释，区域边界独立面积为

$$
A_\partial=\sum_{e\in\partial A}a_{i(e)},\qquad
a_i\ge w_-^{2/3}/\sqrt{g_+}>3/10=:a_-.
\tag{16}
$$

没有由熵定义面积。1221精确身份与(15)联合给

$$
\boxed{D(\sigma\Vert\tau)+D(\tau\Vert\sigma)
=\beta\Delta E_\partial
\le\frac{\beta C_\partial(\beta)}{a_-}A_\partial.}
\tag{17}
$$

在 $5\times5\times5$ 图（125节点）上可直接复用1225已证明的字典记录接口；任意更大 $L$ 的本轮面积不需要额外完整读口认证才具有几何定义。若把有限菜单实际操作推广到其它图，须保持其局部支持隔离及原有限时/高阶矩条件，不能从(17)自动获得单次面积读数。若一次合法面积推断有 $|\widehat A-A_\partial|\le\delta_A$，上界用 $\widehat A+\delta_A$；统计失败概率保持原登记，不称无误差认证。

## 8. 完整互信息：保留中心与强制配对

1221身份为

$$
I(A:B)_{J\sigma J^\dagger}
=H(p_{\rm joint})+2\sum_{e\in\partial A}\mathbb E\log d_{\lambda_e}
  +D(\sigma\Vert\tau).
\tag{18}
$$

每条边只用原商群允许标签；$C_\lambda$ 是同一个 $L_e$ 的本征值。取比较概率
$q_\lambda(t)=d_\lambda^2 e^{-tC_\lambda}/z_{\mathsf G}(t)$，$t>0$。KL非负给

$$
H(p_e)+2\mathbb E\log d_{\lambda_e}
\le t\mathbb E C_{\lambda_e}+\log z_{\mathsf G}(t)
\le t e_\beta/c_T+\log\overline z_{\mathsf G}(t).
\tag{19}
$$

637已给这些完整有限能态的熵有限；也可按有限标签截断与相对熵变分不等式延拓，避免 $\infty-\infty$。用 $H(p_{\rm joint})\le\sum_eH(p_e)$，无需标签独立，得到

$$
\boxed{I(A:B)_{J\sigma J^\dagger}\le
\frac{\beta C_\partial(\beta)+t e_\beta/c_T+
\log\overline z_{\mathsf G}(t)}{a_-}\,A_\partial,\quad t>0.}
\tag{20}
$$

这是保原中心、配对、跨 sector 相干和全部材料的完整扩展互信息上界；不是把其前两项免费抹除后的答案。$\tau$仍是非线性变分比较态，不是已实现的去关联过程或免费功来源。

## 9. 复算内容、误差与停止线

代码分四组：原参数和逐点有理护栏；局部热迹/温度常数；原完整簇的局部关联与多图区域计数；原商群标签的中心/配对 KL 及独立面积误差运输。原结果仅读取并绑定，未重跑旧全谱/矩阵诊断。输出的温度常数用浮点求值，标为解析公式的数值展示，不声称区间外向认证或实际 Gibbs 能源/熵的数值计算；严格结论由上述精确公式和保守有理护栏承担。旧冻结参数本身是所选模型输入，不从数值 RG 误差推断新实验精度。

本轮新增的是原完整无界模型的 $N$ 一致局部热资源和(17)、(20)。仍限定固定格距、正紧几何、指定完整平移 Gibbs 参考和固定材料耦合。没有证明普通热区域熵 $S(A)$ 的面积律、$S=sA$、$\delta S=s\delta A$、UV一致系数、Newton常数、内部热化/温度装置、加速共动钟、引力动力学或阶段结项。$\beta\to\infty$、格距/体积趋零、非均匀或非平衡高能集中态须另审，不自动落入当前结论。
