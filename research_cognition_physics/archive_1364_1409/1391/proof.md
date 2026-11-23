# 一个重下夸克上的两个真空角

## 微扰轨道与 CP

设在紫外有实 CP 基，裸强角为零，真空中除可用规范变换消去的角以外，仅有一个带非零 PQ 荷的标量相位，其他非零真空量 PQ 中性且为实数。精确的微扰 PQ 对称把该相位整圈扫过。CP 共轭真空与原真空在同一轨道，其中存在全部参数和真空量均为实数的代表。任意弱 CP 奇的基不变量既在此轨道上恒定，又须变号，因而等于零。这也覆盖保持 PQ 的高次质量算符。QCD 反常允许共同角改变强角。弱作用具有另一独立显式 CP 源，或者多个已固定的带荷真空量共同提供相对相位时，前件改变。

约束内最少须留下一个不随 PQ 平移的物理真空角。以下两标量、一个重载体给出充分实例。矢量样夸克与标量自发传递 CP 及轴子机制已有具体模型，此处检验同一个质量矩阵上的必要角数、阈值匹配与味流。[Peccei–Quinn](https://doi.org/10.1103/PhysRevLett.38.1440)；[Hashimoto 等](https://arxiv.org/abs/2105.06089)；[Suematsu](https://arxiv.org/abs/2309.04783)。

## 同一物质实例

加入一个与右手下夸克具有相同规范表示的矢量样夸克 $D$。取 $X(D_L)=X(S)=X(T)=1$，$X(D_R)=0$，通常夸克与唯一电弱双态 $H$ 的 PQ 荷为零。实 Yukawa 生成上三行 $A=(m,c)$ 与下行 $(S b_S+T b_T,4S)$。以任意统一质量单位写出见证。上型质量实对角、正且互异，三种上味由各自质量及带电流命名。取 $S=e^{i\phi}$，$T=e^{i(\phi+\delta)}$，$t>0$ 同时放大两个 PQ 真空量。

$$
m=\begin{pmatrix}1&0&0\\1&2&0\\1&1&3\end{pmatrix},
\quad c=(0,1,1/2)^T,\quad b_S=(1,0,1),\quad b_T=(0,1,1),
\qquad
{\cal M}=\begin{pmatrix} A\\ t e^{i\phi} r\end{pmatrix},
\quad r=(b_S+e^{i\delta}b_T,4).
$$

通过 Schur 补直接得到 $\det{\cal M}=6t e^{i\phi}(4-e^{i\delta}/2)$，处处可逆。实上右列允许存在；PQ 禁止重夸克的裸质量，普通电弱质量保持中性。令 $H_d={\cal M}{\cal M}^\dagger$，$P_1,P_2$ 由上型质量投影和带电流拉回下型双态，$P_D$ 为剩余弱单态投影。夸克弱基不变量为

$$
I=\Im\operatorname{tr}(P_1 H_d P_2 H_d P_D H_d)
=\Im[(H_d)_{12}(H_d)_{24}(H_d)_{41}]=-2t^2\sin\delta .
$$

因此即使考虑标量重标记与费米子换基构成的广义 CP，只要 $\sin\delta\ne0$，原真空就不能由保持所有相互作用的广义 CP 固定。轻态四元组取 $J=\Im(V_{11}V_{22}V_{12}^*V_{21}^*)$，轻下态按质量升序，$V$ 为实际 $3\times4$ 带电流。对 $\cos\delta=1/3,\sin\delta>0$，$t=1$ 得 $J=-0.00298844$，真空共轭则符号相反。共同角不改变 $J$。

当 $t\to\infty$，轻左态的平方质量算符趋于 $H_*=AA^\dagger-(Ar^\dagger)(Ar^\dagger)^\dagger/(rr^\dagger)$，由 $H_d-\lambda I$ 的重块 Schur 补即得。这时 $\Im[(H_*)_{12}(H_*)_{23}(H_*)_{31}]=-63\sqrt2/961$，故轻谱非简并且弱 CP 存活。极限 $J=-0.0008514964$，三个轻态的总弱单态份额在 $t=1,10,100$ 时分别为 $0.3250,0.002463,0.00002451$。这些数展示相位供给与重态解耦的差别，质量单位和矩阵尚未用现实质量、CKM 或强 CP 的数据选定。

## 角势、阈值与同一极小

对任意正 $\lambda_S,\lambda_T,\kappa$，一个有界 PQ 与 CP 不变势为

$$
V_0=\lambda_S(|S|^2-v_S^2)^2+\lambda_T(|T|^2-v_T^2)^2+
\kappa\,[\Re(S^\dagger T)-v_S v_T/3]^2 .
$$

在其全部极小，模长分别为 $v_S,v_T>0$，相对角 $\delta=\pm\arccos(1/3)$。共同角自由。取 $v_S=v_T=t$ 返回上述矩阵。一个基本色三重态的手征荷差为一，反常 $N=2T(\mathbf3)=1$；以同一归一和普通电荷 $-1/3$ 计 $E=2N_c(1/3)^2=2/3$。以 Dirac 约定 $\bar\theta=\theta+\arg\det(M_u{\cal M})$，跨重阈值要连同被积掉的质量相位一起匹配。取正质量向量样纯 QCD 的真空能在 $\bar\theta=0\bmod 2\pi$ 达到极小、该处 $\chi>0$。对每个 $\delta$ 可独立取

$$
\phi_*=-\arg(4-e^{i\delta}/2)\pmod{2\pi}.
$$

见证正支 $\phi_*=0.122360745$。这是按 QCD 势求得的真空位置。在不含其它重要显式 PQ 破缺势的有效描述中，任意保 PQ 的复质量修正只要保住非零质量与同一反常，仍容许共同角松弛新的物理强角。一般电弱 CP 奇高阶源可能造成残余偏移。真空存在也未给出宇宙到达极小所需的历史。每个相对角支的 PQ 圆上只有一个该类极小，两个自发 CP 共轭支仍可能形成壁。

$E_{\rm QCD}(0;\delta)$ 可能随奇异谱变化。精确联合实例可将 $V_0$ 指定为匹配后 $\bar\theta=0$ 时的总势，再加 $E_{\rm QCD}(\bar\theta;\delta)-E_{\rm QCD}(0;\delta)$。若只在固定重径向模的低能模型给定锁相势，另有任意保 PQ 的 $\Delta V(\delta)$，则任意联合极小满足

$$
|\cos\delta-1/3|\le
\sqrt{\operatorname{osc}\Delta V/(\kappa v_S^2v_T^2)} .
$$

只需比较其能量与 $\delta_0=\arccos(1/3)$ 的值。当能量起伏小于 $4\kappa v_S^2v_T^2/9$，CP 守恒的 $0,\pi$ 均被排除。取起伏为 $\kappa v_S^2v_T^2/36$ 时 $|\sin\delta|\ge\sqrt3/2$。动力学径向质量有限时还需控制其位移。独审另构造所有夸克质量均与 $\delta$ 无关而弱 CP 随其变化的例子，因而纯 QCD 的质量依赖零点能在该例严格恒定。[独立审阅](independent_review.md)

规范动能取 $|\partial S|^2+|\partial T|^2$。令 $w=v_T^2/(v_S^2+v_T^2)$、$\alpha=\phi+w\delta$，共同模式和相对模式的角动能度量分别为 $F^2=2(v_S^2+v_T^2)$、$F_r^2=2v_S^2v_T^2/(v_S^2+v_T^2)$。设 $\psi(\delta)=\arg(4-e^{i\delta}/2)-w\delta$。在匹配后的精确锁相势、固定径向下，角势 Hessian 为

$$
\chi\begin{pmatrix}1&\psi'\\\psi'&\psi'^2\end{pmatrix}
+\operatorname{diag}(0,K),\qquad
K=2\kappa v_S^2v_T^2\sin^2\delta_0.
$$

行列式为 $\chi K>0$。相对模式足够重时轻本征值趋近 $\chi/F^2$。若 $\chi$ 采用真实轻夸克物理、无额外势且处于该重模极限，得到已有 QCD 轴子质量标定约 $5.7\,\mu{\rm eV}(10^{12}\,{\rm GeV}/F)$。见证的夸克质量未匹配到真实 QCD，不能直接套用此数。举例若 QCD 势取 $\chi(1-\cos\bar\theta)$，额外加 $-\xi^4\cos(k\phi+\beta)$，驻点仅满足 $|\sin\bar\theta|\le |k|\xi^4/\chi$；只有控制修正的能量与分支才可把它转成近零强角界。PQ 精度及额外拓扑扇区仍有独立物理内容。[物理 QCD 标定](https://arxiv.org/abs/1511.02867)

## 同一载体的电弱与轴子味流

在树级匹配、固定相对角和普通双态的纯径向 Higgs 方向，$P=\operatorname{diag}(1,1,1,0)$ 给出 $v\partial_v{\cal M}=P{\cal M}$、$\partial_\alpha{\cal M}=i(1-P){\cal M}$。对角化后 $X=V^\dagger V$，$vY_h=X\widehat m$、$F Y_a=i(1-X)\widehat m$。两者和式为 $vY_h+(F/i)Y_a=\widehat m$。弱单态投影 $C_L=1-X$ 秩一，$|(C_L)_{ds}|^2=(1-X_{dd})(1-X_{ss})$。可通过手征重定相写成 $(\partial_\mu a/F)\bar d_L C_L\gamma^\mu d_L$ 的形式，整体符号依赖重定相约定。额外独立 PQ 带荷物质、不同荷的右手下味、双态与单态径向混合、或相对角在轻轴子中不可忽略，均会改变相应读数。

若在现实介子匹配后该树级味流主导、轴子在实验中稳定且质量可忽略，定义 $q_\mu\langle\pi|\bar d\gamma^\mu s|K\rangle=(m_K^2-m_\pi^2)f_0(0)$（使用相应共轭味流亦同）。两体衰变宽度为

$$
\Gamma(K^+\to\pi^+a)=
\frac{|(C_L)_{ds} f_0(0)|^2}{64\pi F^2}
\frac{(m_K^2-m_\pi^2)^3}{m_K^3}.
$$

E787/E949 给出历史无质量窄线分支比上限 $0.73\times10^{-10}$（90% CL）。用登记的介子中心质量、寿命得 $F/(|(C_L)_{ds} f_0(0)|)>3.4654\times10^{11}$ GeV。若味流为零，实验不约束 $F$。该式保留强子形状因子；替换真实净幅、相对支混合或可见轻粒子需要重新匹配。此实验限与见证任意单位的轻质量未合并拟合。[E787/E949 原始线搜索，第 94 页](https://arxiv.org/pdf/0709.1000#page=94)；[带电 K 资料](https://pdg.lbl.gov/2024/listings/rpp2024-list-K-plus-minus.pdf)；[带电 π 资料](https://pdg.lbl.gov/2024/listings/rpp2024-list-pi-plus-minus.pdf)。
