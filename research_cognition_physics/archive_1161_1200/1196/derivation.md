# 第1196轮参考推导：同一1178父过程的自洽占据准备

2026-10-14。为当前Jacobson入口提供独立选择的候选参考；作为1196参考部分验收，不新增认知公理，不修改既有验收。代码为 [hartree_reference.py](hartree_reference.py)，默认只读输出JSON，`--write`只允许首次创建独立结果。

## 1. 对象与去重

保留[1176完整CAR父对象](../1176/proof.md)及[1178互易项](../1178/proof.md)，固定其 $65^3$ 周期空间、$\ell=16/65$、$m=2$、$\epsilon=3/20$、$\nu=7/20$、$g=1/5$、$\Omega=1/2$、$E_0=0$。背景每复制的源有 $4V$ 个CAR模式，完整Fock及正谱全局基线均保留。这里不改变源正常序，也不将原空真空重定义为填充海。

新检验是：从**原相互作用的Hartree能量驻点**独立选择已知占据准备，再求真实空间子代数的参考。307—311的自由链方法只作工具，不重报发现；1176旧占据 $\operatorname{sgn}(k_x)\alpha_x$ 不是本页负能谱投影。检索当前阶段的proof、编号笔记和准入README未发现本自洽固定点已被验收。

本文选择了背景准备，尚未制造其局域制备器。已知背景可以独立选择，未知原探针及其参考不能据此被复制、后选择或更换读口。

## 2. 完整能量变分先于参考

记 $\langle\cdot\rangle_p=V^{-1}\sum_p$、$k_b=\sin(\ell p_b)/\ell$，$C_p$为源占据协方差，定义

$$
E_m=\langle\operatorname{tr}(C_p\beta m)\rangle_p,\qquad
E_b=\langle\operatorname{tr}(C_p\alpha_b k_b)\rangle_p.
$$

平移均匀Hartree能量（每站、每复制，另保完整全局基线）为

$$
\mathcal E=E_m+\sum_bE_b+z_s(\epsilon E_m+\nu\sum_bE_b)
+g z_t(E_x-E_y)-\Omega(x_s+x_t)-\epsilon\Omega z_sx_t.
\tag{1}
$$

末项来自1178真实互易Hamiltonian，不能删去。按约定以 $\langle\operatorname{tr}(\delta C_p\,\delta\mathcal E/\delta C_p)\rangle_p$ 表示源变分，则

$$
\frac{\delta\mathcal E}{\delta C_p}
=h_p=N\beta m+v_x\alpha_xk_x+v_y\alpha_yk_y+v_z\alpha_zk_z,
\quad N=1+\epsilon z_s,
\quad(v_x,v_y,v_z)=(1+\nu z_s+gz_t,1+\nu z_s-gz_t,1+\nu z_s).
\tag{2}
$$

关系变分同时给

$$
h_s=-\Omega X_s+(\epsilon E_m+\nu\sum_bE_b-\epsilon\Omega x_t)Z_s,
\qquad h_t=-\Omega(1+\epsilon z_s)X_t+g(E_x-E_y)Z_t.
\tag{3}
$$

这三份响应来自同一能量，不另填写来源。有限K已知复制乘积的精确能量期望把式(1)的全部跨复制交互项乘 $(K-1)/K$；Hartree取其极限。父H的全部项及基线原样存在，本文没有把驻点称为有限K精确本征态。

## 3. 各向同性自洽参考与唯一性

取 $z_t=0,x_t=1$，每动量准备完整秩2负能谱投影

$$
C_p(z)=\frac12\left(I-\frac{h_p(z)}{E_p(z)}\right),\quad
N=1+\epsilon z,\quad v=1+\nu z,\quad
E_p=\sqrt{m^2N^2+v^2k^2}.
\tag{4}
$$

$mN>0$，故无零能歧义。每个 $C_p$ 的占据谱恒为 $(0,0,1,1)$；源总粒子数为 $2V$，没有只保正能支或舍弃某个旋量。各向同性给 $E_x=E_y$，因此式(3)的Gt低能Bloch方向确为 $+X$。

源及Gs自洽条件归为

$$
E_m=-2m^2N\langle E_p^{-1}\rangle_p,\quad
E_{\rm hop}=-2v\langle k^2/E_p\rangle_p,
\quad b(z)=\epsilon E_m+\nu E_{\rm hop}-\epsilon\Omega,
\quad F(z)=-\frac{b(z)}{\sqrt{\Omega^2+b(z)^2}},\quad z=F(z).
\tag{5}
$$

直接求导，所有交叉项合并为

$$
b'(z)=-2m^2(\epsilon-\nu)^2\langle k^2/E_p^3\rangle_p\le0,
\qquad F'(z)=-\frac{\Omega^2b'(z)}{(\Omega^2+b(z)^2)^{3/2}}.
\tag{6}
$$

对整个 $[-1,1]$，$N,v$皆正，且

$$
|b|\ge b_*:=2\epsilon\frac{1-\nu}{1+\nu}m+\epsilon\Omega=\frac{131}{360},
\qquad |b'|\le B_*:=\frac{4m(\epsilon-\nu)^2}{3\sqrt3(1-\nu)^3}
=\frac{2560}{6591\sqrt3}.
\tag{7}
$$

第二界用 $s/(m^2+s)^{3/2}\le2/(3\sqrt3m)$。因此

$$
0\le F'(z)\le\frac{\Omega^2B_*}{(\Omega^2+b_*^2)^{3/2}}<\frac6{25}.
\tag{8}
$$

代码将最后不等式平方后用精确有理数核验。$F([-1,1])\subset(0,1)$，压缩定理给唯一各向同性固定点。该结论不是采样推测，不证明所有非均匀或各向异性准备中的全局最低能量。

固定原参数的浮点诊断为

$$
z_*\simeq0.989835554078098,\quad x_s\simeq0.142216651215337,
\quad N_*\simeq1.148475333111715,\quad v_*\simeq1.346442443927334.
$$

$b_*^{\rm root}\simeq-3.480026936435756$，此处导数 $F'\simeq0.0001305224563$。固定点源与其自身 $h_p$ 对易，两种关系自旋也是各自有效场本征态，故该Hartree准备平稳（整体相位保留）。它不是为匹配boost而选的Gibbs态。

## 4. 两个真实空间站的非平凡源区域

选同一复制的源CAR代数，区域为 $(0,0,0)$ 与 $(1,0,0)$ 两站。这里共有8个复CAR模式，区域Fock维数256。使用

$$
(C_A)_{i\alpha,j\beta}=\langle c^\dagger_{j\beta}c_{i\alpha}\rangle
=\frac1V\sum_p e^{ip\cdot(x_i-x_j)}(C_p)_{\alpha\beta}.
\tag{9}
$$

跨站关联的非零性也有解析证据，不仅靠小数：按式(9)的方向约定，

$$
\frac{\operatorname{tr}(\alpha_x C_{01})}{4i}
=\frac{v}{2\ell}\left\langle\frac{\sin^2(\ell p_x)}{E_p}\right\rangle_p>0.
$$

奇项由动量反射消去，有限动量表含非零 $\sin(\ell p_x)$，故严格正。它是同一区域的偶CAR关联，不是免费授予新的实际检测器。

有限求和和8阶谱分解给相关谱约为四份 $0.24445070843193$ 与四份 $0.75554929156807$；跨站块算子范数约 $0.187338020170727$。所以这份参考有真实跨站关联。区域源态高斯，模矩阵由

$$
k_A=\log[(I-C_A)C_A^{-1}]
\tag{10}
$$

独立算出；多体模生成元还含正常归一化常数。代码重建协方差，只作为代数核验，未将这一步当物理boost证明。

### 不依赖浮点谱的严格满秩保证

对任意两站旋量向量 $u=(u_0,u_1)$，其Fourier多项式 $u(p_x)$不依赖 $p_y,p_z$。在固定 $p_x,p_z=0$ 上只保留 $p_y=\pm2\pi16/65$ 两项。令

$$
a_{\max}^2=[m(1+\epsilon)]^2+[(1+\nu)/\ell]^2,\qquad
b_{\min}^2=[(1-\nu)(999/1000)/\ell]^2.
$$

用 $\sin(32\pi/65)=\cos(\pi/130)>999/1000$（由 $\pi^2<10$）及Clifford关系，这对谱投影之和的最小特征值至少

$$
d_*:=\frac{b_{\min}^2}{2(a_{\max}^2+b_{\min}^2)}>0.
$$

离散Parseval给 $\sum_{p_x}\|u(p_x)\|^2=65\|u\|^2$。对正负谱投影完全相同，故

$$
\frac{d_*}{65^2}I\le C_A\le\left(1-\frac{d_*}{65^2}\right)I,
\qquad \lambda_{\min}(\sigma_A)\ge\left(\frac{d_*}{65^2}\right)^8>0.
\tag{11}
$$

代码以精确分数保存这个保守界。数值多体最小本征值约 $1.27505935414\times10^{-5}$，远强于保守式(11)，但不把浮点值提升成已认证的区间下界。

## 5. 有限K能够复用什么

源参考跨站相关，故不能复用1178真空下的空间独立计数改进。应回到1176的**整个复制cluster**，新的互易项给

$$
w\le2V(2\epsilon m+6\nu/\ell+4g/\ell)+2V\epsilon\Omega,
\qquad \mathcal L\le V[2\epsilon(.3+.5)+18\nu/\ell+12g/\ell].
$$

在同类已知复制乘积准备及固定时间窗下，旧计数／Duhamel工具提供

$$
\alpha(T)\le\frac{e^{12wT}-1}{K},\quad
D(\rho_{A,K}(T),\sigma_A)\le\sqrt{\alpha(T)},\quad
\eta_{P,K}\le\frac{4\mathcal LT e^{6wT}}{\sqrt K}.
\tag{12}
$$

这里 $\rho_{A,K}$ 是**精确BG中一个复制的两站源减态**；若原未知探针也加入，另付同父探针反作用比较误差，不能忽略。未知原48码/R的保留属于最后的等距比较，而不是“非线性Hartree映射就是未知量子信道”。此全cluster界极粗；不复制1178的 $K=10^{14}$ 数值证书，不称可实现工程资源。

式(12)不自动控制对数。若另有固定的 $\sigma_A\ge\mu I$、$\|\rho_{A,K}-\sigma_A\|_{\rm op}\le e<\mu$，才可用

$$
\|\log\rho_{A,K}-\log\sigma_A\|_{\rm op}\le e/(\mu-e).
\tag{13}
$$

半迹距给 $e\le2\sqrt\alpha$ 仅为足够条件。式(11)可提供同对象的保守正下界，但本页未选K、时长和精度去签收整份物理模流协议。

## 6. 范围与下一步

1. 得到的是**源CAR子代数**参考及其模矩阵。若把Gs/Gt或原空探针也纳入同一完整区域，Hartree因子仍含纯态；不能因此声称整个父区域忠实，也不能把受限源模流当成所有部门的完整boost。
2. 固定点不是有限K全H的精确基态或Gibbs态，也未证明它是物理真空。可以作为有限时准备候选，实际准备／记录资源仍须交代。
3. 完整能量变分与模能量身份不同。式(2)核对了原同父来源，源参考与独立能源权重的有限比较见[本轮证明](proof.md)；不由两站熵匹配领取面积律、纠缠平衡或引力。
4. 1185匹配来源含 $\epsilon(I-X_t)/2$ 的中央部分；其对应固定点方程与本页 $E_0=0$ 不同，不能直接套本页数值。局部编码运输可复用，不可丢中央能源。
5. 本轮已接[有限boost比较](proof.md)。后续见[NEXT](NEXT.md)，不按目标boost重选态，不强修候选，不将其失败升级为Jacobson路线失败。
