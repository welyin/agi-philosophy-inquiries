# 未编号候选：有限局部 CAR 多副本到同一三维钟与位置运输

2026-10-13。固定参数可行性筛选，尚未准入1175；本文件不是研究结项或新公理。数值部分计算自主 Hartree 比较器，有限父过程由[独立误差证明](../mean_field_bound.md)连接，绝不把比较器数值称为已经模拟了巨大父系统。

## 1. 精确有限父系统与资源

取三维周期格点 $\Lambda=(\mathbb Z/65\mathbb Z)^3$，$V=65^3$，周期 $L=16$，格距 $\ell=16/65$。每个复制 $r=1,\ldots,K$ 包含全空间管内的完整源 CAR Fock 空间（每站4个 Dirac 模式，源质量 $m_S=2$）及每站一个关系自旋。复制不是空间坐标。另有完整探针 CAR Fock 空间，每站4个 Dirac 模式乘2个内部质量标签，$M_P=\operatorname{diag}(0.3,0.5)$。

采用 $\beta=\sigma_z\otimes I$、$\alpha_b=\sigma_x\otimes\sigma_b$。空 CAR 真空约定固定，不将占据制备叫成 Dirac 海。定义完整二次量子化算符

$$
m_i^a=c_{ai}^\dagger\beta M_a c_{ai},\qquad
t_{ij}^a=-\frac{i}{2\ell}c_{ai}^\dagger\alpha_b c_{aj}+\mathrm{h.c.},\qquad
t_i^a=\frac12\sum_{j:j\sim i}t_{ij}^a,\qquad
F_i^a=\epsilon m_i^a+\nu t_i^a.
$$

边项以正坐标方向定义一次，其在相邻两站的半分是同一 Hermitian 算符，不再重复变号。固定 $\epsilon=0.15,\nu=0.35,\Omega=0.5,T=2$。一份 cluster 及两份 cluster 之间的算符为

$$
h=\sum_i(m_i^S+t_i^S)-\Omega\sum_iX_i,
\qquad W_{rs}=\sum_i(Z_{ir}F^S_{is}+F^S_{ir}Z_{is}),
$$

$$
H_{BG,K}=\sum_rh_r+\frac1K\sum_{r<s}W_{rs},\qquad
H_{\rm act}=C_KI+H_{BG,K}+H_{P,0}
+\sum_i F_i^P\overline Z_i,\quad
\overline Z_i=K^{-1}\sum_rZ_{ir}.
$$

明确排除 $r=s$，故独立界式(2)直接适用，不额外忽略自作用。有限源复制 $r$ 的系数是 $1+\epsilon K^{-1}\sum_{s\ne r}Z_{is}$；探针是 $1+\epsilon\overline Z_i$，跳跃系数同理用 $\nu$。这项 $1/K$ 差别在独立证明中已计，不宣称有限 $K$ 每个子系统的系数完全相同。所有跨复制物质算符均为 CAR 偶算符，空间支持仅同站或最近邻边。固定 $K$ 整体有限、自伴、时间不变，整体能量精确守恒；$C_K$ 只把全谱统一平移至正值，并非免费制备或真空能推导。

可使用下列粗而严格的完整 Fock 范数界：

$$
b_S=2\epsilon m_S+6\nu/\ell,\quad
b_P=2\epsilon(0.3+0.5)+12\nu/\ell,\quad
\|W\|\le w=2Vb_S,\quad
\mathcal L=\sum_i\|F_i^P\|\|Z_i\|\le Vb_P.
$$

例如 $C_K=1+\Omega VK+KV[2m_S(1+\epsilon)+6(1+\nu)/\ell]+V[2(0.3+0.5)(1+\epsilon)+12(1+\nu)/\ell]$ 足够。实际与完整比较器同时保留相同的 $e^{-iC_Kt}$，不在准备等距向量界中丢掉一侧的共同相位。所有 $K$ 份物质和自旋、完整负单粒子支、准备能量及末读记录保留。源 $F_i=\epsilon m_i+\nu t_i$ 是所选共同关系路径的完整作用导数；当 $\epsilon\ne\nu$ 时不将它偷称为一个普适总能量密度。整体能量账仍包括全部源、探针、关系自旋及交互。

## 2. 指定背景准备与精确闭合的 Hartree 方程

每份源准备为每站两个 $\beta=+1$ 自旋模式均占据的 Slater 态；每个关系自旋准备 $|+X\rangle$。不同复制之间为固定已知乘积，探针可任意未知且与参考 $R$ 纠缠。这里没有声称未知或跨复制猫态背景也可被同一轨道替代。

平移均匀的 Hartree 一份 cluster 解仍分解为源 Slater 和关系自旋乘积。令 $\mathbf g=(x,y,z)$，$N=1+\epsilon z$、$v=1+\nu z$，每个格点动量 $\mathbf p$ 的 $s_{\mathbf p}=\ell^{-1}\sqrt{\sum_b\sin^2(\ell p_b)}$。两重自旋简并将源协方差精确化成 Bloch 向量 $\mathbf r_s$：

$$
\dot{\mathbf r}_s=2(vs,0,Nm_S)\times\mathbf r_s,
\quad E_m=2m_S\langle r_z\rangle_{\mathbf p},\quad
E_t=2\langle s r_x\rangle_{\mathbf p},
\quad b=\epsilon E_m+\nu E_t,
\quad \dot{\mathbf g}=2(-\Omega,0,b)\times\mathbf g.
$$

初始 $\mathbf r_s=(0,0,1)$、$\mathbf g=(1,0,0)$，故 $b(0)=0.6$、$z''(0)=-4\Omega b=-1.2$。每站平均能量 $NE_m+vE_t-\Omega x$ 守恒。全部 $65^3$ 个源动量按 $|\sin p|$ 与坐标置换精确归并成6545类，权重总数仍是 $V$；不是把源冻在 $p=0$。

## 3. 同一内部钟、实际位置与有限比较器

同一探针 Hamiltonian 的动量块为

$$
h_{P,\mathbf p}(t)=N(t)\beta M_P+v(t)\boldsymbol\alpha\cdot\boldsymbol\ell_{\mathbf p}.
$$

探针准备码是32维：四动量 $\mathcal B=\{0,(1,0,0),(0,1,0),(0,0,1)\}$（动量单位 $2\pi/L$）乘完整4维 Dirac 自旋及2维内部钟。等距映射 $E_{\mathcal B}|p,d,m\rangle=V^{-1/2}\sum_{j\in\Lambda}e^{ip\cdot x_j}c^\dagger_{jdm}|0\rangle$；该码中的所有幅度、内部关联及参考 $R$ 均可未知。等幅空间包、$\chi=(1,0,0,1)/\sqrt2$ 和内部 $|+\rangle$ 仅是预先固定的数值见证。数值同时算整个32维准备空间传播矩阵，没有只保此见证态或正频支。

末读是实际有限格点位置读口加内部 $-Y$ 钟读口；完整位置和钟结果均保留。预先固定的有限汇总效果为 $f_b(j)=[1+\sin(2\pi j_b/65)]/2$，以及 $f_b\otimes P_{-Y}$。它们可由站点细记录后的已知随机后处理得到；现实制备、末读及其资源需要另行声明有限实现误差，不由写下 POVM 自动授予。此筛选只认证终端菜单，中途背景干预后不再免费沿原 Hartree 轨道继续。

独立 $p=0,\beta=+1$ 钟只校准 $\tau=\int_0^T Ndt$，不能单凭末读反推出完整 $N(t),v(t)$。移动探针的内部钟与位置联合结果由同一 $h_{P,\mathbf p}$ 计算，不能用另一静止钟替代。平直对照取同一输入、同一空间标签，传播固有时 $\tau$；因此剥除了单纯统一钟速重标。质量零实验是同一系数规则的另一个探针族，不称与有质量钟同一个载体。

在声明的低动量带上，将 $\ell_{\mathbf p}$ 换成 $\mathbf p$ 的传播误差不超过

$$
\eta_{\rm lat}=T(1+|\nu|)\frac{\ell^2}{6}
\max_{\mathbf p\in\mathcal B}\sqrt{\sum_b p_b^6}.
$$

连续平面波比较器与离散 $f_b$ 的压缩效果矩阵精确相同：有限 Fourier 次数远低于格点 Nyquist 界。格点位置 PVM 的完整 Fourier 运输是有限 cutoff 上的采样位置仪器，不能直接叫连续空间精确点位置 PVM。若运输完整细记录及后态，应使用全 $65^3$ Fourier 模式上的等距，将每个格点投影、指针记录及 Stinespring 映射一起运输；它与初始32维准备码不同。当前连续几何层的直接概率比较只需声明的平滑 $f_b$/钟联合菜单；有限 $K$ 父桥则允许完整离散末读和同一后续 CP 操作。末读后覆盖其它动量时，不再免费领取原低带的连续传播近似。这里没有用不连续半空间端口偷省采样误差，也没有将完整末读后的无限未来投影回初始动量带。

## 4. 几何资格与删弱

连续比较器是经典共形平直 FLRW 字典

$$
ds^2=N(t)^2dt^2-a(t)^2d\mathbf x^2,\qquad a=N/v.
$$

这是由共同完整质量与空间导数两系数解释的有限低动量字典；尚不是 Einstein 动力学。经 $\psi=a^{3/2}\Psi$ 场重标，守恒概率内积 $\int a^3\Psi^\dagger\Psi\,d^3x$ 变成 $\int\psi^\dagger\psi\,d^3x$，自由 Dirac Hamiltonian 为 $N\beta M+(N/a)\boldsymbol\alpha\cdot\mathbf p$。这一成熟字典可从 [Obukhov、Silenko、Teryaev 的一般 Hermitian Hamiltonian，式(2.18)–(2.21)](https://arxiv.org/pdf/1106.0173) 代入 $V=N,W^a{}_b=a\delta^a_b,K=0$ 得到：空间旋量不变量消失，时间连接的对称部分由场重标吸收。这里采用一般时变公式，不把该文随后平稳专化式(2.24)直接外推。当前源不必取得连续物质解释。

初时 $N=v=a=1,\dot a=0$，而 $\ddot a(0)=(\epsilon-\nu)z''(0)=0.24$；平直 FLRW 的 Kretschmann 标量在此为 $12\ddot a^2=0.6912$。因此**已认证的比较器**不是共同参考换表即可消去的平直几何，不能单凭某端口概率改变领取此结论。实际父过程与比较器的有限任务误差须另外使用第5节。

删弱一：令 $\nu=\epsilon$，Hartree 空间均匀时 $a=1$，比较器传播精确等于平直生成元经过共同时间 $\tau$；不是任何空间变化 lapse 都平直。删弱二：源为空 CAR 真空，不含探针的背景 $+X$ 稳定，Hartree 比较器 $N=v=1$。有限 $K$ 实际父系统仍有局部关系涨落和探针反作用，不因这两个删弱就精确化成平直传播；其完整差异只由同一 $\eta_K$ 约束。这两项用同一准备/读口算法核查，不改端口追逐差值。

## 5. 分开登记的误差与当前停止线

1. 有限 $K$：直接用独立界 $\eta_K\le4\mathcal L T e^{6wT}/\sqrt K$。比较背景保留精确 $e^{-itH_{BG,K}}\phi_0^{\otimes K}$，不要求整个背景接近乘积。固定 $\eta_K$ 后只记录足够 $K$ 的 $\log_{10}K$，不优化巨大资源量。
2. 空间格点：上式 $\eta_{\rm lat}$ 与已声明端口 Fourier 压缩关系；源的整个高动量模式族并未被连续近似。
3. 数值积分：当前 RK4 及步长加倍只是可行性诊断；需要严格有限数值误差证书后才可签正式概率界，不用收敛差冒充解析误差。
4. 准备/钟/位置仪器：有限完整实现误差是公开权限输入，尚未由本例造出仪器。全部记录与反作用保持；稀有后选择不享受无分母误差界。

目前不领取三维空间生成、未知背景普遍经典化、GR硬交付完成、正常真空或 Lorentz 连续极限。净检验对象是：同一局部量子父过程的既定来源准备，能否在受控有限任务域给出自主改变的共同钟/空间运输，并排除它只是统一钟速重标。

## 6. 历史复用与净增位置

- [356](../../../../archive_342_369/research_note_356.md) 已有全块/R均场对照；[357](../../../../archive_342_369/research_note_357.md) 已有非交换自主均场，不能再把自主 Bloch 方程记成新结论。
- [525–530入口](../../../../archive_467_530/research_note_525.md) 已有参考自身动力学、校准与完整代价；不给本例实际三维位置资格。
- [758](../../../../archive_742_763/research_note_758.md)、[759](../../../../archive_742_763/research_note_759.md)、[760](../../../../archive_742_763/research_note_760.md)、[761](../../../../archive_742_763/research_note_761.md) 已有有限窗口的经典/量子比较限制；本例同样不抹除背景关联。
- [1174](../../../research_note_1174.md) 是主动局部 lapse 接口；新候选改变 $N/v$，并必须用同一实际移动钟/站点资料和量化比较承担这项新资格。

只有上述具体连接及误差认证成立才值得准入新轮；目前先保存可复算候选，未编号。

## 7. 首次固定参数预览（未认证数值）

[preview.py](preview.py) 已按上述先定参数运行，显式 write-once 保存到 [results_preview.json](results_preview.json)，默认只读复算通过。使用2048步 RK4，4096步仅作为诊断对照。没有改变质量、时间、端口或输入态来追逐差值。

| 项目 | 首次预览 |
|---|---:|
| $a(T)$ | 1.0237655983 |
| 原固定输入，$x$ 位置与 $-Y$ 钟联合概率差 | -0.0020966890 |
| 同一输入，$y$ 位置与 $-Y$ 钟联合概率差 | -0.0025464628 |
| 同一输入，$z$ 联合概率差 | -0.0004499601 |
| 质量零另族，$x$ 位置概率差 | -0.0056962818 |
| 三轴联合效果在完整32维码上的算符距离 | 0.0150246448 |
| $\nu=\epsilon$ 比较器与其固有时平直传播的准备算符差 | $2.63\times10^{-15}$ 以下 |
| 真空源比较器同项 | $3.14\times10^{-15}$ 以下 |

固定 $\chi$ 是 $\alpha_x$ 的正本征态，所以不期待质量零三轴定向平均都非零；三维资格来自同一三轴 Hamiltonian、完整准备码和三轴效果矩阵，而非把一个轴向态包装成各方向都相同的见证。没有按效果最大特征态重新选输入。

首次粗空间界为 $\eta_{\rm lat}=0.0016512218$，平直对照界约0.0011528755；它们之和覆盖了固定 $x/y$ 联合差，故不能用这两个粗数签出稳健联合差。由于此处 $v(t)>0$，同参数更精确的解析界是

$$
\eta_{\rm lat,act}\le\frac{\ell^2p_*^3}{6}\int_0^T v(t)dt,
\qquad
\eta_{\rm lat,flat}\le\frac{\ell^2p_*^3}{6}\tau,
\qquad p_*=2\pi/L.
$$

这不是改变物理或重新拟合，而是用同一自主轨道的已认证积分替代最大值上界；在积分数值认证之前仍仅是待签预算。独立残差证明正在给出该认证。按目前预览，$y$ 联合差有望在此收紧后及公开仪器误差下留下正余量；$x$ 联合差暂不签。充分 $K$ 的粗界极端昂贵：取 $\eta_K=10^{-4}$ 得 $\log_{10}K\ge52275241.074$，只作为数学有限存在性与全部资源账，绝不称为现实可行制备。

## 8. 正式准入采用的公开直和扩展

首次预览保留不改。经主线同意，后续正式模型把 $M_P=\operatorname{diag}(0,0.3,0.5)$ 放在同一个固定父 Hamiltonian 中：每站12个 CAR 模式，准备码为四动量乘4维 Dirac乘3个质量标签，共48维。旧两质量钟实验和质量零实验仅是同一模型的不同准备，而非更换 Hamiltonian。原源、参数、见证和数据均不变，所有跨质量相干也保持；钟效果在后两标签的 $-Y$ 块，质量零标签读位置。

[unified_probe_extension.py](unified_probe_extension.py) 及其 [write-once 结果](results_unified_probe.json) 给48列 Gram 误差 $4.72\times10^{-15}$、旧菜单最大差 $1.12\times10^{-16}$。新的完整 Fock 资源界为

$$
b_P=2\epsilon(0.3+0.5)+18\nu/\ell=25.83375,
\quad\mathcal L=7094593.59375,
\quad \log_{10}K\ge52275241.4213
$$

即可达到原 $\eta_K=10^{-4}$ 目标；$C_K$ 的 probe 跳跃系数相应由12改18。可正式选整数 $K=10^{52275242}$，不在计算机中物化这个数或假装已搭建父体系。此扩展只是同对象合同澄清，不另计科研组。
