# 1175 解析：有限自主背景与共同三维运输

2026-10-13。作者稿，待独立审阅；以下数值结论以保存结果及公开基本浮点算术模型下的残差包络为条件。

## 1. 同一个完整有限父过程

令 $\Lambda=(\mathbb Z/65\mathbb Z)^3$，$V=65^3$，$\ell=16/65$。每份 cluster 含完整源 CAR（每站4个 Dirac 模式，质量2）和每站一个关系 qubit；共 $K$ 份。另有完整探针 CAR，每站12模，内部质量 $M_P=\operatorname{diag}(0,0.3,0.5)$。在空 CAR 真空的固定正规序约定下保留全部 Fock 扇区，不引入 Dirac 海或删除负单粒子支。

取 $\beta=\sigma_z\otimes I_2$、$\alpha_b=\sigma_x\otimes\sigma_b$，固定 $\epsilon=.15,\nu=.35,\Omega=.5,T=2$。记

$$
m_i^a=c_{ai}^\dagger\beta M_a c_{ai},\quad
t_{i,i+\hat b}^a=-\frac{i}{2\ell}c_{ai}^\dagger\alpha_b c_{a,i+\hat b}+\mathrm{h.c.},\quad
t_i^a=\frac12\sum_{e\ni i}t_e^a,\quad F_i^a=\epsilon m_i^a+\nu t_i^a.
\tag{1}
$$

每条边按正坐标方向只定义一次；相邻两站半分同一个 Hermitian 边项。cluster 单体及交换对称偶耦合是

$$
h=\sum_i(m_i^S+t_i^S)-\Omega\sum_iX_i,
\quad W_{rs}=\sum_i(Z_{ir}F_{is}^S+F_{ir}^SZ_{is}),
\quad H_{BG,K}=\sum_rh_r+K^{-1}\sum_{r<s}W_{rs},
\tag{2}
$$

$$
H_K^{\rm act}=C_KI+H_{BG,K}+H_{P,0}
+\sum_iF_i^P\overline Z_i,
\quad H_{P,0}=\sum_i(m_i^P+t_i^P),\quad
\overline Z_i=K^{-1}\sum_rZ_{ir}.
\tag{3}
$$

这是时间不变的完整有限 Hamiltonian；关系、源、探针及交互能全部在同一总能量账。跨复制耦合的空间支持仅同站/最近邻边，在每个 CAR cluster 上都是偶算符。源第$r$份的系数排除自身$r=s$；探针取全部平均。因此有限K源与probe差一个 $1/K$ 项，不宣称它们在微观上精确同系数。

完整 Fock 范数由单粒子正/负谱之和给出：源 $\|m_i^S\|=2m_S$、$\|t_e^S\|=2/\ell$；探针 $\|m_i^P\|=2(0.3+0.5)$、$\|t_e^P\|=6/\ell$。于是

$$
b_S=2\epsilon m_S+6\nu/\ell,\quad
b_P=2\epsilon(0.3+0.5)+18\nu/\ell,\quad
w=2Vb_S\ge\|W\|,\quad \mathcal L=Vb_P\ge\sum_i\|F_i^P\|\|Z_i\|.
\tag{4}
$$

可选共同正谱平移

$$
C_K=1+\Omega VK+KV[2m_S(1+\epsilon)+6(1+\nu)/\ell]
+V[2(0.3+0.5)(1+\epsilon)+18(1+\nu)/\ell].
\tag{5}
$$

实际和比较器都保留同一 $e^{-iC_Kt}$。常数平移不制造正常真空、GR真空能或免费制备。物质的实际来源 $F_i=\partial H/\partial z_i$ 是所选路径的质量/动量通量组合；$\epsilon\ne\nu$ 时不冒充普适总 $T_{nn}$。

## 2. 同源自主 Hartree 与完整背景误差

每份源准备为每站两个 $\beta=+1$ 自旋模均占据的 Slater，关系自旋皆为 $+X$；不同复制初态 $\phi_0^{\otimes K}$，与未知探针/参考独立。Hartree 方程

$$
i\dot\phi=(h+w^\phi)\phi,\quad
w^\phi=(I\otimes\langle\phi|)W(I\otimes|\phi\rangle)
\tag{6}
$$

分成源二次 Hamiltonian 与单站关系自旋 Hamiltonian，因此该指定空间均匀族精确闭合成 Slater协方差与Bloch向量。动量 $\mathbf p$ 对应 $s=\ell^{-1}\sqrt{\sum_b\sin^2(\ell p_b)}$。完整源 Brillouin 区按符号/置换精确归并为6545类，并保整数重数之和 $65^3$。

$$
\begin{aligned}
N&=1+\epsilon z,& v&=1+\nu z,\\
\dot{\mathbf r}_s&=2(vs,0,Nm_S)\times\mathbf r_s,
& E_m&=2m_S\langle r_z\rangle,\quad E_t=2\langle s r_x\rangle,\\
b&=\epsilon E_m+\nu E_t,&
\dot{\mathbf g}&=2(-\Omega,0,b)\times\mathbf g .
\end{aligned}
\tag{7}
$$

初始 $\mathbf r_s=(0,0,1),\mathbf g=(1,0,0)$；每站能量 $NE_m+vE_t-\Omega x$ 守恒。源没有被冻在零动量，背景本身的自主推进也不是外给程序表。

复用已独立审查的[有限复制证明](../_admission/after1174_direction/mean_field_bound.md)，对任意探针/参考输入准备等距 $E_K$，有

$$
\left\|(U_K^{\rm act}(T)-U_K^{\rm cmp}(T))E_K\right\|_{\rm op}
\le\eta_K=\frac{4\mathcal L T e^{6wT}}{\sqrt K},
\tag{8}
$$

其中比较器背景仍是**精确** $e^{-itH_{BG,K}}\phi_0^{\otimes K}$，probe由式(7)的 $N,v$ 推进；不是整体背景乘积近似。式(8)保完整BG后态、源、任意未知probe/R及绝对相位；实际probe反作用没有删除。相同后续完整CP仪器的联合概率总变差不超过 $\eta_K$。但中途改变背景后再套原轨道须另证；广延总能量的绝对差也不由一个小态误差自动趋零。

这里 $w=5015339.0625,\mathcal L=7094593.59375$。取

$$
K=10^{52275242}
\tag{9}
$$

即可确保 $\eta_K\le10^{-4}$；只记录对数预算，不假装已构造或计算此父体系。$K$ 的准备、空间支持及 $O(K)$ 能源全部是巨大但有限的模型资源，不是现实工程可行性证明。

## 3. 原移动钟、三维位置与质量零属于同一个对象

在同一完整父probe里取48维正常单粒子准备码

$$
\mathcal B=\{0,e_x,e_y,e_z\}(2\pi/16),\quad
E_{\mathcal B}|p,d,m\rangle=V^{-1/2}\sum_j e^{ip\cdot x_j}c_{jdm}^\dagger|0\rangle.
\tag{10}
$$

码内所有幅度和跨质量相干均可未知，并可与任意参考纠缠。两质量 $0.3,0.5$ 的钟效果是该子块上的 $P_{-Y}$；质量零标签通过另一准备使用同一Hamiltonian的空间读口。没有另换质量零动力学，也没有只保正能支。

同一比较器动量块为

$$
h_{m,p}(t)=N(t)m\beta+v(t)\boldsymbol\alpha\cdot\boldsymbol\ell_p,
\qquad (\ell_p)_b=\sin(\ell p_b)/\ell.
\tag{11}
$$

完整末读保所有站点及钟记录；预先固定汇总 $f_b(j)=[1+\sin(2\pi j_b/65)]/2$，及 $f_b\otimes P_{-Y}$。这是实际有限格点采样位置菜单。准备/触发/指针/末读资源由明确有限实现合同承担，单支误差上限记 $\eta_{\rm inst}$，不是数学POVM已造出仪器。

固定见证：四模式等幅，Dirac自旋 $\chi=(1,0,0,1)/\sqrt2$，钟两标签 $|+\rangle$；质量零另次准备相同空间/Dirac态。$\chi$ 沿$x$定向，不能声称同一massless态三轴定向均非零。其余未知态由全准备算符证书覆盖，不对每个未知态保证同样大小的留出。

静止$p=0$钟只校准 $\tau=\int_0^TNdt$，不独自识别完整轨道。移动钟×位置由式(11)共同计算。平直关系对照在同一输入、同一采样字典下传播固有时 $\tau$，避免把坐标时 $T$ 不同当作空间变化；它也可由同一父H的真空来源准备在校准时长 $\tau<2$ 的比较器实现，有限父过程仍另付式(8)误差。选择与触发这些有限窗口属于共同内部校准/末读合同，不免费冻结活钟。任何共同参考变换都须同时运输状态、生成元、仪器及记录；其概率精确不变。

## 4. 有限低动量字典及不是换参考的判据

在初始带上，$\sin(\ell p_b)/\ell$ 与$p_b$的差给

$$
d_*\le\ell^2p_*^3/6,\qquad p_*=2\pi/16,
\quad \eta_{\rm lat,act}\le d_*\xi,
\quad \eta_{\rm lat,flat}\le d_*\tau,
\quad\xi=\int_0^Tvdt .
\tag{12}
$$

两系数始终正，故可用积分而非最大值界。几何比较仅低动量probe；源的完整高动量模式不作连续近似。平滑$f_b$及钟联合的压缩效果矩阵在离散求和与连续Fourier积分间精确相同。完整格点PVM运输到全cutoff的采样仪器，不能冒充连续精确点位置PVM；末读后不免费继续原低带传播近似。

连续比较器是

$$
ds^2=N^2dt^2-a^2d\mathbf x^2,\qquad a=N/v.
\tag{13}
$$

它的Dirac作用经 $\psi=a^{3/2}\Psi$ 重标得到 $h=Nm\beta+(N/a)\boldsymbol\alpha\cdot\mathbf p$，可直接复用[一般时变 Hermitian Dirac 公式(2.18)–(2.21)](https://arxiv.org/pdf/1106.0173)，不从平稳专化式外推。这里只是共同有效字典，没有输入或导出Einstein作用。

由式(7)，$b(0)=.6,z''(0)=-1.2,N(0)=v(0)=1,\dot N(0)=\dot v(0)=0$，所以

$$
\ddot a(0)=(\epsilon-\nu)z''(0)=.24,\qquad
R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}\big|_0
=12(.24)^2=.6912.
\tag{14}
$$

这是已指定FLRW比较器的非零曲率标量，不能被共同坐标换表消去。实际父过程能在声明任务上接近该字典仍依赖式(8)/(12)和仪器/数值误差。概率变化本身不等于直接测到了式(14)的全部曲率。

## 5. 删弱与应保留的自由

1. 空间均匀 $\nu=\epsilon$ 时 Hartree $a=1$，比较器所有质量块精确为平直生成元经历同一 $\tau$。有限K父过程仍有涨落，只在式(8)误差内复现这一删弱；空间依赖lapse不享有该结论。
2. 空CAR源真空令没有probe的背景 $+X$ 稳定，比较器 $N=v=1$。实际probe仍反作用，不能把全父Hamiltonian叫作精确平直。
3. 一个共同基变换/坐标字典的运输不改变任何完整任务概率；只改H或只改端口不是此对照。
4. $\epsilon,\nu$、三维图、复制机制、来源准备、质量菜单及有效尺度是输入。源导数不是一般总应力，关系自能不被声称已满足GR普适自耦；本轮不选择Einstein动力学或其传播自由。

## 6. 数值包络与实际合同

首次固定预览及公开直和扩展均保存在准入目录，原参数未优化。正式check绑定这些来源和独立[残差推导](../_admission/after1174_direction/mean_field_numerical_bound.md)、[区间复算](../_admission/after1174_direction/mean_field_numerical_bound.py)。对相同2048步RK4轨道构造逐步三次Hermite重建；原方程二次，故残差最高六次，Bernstein系数凸包给上界。令

$$
\mathcal E^2=|\delta g|^2+2\langle|\delta r|^2\rangle,
\qquad F^2=\langle\nu^2s^2+\epsilon^2m_S^2\rangle .
\tag{15}
$$

相同转动的斜对称项消去，交叉项用有限 $F$ 控制；没有用巨大 $\|W\|$ 的指数控制数值误差。基本浮点加减乘除按IEEE舍入和向外扩张处理，trigonometric输入有独立区间包络，精确十进参数也外扩。这里不宣称已形式验证硬件或数学库实现。

该证明给同轨道积分与probe的保守上界

$$
\delta\tau<3.788\times10^{-5},\quad
\delta\xi<8.837\times10^{-5},\quad
\eta_{\rm num,act}<5.364\times10^{-5},\quad
\eta_{\rm num,flat}<3.381\times10^{-5}.
\tag{16}
$$

质量零块直接用 $\exp[-i\boldsymbol\alpha\cdot\boldsymbol\ell_p\xi]$，其误差由 $p_*(\delta\xi+\delta\tau)$ 控制，比式(16)之和更小。因此包括跨质量相干的48维码可以统一使用 $9\times10^{-5}$ 数值概率预算。向量比较转概率时先规范化数值见证（或极分解准备矩阵），把计算Gram/效果累加的额外缺陷付入 $10^{-9}$ 算术余量；该余量及式(16)之和小于统一预算。末端余量同时覆盖库函数相位与端口求积，不能凭小Gram误差替代。

独立审阅已用解析 Fourier 系数 $0,1/2,\pm i/4$ 核对端口，算符差小于 $5.2\times10^{-17}$；同一平直传播的32阶Taylor展开满足矩阵行范数小于2、尾界 $1.053\times10^{-27}$、基本算术误差不超过 $1.69\times10^{-12}$，与作者库函数传播差小于 $1.7\times10^{-16}$。末端整体基本算术包络可保守付 $10^{-10}$，连同探针自身插值残差 $6.53\times10^{-11}$ 的数倍规范化费用，总数仍小于 $10^{-9}$。具体独立计算绑定于本轮审阅资产，不另假设任意libm相位准确，也不把步长加倍差当作严格界。

四类费用分开相加。取每支 $\eta_K\le10^{-4}$，实际内部准备/有限采样末读合同每支 $\eta_{\rm inst}\le2\times10^{-5}$。它包括校准、触发、完整指针记录及反作用的统一误差，当前为新增有限实际权限输入，未由本例造出仪器。使用经过认证的积分上端点，令

$$
\mathcal B_{\rm pair}=d_*(\widehat\xi+\delta\xi+\widehat\tau+\delta\tau)
+2(10^{-4})+2(2\times10^{-5})+9\times10^{-5}.
\tag{17}
$$

同一仪器合同覆盖两种来源准备与原未知码，不只对选中的见证有效。该菜单是非选择的，完整钟/站点记录保留；罕见选择后的条件分布没有同样的无分母界。

给定式(17)的公开实际合同，固定 $y$ 联合菜单和质量零$x$菜单的差异下界分别为

$$
|\Delta p_{y,-Y}|\ge0.00254646277-\mathcal B_{\rm pair},\qquad
|\Delta p_{x,m=0}|\ge0.00569628185-\mathcal B_{\rm pair}.
\tag{18}
$$

精确存档值见 [results.json](results.json)：$\mathcal B_{\rm pair}=0.002542159115238022$，$y$联合与质量零$x$的正下界分别为 $4.3036578\times10^{-6}$、$0.0031541227$。不把舍入后的展示值用于判定；前者余量很窄，必须保留上述全部条件。首次最大值格点界覆盖联合信号的情况被保留，后续仅用同参数积分改进界。固定$x$联合菜单不领取超过总预算的差；也没有改选完整效果的最大特征态增大信号。

这就给出本轮的有限条件命题：在同一局部量子父过程、指定内部来源族、未知48维码及有限实际读口合同中，原钟/位置/质量零菜单可以共同接到由该来源自主产生的非平直有效关系，同时保留完整BG与资源。它不把该机制或其参数说成现行认知原则唯一强迫的结构，不领取普适GR来源约束或一般几何历史。
