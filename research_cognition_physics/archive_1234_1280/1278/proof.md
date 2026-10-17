# 1278：裸热分配交换驱动同父有效关系几何

2026-10-17。参数先保存于[固定规格](fixed_spec.md)，范围依[准入](admission.md)。本件是公开新有效模型的条件性有限证书，不声称新耦合、准备或精密读口由认知原则唯一选择。复用旧通用均场工具，不搬1175或1272数值证书。

## 1. 同一自治父过程与完整资源

沿1175的三维周期图、完整探针CAR、质量菜单、48维准备码及有限格点读口。每个复制簇、每站改用能源对 $A,B$、关系qubit $G$及保留的热准备纯化。复制数记 $K=10^{600000}$；内部已知程序记 $c=0,1$。令
$$
n_A=|1\rangle\langle1|_A,\quad n_B=|1\rangle\langle1|_B,\quad
D=n_A-n_B,\quad X_{\rm ex}=|10\rangle\langle01|+|01\rangle\langle10|.
$$
固定 $\Delta=\lambda=\kappa=\Omega=1$，取
$$
h_c=\sum_i[\Delta(n_{Ai}+n_{Bi})+\lambda cX_{{\rm ex},i}-\Omega X_{Gi}],
\quad
W_{rs}=\kappa\sum_i(D_{ir}Z_{G,is}+Z_{G,ir}D_{is}),
$$
$$
H_{BG,K,c}=\sum_rh_{c,r}+K^{-1}\sum_{r<s}W_{rs},
\quad
H_{\rm act}=C_KI+H_{BG,K,c}+H_{P,0}+\sum_iF_i^P\overline Z_{Gi}.
\tag{1}
$$
$\overline Z_{Gi}=K^{-1}\sum_rZ_{G,ir}$；
$F_i^P=\epsilon m_i^P+\nu t_i^P$、$\epsilon=3/20,\nu=7/20$沿原探针定义。程序是同一完整Hamiltonian的两个守恒控制块且保留；这里比较已知 $c$，不复制未知程序态。

每站初始
$$
\rho_A=\operatorname{diag}(3/4,1/4),\quad
\rho_B=\operatorname{diag}(3/5,2/5),\quad G=+X.
\tag{2}
$$
两热边缘用有限伙伴纯化且保留；不同簇初始独立，与任意未知探针/参考 $PR$独立。两程序的普通材料、纯化和初始能源相同；控制值不同，不称完整输入相同。纯度、复制、触发、记录均非免费。

$V=65^3,\ell=16/65$，有 $\|h_c\|\le4V,\|W\|\le2V$。原探针全Fock范数给
$$
B_P=V[8/5+18/\ell],\quad
L_P=V[2\epsilon(4/5)+18\nu/\ell],\quad
C_K=1+4VK+V(K-1)+B_P+L_P .
\tag{3}
$$
两程序的完整父Hamiltonian均 $\ge I$，比较过程保同一绝对相位。纯化及程序的简并能源属于声明材料。巨大有限资源和正谱平移不等于可行设备或GR真空源。

式(1)保源、关系、探针反作用和交互能源，自主总能精确守恒。各簇每站的 $n_A+n_B$也守恒，故两裸能源变化相反；$c=0$时各 $n_A,n_B$分别守恒。交互初均能为零不意味着生成元或后来存能可删。准备/终端仪器的控制能源和记录另由有限实现合同承担，不以主过程守恒免费取得。

## 2. 自洽热分配与关系响应

Hartree簇态保完整纯化；有效Hamiltonian分成 $AB$和 $G$两部分，所以指定准备闭合成六变量。单激发块的未归一Bloch向量为 $r=(3/20)u$，且
$$
u(0)=(0,0,-1),\quad g(0)=(1,0,0),\qquad
u'=2(c,0,g_z)\times u,\quad
g'=2(-1,0,(3/20)u_z)\times g.
\tag{4}
$$
两项来自同一交互能源的互易变分。$\lvert u\rvert=\lvert g\rvert=1$，每站
$$
{\cal E}=\Delta(13/20)+c(3/20)u_x-g_x+(3/20)u_zg_z
\tag{5}
$$
守恒，互易存能只计一次。两个程序都有初始热差驱动，必须比较开关差，不把任一非零 $g_z$都归于交换。

$c=0$时 $u_z=-1$恒定；$c=1$时
$$
\Delta E_A^{\rm mf}=\frac{3}{40}[u_{z,1}(T)+1],
\quad\Delta E_B^{\rm mf}=-\Delta E_A^{\rm mf},\quad
g_{z,1}(t)-g_{z,0}(t)=-t^4/10+O(t^6).
\tag{6}
$$
有限结论由下节完整尾界承担。

初片 $N(0)=1$；式(2)相对于已定裸能源 $\Delta n$的参数为
$\beta_A=\log3,\beta_B=\log(3/2)$。后继实际质量钟改变，不能据此声称所有源与交互生成元按该钟同率，或将这两个准备参数当全过程同观察者视界温标。

## 3. 固定短窗的精确有限证书

事前固定 $T=1/8$，不扫描。式(4)各分量的绝对级数由 $y'=4y^2,y(0)=1$控制：每个二次项的系数绝对和不超过4，线性项亦由 $y^2$系数控制。逐阶归纳 $\lvert a_n\rvert\le4^n$，收敛主级数保证所得解满足原ODE；唯一性识别它与真实实解。64阶截断在全窗满足
$$
\left|y_j(t)-\sum_{n=0}^{64}a_{j,n}t^n\right|\le2^{-64},
\quad
\left|\int_0^T y_j-\int_0^T\sum_{n=0}^{64}a_{j,n}t^n\right|\le2^{-67}.
\tag{7}
$$
$y_j$指六变量任一分量，非另设物理场。记 $I_c=\int_0^Tg_{z,c}dt$；有理递推、逐阶能源/长度恒等式及尾给
$$
I_0=0.00019468936427780456\ldots,\quad
I_1=0.00019408102898131251\ldots,
$$
$$
I_1-I_0<-5\times10^{-7},\qquad
\Delta E_A^{\rm mf}>21/10000.
\tag{8}
$$
展示小数不用于判定；[结果](results.json)保存向外包络。全窗结论由解析主级数承担，不由浮点采样外推。

## 4. 完整有限父过程与实际探针

复用[旧有界有限复制证明](../../archive_1161_1200/_admission/after1174_direction/mean_field_bound.md)式(10)、(14)，本轮重新代入 $w=2V,L=L_P$：
$$
\alpha(T)\le e^{12wT}/K,\qquad
\eta_P\le4L_PT e^{6wT}/\sqrt K<2^{-76102}<10^{-12}.
\tag{9}
$$
整数预算只用 $4L_PT<2^{22}$、$6wT=823875/2$、$e<4$及
$\sqrt K=10^{300000}>2^{900000}$，不构造巨大父矩阵。

比较背景始终是**精确**新模型背景联合演化，非整体Hartree产品；实际探针反作用保留。式(9)为指定背景准备上的完整等距误差，对原48维未知探针和任意 $R$统一成立。计数界另给任一簇局部态半迹距 $\le\sqrt\alpha$；加探针反作用误差仍 $<10^{-12}$。因此下述准备误差下真实局部裸能源转移有
$$
\Delta E_{A,c=1}^{\rm act}>1/500,\qquad
\Delta E_B^{\rm act}=-\Delta E_A^{\rm act}.
\tag{10}
$$
这是有界局部能源期望，不将小态误差扩为广延总能源误差。

新增有限准备、正窗触发、终端仪器的每支统一误差合同为
$\eta_{\rm inst}=10^{-10}$，保全部资源、条件后态和记录；其实现须另记完整控制能源。该精度是公开条件输入，未由旧 $2\times10^{-5}$预算或本轮代码制造。若实际误差更粗，须重算签收。式(10)允许初末局部能源各付一次 $10^{-10}$；完整总能守恒仍是自主Hamiltonian的算符事实。

## 5. 原有限格点的独立钟与位置读数

同一Hartree比较器为
$$
h_{m,p}(t)=N(t)m\beta+v(t)\boldsymbol\alpha\cdot\boldsymbol\ell_p,
\quad N=1+\epsilon g_z,\quad v=1+\nu g_z.
\tag{11}
$$
二者恒正。沿[1175原菜单](../../archive_1161_1200/1175/proof.md)及[1179有限格点反演](../../archive_1161_1200/1179/proof.md)，取两次已知见证准备：

- $p=0$、正 $\beta$支、质量 $3/10,1/2$的 $|+\rangle$，原 $P_{-Y}$钟效果。
- 质量零、空间 $(|0\rangle+|p_x\rangle)/\sqrt2$、$\alpha_x=+1$自旋，原完整位置记录与汇总 $f_x(j)=[1+\sin(2\pi j_x/65)]/2$。

同规则、同已知背景的独立试次，不在同个未知样本中途加测量；统一未知码证书仍为式(9)。精确格点给
$$
\tau_c=T+\epsilon I_c,\quad \xi_c=T+\nu I_c,\quad
k_\ell=\sin(2\pi/65)/\ell\in(3/10,2/5),
$$
$$
p_{D,c}=[1+\sin(\tau_c/5)]/2,\qquad
p_{x,c}=1/2+\sin(k_\ell\xi_c)/4.
\tag{12}
$$
不将 $k_\ell$换成连续动量；无连续格点误差费用。位置是既定有限采样菜单，不是连续精确位置仪器。

原单调支上独立反演
$$
\tau(p_D)=5\arcsin(2p_D-1),\qquad
\xi(p_x)=k_\ell^{-1}\arcsin(4p_x-2).
\tag{13}
$$
每项实际概率误差 $\delta=10^{-10}+10^{-12}$。有理三角包络核包括误差后的
$0<2p_D-1<29/1000$、$0<4p_x-2<68/1000$；逆函数Lipschitz常数分别可取11、14。于是
$$
[(\xi-\tau)_1-(\xi-\tau)_0]_{\rm ideal}
=(\nu-\epsilon)(I_1-I_0)<-10^{-7},
$$
$$
[(\xi-\tau)_1-(\xi-\tau)_0]_{\rm actual\ probabilities}
<-9\times10^{-8}.
\tag{14}
$$
这里是真概率证书；经验频率误差另经式(13)运输。

仅式(14)一般不足排除重定时。对本固定基准 $c=0$全窗，
$$
0<d\xi_0/d\tau_0=(1+\nu g_{z,0})/(1+\epsilon g_{z,0})\le M_*:=27/17.
$$
$\tau_1(T)<\tau_0(T)$，唯一 $0<t_*<T$使 $\tau_0(t_*)=\tau_1(T)$。斜率界给
$$
\xi_1(T)-\xi_0(t_*)\le(\nu-M_*\epsilon)(I_1-I_0)
=\frac{19}{170}(I_1-I_0)<0.
\tag{15}
$$
等价的两端诊断
$$
{\cal C}=(\xi_0-\xi_1)-M_*(\tau_0-\tau_1)
>\frac{19}{170}\,5\times10^{-7}
$$
的两支误差不超过 $2(14+11M_*)\delta$，故实际端口有
$\mathcal C>4\times10^{-8}$，且 $\tau_0-\tau_1>0$。它排除沿这份基准钟—位置曲线的共同重定时；$t_*$仅为证明比较点，未授实际中途停机或读口。不声称排除了任意换空间字典、材料或GR规范。

## 6. 热账与退出线

已知初始Gibbs边缘给任意实际后态的精确身份
$$
\Delta S_A=\beta_A\Delta E_A-D(\rho_A^{\rm out}\Vert\rho_A^{\rm in}),\quad
\Delta S_B=\beta_B\Delta E_B-D(\rho_B^{\rm out}\Vert\rho_B^{\rm in}).
\tag{16}
$$
后态不必Gibbs；实际初态若有误差，应使用一般
$\Delta S=\beta\Delta E-\Delta D$，不能忽略初相对熵。纯化、$ABG$关联、程序、探针及记录留在完整后态。式(5)保关系自能与交互存能；两个局部熵不替代全账，保能不推出零熵产。

连续低动量字典 $a=N/v$可沿旧一般Dirac作用识别；本轮严格概率只用有限格点，不认证新连续误差或完整Hilbert $T_{\mu\nu}$。同一声明作用下，裸热分配交换改变后续原生钟/传播，并超出式(15)的基准共同重定时。这是有限同父正向连接，耦合、图、准备和精度仍为输入。

不继承1272钟补偿/源抵销，不宣称Einstein方程、面积熵、视界温标、Clausius或稳定参考的受扰恢复。取得该有限连接即停止扩建。
[check.py](check.py)只计1个作者诊断组：有理ODE及尾、完整有限能量块、资源指数预算和端口反演。默认只读，首次 --write 排他生成[results.json](results.json)。

