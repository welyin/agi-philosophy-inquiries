# 1243证明：有限视界预算到无迹残差及窗口迹闭合

2026-10-15。取$\hbar=c=k_B=1$、号差$-+++$。这是一条有效描述中的条件定理，不以Einstein方程为输入，不采用CT。单位和标架均须由同一操作字典标定。

## 1. 同对象合同

在包含全部比较窗口和连接管的相对紧坐标域$D$上，给定$C^3$ Lorentz度规$g$、$C^1$对称完整来源$T$、光滑正交标架$e_a$，以及常数$s>0$。记$\kappa=2\pi/s$。有效正则性不表示微观必须连续，也不是重新加入三维空间证明已经消去的假设。

在有限覆盖点$p_j$，使用304的九个未来零方向：$k_{\pm i}=e_0\pm e_i$及$k_{ij}=e_0+(e_i+e_j)/\sqrt2$。每个方向对应无旋、无焦散的有限仿射束$\lambda\in[-\ell,0]$，其终截面有正面积$A_0$。所有几何、方向、面积、温标、热流及来源运输均属于同一候选历史；不把独立实验的不同背景拼在一起。

归一末端面积测度记为$d\mu_0$，沿各生成线的面积Jacobian为$F(\lambda,y)$，$F(0,y)=1$。定义

$$
r=R_{ab}k^ak^b,\qquad f=\kappa T_{ab}k^ak^b,\qquad
H=\int d\mu_0\int_{-\ell}^0(-\lambda)Ff\,d\lambda.
$$

公开的有限热—几何合同为

$$
\left|\frac{\Delta A}{A_0}-H\right|\le\delta.
\tag{1}
$$

$\delta$是面积单位的总误差：实际Clausius偏差、面积读出、熵—面积身份偏差、同观察者温标/boost近似、完整热源和运输误差，须经三角界共同进入；已明确单列的几何修正不再重复计入。若原熵单位误差为$\delta_S$，对应项为$\delta_S/(sA_0)$；若$A_0$不精确，须用已认证的正下界及分母误差。不能将$A=S/s$作为面积定义来取得(1)。本轮不证明仪器实现这些精度，也不证明(1)由TR或原公理推出。

每束给出独立统一界（$x=-\lambda\in[0,\ell]$）：

$$
|r-r(p_j,k_A)|\le d_r+L_rx,\quad
|f-f(p_j,k_A)|\le d_f+L_fx,\quad |r|\le R_*,\quad |f|\le F_*,
$$

$$
|\theta|\le b_0+Bx,\qquad
\sigma_{ab}\sigma^{ab}\le s_0^2+D_\sigma x^2.
\tag{2}
$$

$d_r,d_f$包括终截面横向大小和方向/标架运输的预算；剪切平方使用正定屏幕度规。$b_0$控制末端膨胀，$s_0$控制末端剪切。界可由候选的独立有效正则估计提供，不能从有限采样本身自动取得。

## 2. 有限窗口归一：不能混同两种权重

复用316的Raychaudhuri恒等式，改为本轮端点约定，逐生成线有

$$
\frac{\Delta A}{A_0}=\ell\langle\theta(0)\rangle_0+
\int d\mu_0\int_{-\ell}^0(\lambda+\ell)F
\left(r+\sigma^2-\frac{\theta^2}{2}\right)d\lambda.
\tag{3}
$$

这是$F'=\theta F$与$(F\theta)'=F(\theta^2/2-\sigma^2-r)$的一次分部积分。**面积核是$\lambda+\ell$，热核是$-\lambda$；只在主阶常值比较时给同一$\ell^2/2$。** 本证明保留有限差异。

令

$$
z=b_0\ell+B\ell^2/2,\quad M=e^z,\quad a=M-1,
$$

$$
N=s_0^2+b_0^2/2+b_0B\ell/3+(D_\sigma+B^2/2)\ell^2/6.
$$

由(2)有$F\le M,|F-1|\le a$。面积核的归一$x,x^2$矩分别为$\ell/3,\ell^2/6$；热核的归一$x$矩为$2\ell/3$。将$Fr$拆成$r(p_j)+(r-r(p_j))+(F-1)r$，对$f$同样处理，得到

$$
|(R-\kappa T)_{ab}k_A^ak_A^b(p_j)|\le\zeta_j,
$$

$$
\boxed{\zeta_j=
\frac{2\delta}{\ell^2}+\frac{2b_0}{\ell}+d_r+d_f+
\frac{L_r\ell}{3}+\frac{2L_f\ell}{3}
+a(R_*+F_*)+MN.}
\tag{4}
$$

不同束可有不同预算，取共同上界即可。不需要先设$R_{kk}=\kappa T_{kk}$，所以没有以场方程证明自身。

固定未归一误差$\delta>0$时，(4)中的噪声项为$2\delta/\ell^2$；无限缩窗不能制造精确方程。有限有效窗口须同时容纳读口噪声和截断误差。本轮不选定物理最小尺度。

## 3. 从有限覆盖到全部无迹分量

另给独立覆盖合同：九个标架方向的标量$q_A(x)=(R-\kappa T)(k_A,k_A)(x)$在有限覆盖半径$h$内满足共同模连续界$\omega(h)$。覆盖须包含随后连接场的整个支撑和两个窗口，而不是只覆盖中心事件。于是

$$
|q_A(x)|\le\zeta:=\max_j\zeta_j+\omega(h),\qquad x\in D.
\tag{5}
$$

这一步直接沿用1192的“有限覆盖＋独立正则”责任分配。若不提供它，采样外紧支光滑扰动即可逃过全部有限读数，不能使用下文管积分。

定义同一张量

$$
E_{ab}=G_{ab}-\kappa T_{ab}=S_{ab}-\Lambda(x)g_{ab},\qquad
S_{ab}=R_{ab}-\kappa T_{ab}-\frac{R-\kappa T^a{}_a}{4}g_{ab},
$$

$$
\Lambda(x)=\frac{R+\kappa T^a{}_a}{4},\qquad g^{ab}S_{ab}=0.
\tag{6}
$$

复用304/1097九方向反演。为固定所用范数，记$a_i=(q_{+i}+q_{-i})/2$、$b_i=(q_{+i}-q_{-i})/4$，则

$$
S_{00}=\frac14\sum_i a_i,\quad S_{ii}=a_i-S_{00},\quad
S_{0i}=b_i,\quad S_{ij}=q_{ij}-\frac{a_i+a_j}{2}-\sqrt2(b_i+b_j).
$$

在正定辅助标架范数$|S|_e^2=\sum_{ab}|S_{ab}|^2$下，逐项界为$3\zeta/4,5\zeta/4,\zeta/2,(2+\sqrt2)\zeta$。用$2+\sqrt2<7/2$得

$$
|S|_e^2\le\frac{321}{4}\zeta^2<81\zeta^2,
\qquad \varepsilon_S:=9\zeta.
\tag{7}
$$

没有定义Lorentz“正定长度”，也没有把旧九方向代数计为本轮新定理。

## 4. 弱Bianchi控制窗口间的迹差

令$j_b=\nabla^aT_{ab}$，包含所有材料、准备、控制、记录与交互。缩并Bianchi恒等式与(6)给

$$
\partial_b\Lambda=\nabla^aS_{ab}+\kappa j_b.
\tag{8}
$$

**不对(7)微分。** 取两个非负、单位坐标积分窗口$\chi_a,\chi_b$，$\rho=\sqrt{|\det g|}$，以及紧支向量场$V$使

$$
\operatorname{div}_gV=\frac{\chi_a-\chi_b}{\rho},\quad
\Lambda_i=\int\Lambda\chi_i\,d^4x,\quad
W_1=\int|\nabla^aV^b|_e\,d\mathrm{vol}_g.
$$

这是用$\chi_i/\rho$作为度规体积窗口；平均不随意漏掉体积权。分部积分给精确身份

$$
\Lambda_a-\Lambda_b=
\int S_{ab}\nabla^aV^b\,d\mathrm{vol}_g-\kappa\langle j,V\rangle.
\tag{9}
$$

若完整来源弱守恒预算为$|\langle j,V\rangle|\le\varepsilon_j(V)$，则

$$
\boxed{|\Lambda_a-\Lambda_b|\le9\zeta W_1+\kappa\varepsilon_j(V).}
\tag{10}
$$

$W_1$由同一候选几何和选定测试场计算/界定，不是免费选择的小数。若以$\Lambda_a$为待定共同常数，则窗口$b$中、按同一标架平均的完整方程残差满足

$$
\left|\int [G_{ab}+\Lambda_a g_{ab}-\kappa T_{ab}]_e\chi_b\,d^4x\right|_F
\le9\zeta+2[9\zeta W_1+\kappa\varepsilon_j(V)].
\tag{11}
$$

这只控制所声明窗口及比较菜单，不给点态$\nabla\Lambda$界，不自动给初值问题解的距离。共同常数值仍未被选择。

## 5. 一个显式紧支连接与严格有限预算

令$q_h(x)=15(1-x^2/h^2)^2/(16h)$（$|x|\le h$，外为0），$Q=\prod_{i=1}^3q_h(x_i)$，$F_h$为$q_h$的CDF。两窗沿$x^0$相隔$d=b-a>0$：

$$
\chi_a=q_h(x^0-a)Q,\quad \chi_b=q_h(x^0-b)Q,\quad
U=Q[F_h(x^0-a)-F_h(x^0-b)],\quad V^\mu=\delta^\mu_0U/\rho.
$$

$U$紧支于两窗连接管，且$\partial_0U=\chi_a-\chi_b$。由非负性、单位积分和单峰性，

$$
\|U\|_1=d,\qquad\sum_\mu\|\partial_\mu U\|_1
\le2+\frac{45d}{8h}.
\tag{12}
$$

该$C^1$紧支场足够用于弱式。若要求光滑，可在扩大的同一管中共同卷积，导数$L^1$界不增加，并把覆盖扩至新支撑。一般几何下代入$V=U/\rho$，保留$\partial\log\rho$、联络和标架变换。平直标架的校准例中$W_1\le2+45d/(8h)$；此例不作为几何生成模型。

取无量纲预算

$$
\ell=1/10,\ \delta=10^{-6},\ b_0=10^{-6},\ B=1/100,
\ s_0=0,\ D_\sigma=10^{-4},\ L_r=L_f=10^{-3},
$$

$$
d_r=d_f=10^{-5},\ R_*=F_*=1/100,\ \omega(h_{\rm cover})=10^{-5}.
$$

代码用$e^z\le1/(1-z)$（$0\le z<1$）得到完全有理的保守上界，不以浮点指数充当证书。取测试窗$h=1/4,d=1/2$、$W_1\le53/4$、$\kappa\varepsilon_j=10^{-5}$，得到(10)的严格上界小于$0.042$。因此若同一有限窗口比较另有可靠下界$0.05$，这些共同合同不相容，余量超过$0.008$。**这是条件判据，不是已经测得该迹差，亦不证明该预算可由同父仪器达到。** $h$在此为测试窗宽，与覆盖半径有别。

## 6. 删除见证及剩余自由

1. **删除弱守恒：** 在平直几何令$T_{ab}=L(x)\eta_{ab}/\kappa$，则$S=0$，全部零方向残差为零，却有$\Lambda=L(x)$和$j=dL/\kappa$。任意迹变动不能由零方向数据单独排除。此为张量接口反例，不声称其实现完整认知材料动力学。
2. **将小残差直接微分：** 平直标架取$S_{00}=3f(t), S_{ii}=f(t),\Lambda=-3f(t)$，$T=-(S-\Lambda\eta)/\kappa$，则$j=0$。$f=\alpha\sin Nt$时无迹误差与$N$无关而导数随$N$增大；弱式(9)始终成立。不给材料方程或新增物理候选权限。
3. **删除覆盖：** 有限中心点之外的开集内可加平滑紧支扰动。旧1192已给同类见证，本轮只复用，不另制造新模型。
4. **固定噪声无限缩窗：** 相同预算把$\ell$从$0.1$缩到$0.01$时，噪声项从$0.0002$变为$0.02$，(10)的相容上界反而变差。
5. **常数和耦合：** 加同一常数到$\Lambda$不会被(10)确定；$\kappa=2\pi/s$依赖尚未选定的$s$。额外场、高阶修正或不满足(1)的非平衡效应须计入余项，不能被本工具宣布消失。

## 7. 适用结论

新增的是有限热读数、几何聚焦、同对象运输、全无迹张量与弱迹接续的联合误差链。它能筛掉一类在有限窗口内不相容的候选数据；不能独立提供几何熵、热参考、实际Clausius机制、物质公理或共同几何反作用。CT继续是未采用猜想，三项GR整体交付仍待完成。

成熟起点：[Jacobson 1995](https://arxiv.org/abs/gr-qc/9504004)，式(1)—(6)；项目直接复用[304](../../archive_301_341/research_note_304.md)、[316](../../archive_301_341/research_note_316.md)、[1097证明§8](../../archive_1086_1105/1097/proof.md)、[1192](../../archive_1161_1200/research_note_1192.md)及[1233后接口](../../archive_1219_1233/_shared/gravity_reconstruction_interface_after1233.md)。
