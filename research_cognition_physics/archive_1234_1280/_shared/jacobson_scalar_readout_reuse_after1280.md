# 1280之后的标量有限读口：成熟结论接续

2026-10-17。非编号接续，不新增科学轮次、核验组、公理或器件任务。本页给[1280](../1280/proof.md)同一最小无质量标量的一份可选有限校准数学接口；不改其冻结证明。它承接[观察者范围筛查](jacobson_observer_completion_screen_after1280.md)，不把装置自身的完整Einstein解、微观制造链或参考制备工厂设为本阶段必经门槛。真正调用时仍须登记会影响共同预测的宏观有效来源、读数与误差合同。

## 1. 同一参考、同一平坦段及明确紧支撑

1280在$u>-1$有$r=1,\eta=u$，自然参考$\omega_c$在该域恰为Minkowski真空，工作相干均值$f$及其导数为零。因此该域内任何局域场读口的工作态限制与参考限制相同；这不删除原脉冲或跨域记录。

以下取原加速世界线族的$\kappa=1$单位，写

$$
t=\rho\sinh\tau,\qquad x=\rho\cosh\tau,\qquad
u=-\rho e^{-\tau}/\sqrt2,\qquad v=\rho e^\tau/\sqrt2.
\tag{1}
$$

中心$\rho=1,y=z=0$的$\tau$就是原观察者固有时；管内其它$\rho$的lapse为$\rho$，不称整个空间包共用各自固有时。一般$\kappa$下，本文$\beta=2\pi,4\pi$是无量纲boost参数的逆温度，对应中心固有时间的$2\pi/\kappa,4\pi/\kappa$。

令

$$
b(s)=\begin{cases}\exp[-1/(1-s^2)],&|s|<1,\\0,&|s|\ge1,\end{cases}
\qquad I_b=\int_{-1}^{1}b(s)\,ds>0,
$$

$$
h(\tau,\rho,y,z)=\frac{128}{I_b^4}
 b(2\tau-5)b(4(\rho-1))b(4y)b(4z).
\tag{2}
$$

体积元为$dV=\rho\,d\tau\,d\rho\,dy\,dz$，故$\int h\,dV=1$。这是一份明确实、非负、非零的$C_c^\infty$测试，支撑于

$$
2\le\tau\le3,\quad \frac34\le\rho\le\frac54,
\quad |y|,|z|\le\frac14.
\tag{3}
$$

其最小$u$不小于$-5e^{-2}/(4\sqrt2)>-1$，且$v>0$；所以整个有限管位于1280的平坦段。不使用点耦合或尖锐开关。一般“实非零$h$”不足保证可辨，因为$h=\Box g$可以给$\Phi(h)=0$；下面对(2)的正性直接证明排除此例外。

## 2. 两份局域KMS参考的严格协方差差

比较同一平直楔CCR上的标准零均值准自由$\beta$-KMS参考。$\beta=2\pi$恰等于上述$\omega_c$的局部限制；$\beta=4\pi$只作为另一局域校准假设，未赋予相同绝对应力或同一自洽几何。两者在远离视界的紧域均具有Hadamard短程性质，$V_\beta:=W_\beta(h,h)$有限。

实际读取的原始来源是[Moretti–Vanzo, hep-th/9507139v2](https://arxiv.org/pdf/hep-th/9507139)，后发表于Phys. Lett. B375 (1996) 54–59：式(2)给标量热二点函数，式(5)给$\beta=2\pi$与Minkowski核的身份。由这两个明确核相减，令

$$
\cosh\alpha=\frac{\rho^2+\rho'^2+|\mathbf y-\mathbf y'|^2}{2\rho\rho'},
\quad A=\cosh(\alpha/2),\quad B=\cosh((\tau-\tau')/2),
$$

得到光滑实函数

$$
D(X,X'):=W_{2\pi}(X,X')-W_{4\pi}(X,X')
=\frac1{32\pi^2\rho\rho'A(A+B)}>0.
\tag{4}
$$

这是分布奇性相消后的恒等式，不是在类光奇点处直接取绝对值。在(3)上，$\cosh\alpha\le29/9$，所以$A\le\sqrt{19}/3<3/2$；又$B\le\cosh(1/2)<8/7$，最后一个界由$(2n)!\ge2^n$的幂级数上界得到。结合$\rho\rho'\le25/16$及(2)归一，有

$$
\boxed{V_{2\pi}-V_{4\pi}
=\iint h(X)h(X')D(X,X')\,dVdV'
>\frac{14}{2775\pi^2}>0.}
\tag{5}
$$

这已经是固定有限测试的严格差，不需要无限时跃迁率、频谱扫描或新的数值模拟。

## 3. 精确Ramsey响应及有限场注能

作为公开有效读口原语，取探针Pauli算符$Z$在交互绘景恒定，初态$|+\rangle$，并采用

$$
H_I(t)=\lambda Z\int d^3\mathbf x\,h(t,\mathbf x)\Phi(t,\mathbf x),
\qquad \lambda\ne0.
\tag{6}
$$

这要求自由探针Hamiltonian与$Z$对易；前后旋转、加速支撑、开关实现及末记录并未由此自动提供。因自由场CCR为c-number，时间序列指数的二阶交换子仅贡献与探针支无关的相位，高阶交换子为零。因此

$$
U=e^{ic_h\lambda^2}\exp[-i\lambda Z\Phi(h)],
\qquad
p_\beta(+)=\frac{1+\exp[-2\lambda^2V_\beta]}2.
\tag{7}
$$

这是精确Gaussian读口，不依赖无界标量UDW的截断Dyson展开。成熟条件位移/Ramsey接口可见[Tian–Jing–Du (2023), §3.1式(5)—(7)、§3.2式(9)](https://cdn.sciengine.com/doi/pdf/4EADB91F5975453AB235D80F94FD2907)；本页使用(2)的平滑紧支源，不搬用该文环面拓扑或脉冲列。

由(5)，固定$\lambda$给

$$
d:=p_{4\pi}(+)-p_{2\pi}(+)
=\frac12e^{-2\lambda^2V_{2\pi}}
 \bigl(e^{2\lambda^2(V_{2\pi}-V_{4\pi})}-1\bigr)>0.
\tag{8}
$$

它测量参考协方差，不是自动实现非零能隙跃迁的详细平衡测温；在明确两假设菜单内已经有严格有限概率差。

在真实$\omega_c$输入上，令$F_{\rm ret}$为(6)单位源产生的零初值迟滞经典位移（波动方程源符号由(6)固定，不影响以下平方能源）。开关结束后的非选择场边缘是正、负位移的等权混合，故

$$
\Delta E_\Phi=\lambda^2 E[F_{\rm ret}],\qquad
E[F]=\frac12\int_\Sigma\bigl((\partial_tF)^2+|\nabla F|^2\bigr)d^3x,
$$

$$
\boxed{\Delta E_\Phi\le
\frac{\lambda^2}{2}
 \left(\int dt\,\|h(t,\cdot)\|_{L^2(d^3x)}\right)^2<\infty.}
\tag{9}
$$

证明只需波动方程能源恒等式及$|\dot E|\le\|h_t\|_2\sqrt{2E}$；(2)使右侧为明确有限常数。源作用结束时交互项为零，控制来源必须支付这笔场注能。式(9)是有限场转移账，不是整个设备能源上界，也不是任意选择分支的同一上界。

## 4. 另一温度参考不是零来源替代

这里始终采用1280原最小Hilbert应力，不能直接将Moretti–Vanzo式(6)的共形应力式搬给它。可完全由(4)作最小点分裂检查：令$F=\langle\Phi^2\rangle_{4\pi}-\langle\Phi^2\rangle_{2\pi}$，同Minkowski减除下

$$
F=-\frac1{64\pi^2\rho^2},\qquad
\Box F=-\frac1{16\pi^2\rho^4}.
$$

记$D_{ab'}=\nabla_a\nabla_{b'}(W_{4\pi}-W_{2\pi})$。两点核满足质量零方程，故$g^{ab}D_{ab'}|=\Box F/2$。在静止正交时间方向$\hat0=\rho^{-1}\partial_\tau$，(4)给$D_{\hat0\hat0'}|=-1/(512\pi^2\rho^4)$。于是

$$
\boxed{\langle T^{\rm min}_{\hat0\hat0}\rangle_{4\pi}
-\langle T^{\rm min}_{\hat0\hat0}\rangle_{2\pi}
=D_{\hat0\hat0'}|+\frac14\Box F
=-\frac9{512\pi^2\rho^4}\ne0.}
\tag{10}
$$

该差在指定管内有限，且$\beta=2\pi$的绝对参考应力已由1280固定为零。它不证明$\beta=4\pi$能作为相同几何的自洽全局参考；本文只比较其楔内局域限制，不领取跨视界延拓或全球有限总能。

## 5. 重复、来源与停止线

如果另有独立重复输入及相同读口的合同，对(8)两概率使用阈值检验，Hoeffding界给任何$0<\epsilon<1$下有限整数

$$
N\ge\frac{2}{d^2}\log\frac2\epsilon
\tag{11}
$$

足以使错误概率不超过$\epsilon$。单个场受扰后的连续重复不自动独立；式(11)没有提供重置、样本供应、记录维护或实际精度。一轮数学比较已给(8)，无需为本页实现多次实验。

- [1232](../../archive_1219_1233/1232/proof.md)的有限响应、完整来源与误差分账可直接复用；其CAR数值和有界Dyson条件不搬给标量。本页只以线性CCR给同1280材料的精确替代接口。
- 有限管位于$v>0$，因此该管新产生的迟滞扰动不抵达旧$v=0$ cap；这不删除装置已有质量、支撑、控制或记录的来源，也不单凭因果先后证明完整几何误差为零。
- 完整账可采用明确的宏观有效来源及误差合同，不要求先解装置自身Einstein方程。然而(6)的给定源、Ramsey旋转和末读并不自动满足这样的合同；不得把它们当免费控制或把未核预算赋零。
- 本页不改变1280参考，不认证第二参考的同几何实现，不证明实际仪器、热库交换、噪声／反作用共同预算或全目标完成。到此停止器件扩建；只保存成熟结果的有限接续。
