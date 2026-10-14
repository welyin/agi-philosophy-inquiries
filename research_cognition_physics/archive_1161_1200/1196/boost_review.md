# 有限 boost 见证的独立误差审阅

2026-10-14。1196内部非盲独立核验；本独审不额外增加科学核验组。只审固定有效参考候选与有限占据读数，不领取完整1178父过程的参考态、实际仪器或GR交付。

## 1. 对象与归一化

取原声明参数 $z=0.989835554078098$、$N=1+0.15z$、$v=1+0.35z$、$\ell=16/65$、$m=2$，得到 $\mu=mN\ell/v\simeq0.41992380993825135$。横动量为零时，一个等价的两分量CAR块为

$$
h(k)=\frac v\ell[\mu\sigma_z+\sin k\,\sigma_x],\qquad
P_-(k)=\frac12\left[I-\frac{\mu\sigma_z+\sin k\,\sigma_x}{\sqrt{\mu^2+\sin^2k}}\right].
$$

实际有限参考来自65周期动量的负谱投影。取前32站的压缩协方差 $C_A$，没有将区域本身当成半无限空间。在完整三维对象中，这个读口选取整个横向周期截面的零横动量模式，不是三维小球或单点局部仪器。独审另对完整130维实空间周期矩阵取负谱投影，核对Fourier构造；不导入作者代码。

以站点位置 $x_i=(i+1/2)\ell$、键位置 $x_{i,i+1}=(i+1)\ell$ 加权原有效能源，并除以 $v$，候选快度boost为

$$
(B_L)_{ii}=\mu(i+1/2)\sigma_z,\qquad
(B_L)_{i,i+1}=-\frac{\mathrm i(i+1)}2\sigma_x,\quad i=0,\ldots,L-2.
$$

右式分子第一个 $i$ 是虚数单位，括号中的 $i+1$ 是站点指标加一。也即 $B=\sum xh/v$ 的所声明格点离散，不是把参考态先定义为该 $B$ 的Gibbs态。测试函数为 $F(t)=(1+e^{2\pi t})^{-1}$；$2\pi$ 在本页是待检验的几何归一，尚未由实际加速钟推出。

## 2. 远端截断的统一解析界

令 $\psi\in\mathbb C^{2L}$ 是归一化测试模，$b=\langle\psi,B_L\psi\rangle$。将其补零到任意更长的同系数开链 $B_M$，$M>L$。由于只有一条新增跨界键，

$$
r_M^2:=\|(B_M-b)\psi\|^2
=\|(B_L-b)\psi\|^2+\frac{L^2}{4}\|\psi_{L-1}\|^2,
\qquad \langle\psi,(B_M-b)\psi\rangle=0.
$$

其中 $\psi_{L-1}$ 是末站的两分量向量。对精确本征模，第一项为零；代码保留Rayleigh残差，不把浮点本征向量声称为精确本征向量。

标量函数有 $\sup_t|F''(t)|=2\pi^2/(3\sqrt3)$。对 $B_M$ 的谱测度应用带余项的Taylor公式，线性项期望消失，因而

$$
\left|\langle\psi,F(B_M)\psi\rangle-F(b)\right|
\le\frac{\pi^2}{3\sqrt3}\,r_M^2.
$$

该界对所有有限 $M>L$ 一致，完全不依赖逐次加长恰好数值收敛。本页不另声称已证明无限无界boost的自伴性或实际半空间实现。若另给包含这些有限支撑向量的自伴半链实现，同一谱测度证明亦适用。

## 3. 65周期参考与无限链局部参考的alias界

本节的无限链仅为辅助比较对象，不是换掉原65站准备。设 $0<\eta<\operatorname{arsinh}\mu$。复动量条带内

$$
\operatorname{Re}(\mu^2+\sin^2(x+i\eta))
\ge\mu^2-\sinh^2\eta>0,
$$

故负谱符号的解析延拓无分支零点，并可取

$$
M_\eta=\frac12\left[1+
\frac{\mu+\cosh\eta}{\sqrt{\mu^2-\sinh^2\eta}}\right]
$$

作为符号算符范数上界。移位Fourier积分给 $\|C_\infty(r)\|\le M_\eta e^{-\eta|r|}$。绝对可和性允许使用精确周期alias展开 $C_n(r)=\sum_{q\in\mathbb Z}C_\infty(r+qn)$。对于 $|r|<n$，

$$
\|C_n(r)-C_\infty(r)\|
\le E_\eta(r):=
M_\eta\frac{e^{-\eta(n-|r|)}+e^{-\eta(n+|r|)}}{1-e^{-\eta n}}.
$$

令 $a_j=\|\psi_j\|$，占据期望的参考替换误差因此至多

$$
\epsilon_{\rm alias}
=\sum_{i,j=0}^{L-1}a_i E_\eta(i-j)a_j.
$$

采用事前固定 $\eta=0.4$，不靠逐模拟合参考。周期负谱与辅助无限参考的量子占据差受该界约束；它不是对全部无限Fock态的全局迹距离声明。

## 4. 靠近一个切口的读口范围

保留前 $s$ 站，将测试模截断后重新归一为 $\phi_s$，舍弃权重为 $q_s$。两个单模投影的非零特征值为 $\pm\sqrt{q_s}$。对任意 $0\le C\le I$，

$$
|\langle\psi,C\psi\rangle-\langle\phi_s,C\phi_s\rangle|
\le\sqrt{q_s}.
$$

这给占据概率读口的局部化代价。重新归一的前缀模是另一个明确读口，不能声称它与原模完全相同；比较两个模型的差时可保守加两份此界。仅这些effect层面的估计，不完成其实际仪器、后态或后续应力插入的制造与误差认证。

## 5. 独立数值与结论

脚本 [reviewer_boost.py](reviewer_boost.py) 默认只将JSON写到stdout，不保存或覆盖作者results。最低正boost模的独立浮点结果为：

- $b\simeq0.356518434051165$；
- $\langle C_A\rangle\simeq0.104325855286873$；
- $F(b)\simeq0.096209334248420$；
- 差约 $0.008116521038453$；
- 末站权约 $2.8974\times10^{-13}$，统一远端boost截断界约 $1.4089\times10^{-10}$；
- 固定 $\eta=0.4$ 的周期参考alias界约 $7.46\times10^{-10}$；
- 前16站以外的模权约 $9.95\times10^{-8}$。

实际JSON值以默认复算输出为准。代码另核对48、64、96站辅助boost的同一个补零读口，作为解析界的部署检查，不能以数值序列代替上面的证明。固定点、输入小数及矩阵函数的有理区间认证由主线另做；本脚本是独立浮点核验，未将输出小数伪装为向外舍入区间。

**可发布结论：** 在这份固定格点有效参考上，最低正boost模的占据与事前规定的 $2\pi B$ Fermi预测存在有限差异；上述周期参考和远端boost截断效应远不足以解释该差。它不排除格点修正后的模生成元，也不排除有限有效Lorentz/Unruh窗口。测试模的低boost谱值不等同于原物理Hamiltonian的低动量支持，不能据其名称直接领取连续近似。

本页没有证明这份Hartree有效参考就是完整有限K自治父过程的真基态，亦未计入其准备误差；完整父过程若不能提供此参考，结果只保留为候选诊断。区域、参考及几何boost目前是明确声明的数学对象，并非已经认证的全部原仪器。不得将这个有限候选的差升级为Unruh效应、GR或认知纲领的失败。

## 6. 复算

使用既有Python与NumPy，以不写字节码方式运行：

```powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'D:\workspace\AGI的哲学思考\research_cognition_physics\archive_1161_\_admission\after1195_reference_candidate\reviewer_boost.py'
```

原量子概率与辅助边界估计分开输出；没有对 $C_A$ 求对数或裁剪其谱。
