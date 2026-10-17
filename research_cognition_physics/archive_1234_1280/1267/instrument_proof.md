# 1267统计仪器：固定粗化效果、有限环境与后态来源

2026-10-16。按[准入](admission.md)交付一个作者科学组。接[1251 §3](../1251/proof.md)原有限POVM和[1266](../1266/proof.md)同一真实过程；与[几何分项](geometry_proof.md)共享固定读口，最终地位以验收为准。

**主结论：原有限统计效果可选一份有限结果Lüders仪器，明确全部后态和参考，并控制原系统有界来源及记录加权来源。** 数学等距只需要两个有限标签寄存器。它不证明新增参与者的物理能源、压力、世界管或因果权限。第6节的测量—重备仪器只是另一份后态对照，不与主仪器混同。

## 1. 原效果和具名有限后态

七对系统为$\mathcal H_S=L^2(\mathbb R^7)$，$z\in\mathbb R^{14}$。采用原相干态规范：

$$
d\mu_\varepsilon(z)=\frac{d^{14}z}{(2\pi\varepsilon)^7},
\qquad \int |z\rangle\langle z|\,d\mu_\varepsilon(z)=I,
\qquad
|\langle z|w\rangle|^2=e^{-|z-w|^2/(2\varepsilon)}.
$$

相干分析映射$(W\psi)(z)=\langle z|\psi\rangle$是等距。对原事前固定分箱$c(z)$、$\Omega_c=\{z:c(z)=c\}$，令$P_c$为乘以$1_{\Omega_c}$的投影：

$$
E_c=W^\dagger P_cW,\quad 0\le E_c\le I,\quad \sum_{c=1}^D E_c=I.
$$

主仪器选为

$$
\boxed{\mathcal L_c(\rho)=\sqrt{E_c}\rho\sqrt{E_c},
\qquad \Lambda_L=\sum_{c=1}^D\mathcal L_c.}
\tag{1}
$$

它是正常CP仪器，效果**恰为**$E_c$，全部原分箱概率和重复统计不变。对任意正常联合输入$\rho_{SR}$，使用$\sqrt{E_c}\otimes I_R$；非选择参考边缘不变，但系统后态和关联一般改变。不能以来源误差小推断未知态保留或小diamond误差。

完整数学输出可取

$$
V_L\psi=\sum_{c=1}^D
\sqrt{E_c}\psi\otimes|c\rangle_C\otimes|c\rangle_{C'},
\qquad V_L^\dagger V_L=I.
\tag{2}
$$

迹去$C'$即为可见$cq$输出；若未来允许访问$C'$，须保留它。两个标签空间各维数$D$，不要求连续微记录环境。这里$W$只用来证明效果和误差，不被当作另一个实际连续记录器。原$\mathcal H_S$仍无限维，有限结果不等于整个材料Hilbert空间有限维。

式(1)是公开选择的测量后态，不是暗中替换准备。它接收1266同一生成元所得的真实节点态；统计分支此前没有在旗标终点插入另一投影。

### 1.1 三元分箱、固定a边缘与参数范围

在原$\kappa=3$两程序中，$c(z)$是1251对$(a_c,\rho_c,P_c)$、宽$10^{-10}$的三元分箱。几何分项的全$\kappa\in[3,3.001]$族只固定

$$
c_a(z)=q_{10^{-10}}(a_c(z)),\qquad a_c=\chi e^Q.
$$

$\rho_c$含$J_\kappa$，故不可说整个三元分箱跨$\kappa$固定。对$a$边缘可另取式(1)，其结果数$D_a\le D$。**细分与边缘Lüders仪器的非选择通道一般不同**，本稿只给两者同样适用的界，不将其后态等同。

第5节$<2^{-80}$总来源误差只签1266原$\kappa=3$两程序节点族。一般参数族使用自己的测前误差：

$$
\delta_{\rm post}(b)\le\delta_{\rm pre}(b)+q_L(b).
\tag{3}
$$

全参数几何误差$<2^{-40}$不自动成为全部来源$<2^{-80}$。反演固定几何参数不要求本终末仪器恢复输入或继续未声明的未来。

## 2. Fourier交织缺陷：硬分箱无需微分

以下先取实来源符号$b=b_\infty+\widetilde b$，其中$\widetilde b$属原紧支$C^{40}$类。傅里叶约定为

$$
F_j(b)=(2\pi)^{-14}\int|\xi|^j
|\widehat{\widetilde b}(\xi)|\,d^{14}\xi.
$$

常数部分精确保留，$j=0$时算符范数另计$|b_\infty|$。记$B=\operatorname{Op}_\varepsilon(b)$、$T_b=M_bW-WB$。对于字符$u_\xi(z)=e^{i\xi\cdot z}$及幺正$U_\xi=\operatorname{Op}_\varepsilon(u_\xi)$，原Gaussian恒等式给

$$
\operatorname{AW}_\varepsilon(u_\xi)
=e^{-\varepsilon|\xi|^2/4}U_\xi,
$$

$$
(M_{u_\xi}W-WU_\xi)^\dagger(M_{u_\xi}W-WU_\xi)
=2(1-e^{-\varepsilon|\xi|^2/4})I.
\tag{4}
$$

因此由$1-e^{-x}\le x$及Fourier算符积分，

$$
\boxed{\|T_b\|\le\sqrt{\varepsilon/2}\,F_1(b).}
\tag{5}
$$

$M_b$与所有$P_c$对易，故

$$
[E_c,B]=-W^\dagger P_cT_b+T_b^\dagger P_cW,
\qquad
\|[E_c,B]\|\le\sqrt{2\varepsilon}\,F_1(b)=:\delta_b.
\tag{6}
$$

式(4)—(6)对任意可测硬分箱成立；没有对$1_{\Omega_c}$求导，不引入边界平滑、最小谱值或另一分辨率假设。anti-Wick来源可写为$B=\operatorname{Op}_\varepsilon(\beta)$、$\beta=e^{\varepsilon\Delta/4}b$，其$F_j$不增，故同样成立。

## 3. 平方根交换子与非选择来源

对任意$0\le E\le I$和有界自伴$B$，令$\delta=\|[E,B]\|$。强算符积分

$$
\sqrt E=\frac1\pi\int_0^\infty
t^{-1/2}E(E+t)^{-1}\,dt
$$

的交换子在范数中收敛，因为

$$
\|[E(E+t)^{-1},B]\|
\le\min(2\|B\|,\delta/t).
$$

在$a=\delta/\|B\|$处分段，零值情形取极限，得到

$$
\|[\sqrt E,B]\|
\le\frac1\pi(4\|B\|\sqrt a+2\delta/\sqrt a)
=\frac6\pi\sqrt{\|B\|\delta}
<2\sqrt{\|B\|\delta}.
\tag{7}
$$

不要求$E$可逆。由$\sum_cE_c=I$，

$$
\Lambda_L^\dagger(B)-B
=\sum_c\sqrt{E_c}[B,\sqrt{E_c}],
$$

$$
\boxed{
\|\Lambda_L^\dagger(B)-B\|
\le 2D\sqrt{\|B\|\sqrt{2\varepsilon}F_1(b)}
=:q_L(b).
}
\tag{8}
$$

严格号仅用于非零的粗界，不用于零符号。对原菜单统一用

$$
\|B\|\le2^{E_0+1},\quad F_1(b)\le2^{E_0},
\qquad
q_L(b)<4D\varepsilon^{1/4}2^{E_0},\quad E_0=14008.
\tag{9}
$$

额外因子2覆盖原常数加紧支项。纯紧支时可得$3D$，此处不追求该常数。

算符界对任意正常输入、任意参考及有界参考测试均适用，源项与参考测试张量相乘时乘其范数。它控制原系统已声明的有界来源；不控制未完成$\psi^2$、场的全部无界应力，亦不以关系生成元代替物理能源。

## 4. 有限记录与来源的联合界

给定实权$|r_c|\le1$，令$r(z)=r_{c(z)}$。对Weyl来源，

$$
\sum_cr_c\mathcal L_c^\dagger(B)-\operatorname{AW}_\varepsilon(rb)
=
\sum_cr_c\bigl(\sqrt E_cB\sqrt E_c-E_cB\bigr)
-W^\dagger M_rT_b.
\tag{10}
$$

于是

$$
\left\|\sum_cr_c\mathcal L_c^\dagger(\operatorname{Op}_\varepsilon b)
-\operatorname{AW}_\varepsilon(rb)\right\|
\le q_L(b)+\sqrt{\varepsilon/2}F_1(b).
\tag{11}
$$

若后态来源是$\operatorname{AW}_\varepsilon(b)$而目标仍为$\operatorname{AW}_\varepsilon(rb)$，先用$\beta=e^{\varepsilon\Delta/4}b$，再用AW正单收缩性，需额外加

$$
\|\operatorname{AW}_\varepsilon(r(\beta-b))\|
\le\|\beta-b\|_\infty
\le\varepsilon F_2(b)/4.
\tag{12}
$$

所以对两种来源都可用统一界

$$
\boxed{q_{L,r}<8D\varepsilon^{1/4}2^{E_0}}
\qquad(0<\varepsilon\le1,\ D\ge1).
\tag{13}
$$

这里只控制未归一的记录加权来源。不能除以任意小$p_c$后领取各条件支误差，也不把$\operatorname{AW}(rb)$当作两个非对易量子算符的普通乘积。硬分箱始终未被微分。

## 5. 同一已给资源内的有限预算

复用1251的原符号范围，区别它与1249的Hessian常数$M$：

$$
B_{\rm plat}=324\frac{81!}{(40!)^2}2^{40},
\quad L_c=56B_{\rm plat},
\quad A_s=2^{98}(1+L_c)^2.
$$

三元分箱结果数至多

$$
D=(2A_s10^{10}+1)^3<2^{4096}.
\tag{14}
$$

这是明确有限上界；允许包含空箱。代码从阶乘精确重建，核$B_{\rm plat},L_c$与旧账相同，并核$D$实际位数，而不是断言结果数很小。两个数学标签寄存器可各用4096比特的维数上界，不据此推断物理实现成本。

沿原固定资源：

$$
M=2^{3059},\quad E_H=3093,\quad E_0=14008,
\quad E_F=95928+100M,
$$

$$
N=1+\max\left(
\left\lceil\frac{E_F+42}{3}\right\rceil,
\left\lceil\frac{E_H+E_F+42}{6}\right\rceil
\right),\qquad \varepsilon\le2^{-4N}.
$$

由(9)、(13)，

$$
q_L<2^{E_0+4098-N},\qquad q_{L,r}<2^{E_0+4099-N}.
\tag{15}
$$

1266原两程序在指定节点的测前来源误差指数为

$$
e_s=E_0+E_H+15M+22-2N.
$$

精确整数核验给$e_s<-82$以及式(15)两个指数均$<-82$。于是对这两份原准备的同一实际节点态，

$$
|\langle b\rangle_{\rm post}-b_{\rm cl}|
\le2^{e_s}+q_L<2^{-81}<2^{-80}.
\tag{16}
$$

AW正坏域、原有界能源加权尾及1251的$a_c$均使用各自已认证菜单。原硬旗标门未保证的无界共轭矩不会被式(16)修复。一般$\kappa$族使用式(3)，不扩大原测前量词。概率效果精确保持，所以旧统计失败率和新任务系统偏差合同沿原范围保留；本稿没有认证新增仪器的实际系统偏差。

## 6. 明确分开的测量—重备对照

另一份具有同样效果的仪器是

$$
\mathcal I_c(\rho)=\int_{\Omega_c}
\langle z|\rho|z\rangle |z\rangle\langle z|\,d\mu_\varepsilon(z).
\tag{17}
$$

正迹类Bochner积分存在且迹积分不超过$\operatorname{Tr}\rho$；线性扩展为正常CP映射。其参考输出为
$\int_{\Omega_c}|z\rangle\langle z|\otimes\langle z|\rho_{SR}|z\rangle_S\,d\mu_\varepsilon(z)$。
这是公开的重备操作，不是Lüders后态。

记$\mathcal B=\sum_c\mathcal I_c$。原热核身份给

$$
\mathcal B^\dagger(\operatorname{Op}_\varepsilon b)
=\operatorname{Op}_\varepsilon(e^{\varepsilon\Delta/2}b),
\qquad
\mathcal B^\dagger(\operatorname{AW}_\varepsilon b)
=\operatorname{Op}_\varepsilon(e^{3\varepsilon\Delta/4}b).
\tag{18}
$$

两种非选择来源扰动均$\le\varepsilon F_2(b)/2$。记录加权Weyl来源恰为$\operatorname{AW}(r e^{\varepsilon\Delta/4}b)$，相对$\operatorname{AW}(rb)$的误差$\le\varepsilon F_2/4$；AW来源则$\le\varepsilon F_2/2$。这些更小来源界不能替换主仪器的式(8)。

对七模相干输入，Gaussian积分给

$$
\langle z_0|\mathcal B(|z_0\rangle\langle z_0|)|z_0\rangle
=\int|\langle z_0|z\rangle|^4d\mu_\varepsilon(z)
=2^{-7}=1/128.
\tag{19}
$$

故**本测量—重备通道**的半迹距至少$127/128$，不是主Lüders通道的结论。它显示小来源变化与小全态变化完全不同。

其一个等距实现使用连续环境$\mathcal E=L^2(\mathbb R^{14},d\mu_\varepsilon)$：

$$
(V_{\rm MP}\psi)(z)
=|z\rangle_S|c(z)\rangle_C|c(z)\rangle_{C'}\langle z|\psi\rangle .
\tag{20}
$$

积分范数恰为$\|\psi\|$。这种特定实现的连续环境不应成为本轮必须补齐的门槛；主仪器式(2)已经给有限环境替代。两者都未指定记录环境的物理Hamiltonian、来源、准备成本及CC实现，且不能把这些新增来源自动塞进1266原$J+V$抵销。

## 7. 去重、复算与实际范围

- 1249已有相干期望与AW热平滑，用于传播；1251只给效果与概率。
- 1266的$cq$后态属于B旗标，不能替统计后态。
- 历史检索命中1073有限投影测量—重备比较及1131、1219的相干/热指针内容；对象不是本七对有限分箱与原来源连接。一般Lüders或相干测量工具不是新物理定律。

[instrument_check.py](instrument_check.py)登记一个作者科学组，做三类互补核验：

1. 有限Parseval分析框架中直接核效果、平方根、完整等距、交换子和记录加权分解；该有限矩阵诊断只检验一般代数，不冒充原14维连续模型的数值传播。
2. 48点Gauss–Hermite逐坐标核Fourier字符的14维热因子及式(19)；Fraction交错级数包络核字符缺陷和明确有限Fourier菜单的热误差。样本不代替式(4)—(13)的解析全输入定理。
3. 从冻结结果重建$A_s,D,M,E_H,E_0,E_F,N$，以整数核同一资源预算，不计算无法实际展开的$2^N$。

首次 --write 排他创建[instrument_results.json](instrument_results.json)；默认只读重算匹配。绑定本稿、代码、准入和冻结来源；不绑定并行编写的几何文件。未做新ODE、全量子场模拟或实际设备测试。

**净增：原有限统计效果、具名有限环境后态、原有界来源和未归一记录加权来源组成同一数学终末操作。** 剩余是把这份具名原语纳入实际有限世界管时的完整参与者来源、边界和误差登记；可公开选择有限有效输入，不要求无限制造链。该登记不得包含待证明的几何预测作为定义。本轮不证明全量子GR、经验引力或未知态无扰，不新增认知公理，也不修改1266。
