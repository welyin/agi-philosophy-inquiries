# 1268读取分项：同一真实传播中的有限Ramsey尺度记录

2026-10-16。作者分项1组，按[准入](admission.md)；新增采用公理0，CT仍待验证。本文证明一项明确有限块量子读取及共同误差，不把它与[材料分项](material_proof.md)尚待核验的空间来源／联合因果资格拼成已实现的完整装置。

## 1. 同一准备、参数族与有限读口

复用[1267几何证明](../1267/geometry_proof.md)的原七对变量、真实生成元
$\widehat H_\kappa(t)=\operatorname{Op}_\varepsilon(\chi h_\kappa(t))$、已知标准Gaussian初态及完整来源字典。参数族是$B=2,\ \kappa\in[3,3.001]$；另取$\kappa=3,B=1,2$原程序对照。合法初$\pi_Q=-\sqrt{12\rho_0/\kappa}$沿旧规则。单位、杆$k^2=1/10$、$\rho_0=111/100$及物理$\hbar>0$均固定。

采用原有限整体块$V_0=10^{-6}$及已选定物理尺度

$$
s=1000m2^{2N},\quad m=\lceil\max(1,\hbar)\rceil,\qquad
\varepsilon=\frac{\hbar}{V_0s^2}\le2^{-4N}.
\tag{1}
$$

$N$仍是1249/1267保存整数，不另取趋零耦合同时免费保持相位。

新增二级指针$P$及有限末记录，准备$|+x\rangle$，门后读$\sigma_y$。简并基线可声明为原$0.71$材料的有限简并扩张；这是新增有效材料／准备输入，不是原定理已给该物件。数学分项只需指针自由生成元为标量；实际指针、记录和控制的独立来源仍需同对象认证，不能称零资源。

选前置窗：

$$
\tau=10^{-8},\ t_a=2-\tau,\ t_b=2,\ \alpha=1/16,\qquad
g(t)=\frac{140}{\tau}u^3(1-u)^3,\quad u=(t-t_a)/\tau .
\tag{2}
$$

窗外$g=0$，$\int g=1$，$\max g=35/(16\tau)=218750000$，是原$C^2$脉冲类型。整个全参数量子比较留在1267已认证$[0,2]$内，不延用未证明的族后期量子管。

令$\xi_Q$为1249原单$Q$平台，内核$[-1,2]$、支撑$[-2,3]$；新增读口是单坐标乘法算符

$$
F(Q)=\xi_Q(Q)e^Q,\qquad \widehat F=F(\widehat Q).
\tag{3}
$$

未使用全相空间$\chi e^Q$；后者会引入几何动量依赖。沿旧中心$F=a_{\rm cl}$。原平台常数$B_\beta=324[81!/(40!)^2]2^{40}$给
$\|F\|_\infty<27,\ \operatorname{Lip}F\le27(1+B_\beta)<2^{136}$。
左右过渡不重叠，故$\|\xi'_Q\|\le B_\beta$。此有界完成是公开输入，不是局部杆场已自动等同整体算符。

实际数学过程为

$$
i\varepsilon\partial_t\Psi_{\rm r}
=\bigl(\widehat H_\kappa(t)\otimes I+
\varepsilon\alpha g(t)\widehat F\otimes\sigma_z\bigr)\Psi_{\rm r}.
\tag{4}
$$

有界自伴生成元给真实酉传播。原写入、几何运动与全部七对自由度均保留。门强度$\varepsilon\alpha g$和相位$\alpha\int gF$来自同一作用比，不能分别调节。

## 2. 不冻结自由演化的全态界

原统一Gaussian范数引理直接承接1267、1266：

$$
M=2^{3059},\ E_H=3093,\ E_0=14008,\quad
S_*<2^{5M+4},\qquad
d_0:=\|\Psi_0(t)-\Psi_G(t)\|
\le\sqrt\varepsilon\,2^{E_H+15M+16}.
\tag{5}
$$

$\Psi_0(t)$是真实无新增门传播，$\Psi_G$只是二次比较。1267已逐项核整个参数族的旧平台、生成元宽界及来源常数，故(5)统一成立，不能由端点采样推论。

比较向量保留全部原自由演化：

$$
\Phi(t)=\Psi_0(t)\otimes
e^{-i\alpha\sigma_z\int_{t_a}^{t}g(v)a_{\rm cl}(v)\,dv}|+x\rangle .
\tag{6}
$$

积分在门前取0、门后取完整值。乘法Lipschitz界及Gaussian方差给
$\|(\widehat F-a_{\rm cl})\Psi_0\|\le54d_0+\operatorname{Lip}F\sqrt{\varepsilon/2}S_*$。
同一生成元的Duhamel因子$1/\varepsilon$与门项$\varepsilon$抵消，因此

$$
\|\Psi_{\rm r}(t)-\Phi(t)\|
\le\alpha\int|g|\|(\widehat F-a_{\rm cl})\Psi_0\|
<d_R:=\sqrt\varepsilon\,2^{E_H+15M+20}<2^{-100}.
\tag{7}
$$

固定资源整数余量验证最后一步。与Gaussian乘积比较另宽取
$d_{\rm all}=\sqrt\varepsilon\,2^{E_H+15M+21}$。这保留真实读取反作用，没有把小尾再除以$\varepsilon$。

对任意初始$PR$态、且已知Gaussian材料初态与$PR$独立，同一证明用$\|\sigma_z\|=1$成立，允许任意被动参考及其纯化。具体概率另要求指针为$|+x\rangle$；不声称任意未知几何态沿同一经典轨迹。完整末仪器
$\mathcal I_\pm(\rho)=\Pi_\pm^y\rho\Pi_\pm^y\otimes|\pm\rangle\langle\pm|$
保两支及相关性，迹距离收缩运输(7)。不后选、重备或免费复位；非选择读指针不改变原系统边缘，条件后态不冒称无扰。

## 3. 固定$a(2)$标定、有限分离与统计

旧端点给所有本族$a(2)\in(5/2,29/10)$，全前缀膨胀且$H<11/10$。用
$a(t)=a(2)e^{-\int_t^2H}$，短窗内$2<a<3,\ |\dot a|<4$。所以

$$
\bar a=\int ga_{\rm cl},\quad|\bar a-a(2)|<4\tau,\qquad
r(a)=\frac{1+\sin(a/8)}2,\quad p_{\rm cmp}=r(\bar a).
\tag{8}
$$

在$[2,3]$有$1/32<r'(a)=\cos(a/8)/16\le1/16$，由$\cos x\ge1-x^2/2$于$x\le3/8$即可。单任务时间平均概率误差$<\tau/4$，真实数学概率再偏移至多$d_R$。

实际运行另列新的匹配要求
$|p_{\rm physical}-p_{\rm r}|\le\beta_R=10^{-7}$。
它包括本指针准备、钟窗、末读和接口误差；不是1267旧POVM系统偏差的自动继承，也不由本证明置零。全部来源仍须与材料分项为同一过程。

每任务$R=5\cdot10^{14}$次独立重复，统计半宽$\delta_R=10^{-7}$。Hoeffding给失败率$\le2e^{-10}$；至多四个指定任务联合失败率$<8(3/8)^{10}<1/100$，不是不可数全参数族同时数据。成功事件内

$$
|\widehat p-r(a(2))|
<e_p:=2\cdot10^{-7}+\tau/4+2^{-100}.
\tag{9}
$$

旧两程序$\Delta a>683/10^7$，故

$$
\widehat p_{B=1}-\widehat p_{B=2}
>\frac{683}{32\cdot10^7}-\frac{\tau}{2}
 -4\cdot10^{-7}-2^{-99}
>\frac{17}{10^7}>0 .
\tag{10}
$$

旧$\kappa$端点$\Delta a>1/5000$给

$$
\widehat p_{3.001}-\widehat p_3
>\frac1{160000}-\frac{\tau}{2}
 -4\cdot10^{-7}-2^{-99}
>\frac{58}{10^7}>0 .
\tag{11}
$$

固定反演为
$\widehat a=r^{-1}(\operatorname{clip}_{[r(2),r(3)]}\widehat p)$，
$r^{-1}(p)=8\arcsin(2p-1)$只取该已知分支。投影不增误差，因此

$$
|\widehat a-a(2)|<32e_p<\eta_a:=1/150000 .
\tag{12}
$$

无需保留旧$a$读口的$10^{-8}$精度，二结果指针仍给非零几何分离。

## 4. 同态有界来源、正域及候选物理字典

对1249/1267已声明的有界来源、几何、正域尾菜单$b$，用同一Gaussian及$F_0(b)\le2^{E_0+1},F_1(b)\le2^{E_0}$，

$$
|\langle\widehat b\rangle_{\rm r}-b(z_{\rm cl})|
<\sqrt\varepsilon\,2^{E_0+E_H+15M+25}
+\varepsilon2^{E_0}/4<2^{-90}.
\tag{13}
$$

第二项保anti-Wick／Weyl差；正anti-Wick坏域及能源加权尾也受控。非选择指针末读保持这些系统边缘。这不控制未完成无界应力或连续监视路径。

候选直接杆调制是

$$
\Delta\mathcal L=-\sqrt{-g}\,
\frac{\varepsilon\alpha g(\theta)\sigma_z}{3k^2}\mathcal B_\phi,\qquad
\mathcal B_\phi=\sum_Ag^{\mu\nu}\partial_\mu\phi^A\partial_\nu\phi^A .
\tag{14}
$$

它不是带数密度$n$的质量耦合。齐次杆域$\mathcal B_\phi=3k^2a^{-2}$给坐标体积Hamiltonian$\varepsilon\alpha ga\sigma_z$，与(4)内核相同。独立变分给

$$
\Delta\rho=\varepsilon\alpha ga^{-2}\sigma_z,\qquad
\Delta P=-\frac13\varepsilon\alpha ga^{-2}\sigma_z .
\tag{15}
$$

关系钟重建$J_{\rm r}=-H_c-\varepsilon\alpha gF\sigma_z$。相同线性来源映射下，内核显式钟能源与(15)密度可补偿；压力不能删除，$g$端点归零也不保证钟余额恢复。几何动量新增项$-\varepsilon\alpha gF'\sigma_z$必须保留。

在声明的有界量子完成中，可将(15)的比较字典延拓为
$R_F=e^{-3Q}F,\ P_F=-e^{-3Q}F'/3$；
$\|R_F\|\le81,\ \|P_F\|\le27(1+B_\beta)$。
单坐标新增杆源是乘法算符，可用上述上确界。但原完成密度中钟补偿的过渡项是Weyl符号$\chi e^{-3Q}F=\chi\xi_Qe^{-2Q}$，不能以其点值上界冒称算符范数。它含两份原平台和一个复半径1内小于$3^6$的解析因子；1267 §4的Cauchy／Leibniz源类可直接复用；更显式地，$L_c=56B_\beta$时导数界$81r!(2L_c+2)^r$与旧Fourier前因子$C_F=2^{34}(5\cdot256^{13})15^{10}20!$给$F_0\le C_F81(2L_c+2)^{20}<2^{E_0}$，代码作精确整数检查。因此不在管外假设抵销，保守相加得到

$$
\varepsilon\frac{35}{256\tau}
\left[162+27(1+B_\beta)+2^{E_0}\right]
<\varepsilon2^{E_0+26}<2^{-100}.
\tag{16}
$$

这也覆盖钟峰值。外域来源选择是有效完成输入，没有自动成为(14)的完整空间协变实现；原$\chi$过渡代价明确计入，而非按乘法上界删除。

若只审齐次经典$\sigma_z=\pm1$分支，有界力与旧全局Lipschitz界$M$给
$\|\delta z\|\le\varepsilon\alpha\operatorname{Lip}F e^{2M}
<\varepsilon2^{3M+136}$。
合并旧来源导数宽界及门质量，正域裕量损失小于$2^{-100}$，中心不离旧内核。该比较不把量子叠加改成两份各自经典时空。

(13)—(16)只给同一约化候选的有界账本。物理$\theta/\sigma_z$承载、整体块到局部杆的映射、全部材料因果性及真实末读者来源仍由材料审查／具名合同负责。固定杆子块共锥不能代替全模型审查，不拼接不同字典。小归一强度不等于实际总成本为零，仍按(1)恢复体积、时间及能源。

## 5. 参数反演及场域停止线

沿1267全族$\partial_\kappa a(2)>1/5$，

$$
|\kappa-3|\le5(|\widehat a-\widehat a_3|+\eta_a+10^{-16}).
\tag{17}
$$

若接受$|\widehat a-\widehat a_3|\le\eta_a$，则
$|\kappa-3|\le1/15000+5\cdot10^{-16}$。这是条件预测，不是实测；其它物质参数、初态或$\Lambda$同时自由时不能沿用。

在真正未受新增门影响的共同经典前缀$0\le t\le t_a$，原独立来源$\|T\|_{\max}\le111/100$仍给所声明Einstein形状比较族

$$
\|G[g_\kappa]-3T_\kappa\|_{\max}
\le\frac{111}{100}
 \left(\frac1{15000}+5\cdot10^{-16}\right).
\tag{18}
$$

本轮不把(18)延到新读取后的$9/4$历史。门内／门后实际场余项须先将(14)完整控制与来源接入同一模型，再运输(13)、(16)；有界账本不是量子Einstein算符式。原未读取历史的9/4结果继续有效，但不是新历史。

数学指针在$t=2$以后解耦且自由简并，可以延后释放所存相位；不把未来$a$换名为$a(2)$。实际可访问事件、传播、释放控制及来源另按许可合同，不从极短$\tau$推出全块瞬时权限，也不要求完整器件制造链。

## 6. 核验及结论

[代码](readout_check.py)计一组精确整数／有理数核验，绑定旧全族证书，核平台范数、同一$N$余量、固定时标、统计、逆稳定及新增有界来源。旧ODE不重跑，未模拟全部量子波函数。默认只读比较[结果](readout_results.json)，首次--write排他创建。

净增是同一已知准备和真实生成元到具体有限指针、固定$a(2)$与有界来源的共同连接。实际匹配、来源及空间因果资格未由相位证明授予，不能单独宣布GR阶段完成。
