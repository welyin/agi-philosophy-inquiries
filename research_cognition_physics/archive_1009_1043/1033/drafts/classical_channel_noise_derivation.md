# 1033解析稿：经典测量反馈、作用强度和有限噪声门槛

2026-10-08。本文证明一个明确的有限维有效类别内的充要定理。认知原则不预先指定此类别；引力应用另需实际对象映射。结果与成熟经典通道噪声研究相接，不作为新物理定律发表。

## 1. 合同与量词

两方各有一个二态自由度，令 $A=Z\otimes I$、$B=I\otimes Z$，在声明的相互作用表象采用时间齐次生成元

$$
\mathcal L\rho=-ig[AB,\rho]+\gamma_A\mathcal D[A]\rho+
\gamma_B\mathcal D[B]\rho,
\qquad \mathcal D[C]\rho=C\rho C^\dagger-\tfrac12\{C^\dagger C,\rho\}.
\tag{1}
$$

$g\in\mathbb R$、$\gamma_A,\gamma_B\ge0$，均用逆时间单位；物理Hamiltonian为 $\hbar g AB$。此时 $\mathcal D[A]\rho=-[A,[A,\rho]]/2$。采用固定A/B划分、同一准备无关的线性过程、Markov近似和纯局部Z退相干。一般相关噪声、记忆、耗散、时间延迟及未知额外相互作用不在此合同。

经典实现指局部量子仪器、可共享经典随机性及经典通信，不允许预存跨方纠缠、量子传输或隐藏的共同量子门；任意精度极限属于 $\overline{\mathrm{LOCC}}$。此闭包不是声称精确过程有有限步骤或一个统一有限资源上界。

**定理。** 以下三项等价：

1. $e^{t\mathcal L}$ 对每个 $t\ge0$ 保持所有A/B可分输入可分；
2. $e^{t\mathcal L}\in\overline{\mathrm{LOCC}}$ 对每个 $t\ge0$；
3. $\gamma_A\gamma_B\ge g^2$。

“每个t”是这个有效半群的数学分类，不能换成仅在一个末时刻观察不到纠缠；后者明显更弱。例如零噪声、$gT=\pi/2$ 的演化是 $-iZ\otimes Z$，末态操作局部，但中间演化会产生纠缠。物理排除则用下文的**具体有限时间和非零误差裕量**，不要求现实微观连续。

## 2. 必要性：零本征空间中的一阶障碍

取 $\rho_0=|++\rangle\langle++|$，对B作部分转置 $\Gamma$。$\rho_0^\Gamma=\rho_0$ 的零空间包含 $u=|-+\rangle$、$v=|+-\rangle$。直接算式(1)在该二维空间中的导数压缩：

$$
\begin{pmatrix}\langle u|\\\langle v|\end{pmatrix}
(\mathcal L\rho_0)^\Gamma
\begin{pmatrix}|u\rangle&|v\rangle\end{pmatrix}
=\begin{pmatrix}\gamma_A&-ig\\ig&\gamma_B\end{pmatrix}=:M.
\tag{2}
$$

对角非负。若 $\det M=\gamma_A\gamma_B-g^2<0$，取其负本征向量w，则

$$
\langle w|(e^{t\mathcal L}\rho_0)^\Gamma|w\rangle=t\lambda_-+O(t^2)<0
$$

在足够小的正t成立。可分态的部分转置必为正，因此此演化会产生纠缠。这里仅需可分态的PPT必要性，不依赖PPT充分分类。LOCC及其闭包保持可分性，故得到1与2各自蕴涵3。

## 3. 充分性：实际局部Kraus构造

若g=0，式(1)是两份局部退相干，结论直接成立。以下g非零，条件3迫使 $\gamma_A>0$。先引入两个局部测量强度 $\kappa_A,\kappa_B>0$，一向反馈采用

$$
c_A=\sqrt{\kappa_A}A,\quad F_B={g\over2\sqrt{\kappa_A}}B,
$$

反向交换A、B和对应κ。单向的测量反馈生成元为

$$
\mathcal L_{c,F}\rho=-i[\tfrac12(c^\dagger F+Fc),\rho]+\mathcal D[c-iF]\rho.\tag{3}
$$

不只引用式(3)：对一步δ，A处执行两结果测量

$$
M_s=\sqrt{\frac{I+2s\sqrt{\kappa_A\delta}A}{2}},\quad s=\pm1,
$$

把s作为经典消息发至B，B执行 $U_s=\exp[-isg\sqrt\delta B/(2\sqrt{\kappa_A})]$。当 $4\kappa_A\delta\le1$ 时，$\sum_s M_s^\dagger M_s=I$，所有Kraus算子 $U_sM_s$ 均为两方局部算子的乘积。展开并对s求和得 $\Phi_{A\to B}(\delta)=I+\delta\mathcal L_{c_A,F_B}+O(\delta^2)$。平方根在所用小δ邻域解析；半整数奇项由正负s相消。

同样反向执行，两份Hamiltonian各贡献 $gAB/2$；两份耗散中的反对称交叉项相消。合计

$$
\gamma_A=\kappa_A+{g^2\over4\kappa_B},\qquad
\gamma_B=\kappa_B+{g^2\over4\kappa_A}.
\tag{4}
$$

令 $\kappa_A=\gamma_A/2$、$\kappa_B=g^2/(2\gamma_A)$，得到边界率 $(\gamma_A,g^2/\gamma_A)$。再加B处非负的局部退相干 $\gamma_B-g^2/\gamma_A$ 即精确得到式(1)。有限维乘积公式保证N步、$\delta=t/N$ 时这些LOCC过程收敛至 $e^{t\mathcal L}$。故3蕴涵2，也蕴涵1。

这个构造没有把测量、控制器、随机性、时序和经典消息变成整体内部自然产生的免费资源。实际自治实现及其预算仍未认证。构造的价值是证明噪声门槛可达到，而不是要求认知整体就是一个外置反馈装置。

## 4. 完整有限时间解与可见度

计算基的标记为 $(a,b)\in\{\pm1\}^2$，矩阵元满足

$$
\rho_{ab,a'b'}(T)=\rho_{ab,a'b'}(0)
\exp\left[-igT(ab-a'b')-\frac{\gamma_AT}{2}(a-a')^2-
\frac{\gamma_BT}{2}(b-b')^2\right].
\tag{5}
$$

在B的Z本征态上准备A的等权相干态，A的归一化可见度为 $V_A=e^{-2\gamma_AT}$；反向准备给 $V_B=e^{-2\gamma_BT}$。两实验使用相同生成元与时长，局部相位按同一校准字典处理。四路径交互相位的连续标定值为 $\phi=4gT$。由算术—几何均值不等式，经典合同必有

$$
V_AV_B=e^{-2T(\gamma_A+\gamma_B)}\le e^{-4|g|T}=e^{-|\phi|}.\tag{6}
$$

等号要求 $\gamma_A=\gamma_B=|g|$；本轮不增“宇宙选择最小噪声”原则。式(6)是必要条件，单独通过它不能保证乘积门槛或全部合同成立。

φ是依靠同一动力学窗口校准的相位展开值；若只有模 $2\pi$ 的相位，必须对允许的展开分支分别核算，不能任意指定累计相位。实际相干资源和读取权限属于仪器输入，不因形式上可测就默认实现。

这些两次探针统计**不是对所有经典通道的设备无关证伪**。排除仍依赖同一有效生成元对允许未知输入有效、拟合误差已控制；若允许随准备改变机制，或只要求复现这两个已知输入而不承担该生成元，不能套用定理。

## 5. 一个有限时间、有限误差的解析证书

取 $g=1$、$\gamma_A=\gamma_B=1/4$、$T=1/100$。式(2)的负本征值为 $-3/4$。对应的可分态见证是

$$
W=(|w\rangle\langle w|)^\Gamma={I\otimes I-X\otimes X-Y\otimes Z-Z\otimes Y\over4},\quad \|W\|_\infty=1/2.\tag{7}
$$

它只需三组可信局部Pauli相关统计；局部读数记录后再作经典比较，不要求实施非局部投影。可分态满足 $\operatorname{Tr}W\sigma\ge0$。

有 $\|\mathcal L\|_{1\to1}\le2|g|+2\gamma_A+2\gamma_B=3$，且 $\|\rho_0\|_1=1$。用保守的 $\|W\|\le1$ 估余项，$e^x\le1/(1-x)$（$0\le x<1$）给

$$
\operatorname{Tr}W e^{T\mathcal L}\rho_0
\le-{3\over4}{1\over100}+{(3/100)^2\over2(1-3/100)}
=-{273\over38800}<-0.007.
\tag{8}
$$

直接有限矩阵值约 $-0.00746191990$，仅作校准。若真实过程、准备及对齐的共同误差以半迹距离 $\epsilon\le1/1000$ 控制，则式(7)使见证偏差至多ε，仍有

$$
\operatorname{Tr}W\rho_{\rm actual}\le-{1171\over194000}<0.
$$

读数本身另有误差时继续扣除；本文没有认证任何实验的ε。这个证书直接发生在给定有限时间内，不是只有无限精细外推才可见的矛盾。

式(6)也可直接容纳数据区间：令可见度下界为 $v^-_A,v^-_B\ge0$、相位绝对值下界为 $\phi^-\ge0$。若 $v^-_Av^-_B>e^{-\phi^-}$，合同不相容。例如纯示例数据两可见度各 $.99\pm.001$、相位 $.04\pm.001$，则

$$
(.989)^2-\frac1{1+.039}={16267719\over1039000000}>0,
$$

而 $e^{-.039}\le1/(1+.039)$。这些是假设性预算数字，不是已取得的引力实验数据。

## 6. 引力应用是待接桥梁

若两质量均有两条受控路径，并且势在该四维码上可用对角 $V_{ab}$ 表示，去掉各自局部项后，交互率才是

$$
g={V_{00}-V_{01}-V_{10}+V_{11}\over4\hbar}.\tag{9}
$$

采用Newton或945的软核势时，$V_{ab}=-GM_AM_Bu(r_{ab})$；G、核、路径准备、支撑、读取及误差仍是物理输入。[945](../../../archive_935_955/research_note_945.md)只有一端路径相干任务，**没有**给出式(9)需要的双路径交互与新噪声合同。不能把945的约.05相位直接冒作本轮φ。波包移动、泄漏、环境噪声及四路径投影必须另有控制，才可把现实数据接来。

式(1)中的Z噪声与AB对易，故这份二态有效H的平均能量守恒。但它没有动量、支撑与场能账，更没有认证反馈资源、相对论传播、完整应力及几何反作用；不可用这个小模型签收全部957或P981。

## 7. 旧结果、文献和本轮增量

[989](../../../archive_956_989/research_note_989.md)已排除既定电容任务中的LOCC替代，本轮不重称该事实为新定理。新增的是给候选经典引力通道类别一个可达到的精确噪声门槛，并把它传至同一有限相位、两端可见度与误差预算，约束I06与I14的联合采用。

[Kafri–Taylor–Milburn §III](https://arxiv.org/html/1401.0946)以测量反馈研究经典引力通道；[Kafri–Taylor](https://arxiv.org/html/1311.4558)讨论力与噪声的量子关联判据。本文接入其成熟机制，二态必要性及充分性的具体因子、有限见证均由上文独立给出。没有将高斯定理直接套到所有非高斯态，也没有无声采用最小噪声、连续微观时间或任意精确控制公理。

该结果不是所有经典几何描述的排除，不生成G、Einstein作用、空间或完整量子引力；它允许保留经典关系标签和量子实际作用的989分工。跨越该类别需要新的证明，不能仅把变量改名后延伸量词。
