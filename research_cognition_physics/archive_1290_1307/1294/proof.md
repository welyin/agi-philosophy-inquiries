# 1294：有限规范环的局域记录与共同热参考／来源

2026-10-18。按[正式准入](../_shared/round1294_admission.md)固定全部参数；这是新有限图有效候选中的一项共同过程，不是既有公理的普遍反模型。新增公理为零，正式状态以独立审阅及主线验收为准。

## 1. 旧工具、输入与本轮增量

[617 §3—5](../../archive_585_628/research_note_617.md)已给局域带链路CAR跳跃及共同Gauss表示运输；[1021 §3—6](../../archive_1009_1043/research_note_1021.md)已区分完整过程扇区和同一Hamiltonian的热参考；[557](../../archive_554_584/research_note_557.md)、[559](../../archive_554_584/research_note_559.md)要求指针、相互作用及读出注能一起保留。这些成熟工具直接继承。本轮新增的是一份显式局域过程：非收缩环资料通过带荷运输进入当地中性记录，并由同一个完整Hamiltonian计算热参考及声明参数来源；不直接指定Wilson环效果。

比较 $N=2,3$ 的四站有向环 $j=0,1,2,3$（下标模4）。每条链用电基 $|e_j\rangle$，$e_j\in\mathbb Z_N$，且 $X_j|e_j\rangle=|e_j+1\rangle$。环上有一粒正荷CAR物质 $c$；0站有一粒固定占据、无跃迁的负荷物质 $d$。二者质量系数分别 $m_c=m_d=2$。固定一正一负荷是准备／部门输入；负荷无跃迁是有效参数，不是假设无限质量或丢掉它的能源。

2站中性二级记录器的初态为 $|0\rangle_p$，数目 $n_p=|1\rangle\langle1|$。取

$$
\begin{aligned}
K&=\sum_{j=0}^3(c_{j+1}^\dagger X_j^{-1}c_j+\mathrm{h.c.}),\\
H(\lambda)&=m_cN_c+m_dn_d-gK+a_pn_p+r n_{c2}X_p,\\
(m_c,m_d,g,a_p,r)&=(2,2,1,1,1),\qquad
 t_*=2^{-9},\quad\beta_*=2^{-10}.
\tag{1}
\end{aligned}
$$

$X_p$ 是记录器Pauli矩阵，不是链路 $X_j$。所有作用时间无关；没有需要开关而未计账的外部脉冲。裸电能系数在本模型明确取零，但完整 $H$ 非零。本文不借纯Dijkgraaf–Witten的零Hamiltonian声称本模型零源，也不借有隙格点模型的额外操作权限。实际时间归一、固定图和各参数为输入；未构造相对论时钟或几何度规。“非收缩”只指所选四站图的环，未证明连续三维空间的拓扑或四维平坦TQFT资格；实体回路若可缩，其磁通／缺陷及来源的映射另需证明。

在一粒 $c$ 部门，$K$ 的范数不超过2，$n_{c2}X_p$ 范数为1，因此

$$
H\ge I,\qquad \|H-4I\|\le4,\qquad H\le8I.
\tag{2}
$$

全部矩阵有限，正规态、能源及任意固定阶矩均有限。质量常数只在传播计算中提出公共相位，始终保留在物理能源及参数来源中。

## 2. 精确Gauss空间与共同未知输入

Gauss条件为

$$
e_j-e_{j-1}=n_{cj}-\delta_{j0}\pmod N.
\tag{3}
$$

正向跳跃使电场 $e_j$ 降1，恰与两端物质电荷变化配平。式(1)全部项和当地记录读口都保持Gauss。固定 $c$ 在 $x$ 时，全部解恰为

$$
 e_j(x,\ell)=\begin{cases}
 \ell-1,&0\le j<x,\\
 \ell,&x\le j\le3,
 \end{cases}\qquad\ell\in\mathbb Z_N.
\tag{4}
$$

记相应正常物理基为 $|x,\ell\rangle$，连同记录器构成**整个所声明固定电荷部门**，维数 $8N$。它不是任意其它电荷／边界部门的分类。原CAR在单粒子部门的跃迁就是位置矩阵；补偿者仍在全Fock描述中，所有使用的操作均为偶算符。

令 $\phi_q=2\pi q/N$，

$$
 |x,q\rangle=N^{-1/2}\sum_{\ell=0}^{N-1}e^{i\phi_q\ell}|x,\ell\rangle.
\tag{5}
$$

前三条边不变 $\ell$，末边使 $\ell\mapsto\ell-1$，故精确分块为

$$
\begin{aligned}
K_q&=\sum_{j=0}^2(|j+1\rangle\langle j|+\mathrm{h.c.})
  +e^{i\phi_q}|0\rangle\langle3|+e^{-i\phi_q}|3\rangle\langle0|,\\
H_q&=4I-K_q\otimes I+I\otimes n_p+|2\rangle\langle2|\otimes X_p.
\tag{6}
\end{aligned}
$$

分块用于证明，不是许可调用一个全局Wilson测量。共同二级逻辑输入 $Q$ 在两群均编码为

$$
J_N|b\rangle_Q=|x=0,q=b\rangle|0\rangle_p,\qquad b=0,1.
\tag{7}
$$

$J_N^\dagger J_N=I_2$，且 $J_N^\dagger HJ_N=4I_2$、$J_N^\dagger H^2J_N=18I_2$。所以任意未知 $\rho_{QR}$ 的共同准备有相同有限初始能源4及二阶矩18；不以两个任意不同材料准备代替共同输入。环绕量子态本身是登记的准备资源，未声称由局域工厂免费生成。初始当地记录器与 $QR$ 独立。

实际等待 $e^{-itH}J_N\otimes I_R$ 是全等距。若末端调用当地记录器的理想 $Z_p$ 仪器，分支为

$$
 V_{N,b}= (I\otimes|b\rangle\langle b|_p)e^{-itH}J_N,
 \qquad \sum_b V_{N,b}^\dagger V_{N,b}=I_Q.
\tag{8}
$$

该式保所有结果及相应物理后态，并可张量任意 $I_R$。非选择 $R$ 边缘不变，但条件关联可变；记录经典化不保证原未知输入可逆恢复。不同群并没有被宣称拥有同一通道：本轮恰检验它们对同一编码任务的不同响应。

## 3. 全阶有限时间记录证书

置 $\widetilde H_q=H_q-4I$、$P_1=I\otimes|1\rangle\langle1|_p$、$v=|0,0\rangle$。由原局域跃迁直接相乘，

$$
 P_1\widetilde H_q^jv=0\quad(j=0,1,2),\qquad
 P_1\widetilde H_q^3v=(1+e^{-i\phi_q})|2,1\rangle.
\tag{9}
$$

两跳到达站2后才有一次记录翻转，因此不是在初片直接读取全局标签。Taylor尾由式(2)控制：

$$
 \left\|P_1e^{-it\widetilde H_q}v-
 \frac{(-it)^3}{6}(1+e^{-i\phi_q})|2,1\rangle\right\|
 \le \frac{(4t)^4}{24}e^{4t}.
\tag{10}
$$

在 $t=t_*$，用 $e^{1/128}<128/127$ 得右侧严格小于 $t_*^3/40$。因此，记 $p_{N,q}$ 为记录1概率，

$$
\begin{aligned}
p_{N,0}&>\frac{1369}{14400}t_*^6,\\
\frac{289}{14400}t_*^6<p_{3,1}=p_{3,2}
 &<\frac{529}{14400}t_*^6,\\
p_{3,0}-p_{3,1}&>\frac7{120}t_*^6.
\tag{11}
\end{aligned}
$$

其中 $p_{3,1}=p_{3,2}$ 由带补偿相位的环反射推出：令 $P_{13}$ 交换位置1、3，$D_\phi=\operatorname{diag}(1,e^{i\phi},e^{i\phi},e^{i\phi})$，则 $R_\phi=D_\phi P_{13}$ 满足 $R_\phi K_\phi R_\phi^\dagger=K_{-\phi}$，并保持初位置及站2数目；不能泛用复共轭代替真实时间演化。$q=0$ 两群的整个块及初态相同。

对 $N=2,q=1$，$\phi=\pi$。环上二维子空间

$$
\operatorname{span}\{|0\rangle,(|1\rangle-|3\rangle)/\sqrt2\}
\tag{12}
$$

在 $K_\pi$ 下不变，站2数目恒零；记录初态0也不受记录能隙作用。因此 $p_{2,1}(t)=0$ 对所有 $t$ 精确成立。共同逻辑输入 $|1\rangle_Q$ 的跨群差满足

$$
 p_{3,1}(t_*)-p_{2,1}(t_*)>
 \frac{289}{14400}2^{-54}>\frac{2^{-54}}{50}>0.
\tag{13}
$$

式(10)的尾除以 $t^3$ 后随 $t>0$ 增加，故同一界对 $0<t\le t_*$ 全成立。在固定窗口 $[1/1024,1/512]$ 内，跨群差统一严格大于 $2^{-60}/50$；这只是原证书的推论，没有优化参数或增加测试任务。

末端二结果记录分布的总变差就是该概率差。任意扰动若已有本任务输入／仪器误差合同，各分支误差从式(13)相减即可；这里没有声称实际仪器达到这一极小误差。若只改变生成元且 $\|\delta H\|\le\epsilon_H$，有限等待通道对任意伙伴的半diamond距离不超过 $t\epsilon_H$，所以两候选的概率差损失不超过两份此预算之和。这个条件界不引入新的电项或设备。

## 4. 同一全Hamiltonian的参考不是均匀拓扑参考

在完整 $8N$ 维固定电荷部门定义

$$
 \rho_{\beta,N}=\frac{e^{-\beta H_N}}{Z_N},\quad
 Z_q=\operatorname{Tr}_{8}e^{-\beta H_q},\quad
 \bar Z_q=e^{4\beta}Z_q,
 \qquad \Pr_\beta(q)=\frac{\bar Z_q}{\sum_{j=0}^{N-1}\bar Z_j}.
\tag{14}
$$

这是所选完整能源的忠实有限Gibbs比较参考；不是自动热化，也没有要求工作初态等于参考。参考包含带荷粒子、补偿者、记录器和持续相互作用，不能把它换为单独环部门的 $I/N$ 再沿用原账。

对式(6)，长度小于4的闭路径不绕环。第四阶只有四条顺向／逆向环路及两个记录状态贡献环相位，故

$$
\begin{aligned}
\operatorname{Tr}\widetilde H_0^j-
 \operatorname{Tr}\widetilde H_q^j&=0\quad(0\le j\le3),\\
\operatorname{Tr}\widetilde H_0^4-
 \operatorname{Tr}\widetilde H_q^4&=16(1-\cos\phi_q),\\
\bar Z_0-\bar Z_q&=\frac23(1-\cos\phi_q)\beta^4+R_5,\\
|R_5|&\le\frac{16(4\beta)^5}{120}e^{4\beta}.
\tag{15}
\end{aligned}
$$

在 $\beta_*$，$e^{1/256}<2$ 使 $|R_5|<(4/15)\beta_*^4$。$N=3,q=1,2$ 的领先系数为1，$N=2,q=1$ 为 $4/3$，所以两者均有

$$
\bar Z_0-\bar Z_1>\frac{11}{15}\beta_*^4,\qquad
\Pr_{\beta_*}(0)-\Pr_{\beta_*}(1)>
 \frac{11}{240N}\beta_*^4>0.
\tag{16}
$$

后一界使用 $\bar Z_q\le8e^{4\beta_*}<16$。实际未提出质量常数的 $Z_0-Z_1$ 还乘 $e^{-4\beta_*}$，不可漏写。可读取不普遍强迫配分不同，例如等谱扇区仍可有不同多时记录；式(16)是本模型另外核得的结论，绝非以“可读”定义参考改变。

## 5. 完整参数来源、能源守恒与末读账

本轮声明的来源菜单是式(1)同一开参数族的所有五个参数导数：

$$
 S_{m_c}=N_c=I,\quad S_{m_d}=n_d=I,\quad
 S_g=-K,\quad S_{a_p}=n_p,\quad S_r=n_{c2}X_p.
\tag{17}
$$

这些来源在两群均由原物理矩阵定义；不能只保存物质或裸环能源而删 $S_r$。若明确给定参数依赖于有限背景标量 $\gamma$，其响应严格为 $\partial_\gamma H=\sum_\lambda(\partial_\gamma\lambda)S_\lambda$；未给的时空几何耦合、Hilbert应力及其独立变分不由式(17)自动生成。

对任何初态及未知伙伴，完整等待总能守恒。将 $H$ 分成静止、跃迁、记录、相互作用四项，每项真实变化由 $i[H,H_j]$ 给出，四项求和为零；相互作用能源一般不能从记录能源中删去。参考满足

$$
 \partial_\lambda\log Z_N=-\beta\operatorname{Tr}(\rho_{\beta,N}S_\lambda).
\tag{18}
$$

使用同一 $H$ 和同一态求平均，不把关联来源拆成边缘均值的乘积。式(18)是成熟有限热恒等式，本轮没有据此唯一选定全部耦合或群。

对式(8)的声明理想末读，非选择去相干 $\mathcal D_p$ 保所有对角记录项而删除 $n_{c2}X_p$，于是

$$
 \mathcal D_p^*(H)-H=-r n_{c2}X_p,\qquad
 \Delta E_{\rm read}=-r\langle n_{c2}X_p\rangle,
 \quad |\Delta E_{\rm read}|\le |r|\langle n_{c2}\rangle\le |r|.
\tag{19}
$$

该能源变化可以有符号，必须由尚未纳入本 $H$ 的末读装置之已声明工作账承担；不等于最低装置功或实际制造成本。条件后态和记录均保留，不能以非选择均值相等冒称各支无反作用。所选当地指针及其持续作用已经在完整 $H$ 内，准备、计时和末端理想仪器仍为公开有效权限输入。本轮只签这一有限记录形成／声明末读的共同参数源合同，不签外读设备全部微观控制和存储，也不把其未制造提升为原目标必须另建的通用工厂。

## 6. 删除、复算及选择范围

- 删除 $r$ 或 $g$ 后，初始记录0不可能获得本任务的环信息；保全部耦合时式(13)恢复。这是同原生成元的删除检验。
- 保持原链路Hilbert空间、Gauss条件和共同准备，只把末边跃迁项置零，则有效运输图为树，所有保留跃迁边相位可由位置相位变换去除，原固定位置输入及当地记录概率不辨环标签。原环绕资料仍在完整空间中；这里删除的是本菜单访问该资料的闭路运输，未删除共同码。
- 仅测一个直接声明的Wilson效果不替代式(9)—(13)。仅指出热权不均匀也不足以构成本轮净增。
- 固定电荷和初态、群及环、零裸电能、参数和理想读出权限仍是输入。允许任意新电项时式(14)的扇区分解可能失效；需另用全Hamiltonian证明，不能继续引用本式。

[check.py](check.py)只作一个作者核验组：从完整电基枚举Gauss核并构造物理等距，核全部物理空间／Fourier块／来源和未知伙伴分支；用独立的Laurent闭路径算术核三阶记录项及第四阶热项；用有理尾给严格界，小矩阵传播仅为诊断。默认或 `--check` 只读比对 [results.json](results.json)，`--write` 仅首次排他创建。

本轮完成有限图候选中局域运输、当地物理记录与共同热参考／来源的一项绑定。它既不把 $\mathbb Z_2$ 或 $\mathbb Z_3$ 唯一选成自然界规范群，也不证明完整四维拓扑场论、全部认知公理或HG热—几何共同实现。成功后停止扩建本环器件，返回主线结构判断。
