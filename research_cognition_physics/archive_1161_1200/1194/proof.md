# 1194证明：完整来源到有限端点概率的条件桥

2026-10-14。接续1193。本文在同一准备、同一完整演化和同一来源上，将质心界转为固定实际端点任务的概率界。关键控制量是声明输入码上的能源压缩，不要求物理Hamiltonian全谱狭窄，也不执行能量过滤。来源正性、实际CO读口及时间字典仍须明确匹配。

## 1. 净问题与直接复用

[1193](../research_note_1193.md)控制完整应力质心，但均值不等于端点分布。本轮补充

$$
\text{同一正来源的质心与绝对矩}
+\text{输入码能源控制}
+\text{原CO效应匹配}
\Longrightarrow
\text{固定端点任务的有限概率界}.
\tag{1}
$$

以下成果直接复用，不重新认领：

- [182](../../archive_171_189/research_note_182.md)区分线性未归一过滤、归一条件态与确定性过程；[1083](../../archive_1063_1085/research_note_1083.md)保留完整参与权重与归一代表坐标的区别。
- Terno的能源归一化位置POVM见[Moretti，§5.1式(35)、§5.2定理22](https://link.springer.com/article/10.1007/s11005-023-01689-5)。该文在平直时空的有质量实标量单粒子部门处理此构造；本轮不把其正性、因果局域化或实际仪器资格移植给任意量子应力。
- [1117](../../archive_1106_1119/1117/proof.md)区分原CO第一矩、载体质心、概率分布与测后态；[1112](../../archive_1106_1119/1112/proof.md)给原活钟的完整编码和能源保持逻辑过程。

净增是将这些工具接到1193的有限端点缺口，明确能源控制发生在哪个码、如何保参考，以及哪些实际匹配尚未取得。

## 2. 同一准备、完整过程和来源测度

### 2.1 输入码与完整来源

令 $Q$ 为事前声明的有限维未知输入码，$R$ 为任意被动参考。固定准备等距 $J$ 及同一完整演化 $U_t$，记

$$
W_t=U_tJ.
\tag{2}
$$

输出母空间保留全部材料、内部钟、程序、记录、控制、关系承担者和环境。已知混合资源可在准备中纯化，不控制或复制未知 $R$。无限维包络须有共同形式域和尾控制，不能以数值截断替代。

指定1193中的同一体 $B$、同一坐标图和同一来源分配，令 $\widehat{\mathcal E}_{B,t}(\Delta)$ 是区域 $\Delta$ 内的密度化能源来源。连续形式对应同一 $\widehat\tau^{00}_{B,t}$ 的积分。定义

$$
A_t(\Delta)=W_t^\dagger\widehat{\mathcal E}_{B,t}(\Delta)W_t,
\qquad E_t=A_t(\mathbb R^3).
\tag{3}
$$

$A_t$ 必须可数可加，并额外满足

$$
A_t(\Delta)\ge0,\qquad
0<aI_Q\le E_t\le bI_Q<\infty.
\tag{4}
$$

正总能源或母Hamiltonian正谱不能代替区域正性。任意量子场并不自动满足经典局部能量条件；某个声明准备码上的压缩正性可以作为较弱合同，但仍需独立认证。

**能源身份不能只靠同名。** 1193使用的是固定坐标密度矩
$E=\int\sqrt{|g|}T_B^{00}\,d^3x$。曲率或时变背景中，它不自动等于某个物理Hamiltonian、Killing能量或守恒Noether荷。只有同对象字典已经证明
$\widehat{\mathcal E}_{B,t}(\mathbb R^3)=\widehat H_{B,t}$，才可将式(3)写成 $W_t^\dagger\widehat H_{B,t}W_t$。第5节只检查原装置Hamiltonian的压缩，不自动完成这项来源身份。

以下固定时刻并省下标。合同须统一覆盖声明输入族；若仅对受限族成立，结论也限于该族。

### 2.2 数学POVM与原CO效应

由式(4)定义

$$
N(\Delta)=E^{-1/2}A(\Delta)E^{-1/2}.
\tag{5}
$$

它正、可数可加，且 $N(\mathbb R^3)=I$，故为码上POVM。这只是数学参照。能源归一化可能非局域，不能由式(5)领取CC局域装置或实际实施权限。

压缩与归一化一般不交换：

$$
E^{-1/2}W^\dagger\widehat{\mathcal E}(\Delta)WE^{-1/2}
\ne
W^\dagger\widehat H_B^{-1/2}
\widehat{\mathcal E}(\Delta)\widehat H_B^{-1/2}W.
\tag{6}
$$

同样，$\sqrt{W^\dagger\widehat H_BW}$ 一般不等于
$W^\dagger\sqrt{\widehat H_B}W$。后者若无定义，亦不能先行使用。本文不实际施加 $\sqrt{\widehat H_B}$。

对事前固定的二元端点任务，令 $F_{\rm CO}(O)$ 为原实际仪器报告“位于 $O$”的效应，拉回同一输入码。明示匹配合同

$$
\|F_{\rm CO}(O)-N(O)\|_{\rm op}\le\epsilon_{\rm CO}.
\tag{7}
$$

它需由同一实际端口、局部过程、时间、资源和读口校准承担，不由“位置”名称或平移协变获得。本文只使用一个固定二元事件；这不自动给所有细分位置的完整仪器界或系统测后态。

## 3. 对未知输入和参考统一的能源重权界

对任意 $\rho_{QR}$，定义数学比较态

$$
\rho_E=
\frac{(\sqrt E\otimes I_R)\rho_{QR}(\sqrt E\otimes I_R)}
{\operatorname{Tr}[(E\otimes I_R)\rho_{QR}]}.
\tag{8}
$$

分母至少为 $a>0$。该映射一般非仿射，不是确定性CPTP通道，也不作为另一条通道计算diamond距离。本轮不执行过滤，不修改实际未知输入，不删除任何失败分支。

### 3.1 迹距离引理

采用 $D(\rho,\sigma)=\frac12\|\rho-\sigma\|_1$。令

$$
s=\sqrt a,\qquad t=\sqrt b,\qquad
\delta=\frac{t-s}{t+s}.
\tag{9}
$$

则对所有 $\rho_{QR}$，

$$
D(\rho_E,\rho_{QR})\le\delta.
\tag{10}
$$

**证明。** 纯化 $\rho_{QR}$ 为 $|\Psi\rangle_{QRS}$，设
$m=\langle\sqrt E\otimes I\rangle_\Psi$、
$h=\langle E\otimes I\rangle_\Psi$。由谱界

$$
(\sqrt E-sI)(\sqrt E-tI)\le0
\quad\Longrightarrow\quad
h\le(s+t)m-st.
\tag{11}
$$

又

$$
(s+t)^2m^2-4st[(s+t)m-st]=[(s+t)m-2st]^2\ge0,
$$

所以

$$
\frac{m^2}{h}\ge\frac{4st}{(s+t)^2}.
\tag{12}
$$

归一重权后的纯化向量与原向量的重叠为 $m/\sqrt h$；纯态迹距离至多
$\sqrt{1-4st/(s+t)^2}=(t-s)/(t+s)$。部分迹收缩给式(10)，保留任意被动 $R$。证毕。

若 $E=e_*I$，可取紧谱界 $a=b=e_*$，从而 $\rho_E=\rho$、$\delta=0$。这只要求输入码上的能源压缩为标量，不要求母体码不变、定能或母Hamiltonian全谱平坦。

### 3.2 来源概率与POVM概率

正来源定义概率

$$
q_\rho(\Delta)=
\frac{\operatorname{Tr}[A(\Delta)\rho_Q]}
{\operatorname{Tr}[E\rho_Q]}
=\operatorname{Tr}[(N(\Delta)\otimes I_R)\rho_E].
\tag{13}
$$

它一般不等于 $\operatorname{Tr}[N(\Delta)\rho_Q]$，二者事件概率差由式(10)控制为 $\delta$。

例如 $E=\operatorname{diag}(1,4)$、$A(O)=\operatorname{diag}(0,4)$ 时，两本征输入的 $q(O)$ 是0与1，而等概率混合得到 $4/5$，不是 $1/2$。此例只说明182已知的归一化边界，不计新成果。

### 3.3 二元记录与参考的联合统计

令 $\mathcal M_N$ 只输出 $O$／补集的经典记录 $C$；测量系统的后态不在本界的输出中。记 $A_1=A(O)$、$A_0=E-A(O)$。数学来源重权cq态为

$$
\begin{aligned}
\Omega^A_{CR}(\rho)
&=\frac{\sum_{c=0,1}|c\rangle\langle c|
\otimes\operatorname{Tr}_Q[(A_c\otimes I_R)\rho]}
{\operatorname{Tr}(E\rho_Q)}\\
&=(\mathcal M_N\otimes{\rm id}_R)(\rho_E).
\end{aligned}
\tag{14}
$$

故

$$
D\!\left(\Omega^A_{CR}(\rho),
(\mathcal M_N\otimes{\rm id}_R)(\rho)\right)\le\delta.
\tag{15}
$$

$\Omega^A$ 的参考边缘一般不同于 $\rho_R$；这是数学重权，不是已实施的无条件装置或可向 $R$ 发送信号的新操作。

式(7)还给

$$
D\!\left((\mathcal M_{F_{\rm CO}}\otimes{\rm id}_R)(\rho),
(\mathcal M_N\otimes{\rm id}_R)(\rho)\right)\le\epsilon_{\rm CO}.
\tag{16}
$$

证明：令 $D_O=F_{\rm CO}(O)-N(O)$，
$K=\operatorname{Tr}_Q[(D_O\otimes I_R)\rho]$。差态的两块为 $K,-K$，迹距离为 $\|K\|_1$，且

$$
\|K\|_1=
\sup_{\substack{B_R=B_R^\dagger\\\|B_R\|\le1}}
|\operatorname{Tr}[\rho(D_O\otimes B_R)]|
\le\|D_O\|\le\epsilon_{\rm CO}.
\tag{17}
$$

实际二元记录＋参考与式(14)的来源统计相距至多
$\delta+\epsilon_{\rm CO}$。实际完整材料/装置没有被物理删除，只是本界考察其有限输出投影。系统测后态、允许后续、反馈和多时刻联合过程仍未认证。

## 4. 从1193质心界到固定端点概率

定义同一正来源的质心

$$
X_\rho=\int x\,q_\rho(dx).
\tag{18}
$$

它必须与1193同一 $T_B$、准备和演化产生的质心相同，不能另取方便密度替代。事前由共同准备和初始记录固定比较端点 $Y$。若1193仅逐态给出 $Y_\rho$，还须另证到固定 $Y$ 的统一界；不允许获知未知态或输出后再设 $Y=X_\rho$。

假定在声明输入族上统一有

$$
|X_\rho-Y|_\infty\le Z_x,\qquad
\int |x-X_\rho|_\infty\,q_\rho(dx)\le r_{00}.
\tag{19}
$$

若1193的完整绝对应力矩 $r$ 在同一字典中已认证，此处正性成立后有 $r_{00}\le r$。不能从另一份经典模型直接借入数值。

取事前固定 $R_0>Z_x$，并令

$$
O=\{x:|x-Y|_\infty\ge R_0\}.
\tag{20}
$$

$R_0$ 是位置容差，与参考系统 $R$ 区分；本轮统一用开球补集，包含边界。在 $O$ 上
$|x-X_\rho|_\infty\ge R_0-Z_x$。正测度Markov界给

$$
q_\rho(O)\le\frac{r_{00}}{R_0-Z_x}.
\tag{21}
$$

结合式(7)、(10)、(13)，得到

$$
\boxed{
p_{\rm CO}(O\mid\rho_{QR})
\le\min\!\left\{1,\,
\frac{r_{00}}{R_0-Z_x}
+\frac{\sqrt b-\sqrt a}{\sqrt b+\sqrt a}
+\epsilon_{\rm CO}\right\}.
}
\tag{22}
$$

它对合同统一覆盖的输入码和任意参考成立。式(15)—(17)另控制二元记录＋参考的联合统计。这是有限任务概率证书，不是确定性单次端点界或完整仪器证书。

若存在预先记录的不同准备分支，可逐支给固定 $Y_c$ 和预算，再保全部分支记录。不得事后选择成功分支，或把未知态依赖的预测中心称为同一个事前任务。

## 5. 原1112活钟：能源压缩标量不要求停钟

本节只验证原装置Hamiltonian压缩。它没有认证1193的完整空间来源测度，亦未证明其坐标能源矩等于该Hamiltonian。

采用[1112 §2—4](../../archive_1106_1119/1112/proof.md)及[原代码](../../archive_1106_1119/1112/check.py)的Material(1)。保原工作材料正基线：

$$
m_*=17+2\sqrt{19},\quad H_D=m_*I+H_{FK},\quad
H_{FK}=\Gamma(H_c\otimes I_{72})\Gamma^\dagger.
\tag{23}
$$

$H_c$ 是4维钟矩阵，对角元均为 $3/2$，相邻耦合为
$\sqrt3/2,1,\sqrt3/2$；$\operatorname{spec}(H_{FK})=\{0,1,2,3\}$。原运行编码为

$$
J_t|\psi\rangle=
\Gamma\!\left(e^{-it(m_*I+H_c)}|0\rangle_h
\otimes|a_0\rangle\otimes|\psi\rangle\otimes|{\rm vac}\rangle\right).
\tag{24}
$$

保原钟、辅助、未知载荷及接收槽。由
$\langle0|H_c|0\rangle=3/2$、
$\langle0|H_c^2|0\rangle=3$，对任意 $t$ 得

$$
J_t^\dagger H_DJ_t=E_*I,\quad
E_*=\frac{37}{2}+2\sqrt{19}\approx27.217797887081346,
\tag{25}
$$

$$
J_t^\dagger H_D^2J_t=(E_*^2+\tfrac34)I.
\tag{26}
$$

原活钟能源方差为 $3/4$，不是0；完整Hamiltonian谱宽仍为3，钟持续演化。代码压缩没有把母体能源改成常数，亦未令准备码变成定能不变子空间。原能源保持逻辑过程可按1112交织身份复用；体与外界交换时，其局部能源压缩必须重新核验。

因此，“必须让完整母Hamiltonian定能或停止钟才能使归一来源统计保持原未知输入”是过强门槛：在该码上，数学 $\rho_E=\rho$、$\delta=0$。然而物理 $\sqrt{H_D}$ 作用仍可改变完整钟态；本轮不实施它，也不把它等同于码上标量 $\sqrt{E_*}$。

原CO逻辑占据端口不必等于全部carrier或整台装置的能源质心。式(25)只消除码上的归一化障碍，未认证式(3)的空间来源、式(7)的原CO匹配或1193运动合同。

## 6. 固定有限预算

事前固定

$$
r_{00}=\frac1{10000},\quad R_0=\frac1{50},\quad
a=\left(\frac{999}{1000}\right)^2,\quad
b=\left(\frac{1001}{1000}\right)^2,\quad
\epsilon_{\rm CO}=\frac1{5000}.
\tag{27}
$$

$a,b$ 使用公开能源单位，共同比例缩放不改变 $\delta$。这是有效合同预算，不是将第5节装置的 $E_*$ 改成1；两项检验分别登记。

严格复用[1193冻结结果](../1193/results.json)中 finite_budget.position_upper.exact 的22阶含尾位置上界（1193正文式(33)给较粗但同样有效的12阶界，本轮统一采用前者）

$$
Z_x=
\frac{1726118103896327757715686428213275861058247825396308163360516479}
{143065660971439971514979501979910275072000000000000000000000000000000}.
\tag{28}
$$

于是 $\delta=1/1000$、$R_0>Z_x$，式(22)右端为

$$
p_*=\frac{1/10000}{1/50-Z_x}+\frac1{1000}+\frac1{5000}
<0.006204<\frac7{1000}.
\tag{29}
$$

约为 $0.006203018124719993$。事前阈值余量

$$
\frac7{1000}-p_*
>0.0007969.
\tag{30}
$$

可由有理算术直接复算；[check.py](check.py)、[results.json](results.json)另核码压缩、归一化、参考比较及必要删除条件。数字不代表真实仪器已经校准，也不证明全部合同已由同一量子历史同时实现。

## 7. 删除条件与结论边界

### 7.1 正总能源不能替代区域正性

例如 $E=I$，但 $A(O)=\operatorname{diag}(2,-1)$、
$A(O^c)=I-A(O)$。总积分正，区域值却可负或超过1，式(5)不是POVM，Markov推理不成立。这是代数删弱诊断，不是全认知公理或物理应力反模型。

### 7.2 不能按未知态事后定中心

两输入可分别在相距超过 $2R_0$ 的区域集中。逐态取自身质心为 $Y$，都能报告小半径，却没有给同一个预定端点的共同成功率。固定 $Y,R_0$ 与统一输入量词是本定理的组成部分。

### 7.3 效应一致不等于仪器一致

同一POVM可以对应完全不同的系统测后态。式(7)、(16)只控制此次二元记录及被动参考；不认证无扰检测、完整材料后态或全部未来。接回1192还需同历史钟束、原钟的仿射固有时和各次仪器责任。本轮始终使用1193坐标时，未以重参数化获得钟身份。

### 7.4 完整来源与暗部门保留

$A$ 必须来自原完整 $T_B$，不能由Einstein方程倒定义，也不能为取得正性删除改进项、结合能或未计控制。环境和额外暗部门按1193的 $f$ 和总账保留；若它们导致区域非正、矩增长或读口失配，应计入范围，不预先排除其存在。

本轮新增采用公理0、实际器件认证0、GR方程生成0、暗部门身份0。获得的是有限有效合同下从完整来源到端点概率的桥。下一项须实际减少来源正性、原CO匹配或原钟固有时的输入，不继续优化本组概率阈值，也不因归一化数学已处理便宣称整套共同量子几何完成。
