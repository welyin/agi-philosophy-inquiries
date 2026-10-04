# 第719轮：局部关系记录、内部参考与原质量反馈

日期：2026-10-04。接[718](research_note_718.md)、[719入口](719/drafts/causal_accuracy_entry.md)。[代码](719/joint_local_relational_record.py)、[结果](719/joint_local_relational_record_results.json)、[核验](719/research_round_719_checks.json)、[条件账](719/unified_physics_condition_ledger_719.md)。三组核验、二十式；主代理审查，无新增独立代理复核。

## 1. 问题、历史与对象

633及719入口已排除把原分布模Lüders读取当成瞬时因果装置。此次不重复这个反例，而检验一个具体替代：**以原中性物种的另一自旋模式承载内部共享参考，在各节点读取关系量，随后汇集经典记录。**

结果不是原非局部仪器的实现。它能在指定准备条件下重建原模式的统计量，但改变实际后态，耗去所用参考的特定相干，并由原Dirac／跳跃而非Majorana项决定即时能源反馈；Majorana仍会在等待时使系统与参考相关。

|层次|地位|
|---|---|
|认知动机|把共同参考视作同一物质过程中的子系统，测量、准备、等待和资源不能各用一份独立模型|
|继承输入|598原CAR／Gauss／完整H、原两种sterile自旋；706偶局部Kraus；716物理模式归一|
|本轮附加输入|指定系统／参考划分、参考准备、局部三结果投影仪器、内部自旋轴和共同未来中的记录汇集|
|解析新增连接|原模式相干的关系读数、实际参考消耗、同一全质量反馈及等待产生的标定修正|
|数值范围|全64模式系数；条件中性八模256维Fock演化；不是全玻色或完整Gauss Gibbs模拟|
|未完成|这些理想节点读口的自主装置、连续时空局域实现、全尺度过程和动态引力|

524已构造给定背景上的局部**标量**探针，不能将其改名为本轮费米仪器。709已有群平均与逐区删相干的限制，706已有偶局部映射；不把一般超选择、参考激活或局部非信号重记为新定理。

[Szalay等，§7.1—7.4](https://arxiv.org/html/2006.03087)说明费米局部代数与映射需要保留分次结构；[Bartlett—Rudolph—Spekkens，§III.3、VI.4](https://arxiv.org/html/quant-ph/0610030)讨论共享参考作为可消耗资源。这里借用框架，下面的物种、系数、质量响应和联合态均具体映到本项目原对象。没有将光子数U(1)超选择无条件套到含Majorana的模型。

## 2. 不加新物种的局部读取

取两个不同节点A、B。在原ν_R的两种自旋模式中指定

$$
a=\nu_{A\uparrow},\quad b=\nu_{B\uparrow},
\qquad r=\nu_{A\downarrow},\quad s=\nu_{B\downarrow},
\qquad S=(a,b),\quad R_{\rm ref}=(r,s).
\tag{1}
$$

它们都是原规范singlet，但仍满足CAR；不能把跨节点奇算符当成交换变量。选定局部内部轴的关系观测为

$$
X_A=a^\dagger r+r^\dagger a,\qquad
X_B=b^\dagger s+s^\dagger b,\qquad
[X_A,X_B]=0,\quad X_v^3=X_v .
\tag{2}
$$

每个X都为局部偶Gauss算符，其三个谱投影为

$$
P_{v,\pm}=\frac{X_v^2\pm X_v}{2},\qquad P_{v,0}=I-X_v^2,
\quad {\cal L}_v^*(O)=\sum_{\ell=-1}^1P_{v,\ell}OP_{v,\ell}.
\tag{3}
$$

实际联合结果为两个三值标签，K_{ij}=P_{B,j}P_{A,i}。所有结果及其后态均保留，不以诱导效果的平方根替代联合仪器。

任意B区偶观测O_B与A区投影对易，所以

$$
{\cal L}_A^*(O_B)=O_B,\qquad
{\cal L}_A^*{\cal L}_B^*={\cal L}_B^*{\cal L}_A^* .
\tag{4}
$$

这证明给定节点划分上的操作不从A向B发送信号；读数相关性须在真实通信后的共同未来计算。它不证明每个节点内的连续分布模式已能瞬时读取，也不提供装置、时钟、开关和传播光锥。若有有限时间等待，仍用原H，不能冻结跨边相互作用来假称完整因果实现。

## 3. 共享参考使哪些相干可读

先明确一次准备合同：在S与R_ref之间取分次乘积的偶密度算符ρ_S、τ；τ位于参考单粒子部门。记

$$
z_S=\langle a^\dagger b\rangle_{\rho_S},\quad
z_R=\langle r^\dagger s\rangle_\tau,\qquad
\tau_\eta=\frac12
 \begin{pmatrix}1&\eta\\ \eta&1\end{pmatrix}_{|r\rangle,|s\rangle},
\quad 0\le\eta\le1 .
\tag{5}
$$

τ_η是明确的初始共享量子资源，z_R=η/2；没有要求每个节点存在奇宇称期望值。全系统可取固定偶总宇称：S单粒子态与参考单粒子态组成二粒子态。乘以原规范不变光滑紧支撑玻色包即给正常有限能源Gauss准备的一个存在例；准备方法及其成本仍未由此求出。

CAR重排给

$$
\langle X_A X_B\rangle_{\rho_S\otimes_g\tau}
=-2\operatorname{Re}(z_S\overline{z_R})
=-\eta\,\operatorname{Re}z_S .
\tag{6}
$$

参考单粒子假设使异常配对期望为0。全式的负号来自费米重排，普通四qubit张量积会漏掉它。

原目标模式c=√p a+√q b，p+q=1。它的平均占据为

$$
\langle n_c\rangle
=p\langle n_a\rangle+q\langle n_b\rangle
+2\sqrt{pq}\operatorname{Re}z_S
=p\langle n_a\rangle+q\langle n_b\rangle
-\frac{2\sqrt{pq}}{\eta}\langle X_A X_B\rangle,\quad\eta>0 .
\tag{7}
$$

计数项在同分布的另批准备上读取；复相干的另一分量可用指定相位轴，不能同时测量不对易轴却保持原未知样本。η变小使统计重建对误差敏感；本轮不优化资源。

**式(7)是多次准备上的均值重建，不是单副本二值Lüders读取。** 重建权重可以为负或超过1，因此没有自动得到每次试验的合法二值输出概率，更没有原K_σ后态或718的原联合历史。

没有共享参考时，如果辅助态分别对各区域宇称不变，允许的操作仅为局部偶操作和经典通信，则有效任务对局部宇称平均不敏感。用709的群平均方法，两种原相位态满足

$$
{\cal T}_A(|\psi_+\rangle\langle\psi_+|)
={\cal T}_A(|\psi_-\rangle\langle\psi_-|),\quad
|\psi_\pm\rangle=\frac{a^\dagger\pm b^\dagger}{\sqrt2}|0\rangle,\quad
{\cal T}_A(\rho)=\tfrac12(\rho+\Pi_A\rho\Pi_A).
\tag{8}
$$

在上述操作菜单中它们不能区分。共享τ_η不再分别对参考A／B宇称平均不变，故关系记录可恢复部分区分能力。这个菜单限制不是原H的超选择定理：实际跨边量子传播、额外共享资源均会改变前提。

## 4. 同一次读取会改变参考

X_A的零谱子空间是节点(a,r)的偶宇称空间，±1谱子空间各为奇宇称。r或r†连接这两类，因而

$$
P_{A,i}rP_{A,i}=0,\qquad
P_{A,i}r^\dagger P_{A,i}=0,\qquad
{\cal L}_A^*(r^\dagger s)=0 .
\tag{9}
$$

故对任意输入，不论选择还是非选择，每个有非零概率的实际联合结果之后，参考的这份相干都消失：

$$
\langle r^\dagger s\rangle_{\rm after,\;i,j}=0,
\qquad
\langle X_A X_B\rangle_{\rm after}
=\langle X_A X_B\rangle_{\rm before}.
\tag{10}
$$

关系信息转入记录及联合状态，不能把它说成所有量子资源都消失。该τ_η不能作为原样不变的、免费重复使用的参考；恢复、替换或保留更大联合记忆须另算。这里证明的是特定相干消失，不是给出最小热力学擦除功。

在数值p=3/4、η=1时，ψ±分别给〈X_AX_B〉=∓1/2；两边B边缘同为(1/4,1/2,1/4)，共同记录的总变差距离1/2。式(7)重建0.933012701892及0.066987298108。η=.4时联合对比降为.2；η=0时两个联合分布完全相同。每个分支及非选择参考的上述相干均为0。

## 5. 原完整H给出的反馈：Majorana保留，Dirac和跨边项改变

本仪器为配置无关的有限CAR多项式，与原H_b、T_*、W_*对易，保Gauss及623共同域。对X_v的群平均与谱夹断等价：

$$
{\cal L}_v^*(O)=\frac1{2\pi}\int_0^{2\pi}
 e^{i\theta X_v}Oe^{-i\theta X_v}\,d\theta,\qquad
e^{i\theta X_v}\nu_v e^{-i\theta X_v}
=e^{-i\theta\sigma_x}\nu_v .
\tag{11}
$$

原sterile Majorana配对是自旋反对称标量。局部SU(2)满足det=1，所以原κ_v ν_v↑†ν_v↓†及伴随与X_v对易。相反，原Dirac项把ν_v连到未一起旋转的L_v，单个ν算符的平均为0；两个不同节点独立平均同样删掉与读取节点相接的sterile跳跃。

对读取节点集S₀和原声明的同种自旋标量跳跃，完整差量精确为

$$
D_{\rm rel}:={\cal L}_{S_0}^*(H)-H
=-H_{\nu L,S_0}
-H_{\nu{\rm hop},\,e\cap S_0\ne\varnothing},
\qquad {\cal L}_{S_0}^*(H_{\rm Majorana})=H_{\rm Majorana}.
\tag{12}
$$

没有删原动力学，只计算实际非选择读取的状态反馈。若扩展模型另有ν自旋耦合、接触项或来源依赖轴，须重算，不能把(12)当普遍规范原则。

这里用ℋ_v表示原Higgs双重态，以免与本轮关系读口X_v混淆。每个原Dirac块有两个相同奇异值m_v=|Y_ν||ℋ_v|/√F_v。其二次Fock算符范数为2m_v；每条sterile自旋标量边同理给2|t_e|。于是

$$
\|D_{\rm rel}(\phi,g)\|
\le 2|Y_\nu|\sum_{v\in S_0}\frac{|\mathcal H_v|}{\sqrt{F_v}}
+2\sum_{e\cap S_0\ne\varnothing}|t_e|.
\tag{13}
$$

用598／716的F⁻²≤AU+B₀以及|ℋ|≤√(6M)，正常有限能源态有

$$
|\Delta E|
\le 2|Y_\nu|\sqrt{6M}\sum_{v\in S_0}
 (A\langle U_v\rangle+B_0)^{1/4}
+2\sum_{e\cap S_0\ne\varnothing}|t_e| .
\tag{14}
$$

固定图下旧形式下界保证右侧有限。ΔE可以正或负，必须连同装置／控制账处理，不能都称为“消耗”。原完整Gibbs的非负性复用714的被动性；本轮不是能量守恒装置的设计证明。

局部轴固定且不随背景λ改变时，仪器本身来源导数为0；完整准备响应仍须保留：

$$
\partial_\lambda\langle D_{\rm rel}\rangle_{\rho_\beta}
=\langle{\cal L}^*(G)-G\rangle_{\rho_\beta}
-\beta\,\operatorname{Cov}_{\rm KM,\rho_\beta}(G,D_{\rm rel}),
\qquad G=\partial_\lambda H.
\tag{15}
$$

统计重建的p、q和参考η若随几何变化，另按(7)求导；不能因本仪器轴固定就省略716／718的任务和准备变化。能量式中没有即时Majorana差量，不代表Majorana与整个参考过程无关。

## 6. 同一原演化不保持免费的独立参考

对任意联合态定义z_S、z_R及α_S=〈ab〉、α_R=〈rs〉。CAR展开为

$$
X_AX_B=-a^\dagger b^\dagger rs-a^\dagger b\,s^\dagger r
-b^\dagger a\,r^\dagger s-ab\,r^\dagger s^\dagger .
\tag{16}
$$

因此独立准备时的因子化形式及其一般修正为

$$
\langle X_AX_B\rangle
=2\operatorname{Re}(\overline{\alpha_S}\alpha_R-z_S\overline{z_R})
+{\cal C}_{SR},
\quad
{\cal C}_{SR}:=\langle X_AX_B\rangle
-\langle X_AX_B\rangle_{\rho_S\otimes_g\rho_R}.
\tag{17}
$$

C_SR是同一联合态的关联，不是新自由参数。任意态不能靠两个边缘确定它；数值条件二次演化的Gaussian支可由完整协方差复算，不能直接推广到动态玻色耦合后的所有态。

原Majorana已经给出解析的独立参考破坏机制。以原CAR真空乘正常Gauss玻色包χ为初态，κ_A=Y_s s_A/√F_A；在共同核上

$$
\left.\frac{d}{dt}\langle ar\rangle\right|_{0}
=\frac{i}{\hbar}\langle\chi,\kappa_A\chi\rangle .
\tag{18}
$$

选择s_A正侧支撑使其非零。原Dirac和跳跃保粒子数，对真空的这一一阶项无贡献；H_b也与ar对易。任意偶边缘乘积的〈ar〉为0，所以完整原H在小时间内已经把S与参考关联起来。此证明没有关掉其它原项，也没有把配置点态当正常物理准备。

复用718两节点条件中性八模，保全部Dirac、Majorana和跳跃，对声明的初始两粒子乘积态得到：

|t|实际〈X_AX_B〉|由当时两个边缘因子化预测|C_SR|
|---|---:|---:|---:|
|0|−.1303765618|−.1303765618|约0|
|.6|−.0889281224|−.0889659764|.0000378539|
|1.7|−.0030538552|−.0083220591|.0052682040|
|3.1|−.0086239166|−.0138760725|.0052521559|

即使不断更新两个边缘，也不能假装重新获得最初独立参考。三角界|C_SR|≤‖ρ_SR−ρ_S⊗_gρ_R‖₁是成熟估计，不提供该关联的自动闭合方程。

## 7. 核验与下一项

三组检查对应：完整局部记录和共享参考；原质量／跳跃／几何来源；同一原等待下的参考关联。全64模式六组非平凡规范／Higgs配置上，

$$
\max\|h_{\rm measured}-h_{\rm predicted}\|<6.07\times10^{-16},
\qquad
\max\|\Delta_{\rm measured}-\Delta_{\rm original}\|<4.05\times10^{-17}.
\tag{19}
$$

256维条件Gibbs β=.73中，沿718声明的跳跃背景路径：

$$
\Delta E_{\rm cond}=0.0949215073214105,\qquad
\partial_\lambda\Delta E_{\rm cond}=0.0507157391824900.
\tag{20}
$$

末档中央差分误差7.06×10⁻¹¹。矩阵计算只校准原条件部门，不作为完整热参考或连续因果过程的数值模拟。旧空间382—386、425、522—524及699限定反例保留，384已消去的Lipschitz不恢复。

本轮合并C03记录、C15质量、C19参考、C22来源：可读取的统计量、可保留的后态、参考的消耗和同一H中的等待不能分别签收。排除了“以局部关系统计替代整个原Lüders历史”以及“参考在原等待和读取后仍免费不变”的接法；保留实际量子传输、更多联合记忆和不同任务作为替代。

接[720](720/drafts/STATUS.md)：核原质量的自旋结构能否给可保留的联合参考，以及这种结构是否兼容已有空间传播／手征字典；先回查604、614、633、666—667与706，不从新仪器优化转入装置设计，也不因局部SU(2)代数就宣称物理旋转对称或3+1已生成。
