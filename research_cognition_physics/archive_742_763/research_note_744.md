# 第744轮：原读口的精确选择规则与联合历史的占据信号

日期：2026-10-04。接[743](research_note_743.md)及[已冻结入口](744/drafts/research_note_744_working.md)。[代码](744/joint_native_history_readout.py)、[结果](744/joint_native_history_readout_results.json)、[核验](744/research_round_744_checks.json)、[条件账](744/unified_physics_condition_ledger_744.md)。两组复算、十六式；主代理审查，无新增独立代理审查。

## 1. 问题、增量与边界

743证明原量子标量能产生超出确定二次费米演化的关联，但没有证明原占据仪器已经实现。本轮询问：**在同一原Hamiltonian和已经声明的sin s读口下，未知编码输入究竟能留下什么记录？**

解析结果是一个严格选择规则：对称准备时，单次原读口的效果只能读编码相干，任意等待都不分空态与配对态。这个规则不排除两次真实读取的联合历史。保留原五维内部几何、全部质量耦合及标量量子动能的单节点计算，在联合结果的六阶时间导数上发现稳定的占据差异。后者是数值证据，未提升为经过区间认证的存在定理。

|层次|内容|
|---|---|
|认知动机|同一过程须承载未知输入、记录、后态及资源，不能只比较已知态的均值|
|继承输入|598原有限图H、623共同域、原物种与耦合、717的sin s效果及平方根仪器|
|准备输入|解析部分取全局singlet反射对称的正常Gauss准备；数值部分另取明确的一节点Gaussian波函数|
|解析证明|原全H联合对称性；未知输入效果选择规则；真实联合历史及来源表达|
|数值验证|原一节点H的六阶时间jet；没有把标量冻结成经典背景|
|未完成|非零信号的严格符号证书、有限时间窗口、多节点推广、占据Lüders仪器、末读自治、连续及动态引力|

H⁵是内部标量目标，绝不是本轮推出的五维物理空间。空间382—386、425、522—523的条件性成果和604、649／699的边界直接保持。

## 2. 实际未知输入与完整输出

取一个节点v的原sterile模a、b，原Majorana系数Y_s。ψ是正规化、Gauss不变且在所有节点同时singlet反射下不变的玻色准备。编码等距嵌入为

$$
|0_L\rangle=\Omega,\quad
|1_L\rangle=e^{i\arg Y_s}a_v^\dagger b_v^\dagger\Omega,\qquad
V|j\rangle=\psi\otimes|j_L\rangle .
\tag{1}
$$

未知编码态ρ允许与被动参考纠缠。这里“Gauss”指规范约束，不能与后文Gaussian波函数的名称混淆。实际效果和Kraus算符完全沿用原读口：

$$
E_r=\tfrac12+\tfrac r4\sin s_v,\qquad L_r=\sqrt{E_r},\quad r=\pm1,\qquad
\tfrac14I\le E_r\le\tfrac34I,\quad \sum_r L_r^\dagger L_r=I.
\tag{2}
$$

U(t)=exp(−itH_F/ℏ)。原全物理系统作为输出保留，不投影回二能级编码：

$$
\mathcal I_{r,t}(\rho)=L_rU(t)V\rho V^\dagger U(t)^\dagger L_r,\qquad
M_r(t)=V^\dagger U(t)^\dagger E_rU(t)V .
\tag{3}
$$

每个分支完全正，求和保迹，且保Gauss；作用被动参考时张量恒等映射。准备V和最后L_r仍是声明的操作输入，这个表达没有制造自治探测装置。

## 3. 原完整H的精确对称性

S同时反射全部节点的singlet坐标，保持Higgs和规范链路。N是全CAR数算符。令

$$
\Theta=S e^{i\pi N/2},\qquad \Theta^2=(-1)^N,\qquad
\Theta H_F\Theta^\dagger=H_F .
\tag{4}
$$

证明：原目标度量、测度和动能保S；原位势依赖singlet平方，全局同时反射保原测地边势。Dirac和跳跃保持N且保S。Majorana项含x₅和两个创建／湮灭算符，两个变换各给负号，合起来不变。原自伴闭包继承该对称性。不是逐节点独立反射，也不是连续手征测度的反常结论。

准备条件和原效果给

$$
\Theta V=VZ,\qquad \Theta E_+\Theta^\dagger=E_-,\qquad
ZM_+(t)Z=I-M_+(t).
\tag{5}
$$

在编码Pauli基展开，由(2)(5)精确推出

$$
M_+(t)=\tfrac12I+b_x(t)X+b_y(t)Y,\qquad
b_x(t)^2+b_y(t)^2\le\tfrac1{16}.
\tag{6}
$$

因此两个编码占据输入的单次加号概率都等于1/2，任意t成立。没有证明两个条件后态相同，没有排除其他准备、其他观测、联合历史或完整输出中的信息。

717的实际标量力直接复用；将它压到当前编码得到相干读口的非零首项：

$$
M_+(t)=\tfrac12I-
\frac{|Y_s|}{8w}\langle\sqrt F\cos s_v\rangle_\psi\,t^2 X+O(t^3)
\quad(\hbar=1).
\tag{7}
$$

原743正常紧支撑包与w=1时，X系数约−0.05019934720188815。此处复用旧力公式，不把它计为新发现。原全部32模质量、势、边势和度量的对称校准误差为0。

## 4. 联合读取为什么不是单次结论的重复

任意真实有限历史取时间有序Kraus乘积K_r，保留全部标签及最终系统。由每段U与Θ对易、ΘL_rΘ†=L_−r得

$$
\Theta K_{\boldsymbol r}\Theta^\dagger=K_{-\boldsymbol r},\qquad
ZV^\dagger K_{\boldsymbol r}^\dagger K_{\boldsymbol r}VZ
=V^\dagger K_{-\boldsymbol r}^\dagger K_{-\boldsymbol r}V .
\tag{8}
$$

对所有标签同时翻号为奇的报告只能读X、Y；为偶的报告允许I、Z。允许Z不代表Z系数已经非零。

本轮先取：等待t，随后进行两次原平方根读取，中间等待为零。若第二个结果为q，则

$$
K_{qr}(t)=L_qL_rU(t)V,\quad
\Omega_t(\rho)=\bigoplus_{q,r}K_{qr}(t)\rho K_{qr}(t)^\dagger,\quad
\sum_{q,r}K_{qr}^\dagger K_{qr}=I.
\tag{9}
$$

保留中间反作用后，乘积标签qr的期望对应

$$
A=\sum_{q,r}qr\,E_qE_r=\tfrac14\sin^2s_v,\qquad
P_{\rm even}=\tfrac12(I+A),\qquad 0\le A\le\tfrac14I.
\tag{10}
$$

这里只在计算概率时合并对易效果。完整仪器仍是两次L的真实乘积；不能替换成sqrt(P_even)后声称后态不变。关于诱导效果与先后仪器的成熟框架，参见[Fewster—Verch第3.2—3.3节](https://arxiv.org/html/1810.06512v3)。该文不替本项目证明原耦合的实际占据信号或自治终端。

## 5. 数值对象：原全位点H，不是玻色背景近似

本节明确取一顶点、无边图。所有五个原标量坐标、全部32CAR模式与原Dirac／Majorana耦合保留；无边是图的选择，不是在多节点方程中偷偷丢弃跳跃。另取α=1的正常Gauss准备：

$$
\psi_\alpha(x)=Z_\alpha^{-1/2}e^{-\alpha|x|^2/2},\qquad
d\mu_{\mathcal K}=\frac{d^5x}{\sqrt{1+|x|^2/6}},\qquad
H=-\frac{\Delta_{\mathcal K}}{2w}+wU+\sum_{i=1}^5x_iO_i .
\tag{11}
$$

该准备不是743紧支撑波包，也不是730连续Hadamard参考。U使用623原多项式，O_i由原完整质量矩阵提取。所有有限H矩存在；没有假设任意t的Taylor级数收敛。

用稀疏多项式乘Gaussian表示H^kψ。令D=Σx_i∂_i，原Laplacian及Gaussian共轭为

$$
\Delta_{\mathcal K}=\Delta_{\mathbb R^5}+\tfrac16(D^2+4D),\qquad
\widetilde D=D-\alpha|x|^2,\qquad
\widetilde\Delta_{\mathcal K}=\Delta_{\mathbb R^5}-2\alpha D+\alpha^2|x|^2-5\alpha
+\tfrac16(\widetilde D^2+4\widetilde D).
\tag{12}
$$

对j=0,1与Ψ_j=ψα⊗|j_L>，有弱时间导数

$$
m_j^{(n)}(0)=i^n\sum_{k=0}^n(-1)^k{n\choose k}
\langle H^{n-k}\Psi_j,A H^k\Psi_j\rangle,\qquad
m_j(t)=\langle U(t)\Psi_j,A U(t)\Psi_j\rangle .
\tag{13}
$$

先完成全部五方向微分和实际CAR作用，再借Gauss不变性对标量内积采用Higgs径向代表。没有先冻结角变量。取u=αr²、v=√αx₅，求积权重和A为

$$
d\nu_\alpha\ \propto\
\frac{u e^{-u-v^2}\,du\,dv}{\sqrt{1+(u+v^2)/(6\alpha)}},\qquad
A=\tfrac14\sin^2\left(\sqrt{\frac{2}{1+(r^2+x_5^2)/6}}\,x_5\right).
\tag{14}
$$

## 6. 可复算结果及其准确强度

|Laguerre／Hermite阶数|空态六阶导数|配对态六阶导数|配对−空态|
|---|---:|---:|---:|
|48|−7.592298838157883|−8.142349609604235|−0.5500507714463518|
|64|−7.592298838261713|−8.142349609715545|−0.5500507714538321|
|80|−7.592298838255378|−8.142349609709644|−0.5500507714542655|

0至4阶差数值为零到检查精度，六阶恒等观测残差小于3.2×10⁻¹¹。H^0至H^6的稀疏项数，空态末端27089、配对末端55053；所有可达CAR态通过原作用生成。配对相位不影响这两个输入的概率。

六阶导数给六阶Taylor项的数值系数：

$$
\frac{m_1^{(6)}(0)-m_0^{(6)}(0)}{6!}\simeq-0.000763959404797591,\qquad
\frac{P_{{\rm even},1}^{(6)}(0)-P_{{\rm even},0}^{(6)}(0)}{6!}
\simeq-0.0003819797023987955.
\tag{15}
$$

不能把这些数字写成已证的有限时间判别概率：未解析证明所有低阶差消失，未认证六阶符号，未控制余项。代码有浮点运算和10⁻¹²稀疏系数舍弃；阶数收敛与恒等残差是校准，不是严格误差证书。第五阶未列入0至4阶检查，故尤其不能称六阶是已证的首个非零阶。

本轮的科学进展是：原全H的严格单次限制已确定；实际联合仪器有可复算的非零候选信号，值得继续做严格判定。既不宣称已实现理想n_f读取，也不把信号小自动解释成不能读取。

## 7. 后态、来源及接续条件

读取的完整分支作用在同一个物理状态上，任一原来源J的非选择期望必须使用实际后态。对外部几何参数γ的差分也必须保留全部因子：

$$
\langle J_\gamma\rangle_{\rm out}
=\sum_{q,r}\operatorname{Tr}\!\left(J_\gamma K_{qr,\gamma}\rho K_{qr,\gamma}^\dagger\right),\qquad
\partial_\gamma\langle J_\gamma\rangle_{\rm out}
=\sum_{q,r}\operatorname{Tr}\!\left[\partial_\gamma
(J_\gamma K_{qr,\gamma}\rho K_{qr,\gamma}^\dagger)\right].
\tag{16}
$$

若输入准备ρ也变，导数同样作用于ρ。623—625、704和718已有域、历史与二阶来源合同，本轮不重证。式(16)是来源的共同归属，不是新的Einstein解；读取注能、准备和终端实现仍不能省。

下一项[745入口](745/drafts/STATUS.md)：优先把六阶候选信号变成带可审计误差的严格符号或有限时界，或找出导致数值信号误判的结构；随后接正间隔和连通图。避免无止境提高求积阶数。原目标、连续参考、动态几何和空间既有结论保持。
