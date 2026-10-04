# 第633轮：连续记录、态正则性与同一物质的几何来源

日期：2026-10-01。接[632](research_note_632.md)，回查[524](../archive_467_530/research_note_524.md)、[592](../archive_585_628/research_note_592.md)、[598](../archive_585_628/research_note_598.md)、[602](../archive_585_628/research_note_602.md)及[624](../archive_585_628/research_note_624.md)。[代码](633/joint_continuum_record_sources.py)、[结果](633/joint_continuum_record_sources_results.json)、[核验](633/research_round_633_checks.json)、[条件账](633/unified_physics_condition_ledger_633.md)。四组检查、十八式；主代理审查，无新增独立代理审查。

## 1. 本轮合并什么

630—632连接了同一连续真空的谱、应力与局部项，尚未说明发生记录后是否仍有可用的连续态和来源。本轮在**同一原物质的自由、恒定背景分支**中，定义原sterile场的一个光滑局域模，保留原Dirac／Majorana混合，联立计算二值记录、后态的短距离正则性、有限激发能、质量来源和应力迹。

一般CAR、CP和Hadamard封闭性质均属成熟工具；新增的是原混合质量和当前几何来源上的具体共同映射。524已用Fewster—Verch方法研究另一标量探针；592／624已处理原有限图仪器。这里不把这些结果重新编号，也不宣称已把原图仪器映入连续理论。

|层次|地位|
|---|---|
|认知动机|发生记录后，态、物质与几何来源仍应来自同一过程|
|继承对象|598一整代物质、原复Y、原F；632树背景与同一固定真空减除|
|附加输入|给定平直3+1、canonical自由费米场、光滑Cauchy涂抹；声明的二值Lüders操作|
|解析连接|后态的平滑两点差、各阶能源矩、原混合谱所固定的注能和来源|
|明确冲突|整个分布模的瞬时读取可在其两个子区域之间传信，不能直接充当有限传播的真实过程|
|未完成|自主装置、能源供给、有限时因果实现、相互作用态、原图映射及动态引力|

局域代数归属与装置可因果实现分开验收。后者是条件，不在本轮转入认知硬件设计。

## 2. 成熟结果与原对象映射

[Sanders，The locally covariant Dirac field](https://arxiv.org/html/0911.1304)，命题4.8给局域可观测代数的因果性和time-slice性质，命题4.12讨论Hadamard态对代数操作的封闭；§4.3连接应力与几何变化。本文只借用这些自由场工具。原有限个混合中性自由场可用恒定Takagi变换化为两个Majorana场；下面也直接证明本特定光滑操作的正则性，不将该文当成相互作用SM或装置定理。

取632原树值，场排列为(ν_L↑,ν_L↓,ν_R↑,ν_R↓)，保留原代码的复质量与配对块：

$$
\phi_\star=(0,\sqrt{u_1},0,0,\sqrt{u_2}),\quad
B(k)=\begin{pmatrix}K(k)+h&\Delta\\
\Delta^\dagger&K(k)^T-h^T\end{pmatrix},\quad
K(k)=\operatorname{diag}(-\sigma\cdot k,\sigma\cdot k).
\tag{1}
$$

h、Δ直接取598的索引(24,25,30,31)。下方动能符号来自Nambu中的c†(−k)，不能误用−K(k)ᵀ。B(−k)与B(k)满足粒子—空穴关系；B₀=B(0)与动能反对易。其余带电物种仍在同一自由理论中，当前操作对它们为恒等，没有删去物种或改变反常计数。

令D=|Yν|√u₁/√F，M_R=|Ys|√u₂/√F。sterile原场在两个质量部门的权重为

$$
m_\pm=\frac{\sqrt{M_R^2+4D^2}\pm M_R}{2},\qquad
w_\pm=\frac12\left(1\pm\frac{M_R}{\sqrt{M_R^2+4D^2}}\right),\quad
w_-+w_+=1.
\tag{2}
$$

这不是把原sterile场替换成独立质量本征场；实际读取仍用原场，w保留混合。数值m₋=0.148500538672、m₊=0.277129738518，w₋=0.348895618170、w₊=0.651104381830。原8维BdG对24组动量的独立核对最大误差2.67×10⁻¹⁵。

## 3. 同一记录对态的要求

取固定Cauchy面上的f∈C∞₀(R³,C²)，∥f∥₂=1，涂抹原ν_R。由原规范表示它为G-singlet；虽c是奇场，下列仪器为偶：

$$
c=\nu_R(f),\quad \{c,c^\dagger\}=1,\quad c^2=0,\qquad
n=c^\dagger c=n^2,\quad L_1=n,\quad L_0=1-n.
\tag{3}
$$

空间支持位于一个有界区域时，借time-slice可把它放入包住该数据的时空区域代数。它没有读取一个外加无质量探针。对原自由真空Ω，令p=〈n〉；0<p<1时

$$
p_1=p,\quad p_0=1-p,\qquad
\omega_r(A)=\frac{\langle\Omega,L_r A L_r\Omega\rangle}{p_r},\qquad
\bar\omega(A)=\tfrac12\{\omega(A)+\omega(R A R)\},\quad R=1-2n.
\tag{4}
$$

CP及归一来自投影代数。与c、c†都分离的偶可观测A与L_r对易，所以其非选择期望不变；选择性条件期望可以改变原关联，这不是可控的远程传信。

自由真空作Wick展开时，插入两个L_r后的两点函数，除原两点函数外仅含有限多个经f涂抹的真空收缩之积：

$$
S_r(x,y)-S_\Omega(x,y)
 =\sum_{a,b=1}^{N_r}d^{(r)}_{ab}u_a(x)v_b(y),\qquad
S_r-S_\Omega\in C^\infty.
\tag{5}
$$

理由是原恒定质量自由算子的正负频投影为有界符号矩阵，乘上f的快速下降Fourier变换仍具有全部加权L²矩；自由演化保这些矩，Sobolev嵌入使各收缩及其时空导数光滑。Wick项数有限，p_r>0仅带来有限归一因子。因此后态保持同一Hadamard奇性。没有假设后态仍是真空或一般Gaussian态。

以同一真空正规序定义H_vac≥0。L_rΩ至多含真空与两个自由准粒子，各波函数有全部能源矩，故

$$
L_r\Omega\in\bigcap_{j\geq1}D(H_{\rm vac}^j),\qquad
\Delta\langle T_{\mu\nu}(x)\rangle_r
 =\left.\mathcal D_{\mu\nu}(S_r-S_\Omega)(x,y)\right|_{y=x}.
\tag{6}
$$

右侧有限且光滑；质量来源同理。这是本光滑操作的充分条件，不能泛化为任何有界局域操作都具有有限能源。有限的是相对于同一固定减除的激发能，不是未经重整化的绝对真空能。

## 4. 记录、注能与共同来源的实际计算

令e为Nambu空间中f占据原sterile粒子位置的向量，C为反线性粒子—空穴映射，P₋=1_{B<0}。Q投影到e、Ce张成的二维子空间；R在一粒子空间为1−2Q。非选择协方差差为

$$
Q=|e\rangle\langle e|+|Ce\rangle\langle Ce|,\qquad
\Delta S=-QP_--P_-Q+2QP_-Q,\qquad
\Delta\langle J\rangle=\tfrac12\operatorname{Tr}(J\,\Delta S).
\tag{7}
$$

J为同一二次可观测的BdG来源矩阵；常数项在状态差中抵消，1/2是Nambu重复计数。ΔS秩至多4，且其核光滑。能量、标量和几何来源使用同一个ΔS，不能各选一个“读后参考态”。

一般f下，设e₊=〈e,BP₊e〉、e₋=−〈e,BP₋e〉。两者非负。把nΩ写成pΩ+χ₂，由CAR，两个准粒子波函数正交，范数平方为p与1−p，于是

$$
\|\chi_2\|^2=p(1-p),\quad
\mathcal E_\chi=p e_+ +(1-p)e_-,\quad
E_1=\frac{\mathcal E_\chi}{p},\quad
E_0=\frac{\mathcal E_\chi}{1-p},\quad
\bar E=2\mathcal E_\chi .
\tag{8}
$$

这里H_vacΩ=0使真空—双粒子交叉能源项消失；其他来源的选择性期望一般仍有交叉项，不应套用E_r公式。

取实径向包乘固定spinor。角平均使各质量的σ·k项为零，p=1/2。ν_R本身不是一个质量本征场，其完整谱给

$$
p=\tfrac12,\qquad E_0=E_1=\bar E
 =\sum_{\alpha=\pm}w_\alpha\int d^3k\,|\widehat f(k)|^2
       \sqrt{k^2+m_\alpha^2}.
\tag{9}
$$

以下量是理想更新后同一非选择态的全空间积分，不是局部密度恒定或已经求出的引力场。O是均匀质量缩放源，P_av是空间应力迹的三分之一：

$$
\Delta O=\sum_\alpha w_\alpha\int d^3k\,|\widehat f|^2
          \frac{m_\alpha^2}{\sqrt{k^2+m_\alpha^2}},\qquad
\Delta P_{\rm av}=\frac13\sum_\alpha w_\alpha\int d^3k\,|\widehat f|^2
          \frac{k^2}{\sqrt{k^2+m_\alpha^2}},\qquad
\bar E-3\Delta P_{\rm av}=\Delta O .
\tag{10}
$$

式(7)用J=B₀及J=B−B₀分别得到(10)；没有把“记录能源”直接宣称为独立的几何项。在同一径向q参数化上，m′α=bmα，因此

$$
\Delta J_q=b\,\Delta O,\qquad
b=\frac1{\sqrt6}\coth\frac{q_\star}{\sqrt6}
  =1.64312245702040,\qquad
\int d^3x\,\Delta T^\mu{}_\mu=-\Delta O
\quad(\eta=(-,+,+,+)).
\tag{11}
$$

本轮在632原树背景计算，所以b不同于630的q=0.27诊断点。这一差异已由原F与场值决定。固定背景上的态无关异常项和局部减除在状态差中抵消；不由(11)推出完整曲背景Ward身份。

632保持平直树背景的有限匹配系数继续固定。若记录之后又重选一个真空或调一个常数来抹去(9)，便更换了模型。闭合整体且开关前后相互作用能为零时，装置与控制部门至少须共同满足

$$
\Delta E_{\rm field}+\Delta E_{\rm apparatus+control}=0.
\tag{12}
$$

本轮只给ΔE_field，没有完成供能、记录稳定性或全部熵成本；它们不能因有限注能就被签收。

## 5. 真正紧支撑的连续样本与误差范围

没有用Gaussian的近似局域代替紧支撑。本数值样本为

$$
f_\ell(x)=\ell^{-3/2}N
 \begin{cases}\exp[-1/(1-|x/\ell|^2)],&|x|<\ell,\\0,&|x|\geq\ell,\end{cases}
 \ \chi_\uparrow,\qquad
\int|f_\ell|^2=1.
\tag{13}
$$

用单位化Fourier约定进行径向积分。基础包N=3.22576097031，K₂=∥∇f₁∥²≈12.3017353737，K₄=∥Δf₁∥²≈543.338856442。原谱给解析界

$$
\sum_\alpha w_\alpha m_\alpha\leq\bar E_\ell
 \leq\sqrt{K_2/\ell^2+D^2+M_R^2},\qquad
E_{\text{tail},\,|q|>Q}
 \leq\frac{K_4}{\ell Q^3}
       \sqrt{1+\frac{\ell^2m_+^2}{Q^2}} .
\tag{14}
$$

尾界由q⁴加权范数控制，存在性证明不使用截断。数值在q≤180积分；尾界系数K₄由数值求得，未作区间认证，所以不把网格差声称为严格总误差。加密两组求积的变化≤8.44×10⁻¹⁴；第一行解析尾界的数值估计为9.32×10⁻⁵，明显比网格差保守。

|支持半径ℓ|E₀=E₁=平均注能|质量源ΔO|径向源ΔJq|平均积分压力|
|---|---:|---:|---:|---:|
|1|3.141211148765|0.023200355096|0.038121024470|1.039336931223|
|2|1.587654563428|0.044648912625|0.073363631015|0.514335216935|
|4|0.825609486425|0.081033493319|0.133147952644|0.248191997702|

原参数单位为模型诊断单位，不是实测质量、米或焦耳。三行只是同一连接的尺度样本，不是寻找最优仪器。原四中性CAR模式的16维Fock核对给协方差／来源最大差5.89×10⁻¹⁶；原s来源在零动量样本增加0.0824734162447。有限Fock测试只核代数与原质量，连续正则性由第3节证明。

## 6. 不能把局域代数中的投影直接当成瞬时物理记录

[Fewster—Verch，Quantum fields and local measurements](https://arxiv.org/html/1810.06512)，§3以时空局域耦合定义仪器，并检验因果有序组合。524已使用其方法。我们继承这一要求，给当前原sterile模的一个具体排除例，不把一般“测量须因果”当新定理。

设a=ν_R(f_A)、b=ν_R(f_B)，f_A、f_B光滑且空间支持分离。它们属于同一原场，满足两个独立CAR模式的代数。令

$$
c=\cos\theta\,a+\sin\theta\,b,\quad n_c=c^\dagger c,\quad
\Phi^*(A)=n_c A n_c+(1-n_c)A(1-n_c),\qquad \theta=\pi/6.
\tag{15}
$$

两区域的偶可观测互相对易，Φ也保持整个支持外的观测。但其对B区占据数的作用包含跨区相干：

$$
\Phi^*(n_b)=2\cos^2\theta\sin^2\theta\,n_a+
 (1-2\cos^2\theta\sin^2\theta)n_b
 -\cos\theta\sin\theta\cos(2\theta)(a^\dagger b+b^\dagger a).
\tag{16}
$$

选两个局域CAR模式中固定奇宇称的单粒子态，局域偶相位U_A=1−2n_a把二者互换至无关总相位：

$$
|\psi_\pm\rangle=(a^\dagger\pm b^\dagger)|00\rangle/\sqrt2,\qquad
\langle n_b\rangle_+=\langle n_b\rangle_-=\tfrac12,\qquad
\langle\Phi^*(n_b)\rangle_\pm=\tfrac12\mp\frac{\sqrt3}{8}.
\tag{17}
$$

故若把整个Φ安排为相邻瞬间、允许B区在A区信号到达前读取n_b，A的可控相位产生√3/4≈0.433012701892的可读差异。该反例针对**无传播时长的分布式实现**；完成有限时跨区通信后实施一个联合操作不受此结论排除。

本反例初态不是前面计算的真空。可从真空先取(1−n_a)(1−n_b)分支，再用(a†±b†)生成这些模式态，复用第3节的光滑多项式论证。为直接核清分支非零，取两个分离平移的相同实径向包；p_a=p_b=1/2。Wick给空分支概率为1/4−|〈a†b〉|²+|〈ab〉|²。原真空的动量占据权在开集上为正，两个平移包的相位不是常数，严格Cauchy—Schwarz给|〈a†b〉|<1/2，故空分支严格非零。这样得到同一场中的Hadamard、各阶有限能源反例态；没有仅凭有限维态自动宣称连续态可用。该准备不是免费或瞬时的，所排除的是随后Φ的超光速实现。有限CAR矩阵只核(16)—(17)，未模拟真实探针。

因此目前能共同签收与仍独立的条件为

$$
\text{光滑原场记录}\ \Longrightarrow\
 \{\text{正态、同一短距结构、有限激发能、共同来源}\},
\qquad
\text{上述条件}\ \not\Longrightarrow\
 \text{分布式仪器可瞬时因果实现}.
\tag{18}
$$

## 7. 条件账与下一步

本轮把C03／C15／C17／C19／C22的一个连续自由场接口接通：同一原混合质量决定记录后的状态及物质／几何来源。原真空的局部匹配不能在每次记录后独立重置。

C02／C05又共同限制C03：局域代数归属、正态和有限能源仍不足以签收真实时空中的仪器组合。新的具体反例排除当前瞬时分布模读取接法，保留有限时因果过程。这不是整个统一计划或全部局域测量的反证。

[634入口](634/drafts/STATUS.md)继续整合记录、连续来源与几何约束：优先区分逐记录条件期望和无条件物理来源，核跨记录界面的守恒及共同状态边界条件。复用600—601、620、624及524的真实过程条件；不继续优化本波包或设计认知硬件。原图映射、相互作用、一般曲背景、内部供能及GR仍开放。
