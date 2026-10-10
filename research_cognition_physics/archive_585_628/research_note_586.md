# 第586轮：共同记录的压缩、后继量子过程与同源几何响应

日期：2026-10-01。接[正向成果回顾](585/cognition_forward_bridge_review_585.md)及[585](research_note_585.md)。五组复算，主代理推导及核验；未取得新的独立代理审查。[代码](586/joint_record_source_compression.py)、[结果](586/joint_record_source_compression_results.json)、[正式核验](586/research_round_586_checks.json)。

## 1. 要检验的共同映射

**问题：** 当前物质模型中的若干读取只保留一个共同摘要时，能否把它们替换为“直接读取该摘要”，同时保持后继物质过程与几何来源？

**结论：一般不能。** 对577原二元读取连续使用两次，只报告两个结果的奇偶；与直接使用同一奇偶效果的平方根仪器相比，它们在所有未知输入上的即时报告概率完全相同，却有不同的完整后态。在原574完整Hamiltonian下，原有二元读口能在随后充分小的正时间区分二者；同时局部能源和空间权响应不同。无需新增物质字段，也未改原Hamiltonian。

还给出一个限定的正向判据：在本轮的正、光滑、配置乘法平方根仪器类中，**被省略的细记录在给定粗记录后不再依赖配置**，恰是粗化后仪器等于直接粗读、并对所有允许来源保持能源的条件。保留完整粗化CP映射则始终合法；不要求所有摘要都使用平方根仪器。

|层次|本轮地位|
|---|---|
|认知动机|共同认可的摘要应能够继续使用；不能省去仍改变后继过程的关系或资源。|
|继承|574原曲目标、物质、规范群、图及量子化；577末读；584正规Gauss包；585局部来源定义。|
|额外操作输入|理想末读可以有限次组合，结果可分组并隔离；直接粗读是待比较的另一种指定仪器。|
|解析增量|原物质中的粗记录—后态—能源相容判据；原二次读取的严格Gauss反例；同H有限时间可读见证与几何响应。|
|未完成|原H自主制备仪器、存档装置与内部控制，宏观尺度映射、量子引力约束及连续极限。|

### 去重及成熟工具

[506](../archive_467_530/research_note_506.md)已证明角色细读与集体读的过程不同；[557](../archive_554_584/research_note_557.md)已有指定Gaussian仪器的注能形式；[577](../archive_554_584/research_note_577.md)已给原末读的Gauss和形式域合法性。它们直接复用，不将一般“POVM不能决定后态”当新发现。本轮新增的是**当前完整曲目标物质中，同一粗记录的后继读数与几何来源不能共同替换，以及替换成功的明确条件**。

[Gudder，§2、Theorem 2.6](https://arxiv.org/pdf/2109.07019)区分仪器、其测量的效果和后处理。[Chiribella–Yang，coherent coarse-graining定义](https://hub.hku.hk/bitstream/10722/245818/1/content.pdf)明确讨论将一个操作换为其效果平方根操作；这与普通结果合并不同。本轮直接给所需核函数证明，不借此假定新仪器的物理实现。[Radaelli等，Appendix B](https://arxiv.org/html/2206.00463)给经典Fisher信息的条件分解；下文采用配置依赖概率的有限和版本，不把它当量子Fisher信息的一般链式法则。

## 2. 固定同一个物质对象

不改574的物理空间与自伴算符。记节点B的singlet为s，Higgs长度为h，M=2；全图配置空间及节点动能为

$$
\mathcal Q=(\mathbb H^5)^V\times G^E,\quad
\mathcal H_{\rm phys}=P_{\mathcal G}L^2(\mathcal Q,d\mu),\quad
H=-\frac{\hbar^2}{2}\Delta_{\mathsf G}+W,\quad
T_B=-\frac{\hbar^2}{2w_B}\Delta_{\mathcal K},\quad
F=M-\frac{h^2+s^2}{6},\quad
\mathcal K^{-1}=F\left(I-\frac{\phi\phi^T}{6M}\right),\quad
w_B=\epsilon^3\psi_B^6.
\tag{1}
$$

原W保留全部现场势、测地边势和规范磁势；链路动能也保留。几何权仍固定，Gauss是内部规范约束。B指原节点及其访问代数，不假定物理Hilbert空间按规范区域自由分解。

先取节点B上的有限组光滑、严格正、规范不变配置效果，和式恒为1；为取得所有有限能量态的统一形式域结论，要求各效果有正常数下界、曲度量梯度一致有界：

$$
e_r(\phi_B)>0,\qquad \sum_r e_r=1,\qquad
L_r=\sqrt{e_r},\qquad
\mathcal I_r(\rho)=L_r\rho L_r,\qquad
\sum_r\mathcal I_r\text{ 为CPTP通道},\qquad
[L_r,P_{\mathcal G}]=0.
\tag{2}
$$

任意旧参考R上作恒等扩张，所有分支仍保Gauss。每个结果的平方根是本轮指定的仪器，不是由效果唯一推出。

## 3. 一个形式恒等式把记录与原动能接起来

定义配置参数的经典Fisher张量$I_{ab}=\sum_r(\partial_a e_r)(\partial_b e_r)/e_r$。它描述这份仪器的配置分辨性，不是已涌现的空间度规，也不是未知量子态的Fisher张量。

在原共同形式核上展开$d(L_r\Psi)$，由$\sum L_r^2=1$、$\sum L_r dL_r=0$得到

$$
\Delta E_{\mathcal I}(\rho)
=\sum_r\operatorname{Tr}(H L_r\rho L_r)-\operatorname{Tr}(H\rho)
=\frac{\hbar^2}{2w_B}
\left\langle\sum_r\mathcal K^{-1}(dL_r,dL_r)\right\rangle_\rho
=\frac{\hbar^2}{8w_B}\langle \mathcal K^{ab}I_{ab}\rangle_\rho.
\tag{3}
$$

交叉项在求和后精确消失，即使原态具有非零概率流或与参考关联也成立。所有配置乘法势及其它因子的动能平均不变。梯度界使各L_r保持原能量形式域；由闭形式及正常态分解，式(3)延拓到任意有限平均能源的物理态，无需假设属于D(H²)。这里只是**物质源注能**，不等于装置总耗散、热力学最小功或所有仪器的普适下界。

## 4. 粗记录相同，何时完整过程也相同

按固定分组a合并细结果r，记$p_a=\sum_{r\in a}e_r$、$q_{r|a}=e_r/p_a$。对这些正函数，有限和求导直接给

$$
I^{\rm fine}_{ij}-I^{\rm coarse}_{ij}
=\sum_a p_a\sum_{r\in a}
\frac{\partial_iq_{r|a}\,\partial_jq_{r|a}}{q_{r|a}}
\succeq0.
\tag{4}
$$

证明：写$de_r=q_{r|a}dp_a+p_a dq_{r|a}$代入；两交叉项因$\sum_r dq_{r|a}=0$消失。这不是对实际量子态的经典隐变量假设，只是对配置乘法效果的函数恒等式。

两种具有同一报告概率的仪器分别是

$$
\mathcal J_a^{\rm fine}(\rho)=\sum_{r\in a}L_r\rho L_r,\qquad
\mathcal J_a^{\rm direct}(\rho)=\sqrt{p_a}\rho\sqrt{p_a},\qquad
\operatorname{Tr}\mathcal J_a^{\rm fine}(\rho)
=\operatorname{Tr}\mathcal J_a^{\rm direct}(\rho)
=\operatorname{Tr}(p_a\rho).
\tag{5}
$$

在本轮连通配置域、严格正光滑效果及完整Gauss来源域内，以下三项等价：全部粗分支通道相同；所有条件概率q在配置上为常数；两种非选择仪器对每个允许有限能源来源有相同注能。

$$
\mathcal J_a^{\rm fine}=\mathcal J_a^{\rm direct}\ (\forall a)
\quad\Longleftrightarrow\quad
q_{r|a}(\phi_B)=c_{r|a}\ (\forall a,r)
\quad\Longleftrightarrow\quad
\Delta E_{\rm fine}(\rho)-\Delta E_{\rm direct}(\rho)=0\ (\forall\rho).
\tag{6}
$$

**证明。** 两个配置点x、y间，细分支核与粗分支核之比为$\sum_{r\in a}\sqrt{q_{r|a}(x)q_{r|a}(y)}$。它不超过1，等号恰要求两个归一正向量$\sqrt{q(x)}$、$\sqrt{q(y)}$相同。用跨两个规范轨道邻域的平滑不变波包检验核，得到对完整物理态域的必要性；常数条件反过来逐Kraus直接保证通道相同。能源条件由式(3)—(4)及$\mathcal K^{-1}>0$得到：若差对所有不变紧支撑源为零，各条件概率梯度必须逐点为零；连通性将它提升为常数。

若只限制到特殊码、非连通域或一种固定态，则全域必要性不能照搬；尤其单个态的能源差为零未必表示所有态上的仪器等价。允许保留完整$\mathcal J_a^{\rm fine}$时，粗报告本来就是合法过程，不要求换为$\mathcal J_a^{\rm direct}$。

## 5. 只用577原读口的具体反例

同一B上连续做两次原理想读口，中间等待取0。该选择只是有限仪器复合的明确测试，不宣称真实器件零时长。令$u=\sin s/2$，两个细结果r,t各为±1，报告a=rt：

$$
e_r=\frac{1+ru}{2},\qquad
L_{rt}=\sqrt{e_re_t},\qquad
p_a=\frac{1+a u^2}{2},\qquad M_a=\sqrt{p_a},\qquad
\mathcal J_a^{\rm fine}=\sum_{rt=a}L_{rt}(\cdot)L_{rt},\quad
\mathcal J_a^{\rm direct}=M_a(\cdot)M_a.
\tag{7}
$$

细读可由两份原末读复合实现；直接奇偶读是另给的比较仪器，未声称原H自动产生它。两者报告a的完整概率相同；奇分支实际上相同，偶分支的两个条件概率依赖s，所以不相同。

将导数相对于s的平方和记为A，直接求导或用式(4)给

$$
\Delta A(s):=\sum_{r,t}|L'_{rt}|^2-\sum_a|M'_a|^2
=\frac{(u')^2}{2(1+u^2)}
=\frac{\cos^2s}{8(1+\sin^2s/4)}.
\tag{8}
$$

因此，**完全相同的粗报告不能消除源中的额外反作用**：

$$
\delta E:=E_{\rm after,fine}-E_{\rm after,direct}
=\frac{\hbar^2}{2w_B}\langle \mathcal K^{ss}\Delta A\rangle,
\qquad
\mathcal K^{ss}=F\left(1-\frac{s^2}{6M}\right).
\tag{9}
$$

在$0<s<\pi/2$的非零紧支撑Gauss来源上严格为正。将细记录放在整体内部但限制当前访问，不会物理撤销这份反作用；只有另证相干恢复／反馈过程才可能改变它。

## 6. 同一原Hamiltonian、同一旧读口能够看到差别

不能仅以无界能源差代替实际读数见证。沿577取$b=\sin s_B$、$D=(i/\hbar)[H,\cdot]$，将两种非选择通道记为$\Phi_f,\Phi_c$。所有配置乘法量在两通道下保持；一阶微分项因$\sum LdL=0$也相同。

原完整H给$D^2b=C_B+f_B$，其中f_B是含完整势的乘法项，C_B仅对B微分。C_B的二阶系数为$w_B^{-2}\mathcal K^{ia}\mathcal K^{jb}\nabla_a\nabla_b b$；其余一阶和零阶项在通道差中抵消。于是得到完整算符核上的恒等式

$$
(\Phi_f^\dagger-\Phi_c^\dagger)b=0,\qquad
(\Phi_f^\dagger-\Phi_c^\dagger)Db=0,\qquad
(\Phi_f^\dagger-\Phi_c^\dagger)D^2b
=\frac{\hbar^2}{w_B^2}\Delta A(s)\,
\operatorname{Hess}_{\mathcal K}b(\nabla_{\mathcal K}s,\nabla_{\mathcal K}s).
\tag{10}
$$

这里保留了全部其它节点、测地边势、链路动能与磁项：它们在此特定通道差中或不出现，或以配置乘法出现而抵消，并非从H删除。二阶微分算子乘法展开给式(10)，不取经典或小ℏ极限。

记$a_s=1-s^2/(6M)$。由$\nabla s(\mathcal K^{ss})=-2F^2s a_s/(3M)$和Hessian链式法则，

$$
\operatorname{Hess}_{\mathcal K}(\sin s)(\nabla s,\nabla s)
=-F^2\left[a_s^2\sin s+\frac{s a_s\cos s}{3M}\right]<0
\qquad(0<s<\pi/2).
\tag{11}
$$

采用584相同紧支撑径向包，$h\in[.42,.62]$、$s\in[.25,.55]$，其它节点取同类包，规范链路取Haar常函数。该源正规、严格Gauss，并属于各固定次H幂域；允许附加不变相位$e^{i\eta hs}$，因为支持远离h=0。读取的光滑乘法保紧支撑及各幂域。

令$\rho_f=\Phi_f(\rho)$、$\rho_c=\Phi_c(\rho)$，后继都只按原H等待，再使用原$E_+=1/2+b/4$末读。存在充分小的正时间，使两概率不同，且

$$
P_f(+;t)-P_c(+;t)
=\frac{\hbar^2t^2}{8w_B^2}
\left\langle\Delta A\,
\operatorname{Hess}_{\mathcal K}b(\nabla s,\nabla s)\right\rangle_\rho
+o(t^2)<0.
\tag{12}
$$

有限时间存在性来自C²期望、相同的零一阶以及严格负的二阶系数，不靠外推有限网格时间模拟。未给出对尺度、源族或ℏ统一的可读时间窗。原奇偶报告也可保留；非选择未来概率已不同，故完整“奇偶—等待—再读”联合通道不可能相同。

## 7. 同一差异进入局部几何来源

原配置分布不变，故现场、测地边、磁势期望不变；只有B动能变化。在585的端点／面顶点正分配中，$\delta\langle h_i\rangle=\delta_{iB}\delta E$。对于该轮邻居平均Pα和相同几何校准，

$$
\delta\langle H_\alpha[N]\rangle=(P_\alpha N)_B\delta E,
\qquad
\left.\partial_{\psi_B}\delta E\right|_{\rho,\epsilon}
=-\frac6{\psi_B}\delta E.
\tag{13}
$$

固定态偏导使用与ψ无关的原配置测度字典；不能解释为改变几何时重新选择态后的全导数。前式也不能从单位lapse唯一决定α，585的非唯一性仍保留。得到的是同一物质操作对**已定义几何来源的响应**，不是完整协变应力或动态Einstein解。

## 8. 有限历史合法，但装置仍需内部实现

以原H等待与原二元读取串接，保留每个结果，任意参考恒等扩张。各分支为有界Kraus算符，完备性逐步给完整记录通道；Gauss逐分支保持。单次原读的$\sum_r|L'_r|^2=\cos^2s/[16(1-\sin^2s/4)]\le1/16$，且$\mathcal K^{ss}\le M$，所以对固定有限图、最小节点体积$w_{\min}>0$，

$$
K_{\boldsymbol r}=L_{r_k}^{(B_k)}U_{t_k}\cdots L_{r_1}^{(B_1)}U_{t_1},
\quad \sum_{\boldsymbol r}K_{\boldsymbol r}^\dagger K_{\boldsymbol r}=I,
\qquad
\mathbb E[E_{\rm source,after\ k}]
\le E_0+\frac{k\hbar^2M}{32w_{\min}}.
\tag{14}
$$

在结果控制下一次节点和等待时，按条件态再平均仍成立，前提是有限步数、相同H及同一体积下界。各结果概率至少1/4，因此每个有限历史的条件能量也有限，但不主张与历史长度无关的逐分支界。该量是源能源，不包含实现、保存、隔离记录及控制器的全部代价。

这一步使“单次允许末读”成为可组合的**指定有限操作接口**；没有证明这些理想仪器已由原H自行运行。不能把记录寄存器的数学存在当成新增自由度零成本，也不将式(14)的固定图界推广到$w_{\min}\to0$的无限细化。

## 9. 可复算结果与核验范围

五组检查均通过：原两读奇偶仪器；多变量Fisher分解和无损条件正对照；原紧支撑Gauss包的曲动能；完整五配置坐标的嵌套微分对易子；同源空间权及lapse响应。

取ℏ=.7、w_B=.8、包内不变相位η=.23。48、80、128阶独立求积趋于相同结果：

|量|两次细读、仅报告奇偶|直接平方根奇偶读|
|---|---:|---:|
|即时偶报告概率|0.51918059089555|0.51918059089555|
|源注能|0.06403886186800|0.00467852164738|
|读取后B动能|263.5666107496204|263.5072504093998|

差能源为0.05936034022062。独立展开各读后波函数梯度与式(9)的差约5.63×10⁻¹⁴；允许非零原概率流，不依赖实态交叉项恰为零。其它原H项保持由完整算符证明承担，表中未把单节点动能当全图总能量。

后继原＋记录概率差的t²系数为−0.01578008682333。五坐标局部嵌套对易子诊断在步长.024、.012、.006的误差约4.01×10⁻⁴、1.00×10⁻⁴、2.65×10⁻⁵；局部有限差分不是完整图传播模拟，有限时间结论由式(10)—(12)证明。

在ε=.8、ψ_B=1.13时，同源差能源为0.04454979875422，空间权偏导−0.23654760400469。背景权中心差分误差按二阶缩小；585局部lapse分配的对偶求和误差为零。没有图像检验。

## 10. 对联合目标的实际推进

[条件账](586/unified_physics_condition_ledger_586.md)将C01操作、C03记录、C20尺度压缩、C21物质读者和C22反作用接入同一原模型。这次排除的是：**只保持共同报告的效果概率，就用直接平方根粗读替代真实细过程，并继续沿用同一后态、未来记录和几何来源。**

可行替代很明确：保留真正的粗化CP映射及隐藏关联；或在当前仪器类内证明式(6)；或另行实现并记账相干恢复。没有证明所有粗粒化都失败，也不把这个有限接口当成统一理论完成。

下一步回到同一物质的局部能源流及切向来源：完整过程已明确必须保留哪些仪器反作用，不能再以摘要效果替代它。原586能源流[候选原稿](586/drafts/STATUS.md)保留，后续另起编号，先检验原H的流与几何生成元，不继续扫描奇偶读数精度或弱测参数。
