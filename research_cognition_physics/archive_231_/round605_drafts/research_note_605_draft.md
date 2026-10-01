# 第605轮：手征配置限制、Gauss热态与共同来源域

日期：2026-10-01。接[604](research_note_604.md)，复用[569](research_note_569.md)全局商群、[598](research_note_598.md)有限CAR和[603](research_note_603.md)完整热迹。[代码](joint_admissible_gauss_domain.py)、[结果](joint_admissible_gauss_domain_results.json)、[核验](research_round_605_checks.json)、[条件账](unified_physics_condition_ledger_605.md)。四组复算、十四式；主代理审查，无新增独立代理审查。

## 1. 明确接口，而非以一个理论名称代替连接

604排除中心差分的一条物种接法。本轮先核成熟手征路线所需的配置条件，是否与原Gauss热态及能源来源共存；不重复节点复制，也不开始主体装置设计。

[Hernández、Jansen与Lüscher](https://arxiv.org/abs/hep-lat/9808010)证明Neuberger算符在充分光滑规范背景上的指数局部性。[Lüscher的具体表述](https://arxiv.org/pdf/hep-lat/9909150)第2节给出一个充分条件：所有面满足表示中的算符范数界，ε<1/30；随后还需处理随规范场变化的手征投影与测度。这是特定Euclidean格点框架的充分条件，不是所有手征理论的必要条件。

**本轮新增连接：原完整有限图的忠实Gauss热态不能精确满足排除非空开放配置区的硬限制；可以改变闭形式定义域，使空间限制与Gauss、热迹及原读取来源同时成立，但这会改变动力学和参考态。** 直接筛选旧热态一般不等价，甚至可能注入无限动能。把同样截断直接施加给时间面权重，又须另查转移矩阵正性。

|层次|地位|
|---|---|
|认知动机|共同物质、态、允许操作与来源必须在同一个实际配置域中成立|
|继承输入|固定有限图、原曲目标及非线性束缚势、商群、正电动能、有限CAR、正给定几何|
|新增分支输入|一个规范不变的空间允许域Ω、硬边界或有限惩罚参数λ；不把ε当成认知公理|
|解析结果|全支撑障碍、受限闭形式的Gauss与热迹、软惩罚极限和概率界、原仪器来源保留|
|数值诊断|原表示与磁势；一维规范转子中的惩罚和投影；时间单面权重的负Fourier模|
|未完成|四维Euclidean手征测度、完整Hamiltonian对应、目标手征谱、Lorentz及GR连续匹配|

原[605候选](round605_drafts/STATUS.md)保留，以下是其配置域接口的实际结果。

## 2. 同一商群上的允许配置

令Q_g=G^E，G为569的Z6商群，R取598全部单粒子表示的直和；自旋重复不改变范数。图包含至少一个普通非退化面。定义

$$
r_p(U)=\|I-R(U_p)\|,\qquad
\Omega_\epsilon=\{U\in Q_g:\max_p r_p(U)<\epsilon\},\qquad
E_\Omega=1_{\Omega_\epsilon},\quad E_{\rm bad}=I-E_\Omega.
\tag{1}
$$

R下降到商群，故该定义不依赖链路lift；面回路在基点共轭下范数不变，EΩ与Gauss投影交换。Ω是含平坦配置邻域的开放集。对ε=.02，取单面纯超荷回路z=exp(iπ/12)，R中电荷−6块给范数√2；因此坏区也含非空开放集。原569磁势在该点有限，数值约7.85481308。

限制作用于空间链路。不能据此宣称已经满足四维Euclidean所有时空面上的局部性定理；本轮还没有四维Dirac算符。

## 3. 原完整热态不能把坏区精确删除

记H为598—603原完整物理Hamiltonian，β>0。603保证Z_phys有限，e^(−βH)在Gauss空间上无核。E_bad在该空间非零：可取坏区内的光滑规范不变面函数，乘规范不变的标量紧支撑径向函数与费米空态。这个向量存在不要求费米空态本身是H的不变部门。

因此

$$
p_{\rm bad}(\beta)
=\operatorname{Tr}_{\rm phys}(\rho_\beta E_{\rm bad})
=\frac{\|E_{\rm bad}e^{-\beta H/2}\|_{\rm HS}^2}{Z_{\rm phys}}
>0.
\tag{2}
$$

严格性来自热算符稠密值域与非零投影，而不是经典Boltzmann权或无符号路径积分；允许原矩阵质量和配对。

即使加有限惩罚，也有

$$
H_\lambda=H+\lambda E_{\rm bad},\qquad 0\le\lambda<\infty,\qquad
D(q_\lambda)=D(q),\qquad p_{{\rm bad},\lambda}>0.
\tag{3}
$$

该项有界、非负、规范不变，故不破坏原半有界自伴性及热迹存在；但不产生精确零概率。它的能源与全部H一起进入全局lapse来源，不能视为免费删除。此结论只否定“原模型已自动满足硬允许域”，不证明所有不满足充分光滑界的配置都导致非局部算符。

## 4. 改变定义域可以保留正量子过程和热态

为精确施加空间限制，定义新的物理Hilbert空间及闭形式：

$$
\mathcal H_\Omega=E_\Omega\mathcal H_{\rm phys},\qquad
D(q_\Omega)=D(q)\cap\mathcal H_\Omega,\qquad
q_\Omega=q|_{D(q_\Omega)}.
\tag{4}
$$

这里使用零延拓定义域：波函数在坏区为零，并仍具有原有限形式能源。它在形式范数下闭合，在HΩ中稠密；后者由Ω内部紧支撑光滑函数及紧群平均得到。由闭形式表示定理得到HΩ自伴且与原H同下界。光滑／Lipschitz边界时就是通常Dirichlet条件；一般边界以(4)定义，不擅自认定任意Sobolev空间等同。

Gauss保持、原势对CAR的相对形式界均在该子域继续成立。将HΩ本征值按重数排列，min–max所用测试子空间变少，因此

$$
E_k(H_\Omega)\ge E_k(H),\qquad
0<Z_\Omega(\beta)=\operatorname{Tr}_{\mathcal H_\Omega}e^{-\beta H_\Omega}
\le Z_{\rm phys}(\beta)<\infty.
\tag{5}
$$

紧预解来自原形式域紧嵌入的限制。Ω内仍有无限多个标量径向函数，物理空间不是零空间。全部能量矩随之有限。**这是一份明确改变了模型的可用分支，不是原热态的同义表示。**

若Ω及ε不依赖给定几何参数θ，则形式域对θ固定，原来源G=∂θq限制后仍满足603的界：

$$
G_\Omega=A_\Omega^{1/2}\mathsf B_\Omega A_\Omega^{1/2},\quad
\|\mathsf B_\Omega\|\le L,\quad A_\Omega=H_\Omega+c\ge1,\qquad
0\le\chi^R_\Omega(0^+)\le\chi^E_\Omega
\le L^2\bigl[4\langle A_\Omega\rangle+5\beta\langle A_\Omega^2\rangle\bigr].
\tag{6}
$$

共同界从原形式范数限制得到；均值和谱必须使用新热态。若让ε依赖θ，边界随参数移动，本证明不覆盖形状导数和边界力；不能继续套用旧G而遗漏新增来源。

原L±(s_B)只作用标量，保持支持且保形式域，故

$$
\mathcal I^*(H_\Omega)-H_\Omega=D_B|_{\mathcal H_\Omega},\qquad
\mathcal I^*(G_\Omega)-G_\Omega=(\partial_\theta D_B)|_{\mathcal H_\Omega}.
\tag{7}
$$

同理有限λ惩罚与L±对易，也保持原注能身份。这把限制、Gauss、读取和来源接在一起；不保证全部旧仪器或任意锐投影都保新域。

## 5. 有限惩罚到硬域的受控连接

qλ随λ单调增加。有限极限域恰为D(q)∩ker E_bad，即(4)。因此在原Hilbert空间上，

$$
e^{-\beta H_\lambda}\longrightarrow
E_\Omega e^{-\beta H_\Omega}E_\Omega,\qquad
Z_\lambda(\beta)\downarrow Z_\Omega(\beta),\qquad
\rho_\lambda\longrightarrow\rho_\Omega
\quad\text{以迹范数}.
\tag{8}
$$

说明所用极限：先移位使H≥1。闭形式单调极限给广义强预解及热半群强收敛；极限在HΩ正交补为零。对每个k，E_k(Hλ)单调且不超过E_k(HΩ)。原紧形式嵌入使相应低能向量有强收敛子列，而λ‖E_badψ‖²有界强迫极限进入Ω；min–max遂给本征值极限等于E_k(HΩ)。由原可和热权支配逐项求和，得到迹收敛；正热算符的强收敛加迹收敛给迹范数收敛。没有错误使用“算符指数对任意算符序单调”。

设Fλ=−β⁻¹log Zλ。Gibbs变分原理对ρλ与新ρΩ分别比较，给

$$
0<p_{{\rm bad},\lambda}
=\partial_\lambda F_\lambda
\le \frac{F_\Omega-F_0}{\lambda}\quad(\lambda>0),\qquad
F_0\le F_\lambda\le F_\Omega.
\tag{9}
$$

证明不需要H与E_bad对易：Fλ=F₀[ρλ]+λp_bad≥F₀+λp_bad，且ρΩ为λ族合法试探态并使惩罚均值为零。常数依赖原模型、β、Ω；不是全图细化一致误差。该界控制坏区概率，不单凭它推出所有无界来源的误差。

## 6. 筛选旧热态并不实现新动力学

一次硬筛选旧热态得到

$$
\rho_{\rm sel}
=\frac{E_\Omega e^{-\beta H}E_\Omega}
{\operatorname{Tr}(E_\Omega e^{-\beta H})}
\quad\not\equiv\quad
\frac{e^{-\beta H_\Omega}}{Z_\Omega}.
\tag{10}
$$

左边原演化允许离开Ω再返回；右边是整段演化都满足新边界。它们一般不同，且锐乘法不一定保原动能形式域。

取明确诊断：一维U(1)物理回路转子，H=−∂θ²，θ∈S¹，Ω=(−θ₀,θ₀)。这是域机制的可解例子，不是原全非线性模型的配分函数。转子Gibbs对角位置密度均匀，所以

$$
p_\Omega=\theta_0/\pi,\qquad
\rho_{\rm sel}=\frac1{p_\Omega}\sum_{n\in\mathbb Z}
\frac{e^{-\beta n^2}}{Z}\,
|E_\Omega n\rangle\langle E_\Omega n|.
\tag{11}
$$

每个被截断平面波在端点跳跃，不属于圆上的H¹形式域；n=0的权已严格正。正性保证各项不能相消，故

$$
\operatorname{Tr}(H\rho_{\rm sel})=\infty,\qquad
\operatorname{Tr}(H_\Omega\rho_\Omega)<\infty.
\tag{12}
$$

前者是解析域结论，不由有限网格发散趋势代替证明。θ₀=π/3、β=.7的网格核验中，锐筛选动能随N=96/192/384为15.75/30.35/59.53，受限热态动能约2.3085/2.3092/2.3094。

## 7. 不能把空间限制直接搬成所有时间面权重

[Creutz，Positivity and topology in lattice gauge theory](https://arxiv.org/pdf/hep-lat/0409017)明确假设单面实非负权重、单位元附近解析、其他处逐段光滑；在该类中，禁止一块开放面配置区与正转移矩阵冲突。不是针对所有Hamiltonian或所有手征正则化的普遍不可能性定理。

其U(1)必要条件可直接复算。时间单面核T(θ′,θ)=W(θ−θ′)若为正算符，Fourier本征值须满足

$$
\widehat W_n=\frac1{2\pi}\int_{-\pi}^{\pi}
W(\phi)e^{in\phi}\,d\phi\ge0
\quad\text{对全部整数 }n.
\tag{13}
$$

具体截断Wilson型权重给

$$
W(\phi)=e^{\kappa(\cos\phi-1)}1_{|\phi|<\theta_0},\qquad
\kappa=1.5,\quad\theta_0=.4,\qquad
\widehat W_9\simeq-0.0121367421<0.
\tag{14}
$$

原群包含这样的U(1)方向，故不能只凭W逐点非负便签收转移矩阵正性。此例旨在排除“把硬截断随手加到时间面上”这条修补接法，不是在重新发现Creutz定理。

第4节的受限空间Hamiltonian仍有正热算符e^(−βHΩ)，并不矛盾：它没有要求所有时空面共用上述局部单面核，也没有因此建立Euclidean GW测度或旋转／Lorentz一致性。恢复这几个性质仍是联合条件。

## 8. 复算与条件压缩

|组|内容|结果|
|---|---|---|
|1|598原表示的允许域在商lift及规范共轭下不变；569坏点磁势|12组最大误差1.11×10⁻¹⁵；坏点势有限|
|2|含.4(1−cosθ)势的N=120转子，β=.8，λ=0/5/50/500/5000|坏区概率始终正，下降至2.86×10⁻⁶；满足(9)；Fλ趋近FΩ≈2.27107|
|3|锐筛选与新定义域的差别|N=120态迹距≈.378219；细化锐筛选能源发散、新热能源收敛|
|4|时间面权重Fourier模|首个明显负模n=9；256与512点求积差4.17×10⁻¹⁶|

三、四组分别是域和正性诊断；完整原模型的存在与支撑结论由第2—5节证明。没有用转子成功代替标准模型完成。

本轮条件合并得到三条分支：保留原全配置域，另找不要求硬限制的局部性控制；改变空间闭形式域并保留既有Gauss／热来源结论；或接受某种Euclidean截止处正性限制，再单独证明物理极限。后两项都不能隐瞒其新增输入。

下一步[606](round606_drafts/STATUS.md)核随规范场变化的手征子空间与原电动能是否相容：把固定CAR纤维替换为变化投影时，联络和几何补项是否被迫进入同一能源／来源。先复用已有几何和反作用结果，再用明确投影族检验；不把一般投影恒等式单独冒充物质统一证明。整体目标保持，认知实现后置。
