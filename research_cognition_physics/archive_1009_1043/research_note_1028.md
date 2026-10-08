# 第1028轮：完整物质源到单一活跃引力理想的条件桥梁

2026-10-08。接[1027下一项](1027/NEXT.md)。本轮将1027的完整带荷源与1018二阶相容性接到同一平直父作用，新增一组合法分支及规范不变量校准，累计3805；新认知公理0。不是新发现1018的耦合代数式，也不是全部认知输入已被消除。

## 1. 问题与结论

[1027](research_note_1027.md)已在明确源类中锁定P981全部带荷部门的相对来源，但[1018](research_note_1018.md)的二阶证明使用正则标量物质。[1024](research_note_1024.md)正确保留了这条完整物质桥梁，不能用任意Higgs数据令规范场为零来代替证明。

本轮给一条较小的真实连接：完整P981物质作用存在规范流恒零的**固定内部方向、实径向Higgs局部解**。对其规范不变量 $O=\Phi^\dagger\Phi$，标准内部规范补偿恒为零；1018的首阶反作用及Killing参数论证可以实际接入。

在第2节完整一阶最小数据下，任意正规局域完成必须满足

$$D_{ab}:=q_aq_b-\delta_{ab}\kappa_aq_a=0.\tag{1}$$

因此其最小一阶动态源权只能全为零，或只在一个非自由引力理想上非零，并有 $q_r=\kappa_r$。允许分支存在一个经典最小完成。**完整一阶代表、自由自旋2、共同Lorentz主部仍是条件**；本轮没有证明所有第一阶耦合都等价于此代表，也不把一阶零权误称为所有可能高阶交互均为零。

## 2. 同一扩大父族与完整一阶输入

物质部分是[P981](../archive_956_989/981/drafts/common_parent_contract_v1.md)的平直领先规范—Higgs—Weyl作用，保无导数Yukawa、所声明Weinberg项及实际Higgs势。可纳入[B993](../archive_990_1008/993/common_candidate_v1.md)中性标量；为把它也纳入同一 $T_{\rm full}$，使用1027所需的非零门户。零门户分离边界另列。

规范、物质相互作用不随形式参数 $\varepsilon$ 关闭。引入同背景、正号Pauli–Fierz字段 $h^a$，纯引力理想的三次系数为 $A^c{}_{ab}=\delta_{ab}\delta_{ac}\kappa_a$，沿用[358](../archive_342_369/research_note_358.md)与1018的归一。写

$$S=S_0+\varepsilon S_1+\varepsilon^2S_2+\cdots,\qquad S_0=\sum_a S_{\rm PF}[h^a]+S_m.\tag{2}$$

采用完整最小一阶Noether代表：

- $S_1$ 为各理想的Einstein三次项及 $\tfrac12\sum_aq_a\int h^a_{\mu\nu}T_{\rm full}^{\mu\nu}$。
- 纯自旋2零阶作用 $\delta_0h^a=2\partial_{(\mu}\xi^a_{\nu)}$，物质上的纯自旋2零阶作用为零。
- $\delta_1h^a=\kappa_a\mathcal L_{\xi^a}h^a$；物质采用 $\sum_aq_a\mathcal L_{\xi^a}$ 的标准最小张量／旋量提升及所需内部补偿，因此 $\delta_1O=\sum_aq_a\xi^{a\mu}\partial_\mu O$。
- 两个纯自旋2参数的一阶括号，其**完整引力分量**为 $B_1^a=\kappa_a[\xi^a,\zeta^a]$，没有另留作用于 $O$ 的未知一阶参数分量。
- 内部Yang–Mills与局部Lorentz行动保持标准，所保留各阶均恒湮灭 $O$；允许场依赖的内部补偿。普通Lie与规范协变Lie的内部差别也恒湮灭 $O$。

在metric-PF及对称vierbein写法中，局部Lorentz项是旋量提升所需的局部标架补偿；未定规vierbein写法中则为既有标架冗余，不是新添一套独立传播字段或自由规范内容。

这比只规定 $hT$ 顶点更强。可用辅助 $g_{\rm aux}=\eta+\varepsilon\sum_aq_ah^a$ 的最小物质协变化**仅展开到一阶**来定义该代表；任意 $q_a$ 此时均允许，尚未选择唯一活跃理想。不能把此辅助写法无审计延成任意 $q$ 的完整多度规理论。

本轮取1027源类的最小代表 $I_a=0$、额外独立常量 $C_a=0$；原物质势内的真空常数保持。没有证明全部非最小改进和独立tadpole完成均属此类。允许任意未知 $S_2,\delta_2,B_2$ 等、允许在壳闭合，但须在 $\varepsilon=0$ 正规、局域、每阶有限导数且不增加自由规范内容。采用自由完整代数离壳闭合代表 $M_0=0$，不假定 $M_1=0$。

## 3. 合法径向分支为何属于完整原作用

取一个局部空间坐标 $y$，令

$$A=0,\qquad\psi=0,\qquad\sigma=0,\qquad \Phi=\frac1{\sqrt2}\begin{pmatrix}0\\f(y)\end{pmatrix},\qquad f\in\mathbb R.\tag{3}$$

这里零规范场只用于证明存在一个特殊完整解，并非任意Higgs解的截断规则。对每个Hermitian内部生成元 $T_A$，Higgs电流正比

$$i\bigl(\Phi^\dagger T_A\partial^\mu\Phi-(\partial^\mu\Phi)^\dagger T_A\Phi\bigr)=0,\tag{4}$$

因为实 $f$ 与导数沿同一固定内部向量；色作用则本来平凡。因此全部规范方程在 $F=0$ 下成立，Gauss也成立。

Higgs规范不变势在该分支成为 $U_{\rm rad}(f)$。完整Higgs方程沿径向化为

$$f''(y)=U_{\rm rad}'(f(y));\tag{5}$$

其它实Higgs方向的方程为零。标准 $U_{\rm rad}=\lambda_H(f^2-v_H^2)^2/4+V_0$ 给 $f''=\lambda_Hf(f^2-v_H^2)$。Yukawa和Weinberg项的物质变分在全部费米字段为零时成立，且没有额外Higgs源。B993势及门户关于 $\sigma$ 为偶，其一阶变分在 $\sigma=0$ 为零。

光滑径向势使任意有限初值 $f(0),f'(0)$ 有局部真解；可选 $f(0)f'(0)\ne0$。这只是本类别内的合法数学解，不声称是自然宇宙选中的状态、稳定真空或有限能全球解。

于是 $O=f^2/2$ 且 $\partial_yO=ff'\ne0$。这里O是用于经典Noether与形式形变测试的局部**内部**规范不变量，不是完整引力的微分同胚不变Dirac可观测量；它恰要按Lie作用变化。没有构造重整化复合算符或量子测量。证明不要求任意字段分量都独立活跃。若改用旋转相位的Higgs且仍设 $A=0$，一般会有非零超荷流；这正是1024原反例所禁止的更广截断。

## 4. 来源与首阶反作用保持在同一对象内

式(3)—(5)满足**全部**平直物质方程，其完整应力守恒。仅有物质非零而令全部引力场始终为零，一般不满足耦合方程，不能如此消去开放代数项。

直接复用[1018§4](research_note_1018.md)：对各守恒源 $q_aT[\Phi_0]$，局部Pauli–Fierz方程存在满足harmonic约束的响应 $h_1^a$。旧引理用初片Poisson资料使 $H_\nu=\dot H_\nu=0$，再由守恒保证约束传播；它不依赖物质是否只是单标量，只依赖当前源确实守恒。

取 $h^a=\varepsilon h_1^a$，全部物质保持上节解，则物质一阶耦合 $\varepsilon hT=O(\varepsilon^2)$，引力首阶源由 $h_1$ 抵消。因此全部Euler残差及任意固定阶局部导数满足

$$E=O(\varepsilon^2).\tag{6}$$

未知 $S_2$ 本身也从该阶开始。对正规开放代数算符 $M=\varepsilon M_1+\cdots$，有 $ME=O(\varepsilon^3)$。这消去了二阶测试中的方程平凡项，而没有另设 $M_1=0$。有限Taylor校准只核局部系数，实际光滑解与反作用的阶数由解析引理承担。

## 5. 规范不变量上的二阶障碍

在测试邻域选择平直Killing参数，所以纯自旋2 $\delta_0h$ 及所有导数为零，$\delta_0$ 对物质也为零。允许在邻域外乘紧支函数，局域有限导数身份不变。

故 $\delta_0\delta_2O$ 项不进入测试；$B_2$ 的自旋2部分通过 $\delta_0O=0$ 作用，也消失。内部Yang–Mills／局部Lorentz参数，无论是否依赖字段，均对 $O$ 恒为零。**须保第2节关于内部行动各阶的限定**，否则不能只凭零阶规范不变性忽略 $B_{1,\rm int}\delta_{1,\rm int}O$。式(6)另消去开放项。

沿1018的交换子顺序约定，唯一剩余条件为

$$\sum_{a,b}D_{ab}\mathcal L_{[\xi^a,\zeta^b]}O_0=0.\tag{7}$$

令仅 $a,b$ 槽活跃，选

$$\xi=\partial_x,\qquad \zeta=x\partial_y-y\partial_x,\qquad[\xi,\zeta]=\partial_y.\tag{8}$$

二者为Killing，且 $\mathcal L_{[\xi,\zeta]}O_0=ff'\ne0$。因此每个 $D_{ab}=0$。所测量的代数身份必须对每个合法局部解成立，一个真实父解已经足以提供必要条件；没有把整个父作用换成独立标量模型。

## 6. 条件分类、存在完成与未被排除的分支

对 $a\ne b$，式(1)为 $q_aq_b=0$；对角为 $q_a(q_a-\kappa_a)=0$。故

$$q=0\quad\text{或}\quad q_r=\kappa_r\ne0,\quad q_{a\ne r}=0.\tag{9}$$

给定允许非零分支，可把全部原规范—Higgs—Weyl物质及所声明Weinberg项最小协变化到 $g_r=\eta+\varepsilon\kappa_rh^r$，采用该理想的Einstein完成；其余引力理想保独立Einstein或自由Pauli–Fierz作用。这给一个正规**经典**完成的存在例，包含完整原物质作用而非只含径向截面。$q=0$ 则允许物质脱离这些引力部门。

这个存在例不选择绝对耦合、唯一完成、量子反常处方或真实宇宙态，也不以预置Einstein项证明认知生成GR。中性门户为零时，1027原本允许独立来源；不同已分离物质块可分别接不同理想，不被这里排除。完全脱离的引力字段、massive spin2、额外第一阶结构和改变自由规范内容的理论也不在结论中。

## 7. 本轮校准与一般证明的分工

[科学代码](1028/full_parent_gravity_bridge.py)与[结果](1028/full_parent_gravity_bridge_results.json)针对本轮新增映射核验：径向ODE的有限Taylor系数、完整四实Higgs方程、所有Higgs规范流与Gauss、规范不变量的逐次Killing变分，以及旋转相位但删规范场的失败对照。分支上的零费米源、偶中性场和门户边界另按完整菜单检查。

校准参数为 $\lambda_H=2/5$、$v_H^2=9/4$、$f(0)=5/4$、$f'(0)=2/7$。12阶Taylor给径向方程0—10阶、应力散度0—10阶以及四项电弱流0—11阶精确零；色流由表示平凡而为零。零费米分支由每个保留方程项至少一次含费米变量的解析次数判断承担，不冒称运行了完整Grassmann场求解器。

独立变化量 $\partial_yO(0)=5/14$，双重态复合量的四个乘积变分项给同一个Killing括号值。双活跃坏权样例的闭合残差为 $\tfrac5{14}\begin{pmatrix}-1&1\\1&-2\end{pmatrix}$；允许样例为零。相位变化率 $3/5$ 而错误删除规范连接时，除去各自耦合的超荷流、弱三流分别为 $-15/32$、$15/32$；补回对应纯规范连接后归零。这个对照检验截断合法性，不给出一个新的相位演化模型。

精确有理数承担有限代数与Taylor残差；局部ODE真解、全部有限导数计阶和Pauli–Fierz约束传播由第3—4节解析证明承担。没有把Taylor多项式当全局解，也不重跑1018大量耦合扫描或重新制造已知障碍矩阵。

默认复算不改保存结果：

```powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'research_cognition_physics/archive_1009_/1028/full_parent_gravity_bridge.py'
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'research_cognition_physics/archive_1009_/1028/verify_round1028.py'
```

[独立审阅](1028/review.md)专门核内部规范括号、首阶反作用和原父分支；[收据](1028/research_round_1028_checks.json)冻结本轮及明示旧输入，不冻结活动导航。

## 8. 与成熟研究、旧轮次的关系

[Boulanger等§9](https://arxiv.org/html/hep-th/0007220#S9)的显式物质论证针对标量，原文也没有声称覆盖所有更复杂物质。1018已核其相关限制并直接证明指定完整一阶数据下的多标量条件。本轮不把该文当成所有SM耦合已经分类的证据。

[326](../archive_301_341/research_note_326.md)为软散射接口，732／735为已给共同几何下的Ward与限定量子处方，均不是当前扩大父族的二阶接线。[去重审计](1028/selection_audit.md)列明增量：1027先锁完整来源，再以一个合法完整父解和规范不变量运输1018必要条件。

首阶反作用引理、Killing障碍形式和式(9)支持分类全部复用旧结果，不称新定理发现。新增的是其真实适用对象与证书；这补了一个此前明确留空的连接。

## 9. 依赖账的实际变化

1024的 scalar_source_to_parent 接口在1027指定领先源类中获得分类；gravity_ideal_to_parent 现在在**最小完整一阶代表**中有明确运输，强度为条件性共同对象映射。二者均不等于完整父量子过程已经实现，也没有关闭认知选择来源。

认知动机仍是同一过程的源、变换和反作用相容。独立物理输入仍包括共同Lorentz主部、字段及相互作用、健康自由自旋2理想、局域正规形变和完整最小代表。1027删除相对部门权，本轮进一步把其整体物质耦合与活跃理想的三次耦合相联；绝对数值及零耦合分支仍留存。

[依赖更新v0.17](1028/input_dependency_update_v0_17.md)明确这个变化。新增校准1组，累计3805，新认知公理0；有限有效域保持，没有要求先构造任意高能或无限精细完成。

## 10. 下一步与停止线

停止该桥的Killing样本、额外径向解和多引力候选修理。下一步先审“完整最小一阶代表”到底哪些部分可由更宽作用类别推出，哪些仍属独立采用；必须有真实的新选择命题才能开轮。若只是把现有Einstein最小协变化再展开，直接复用，不另计研究。

此外，非最小曲率系数、真空项和量子共同过程保留原登记；[983](../archive_956_989/research_note_983.md)早已说明共形值属于选择，不能重复把它当新缺口。详见[后续准入](1028/NEXT.md)。整体生成目标继续，当前桥完成不等于阶段统一完成。
