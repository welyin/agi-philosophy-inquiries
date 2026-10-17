# Witten世界线代数：原文身份及有限任务边界

2026-10-17。非编号文献核查，不新增模型、代码或科学轮次；不改变现有目标或参考。

## 1. 编号与查重

本次原文是 Edward Witten，[A Background-Independent Algebra in Quantum Gravity，2308.03663v3](https://arxiv.org/html/2308.03663v3)。[2309.15897](https://arxiv.org/abs/2309.15897)对应另一篇 Generalized Black Hole Entropy is von Neumann Entropy，不能互换引用。

已对研究目录检索编号、标题、worldline algebra、世界线代数及Hartle–Hawking。未找到2308.03663的既有专门范围审核；1004的Hartle–Hawking黑洞参考不是本次普适无边界态命题。此检索不意味着没有任何间接引用，更不意味着文献结果是本研究的新发现。

## 2. 原文定位摘要

Primary-source findings, limited to §§2–4:

> 本节为概述，不是原文引语。

Section 2 introduces a minimal geodesic observer with Hamiltonian $m+q$, $q\ge0$, a worldline action, and constrained field operators. Background independence concerns an OPE algebra; its Hilbert representation requires a spacetime and worldline. Remark 2 limits the clock interpretation. Remark 3 requires an additional frame for tensor observables and omits its compact rotation sector from the explicit treatment.

Sections 3.1–3.2 establish the tracial property of the chosen de Sitter reference. Its completion and entropy concern the stated perturbative representation. Section 3.3 permits an upper clock-energy bound; this is not a finite-time instrument construction.

Section 4.1 relates the semiclassical de Sitter discussion to generalized entropy and an observer contribution, retaining normalization and approximation qualifications. Section 4.2 explicitly presents the universal Hartle–Hawking state and universal maximum-entropy interpretation as hypotheses. Neither a universal theorem nor a completed finite apparatus follows from that proposal. [原文§§2—4](https://arxiv.org/html/2308.03663v3)

## 3. 最小作用确实给来源，但“完整”须注明对象

不能把本论文说成只有抽象代数、完全没有物质来源。对其最小观察者作用作通常的度规变分，在给定光滑背景和世界线的范围得到分布来源

$$
T_{\rm obs}^{\mu\nu}(x)
=\int d\tau\,(m+q)\,U^\mu U^\nu
\frac{\delta^{(4)}(x-\gamma(\tau))}{\sqrt{-g(x)}}.
\tag{1}
$$

这是对原文作用的直接变分核对，不是新增材料。若 $W_{\mu\nu}$ 是紧时间支持的平滑测试，世界线段上有界，且原准备有有限 $\langle q\rangle$，则

$$
\left|\langle T_{\rm obs}[W]\rangle\right|
\le (m+\langle q\rangle)
\int_{\gamma\cap\operatorname{supp}W}
|W_{\mu\nu}U^\mu U^\nu|\,d\tau<\infty.
\tag{2}
$$

因此，最小观察者的有限涂抹源不存在“因采用世界线代数就必然无限”的障碍。式(2)只核此给定背景中的原点粒子来源，既不求解其自场，也不证明自引力处方、读出控制、支持结构和新记录的完整来源。有限一次矩不是后者的替代证书。

对当前[有界参考查询](../bounded_reference_query_after1270/observer_finite_query_gate_after1270.md)，$q$的两份准备及有界效果本身已有正常性和有限矩。这项旧代数结果与(2)相容；但它没有把[读出映射审计](../reference_query_qualification_after1270/reference_query_readout_mapping_after1270.md)缺少的控制与约束来源身份自动补上。不能一边沿用无新增相互作用的最小作用，一边把实际查询的驱动、旗标及后态认作已由该作用实现。

## 4. 与2026全荷静态patch构造不是无损替换关系

已有[2601.07910范围审核](../quantum_identity_gate_after1270/realistic_observer_2026_screen_after1270.md)核定：它保留同微扰阶的角向边界荷，观察者的旋转与形变能源进入完整约束、参考和熵。相关原文入口为[§§2.6—2.7、5.1—5.2](https://arxiv.org/html/2601.07910v1)。这与本次最小作用的用途不同。

本项目的独立匹配判断是：

- **不能按名称直接替换。** 世界线、正能钟、有限迹这些共同名称，不给两个物理代数、完整来源和同一参考的等同映射。
- **限制查询不等于删除来源。** 只查询 $e^{-\beta q}$ 可以限制菜单；若仍要求是原全荷观察者，就须保其辅助荷准备及对总能源、面积荷的贡献。菜单看不见某部门不证明该部门源为零。
- **代数生成与实际访问不同。** 即使一个目标量可由世界线代数表达，也尚未给该表达在当前材料中的有限时间、精度、后态和控制来源。已有1270有限边界原语可作模板，但其原SC码、FLRW来源和钟不能直接迁移。
- **不预设必须全模式制造。** 项目可以采用明确有限有效输入；这里仅判定所查文献没有给出保留本次全部身份的现成替换。若改用受限最小观察者，将是另一份公开限定的候选，不能声称仍保原全荷参考与几何熵。

这里没有证明两个构造不能存在有用映射。当前缺的是映射本身：同一准备和有界效果的限制/运输，完整能源与面积荷的对应，以及新增读口的来源和误差。不能用“普适最大熵”猜想跳过它，也不要求先实现完整量子引力或所有微观装置。

## 5. 判定与停止

**PASS_AS_RESTRICTED_ALGEBRAIC_REFERENCE；NOT_A_COMPLETE_FINITE_TASK_CERTIFICATE。**

可以保留的是：原理层面的世界线关系描述、给定最小作用的来源，以及特定静态patch参考下的代数与熵接口。它不能直接替换2026无限辅助部门同时保当前全荷来源和参考，更未补出本项目缺少的实际查询过程。

因此，本次不推荐因“世界线代数更简单”重开一个观察者制造支线。已有有界查询和读出映射的停止线保持；没有新实验、参数扫描或编号。Hartle–Hawking普适最大熵继续只记为文献猜想。
