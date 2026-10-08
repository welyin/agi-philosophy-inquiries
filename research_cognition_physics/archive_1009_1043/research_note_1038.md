# 第1038轮：运动参考之间的方向接口需要哪些联合资料

2026-10-08。基线1035；与1036、1037独立。新增完整校准1组，新采用认知公理0，累计由主代理统一结算。

## 1. 本轮结论

**同一个方向qubit不能在任意运动参考之间独立运输，即使另外保存完整动量边缘，也可能丢掉决定读数的关联。** 在明确采用的有质量、自旋1/2、正能自由单粒子表示中，本轮给正常有限能波包：两个输入的spin边缘和完整动量边缘完全相同，换到同一惯性参考后，指定方向效果的概率差仍至少

$$\boxed{783/1700>0.46.}$$

因此，仅据这两个边缘预测，至少一份输入的误差不小于783/3400。保留动量分辨的条件spin资料可准确比较此菜单；有限分箱另有保任意旧参考的误差界。完整联合过程没有矛盾。

本轮还分类所有boost方向效果生成的数学von Neumann代数。**忠实保该代数是比逐项单次方向均值更强的合同**，不是已经从认知原则推出的权限。

## 2. 去重与输入

已读总导航、[1035](research_note_1035.md)、其NEXT/准入、[1009输入账](1009/input_dependency_ledger_v0_1.md)及[957假说](../archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md)。完整去重见[选题审计](1038/selection_audit.md)。

|旧结果|直接复用及边界|
|---|---|
|[382](../archive_370_428/research_note_382.md)—[384](../archive_370_428/research_note_384.md)|条件性三维针对真实方向壳与同一qubit效果；不重加384已消去的下界Lipschitz条件|
|[386](../archive_370_428/research_note_386.md)或[425](../archive_370_428/research_note_425.md)|真实邻域/光滑坐标的两条替代桥；不将这里的动量坐标冒称涌现位置|
|[522](../archive_467_530/research_note_522.md)—[523](../archive_467_530/research_note_523.md)|热参考与方向仪器已有共同条件实现；没有认证全部运动主体的独立qubit摘要|
|[720](../archive_702_741/research_note_720.md)|固定自由H的FW自旋与读口/来源运输；不同于本轮惯性参考比较|
|[923](../archive_923_934/research_note_923.md)|方向报告不等于全部未来资料；本轮给另一个明确物理表示、正常波包及任务代数分类|
|369、305、312、317、639、837、1020、963|模结构、热参考与钟标定的条件工具不重做|

**认知动机：** A4/A5/A7要求保同一任务的关联与比较字典。**物理采用：** 3+1 Lorentz结构、m>0、自旋1/2自由正能表示、canonical spin及声明的效果菜单。**未采用：** 把被动boost当作可任意执行的加速动作，或把canonical效果直接叫作现实Stern–Gerlach仪器。

## 3. 给定表示后，比较律由同一对象固定

取c=ħ=1、E(p)=sqrt(m²+|p|²)、dμ=d³p/(2E)，

$$
\mathcal H=L^2(\mathbb R^3,d\mu)\otimes\mathbb C^2,
\quad B(p)=\frac{(E+m)I+p\cdot\sigma}{\sqrt{2m(E+m)}},
\quad W(A,p)=B(\Lambda p)^{-1}AB(p).
\tag{1}
$$

U_Aψ(Λp)=W(A,p)ψ(p)是酉变换，且

$$W(A_2A_1,p)=W(A_2,\Lambda_1p)W(A_1,p).\tag{2}$$

无需为方向另拟独立旋转规则。纯旋转的W与p无关；一般boost随p变化。[Polyzou—Glöckle—Witała，§4、5.1](https://arxiv.org/html/1208.5840)提供成熟canonical-spin表示；[Peres—Scudo—Terno，式(9)及后文](https://arxiv.org/html/quant-ph/0203033)已有约化spin非自主协变结论，不作为本项目发现。

对单位轴n、rapidity ξ，设t=tanh(ξ/2)、v=p/(E+m)，直接相乘给

$$
W(p)=\frac{(1+t n\cdot v)I+i t(n\times v)\cdot\sigma}
{\sqrt{1+2t n\cdot v+t^2|v|^2}}.
\tag{3}
$$

详证、作用域及完整计算在[解析稿](1038/drafts/moving_direction_transport_derivation.md)。代码由独立SL(2,C)构造核式(3)及组合律。

## 4. 自主qubit运输的充要条件

在给定动量域上允许全部正常联合态时，存在仅依赖spin边缘的固定预测映射，满足

$$\operatorname{Tr}_p(U_A\rho U_A^\dagger)
=T(\operatorname{Tr}_p\rho)\quad(\forall\rho),\tag{4}$$

**当且仅当Ad_W(p)几乎处处相同。** 必要性不要求T线性：对|f><f|⊗τ，左侧是∫|f|²Ad_W(τ)dμ。其对所有f相同，迫使积分核几乎处处恒定；用张成Herm(2)的有限组τ完成证明。充分性给同一酉CPTP通道。

若先给已知、与spin独立的固定动量分布，可得到该分布加权的CPTP通道；它不具有式(4)的任意输入量词。下例进一步排除“分别保存两个边缘就足够”。

## 5. 正常有限波包的严格概率差

取m=1、x向boost ξ=ln3，动量中心p±=(0,0,±4/3)。中心处

$$W_\pm=(4I\mp i\sigma_y)/\sqrt{17}.\tag{5}$$

取任意归一C∞波包f+支持在p+半径ε=1/100的动量球中，f−(p)=f+(−p)。具体紧支指数轮廓见解析稿§4。比较

$$
\rho_{a,b}=\tfrac12|f_+\rangle\langle f_+|\otimes P_z^{\pm}
          +\tfrac12|f_-\rangle\langle f_-|\otimes P_z^{\mp}.
\tag{6}
$$

两spin边缘都是I/2，完整动量边缘也相同。它们是正常rank-2态，所有动量矩及能量有限；不使用精确动量本征态。位置波包没有紧支，未认证局域制备。

将包内W替换成中心W的比较基准，换参考后的Px+概率分别为25/34和9/34，差8/17。有限宽度不能忽略：由p↦v的导数及式(3)的归一化四元向量直接得全局界

$$
\|W(p)-W(p_0)\|\le L_{m,\xi}|p-p_0|,
\qquad L_{m,\xi}=\frac{|\tanh(\xi/2)|}{2m[1-|\tanh(\xi/2)|]}.
\tag{7}
$$

当前L=1/2。每份输出概率误差≤Lε，故真实包的差至少8/17−1/100=783/1700。两边缘相同的任何预测器最坏误差至少该差的一半。这个界由解析证明承担，正常波包的数值积分只校准实现。

## 6. 完整菜单需要什么资料

正常态限制到动量对角代数后由可积正矩阵Fρ(p)表示。任意boost方向报告为

$$\int d\mu\,\operatorname{tr}[F_\rho(p)W(A,p)^\dagger E W(A,p)].\tag{8}$$

因此保Fρ足以正确运输本菜单；只存其迹与积分两个边缘不够。不同动量之间的相干不被这份菜单看到，但可能影响位置、局域交互或其它后继任务，不能从完整态删除。

若另要求忠实保所有这些效果**生成的von Neumann代数**，则

$$\mathcal M=L^\infty(\mathbb R^3,d\mu)\bar\otimes M_2.\tag{9}$$

证明要点：各效果都是乘法矩阵，常Pauli在其中；第i轴boost在零点的算子范数导数为D_i=Σ_{j≠i}v_jσ_j。于是{σ_j,D_i}/2=v_jI。p↔v是Borel双射，三个v_j的联合谱演算生成全部动量乘法函数；再乘Pauli完成式(9)。

式(9)只分类数学代数闭包。证明使用乘积及弱极限，**没有证明逐项单次统计自动提供这些任务，或顺序仪器可由真实材料实现**。也未将连续统计限制当成原Hilbert空间上精确、无扰、正常的动量退相干装置。

## 7. 有限任务与物理读口

给有限动量箱、代表p_k，域内|p−p_k|≤r，保联合矩阵F_k=∫BkFρ并对尾箱保迹。若尾权重≤η，则

$$\tfrac12\|\tau_A-\widehat\tau_A\|_1\le L_{m,\xi}r+\eta.\tag{10}$$

将旧参考腿一起保留时，同界有效。有限任务因此不必先取得无限精度动量记录。真实分箱仪器、参考准备与资源代价仍需另接；本文没有优化它们。

同一实验的状态与效果共同运输时概率不变。式(6)是两种不同准备在同一个运动参考菜单中的比较，不是“同一事件不同观察者得到矛盾概率”。canonical spin也不能未经模型识别直接叫现实电磁spin读口；[Saldanha—Vedral](https://arxiv.org/abs/1111.7145)与[量子参考系操作研究](https://arxiv.org/html/1811.08228)提醒必须说明实际测量耦合。

## 8. 复算结果与核验

[代码](1038/moving_direction_transport.py) · [结果](1038/moving_direction_transport_results.json) · [核验器](1038/verify_round1038.py) · [核验收据](1038/research_round_1038_checks.json) · [审阅范围](1038/review.md)

|检验|结果/含义|
|---|---|
|36份完整诱导表示与闭式交叉|质量壳、酉性、组合及共同字典误差在浮点精度内；差分导数最大残差3.68×10⁻¹¹|
|正常动量包两套求积|细网格概率差0.470585852807；解析保证≥783/1700|
|精确四动量块|有理Lagrange投影，数学代数维16；不以样本证明无限闭包|
|任意参考的有限箱校准|3维参考、60维联合输入，输出半迹距离约0.000132807，低于解析预算|

默认复算只读，`--write`仅允许结果首次创建。执行：

```powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 research_cognition_physics/archive_1009_/1038/moving_direction_transport.py
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 research_cognition_physics/archive_1009_/1038/verify_round1038.py
```

## 9. 对生成目标的净变化

给定完整物理表示以后，方向比较律与所需联合资料不再能独立任意配置；强裸qubit运输被有限预测反例排除。该结论补旧固定参考方向接口的适用边界，未推导时空维数、Lorentz群、质量、规范物种或引力。

这是自由一粒子条件工具，尚未接成完整相互作用P981、量子来源/反作用或自主测量过程；不把高层采用条件重新称作认知公理。见[本地依赖账](1038/input_dependency_update.md)及[下一步](1038/NEXT.md)。完成本轮后停止该波包/boost候选扩建。
