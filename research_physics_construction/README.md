# 认知物理构造路线

建立日期：2026-09-17。当前版本：**C13有限径向多模与单模误差，第13轮完成**；C0—C12保留。

本路线从已有研究继续向量子理论与广义相对论的统一模型推进。方法改为：**固定可修订的工作假设，构造一个模型，同时检验多类现象；让具体成功、失败和新数据反过来约束前提。** 不再把证明每项假设唯一必然，作为每一步构造的开工条件。

**用户确认的当前约束：只论证充分性，必要性只记录、不论证。** 先从明确认知模型建立通向量子理论与广义相对论的完整充分性链条，再回头讨论前提为什么不可替代。必要性事项统一登记在[暂存清单](NECESSITY_BACKLOG.md)，不因局部失败自动重开；内部一致性、物理符合性与输入—结论核对仍是必需验收。

旧路线的216轮认知研究与8轮信息几何研究完整继承为证据库，不删除、不复制覆盖，也不自动继续其最后的公理追问。新路线独立编号，从00开始，不使用旧编号217。

## 从这里开始

1. [起点整合](FOUNDATIONS.md)：全部研究阶段怎样进入新路线，哪些是约束、工具、工作假设与未决问题。
2. [C0模型约定](MODEL_CONTRACT.md)：同一个模型的状态、操作、记录、反馈、组合与资源边界。
3. [逐轮继承索引](INHERITED_INDEX.md)：224篇编号笔记的原结论、条件和源文件链接。
4. [当前状态与路线图](RESEARCH_STATE.md)：下一轮的具体任务和停止条件。
5. [第00轮整合记录](research_note_00.md)：可运行基线、结果及其物理解释边界。
6. [C1传播约定](C1_CONTRACT.md)与[第01轮](research_note_01.md)：同一局部动力学中的传播、干涉和消息延迟反馈。
7. [C2钟与标定约定](C2_CONTRACT.md)与[第02轮](research_note_02.md)：两个有限钟、计数区间、读数回传和经标定反馈。
8. [C3共享背景约定](C3_CONTRACT.md)与[第03轮](research_note_03.md)：非均匀链路、两类已标定探针、共同参数见证及未参与拟合读数的预测。
9. [C4长链约定](C4_CONTRACT.md)与[第04轮](research_note_04.md)：从同一中点交换规则得到双向变速包络，并给连续单元记录的二阶误差界。
10. [C5质量约定](C5_CONTRACT.md)与[第05轮](research_note_05.md)：局部交错质量耦合两个载波，给连续参考截断、变质量传播及平滑终端读数保证。
11. [C6量子背景约定](C6_CONTRACT.md)与[第06轮](research_note_06.md)：同一有质量链、两个钟和一个集体背景模联合演化，给双向响应、准备一致性、能量账及模截断记录预算。
12. [C7局域背景约定](C7_CONTRACT.md)与[第07轮](research_note_07.md)：相邻两模只耦合附近数据及本站钟，隔离背景传递路径，并给两模截断控制的远端源记录差异。
13. [C8背景场约定](C8_CONTRACT.md)与[第08轮](research_note_08.md)：显式采用梯度振子生成元，得到标量有效方程与有限读数尺度；条件高斯联合态另行检验冻结源和钟。
14. [C9弱场度规约定](C9_CONTRACT.md)与[第09轮](research_note_09.md)：采用明确四维二次张量作用量，构造守恒孤立源、张量潮汐及同度规双物种和钟记录，保留非线性残差。
15. [C10流体度规约定](C10_CONTRACT.md)与[第10轮](research_note_10.md)：在声明的非线性引力与因果EOS下求压力支撑的球体，核对协变平衡、真空匹配和完整度规探针记录。
16. [C11径向扰动约定](C11_CONTRACT.md)与[第11轮](research_note_11.md)：同一EOS球的径向线性稳定性有独立参数能量证书，数值模谱、移动表面守恒和时变钟读数分别核对。
17. [C12径向量子模约定](C12_CONTRACT.md)与[第12轮](research_note_12.md)：由物理惯性和量子尺度检查出发，明确宏观同形源族，构造模、数据与双钟的自主联合能量交换。
18. [C13多模约定](C13_CONTRACT.md)与[第13轮](research_note_13.md)：前两个泛音的真空源响应、有限三模联合演化，以及对C12原有记录的单模遗漏误差。

## 已经落地

- [继承清单](inherited_results.json)：224篇编号笔记、705份源材料，包含223份Python文件、252份JSON；旧可变导航文件不参与哈希锁定。逐轮结论与适用边界均保留，历史结果不是本轮重新证明。
- [统一基线](baseline.py)及[结果](baseline_results.json)：一个四槽联合状态引擎，贯通干涉、Bell关联、公开记录、预测反馈、退相干对照及声明网络上的传递。
- 原复合整体和共同J实表示逐操作对应；局部实摘要作为受限接口，未实施“删除虚部”通道。
- 本轮不声称从认知推导了Born规则、量子必然性、物理时间、洛伦兹对称或引力。

**第01轮新增**：[传播模块](propagation.py)与[结果](propagation_results.json)。同一个局部交换生成元实现端点传播与干涉；记录经M→B→A两跳送达，控制器在等待期间更新模型。延迟测试中到达时策略为25/32，过时动作为23/32；一阶近似有解析误差上界并足以保持示例中的行动选择。14项新检查通过，物理解释按工作输入与充分性结果分别登记。

**第02轮新增**：[有限钟标定模块](clock_calibration.py)与[结果](clock_calibration_results.json)。在声明的同步准备、稳定参数和相位窗口内，152,020次模拟新试验给出两端相容的相对刻度区间；解析并合界保证同时覆盖概率至少99%。估计器只接收计数，标定行动码和新任务记录均计传输延迟，后续任务成功率约0.958732。13项新检查通过；同步、钟动力学和工作窗口是输入，不是已推导的物理时空。

**第03轮新增**：[共同背景模块](shared_background.py)与[结果](shared_background_results.json)。339,056次模拟新试验使两类探针的背景区间同时获得至少99%的条件覆盖保证；计数给出的同一组参数可以放回六槽引擎运行，全部24个均值残差小于0.00385。端点观测量未用于拟合，预测区间还指导一次新反馈，成功率约0.997417。几何是明确插值与局部色散标尺下的候选表示，共享耦合是输入，尚未建立连续场或引力方程。

**第04轮新增**：[包络极限模块](envelope_limit.py)与[结果](envelope_limit_results.json)。在新声明的平滑周期背景上，C3同一中点耦合导出双向方程中的v'/2修正；Taylor与Duhamel给固定平滑准备、有限时间的O(a²)误差。256节点、参数2时连续单元记录TV上界约0.002816；128节点的包络预测分类成功率0.727005，格点结果0.726909，读数收齐前禁止判断。14项新检查通过。结果是静态背景上的无质量单粒子有效传播，未给质量、反作用或Einstein动力学。

**第05轮新增**：[有质量模块](massive_envelope.py)与[结果](massive_envelope_results.json)。同一链的交错局部项精确交换载波，得到带m ell sigma_x的有效方程；质量与钟率耦合是工作输入。有限Fourier参考具有独立截断界，变系数正则性和格点残差给有限时间误差。512节点平滑终端任务的预测正确率0.590022、格点结果0.590018，总误差上界约0.023892；14项新检查通过。原始微观位置干涉不被隐藏为连续密度，背景反作用和Einstein方程仍未构造。

**第06轮新增**：[量子背景模块](quantum_chain_background.py)与[结果](quantum_chain_background_results.json)。在同类16节点链上保留两个相干钟和一个10级集体模，单一哈密顿量实现数据与背景、钟的双向响应；自由总能量及准备混合一致性保持。T=2的模截断记录TV预算约0.00012449，背景非真空概率扣除该预算后仍至少约0.00264083。19个原始结果经72次bit-hop在参数92收齐，再发布联合输出。13项新检查、路线共97项通过。集体模尚非局域引力场，C5静态连续误差界不能自动沿用，Einstein动力学仍未得到。

**第07轮新增**：[局域背景模块](local_quantum_background.py)与[结果](local_quantum_background_results.json)。8节点同类链的0、1节点各有一个模与钟，模只调制本站质量、钟及相邻边，两模间有明确背景交换。暂停数据跃迁的诊断隔离源到另一处钟的背景路径；每模10级时，概率差约0.00030433，扣除两次准备的预算后仍至少约0.00021675。8级预算不足的结果保留。默认跃迁开启时仍有双向作用及完整能量账；12个原始结果经18bit-hop在参数47收齐。16项新增、路线共113项通过。两节点片段尚非连续场或Einstein动力学。

**第08轮新增**：[背景场模块](background_field_limit.py)与[结果](background_field_limit_results.json)。显式将数守恒交换改为动量平方与邻点位移差平方，固定平滑模式得到标量Klein–Gordon型方程和二阶空间界。有限分箱、整数权重与饱和预算一并计入；协同精化至128节点的完整记录误差界约0.00000746836，原始3202位、102401bit-hop，未隐藏精度代价。另一有限冻结源任务用完整条件高斯态保留钟相位、准备一致性和能量，无Fock截断。17项新增、路线共130项通过。自由场连续界不覆盖有源过程；新生成元为输入，尚无度规约束或Einstein方程。

**第09轮新增**：[弱场度规模块](weak_metric_response.py)与[结果](weak_metric_response_results.json)。在明示的3+1维二次度规作用量下，变分得到线性化Einstein方程，核对规范、Bianchi及正高斯孤立源；不将作用量输入称为认知推导。两种质量的原有链与两个钟使用同一源度规，归一化近远钟率比约0.99258044。每次34个原始位、272bit-hop，参数175.656389收齐。16项新增、路线共146项通过。完整非线性残差仍非零；源支持应力、量子源与完整引力闭合尚未完成。

**第10轮新增**：[非线性流体模块](nonlinear_fluid_metric.py)与[结果](nonlinear_fluid_metric_results.json)。先用均匀密度解析球核对完整曲率与积分器，再用声速平方1/3的EOS求出压力支撑的非线性平衡球，R约4.134864、M约3.130386，表面p=0并接真空。独立焓–钟率关系、密集输出差分和精化检查与场方程相容；不把代回残差当全局误差证书。完整面积半径度规给近远钟率比约0.92287658，每试验34原始位、272bit-hop，参数212.210453收齐。17项新增、路线共163项通过。非线性作用量和EOS仍为输入，径向稳定性与量子源闭合待做。

**第11轮新增**：[径向模模块](radial_fluid_modes.py)与[结果](radial_fluid_modes_results.json)。保持C10参数不变，正确采用有限表面密度下的自由压力边界，独立参数能量界和精确有理比较给所有允许径向线性模omega²>0.015；正Ritz值不单独当作证明。最低模omega²约0.1256416728，与射击、节点和精化结果一致。粒子数和ADM质量的移动表面项保留；时变探针能量变化与做功积分相符，完整记录相对静态TV约0.00019715。17项新增、路线共180项通过。结论不覆盖非径向或非线性稳定，模耗尽与量子探针反作用仍未闭合。

**第12轮新增**：[径向量子反作用模块](quantum_radial_backreaction.py)与[结果](quantum_radial_backreaction_results.json)。物理惯性为4pi乘C11积分，并用固定粒子数ADM动能变化独立核对。原尺寸真空约93.3%落在线性窗口外，结果保留；明确放大1024倍的物理同形族后，构造自主单模–数据–钟H，总能量保持且模有双向响应。D20相同记录的Fock预算约1.65e-9，任一时刻幅度窗口外概率上界约0.007302，二者不混同非线性误差。17项新增、路线共197项通过。其他模、完整二阶几何及认知到动力学输入的连接未完成。

**第13轮新增**：[径向多模模块](radial_multimode_reduction.py)与[结果](radial_multimode_reduction_results.json)。保持C12物理尺度，加入两个独立真空泛音，核对惯性正交及非对易源的完整三模演化。对新增模取迹后，C12原有记录的有限遗漏上界约1.30661e-9，实际精化差约1.75472e-10；完整态与保留记录的误差界分开，不能只凭频率高就删除泛音。三模幅度窗口及Fock预算另记，保留513点预算不足与1025点精化结果。17项新增、路线共214项通过。无限高阶模尾部、非线性几何及认知动力学连接仍未完成。

第00轮检查见[验证记录](baseline_checks.json)：新路线16项加29项直接相关旧检查，共45项通过；不是旧研究全量重验。

## 复算

沿用现有Python与NumPy；C10的TOV积分另需SciPy，当前环境已安装。清单及文档工具使用已有pypandoc，不需重新建立环境。

```powershell
python -B -X utf8 research_physics_construction/inheritance.py --verify
python -B -X utf8 research_physics_construction/baseline.py --write-results
python -B -X utf8 research_physics_construction/propagation.py --write-results
python -B -X utf8 research_physics_construction/clock_calibration.py --write-results
python -B -X utf8 research_physics_construction/shared_background.py --write-results
python -B -X utf8 research_physics_construction/envelope_limit.py --write-results
python -B -X utf8 research_physics_construction/massive_envelope.py --write-results
python -B -X utf8 research_physics_construction/quantum_chain_background.py --write-results
python -B -X utf8 research_physics_construction/local_quantum_background.py --write-results
python -B -X utf8 research_physics_construction/background_field_limit.py --write-results
python -B -X utf8 research_physics_construction/weak_metric_response.py --write-results
python -B -X utf8 research_physics_construction/nonlinear_fluid_metric.py --write-results
python -B -X utf8 research_physics_construction/radial_fluid_modes.py --write-results
python -B -X utf8 research_physics_construction/quantum_radial_backreaction.py --write-results
python -B -X utf8 research_physics_construction/radial_multimode_reduction.py --write-results
python -B -X utf8 -m unittest discover -s research_physics_construction -p "*.py"
```

清单完整性检查会拒绝静默改变旧源、轮次范围或继承条件；它不是证明验证器。后续有必要更新起点时应建立明确的新版本与变更记录，而不是把新结果追写成C0原本就有的事实。

## 本路线的工作规则

- 只开展充分性构造与验证；必要性仅登记，当前不开展证明、替代模型排除或相应扫描。
- 所有模型共用明确的状态、混合、控制和组合合同；不同任务允许不同已知准备，但不能悄悄换理论。
- 建模输入与推导结果分开。复规则若已作为输入，复算出标准量子预测只能算内部一致性或工程验收。
- 每轮至少给一个可检验增量、复算结果及对统一模型的影响；参数扫描不自动拆成多轮。
- 不为了得到预期结论而隐藏参考、丢弃失败记录、改变输入独立性或免费重置资源。
- 历史反例继续约束新模型；阶段性解释缺口不自动阻止构造，违反概率一致性或已承诺任务的具体错误则必须修正。
- 不创建后台研究任务，不自动执行Git操作。用户的新方向优先于旧笔记中的“下一步”。

当前最重要的边界：**采用“整体复等价、局部允许实接口”作为工作假设，不等于已经证明世界只能如此。**