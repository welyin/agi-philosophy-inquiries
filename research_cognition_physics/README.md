# 认知物理研究：阶段论文与物理生成

更新：2026-09-22。**第一阶段1—222轮、第二阶段223—230轮已结项；第三阶段完成231—232轮。**

## 当前入口

- **[第二阶段整合论文 v1.1](可组合认知结构与有限维量子理论_阶段论文.md)**：229—230轮补充已合入定义、综合定理、证明与实验误差范围。结项范围是有限维量子操作规则与连续可逆动力学形式的条件性重建。
- **[第三阶段研究目录](archive_231_/README.md)**、[用户的物理生成猜想](archive_231_/physical_generation_conjecture.md)。
- **[第231轮：事件次序审计](archive_231_/research_note_231.md)**：固定操作下影响不传递；在明示无环线路中，可能影响等于可达性；可达性仍不唯一确定钟尺。8项检查通过。
- **[第232轮：记录与组合一致性](archive_231_/research_note_232.md)**：相容记录不迫使固定全局次序；已知单次一致模型在双副本联合操作下失败。8项检查与具体解析反例。下一步审计合法端口合并。
- [研究方向](research_direction.md)、[当前状态及下一步](RESEARCH_STATE.md)。

## 已完成阶段

| 阶段 | 成果 | 原始证据 |
|---|---|---|
| 1—222 | [复量子状态空间论文](可组合认知结构与复量子状态空间_阶段论文.md)：F＋U＋C＋P工作定义、P⇔L及复矩阵状态锥的定理衔接 | [档案](archive_001_222/README.md)、[逐轮索引](archive_001_222/ROUND_INDEX.md)、[完整性核验](archive_001_222/verify_stage1.py) |
| 223—230 | [有限维量子理论论文](可组合认知结构与有限维量子理论_阶段论文.md)：相容张量表示、CP仪器、操作完备性、连续可逆动力学形式及有限实验误差 | [档案](archive_223_230/README.md)、[当前核验入口](archive_223_230/verify_archive.py)、[原v1.0版本](archive_223_230/paper_versions/finite_quantum_v1.0.md) |

第二阶段原合同给全部有限仪器任意逼近；某个对象的非恒定连续可逆路径R_seed足以给精确权限。另加Time给薛定谔形式，非平凡Time还提供种子。原合同是否自动提供该路径仍开放。未选定具体Hamiltonian、物理钟尺、无限维极限、场论或引力。

按用户新指示，旧archive_223_已重命名为archive_223_230；旧论文和补充稿保留作历史版本，当前阅读以整合论文为准。历史笔记的待办按当时语境理解。

## 工作约定与目录

以选定数学条件定义本文的可组合认知结构，暂不重开定义的普适认知解释。区分认知动机、新增输入、解析证明、数值验证和物理解释；引用成熟定理时逐项核验前提。相似结构、共同维数或一个可实现模型不足以证明理论等价或必然性。

顶层文件只保留本README、research_direction.md、RESEARCH_STATE.md与阶段论文。所有代码、结果和过程笔记位于阶段子目录。第一阶段过程目录是archive_001_222/research_process；独立[物理构造](../research_physics_construction/README.md)与[信息几何](../research_information_geometry/research_note_01.md)不再复制进本档案。

## 复算

在项目根目录使用既有Python与NumPy运行时：

    python -B -X utf8 research_cognition_physics/archive_001_222/verify_stage1.py
    python -B -X utf8 research_cognition_physics/archive_223_230/verify_archive.py --run-tests
    python -B -X utf8 research_cognition_physics/archive_231_/event_order_audit.py
    python -B -X utf8 research_cognition_physics/archive_231_/verify_round231.py
    python -B -X utf8 research_cognition_physics/archive_231_/verify_round232.py

核验不覆写旧科学结果。第二阶段70项检查与新阶段两轮各8项检查分别计数；当前公式渲染与目录迁移记录见各阶段索引。
