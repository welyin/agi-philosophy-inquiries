# 认知物理研究：阶段论文与后续起点

更新：2026-09-20。第1—222轮已告一段落，按用户要求归档；本次整理不另计第223轮。

## 当前成果

- [阶段论文：可组合认知结构与复量子状态空间](可组合认知结构与复量子状态空间_阶段论文.md)：工作定义、定理链、主题综合、全部222轮索引及后续任务。
- [原始研究档案](archive_001_222/README.md)：笔记、代码、结果、严格证书与219轮修订历史。
- [逐轮结论与范围](archive_001_222/ROUND_INDEX.md)、[原文件校验清单](archive_001_222/ARCHIVE_MANIFEST.json)。
- [迁移与复算检查](archive_001_222/stage1_archive_checks.json)、[归档代码回归结果](archive_001_222/stage1_test_results.json)。

## 已接受的工作方向

以有限维因果凸操作框架F、同操作维数的普遍系综导引U、连续可逆纯态传递C、完整动作本地校准可迁移P，定义本文研究的可组合认知结构。在该框架内P等价于局部层析L，既有重建定理给出复矩阵状态锥。

这是一套明确的充分性路线与工作定义。定义的普适认知解释暂不作为继续工作的门槛；不把全部旧模型的条件同时采用，也不宣称已经得到全部物理动力学。

下一阶段采用定理衔接路线：核验成熟定理的全部前提，先从状态结构走向量子过程与连续动力学，再推进场论和Jacobson等引力路线；弦论作为可选比较方向。见[当前研究方向](research_direction.md)和[当前状态](RESEARCH_STATE.md)。历史笔记中的“下一步”记录的是当时计划，不覆盖本次决定。

归档的研究过程目录为`archive_001_222/research_process/`。按用户要求，已删除归档内的物理构造、信息几何两个副本目录，相关引用直接指向项目根目录的独立研究资料，详见[归档范围](archive_001_222/README.md#研究过程与归档范围)。

顶层文件只保留本README、`research_direction.md`、`RESEARCH_STATE.md`及阶段论文。测试入口、校验脚本与JSON记录已移入`archive_001_222/`；后续研究代码、结果和过程记录放入相应阶段子目录。

## 复算

在项目根目录，使用已有Python与NumPy运行时：

```powershell
python research_cognition_physics/archive_001_222/run_stage1.py local_calibration_contract
python research_cognition_physics/archive_001_222/run_stage1.py --all
python research_cognition_physics/archive_001_222/verify_stage1.py
```

入口只执行单元测试，不调用历史结果写入器。需要修改或重新生成实验时，从档案复制到新的研究文件后再执行。另一个[物理构造目录](../research_physics_construction/README.md)保留独立的C0—C13合同，未合并为这222轮的结论。
