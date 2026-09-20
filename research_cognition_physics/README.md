# 认知物理研究：阶段论文与后续起点

更新：2026-09-20。第1—222轮已归档；第二阶段在`archive_223_/`推进，已完成第226轮。

## 当前成果

- [第二阶段入口](archive_223_/README.md)、[第226轮：纯导引资源与辅助测量权限](archive_223_/research_note_226.md)：同维辅助实现满秩态的全部有限系综，迫使纯满Schmidt秩资源及该辅助的全部POVM权限。8项检查通过。
- [第225轮：全局取向与跨系统完全正性](archive_223_/research_note_225.md)：由三体正边缘统一全部矩阵取向，将CP结论推广到所有实际A→B过程；尚不授予全部CP实现权限。9项检查通过。
- [第224轮：张量复合与完全正性](archive_223_/research_note_224.md)：合法边缘与维数饱和识别二元张量结构，推出实际A→A过程完全正。11项检查通过。
- [第223轮：从复状态空间到连续可逆动力学](archive_223_/research_note_223.md)：连续可逆过程属于酉分支；加入明示的连续齐次时间群合同后得到交换子生成元。10项检查通过。
- [阶段论文：可组合认知结构与复量子状态空间](可组合认知结构与复量子状态空间_阶段论文.md)：工作定义、定理链、主题综合、全部222轮索引及后续任务。
- [原始研究档案](archive_001_222/README.md)：笔记、代码、结果、严格证书与219轮修订历史。
- [逐轮结论与范围](archive_001_222/ROUND_INDEX.md)、[原文件校验清单](archive_001_222/ARCHIVE_MANIFEST.json)。
- [迁移与复算检查](archive_001_222/stage1_archive_checks.json)、[归档代码回归结果](archive_001_222/stage1_test_results.json)。

## 已接受的工作方向

以有限维因果凸操作框架F、同操作维数的普遍系综导引U、连续可逆纯态传递C、完整动作本地校准可迁移P，定义本文研究的可组合认知结构。在该框架内P等价于局部层析L，既有重建定理给出复矩阵状态锥。

这是一套明确的充分性路线与工作定义。定义的普适认知解释暂不作为继续工作的门槛；不把全部旧模型的条件同时采用，也不宣称已经得到全部物理动力学。

第二阶段采用定理衔接路线：核验成熟定理的全部前提，从状态结构走向量子过程与连续动力学。223—225轮已给可逆动力学形式、相容复合及全部实际过程的CP表示，226轮补上指定辅助的全部POVM权限；下一步核验实际可逆控制能否支持全部量子操作，随后推进场论和Jacobson等引力路线。弦论作为可选比较方向。见[当前研究方向](research_direction.md)和[当前状态](RESEARCH_STATE.md)。历史笔记中的“下一步”记录的是当时计划，不覆盖当前决定。

归档的研究过程目录为`archive_001_222/research_process/`。按用户要求，已删除归档内的物理构造、信息几何两个副本目录，相关引用直接指向项目根目录的独立研究资料，详见[归档范围](archive_001_222/README.md#研究过程与归档范围)。

顶层文件只保留本README、`research_direction.md`、`RESEARCH_STATE.md`及阶段论文。测试入口、校验脚本与JSON记录已移入`archive_001_222/`；后续研究代码、结果和过程记录放入相应阶段子目录。

## 复算

在项目根目录，使用已有Python与NumPy运行时：

```powershell
python research_cognition_physics/archive_001_222/run_stage1.py local_calibration_contract
python research_cognition_physics/archive_001_222/run_stage1.py --all
python research_cognition_physics/archive_001_222/verify_stage1.py
python -B -X utf8 research_cognition_physics/archive_223_/reversible_dynamics_bridge.py
python -B -X utf8 research_cognition_physics/archive_223_/tensor_process_bridge.py
python -B -X utf8 research_cognition_physics/archive_223_/global_orientation_bridge.py
python -B -X utf8 research_cognition_physics/archive_223_/steering_permission_bridge.py
```

入口只执行单元测试，不调用历史结果写入器。需要修改或重新生成实验时，从档案复制到新的研究文件后再执行。另一个[物理构造目录](../research_physics_construction/README.md)保留独立的C0—C13合同，未合并为这222轮的结论。
