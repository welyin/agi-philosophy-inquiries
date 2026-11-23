# 原第16论文撤回与中间成果归位

2026-10-19，按用户要求将原根目录第16阶段论文移到[中间成果](../../archive_1323_1352/未知参考交接的作用选择与动态来源_中间成果.md)，取消编号及完整阶段结项定位。实际研究结果仍保留；没有新科研轮次、公理或应用目标。

当前引用及研究状态已修正。旧验收JSON、旧程序与应用目标完成记录保留为当时的窄工作包证据，里面的旧路径和旧哈希不被篡改为新的历史。原“阶段已完成”不再代表现行研究分类；以[当前研究目录](../../archive_1323_1352/README.md)及[长程目标建议](../../archive_1323_1352/_shared/next_long_horizon_proposal.md)为准。

[迁移清单](manifest.json)仅登记移动和可逆文字修改。原论文不在根目录保留副本，也不生成ZIP。历史复核时读取实际新位置，验证当前内容，再在内存恢复迁移前字节。

当前只读入口为[scripts/run_research_interaction_rules_current.py](../../../scripts/run_research_interaction_rules_current.py)：

```text
python -B -X utf8 scripts/run_research_interaction_rules_current.py --verify-reclassification
python -B -X utf8 scripts/run_research_interaction_rules_current.py --verify-navigation
python -B -X utf8 scripts/run_research_interaction_rules_current.py --script research_cognition_physics/archive_1323_/_shared/check_stage_joint.py -- --check
```

它不恢复磁盘旧论文、不写研究结果，也不重新发布旧结项。涉及历史论文／说明哈希的旧核查通过此入口执行；只涉及未改科学文件的逐轮程序仍可直接按原环境运行。
