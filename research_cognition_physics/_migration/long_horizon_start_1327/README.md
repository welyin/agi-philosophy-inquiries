# 长程目标启动后的当前导航层

2026-10-19。用户明确要求启动新目标；[目标正文](../../archive_1323_/_shared/long_horizon_goal.md)及[启动工具回执](../../archive_1323_/_shared/long_horizon_goal_start_receipt.json)已经保存。本次仅在五份当前导航标题后增加启动说明，科学轮次0、新公理0，下一可用编号1327。

[清单](manifest.json)登记可逆插入、前后哈希及原迁移清单哈希。原文不覆盖，不维护导航副本或ZIP。历史复核先在内存撤销此次插入，再恢复[上次撤回论文的迁移](../reclassify_1326_interim/README.md)，然后进入原科学证据层。旧验收的狭窄工作包含义不变，不自动认定新目标完成。

只读入口保持不变：

```text
python -B -X utf8 scripts/run_research_interaction_rules_current.py --verify-reclassification
python -B -X utf8 scripts/run_research_interaction_rules_current.py --verify-navigation
python -B -X utf8 scripts/run_research_interaction_rules_current.py --script research_cognition_physics/archive_1323_/_shared/check_stage_joint.py -- --check
```

[核验记录](verification.json)只检查导航、哈希和已有科学证据，不重复科学矩阵实验。启动回执的 active 是启动时事实，非未来永远运行的断言。
