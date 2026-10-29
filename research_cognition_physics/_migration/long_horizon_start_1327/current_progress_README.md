# 当前进展区：单份维护，不新增导航副本

五份导航仅在首标题后保留一块 `CURRENT_RESEARCH_PROGRESS_BEGIN/END` 标记区。[当前进展清单](current_progress_manifest.json)绑定整文、区块、剥离后基底和实际科学证据。每次发布替换这同一块和这同一份清单，旧科学报告与阶段快照承载历史；不复制整份导航，不新增逐轮迁移链。

两份旧迁移清单保持冻结。只读兼容器先校验进展证据，剥离本区，再恢复目标启动层、撤稿层及旧科学核验层。任何区外改动、标记重复、证据哈希不符均拒绝；本区不能改变旧科学字节。当前进度的有效性来自已保存验收，哈希检查自身不是重新证明科学结论。

复核入口仍为 `scripts/run_research_interaction_rules_current.py` 的 `--verify-reclassification`、`--verify-navigation` 和旧 `check_stage_joint.py --check`。本轮核验见[当前复核记录](current_progress_verification.json)。更新时应先确认新的验收证据，再替换动态区并重算同一清单、同一复核记录；不得改写两份冻结的 `manifest.json`。
