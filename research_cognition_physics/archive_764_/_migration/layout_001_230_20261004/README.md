# 001—230轮目录整理与末阶段更名

2026-10-04。本次是文件组织维护，不增加研究轮次或科学结论。

## 变更

- 001—222进一步拆分为11个[主题阶段](../../../archive_223_230/001—230阶段总览.md)；223—230仍保留原边界。230篇正式报告放相应阶段根目录，配套材料归入数字目录。
- `archive_764_775` 更名为 `archive_764_`，后续仍在此研究。下一轮入口为 [776/drafts/STATUS.md](../../776/drafts/STATUS.md)。
- 共用背景与原导航单列；219修订前版本、论文旧版和原核验保留，重复导航副本与迁移ZIP已清理。

[1045份早期文件的迁移清单](manifest.json) · [快照哈希](snapshot.json) · [全部775轮及链接检查](layout_checks.json) · [早期全量单元测试](replay_checks.json)

## 字节与路径合同

代码、结果和历史收据保持原字节。阅读版Markdown只更新链接目标；精确替换位置记入清单。原导航的完整改写另存前版本和变更收据。冻结稿、原哈希清单中的旧路径按写作时布局理解。

231—775的[当前迁移清单](../manifest.json)已同步目的路径，ZIP仅是当时的备份，现已由直接读取当前文件的[复算入口](../README.md)取代。本次不重新验证数学论证，也不把两条已有整链审计限制改成通过，详见[此前复算范围](../README.md)。

## 复算

从项目根目录用既有Python/NumPy运行：

```powershell
python -B -X utf8 research_cognition_physics/archive_764_/_migration/layout_001_230_20261004/replay_early.py
python -B -X utf8 research_cognition_physics/archive_764_/_migration/layout_001_230_20261004/replay_early.py --stage 1 --module local_calibration_contract
python -B -X utf8 research_cognition_physics/archive_764_/_migration/replay.py --script verify_round775.py
```

早期入口通过迁移清单直接读取现有阶段目录的代码，只运行原单元测试，不解压ZIP、不恢复旧目录、不调用结果写入器。默认检查第一阶段原218个测试模块的1915项测试与第二阶段八轮70项测试；SciPy复用项目已有 `.research_runtime`。可用 `--stage 1` 或 `--stage 2` 限定范围，`--report` 只允许创建新的收据。

历史发布核验要求当时的完整导航和论文哈希；本次单元测试与当前路径完整性检查分别记录，不把历史发布合同默默套到新阅读版。

## 追加分期的追溯

[分期快照哈希](phase_split_snapshot.json)、[逐文件移动计划](phase_split_plan.json)和[链接变更](phase_split_link_changes.json)记录本轮用户追加要求。重复整包和导航副本已列入清理；保留迁移清单、哈希及核验收据。
