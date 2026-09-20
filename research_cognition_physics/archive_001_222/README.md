# 第1—222轮研究档案

归档日期：2026-09-20。按用户决定结束本阶段，转为公理化工作定义及后续过程研究。

## 阅读入口

- [新的阶段论文](../可组合认知结构与复量子状态空间_阶段论文.md)。
- [全部222轮：结论与适用范围](ROUND_INDEX.md)。
- [原研究索引](research_process/README.md)、[222轮最终笔记](research_process/research_note_222.md)。
- [219轮修订稿](research_process/research_note_219.md)、[219轮原稿归档清单](research_process/history/round219_v1/archive_manifest.json)。
- [文件哈希清单](ARCHIVE_MANIFEST.json)。

## 研究过程与归档范围

`research_process/`表示本阶段的研究过程，保存721个研究文件，含222篇编号笔记、221份主目录Python脚本、结果、账本、原导航和219轮历史版本。既有缓存随目录保存，不计科研证据。

2026-09-20按用户明确要求，删除归档内的`research_physics_construction/`与`research_information_geometry/`两个副本目录，共102份引用副本；原`research_cognition_physics/`重命名为`research_process/`。两条独立路线仍使用项目根目录的[物理构造资料](../../research_physics_construction/README.md)与[信息几何资料](../../research_information_geometry/research_note_01.md)，不纳入本档案。

## 文件布局

| 位置 | 归属与用途 |
|---|---|
| `research_process/` | 第1—222轮研究过程与结果 |
| 外层3份讨论／架构文档 | 其余外部阅读快照，不计研究轮次 |
| `migration_backups/` | 迁移前文件备份及4份链接调整前的研究原文，供追溯 |
| `ROUND_INDEX.md` / `ROUND_INDEX.json` | 指向新目录的222轮索引 |
| `ARCHIVE_MANIFEST.json` / `PATH_MIGRATION.json` | 文件校验、已删除副本的元数据及路径迁移记录 |
| `run_stage1.py` | 运行历史单元测试；不执行结果写入器 |
| `verify_stage1.py` | 检查原始证据哈希、222轮覆盖、索引与本地链接 |
| `stage1_test_results.json` | 阶段归档时的完整数值回归记录：1915项测试 |
| `stage1_archive_checks.json` | 最近一次归档完整性与链接检查 |
| `stage1_compatibility_checks.json` | 初次迁移时与独立构造路线的兼容性检查记录 |
| `stage1_layout_checks.json` | 删除两个副本目录及改名时的24项检查记录 |
| `paper_math_review_checks.json` | 阶段论文的公式语法、编号、内容保持及渲染审核记录 |

717份研究文件保持原字节；4份Markdown只调整指向两条独立路线的链接，其原文保存在`migration_backups/process_link_originals/`。代码、实验结果、数学论证与历史结论未改写。清单保留原始哈希，并为链接调整后的版本单列当前哈希和可复核的替换记录。

旧笔记中的研究方向与“下一步”按历史原样保留。本次之后以[当前研究状态](../RESEARCH_STATE.md)为准。历史证书中的逻辑路径和旧哈希保持不变，复算及校验入口负责解析到新目录，并在需要时验证链接调整前的原文。

辅助脚本及检查记录集中放在本归档目录；上一层仅保留三个导航文档与阶段论文。这些辅助文件不属于新增研究轮次，已有测试报告记录的是当时实际执行的检查。

## 复算与保护

在项目根目录使用[运行入口](run_stage1.py)执行历史单元测试，它补齐同目录模块及现有`.research_runtime`依赖搜索路径，不执行结果写入器：

```powershell
python research_cognition_physics/archive_001_222/run_stage1.py local_calibration_contract
python research_cognition_physics/archive_001_222/run_stage1.py --all
python research_cognition_physics/archive_001_222/verify_stage1.py
```

实验修改和重新生成结果应放入新的轮次文件，避免覆盖归档材料。哈希检查验证原始证据及明示的链接调整，不替代数学证明。历史中已存在的失效链接单独列入清单，迁移不得新增失效链接。
