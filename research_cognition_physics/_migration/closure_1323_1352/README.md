# 第1323—1352轮闭档复算

在项目根目录用现有的 Python 和 NumPy 执行以下命令。第一条检查改名、论文、审阅及历轮证据，第二条重新运行30轮验收中登记的科学程序。

```powershell
python -B -X utf8 scripts/run_research_interaction_rules_current.py --verify-stage-closure
python -B -X utf8 scripts/run_research_interaction_rules_current.py --replay-stage-science
```

历次验收中写入的原目录名按相对路径映射到同一份[正式档案](../../archive_1323_1352/README.md)。闭档清单绑定了改名前414个文件的摘要；对改过正文的八篇 Markdown，另外记录新、旧摘要和恢复旧文所需的差异。复算在内存里重建验收时的原始字节。[闭档清单](manifest.json)

`--verify-reclassification` 和 `--verify-navigation` 仍检验它们各自的历史边界。`--verify-current-progress` 可核对第1352轮结束时已冻结的进展证据；其中的目标状态是该轮发布当时的记录。闭档后的论文及独立评审另行绑定。[本次结项成果](products.json)

可通过原入口单独复算，例如以下命令核对1352轮独立作者的结果。

```powershell
python -B -X utf8 scripts/run_research_interaction_rules_current.py --script research_cognition_physics/archive_1323_1352/1352/independent_check.py -- --check
```

