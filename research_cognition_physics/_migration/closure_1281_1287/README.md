# 1281—1287归档与1288启动的可核验迁移

2026-10-18按用户要求，将`archive_1281_`收尾为[archive_1281_1287](../../archive_1281_1287/README.md)，下一阶段进入[archive_1288_1289](../../archive_1288_1289/README.md)。本记录属于归档维护，不是新的科学轮次。

## 保留范围

- 原阶段100件文件全部移入新路径，包括7轮报告、代码、结果和独立审阅。
- 48件Python／JSON科学文件逐字节保留。Markdown仅修改迁移路径和归档导航；原字节可根据编辑清单恢复。
- 论文11更新目录引用；四份研究根导航先登记结项，再以单独启动层登记1288。历史阶段结论不扩大。
- 原目录不保留影子副本，不创建ZIP或符号链接；新复算入口只在内存映射旧路径。

## 证据

- [迁移前独立库存](independent_before.json)：逐文件路径、大小及SHA256。
- [最终迁移清单](manifest.json)：100件迁移文件、5件根文档；原始／当前字节摘要与可逆编辑。
- [独立归档审计](independent_audit.md)：逐项保全、反向恢复及链接检查。
- [复算验证](runtime_validation.json)：新入口和代表性历史作者／独立检查。
- [1288启动层](../structure_selection_active_1288/layer.json)：根导航单次前缀插入及封存依赖，不改科学档案。

## 当前复算方式

从仓库根目录，使用现有Python运行以下入口（无需安装新依赖）：

```text
python -B -X utf8 scripts/run_research_structure_selection_active.py --verify-navigation
python -B -X utf8 scripts/run_research_structure_selection_active.py --verify-closure
python -B -X utf8 scripts/run_research_structure_selection_active.py --verify-matter-navigation
python -B -X utf8 scripts/run_research_structure_selection_active.py --script research_cognition_physics/archive_1281_1287/1287/check.py
python -B -X utf8 scripts/run_research_structure_selection_active.py --script research_cognition_physics/archive_1281_1287/1287/independent_check.py
```

历史逻辑路径`archive_1281_`也可通过该入口使用。新阶段已经修改根导航，因此最终入口使用`structure_selection_active`；单独的旧冻结入口不负责恢复后加导航。

封存层分别撤销1288启动前缀、1281—1287迁移编辑、此前导航层，以原摘要核对历史材料；不是忽略摘要或把新文件冒充原文件。科学检查运行于禁止写研究档案的模式。

## 原有验收范围

原阶段的`COMPLETE_WITH_SCOPE`及论文11保持不变。结项记录已说明的动态台账快照差异，以及继承层中论文08的既有外部编辑例外，仍按原范围处理；本次归档不重写这些历史摘要，也不声称消除了所有既有历史例外。

目录因Windows占用曾拒绝重命名；用户关闭相关标签后，原地重命名成功。未通过复制科学材料绕过占用。
