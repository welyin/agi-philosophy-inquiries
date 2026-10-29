# 1313—1322 归档迁移

2026-10-19。已将 archive_1313_ 原地重命名为 [archive_1313_1322](../../archive_1313_1322/README.md)。本次只整理已结项阶段，不新增研究轮次、公理或科学结论。

## 保留与核验

- 209 份原文件全部迁移；119 份科学 Python／JSON 逐字节不变。
- 214 条迁移记录包含全部阶段文件和根四导航、论文 15。仅 10 份 Markdown 修改路径或归档导航，原字节可精确恢复。
- 十轮的 91 条验收绑定和阶段的 105 条绑定均核对原始字节；继承支持文件保持冻结。
- 代表复算：1322 作者程序通过 155 项、独立程序通过 216 项，1321 作者程序通过 176 项；均只读，旧逻辑路径及新物理路径都可使用。
- 未用 ZIP、目录副本或隐藏旧目录。新阶段 archive_1323_ 及其目标不属于本迁移清单。

[迁移清单](manifest.json) · [原始盘点](before.json) · [验收记录](verification.json) · [只读核验程序](verify.py)

## 历史复算入口

使用现有 Python 运行 scripts/run_research_interaction_generation_closed.py，以 --script 指定本阶段旧路径或新路径，再传 -- --check。入口在内存恢复冻结文件的历史视图，不修改科学脚本或结果。

新阶段接入后，以 scripts/run_research_interaction_rules_active.py 为最外层入口。--verify-closure 检查本次归档，--verify-interaction-navigation 检查恢复后的 1313 历史导航；--verify-navigation 检查最外层当前导航。verify.py 在新导航层存在时自动穿过该层，默认只读；--write-new 仅用于首次保存验收记录。

## 可逆与接入边界

manifest.json 保存每项原路径、现路径、原／现 SHA-256 与 Markdown 精确替换区间。ClosedInteractionGenerationLayout.before_interaction_closure() 恢复迁移前字节；before_interaction_navigation() 进一步撤销旧阶段导航前缀。后继只在根四文件标题之后添加可逆前缀，以本次清单的 current_sha256 为其 original_sha256；不得改写旧导航层或科学证据。

本次验收限于该迁移、完整冻结绑定与已列代表复算。既往论文 08 的历史全布局差异仍按原审计保留，不在本次修订，也不据此宣称所有历史全布局检查均通过。
