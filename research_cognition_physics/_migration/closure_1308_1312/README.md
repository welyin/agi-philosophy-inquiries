# 1308—1312归档迁移

2026-10-18。已将研究目录从 archive_1308_ 实际重命名为 [archive_1308_1312](../../archive_1308_1312/README.md)。本次只整理目录，不新增科学结论、研究轮次或公理。

## 保留与核验

- 84份原文件全部迁移；42份 Python／JSON 逐字节不变。
- 89条迁移记录包括84份阶段文件和5份外部Markdown。仅6份Markdown有可逆修改：阶段README、根目录四份导航和论文14。
- 5轮的48项验收绑定和阶段验收绑定按原始字节恢复核验；153份继承运行支持文件保持不变。
- 代表复算覆盖1312作者及独立程序、1311跨轮导入以及旧1307程序。旧路径和新路径均可进入同一实际文件。
- 未使用ZIP、目录复制或隐藏历史目录；新阶段1313文件不属于本清单，由主线维护。

[完整清单](manifest.json) · [原文件盘点](before.json) · [核验记录](verification.json) · [可复核核验程序](verify.py)

## 复算入口

使用现有Python和 scripts/run_research_conditional_selection_closed.py，以 --script 指定本阶段旧路径或现路径，再传 -- --check。该只读入口在内存撤销归档路径修改后执行冻结代码；原始科学文件不改写。

新阶段导航接入后，使用 scripts/run_research_interaction_generation_active.py 作为最外层入口；其 --verify-closure 验证本次归档，--verify-conditional-navigation 验证1308导航。迁移核验程序 verify.py 会在1313导航层存在时自动通过该层读取原字节，默认只读；--write-new 仅用于首次保存核验记录。

## 可逆边界

manifest.json 每条记录含原路径、现路径、原／现SHA-256及精确Markdown替换区间。闭包类 ClosedConditionalSelectionLayout 的 before_conditional_closure() 恢复迁移前字节；before_conditional_selection_navigation() 进一步撤销原1308导航。后继导航只允许在根四文件标题后添加一个可逆前缀，并以本清单的现哈希作为其原边界。

本次仅验收当前迁移和已列代表复算。既往论文08的历史全布局差异仍按原审计记录保留，不宣称所有旧层的全布局检查全部通过，也未借本次迁移修改该论文。
