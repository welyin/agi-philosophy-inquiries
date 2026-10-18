# 1313起当前导航与历史复算

本层只记录研究根目录README、research_direction、RESEARCH_STATE及ROADMAP在原标题后的一处新增导航。历史正文和1308—1312归档保持；导航修改可在内存精确撤销。

- 当前阶段：[archive_1313_](../../archive_1313_/README.md)。
- 前一归档：[1308—1312迁移](../closure_1308_1312/README.md)。
- 本层清单：`layer.json`，由新阶段`_shared/publish_start.py`首次发布。
- 当前只读复算入口：仓库`scripts/run_research_interaction_generation_active.py`。

`--verify-navigation`只核当前四份导航及继承支持；`--verify-closure`核1308—1312归档；`--verify-conditional-navigation`核上一活动层。`--script`可指定历史旧路径或当前路径。新阶段科学脚本直接使用当前文件。复算复用现有Python，并使用`-B -X utf8`。

后续更新四份根导航时，先用`relocation_original_bytes`撤销本层插入，保持其原边界不变，再替换这一处插入内容并更新本层现哈希。应在写入前核当前文件和清单未被并发修改。不要重复运行首次发布器，不要直接改旧清单、冻结代码或保存结果来适配新目录。

历史层全面验证与本次局部检查范围不同；本层不声明已消除既存的全部历史布局差异。目标启动不计科学轮次。
