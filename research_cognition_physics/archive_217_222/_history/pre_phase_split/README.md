# 第1—222轮研究档案（已结项）

[研究总目录](../README.md) · [阶段论文](../可组合认知结构与复量子状态空间_阶段论文.md) · [逐轮结论与范围](ROUND_INDEX.md) · [配套文件索引](文件索引.md)

## 阶段范围

本阶段建立并审查从认知操作要求到复量子状态空间的条件性重建链；阶段论文保留明确的工作定义、数学前提和适用范围。不把工作定义的采用等同于证明所有认知系统或所有物理现象都满足这些前提。

正式报告为本目录的 `research_note_01.md` 至 `research_note_222.md`，保留原有编号写法。代码、结果、检查及补充材料按轮次放在 `1/` 至 `222/` 中。共用文件按首次可确认的轮次归属，后续轮次通过链接引用。

- [第219轮修订稿](research_note_219.md)及[原稿和清单](219/history/round219_v1/archive_manifest.json)均保留。
- [第222轮报告](research_note_222.md)接至[223—230轮有限维量子理论阶段](../archive_223_230/README.md)。
- [原研究导航](_history/indexes/README.md)和[原归档导航](222/archive_closure/README_before_round_layout.md)按历史语境阅读；当前方向以[研究状态](../RESEARCH_STATE.md)为准。

## 目录约定

|位置|内容|
|---|---|
|`research_note_编号.md`|正式编号研究报告|
|`1/` … `222/`|对应轮次的代码、结果和核验；联合检查放在文件名中最早一轮|
|`219/history/`|219轮修订前的冻结版本|
|`222/archive_closure/`|阶段论文维护、原索引、原清单、历史归档核验和迁移备份|
|`_shared/`|没有单一轮次归属的背景讨论及共用依赖说明|
|`_history/`|原导航和保留的运行缓存；缓存不计科研证据|

原 `research_process/` 的材料已实际迁移。独立的物理构造、信息几何路线仍在项目根目录，不恢复此前已删除的两套归档副本；旧迁移备份仅作为已有历史原件保存。

## 复算与原件保护

代码、JSON结果及历史收据保持原字节；阅读版Markdown仅调整链接目标。原件与路径变化见[本次迁移说明](../_migration/layout_001_230_20261004/README.md)。旧脚本使用固定目录与跨轮导入，统一通过隔离复算入口恢复原布局：

```powershell
python -B -X utf8 research_cognition_physics/_migration/layout_001_230_20261004/replay_early.py --stage 1
```

命令从项目根目录运行，使用既有Python/NumPy环境；可加 `--module local_calibration_contract` 只检查指定模块。[本次复算记录](../_migration/layout_001_230_20261004/replay_checks.json)与[原1915项回归记录](222/archive_closure/stage1_test_results.json)分别保留，不把目录整理计为新研究轮次。
