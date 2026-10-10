# 1086—1105阶段结项迁移记录

本次将 `research_cognition_physics/archive_1086_` 实际更名为 `archive_1086_1105`，不保留旧目录壳、副本树或 ZIP。阶段结项是公理化、论文综合与文件归档，不是完整物理狭义相对论或更强目标完成；没有新研究轮次、科学实验或校准增量。

## 内容与边界

- `preparation_snapshot.json`：迁移前原235份阶段资产、7份相关外部Markdown及冻结支持文件的284项哈希基线，不含正文副本。
- `manifest.json`：实际迁移的243份阶段资产及7份外部文档，共250项；其中11份Markdown有可逆修改。科学PY/JSON不改字节，新增采用与审阅文件单独计入。
- `close_stage.py`：一次性迁移驱动及事后核验。源、目标必须是项目研究目录下的指定兄弟目录，迁移前再次校验全部基线。
- `preflight_checks.json`：迁移前新读取边界与所有旧层兼容性的检查记录；临时预检清单已删除。
- `replay_checks.json`：迁移、全部旧层检查和1097、1100、1101、1102、1105默认复算的实际输出。

唯一允许迁移后补入的文件是阶段目录 `_shared/stage_paper_closure/publication_checks.json`，由出版审核最后写入；它不是旧科学资产。校验器检查原盘点完整且没有其它额外文件，并单独报告这份收据的当前哈希。

## 复算入口

使用新的最外层入口，路径可以采用实际的新阶段目录：

```powershell
python -B -X utf8 scripts/run_research_sr_archive.py --verify-layout
python -B -X utf8 scripts/run_research_sr_archive.py --script research_cognition_physics/archive_1086_1105/1105/check.py
python -B -X utf8 research_cognition_physics/_migration/sr_closure_1086_1105/close_stage.py --verify
```

本次实际复算复用了既有 Python 运行时。新入口先验证真实文件，内存中撤销本次导航和路径修改，再依次交给冻结的active、公式、历史导航及旧归档层。1086—1105源文件获得迁移前逻辑路径与内容视图；更早的实验继续获得各自既有历史视图。旧程序与旧清单中的历史目录名保留作为逻辑键，不依赖磁盘上的旧目录。

冻结的旧 `run_research_active.py` 和active清单没有改写。请在迁移后的工作区使用上述新入口，不直接运行只认识旧物理目录的历史包装器。复算禁止研究资产写入和ZIP读取，不运行发布脚本。

## 核验范围

迁移后144份PY/JSON与各自迁移前字节一致；其中142份属于原235项阶段基线，另2份为本次新采用记录。全层校验还验证旧6949份PY/JSON、8340个公式体，以及旧布局实际存在的7146份资产；2819份按照原清单规则标为可选本地资产，未将缺省项虚报为已复算。代表轮次只重新比较原结果，不增加科学证据或把有限示例升级为实际物理实现。

阶段正文的旧“当前”“候选”“下一步”作为历史保留。六份导航顶端的新结项栏优先说明采用范围、论文入口、后继责任及目标仍暂停的状态。
