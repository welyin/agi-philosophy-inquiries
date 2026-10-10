# 迁移与当前目录复算


## 2026-10-09：1044—1062归档

[19轮归档、路径核验与只读复算](layout_1044_1062_20261009/README.md)。新入口`scripts/run_research_recent.py`承接1009关闭层及本次两阶段路径；不恢复旧目录、不复制科学材料。

## 2026-10-08：764—1008再次分期

旧迁移元数据已集中到本目录，科学文件由当前物理位置读取。[本次10阶段划分、完整核验及复算说明](layout_764_1008_20261008/README.md)。新增一层精确路径与文本变更映射，保留旧清单和旧收据，不创建旧目录或ZIP。776—1008也通过统一入口运行，例如：

```powershell
python -B -X utf8 scripts/run_research_current.py --script archive_990_1008/1008/overall_completion_audit.py
```

001—775轮正式报告和配套材料已实际放入各阶段目录。目录位置由当前迁移清单记录，科学代码与结果保持原字节。

## 当前复算入口（2026-10-05修正）

旧入口曾通过解压整包ZIP恢复旧目录，这项兼容绕路已移除。现在入口根据迁移清单直接加载各阶段目录中的代码和数据，不解压ZIP、不创建旧布局、不复制整套研究。运行器禁止读取ZIP及`navigation_before_*`，并禁止写入研究证据。

```powershell
python -B -X utf8 scripts/run_research_current.py --stage 1
python -B -X utf8 scripts/run_research_current.py --stage 2
python -B -X utf8 scripts/run_research_current.py --script joint_wick_source_anomaly.py
python -B -X utf8 research_cognition_physics/_migration/replay.py --suite
```

从项目根目录运行。早期的[兼容入口](layout_001_230_20261004/replay_early.py)也已转接当前目录运行器；`--script`既支持旧脚本名，也支持当前相对路径。跨轮次导入和数据引用由[统一路径解析器](../../scripts/research_layout.py)定位到物理文件，不要求逐个修改冻结算法。

## 如何核对原收据

[231—775迁移清单](manifest.json)和[001—230迁移清单](layout_001_230_20261004/manifest.json)记录原路径、当前路径、原哈希、当前哈希和Markdown链接替换。

运行器先检查当前文件的真实字节。科学Python与JSON直接使用当前文件。对于历史收据中的Markdown原哈希，只在内存中逆转已记录的链接替换，再计算并核对原哈希；不改正文、不替换磁盘上的新链接、不直接返回清单中的哈希充当验证。当前目录核验也直接用这些文件和替换记录，无需ZIP。

## 重复备份与历史发布核验

README、方向和状态只在正式位置维护一份，修改过程交给Git。`navigation_before_*`、整理前ZIP及被当前清单接替的大型清单副本是重复维护材料，已列入清理；今后不再按轮次复制整套导航。

旧`postcheck_round*.py`还要求“发布当时README逐字节相同”，这种检查针对当时发布状态，不适用于已经重新编排的现行目录。旧收据保留为历史记录；当前入口分别核验科学复算和当前布局，不把删除导航快照后的旧发布检查宣称为通过。

个别全链来源检查会读取第660轮原始PDF；它是外部研究资料，可按[来源收据](../archive_653_701/660/drafts/source_receipt.json)获取并核对哈希。常规科学复算套件不需要该PDF。

## 已有核验与范围

修正前的[首次整链复算](replay_checks.json)保留425字段缺失和523工作区外Spark引用两项旧限制，不改写其失败记录。原单元测试收据也保留。当前验证结果见[当前目录检查](current_layout_20261005_checks.json)。未重跑全部545轮，也未重审全部数学证明。

[研究总目录](../README.md) · [阶段成果](../_shared/notes/阶段成果总览.md) · [版本管理说明](../../scripts/GIT_STORAGE.md)
