# 1044—1062轮归档与当前路径复算

2026-10-09。整理既有19轮，科学、经验和认知公理增量均为0；最新正式1062，累计3838。研究目标仍暂停，完整ROADMAP没有结项。

## 实际目录

|原目录|当前目录|科学状态|
|---|---|---|
|archive_1044_|[archive_1044_1045](../../archive_1044_1045/README.md)|1044—1045的原有限目标已结项|
|archive_1046_|[archive_1046_1062](../../archive_1046_1062/README.md)|1046—1062进展快照；整体目标未完成|

两个目录共595份既有文件全部迁移；正式报告仍在阶段根目录，证明、代码、结果及审阅仍在各轮编号目录。跨轮`_shared`和准入`_admission`不拆散、不复制。第二阶段README另给三组阅读次序，保留跨组依赖。

氢谱共同匹配保留作者稿及暂停审阅检查点，**尚未独立终签**；目录整理不将其计为已验收采用或1063。

## 历史内容与当前导航

[迁移清单](manifest.json)登记604项文件映射与前后版本，包括595个移动文件和9份外部引用文档。科学Python、JSON结果及冻结验收收据保持原字节；37份Markdown的路径或导航变动登记为可逆文本补丁。

两阶段README和索引改成单一当前入口，移除“最新1060”等过时导航；研究根目录同步归档状态。旧导航文字可由清单中的差分在内存中恢复并核对原哈希，不另外复制导航快照、整套目录或ZIP。项目根`README.md`和既有四份路径运行工具保持原字节。

1044目录直接重命名完成。Windows拒绝整体重命名1046目录后，逐个直接子项执行同盘移动，再删除已核为空的源目录；没有复制或遗漏科学文件，没有保留旧目录壳。此操作只改变物理位置。

## 复算入口

在项目根目录，使用已有Python运行：

```powershell
python -B -X utf8 scripts/organize_research_1044_1062.py --verify
python -B -X utf8 scripts/run_research_recent.py --verify-layout
python -B -X utf8 scripts/run_research_recent.py --script archive_1044_1045/_shared/verify_stage.py
python -B -X utf8 scripts/run_research_recent.py --script archive_1046_1062/_shared/verify_overall_evidence_after1062.py
python -B -X utf8 scripts/run_research_recent.py --script archive_1046_1062/1062/verify_round1062.py
```

若系统`python`指向不可用别名，请沿用已有解释器：
`C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe`。

[新只读入口](../../../scripts/run_research_recent.py)复用旧解析器并增加本次映射层。它先读取当前真实文件、核对当前哈希，再逆转已登记的Markdown改动并核旧哈希；不以清单里的摘要代替真实读文件。遇到更早的1009或764归档，会依次应用各自映射，历史子进程也通过新入口运行。科学运行禁止写入证据文件。

历史冻结脚本直接运行可能仍使用旧目录，应使用上述入口。旧脚本记录的是当时结果和假设；当前目录检查通过不等于重审全部数学或取得新物理证据，也不把原有环境相关浮点限制、旧失败或未终审材料改写为通过。

## 核验记录

本次核文件逐一对应、原字节可恢复、现行链接、编号连续性及代表性冻结入口；[总核验记录](checks.json)与[历史入口执行记录](replay_checks.json)保存实际结果。清单和科学结果均不覆盖。需要以后再改冻结正文或目录时，应新增明确维护层，不能静默跳过原哈希。

历史入口同时保留各时期的目录名称和相对路径写法；新运行器修正了将1044之后的旧阶段引用过早改名、将词法`..`过早折叠的兼容问题，未修改原收据。早期764目录层的三个根导航在本次整理前已随研究演化，其旧快照差异单列，不声称重新验证全部更早导航。数值线程环境沿用调用者设置。

[研究总目录](../../README.md) · [研究状态](../../RESEARCH_STATE.md) · [原路线图](../../ROADMAP.md)
