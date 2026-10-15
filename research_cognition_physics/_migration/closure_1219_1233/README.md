# 1219—1233：归档收尾与广义相对论研究准备

2026-10-15。目录已从`archive_1219_`实际重命名为[archive_1219_1233](../../archive_1219_1233/README.md)。本次是目录维护和后继准备，没有新增科学轮次或公理。

## 归档内容

- 15轮编号报告、231个目录文件完整保留。
- 35份Python、102份JSON及原目标文本保持字节不变，验收哈希没有改写。
- 8份Markdown进行了可逆的链接或当前导航调整，涵盖阶段README、索引、根目录4份导航及09论文链接等。
- [迁移清单](manifest.json)共236项，登记原路径、现路径、前后SHA256及逆向编辑；133份既有支持文件保持不变。

文件只存在于当前目标目录；不建立旧目录外壳、重复资料目录或ZIP。[迁移程序](../../../scripts/close_research_1219_1233.py)先验证路径、库存和并发修改，再执行重命名。旧科学文档中的原命令与原回执路径属于历史记录，不批量重写其结果。

## 当前复算入口

在项目根目录通过[归档入口](../../../scripts/run_research_thermal_geometry_closed.py)复算：

```powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 scripts/run_research_thermal_geometry_closed.py --verify-closure
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 scripts/run_research_thermal_geometry_closed.py --verify-navigation
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 scripts/run_research_thermal_geometry_closed.py --script research_cognition_physics/archive_1219_1233/1233/check.py
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 scripts/run_research_thermal_geometry_closed.py --script research_cognition_physics/archive_1219_1233/1233/independent_check.py
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 scripts/run_research_thermal_geometry_closed.py --script research_cognition_physics/archive_1219_1233/_shared/check_stage_completion.py
```

入口同时接受原目录名与完整目录名。读取当前唯一文件，核验当前字节，在内存中逆转本次文档编辑并恢复原逻辑路径，再核对原证据；不把保存的哈希当作文件内容。归档复算只读，不允许改写原验收、结果或历史导航。直接运行旧导航更新器或旧带写入参数的脚本不再是维护入口。

## 检查记录及范围

最终检查见[verification.json](verification.json)、[复算记录](runtime_checks.json)、[链接检查](link_checks.json)与[独立审阅](independent_review.md)。当前归档核验和历史阶段科学验收分别登记，不互相替代。

此前08论文已有一次独立的外部格式编辑，旧1161—1218全量冻结检查因此保留`EXTERNAL_EDIT`状态；更早全量检查也有已记录的猜想文档变更。两者不由本次迁移造成。本次保留这些文件、旧哈希和例外证据；不把当前归档通过说成所有历史版本全量一致。09论文仅更新引用，原科学文本可按迁移清单逐字恢复。

## 后继研究

下一阶段可以直接进入**广义相对论的共同几何响应与动力学重建**，从1234起。首要问题是完整热流和独立几何读口之间，何种认知要求能够选择共同几何响应，随后连接成熟引力重建，而非重复一般热学或持续修补某台设备。

完整拟定范围、旧成果复用、首项问题和验收要求见[下一阶段准备](next_stage_plan.md)。本次仅准备，未开启新应用目标。1219—1233阶段原目标已完成；广义相对论目标仍待后续研究。
