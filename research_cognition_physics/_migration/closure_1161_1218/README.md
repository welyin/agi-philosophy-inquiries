# 1161—1200与1201—1218归档交接

本次按用户要求，将两个活动目录实际归档为：

|轮次|归档目录|研究状态|
|---|---|---|
|1161—1200|[archive_1161_1200](../../archive_1161_1200/README.md)|40轮、130组验收；GR研究记录保留，完整GR目标未完成。|
|1201—1218|[archive_1201_1218](../../archive_1201_1218/README.md)|18轮、36组验收；热力学四律在声明的操作及资源范围内结项。|

本次是目录和引用维护，不增加科学轮次，不改变公理、证明、模型输入或目标状态。下一研究轮从1219起；尚未创建新研究目录或启用新目标。

共实际移动984个文件，其中1161—1200目录717个、1201—1218目录267个。58份编号报告与58份验收齐全；507份Python／JSON保持字节不变。24份Markdown发生可逆的链接或管理说明调整，包括8份目录外文档。

## 保存方式

[manifest.json](manifest.json)登记每个文件的原路径、当前路径和前后SHA256，以及文档引用修改的逆向差分。文件只在当前目标目录中保存一份；不建立旧目录外壳、影子目录、符号链接、ZIP或第二份科学资料。

Python、JSON、目标原文及结果保留原字节。Markdown的当前链接、归档提示和08论文复算命令通过可逆修改更新；复算时先核当前字节，再恢复历史字节以匹配旧验收。原回执中的路径和哈希代表原时点，不通过重写它们制造“验收通过”。

## 当前只读复算入口

在项目根目录使用[统一入口](../../../scripts/run_research_thermodynamics_closed.py)：

```powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'scripts/run_research_thermodynamics_closed.py' --verify-closure
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'scripts/run_research_thermodynamics_closed.py' --verify-layout
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'scripts/run_research_thermodynamics_closed.py' --script 'research_cognition_physics/archive_1201_1218/1218/check.py'
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'scripts/run_research_thermodynamics_closed.py' --script 'research_cognition_physics/archive_1201_1218/_shared/after1218_common_prediction_check.py'
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'scripts/run_research_thermodynamics_closed.py' --script 'research_cognition_physics/archive_1201_1218/_shared/finalize_four_laws_stage.py' -- --verify
```

`--script`同时接受原目录名和完整轮次目录名。历史脚本按原逻辑路径与原文档视图运行，实际读取当前唯一文件。路径代理只用于受控研究模块；不全局替换`pathlib`。复算入口禁止向研究证据写入，不允许`--write`重建原结果。

旧`archive_1161_1200/_shared/update_navigation.py`是历史阶段的导航维护器。归档后仅通过上述入口执行`--check`核验原状态，不再直接执行`--apply`。下一阶段的导航应在当前归档边界上建立新的管理层，保留本次与更早迁移证据。

## 下一阶段的研究交接

继承[第08篇阶段论文](../../08_认知热学角色与热力学四律_阶段论文.md)、[逐项目标验收](../../archive_1201_1218/_shared/four_laws_goal_completion_audit.md)与[实际过程对象表](../../archive_1201_1218/_shared/current_process_geometry_interfaces.md)。后继方向是同一实际过程的热力学—几何连接，首先确认：

1. 实际钟、三维位置、信号和区域边界属于同一过程，不能把内部控制坐标直接当作物理空间。
2. 热参考、区域能流及观察者的演化参数如何接续；普通热平衡本身不自动给出Unruh关系。
3. 面积响应具有独立的操作或几何来源；不能用`A=S/s`定义出所要证明的面积律。
4. 在受控有效范围中，完整来源、反作用与几何动力学共同成立；现有四律不等于Einstein方程。

复用既有空间与跨尺度接口，特别是382—386、425和522—523；具体候选、权限、误差、资源增长及失败判据应在新目标中登记。本次不选定新的科学模型，也不把旧GR阶段的未完成项改写为已完成。

## 核验记录

迁移执行器为[close_research_1161_1218.py](../../../scripts/close_research_1161_1218.py)。最终迁移、路径、冻结哈希及代表性复算结果保存在[verification.json](verification.json)。管理资料不计科学轮次。

当前归档核验、原GR导航核验、1218作者与独立程序、四律共同预测和原结项核验均已通过。[链接检查](link_checks.json)覆盖477份文档、9,074处链接，无缺失或仍指向两旧目录的链接目标；[运行记录](runtime_checks.json)保存各项原输出。

**旧全仓检查的已知范围：** `--verify-layout`继续进入更早公式规范化层时，会报告`猜想/2029-9-29.md`与该旧层的冻结版本不同。逆向恢复本次迁移后，读到的文档哈希与本次迁移前完全一致，仍不同于那份更早的公式检查清单，证明差异早于本次归档。这里保留后来更新的猜想文档及旧验收，不刷新旧哈希。该历史全量检查没有通过；不能把上述当前归档与代表性复算的通过表述为全仓全部历史版本均一致。
