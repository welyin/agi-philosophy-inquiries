# 1234—1280目录收尾与复算

2026-10-17。对应[研究档案](../../archive_1234_1280/README.md)和[阶段论文10](../../10_认知操作公理与有效广义相对论_阶段论文.md)。

## 整理决策

保留一个 `archive_1234_1280`。47轮共同服务于同一结项成果，父模型、共享来源账和后续资格检验跨越多轮；拆分会割裂证据链。阶段README按1234—1244、1245—1252、1253—1270、1271—1278、1279—1280及收尾接续五段导航。这些只是阅读分组。

实际移动原1308个文件；报告在阶段根目录，代码、结果和各轮证书在数字子目录，跨轮材料在 `_shared` 和 `_admission`。不复制共享稿，不创建ZIP、旧目录副本或路径替身。

## 证据与路径

[manifest.json](manifest.json)记录实际文件的原路径、新路径、前后SHA256与可逆Markdown编辑。科学Python与JSON字节不变；文档只更改引用和归档导航。历史冻结校验通过只读入口恢复旧逻辑路径和修改前的Markdown字节，不改写摘要函数，也不替换科学结果。

历史报告中的命令保留原时点。复算时使用以下新入口；`--script`同时接受旧路径与新路径。不要运行旧出版脚本的写入选项。

```powershell
$python = 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
& $python -B -X utf8 scripts/run_research_gr_response_closed.py --verify-closure
& $python -B -X utf8 scripts/run_research_gr_response_closed.py --verify-navigation
& $python -B -X utf8 scripts/run_research_gr_response_closed.py --script research_cognition_physics/archive_1234_1280/1280/check.py
& $python -B -X utf8 scripts/run_research_gr_response_closed.py --script research_cognition_physics/archive_1234_1280/_admission/jacobson_effective_joint_after1280/independent_check.py
```

以上命令在项目根目录运行；从其它位置运行时给入口和脚本完整路径。运行器只读，禁止写研究文件。

`--verify-closure`只验本次1234—1280归档；`--verify-navigation`核恢复后的导航层。更早的paper08外部编辑例外仍按原边界保留；本次不重写它的验收，也不把本次归档通过称为全部历史问题已消除。

## 收尾验收

- [独立迁移审计](independent_migration_audit.json)：47份报告和47个轮次目录齐全；1308个迁入文件与清单一致，678个Python／JSON字节不变，13个Markdown编辑可逆；702个新路径链接已核对。
- [复算记录](runtime_verification.json)：本次闭包、恢复后的导航、1280独审、联合合同独审、1234嵌套历史复算，以及跨工作目录运行，6项均通过。
- [运行器代理审阅](agent_runtime_review.json)：另外覆盖1238嵌套读取器、1245子进程和跨轮根相对哈希，并确认旧paper08例外仍显式报出。

独立审阅进一步确认只读路径适配和写入保护；复算资格限现有冻结入口及已覆盖调用形式，不能据此保证所有旧脚本绕过入口直接运行。不会向已冻结阶段塞入导航备份。本次无新增科学轮次、公理、目标或实验。
