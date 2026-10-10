# 1063—1085阶段：实际目录归档

2026-10-10。按用户要求，将实际目录 `archive_1063_` 一次改名为 `archive_1063_1085`。416份文件全部原位移动，没有旧目录壳、复制目录、junction或ZIP；没有新增科学轮次。本阶段仍为22组科学校准，累计3860，最后正式报告1085为计0综合轮。

[逐文件迁移清单](manifest.json)覆盖416份移动文件及6份外部引用文档。17份Markdown仅修改阶段目录或复算入口路径，保存精确可逆字符补丁；80份Python和122份JSON保留原字节。旧科学收据、旧运行器、旧导航清单及上一阶段迁移回执均未重写。

## 只读复算

从项目根目录运行：

```powershell
python -B -X utf8 scripts/run_research_spatial_archive.py --verify-layout
python -B -X utf8 scripts/run_research_spatial_archive.py --script archive_1063_1085/1085/verify.py -- --recompute
python -B -X utf8 scripts/run_research_spatial_archive.py --script archive_1046_1062/_shared/verify_overall_evidence_after1062.py
```

第二条保留原1085入口，实际执行作者及独立模型的两套有限校准。原入口输出的 `goal_status=active` 是冻结历史字段；目标后来完成的依据仍是[实际结项回执](../../archive_1063_1085/_shared/spatial_goal_closure/goal_completion.json)，本次迁移不改历史字段。

旧逻辑路径亦可传入 `--script`。直接执行阶段历史脚本、或单独使用旧 `run_research_live.py`，仍可能按旧目录查依赖；本次新入口在进程内定位移动后的导航清单并保持原脚本的路径字典语义。无需恢复旧物理目录。

## 验证与层次

[实际检查结果](checks.json)登记哈希、链接、原验证器输出及运行器身份。

读取顺序是实际文件、撤销本次目录引用补丁、再进入原导航和旧归档层。每一步先校验当前字节与长度，再反演补丁并校验原哈希。1063—1085脚本使用本阶段迁移前视图；更早脚本仍使用它们各自的历史视图。没有将保存的哈希当作文件内容或直接当作校验通过。

新运行器的 `SpatialLayout.physical_bytes(path)` 是最外层读取点。以后若增加新的导航或数学格式维护，外层须先验证真实当前文件并精确反演到本次清单的字节，然后交回此层；不得重写本迁移清单、旧收据或科学程序。新增阶段的资料及维护层另行保存。

本次没有转换公式分隔符、重写论文观点或改变科学结论。更早764迁移记录中的三份已演化根导航差异仍单列保留，不声称本次恢复了那份无关历史快照。
