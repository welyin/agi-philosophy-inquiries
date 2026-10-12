# 1120—1158合并迁移与历史复算

2026-10-12。按用户要求，将1120—1156和1157—1158实际合并至[archive_1120_1158](../../archive_1120_1158/README.md)。旧目录`archive_1120_`和`archive_1157_1158`不再保留。

## 实际布局

- 39份编号报告放在合并目录根部；1120至1158的39个编号目录分别保留各轮代码、证明、结果与审阅。
- 两段原共享材料、README和索引分别移至`_shared/1120_1156/`、`_shared/1157_1158/`，不覆盖同名文件。
- 新根README和统一轮次索引提供当前入口。两段原目标、结项范围及科学结论不变。

[单份清单](manifest.json)记录597项实际移动及6份外部Markdown的引用更新；247份科学PY／JSON保持原字节。Markdown引用和导航编辑可在内存中逆向恢复到迁移前的精确字节，以便核对旧审阅哈希。旧迁移清单、旧运行时和活动导航基线不改写。

一份原有且已被Git忽略的Python字节码缓存随目录原样移入，标为可选的本地附属文件：存在时核对原字节，缺失时也不阻断复算；不加载、不作为科研依据。没有创建ZIP、旧目录壳、重复研究文件或导航副本。

## 当前只读入口

从项目根目录使用已有Python与NumPy运行时：

```powershell
python -B -X utf8 scripts/run_research_combined.py --verify-merge
python -B -X utf8 scripts/run_research_combined.py --verify-layout
python -B -X utf8 scripts/run_research_combined.py --script research_cognition_physics/archive_1120_1158/1156/check.py
python -B -X utf8 scripts/run_research_combined.py --script research_cognition_physics/archive_1120_1158/1157/check.py
python -B -X utf8 scripts/run_research_combined.py --script research_cognition_physics/archive_1120_1158/1158/check.py
python -B -X utf8 scripts/run_research_combined.py --script research_cognition_physics/archive_1120_1158/_shared/1120_1156/stage_closure/verify_assets.py
python -B -X utf8 scripts/run_research_combined.py --script research_cognition_physics/archive_1120_1158/_shared/1157_1158/verify_stage.py -- --require-closure
```

新入口同时接受原逻辑脚本路径和当前物理路径。每个冻结模块保有原有`__file__`及对应历史视图，再解析到唯一现存文件，因此原父目录层数、哈希和跨阶段导入可以复用。禁止复算时写研究资产；代码没有被改成忽略哈希或读取替代结果。

旧文档中的直接执行方式、旧运行时入口属于当时记录；合并后需要验证这些冻结代码时，请使用上述最外层入口。新层复用此前的历史恢复层，不重新生成早期档案。

## 核验

清单、文件唯一性、字节恢复、原目录消失、本地引用、两段结项校验、1156—1158默认只读复算及既有全布局检查均通过。另实际测试了移开可选缓存后核验仍通过，并恢复原缓存字节。

[执行记录](verification.json)区分主线资产核验、运行时代理执行及最终缓存测试；[独立审计](independent_review.md)核对迁移、原审阅哈希与实际引用。本次不新增科研轮次或数值实验组。
