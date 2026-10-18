# 1288—1289与1290—1307归档

2026-10-18，按用户要求完成两阶段目录收尾：

- [archive_1288_1289](../../archive_1288_1289/README.md)：2轮，对应论文12。
- [archive_1290_1307](../../archive_1290_1307/README.md)：18轮，对应论文13。

原377件文件全部原地迁入，193件Python／JSON逐字节保留。27件Markdown只改迁移路径或归档导航，原字节可从唯一[迁移清单](manifest.json)反向恢复；另有8件外部Markdown纳入本次记录。旧目录已移除，没有复制档案、ZIP、影子目录或符号链接。原先已有的缓存随目录保留为不参与复算的普通库存，本次不生成缓存。

本次不增加科学轮次或公理，不扩大既有阶段验收范围。报告内“待审”等字样保留当时状态，最新结论以各轮正式验收及阶段README为准。

## 当前复算入口

从项目根目录使用已有Python环境运行：

```text
python -B -X utf8 scripts/run_research_structure_selection_closed.py --verify-closure
python -B -X utf8 scripts/run_research_structure_selection_closed.py --verify-navigation
python -B -X utf8 scripts/run_research_structure_selection_closed.py --script research_cognition_physics/archive_1288_1289/1288/independent_check.py
python -B -X utf8 scripts/run_research_structure_selection_closed.py --script research_cognition_physics/archive_1290_1307/1307/check.py -- --check
python -B -X utf8 scripts/run_research_structure_selection_closed.py --script research_cognition_physics/archive_1290_1307/1307/independent_check.py -- --check
```

入口同时接受旧逻辑路径和新物理路径，在内存撤销有记录的Markdown修改后核对旧摘要，不忽略哈希，不改旧科学代码。历史写结果模式被禁止。

若后续1308阶段已添加新的根导航，使用`run_research_conditional_selection_active.py`作为最外层入口；上述参数和脚本路径不变，旧入口仍为冻结依赖。

## 验证范围

目录清单、385件条目的原始／当前摘要、149件冻结支持，以及1307作者／独立、1288独立和1290—1307联合复核已通过。本次检查由主线执行；代理完成入口初稿后遇模型容量错误，不将尚未完成的代理验收声称为独立通过。

全部旧层的`--verify-layout`仍会报告论文08在更早阶段已登记的外部编辑差异，见[前阶段审计](../closure_1281_1287/independent_audit.md)。本次没有改动论文08，也不把这项既有差异说成本次迁移造成或已经修复。

本次还核定235项原轮次验收绑定和153处有效迁移链接。论文12的历史`phase_paper_draft.md`原本采用拟发布于研究根目录的路径写法，19处链接在草稿存放位置原已不可达；正式论文链接有效。这些原有草稿位置问题在核验中单列，不算新增断链，也没有重写冻结草稿的历史内容。

[迁移前库存](before.json) · [迁移程序](migrate.py) · [核验记录](verification.json)
