# 1159—1160完成阶段的可逆归档

本迁移把`archive_1159_`实际重命名为`archive_1159_1160`。目录深度不变，不保留旧壳目录、不复制科学资料、不读取ZIP。新阶段从独立的`archive_1161_`开始。

## 冻结边界

[manifest.json](manifest.json)逐项记录41个实移文件及5个原地外部文档，共46条。8份Markdown的当前引用或导航说明使用可逆文本差分；16份Python／JSON文件保持字节不变。80份既有运行入口与迁移证据的SHA保持。基线来自已完成目标的真实回执及其资产清单，保留已有未提交研究成果，不依赖Git提交状态。

历史目标原文、历史复算命令和JSON路径仍保留其原记录；不能为消除旧路径字符串而刷新验收SHA。07论文的当前链接指向新目录，历史核查时恢复原链接字节。

## 当前只读复算入口

在项目根目录，使用既有Python运行[只读入口](../../../scripts/run_research_gr_active.py)，即`scripts/run_research_gr_active.py`，其参数如下：

```text
--verify-closure
--verify-layout
--verify-navigation
--script research_cognition_physics/archive_1159_1160/1159/check.py
--script research_cognition_physics/archive_1159_1160/1160/check.py
--script research_cognition_physics/archive_1159_1160/_shared/verify_stage.py
```

`--script`也接受历史`archive_1159_`路径。每份旧脚本得到各自原逻辑`__file__`、原文档字节和原文件枚举；后继1161脚本读取实际当前文件。路径代理只在受控脚本模块内生效，不全局替换pathlib。研究目录在复算期间只读。

旧阶段资产校验仍与原`stage_assets_check.json`比较：26份Markdown、246条链接和38个绑定资产。这里没有通过修改旧回执来适配改名。

## 1161活动导航

新导航的单份正文为[archive_1161_/_shared/navigation_notice.md](../../archive_1161_1200/_shared/navigation_notice.md)，正文链接以`research_cognition_physics`为基准。[新维护器](../../archive_1161_1200/_shared/update_navigation.py)的`--check`只读预览，`--apply`替换四根导航顶部唯一的新阶段栏目。

活动差分另存于`_migration/gr_active_1161/layer.json`，可逆还原本次归档边界，再连接冻结的`local_dynamics_active_1159`及此前各层。新notice可按实际研究进度维护；不得直接运行已归档1159维护器或手改冻结manifest/layer。新导航维护器需要直接执行，不能通过只读复算入口执行`--apply`。

## 核验

迁移及新导航的实际检查记录见[verification.json](verification.json)。完整历史层只在本次迁移实现稳定后核验；后继纯notice更新核验当前外层即可，旧层证据保留原适用范围。本文和迁移脚本属于行政资料，不计科学轮次、公理或经验校准。
