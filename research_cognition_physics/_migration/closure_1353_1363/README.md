# 第1353—1363轮目录收尾

第17篇阶段论文的11轮记录现集中在[归档目录](../../archive_1353_1363/README.md)。[逐轮验收](../../archive_1353_1363/_shared/stage_acceptance.json)和全部科研代码与结果保留原字节。第17篇论文、研究目标及第1478轮准入的改动仅涉及目录名称。档案说明补入了第18篇阶段论文的入口。

[迁移清单](manifest.json)登记原目录119个文件的 SHA-256，以及6个确实改动文件的新旧摘要和恢复方法。三个只改路径的文档可直接逆替换，其余文件以字节位置登记反向差异。原来的11份逐轮验收清单、阶段验收及启动记录都保留历史路径；[读取与核验程序](../../../scripts/research_sm_1353_closure.py)将这些路径映射到现有文件并验证原摘要。

从项目根目录，使用装有 NumPy 的既有 Python 运行以下命令。第一条复核归档和外部改动的原始字节；第二条独立重跑11轮的23项计算；第三条还检查与前后相邻阶段的导航哈希。

```powershell
python -B scripts/research_sm_1353_closure.py
python -B research_cognition_physics/archive_1353_1363/_shared/verify_stage.py --replay
python -B scripts/research_navigation_1353.py
```
