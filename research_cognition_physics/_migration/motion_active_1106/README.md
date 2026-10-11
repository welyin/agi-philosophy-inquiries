# 1106起活动导航层

本层维护研究根目录README、research_direction、RESEARCH_STATE、ROADMAP的当前导航，以及当前洛伦兹阶段论文的认知公理体系措辞。1086—1105档案、旧运行器、旧清单和原论文发布回执保持原字节；历史复算通过同一层的可逆差分恢复发布时论文。

`layer.json`目前保存5项可逆Markdown差分及前后哈希，不存整份文档副本；`start_navigation.py`为本次一次性插入脚本，不可重复执行。`checks.json`保留启动时的旧层兼容性校验；本次措辞修订及当前校验见[修订记录](../../archive_1106_/_shared/axiom_wording_revision/README.md)。导航只是启动记录，尚不计为新正式研究报告或科学校准。

迁移后的历史复算使用：

```powershell
python -B -X utf8 scripts/run_research_motion_active.py --verify-layout
python -B -X utf8 scripts/run_research_motion_active.py --script research_cognition_physics/archive_1086_1105/1105/check.py
```

新入口先验证真实活动文档并在内存中恢复结项时字节，再调用冻结的SR归档、active、公式与更早归档验证逻辑。没有改写旧哈希，也没有复制目录或读取ZIP。直接使用旧最外层入口核验新导航会报告版本差异，应以上述入口为准。

该入口用于历史实验复算，未替新阶段代码定义科学权限或运行模型。以后更新列入本层的活动Markdown时，应继续登记与当前层起始字节之间的可逆差分，不修改任何更早层；不要为每次文字更新再复制整份历史文件。
