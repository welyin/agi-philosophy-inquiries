# Markdown公式格式统一与历史复算入口

本次是用户明确授权的格式维护，不是新科学轮次。全项目Markdown统一使用独立行的块级双美元定界符和单美元行内定界符。既有公式正文、非数学文字、代码、旧Python/JSON和历史科学结论保持原字节或可精确逆恢复；没有原文备份、ZIP或第二份历史目录。

## 本次结果

- 扫描3488份Markdown，修改654份。
- 转换3428个旧块级公式、4027个旧行内公式。
- 另对885个原有双美元公式补齐独占行或段落边界。
- 8340份涉及的公式正文、全部差分逆哈希及重复执行幂等核验通过。
- 6949份既有Python/JSON逐字节未改；7份纯导航快照被排除并保留。
- 代码围栏、缩进代码、行内/跨行代码、HTML注释、转义以及TeX换行间距被保护；没有替换裸方括号。
- 后扫描无剩余可转换公式、无未配对公式异常。唯一正文中的双美元符号说明保留为普通文字，并在清单记录。

唯一可逆格式差分是[manifest.json](manifest.json)。其中只保存定界符/必要边界差分、公式体摘要和输入/输出哈希，不保存全文副本。[格式核验](format_verification.json)及[完整集成回执](integration_checks.json)分别记录结果。

## 现在使用的只读复算入口

旧入口源码不改。经过公式维护的文档由新最外层入口提供原始读取视图：

```powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 scripts/run_research_formatted.py --verify-layout

& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 scripts/run_research_formatted.py --script research_cognition_physics/archive_1063_1085/1085/verify.py -- --recompute
```

读取顺序为：实际文档 → 撤销本次公式差分 → 撤销1086新阶段导航差分 → 撤销1063—1085目录迁移 → 继续既有导航/更早布局层 → 原冻结字节。每一层先验证当前哈希，再验证恢复后的原哈希；子进程同样通过新最外层入口。

实际文档仍显示统一后的公式。复算只是在内存中读取历史视图，不回写文档，也不改历史哈希“迎合”新格式。

## 自身维护检查

```powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 scripts/normalize_markdown_math.py verify
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 scripts/normalize_markdown_math.py plan
```

本次后计划为零修改。应用模式拒绝覆盖已存在的格式层；以后若进行其他维护，应建立有明确授权的新差分层，不覆盖这份回执。

## 验证范围

最新集成核验通过422项空间阶段迁移入口、604项近期迁移入口以及7146项更早源布局检查。1085的14份资产、5份冻结来源、21条本地链接，以及作者和独立模型的实际重算全部通过。原有三份更早导航时点差异仍由其原迁移记录解释，没有借本次格式维护覆盖。

转换器使用保守扫描而非完整Markdown AST；保护不明确区域、不把货币符号或代码中的美元认作公式。不进行图像检验，也不承诺每一种Markdown阅读器都启用了数学渲染。全部数学正文等价性由实际差分、保护区检查、公式体摘要、逆字节哈希和幂等检查承担。

## 审计时点与后续正常研究

清单中的全量Markdown输入库存用于证明**本次维护时点**只改了计划内容。它不是把所有活动文档永久冻结：根目录当前导航及1086阶段非候选文档仍可按正常研究更新；未来全量核验若发现这类变更，应按其自己的导航/阶段差分解释，不能称为历史科学资产损坏。读取冻结历史时，真正必须严格恢复的是被转换条目的原始字节及旧科学资产；后续若升级维护入口，应沿外层追加明确差分，不重写本次manifest、旧Python/JSON或历史哈希。

旧1062默认整体证据入口也已通过新增最外层入口实跑，证明格式逆恢复穿透更深历史层；完整命令和回执追加在integration_checks.json。
