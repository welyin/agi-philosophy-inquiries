# 第二阶段档案：第223—230轮（已结项）

更新：2026-09-21。按用户指示从archive_223_重命名为archive_223_230；所有223—230轮科学证据保持原字节。

## 当前论文与历史版本

**[整合论文v1.1](../可组合认知结构与有限维量子理论_阶段论文.md)**已把229—230轮结果直接纳入定义、主定理、证明和有限实验范围。结项范围为有限维量子操作规则与连续可逆动力学形式的条件性重建。

[旧论文v1.0](paper_versions/finite_quantum_v1.0.md)、[原补充稿](STAGE2_ADDENDUM.md)、[旧档案索引](paper_versions/archive_readme_before_rename.md)保留原文。旧文件内的停止指令和路径按写作时语境解释；用户已经补充猜想，后续进入[第三阶段](../archive_231_/README.md)。

第229轮临时历史索引中的旧archive_223_链接也保留原文，当前目标为本索引。核验报告单独列出这一旧链接映射；它不属于当前导航链接。

## 逐轮证据

| 轮次 | 结果 | 范围 | 实现与检查 |
|---|---|---|---|
| [223](research_note_223.md) | 从态集仿射双射接入酉／转置酉分类；恒等连续分支为酉；连续齐次可逆时间群有−i[K,·]生成元 | Time为明示新增合同；具体K、钟尺、ℏ及全部CPTP仪器未选定 | [代码](reversible_dynamics_bridge.py)、[结果](reversible_dynamics_bridge_results.json)、[轮次检查](research_round_223_checks.json)；10项测试通过 |
| [224](research_note_224.md) | 合法边缘与维数饱和给二元张量识别；实际A→A过程完全正 | 不新增全效果权限；逐对取向尚待统一，全部CP实现权限未推出 | [代码](tensor_process_bridge.py)、[结果](tensor_process_bridge_results.json)、[轮次检查](research_round_224_checks.json)；11项测试通过 |
| [225](research_note_225.md) | 三体正边缘约束相对取向；全局相容表示中的全部实际A→B过程为CP | 仅剩共同转置自由；不推出全部CP实现权限或所有整数维系统 | [代码](global_orientation_bridge.py)、[结果](global_orientation_bridge_results.json)、[轮次检查](research_round_225_checks.json)；9项测试通过 |
| [226](research_note_226.md) | 同维辅助对满秩态的全部系综导引，迫使纯满Schmidt秩资源与辅助全部有限POVM | 权限属于指定辅助；目标侧全部测量、全部仪器及实际可逆纯化唯一性仍待核验 | [代码](steering_permission_bridge.py)、[结果](steering_permission_bridge_results.json)、[轮次检查](research_round_226_checks.json)；8项测试通过 |
| [227](research_note_227.md) | 独立组合排除受限辛控制群，双副本可逆群稠密于全部酉通道群 | C_path或闭群分支给精确全酉权限；不默许极限闭包，不给控制成本 | [代码](reversible_control_bridge.py)、[结果](reversible_control_bridge_results.json)、[轮次检查](research_round_227_checks.json)；8项测试通过 |
| [228](research_note_228.md) | 显式酉扩张实现全部有限CP仪器，并转移全部POVM权限 | 原合同给diamond范数任意逼近；加C_path或闭群给精确实现；只针对已有系统类型 | [代码](instrument_completion_bridge.py)、[结果](instrument_completion_bridge_results.json)、[轮次检查](research_round_228_checks.json)；9项测试通过 |
| [229](research_note_229.md) | 一个非恒定连续可逆操作路径，借助组合和简单Lie代数，推出所有对象的全酉及全仪器权限 | 非平凡Time足够；没有从原合同证明种子自动存在 | [代码](continuous_seed_bridge.py)、[结果](continuous_seed_bridge_results.json)、[轮次检查](research_round_229_checks.json)；8项测试通过 |
| [230](research_note_230.md) | 有限自适应实验的最终误差不超过各步骤diamond误差之和；后选择必须保留成功率分母 | 有限正容差下的实验逼近，不等于精确集合相同、相同成本或无限极限等价 | [代码](finite_protocol_closure.py)、[结果](finite_protocol_closure_results.json)、[轮次检查](research_round_230_checks.json)；7项测试通过 |

## 完整性与复算

[原清单](STAGE2_MANIFEST.json)冻结24份科学文件，[补充清单](STAGE2_CLOSURE_ADDENDUM.json)冻结另外8份；八轮共70项检查。两份清单均不改写。[迁移前快照](PRE_RENAME_FILES.json)记录47个原文件与版本备份；原README哈希由备份承接。

当前入口为[verify_archive.py](verify_archive.py)：核对迁移、旧证据及旧论文备份，核对当前整合来源与公式审核，并可复跑70项测试。旧verify_stage2及verify_stage2_completion脚本保留原哈希合同；不能直接用于已经更新为v1.1的顶层论文。

    python -B -X utf8 research_cognition_physics/archive_223_230/verify_archive.py
    python -B -X utf8 research_cognition_physics/archive_223_230/verify_archive.py --run-tests

结果见[integration_checks.json](integration_checks.json)。整合公式检查见[integrated_math_checks.json](integrated_math_checks.json)，复核脚本为[review_integrated_math.mjs](review_integrated_math.mjs)。旧[收尾核验](stage2_completion_checks.json)与[公式报告](completion_math_checks.json)仍是当时版本记录，不重写。

[整合生成脚本](integrate_stage2_paper.py)可从v1.0复算正文变更；当前核验检查生成内容与v1.1完全一致。核验不覆盖历史科学结果。

## 保留的开放问题

原F＋U＋C＋P是否自动保证某处存在非恒定可逆连续路径仍开放。229轮以该种子为充分条件给精确操作，Time非平凡时自动提供；230轮确认原合同的有限精度实验范围。尚未选择自然Hamiltonian、物理钟尺、无限维极限或引力。

本阶段已结束，不把这些后续问题改写为223—230轮尚未结项。当前方向见[研究方向](../research_direction.md)。
