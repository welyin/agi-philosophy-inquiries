# 从认知闭环到物理：研究索引

更新：2026-09-15，已推进至第二十四轮。

## 项目定位

项目的总问题是：认知的记录、预测、行动、更新和主体组合，能否在不预设量子态与物理时空的前提下，约束出物理规律？

- [SoCA 工程方案](../认知联合体架构-AGI工程方案.md)提供认知组织原则。
- [研究方向讨论](../研究方向讨论：从认知架构到物理世界.md)提出从认知结构出发的转向。
- [当前研究方向](research_direction.md)将哲学动机、建模输入、数学结果和物理对应分开，是本目录的主线约定。
- [信息几何研究库](../research_information_geometry/research_note_01.md)保留八轮既有物理模型的兼容性检查与反例，作为后续验收工具。
- [第八篇假说稿](../article8_hypothesis.md)是历史论文稿，不随每轮模型试验自动改写。

当前所有可行认知模型仍有明确经典实现。它们的进展在于约束状态、更新和接口的自洽性；尚未导出量子理论或引力。

## 研究链条

| 记录 | 解决的问题 | 已知边界 | 复算脚本 |
|:---|:---|:---|:---|
| [起点](research_direction.md) | 实现预测—行动—反馈闭环 | 闭环本身允许经典实现 | [cognitive_loop.py](cognitive_loop.py) |
| [01](research_note_01.md) | 从未来可区分性定义预测状态，检验主体交流 | 本地等价不保证组合等价 | [predictive_states.py](predictive_states.py) |
| [02](research_note_02.md) | 区分静默更新、统计可逆和有备份恢复 | 连续性与非交换性仍不排除经典模型 | [reversible_inference.py](reversible_inference.py) |
| [03](research_note_03.md) | 完整读取容量及复制权限 | 容量有限不等于状态有限维 | [operational_capacity.py](operational_capacity.py) |
| [04](research_note_04.md) | 精确预测闭包与受控近似 | 有限工作记忆不能直接保证全协议精确闭包 | [predictive_closure.py](predictive_closure.py) |
| [05](research_note_05.md) | 有界读取与有限平移空间产生有限频率形式 | 单频、对比度、更新和组合仍需独立指定 | [bounded_responses.py](bounded_responses.py) |
| [06](research_note_06.md) | 弱对比下实现部分保留与精确闭包 | 摘要正性不足以保证实际核；最优性未定 | [outcome_updates.py](outcome_updates.py) |
| [07](research_note_07.md) | 用后续读取衡量旧区别，构造改进，证明受限最优与非唯一性 | 一般核只有上下界；最优保留未选出唯一更新 | [retention_tradeoff.py](retention_tradeoff.py) |
| [08](research_note_08.md) | 检验试图饱和一般上界的二原子核 | 摘要与点态核均合法，完整输出仍违反原准备集合 | [preparation_obstruction.py](preparation_obstruction.py) |
| [09](research_note_09.md) | 求已有原语随机混合的最优保留曲线，恢复零强度恒等更新 | 只证明混合类别最优；仍有一阶弱扰动与更新非唯一性 | [randomized_retention.py](randomized_retention.py) |
| [10](research_note_10.md) | 检验重复弱读取的完整记录与累积扰动 | 公开信息可趋零而状态仍改变；公开执行模式会改变结论 | [weak_readout_limit.py](weak_readout_limit.py) |
| [11](research_note_11.md) | 推出公开模式的条件跳跃过程，计算全部事件标记的信息 | 明确经典实现；事件参数不是物理时间 | [marked_event_limit.py](marked_event_limit.py) |
| [12](research_note_12.md) | 证明同一圆盘接口上，记录不能提供通用确定性恢复 | 反向推断与通道恢复不同；不包含输入相关备份 | [record_recovery.py](record_recovery.py) |
| [13](research_note_13.md) | 取得一般核纵向保留锐界，排除固定模型的全方向二阶扰动 | 纵向目标与旧横向目标不同；依赖固定准备宽度 | [longitudinal_bound.py](longitudinal_bound.py) |
| [14](research_note_14.md) | 缩小准备宽度，构造有累计信息的二阶弱扰动经典模型 | 改变了准备规则；本轮仅核对均值、记录矩和局部条件矩 | [shrinking_preparations.py](shrinking_preparations.py) |
| [15](research_note_15.md) | 完成固定方向的联合记录—状态路径极限证明 | 有限操作时段、弱收敛；没有总变差收敛或物理时间结论 | [joint_diffusion_limit.py](joint_diffusion_limit.py) |
| [16](research_note_16.md) | 由公开记录的区段相关性识别准备噪声与读取强度 | 总体可识别性不是单次精确估计；仍有无量纲自由参数 | [record_identifiability.py](record_identifiability.py) |
| [17](research_note_17.md) | 证明全部固定方向记录等价于二状态隐藏 Markov 模型 | 不能扩展为任意旋转探针的二状态实现；经典解释仍完整 | [two_state_filter.py](two_state_filter.py) |
| [18](research_note_18.md) | 多方向记录排除共同二状态实现，确定线性预测维数为 3 | 预测维数不等于隐藏标签数；未排除更大的经典模型 | [multidirectional_records.py](multidirectional_records.py) |
| [19](research_note_19.md) | 无限次精确无理角平移排除所有有限 Markov 实现 | 连续经典实现仍成立；有限记录前缀有有限实现 | [finite_rotation_obstruction.py](finite_rotation_obstruction.py) |
| [20](research_note_20.md) | 22 标签精确实现 fresh 多方向读取，并给出有限平移的完整记录误差预算 | 构造上界不是全局最小标签数；未要求操作等价准备共享隐藏编码 | [finite_hidden_readout.py](finite_hidden_readout.py) |
| [21](research_note_21.md) | 四准备见证检验准备非情境性，并量化必需的隐藏准备差异 | 新增条件不是既有认知原则的推论；原受限经典模型仍成立 | [preparation_contextuality.py](preparation_contextuality.py) |
| [22](research_note_22.md) | 整个圆盘单次读取的准备非情境可见度阈值恰为 2/π | 只实现一次最终读取；未实现 fresh 选择性更新 | [affine_preparation_threshold.py](affine_preparation_threshold.py) |
| [23](research_note_23.md) | 对比度 0.5 的五次顺序读取排除任意准备非情境扩展 | 一次读取却有四标签实现；多数后处理未必最优 | [sequential_contextuality.py](sequential_contextuality.py) |
| **[24](research_note_24.md)** | **构造低强度下完整的准备非情境连续更新核，保持所有有限自适应协议** | **η≤0.00174982 是充分区间；一般最优强度边界仍未确定** | **[noncontextual_update_kernel.py](noncontextual_update_kernel.py)** |

## 最新结论与下一步

最近继续完成四轮：

1. **第二十一轮明确新增条件的作用**：“操作等价准备具有相同隐藏分布”称为准备非情境性。加入它并尊重随机混合，四准备猜测任务的分数至多 3/4；原强读取给出约 **0.849882**，要求至少 **0.399528** 的隐藏奇偶分布距离。旧经典模型保留不可访问的准备差异，仍可实现这些统计；新增条件并非既有认知规则的推论。
2. **第二十二轮求出单次阈值**：整个圆盘的准备非情境单次读取实现恰需可见度 **β≤2/π**，原对比度对应 **η≤0.6432999034**。给出达到边界的连续模型；有限标签无法恰好达到边界。尚不能据此断言顺序更新可行。
3. **第二十三轮发现顺序障碍**：η=0.5 的所有方向单次读取有四标签准备非情境模型；但同一对象连续五次 fresh 读取的多数猜测分数为 **0.7717466112>3/4**，排除任意准备非情境的顺序扩展。结论有完整串递推、独立解析式及有理数余量验证。
4. **第二十四轮得到完整低强度构造**：固定原准备宽度，对 **η≤0.001749820194** 构造非负连续隐藏角更新核，逐分支保持整个唯一仿射准备分布，精确复现任意有限自适应读取与静默旋转。这个充分区间与 η=0.5 的排除之间仍有待解边界。

**下一轮主任务**：改进第二十四轮显式核的非负性证书或构造，与顺序完整记录的排除界比较，收窄准备非情境完整仪器的可行强度范围。区分当前核的极限与所有实现的极限，不能靠未违反的有限扫描宣称模型存在。准备非情境性的认知来源、有限标签与有限精度成本仍保留为独立问题。

**保留的未解问题**：第七轮固定 $\eta=0.5$ 的一般周期核横向最优保留区间仍是 **[0.8570324548, 0.8660254038]**。第十三轮解决的是另一纵向任务；第十四轮改变了准备集合。当前各轮仍有完整经典解释。

**接续方式**：已按用户要求在当前任务设置自动续研，每 30 分钟读取最新进度并接续。仅在有实质结果、失败或需要用户决定时通知；手动研究也以本索引的最新未完成问题为入口。

## 复算

在项目根目录执行：

```text
python -X utf8 research_cognition_physics/preparation_contextuality.py --write-results
python -X utf8 research_cognition_physics/affine_preparation_threshold.py --write-results
python -X utf8 research_cognition_physics/sequential_contextuality.py --write-results
python -X utf8 research_cognition_physics/noncontextual_update_kernel.py --write-results
python -X utf8 -m unittest discover -s research_cognition_physics -p "*.py"
```

依赖 Python 和 NumPy。本轮验证版本为 Python 3.12.14、NumPy 2.3.5；第二十一至二十四轮分别新增 10、10、15、10 项检查，目录全部 **324 项检查通过**。`--write-results` 先运行对应轮次检查，通过后再写入 JSON；目录回归命令不重写旧结果。

本机若 `python` 指向不可用的 Windows 应用别名，可在 PowerShell 使用本轮已有运行时：

```powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -X utf8 research_cognition_physics/noncontextual_update_kernel.py --write-results
```

这条绝对路径只是本次环境的运行入口，不是其他机器的安装要求。
