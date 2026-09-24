"""Guarded navigation update for round 428 and the new internal-rule proposal."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第428轮完成：**' not in texts[BASE/'RESEARCH_STATE.md'], 'Already integrated'

summary = ('**第428轮完成：** [统一内生演化与内部随机终止]({prefix}research_note_428.md)接续用户新原则候选：'
           '选择与触发由同一内部规则产生；原则的存在要求不等于已找到具体演化律。可数首次终止分支归约为固定通道，'
           '故共同可分正常程序、固定返回接口下，随机重试与增长仍不能实现全部精确酉且几乎必然有限终止。'
           '这些实施条件不是新原则本身；复用已有相位程序核清平均消费与预备容量。6项检查，累计1996项；'
           '595份编号科学文件、630份保护证据。具体完整规则、实际三维及GR仍未得到。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**427后范围复核（不增轮次）：**' in line)
    historical = old.replace('**427后范围复核（不增轮次）：**', '**427后范围复核（历史，不增轮次）：**')
    historical = historical.replace('当前只作概念澄清', '当时只作概念澄清').replace('科学状态保持427／1990', '当时科学状态为427／1990')
    insertion = ('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix)
    texts[path] = texts[path].replace(old, historical+'\n\n'+insertion, 1)
    texts[path] = texts[path].replace('231—427轮', '231—428轮')

next_step = ('**下一项：统一内生规则怎样实现原有权限和实际关系。** 用户允许探索统一内生演化原则，'
             '并追问所有精确操作的具体生成规则。428已区分原则、额外实施条件和未找到的完整演化律；'
             '不能把相位反馈实例或占位符T当作答案，也不能用该实例支持的操作偷换229的全部精确权限。'
             '后继先核对连续目标是否必须进入共同正常程序、返回接口为何成立，再构造同一内部过程中的'
             '目标、触发及关系操作。可分性、固定返回、几乎必然有限终止均未自动升为认知原则。'
             '空间主线保留实际关系读取、改变与跨基点比较；不继续相位程序参数扫描或重复不可编程定理。')
old_heading = '**下一项：有限测试怎样读取并改变实际定位关系。**'
for path in (BASE/'research_direction.md', HERE/'README.md'):
    assert old_heading in texts[path]
    texts[path] = texts[path].replace(old_heading, next_step+'\n\n**427后的空间接口问题（继续保留）：**', 1)

direction = BASE/'research_direction.md'
texts[direction] = texts[direction].replace('最新科学轮次与检查数为427／1990', '最新科学轮次与检查数为428／1996')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新427轮及累计1990项见本文开头', '最新428轮及累计1996项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[428科学核验](archive_231_/verify_internal_termination_round.py)、[428整合核验](archive_231_/verify_round428_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计1990项；旧科学证据保持原字节', '第三阶段累计1996项；旧科学证据保持原字节')

archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成427轮', '当前完成428轮').replace('当前427不作为预定终点', '当前428不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [427：'))
texts[archive] = texts[archive].replace(row, row+'\n| [428：统一内生演化与随机终止](research_note_428.md) | 首次终止CP求和接精确编程边界；候选原则与共同正常程序等额外条件分开 | [代码](internal_termination_audit.py)、[结果](internal_termination_audit_results.json)、[核验](research_round_428_checks.json)；6项检查；未找到完整内生规则或实际三维 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[428整合核验](verify_round428_integration.py)、[整合记录](round428_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前1990项科学检查；627份保护证据，其中本阶段编号科学文件592份',
    '当前1996项科学检查；630份保护证据，其中本阶段编号科学文件595份')
texts[archive] += '\n428以427为科学基线，只读复算：`python -B -X utf8 research_cognition_physics/archive_231_/verify_internal_termination_round.py`。6项检查、11个编号公式；已知概率程序明确归属原文，内部终止归约的附加条件单独列出。\n'

audit = HERE/'spatial_premise_closure_audit.md'
assert '## 56.' not in texts[audit]
texts[audit] += '''
## 56. 第428轮：新原则候选与具体内生规则的区分

2026-09-24。用户在第55节澄清以后明确表示：不要排除所有触发与选择由同一套底层内部规则产生成为新增认知原则，只要符合认知直觉。随后追问，在所有精确操作都由同一规则自行完成的要求下，这套规则究竟是什么、是否就是新的认知原则。因此第55节的“仅有选项疑问”不再是最新方向；不再等待原先二选一。

### 56.1 新增研究要求的准确地位

[428第1节](research_note_428.md)提出统一内生演化原则候选E：动作、选择和触发来自整体内部状态、记忆、目标和交互，同一规则允许随机性、设备增长及新任务。整体内部的伙伴可提供任务，所有相关信息与资源计账。它不要求继续推导产生基本规律的更低规律，也不暗加固定有限宇宙或固定有限Hamiltonian。

“存在统一内生规则”是一条候选原则；某个具体更新核、通道或H如果是选定而非推出，另算模型假设。占位符T并不是已找到具体律。用户的允许探索不等于已经证明该要求可与229全部精确权限、正常程序和有限内部完成同时成立，也不等于接受本轮为分析选用的每条技术条件。

### 56.2 本轮新增终止规约及边界

对同一内部规则的可数首次终止历史，正常CP分支和为迹不增CP映射；数学上补齐缺失的迹即得固定CPTP处理器。若所有目标程序都几乎必然有限步完成精确酉，补全项在这些程序上为零。精确编程定理于是要求不同目标程序支持正交。共同可分正常程序中至多容纳可数目标，故不能覆盖全部相位酉或全酉。

证明允许随机重试和增长内存，不要求共同截止时间、有限平均运行时间、固定H或能量条件。共同可分正常初态、全部目标依赖的初态承载、可数内部终止和固定返回接口，均是明确的额外实施合同。数学补全不要求有限检测永不终止；本结论不否定E的所有形式，也不是完整认知模型或三维反证。

345、406、409的既有范围与本轮有别；235／254的有限过程相容性、420的指定端点编译也未包含这一内部停止量词。两位独立审计均确认增量。Nielsen–Chuang工具和Vidal–Masanes–Cirac倍相位方案均明确引用原研究，不声称新基本定理或新程序算法。

### 56.3 完整代码证据与下一项

以目标无关受控线路复算概率相位程序；全部首次成功分支、失败通道和未知参考关联均检查。平均消耗趋2但N片预备容量不随之变成2；两个目标的完整前缀重叠趋零，按需生成也要核来源。正常程序假设以外的无限乘积表示保持开放，不直接当作解决方案。

6项新检查、11个编号公式、24组仪器实例通过；累计1996项，595份编号科学文件、630份保护证据。最终独立代理已用既有Python只读复算完整JSON，旧科学文件及历史检查不改写，不作图像验证。

下一项须构造或约束真正的内部演化机制，先核共同正常程序和返回接口的来源，不能重复不可编程范围审计来凑轮次。实际空间仍需同一关系的读取、改变及跨基点比较；相位控制实例不等于位置生成。三维阶段与GR目标保持未完成，本轮有实质新增量词结果，目标保持活动。
'''

prepared = {}
for path, content in texts.items():
    newline = '\r\n' if b'\r\n' in original[path] else '\n'
    prepared[path] = (content.rstrip()+'\n').replace('\n', newline).encode('utf-8')
for path in PATHS:
    assert path.read_bytes() == original[path], f'Concurrent edit: {path}'
for path in PATHS:
    assert path.read_bytes() == original[path], f'Concurrent edit: {path}'
    path.write_bytes(prepared[path])
    print(path.relative_to(ROOT))
