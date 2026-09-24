"""Guarded navigation update: pair rule and continuous-evolution clarification."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第429轮完成：**' not in texts[BASE/'RESEARCH_STATE.md'], 'Already integrated'
summary = ('**第429轮完成：** [相同独立状态与具体交换规则]({prefix}research_note_429.md)'
           '在明示的同维、封闭可逆两体模型中，相同独立混态各自不动恰好选出部分交换；'
           '连续齐次生成元为aI＋bSWAP。纯态条件不足、相关输入条件过强分别给出精确边界。'
           '新功能条件尚未升为公理；6项检查，累计2002项；598份编号科学文件、633份保护证据。'
           '按用户最新提示，区分连续自然演化与任意精确程序实现，后者不再作为空间生成的先决门槛。'
           '用户已要求沿交换规则继续，优先检验多主体内部关系；具体完整宇宙律、三维与GR仍开放。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第428轮完成：**' in line)
    insertion = ('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix)
    texts[path] = texts[path].replace(old, old+'\n\n'+insertion, 1)
    texts[path] = texts[path].replace('231—428轮', '231—429轮')

next_step = ('**下一项：沿429交换规则研究连续自然演化中的内部关系。** 用户明确要求先探索完429方向。'
             '先由已选两体规则检验多主体关系自由度、可读比较和持续演化，再核其与实际坐标的桥梁。'
             '统一内生演化不自动要求宇宙充当任意精确程序处理器；229的操作定理及428的限定实现障碍保留，'
             '不将前者直接赋予无控制的自然流，也不将后者设成空间生成门槛。连续性、可逆性、齐次性各自列明，'
             '不暗加均匀微观钟或离散停机。429的相同独立状态不动仍为候选输入，强度及主体参与关系未被选出。'
             '复用343、250—255、398、411；不重扫预定图维数或精确程序参数。')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：统一内生规则怎样实现原有权限和实际关系。**'))
    texts[path] = texts[path].replace(old, next_step+'\n\n**428后实施问题（历史，限于精确处理器解释）：** '+old.split('** ', 1)[1], 1)
direction = BASE/'research_direction.md'
texts[direction] = texts[direction].replace('最新科学轮次与检查数为428／1996', '最新科学轮次与检查数为429／2002')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新428轮及累计1996项见本文开头', '最新429轮及累计2002项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[429科学核验](archive_231_/verify_agreement_exchange_round.py)、[429整合核验](archive_231_/verify_round429_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计1996项；旧科学证据保持原字节', '第三阶段累计2002项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成428轮', '当前完成429轮').replace('当前428不作为预定终点', '当前429不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [428：'))
texts[archive] = texts[archive].replace(row, row+'\n| [429：相同独立状态与交换规则](research_note_429.md) | 完整混态不动选部分交换；纯态与关联量词的边界；自然连续流和精确控制分层 | [代码](agreement_exchange_audit.py)、[结果](agreement_exchange_audit_results.json)、[核验](research_round_429_checks.json)；6项检查；关系与三维来源仍开放 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[429整合核验](verify_round429_integration.py)、[整合记录](round429_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前1996项科学检查；630份保护证据，其中本阶段编号科学文件595份',
    '当前2002项科学检查；633份保护证据，其中本阶段编号科学文件598份')
texts[archive] += '\n429以428为科学基线，只读复算：`python -B -X utf8 research_cognition_physics/archive_231_/verify_agreement_exchange_round.py`。6项检查、11个编号公式；额外认同条件与连续演化语义分别列明。\n'
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 57.' not in texts[audit]
texts[audit] += '''
## 57. 第429轮：具体交换规则与连续自然演化的范围修正

2026-09-24。前一轮428有实质结果；本轮提出并验证一个具体规则选择，非连续无进展状态。用户进一步提醒宇宙可能只需连续演化、不需执行任意指定精确动作，并明确要求沿429交互方式继续探索。

### 57.1 新结果与附加输入

[429](research_note_429.md)从同维、已识别、封闭可逆两体和所有独立同态边缘不动条件出发，以互信息为零、相邻混态差及已知二副本Schur–Weyl对易代数推出部分交换。连续齐次生成元为aI＋bSWAP，非零b及参与关系未被选出。新功能条件不是原FUCP或“持续认同”自动包含的公理；Ziman等的均匀化和交换机制明确归属原研究。

同一固定规则响应不同内部伙伴态，并保留完整未知参考信息于整体。qutrit有理例说明只保持相同纯态不够，效果差精确1／8；有关联的相同边缘例说明更强条件会使该连续交互平凡。单次局部响应的瞬时Hamiltonian不能升级为任意有限精确控制。

### 57.2 去重与验证

取消抽象类型骨架、未知输入无扰装配和动态返回对齐三种重复候选，分别复用342、177、428。裸全酉协变更新的分类另加无内部标记条件，不能替代E，亦未登记新轮。交换图维数、局部重置、平均场分别由343、390、350／356覆盖，不重复实验。

独立代理完成6项检查的只读全JSON复算，审查混态证明、参考扩张、关联边界和程序量词；其所指出的三处范围措辞已修正。科学核验另复算全部新结果并检查旧630份证据哈希。新增11个编号公式；累计2002项、598份编号科学文件、633份保护证据。未执行图像检查，历史科学文件和旧检查快照未改写。

### 57.3 用户连续演化提示的落实

统一内生要求应首先约束整体变化、内部控制和记录的共同来源，而非默认宇宙是一台任意精确通用处理器。连续与离散是一组区别，精确规律与近似控制是另一组区别；连续性本身既未被证实为宇宙底层性质，也不自动提供全酉控制。

回读229第5—6节，其定理使用完整操作选择、组合及实际连续可逆种子，并已声明不提供具体H和精度成本。故不改写229、不把已证明权限无条件解释成有限自主程序，也不把任意有限误差近似偷换成原精确结论。若以后采用有效／理想极限语义，须重核相关操作前提。230的有限实验误差工具保留复用。

428仍限制其精确程序合同，连续实现若满足该合同亦受限制；它不否定所有连续自然律，也不再作为空间生成必须先解决的门槛。429可直接解释为持续相互作用，θ是该模型演化参数，不必是外部命令；真实钟、参与关系和初态仍待来源。

### 57.4 下一项

按用户最新决定沿429探索：先考察若干主体通过同一交换机制形成什么内部可观测关系，以及持续演化怎样改变这些关系。优先对接既有交换相互作用、关系自由度及内部记录结果，不另加全球离散调度或停机。实际空间仍需关系的读取、改变、跨基点一致性；三维和GR均未完成，不触发归档补名与阶段结项。
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
