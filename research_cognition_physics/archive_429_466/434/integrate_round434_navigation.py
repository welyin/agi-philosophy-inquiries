"""Guarded navigation update for the continuous encoded-exchange response."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第434轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第434轮完成：** [持续交换的编码条件响应]({prefix}research_note_434.md)'
           '复用四体单态编码，以八个原始qubit上的固定两体交换产生关系读者的低能条件响应；'
           '给任意编码输入及参考的统一动力学误差界，有限时间窗的概率差严格大于0.729。'
           '不需新增基本翻转或受控作用，但分组、强弱耦合及码态准备仍为输入，尚非433完整网络。'
           '6项检查，累计2032项；613份编号科学文件、648份保护证据。三维与GR仍开放。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第433轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix), 1)
    texts[path] = texts[path].replace('231—433轮', '231—434轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：关系载体与更新能否来自已有交换结构。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：连续交换关系块的跨块组合。** 434已给一个条件响应原语的持续两体交换实现和误差界，'
        '不需逐门编译；分组、能隙、弱耦合及码态准备仍属输入。后继检验多条关系同时参与时，'
        '二阶虚路径会不会产生额外跨关系作用，能否保持共享数据及相干关系；'
        '不以逐对隔离或抽象普适性冒充433完整网络。不继续扫描434的ε，也不预置三维。'
        '连续自然演化及内部资源记账保持，三维与GR仍开放。\n\n**433后原语任务（已有有限有效实现）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为433／2026', '最新科学轮次与检查数为434／2032')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新433轮及累计2026项见本文开头', '最新434轮及累计2032项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[434科学核验](archive_231_/verify_encoded_exchange_response_round.py)、[434整合核验](archive_231_/verify_round434_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2026项；旧科学证据保持原字节', '第三阶段累计2032项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成433轮', '当前完成434轮').replace('当前433不作为预定终点', '当前434不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [433：'))
texts[archive] = texts[archive].replace(row, row+'\n| [434：持续交换的编码条件响应](research_note_434.md) | 八体交换有效条件作用；全未知输入误差界；无需精确停机的响应窗口 | [代码](encoded_exchange_response_audit.py)、[结果](encoded_exchange_response_audit_results.json)、[核验](research_round_434_checks.json)；6项检查；分组及能隙仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[434整合核验](verify_round434_integration.py)、[整合记录](round434_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2026项科学检查；645份保护证据，其中本阶段编号科学文件610份', '当前2032项科学检查；648份保护证据，其中本阶段编号科学文件613份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 62.' not in texts[audit]
texts[audit] += '''
## 62. 第434轮：新增基本作用可在有限有效模型中削去，组合来源仍待核

2026-09-24。接续433回读导航、笔记及结果，查旧文后对接Bacon等四体单态编码、Weinstein等编码耦合及Schrieffer–Wolff方法。原文工具按既有工作归属；未将量子点布局或低温噪声结论引入三维生成。

两个四体码块以强块内交换保持低能关系自由度，两条弱跨块交换通过二阶虚路径产生非零ZZ；块内固定交换小修正给关系更新与偏置。八个原始单元全部只用两体SWAP，14条见证权重为正、始终开启。单跨块路径只有二阶标量，双路径产生条件作用。用非归一化数学交织映射直接给Duhamel误差界，覆盖全部未知编码态和任意内部参考；不要求额外制备微扰基态或施加绝热开关。

ε＝1／512时，同一空白读者对两种数据关系的激活概率差约0.973608，解析界保证整个逻辑时间窗差大于0.729。低能有效模型越准，作用所需模型时间按ε的平方倒数增长。未宣称永久记录、物理时钟、码外任意输入保护或完整宏观读出。

本轮消去的是有限例中的新基本翻转／条件作用类型；并未消去四体分组、参与表、耦合层级及码态准备。也没有实现433全部三边共享数据结构。432等强全局冻结定理保持；此例使用明确非均匀权重。后继先核同时参与的多关系组合及虚路径串扰，不继续扫ε或直接指定维数。

6项检查、12个编号公式，累计2032项；613份编号科学文件、648份保护证据。433式(1)与式(5)两处漏写TeX反斜杠在434第7节列明勘误，不修改旧科学文件或哈希。新结果只读复算，无图像检验或独立代理终审声明。目标活动；物理位置、三维与GR尚未生成，未触发空间阶段结项。
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
