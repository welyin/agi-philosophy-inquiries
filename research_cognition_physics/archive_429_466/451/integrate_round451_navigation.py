"""Guarded navigation update for round 451 and the user's SoCA organization question."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第451轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = '**第451轮完成：** [多主体比较字典的相容条件与SoCA接口边界](%sresearch_note_451.md)在明确的最小地址载体及逐对逐项合同下，证明置换字典恰对应2的幂主体数，相干字典恰对应实Hadamard存在性；给12主体相干、6主体扩容正例及逐对相容不等于全网络对称的反例。6项检查，累计2139项；664份编号科学文件、699份保护证据。按用户SoCA提醒，后继转核复量子有效主体的组织与保持：实际接口闭合、未知信息保留和再次组合；不把治理动机、底层复数或具体分层直接当成必要性证明。三维未生成。'
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第450轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—450轮', '231—451轮')

next_item = '**下一项：宏观复量子主体的组织与保持机制。** 用户说明SoCA源于社会治理与防止权力异化，并问多个认知单元组成宏观复量子整体需要哪些结构。先定义整体实际可访问的准备、读出、操作及再次组合接口，复用221、430和450核其在真实后继交互中的闭合与未知信息保持；仅找到复矩阵子块或底层复数不算新主体形成。SoCA保真、模块边界、按需共享、闭环及内部资源账作为功能候选，分别检验充分性与必要性，不先认定具体分层唯一。449接收—处理—关系重组的反馈保留为具体候选；完整字典不先作为普遍必需。先去重194—200、250—258及433—449，不重做经典闭环、一般模拟、抽象核空间定理或Hadamard阶数扫描。仍沿429自然交互路线；三维及GR目标未完成。'
for path in (BASE/'research_direction.md', BASE/'RESEARCH_STATE.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：多主体比较识别与自然交互的整体相容性。**'))
    texts[path] = texts[path].replace(old, next_item+
        '\n\n**450后比较识别任务（451已给指定合同的充要条件与范围）：** '+old.split('** ', 1)[1], 1)

texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace(
    '最新科学轮次与检查数为450／2133', '最新科学轮次与检查数为451／2139')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新450轮及累计2133项见本文开头', '最新451轮及累计2139项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[451科学核验](archive_231_/verify_global_comparison_frames_round.py)、[451整合核验](archive_231_/verify_round451_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2133项；旧科学证据保持原字节', '第三阶段累计2139项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成450轮', '当前完成451轮').replace('当前450不作为预定终点', '当前451不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [450：'))
texts[archive] = texts[archive].replace(row, row+
    '\n| [451：多主体比较字典的相容条件与SoCA接口边界](research_note_451.md) | 对易对合与Hadamard；12主体相干；6主体扩容；逐对与全网络区别 | [代码](global_comparison_frames_audit.py)、[结果](global_comparison_frames_audit_results.json)、[核验](research_round_451_checks.json)；6项检查；后继核宏观有效主体的组织与保持 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[451整合核验](verify_round451_integration.py)、[整合记录](round451_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace(
    '当前2133项科学检查；696份保护证据，其中本阶段编号科学文件661份',
    '当前2139项科学检查；699份保护证据，其中本阶段编号科学文件664份')

audit = HERE/'spatial_premise_closure_audit.md'
assert '## 79.' not in texts[audit]
texts[audit] += '''
## 79. 第451轮：全体比较字典的条件及SoCA组织问题

2026-09-24。接续450的固定识别来源。在N主体、恰N标记＋blank、每主体A/D两寄存器上，设F_ij＝U_j U_i†、U_0＝I，保持共同blank且同一F作用A/D。逐个隔离主体对，要求原始packet交换、互指记录、主动记录及均匀内部交换分别与扭转交换T_ij相容。完整packet相容与blank迫使F_ij²＝I；互指记录迫使F_ij|j>落在|i>射线。于是Ui为彼此对易的Hermitian对合，Ui|0>为第i个正交标记乘相位；反向亦成立。

置换字典的生成群是忠实、自由、传递的阿贝尔2群，故N＝2^m，XOR给反向构造。相干酉字典共同对角化给H diag(w) H^T＝I_N，方阵性迫使w＝1/N，得到实Hadamard HH^T＝NI。任意第一行全正的Hadamard由Ui＝H diag(H_i) H^T/N反向给字典；这是存在性充要条件，不是完整等价类分类。复用Paley12构造，九项幅度±1/3的列证明相干字典超出置换。未声称解决全部Hadamard阶数。

最小N6不可行但加到M8标记可行，原完整主体维数49增至81；显式等距嵌入保持未知态及任意参考，但比较可把旧标记送入备用能级，未给自主增长／新主体加入。N4全D空白见证(0,2,0,0)经T01变(3,1,0,0)，h_R能量4变2，排除把逐对隔离合同升为全网络交换对称。地址比特数不等于空间维数。

用户两次SoCA提示已落实：回读根工程方案及437已有资源审计，补全“完整主体、接口保真、按需共享、临时组队、预算、闭环”与当前结果／缺项对照。完整全体字典是本轮额外合同，SoCA原文没有要求它。用户说明该设计起于治理与防权力异化，进而提出宏观复量子主体需要哪些组织条件。后继按此方向研究有效主体的形成、保持与再次组合：先声明真正准备／读出／操作接口，再核同一真实交互下的闭合、未知信息及伙伴扩张。221操作闭合、430关系C⊕M2和450组合合同直接复用；单独找到复矩阵块不算新成果。保真、保留差异及协调不越界可启发候选约束，不直接推出相干、纠缠或SoCA具体层级的唯一必要性。经典工程功能与量子结构的区别、未知态不可免费复制及用户内部资源约束保留。

6项检查、13公式，累计2139项；664份编号科学文件、699份保护证据。独立只读审核修正ST记号后通过，旧科学文件与结果保持原字节，无图像。没有从429选出新势、没有实际生成比较字典、没有完成SoCA闭环或三维；完整GR目标仍开放。
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

