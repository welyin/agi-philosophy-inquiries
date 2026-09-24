"""Guarded navigation update for autonomous distance deformation on trees."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第442轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第442轮完成：** [分叉树的自主距离变形与完整量子界](%sresearch_note_442.md)'
           '在440未改动的自主H下，树形可真正改变，不需共同固定链；精确计数叶对距离增减，'
           '得到规模一致的加权范数和概率界，覆盖任意树相干、未知数据及参考。'
           '六主体给2／9与8／9相干区别及所有未知数据的正概率证书。'
           '6项检查，累计2083项；637份编号科学文件、672份保护证据。'
           '初始树与叶资格仍为输入；距离变形尚非完整数据信号界，三维与GR未完成。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第441轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—441轮', '231—442轮')
next_item = (
    '**下一项：分叉重组中的实际数据信号。** 442消去了距离变形证明中共同固定链的需求，'
    '但叶对变近的概率不是消息送达概率，数据可经中继传播。'
    '在同一自主树模型内，检验按实际关系定义的局部算符／通道截断能否覆盖未知图相干及参考，'
    '并接上主体实际读数。复用343、347、441及既有图动力学定理，不能倒置其局域性前提；'
    '不重跑距离枚举或把可编排换边路线当作自主行为。若仍需外加固定几何，记录准确缺口并转回关系载体来源。'
    '初始图选择、内部测距、三维和新主体接入仍开放。')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：传播桥梁能否离开预置路径次序。**'))
    texts[path] = texts[path].replace(old, next_item+
        '\n\n**441后一般图任务（已补分叉树距离变形，完整信号仍开放）：** '+old.split('** ', 1)[1], 1)
state = BASE/'RESEARCH_STATE.md'
old = next(line for line in texts[state].splitlines() if line.startswith('**下一项：实际关系的位移来源。**'))
texts[state] = texts[state].replace(old, next_item+
    '\n\n**421—422后方向（历史；最新进展见442）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为441／2077', '最新科学轮次与检查数为442／2083')
texts[state] = texts[state].replace('最新441轮及累计2077项见本文开头', '最新442轮及累计2083项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[442科学核验](archive_231_/verify_branching_tree_distance_round.py)、[442整合核验](archive_231_/verify_round442_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2077项；旧科学证据保持原字节', '第三阶段累计2083项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成441轮', '当前完成442轮').replace('当前441不作为预定终点', '当前442不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [441：'))
texts[archive] = texts[archive].replace(row, row+'\n| [442：分叉树的自主距离变形与完整量子界](research_note_442.md) | 无共同固定链的距离变形；精确转换计数；未知树相干、数据与参考界 | [代码](branching_tree_distance_audit.py)、[结果](branching_tree_distance_audit_results.json)、[核验](research_round_442_checks.json)；6项检查；树部门和叶身份仍为输入，尚非数据信号定理 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[442整合核验](verify_round442_integration.py)、[整合记录](round442_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2077项科学检查；669份保护证据，其中本阶段编号科学文件634份', '当前2083项科学检查；672份保护证据，其中本阶段编号科学文件637份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 70.' not in texts[audit]
texts[audit] += '''
## 70. 第442轮：自主分叉树的距离变形，不再需要共同固定链

2026-09-24。接441，保持440纯换边Hamiltonian不变，在固定有限N、有界度树、两个指定叶标签的严格不变部门内研究。树上的NNI及单次叶距改变至多1是成熟图论事实；本轮接上全部顶点标签、量子转换幅度和未知数据。没有把Arrighi–Martiel预设精确因果性的表示定理反用作因果性来源，也未移植到439单边增删合并模型。

精确计数：叶距d下降的转换至多2(C−2)(d−2)，上升至多2(C−1)(C−2)(d−1)。D=d−1、正反权重D及D^−1使非Hermitian部分范数≤B=2C(C−2)|kappa|，所有图对角数据作用相消。于是U†D^−2U≤exp(2B|t|)D^−2，覆盖任意树相干、未知数据和内部参考；初态支持d≥d0时，缩到d≤r的概率≤min(1,[(r−1)/(d0−1)]²exp(2B|t|))。不是先测树再作经典平均。

3413棵指定叶对的3至7点树逐转换核验，至6点与440生成元对照。六主体384维完整块核算符界；对称未知数据的确定远图与均匀相干远图在pi/6分别给2/9和8/9近距概率。一般未知数据在1/1000时有严格正概率下界181413961/98208100000000，全部参考保留。d0=101到r=2的全规模证书在1/100时给概率≤1/8800。

已有flip可按额外指定层次把梳树主干逐层减半，18点例直径从9变6，故不是固定图改标签。该可编排路线需选层与控制，未证明均匀自主F按此运行，不把它当自主快速通信。441的J=0反例继续提醒：图变化不等于主体间数据交换。

实际推进是消去距离变形证明中的共同路径模板；初始树、容量、叶资格、关系载体和换边项仍为输入。全图距离分布不是已实施的内部测距，亦非消息传输概率。下一项核同一自主分叉模型的完整数据读数及关系局域截断，不能用额外三维图案替代来源问题。6项检查、14公式，累计2083项、637份编号科学文件、672份保护证据；历史字节保留，三维与GR目标继续。
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
