"""Guarded navigation update for continuous three-body exchange relations."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第430轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第430轮完成：** [三方交换关系与连续读出]({prefix}research_note_430.md)'
           '接入已有交换编码，三个qubit首次形成非交换M₂关系块；给所有完整两体边缘相同、'
           '同一固定交换律下未来两体概率差为√3／2的具体见证。完全等强三角耦合则冻结关系态。'
           '6项检查，累计2008项；601份编号科学文件、636份保护证据。'
           '未把关系Bloch分量认作位置；下一项用同一交换机制构造内部读者与记录。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第429轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix), 1)
    texts[path] = texts[path].replace('231—429轮', '231—430轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：沿429交换规则研究连续自然演化中的内部关系。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：同一交换机制怎样产生内部记录。** 430已经找到三方关系qubit及持续演化的两体可辨差异。'
        '后继加入明确的内部读者和空白记忆，检验固定交换律能否把关系写到读者，并计入初态、耦合、扰动和参考。'
        '优先建立有限可读窗口，不擅加永久存储或精确停机。250的CZ记录不能直接替代该构造；'
        '关系与实际坐标的来源仍未闭合。\n\n**429后方向（继续保留）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为429／2002', '最新科学轮次与检查数为430／2008')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新429轮及累计2002项见本文开头', '最新430轮及累计2008项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[430科学核验](archive_231_/verify_exchange_relation_round.py)、[430整合核验](archive_231_/verify_round430_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2002项；旧科学证据保持原字节', '第三阶段累计2008项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成429轮', '当前完成430轮').replace('当前429不作为预定终点', '当前430不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [429：'))
texts[archive] = texts[archive].replace(row, row+'\n| [430：三方交换关系与持续读出](research_note_430.md) | 最小非交换关系块；同两体边缘而未来不同；等强冻结边界 | [代码](exchange_relation_audit.py)、[结果](exchange_relation_audit_results.json)、[核验](research_round_430_checks.json)；6项检查；内部记录器仍待生成 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[430整合核验](verify_round430_integration.py)、[整合记录](round430_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2002项科学检查；633份保护证据，其中本阶段编号科学文件598份', '当前2008项科学检查；636份保护证据，其中本阶段编号科学文件601份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 58.' not in texts[audit]
texts[audit] += '''
## 58. 第430轮：交换关系中的量子载体与固定连续响应

2026-09-24。接续用户要求沿429深入。本轮对接DiVincenzo等的三qubit交换编码和Bacon等的编码普适性，不将成熟编码重新命名为原创理论。与250的CZ指针、411的一般源作用代数、343的预定图几何去重后，集中于具体交换族中完整两体摘要遗漏的关系及其自然响应。

两个qubit的共同不变关系代数为C⊕C；三个为C⊕M₂。由交换算符本身构造P及关系X、Y、Z；Y来自两个关系的不对易性。完全相同的三组完整两体密度矩阵不能决定三体关系态。两个正交关系态在同一H=S₁₂＋S₂₃下，未来二体singlet概率差达到√3／2，另给有限演化读数重建第三分量的解析式。未知规范、逻辑及参考关联按同一整体酉保留。

三条耦合若完全等强，H在关系块为标量，关系不动。故配对强度、主体识别、初态和两体相加规则仍是明示输入；不能从作用形式相同推得空间、坐标或自动结构选择。“最少三个主体”也不等于物理三维。

6项新检查和10个编号公式通过。累计2008项；601份编号科学文件、636份保护证据。核验复算新JSON并保护全部旧证据；本轮没有独立代理终审，不将429的代理意见迁移为430已审。没有图像检查或旧实验重跑。

下一项从概率可辨性走向实际内部载体：只用同一交换机制，加上明确记账的读者初态与耦合，让关系信息进入读者，检查有限读数窗口、扰动及参考。若需外控关闭作用或新的测量相互作用，应列明，不能藏入“读取”一词。目标保持活动；空间和GR的来源未闭合。
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
