"""Guarded navigation update for round 448 common-carrier and protection audit."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第448轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = '**第448轮完成：** [伙伴重组与内部处理的共同载体](%sresearch_note_448.md)将完整主体交换与内部交换接入同一单标记空间，精确保留未知伙伴、逻辑和参考并追踪实际读数。四主体给λ＝g＝1、t＝1／8下对每个Δ>0的非互指概率>1／40反例，说明活跃端势不能统一保护未知逻辑；完整记录势有原始空间正性并恢复445保护，但均匀内部交换精确解耦，失去447条件逻辑作用。6项检查，累计2120项；655份编号科学文件、690份保护证据。两种势及接触仍为输入，未生成三维。'
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第447轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—447轮', '231—448轮')
next_item = '**下一项：记录保护与条件处理的兼容条件及作用来源。** 448已建立完整伙伴／逻辑共同载体；活跃端势的全Δ反例和完整记录势的精确解耦说明两种功能不能混用。先核同一载体中什么对角能量条件保护全部未知逻辑又保留伙伴决定的响应；若需独立作用，交代系数与整体内部代价，不能藏入“稳定性”。不扫描Δ，不把暂时非互指禁令增为公理，不重复434—436通用模拟。只在新增兼容性定理、严格反例或可检验来源约束时新增轮次。实际接触、定位、三维及新主体接入仍开放。'
for path in (BASE/'research_direction.md', BASE/'RESEARCH_STATE.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：伙伴变化与未知逻辑处理怎样共用载体。**'))
    texts[path] = texts[path].replace(old, next_item+
        '\n\n**447后共同载体任务（448已完成兼容表示与保护边界）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace(
    '最新科学轮次与检查数为447／2114', '最新科学轮次与检查数为448／2120')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新447轮及累计2114项见本文开头', '最新448轮及累计2120项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[448科学核验](archive_231_/verify_joint_marker_dynamics_round.py)、[448整合核验](archive_231_/verify_round448_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2114项；旧科学证据保持原字节', '第三阶段累计2120项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成447轮', '当前完成448轮').replace('当前447不作为预定终点', '当前448不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [447：'))
texts[archive] = texts[archive].replace(row, row+
    '\n| [448：伙伴重组与内部处理的共同载体](research_note_448.md) | 完整共同载体；实际读数；全Δ保护反例；完整记录势正性及精确解耦 | [代码](joint_marker_dynamics_audit.py)、[结果](joint_marker_dynamics_audit_results.json)、[核验](research_round_448_checks.json)；6项检查；保护与处理的作用来源仍开放 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[448整合核验](verify_round448_integration.py)、[整合记录](round448_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace(
    '当前2114项科学检查；687份保护证据，其中本阶段编号科学文件652份',
    '当前2120项科学检查；690份保护证据，其中本阶段编号科学文件655份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 76.' not in texts[audit]
texts[audit] += '''
## 76. 第448轮：共同单标记载体、全Δ反例与完整记录势

2026-09-24。完成447提出的伙伴重组／未知逻辑共用载体任务：每主体(j,⊥)或(⊥,j)，标记π从匹配扩展到全部置换，保留自指、长环及所有非互指中间态。完整双寄存器主体SWAP与内部SWAP严格闭合N!2^N代码；N＝4为384维，含48维匹配／逻辑空间。半数单跨寄存器交换离开此代码，未据此禁止虚过程。

受控Tπ同时变换作用和实际主体空值读数。完整主体交换变为标记Γ，均匀内部交换仍为ΣX，活跃端互指势因互指对无序集合不变而保留。g＝0精确恢复447；共同载体保留未知伙伴／逻辑／参考，但不等于互指子空间始终不变。

N＝4、λ＝g＝1、t＝1／8，任一指定匹配基态乘|1111>给全Δ>0的反例。去掉内部L仅作比较，此部门h_A＝4I、标记自由Γ；Γ的精确多项式与二／四阶矩给奇置换概率sin²(6t)／12＋3sin²(2t)／4。恢复非零L后，压缩到全1的Dyson奇数阶为零，误差≤cosh(1／2)−1≤6／47。有限有理sin下界得非互指且全1概率>582169／22090000>1／40。它否定随Δ增大对全部未知逻辑统一保护的主张，不把非互指视作信息遗失或禁态；不是扫描趋势或全输入泄漏下界。

完整记录势采用互斥R_i(j)=|j,⊥><j,⊥|+|⊥,j><⊥,j|，原始空间每主体至多一条互指边，故h_R非负且偶数N能隙2；未藏入全局代码投影。在同一载体中C†H_RC=(Δh_R+gΓ)⊗I＋λI⊗ΣX，均匀λ是范围条件。445全时间、未知代码／参考保护界因此适用，但447条件逻辑ZZ消失；实际主体运输通信仍可存在。

两种势分别承担完整记录保护与伙伴条件处理，并非同一个未说明的“稳定性”要求。相加可列候选，但属于额外模型输入，不等于429已经生成，亦不为此重复436通用模拟。下一项只核新的兼容性条件、独立作用代价及可检验来源约束，不扫描罚能。

6项新检查、12公式，累计2120项；655份编号科学文件、690份保护证据，旧科学文件保持原字节。独立只读审核通过，无图像。实际接触、定位、三维、GR仍开放，阶段未结项。
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
