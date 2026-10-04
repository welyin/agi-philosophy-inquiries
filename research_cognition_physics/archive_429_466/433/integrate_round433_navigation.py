"""Guarded update for the exchange-participation source audit."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第433轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第433轮完成：** [动态量子参与及新增原语]({prefix}research_note_433.md)'
           '对接Quantum Graphity与演化图物质模型，在明确扩展中以相同参数、相同空白关系态产生依数据关系而异的激活概率。'
           '完整联合律线性且保持未知参考；也证明关系翻转与受控交换不是六个原始qubit的429交换原语直接给出的。'
           '6项检查，累计2026项；610份编号科学文件、645份保护证据。'
           '关系载体架构和更新仍为新增输入，未导出稀疏几何、三维或GR。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第432轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix), 1)
    texts[path] = texts[path].replace('231—432轮', '231—433轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：由真实内部关系变量产生交换参与条件。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：关系载体与更新能否来自已有交换结构。** 433已对接动态图工具，给同一线性整体律中的真实参与响应，'
        '同时核出新增的关系qubit、翻转及受控交换。后继先检验能否由430现有关系块承担这些变量，'
        '并在429交换机制内产生其更新；不能把编码普适性、另加边标签或逐门控制表当实际来源。'
        '若需新功能条件，应明确其最小内容与认知动机。暂不套入优选度数、相变或参数扫描来指定三维。'
        '连续自然演化及内部资源记账保持，空间与GR仍开放。\n\n**432后动态图任务（已完成有限实现与新增原语审计）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为432／2020', '最新科学轮次与检查数为433／2026')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新432轮及累计2020项见本文开头', '最新433轮及累计2026项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[433科学核验](archive_231_/verify_quantum_participation_round.py)、[433整合核验](archive_231_/verify_round433_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2020项；旧科学证据保持原字节', '第三阶段累计2026项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成432轮', '当前完成433轮').replace('当前432不作为预定终点', '当前433不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [432：'))
texts[archive] = texts[archive].replace(row, row+'\n| [433：动态参与与新增原语](research_note_433.md) | 同一空白关系态的状态依赖激活；线性整体过程；原始交换实现边界 | [代码](quantum_participation_audit.py)、[结果](quantum_participation_audit_results.json)、[核验](research_round_433_checks.json)；6项检查；关系载体及其更新仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[433整合核验](verify_round433_integration.py)、[整合记录](round433_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2020项科学检查；642份保护证据，其中本阶段编号科学文件607份', '当前2026项科学检查；645份保护证据，其中本阶段编号科学文件610份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 61.' not in texts[audit]
texts[audit] += '''
## 61. 第433轮：动态图工具能实现参与响应，但没有免费给出新原语

2026-09-24。前一轮432有实质来源边界。本轮回读导航及432，回查350、356—357、398、362和365，核对Quantum Graphity 2008与演化图Bose–Hubbard模型2009原文。前者的优选度数／环项、后者的边变量／转化及短路程限制均是输入，不能直接提升为认知到三维的证明。

给三个数据qubit和三个关系qubit，所有对使用相同Ω、Δ、g，关系初态全0。固定H包含关系翻转、关系偏置及n_e控制的数据SWAP。完整重标与数据共同换轴协变；空白关系块中的二阶、四阶算符系数为2I和−12I−4S_e，精确整数核验。数据初始关系不同使激活概率在四阶出现差异，同一过程也反作用于数据；不是预装非均匀数字边表，也不是将未知态期望值直接代入Hamiltonian。完整过程仿射并保持未知参考信息。

新增资源没有消失：每个潜在对已有关系载体与归属，空白态、翻转、偏置及受控交换都要来源。六个原始qubit的全部SWAP保持共同Z，但新H不保持，并含三体Pauli项，故不能说已经由429裸交换直接得到。可能的关系编码扩展仍开放，不以抽象编码普适性替代实际实现。去掉Ω时全部边位守恒、全0关系使数据不动，显示新增更新项承担了实质任务。

6项检查、10个编号公式。累计2026项、610份编号科学文件、645份保护证据；新结果只读复算，历史哈希保持，无图像检查或独立代理终审声明。此轮为明确附加模型中的正实现及原语审计，未证明新认知公理、稳定几何、三维或GR。

后继沿429继续核关系变量能否由430已有关系自由度承担及其更新来源。若必须增加新的功能条件，应清楚列出；不直接选择度数或温度相来获得目标维数。目标活动，当前未触发空间阶段结项。
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
