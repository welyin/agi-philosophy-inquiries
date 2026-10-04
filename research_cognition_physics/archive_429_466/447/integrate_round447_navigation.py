"""Guarded navigation update for internal swaps and exact distributed-record communication."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第447轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = '**第447轮完成：** [双寄存器伙伴记录与内部交换](%sresearch_note_447.md)每主体一个伙伴标记在两寄存器间交换，精确保留可恢复身份；同一互指势在单标记代码中给按记录配对的量子作用，覆盖未知匹配、逻辑态及参考。固定空值效果在[0.995,1.005]读出两消息的概率差严格大于0.1，无低能近似或相位跟踪。同时核低能相位码的读出边界及任意原始数据四阶改图反例。8项检查，累计2114项；652份编号科学文件、687份保护证据。互指势仍是跨主体输入，伙伴改变尚未与本处理代码合并，未生成三维。'
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第446轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—446轮', '231—447轮')
next_item = '**下一项：伙伴变化与未知逻辑处理怎样共用载体。** 447已给单标记代码的精确处理与固定读数，但可恢复伙伴M严格守恒；445的改配对和446的整块交换不能直接与该代码相加。先核哪些交换在当前单标记空间中闭合，恢复出的伙伴、内部逻辑和实际读数分别怎样变化；覆盖未知匹配／逻辑／参考，明确临时非法关系及作用代价。不继续扫描低能相位、提高Taylor阶数或重述436普适模拟。若只是重新命名主体或另写条件项，记录缺口不另开轮。互指势来源、实际测距、三维及新主体接入仍开放。'
for path in (BASE/'research_direction.md', BASE/'RESEARCH_STATE.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：内部子结构交换能否突破标签随载守恒。**'))
    texts[path] = texts[path].replace(old, next_item+
        '\n\n**446后子结构任务（447已补精确伙伴处理与固定读数）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace(
    '最新科学轮次与检查数为446／2106', '最新科学轮次与检查数为447／2114')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新446轮及累计2106项见本文开头', '最新447轮及累计2114项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[447科学核验](archive_231_/verify_internal_partner_response_round.py)、[447整合核验](archive_231_/verify_round447_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2106项；旧科学证据保持原字节', '第三阶段累计2114项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成446轮', '当前完成447轮').replace('当前446不作为预定终点', '当前447不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [446：'))
texts[archive] = texts[archive].replace(row, row+
    '\n| [447：双寄存器伙伴记录与内部交换](research_note_447.md) | 精确单标记代码；配对处理；固定空值信号；低能读出边界；原始数据反例 | [代码](internal_partner_response_audit.py)、[结果](internal_partner_response_audit_results.json)、[核验](research_round_447_checks.json)；8项检查；互指势与编码仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[447整合核验](verify_round447_integration.py)、[整合记录](round447_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace(
    '当前2106项科学检查；684份保护证据，其中本阶段编号科学文件649份',
    '当前2114项科学检查；687份保护证据，其中本阶段编号科学文件652份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 75.' not in texts[audit]
texts[audit] += '\n## 75. 第447轮：共同伙伴记录、内部交换与固定读数的精确接口\n\n2026-09-24。接续446子结构交换问题，复查430／434已有关系处理和436通用静态模拟，不重复首次条件响应或模拟存在性。沿445互指势h0，扩展共同空符号，每主体两份(N＋1)维寄存器，只加入主体内部SWAP。跨主体物理作用仍来自已输入的h0，不宣称从纯局部交换产生通信。\n\n主代码每主体只含一个伙伴j和一个⊥，逻辑0＝|j,⊥>、1＝|⊥,j>。两个寄存器联合精确保留伙伴M；单份A可暂为空。这不是把暂时A不互指当成遗忘。内部SWAP恰为逻辑X，A端空值效果恰为逻辑1投影；同一多重集合内可变化，突破446标签随载限制。全部未知匹配／逻辑／参考有精确等距交织，四主体三匹配48维代码在243维扩展空间中完全不泄漏。\n\n每个已记录pair的精确H＝2Δ(I−|00><00|)＋λ(X_i＋X_j)，含由原互指势直接投影而来的ZZ项，并非四阶生成。固定M时按pair张量分解；保持代码的局部TP干预不能向别的pair发信。M相干时保留完整受控量子演化，不把关系改成未声明的经典配对。\n\nΔ＝λ＝1，源0／1、目标0，始终读同一个目标A端空值效果。32阶Gaussian有理Taylor、||H||≤4及exp4<81给t＝1概率差>11／100；||[H,E]||＝1给差的Lipschitz常数2，从而[199／200,201／200]整个窗口严格>1／10。无需低能、后选或相位跟踪；准备、计时窗口和实际仪器仍为输入，没有空间钟或完整自治读出器的结论。\n\n保留另一低能码|j,j>、|j,⊥>及激发|⊥,j>作范围对照。精确能谱给连通相位χ<0，首项−λ^4／(4Δ³)，裸码每pair全时间误差≤2|ε|、全N≤N|ε|。其大迹距离需跨局部守恒部门的相干准备和相位读取，且最优效果会随时间转动；不能以此假称同一固定读出器。主代码避免的是这项接口缺口。\n\n任意原始D态不享受固定伙伴合同：A=M=(01)(23)、D=M′=(02)(13)，四个内部SWAP将两者对换，24条虚路径给跨匹配系数−(5／2)λ^4／Δ³。说明代码准备不是可省假设。主代码虽能处理与通信，M严格守恒，尚未结合445动态重配。\n\n8项检查、16公式，累计2114项；652份编号科学文件、687份保护证据。旧科学文件保持原字节，无图像。下一项核伙伴更新与未知逻辑能否共用当前单标记载体，先审计交换闭合与真实读数，不直接拼接不同代码、不继续相位调参。互指势来源、三维、GR继续开放；阶段未结项。\n'
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
