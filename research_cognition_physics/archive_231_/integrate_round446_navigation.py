"""Guarded navigation update for whole-subject exchange and fixed-subject data reception."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第446轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = '**第446轮完成：** [完整主体交换与真实接收](%sresearch_note_446.md)用同一个完整地址＋数据SWAP同时更新伙伴与搬运数据，给精确条件置换及全部未知联合态／参考误差。标签关联数据代数严格守恒，完整地址回路的数据块只有标量恒等；四主体跨旧分量实际接收信号严格大于0.35，旧伙伴源权重始终不超过9／4096。互指关系尚非通信邻接。5项检查，累计2106项；649份编号科学文件、684份保护证据。罚能与接触权限仍为输入，未生成三维。'
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第445轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—445轮', '231—446轮')
next_item = '**下一项：内部子结构交换能否突破标签随载守恒。** 446已给完整主体交换的真实数据接口及结构性边界；不继续调该作用族的标量权重、窗口或虚过程阶数。先回查430关系代数、434持续响应及436模拟结论，区分已有内部比较／处理与尚未由内生关系选择的接触权限。候选须具体打破446标签关联数据守恒，同时交代伙伴一致性及未知输入；若只另写任意条件数据项或重述普适编译，不另开轮。互指罚能来源、参与选择、实际测距、三维及新主体接入仍开放。'
for path in (BASE/'research_direction.md', BASE/'RESEARCH_STATE.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：完整主体交换怎样同时搬运数据与改变关系。**'))
    texts[path] = texts[path].replace(old, next_item+
        '\n\n**445后完整主体任务（446已补真实数据接口及随载边界）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace(
    '最新科学轮次与检查数为445／2101', '最新科学轮次与检查数为446／2106')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新445轮及累计2101项见本文开头', '最新446轮及累计2106项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[446科学核验](archive_231_/verify_whole_subject_exchange_round.py)、[446整合核验](archive_231_/verify_round446_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2101项；旧科学证据保持原字节', '第三阶段累计2106项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成445轮', '当前完成446轮').replace('当前445不作为预定终点', '当前446不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [445：'))
texts[archive] = texts[archive].replace(row, row+
    '\n| [446：完整主体交换与真实接收](research_note_446.md) | 精确条件置换；标签数据守恒；未知联合接口；跨分量实际信号；旧伙伴禁则 | [代码](whole_subject_exchange_audit.py)、[结果](whole_subject_exchange_audit_results.json)、[核验](research_round_446_checks.json)；5项检查；罚能和潜在接触仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[446整合核验](verify_round446_integration.py)、[整合记录](round446_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace(
    '当前2101项科学检查；681份保护证据，其中本阶段编号科学文件646份',
    '当前2106项科学检查；684份保护证据，其中本阶段编号科学文件649份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 74.' not in texts[audit]
texts[audit] += '\n## 74. 第446轮：完整主体交换的实际数据通道与随载守恒\n\n2026-09-24。承接445的完整数据接口。440／441已有动态标签与真实读数，但另输入沿占用边的数据SWAP；本轮没有该项，也不输入固定链。把445每个地址SWAP改成完整主体地址＋数据SWAP，其他互指罚能与全对权重保持。429只支持交换项形式，罚能、分组与潜在接触仍为输入。\n\n在每地址身份一次的精确不变部门，令Tπ按当前持有者排列数据，C＝Σ|π><π|⊗Tπ。每个完整SWAP与地址SWAP严格共轭，故H*＝C(H445⊗I)C†。此式对任意固定逐对权重及任意地址对角势仍成立。任意未知代码／数据／参考直接继承445算符误差与泄漏界，无需新的准备或扫描。\n\n固定主体b的数据O在变换后读取π(b)标签；随标签j的B_j(O)则与H严格对易。标签数据的完整矩阵代数保持，任意同端点地址路径给同一数据置换，闭合地址回路数据块为标量恒等。因此不能只通过调标量接触权重摆脱随载结构。B_j对易性也适用重复地址原空间，但独立单因子代数及C描述只用于每身份一次部门。\n\n已知初始匹配M0时，任意关联数据／参考的实际接收边缘为Σ_a q_ba ρ_aR，q_ba按π(b)＝M0(a)的全部地址结果加总。地址测量概率与数据无关，完整地址边缘态并非如此；任意未知图相干不能改作同一经典混合。初始源a的旧伙伴b＝M0(a)要求最终b自指，合法互指代码中不可能，所以该来源权重≤码外概率。\n\n四主体初始01、23、源0、接收3，在445同一ε＝1／128与τ∈[1.9,2.1]窗口，用既有严格分数证书得到实际接收信号>7／20；旧伙伴1对所有时间≤9／4096。实际零／一消息的固定接收者迹距离就是此权重，不是图变近概率；未知量子消息给有噪通道而非完美传送。互指因此尚非通信邻接，没有生成空间局域性或三维。\n\n5项新检查：精确完整SWAP共轭、联合二阶与误差、主体读数及标签守恒、未知数据／参考通道、实际信号与旧伙伴禁则。384维整数与谱证书、独立审查通过。累计2106项、649份编号科学文件、684份保护证据。下一项回查430／434内部子结构交换，候选必须明确突破随载守恒并交代关系一致性，不能重述436通用模拟或另加无来源条件数据项；阶段未结项。\n'
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
