"""Guarded navigation update for conserved endpoint resources and regrouping."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第438轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第438轮完成：** [共享端点资源与自主重组]({prefix}research_note_438.md)'
           '对接已有边／资源转换机制，局部占用荷守恒使每个可测配置满足容量上限，'
           '三主体持续演化有精确2／9换伙伴概率；完整未知联合信息保留。'
           '同时证明空边仍有全对三阶数据影响，故占用稀疏不等于传播局域性。'
           '7项检查，累计2057项；625份编号科学文件、660份保护证据。'
           '配额、潜在载体及速率仍为输入；后继核建联资格，三维与GR未完成。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第437轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix), 1)
    texts[path] = texts[path].replace('231—437轮', '231—438轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：参与关系是否消耗同一份端点资源。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：建联资格怎样产生可持续的传播局域性。** 438已给共享占用荷、有限参与及自主重组的明确正例，'
        '并以全对三阶数据响应证明占用快照不是完整因果作用图；不再扫描容量与三节点概率。'
        '先核Hamma等已有短路径建联条件的启动和重组资格：是否依赖预置骨架、是否使空图或匹配部门冻结，'
        '以及是否有不预置空间的内部相遇接口。复用255—258协商账和436静态模拟定理，'
        '分别追踪潜在作用、平均边、每配置占用与实际可持续影响；不把控制通道或容量3隐去为免费输入。'
        '保留共享配额的正结果，建联范围、结构稳定和实际三维仍开放。'
        '\n\n**437后共享端点资源任务（已完成模型与边界）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为437／2050', '最新科学轮次与检查数为438／2057')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新437轮及累计2050项见本文开头', '最新438轮及累计2057项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[438科学核验](archive_231_/verify_endpoint_resource_exchange_round.py)、[438整合核验](archive_231_/verify_round438_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2050项；旧科学证据保持原字节', '第三阶段累计2057项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成437轮', '当前完成438轮').replace('当前437不作为预定终点', '当前438不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [437：'))
texts[archive] = texts[archive].replace(row, row+'\n| [438：共享端点资源与自主重组](research_note_438.md) | 占用荷给每配置容量；精确自主换伙伴；全对三阶响应揭示传播边界 | [代码](endpoint_resource_exchange_audit.py)、[结果](endpoint_resource_exchange_audit_results.json)、[核验](research_round_438_checks.json)；7项检查；资源律与建联资格仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[438整合核验](verify_round438_integration.py)、[整合记录](round438_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2050项科学检查；657份保护证据，其中本阶段编号科学文件622份', '当前2057项科学检查；660份保护证据，其中本阶段编号科学文件625份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 66.' not in texts[audit]
texts[audit] += '''
## 66. 第438轮：共享配额兼容重组，占用图仍不等于完整传播图

2026-09-24。接437与SoCA共享容量动机，对接Hamma等0911.5075式(16)的已有边／端点资源转化；删除资源跳跃、加入条件数据SWAP，明确这些是模型选择。各端点Q_i=f_i+Σn_ij精确守恒，初态Q_i=C使所有可测关系配置degree≤C。C=1部门精确等距到全部匹配，原9-qubit矩阵的32列整数交织成立；C=2另核资源部门。共同参数与对称资源初态还保相同独立数据的各自边缘不动，但不等于429已经选择新增规则。

在对称数据部门，三个单边状态与空图构成固定星；从12关系出发于π/2转为13、23的概率各2/9，原边5/9，空图0，未知对称数据及参考保持。一般两匹配间最低阶振幅由合法最短路径数乘数据恒等给出，故可重组而无需外部开关；这不是指定伙伴控制或不可逆选图。一般数据会与关系相关，只有完整联合未知信息可逆。

严格区分每配置容量和平均边支撑：后者仍可完整。显式潜在关系硬件仍为N(N−1)/2；完整带标签匹配接口需Theta(N log N)总量子位，但单条对称轨迹仅需floor(N/2)+1维，已分别核算。端点转换生成元在配额部门的范数恰为|Omega|sqrt(N−1)，不能把一次至多一位伙伴当作速率和硬件都规模无关。配额是占用荷，不是已推导的能量、功耗或熵。

进一步给空关系初态的数据通道三阶项−i Omega²J t³[ΣSWAP,rho]/3；已有容量不阻止任意潜在主体对的短时影响。三主体1/100时接收Y差有严格有理下界13/27000000；复用430关系方法加内部参考后，共同换轴不变的SWAP读数保留一半信号，准备、隔离和读取仍计账。没有物理距离或光速输入，不称为超光速。该结果把下一问题落到建联资格与传播稳定性，而非继续模拟有限参与的正例。

后继核成熟短路径条件是否可启动、重组，还是依赖预先骨架或冻结匹配；空间的来源仍需补齐。7项检查、15个公式；累计2057项、625份编号科学文件、660份保护证据。只读独立审查与本地证书均已完成，旧科学字节不变，目标活动，空间与GR未结项。
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
