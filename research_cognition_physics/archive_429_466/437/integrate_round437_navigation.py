"""Guarded navigation update for the participation-capacity audit."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第437轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第437轮完成：** [关系容量与参与稀疏性]({prefix}research_note_437.md)'
           '复用singlet单配性，证明强相关预算不等于参与数上限；'
           '433任意规模、任意输入的短时参与度有线性下界，对称数据部门给精确稠密图分布。'
           '32主体至多三邻居概率小于10⁻⁷⁰；潜在边硬件及准备能量按N²计账。'
           '6项检查，累计2050项；622份编号科学文件、657份保护证据。'
           '只否定当前规则自动稀疏，后继检验共享端点资源，三维与GR未完成。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第436轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix), 1)
    texts[path] = texts[path].replace('231—436轮', '231—437轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：从交换模拟能力回到内部结构选择。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：参与关系是否消耗同一份端点资源。** 437已证明433原规则的激活量不受数据qubit单配性限制，'
        '并给保留未知对称数据及参考的精确稠密反例；无需再扫描其规模或等待稀疏相。'
        '检验在同一固定连续律中，建立关系占用两端内部资源、解除关系归还的候选机制，'
        '能否同时给参与上限、未知信息保留及重组。共享账、容量常数、潜在边硬件和耦合速率分别核算；'
        '任何新增守恒量仍须标为建模输入，不预设容量3，不以有限度数代替三维。'
        '436的有限模型静态交换模拟结论直接复用，不重复编译通用模拟器。三维与GR仍开放。'
        '\n\n**436后结构选择任务（有限容量候选已核）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为436／2044', '最新科学轮次与检查数为437／2050')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新436轮及累计2044项见本文开头', '最新437轮及累计2050项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[437科学核验](archive_231_/verify_exchange_participation_capacity_round.py)、[437整合核验](archive_231_/verify_round437_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2044项；旧科学证据保持原字节', '第三阶段累计2050项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成436轮', '当前完成437轮').replace('当前436不作为预定终点', '当前437不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [436：'))
texts[archive] = texts[archive].replace(row, row+'\n| [437：关系容量与参与稀疏性](research_note_437.md) | 强相关预算与参与量区别；全规模激活下界；精确稠密反例及内部资源账 | [代码](exchange_participation_capacity_audit.py)、[结果](exchange_participation_capacity_audit_results.json)、[核验](research_round_437_checks.json)；6项检查；共享端点资源待检验，三维仍开放 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[437整合核验](verify_round437_integration.py)、[整合记录](round437_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2044项科学检查；654份保护证据，其中本阶段编号科学文件619份', '当前2050项科学检查；657份保护证据，其中本阶段编号科学文件622份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 65.' not in texts[audit]
texts[audit] += '''
## 65. 第437轮：强关系预算并非参与上限，稀疏性需要真正共享的容量

2026-09-24。接436回查SoCA的有限工作空间、稀疏激活、代谢与临时组队；工程设计给功能动机，不直接指定量子资源律。复用235的单配性来源与282的三邻居输入边界，避免重报不可克隆或图度数不决定维数。

标准星块谱界给Σ(2p_singlet−1)₊≤1；固定强度阈值才给伙伴数上限，且可有任意多逐渐变弱的纠缠伙伴。这不约束433另一批关系qubit的激活。针对该模型，局部Duhamel比较消去巨大环境范数，任意N、任意未知数据及参考均有t＝1／2时每边激活≥121／2304，平均度随N增长。

完全对称数据部门还给全时间精确因子分解，未知数据和参考原样保留；在π／(2√2)读取关系表得到G(N,1/2)分布。32主体最大度≤3的概率由精确二项式尾界限制到约1.015×10⁻⁸²。大规模结论来自解析恒等式和整数证书，未冒充运行528-qubit矩阵。436的全状态模拟误差同样保留这个稠密事件，故通用交换模拟能力不能自动筛选空间。

内部资源账揭示：模型已预置N(N−1)／2个关系载体，空白准备相对对称部门基态具有(√2−1)N(N−1)／2能量。总平均能量虽守恒为零，正占据代价可被翻转项负期望抵消；没有外部泵浦，也没有每主体有限消费份额。当前自动稀疏候选失败，并不排除增加明示共享资源后的正向机制。

下一项具体检验建边占用两端内部资源、解边归还的固定连续律，追踪未知数据、重组、容量及潜在硬件账。新增守恒量或容量都是待解释模型输入，不预设3；有限度数不等于空间三维。6项检查、12个公式；累计2050项、622份编号科学文件、657份保护证据。旧科学字节不变；目标活动，未触发阶段结项。
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
