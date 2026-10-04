"""Guarded navigation update for short-path gates and resource-limited rewiring."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第439轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第439轮完成：** [短路径建联与有限容量重组](%sresearch_note_439.md)'
           '接入已有二步路径门，证明固定连通分区内严格无信号，端点转换范数有不随N增长的上界。'
           '同时完整分类C＝2部门：仅三顶点成分可改占用图，长链与长环图固定但数据仍交换。'
           '7项检查，累计2064项；628份编号科学文件、663份保护证据。'
           '初始分区及建联律仍为输入；后继核逐端点守恒的同时换边，三维与GR未完成。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第438轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—438轮', '231—439轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：建联资格怎样产生可持续的传播局域性。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：没有临时空槽时能否自主交换关系。** 439已证明短路径门给固定分区无信号及规模无关端点转换范数，'
        '也完整分类C＝2的三点可变、长链长环图冻结部门；不再扫描容量或重复星模型概率。'
        '先查已有同时撤换多条边、逐端点保占用的边交换机制，核在内点满槽的长链或饱和网络上是否仍可重组，'
        '并分别记明新生成元、未知数据、初图不变量和传播边界。复用436静态模拟桥但不将可模拟性等同于429已选择规则。'
        '初始连通分区、新主体接入、实际测距及三维仍开放；不预填C＝3或预置空间骨架。'
        '\n\n**438后短路径任务（已完成分类与资源界）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为438／2057', '最新科学轮次与检查数为439／2064')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新438轮及累计2057项见本文开头', '最新439轮及累计2064项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[439科学核验](archive_231_/verify_short_path_exchange_round.py)、[439整合核验](archive_231_/verify_round439_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2057项；旧科学证据保持原字节', '第三阶段累计2064项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成438轮', '当前完成439轮').replace('当前438不作为预定终点', '当前439不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [438：'))
texts[archive] = texts[archive].replace(row, row+'\n| [439：短路径建联与有限容量重组](research_note_439.md) | 固定分区无信号；端点范数规模一致；C＝2全部图部门分类 | [代码](short_path_exchange_audit.py)、[结果](short_path_exchange_audit_results.json)、[核验](research_round_439_checks.json)；7项检查；建联律与初始分区仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[439整合核验](verify_round439_integration.py)、[整合记录](round439_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2057项科学检查；660份保护证据，其中本阶段编号科学文件625份', '当前2064项科学检查；663份保护证据，其中本阶段编号科学文件628份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 67.' not in texts[audit]
texts[audit] += '''
## 67. 第439轮：短路径门给端点速率界，也暴露有限槽位下的重组障碍

2026-09-24。沿429交换方向，接438共享端点资源，核Hamma等0911.5075式(17)—(18)的L＝2路径门。门是共同邻居计数，多路径时不是投影；原文完整模型含资源跳跃，本轮只在438无跳跃、局部Q_i=C部门加入此门。N＝3、C＝2的216维资源关系空间到8维合法图基有精确整数交织，完整未知数据及参考由张量延伸保留。

加边已有二步路、删边仍留替代路，故连通分区投影守恒。在一个固定分区部门，合法空间按分量张量分解，H是各分量H之和；任意初始纠缠与参考、留在部门的局部迹保持操作均不能跨分量发信号。这是438空边全对影响的明确修复，但同时空图不启动、不同旧分量不能合并。它复用257已识别的接续边界，不将一般路由不可能性重新包装。

同一门还给端点转换范数≤|Omega|C(C−1)，由对称配置矩阵行和和二步路径数恒等式直接证明；固定C下规模一致。C＝2、N≥3的全合法部门范数精确为sqrt(2)|Omega|。该积极资源结果补438的sqrt(N−1)增长，尚不是实际度量光锥。

C＝1所有占用守恒，但相干叠加与数据关联仍可演化。C＝2每连通成分只能为路径或环；只有三顶点P3与K3互换，其他图固定。若初图有r个三顶点成分，非零Omega时可达配置是r颗三叶星的笛卡尔积、共4^r个；此数不是单条酉轨迹维数。N＝6全部1858个合法图分成1378个单配置、80个四配置、10个十六配置部门。未知数据不阻断最短路径最低阶恒等振幅；谱例只复用438星方法，不另为同样的2/9概率开轮次。

任意C的饱和无三角正则图结构冻结，数据交换仍活跃。反面范围由容量充足例约束：每分量s点且C≥s−1时，可经三角补边到完全图再回任意连通目标，此为源文献已有可达性结果的容量接口。不能据小容量限制否定所有局域模型、原文完整跳跃模型或认知原则。

后继先查成熟的多边同时交换，核其能否绕过临时空槽而保逐端点占用，以及额外规则与当前交换原语的关系；不静默增加基本控制或将编码可模拟性当自然选择。7项检查、15个公式；累计2064项、628份编号科学文件、663份保护证据。独立只读审查及本地复算完成，历史科学字节不变，空间与GR目标仍活动。
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
