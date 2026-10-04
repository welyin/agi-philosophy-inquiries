"""Guarded navigation update for exchange-block composition and triangle terms."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第435轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第435轮完成：** [交换关系块的组合与环路作用]({prefix}research_note_435.md)'
           '证明任意有限单态块图的二阶不同边交叉项为零，允许共享块上的非对易作用同时存在；'
           '指定三角在三阶产生精确系数−1／24的ZZZ项，并给覆盖未知编码态及参考的动力学误差界。'
           '整数证书与6项检查通过，累计2038项；616份编号科学文件、651份保护证据。'
           '连接图及端口权重仍为输入，未实现433完整网络或生成三维。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第434轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix), 1)
    texts[path] = texts[path].replace('231—434轮', '231—435轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：连续交换关系块的跨块组合。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：从有效集体作用到同时工作的动态参与网络。** 435已给二阶组合律和三阶环路作用，'
        '不需把多条边逐次隔离；但二阶生成元有双线性秩限制，三阶ZZZ也未自动补齐433的全部受控SWAP。'
        '后继核能否保持未知共享数据与相干关系，构成同时工作的参与更新；'
        '新增辅助块、指定作用图、准备与不同阶次的误差必须列清，不能恢复未记账控制表。'
        '图及端口权重的自然来源仍开放，不预置目标维数。\n\n**434后组合任务（已完成有限阶组合及环路证书）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为434／2032', '最新科学轮次与检查数为435／2038')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新434轮及累计2032项见本文开头', '最新435轮及累计2038项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[435科学核验](archive_231_/verify_exchange_block_composition_round.py)、[435整合核验](archive_231_/verify_round435_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2032项；旧科学证据保持原字节', '第三阶段累计2038项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成434轮', '当前完成435轮').replace('当前434不作为预定终点', '当前435不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [434：'))
texts[archive] = texts[archive].replace(row, row+'\n| [435：关系块组合与三阶环路](research_note_435.md) | 二阶共享块组合；三阶ZZZ整数证书；任意编码输入动力学误差界 | [代码](exchange_block_composition_audit.py)、[结果](exchange_block_composition_audit_results.json)、[核验](research_round_435_checks.json)；6项检查；作用图及端口仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[435整合核验](verify_round435_integration.py)、[整合记录](round435_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2032项科学检查；648份保护证据，其中本阶段编号科学文件613份', '当前2038项科学检查；651份保护证据，其中本阶段编号科学文件616份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 63.' not in texts[audit]
texts[audit] += '''
## 63. 第435轮：共享块可组合，高阶环路并非独立二块作用

2026-09-24。先读434及导航结果，复用已有低能相互作用支撑工具，直接证明指定单态码中的更强选择规则。任意有限块图、任意实端口权重下，一次跨块作用只激发两个端点至能量4；不同块对的激发部门正交，故二阶交叉项严格为零，包括共享一个块。二阶作用可相加但不必对易；真实固定总Hamiltonian同时作用于任意未知码态。

能力范围同样列清：在固定实码基中二阶只含I、X、Z，两逻辑比特的双线性系数秩至多2，而逻辑SWAP为3；这不是对所有高阶编码实施的不可能性。完全均匀跨块权重在单态码上无作用，端口选择仍是输入。

指定三块环路以三条原始跨块交换得到二阶纯标量，三阶含ZZZ系数−1／24。27个有序词分解为3个同边立方、6个三角和18个零项；另用整数码列、整数置换及Fraction独立核对系数。三阶非归一化交织映射给精确残差多项式，在τ＝ε³t固定时误差为O(ε)，覆盖所有码输入及内部参考。ε＝1／512、τ＝3π的解析界约0.051476；不声称做了极长时间4096维谱指数模拟。

本轮证明简单交换能够形成有效三块作用，同时指出逐对计算会漏掉的高阶影响。没有生成块划分或作用图，也没有自动实现433完整受控SWAP网络；后继核动态参与的同时组合与资源账，不将给定图或关系空间的维数认作物理三维。

6项检查、12个编号公式；累计2038项、616份编号科学文件、651份保护证据。434及更早证据保持原字节，新结果只读复算，无图像检查或独立代理终审声明。目标活动；空间阶段尚未结项，GR未完成。
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
