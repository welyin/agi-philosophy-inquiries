"""Guarded navigation update for the full-network static-exchange bridge."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第436轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第436轮完成：** [完整网络的静态交换模拟桥]({prefix}research_note_436.md)'
           '核验已有低能局部模拟定理可覆盖433全网络，含任意未知输入及参考；'
           '具体构造六辅助首层，核全64列及误差，并用有理证书保证模拟后有限窗参与差异大于0.0114998。'
           '最终Heisenberg清单未编译；分组、目标依赖编码及非均匀权重仍为输入。'
           '6项检查，累计2044项；619份编号科学文件、654份保护证据。后继回到结构选择来源，三维与GR未完成。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第435轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix), 1)
    texts[path] = texts[path].replace('231—435轮', '231—436轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：从有效集体作用到同时工作的动态参与网络。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：从交换模拟能力回到内部结构选择。** 436已核成熟静态Hamiltonian模拟定理覆盖433完整有限网络，'
        '并具体检查首层和未知参考误差；不再为每个目标高阶项重做实现能力证明。'
        '后继回查SoCA功能定义及429的前提，检验哪些内部功能条件能选择主体分组、稳定关系块或参与模式。'
        '目标依赖的辅助、能隙、非均匀权重和初始编码仍要来源；'
        '若无法消去，给独立性证据或可检验的新候选。不得把可编译某种几何解释为必然采用该几何，'
        '也不将最终交换清单尚未编译说成现有定理没有覆盖有限模拟存在性。连续交换主线保持，三维与GR仍开放。'
        '\n\n**435后全网络任务（静态模拟存在性已对接）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为435／2038', '最新科学轮次与检查数为436／2044')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新435轮及累计2038项见本文开头', '最新436轮及累计2044项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[436科学核验](archive_231_/verify_static_exchange_network_bridge_round.py)、[436整合核验](archive_231_/verify_round436_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2038项；旧科学证据保持原字节', '第三阶段累计2044项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成435轮', '当前完成436轮').replace('当前435不作为预定终点', '当前436不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [435：'))
texts[archive] = texts[archive].replace(row, row+'\n| [436：全网络静态交换模拟](research_note_436.md) | 已有定理覆盖全目标；显式六辅助首层；完整未知参考误差与严格响应预算 | [代码](static_exchange_network_bridge_audit.py)、[结果](static_exchange_network_bridge_audit_results.json)、[核验](research_round_436_checks.json)；6项检查；最终交换清单未编译，结构选择仍开放 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[436整合核验](verify_round436_integration.py)、[整合记录](round436_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2038项科学检查；651份保护证据，其中本阶段编号科学文件616份', '当前2044项科学检查；654份保护证据，其中本阶段编号科学文件619份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 64.' not in texts[audit]
texts[audit] += '''
## 64. 第436轮：有限全网络的静态模拟已有覆盖，研究关口回到结构选择

2026-09-24。回读435和导航，追溯367—368对局部编码定义的引用，核Cubitt–Montanaro–Piddock完整稿Definition 23、Lemmas 25／39／40、Theorems 41／42及动力学界。此链是整个固定Hamiltonian的低能局部模拟，并非逐门控制的抽象普适性；适用于433完整六qubit三体实Hamiltonian。Heisenberg＝2SWAP−I，最终作用类型仍在429内。

本地代码完整展开25个Pauli项，保留6个成对Y项；逐项配置6个内部辅助构成12-qubit无Y首层，给全部作用清单、64列算符核验及保守Duhamel误差界。首层构造不等于最终Heisenberg清单已经编译，后者存在性依赖已核文献定理。普通复线性等距编码保留所有未知输入和参考，不采用反线性态共轭、输入估计或码外后选择。

对433原初态的参与差异追加60阶整数／有理证书：模型时刻1／2严格大于7／500。结合导数界、全链η＝ε＝10⁻⁴的静态模拟合同，窗口[0.499,0.501]的模拟差异严格大于0.0114998。此为定理条件下的完整网络信号保证；没有数值运行最终交换硬件，抽象4维合同例明确只作查错。

关闭的是有限目标的交换模拟存在性，不是自然架构来源。分组、辅助准备、能隙、目标依赖编码及非均匀实权重仍是输入；没有证明同一有限硬件可精确完成一切任务，更没有推出三维或GR。以后不为每个高阶项重复编译模拟器来充当生成进展；回到SoCA和429候选功能条件是否能选择实际分组／参与结构，以及这些输入是否独立。

6项检查、12个编号公式；累计2044项、619份编号科学文件、654份保护证据。旧科学文件和旧结果保持原字节，新结果只读复算；无图像检查或独立代理终审声明。目标活动，未触发空间阶段结项。
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
