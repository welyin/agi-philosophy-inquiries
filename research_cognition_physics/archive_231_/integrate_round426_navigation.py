"""Guarded navigation update for round 426; frozen evidence remains read-only."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = BASE.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第426轮完成' not in texts[BASE/'RESEARCH_STATE.md'], 'Already integrated'

summary = ('**第426轮完成：** [有限控制预算与操作闭性]({prefix}research_note_426.md)将415的紧预算覆盖接到Baire与NSS，'
           '在统一操作拓扑闭合及精确有限预算覆盖下得到Lie群和连续最短成本，无需半幅输入。'
           '同一两个有界原语给可任意逼近、却无有限精确预算的弱块目标，说明闭性不能免费补入。'
           '8项检查，累计1985项；589份编号科学文件、624份保护证据。实际位置与三维仍未生成。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    marker = next(line for line in texts[path].splitlines() if '**第425轮完成：**' in line)
    insertion = ('> ' if path == ROOT/'README.md' else '') + summary.format(prefix=prefix)
    texts[path] = texts[path].replace(marker, marker+'\n\n'+insertion, 1)
    texts[path] = texts[path].replace('231—425轮', '231—426轮')

next_step = ('**下一项：认知实施要求中的完备性与实际关系。** 426从有限有界控制及闭合的精确可达群得到光滑操作结构，'
             '但统一操作拓扑、闭性及实际位置接口仍是输入。接续检查认知合同要求哪种误差控制和可实现性闭合，'
             '并核对它是否作用于主体间的定位关系；不把每个极限都有限成本可达升级为已接受公理。'
             '有限载体全酉Lie性质、旧端点商及415紧性直接复用；不继续脉冲优化或预算扫描。'
             '局部图册保留，像与透镜暂缓；三维及GR仍未完成。')
direction = BASE/'research_direction.md'
texts[direction] = texts[direction].replace('最新科学轮次与检查数为425／1977', '最新科学轮次与检查数为426／1985')
old_heading = '**下一项：实际关系上的重放与半幅一致性。**'
assert old_heading in texts[direction]
texts[direction] = texts[direction].replace(old_heading,
    next_step+'\n\n**425后的历史问题（半幅分支保留）：**', 1)
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新425轮及累计1977项见本文开头', '最新426轮及累计1985项见本文开头')

readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[426科学核验](archive_231_/verify_closed_control_round.py)、[426整合核验](archive_231_/verify_round426_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计1977项；旧科学证据保持原字节', '第三阶段累计1985项；旧科学证据保持原字节')

archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成425轮', '当前完成426轮').replace('当前425不作为预定终点', '当前426不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [425：'))
texts[archive] = texts[archive].replace(row, row+'\n| [426：有限控制预算与操作闭性](research_note_426.md) | 紧预算覆盖与Baire给局部紧，完整操作范数给NSS；弱块目标可一致逼近却无有限精确预算 | [代码](closed_control_lie_audit.py)、[结果](closed_control_lie_audit_results.json)、[核验](research_round_426_checks.json)；8项检查；未把操作维数当空间 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[426整合核验](verify_round426_integration.py)、[整合记录](round426_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前1977项科学检查；621份保护证据，其中本阶段编号科学文件586份',
    '当前1985项科学检查；624份保护证据，其中本阶段编号科学文件589份')
assert old_heading in texts[archive]
texts[archive] = texts[archive].replace(old_heading,
    next_step+'\n\n**425后的历史问题（半幅分支保留）：**', 1)
texts[archive] += '\n426以425为科学基线，只读复算：`python -B -X utf8 research_cognition_physics/archive_231_/verify_closed_control_round.py`。8项检查、16个编号公式；操作闭性不自动来自内部资源计账。\n'

audit = HERE/'spatial_premise_closure_audit.md'
assert '## 53.' not in texts[audit]
texts[audit] += '''
## 53. 第426轮：从实际控制预算接到光滑性，保留闭性的代价

2026-09-24。上一目标回合完成425的证明、代码、结果与冻结整合，记为实质进展。本回合读取导航及425后，独立检查过程细分与有限量子接口。程序周期、摘要自主更新及有限对象全酉Lie性质已被旧轮覆盖，取消重复候选；实际有限预算覆盖与操作拓扑的闭性有不同连接，形成426。

### 53.1 实际推进

[426](research_note_426.md)复用415的紧预算K_T，不重跑紧性实验。在同一个完整操作范数拓扑内，全部实际有限预算端点G=∪K_n若为Baire空间，某K_n有内点，便给局部紧。Banach代数范数的倍增不等式直接排除小子群，无须425的半幅输入；Hilbert第五问题给有限维Lie群。精确有限预算覆盖又使生成元的Lie代数张满整个群，以短共轭端点图和逆函数定理推出同一最短控制成本连续。

这套合同内Baire、局部紧与操作闭性等价，不将Baire宣传为免费或严格更弱要求。完整酉通道的diamond完备性已核，不与有限能量或逐探针拓扑混淆。结论首先属于实际完整操作群，尚须真实位置挠子或相应作用商接口才能接空间图册；群维数1、3、8均可出现。

### 53.2 一个不能被免费加入的极限操作

同一两个固定有界生成元H_X=⊕4⁻ʲσ_x、H_Z=⊕4⁻ʲσ_z，允许全部415控制。重复括号和成熟Li–Khaneja奇多项式方法，结合本文的Bernstein全区间误差界，证明目标⊕exp(−i2⁻ʲσ_z)在完整操作范数闭包内。它仍无法在任何有限预算下精确完成：第j块要求T≥2ʲ−2⁻ʲ/6。全部预算T≥1的协议更有统一通道误差下界1/(12T)，但不排除任意精度逼近，也不声称该界最优。

这里的∞量词由解析证明承担，有限截断Vandermonde及脉冲数值只核公式。无限载体、控制原语、权重、任意可测控制及预算均为明示模型输入，不称为完整认知反模型。该边界说明闭包新增端点未必属于任何旧K_T，而不是重复415的近回返参数样本。

### 53.3 交付与范围

8项新检查、16个编号公式通过，第三阶段累计1985项；589份编号科学文件、624份保护证据。两个独立代理分别终审通道拓扑／成本连续性和无限块闭包／资源下界，其中后者实际只读复算并核对正式JSON。旧科学文件、结果与核验保持原字节，未做文章图像检查。

后继核查认知定义真正要求哪种完备性与统一误差控制，以及这些操作是否改变同一个实际定位关系。不得将‘所有理想极限必须有限成本精确实现’新增为用户已接受公理。没有新的来源证明时不扩展脉冲最优控制、弱块参数或预算扫描。三维阶段仍未结项，完整GR目标保持未完成；342—343原形式不足性与强化内部实施范围继续分开。
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
