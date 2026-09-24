"""Guarded navigation update for round 427; frozen evidence remains read-only."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = BASE.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第427轮完成' not in texts[BASE/'RESEARCH_STATE.md'], 'Already integrated'

summary = ('**第427轮完成：** [有限概率测试与完整操作误差]({prefix}research_note_427.md)证明紧预算集合上有限测试可认证统一误差，'
           '测试拓扑另满足Baire时恢复426的统一拓扑及Lie结论。弱块源给明确菜单、无限尾界和完整统计资源账；'
           '取消预算时同一菜单可漏掉最大差异。5项检查，累计1990项；592份编号科学文件、627份保护证据。'
           '相干测试权限、Baire及实际位置接口仍为输入，尚未生成三维。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    marker = next(line for line in texts[path].splitlines() if '**第426轮完成：**' in line)
    insertion = ('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix)
    texts[path] = texts[path].replace(marker, marker+'\n\n'+insertion, 1)
    texts[path] = texts[path].replace('231—426轮', '231—427轮')

next_step = ('**下一项：有限测试怎样读取并改变实际定位关系。** 427补上有预算承诺的有限概率认证，并将统一拓扑的来源改写为'
             '明确的测试族与Baire条件；这些条件仍未由认知原则推出。后继必须在同一内部过程模型中交代实际关系的读取、'
             '改变和跨基点比较，不能只继续加测块数、缩小误差或优化脉冲。一般端点商、有限维层析及已有坐标算法直接复用。'
             '用户四种变化性质保持为认知动机，不等同完整位移合同；三维与GR仍未完成，像和透镜暂缓。')
old_heading = '**下一项：认知实施要求中的完备性与实际关系。**'
for path in (BASE/'research_direction.md', HERE/'README.md'):
    assert old_heading in texts[path]
    texts[path] = texts[path].replace(old_heading,
        next_step+'\n\n**426后的历史问题（完备性分支保留）：**', 1)

direction = BASE/'research_direction.md'
texts[direction] = texts[direction].replace('最新科学轮次与检查数为426／1985', '最新科学轮次与检查数为427／1990')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新426轮及累计1985项见本文开头', '最新427轮及累计1990项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[427科学核验](archive_231_/verify_finite_probe_round.py)、[427整合核验](archive_231_/verify_round427_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计1985项；旧科学证据保持原字节', '第三阶段累计1990项；旧科学证据保持原字节')

archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成426轮', '当前完成427轮').replace('当前426不作为预定终点', '当前427不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [426：'))
texts[archive] = texts[archive].replace(row, row+'\n| [427：有限概率测试与完整操作误差](research_note_427.md) | 紧预算上的有限菜单与无限尾界；测试Baire恢复统一拓扑，额外权限及成本明示 | [代码](finite_probe_topology_audit.py)、[结果](finite_probe_topology_audit_results.json)、[核验](research_round_427_checks.json)；5项检查；未生成实际三维 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[427整合核验](verify_round427_integration.py)、[整合记录](round427_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前1985项科学检查；624份保护证据，其中本阶段编号科学文件589份',
    '当前1990项科学检查；627份保护证据，其中本阶段编号科学文件592份')
texts[archive] += '\n427以426为科学基线，只读复算：`python -B -X utf8 research_cognition_physics/archive_231_/verify_finite_probe_round.py`。5项检查、11个编号公式；预算承诺、测试权限与Baire条件分别明示。\n'

audit = HERE/'spatial_premise_closure_audit.md'
assert '## 54.' not in texts[audit]
texts[audit] += '''
## 54. 第427轮：让有限读取在明确预算下认证完整操作

2026-09-24。回读导航、425—426和四种认知变化性质的已有登记后，不重新检验原猜想。独立范围审计确认：原形式合同不足性342—343不撤销，但固定单一H、永久输出、有限初始宇宙及全部理想极限有限成本实现均非用户已接受原则。420尚缺完整过程接口相容性，不能仅靠扩大输出对象自动补齐。

### 54.1 新的测试到拓扑连接

[427](research_note_427.md)保留415的真实紧预算集合，以合法有限实验概率定义分离测试拓扑。每个预算集合上，紧到Hausdorff的恒等映射为同胚；进一步用概率差的有限开覆盖，取得给定精度的有限认证菜单。这个结论无需Baire。

若全部精确有限预算端点在测试拓扑下另满足Baire，可数紧覆盖给一个预算集合的开内点，固定前后复合把局部拓扑一致扩展到全群。于是恢复426统一操作拓扑及条件性Lie结论。Baire不是认知原则的新推论，测试权限和精确预算覆盖均未省去；不存在无条件三维结论。

### 54.2 可复算的真实菜单及代价

在426弱块源增加静止相干参考，明确给跨块制备和二元效果权限。每块四项概率重建SU(2)第一列，前N块概率差加已知控制预算的无限尾界，给覆盖全部未知输入和参考的半diamond上界。块内层析无法读出相对参考相位，原XZ控制也不能自己产生这些跨块测试。

预算B=4的例子用16个设置给完整误差上界约0.0101456。若每个概率精度0.003、置信失败不超过0.01，保守独立抽样账需15580544次源调用；另外计算准备、读出、重置和存储，不将其并入单次预算4。代码未伪造已做的物理实验。无共同预算时，同一单Z脉冲族可令整个有限菜单与恒等一致、完整误差却为1，具体说明预算为何必要。

### 54.3 交付与仍然缺失的来源

5项新检查、11个编号公式通过；累计1990项、592份编号科学文件、627份保护证据。独立代理核对证明和实际效果并只读复算正式JSON；所有旧科学文件、结果和检查快照保持原字节，没有图像检查。

“连续变化、需要资源、幅度可调、方向可选”支持变化几何的研究，仍不独自选出三个独立方向或完整位移合同。后继集中于主体间真实关系的可实施读取、改变及跨基点比较，不继续弱块参数扫描或重复一般拓扑区别。

内部实施的在线量词仍未被正式追加。若以后采用，正确候选是同一历史事后接纳新任务、按统一规则增长已计账设备并保留旧未知关联，不必要求固定单一有限装置或永久保存所有输出。这个待定量词不阻断正向关系研究。空间阶段收尾和完整GR仍未完成，本回合有实质新结果，目标保持活动。
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
