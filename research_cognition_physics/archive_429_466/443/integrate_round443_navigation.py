"""Guarded navigation update for actual signalling with coherent tree inputs."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第443轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第443轮完成：** [合法邻域通道与分叉图实际数据信号](%sresearch_note_443.md)'
           '保持440—442自主H，在六树部门接上主体数据编码与实际接收读数。'
           '未知远图相干及任意内部参考下，t＝1的接收迹距离有严格下界10⁻⁴；'
           '指定相干态可消去三阶振幅，信号概率由六阶改为十阶。'
           '已接入2024量子网络文献对旧偏迹的修正。'
           '6项检查，累计2089项；640份编号科学文件、675份保护证据。'
           '全零数据准备、编码和读出仍为合同；未生成三维或完成一般传播界。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第442轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—442轮', '231—443轮')
next_item = (
    '**下一项：可变伙伴的内部关系载体。** 443已在不改自主H的有限分叉模型内给实际通信证书；'
    '不继续优化六树概率、等待时刻或Taylor阶数。先复查433、437—440及既有端口工作，'
    '检验仅随主体及实际端口规模增长的内部载体，能否保持未知伙伴身份及叠加，'
    '并承载429交换和现有关系重组。438已给完整带标签匹配的N log N记忆量级，不能重算作新发现；'
    '新的问题是物理因子化与自主作用项，而非仅做部门维数压缩。'
    '编码若把关系操作变成全局多体作用，须计入访问与实施成本。'
    '一般分叉传播界保留为技术问题；初始关系选择、测距、三维及新主体接入仍开放。')
for path in (BASE/'research_direction.md', BASE/'RESEARCH_STATE.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：分叉重组中的实际数据信号。**'))
    texts[path] = texts[path].replace(old, next_item+
        '\n\n**442后通信任务（443已补有限全相干信号，一般界仍开放）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace(
    '最新科学轮次与检查数为442／2083', '最新科学轮次与检查数为443／2089')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新442轮及累计2083项见本文开头', '最新443轮及累计2089项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[443科学核验](archive_231_/verify_quantum_neighborhood_signal_round.py)、[443整合核验](archive_231_/verify_round443_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2083项；旧科学证据保持原字节', '第三阶段累计2089项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成442轮', '当前完成443轮').replace('当前442不作为预定终点', '当前443不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [442：'))
texts[archive] = texts[archive].replace(row, row+
    '\n| [443：合法邻域通道与分叉图实际数据信号](research_note_443.md) | 主体读出相容；相干改变信号阶次；任意远图及参考的有理通信下界 | [代码](quantum_neighborhood_signal_audit.py)、[结果](quantum_neighborhood_signal_audit_results.json)、[核验](research_round_443_checks.json)；6项检查；全零数据准备及读出为输入，非一般传播锥 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[443整合核验](verify_round443_integration.py)、[整合记录](round443_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace(
    '当前2083项科学检查；672份保护证据，其中本阶段编号科学文件637份',
    '当前2089项科学检查；675份保护证据，其中本阶段编号科学文件640份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 71.' not in texts[audit]
texts[audit] += '''
## 71. 第443轮：合法局部通道与完整有限数据信号

2026-09-24。继续429交互方向，接442“图距离改变不等于数据已送达”的缺口，未改440纯换边H。先对接Arrighi、Durbec、Wilson的Quantum networks theory（Quantum 8,1508,2024）：其第15页明确修正2017旧图偏迹不保正的问题。本文六树中的两图相干例复现最小特征值−1／2，是已知缺陷的当前模型核验，不称新发现；此前具体科学证明未使用旧偏迹作为通道。

对单个基态先拆分数据及全部边位，保留因子身份与空占用，得到等距V与普通偏迹通道。其后读出中心b数据，严格等于原始b因子的偏迹，覆盖未知图、数据与任意参考。384维部门全部147456个矩阵单位通过；动态图球的嵌套、图球搜索和路由的内部实施未由此推出。主体实际读数直接作用于固定b数据因子。

真实通信采用明确准备合同：数据全零，远图四维部门及其参考完全未知，发送者在0号主体作I或X编码，5号主体读激发。精确单激发部门36维，接收概率p就是两编码接收数据态的迹距离；含参考仍相等。未知远图不是先测量为经典图。J＝kappa＝1时三阶接收块Gram有22、4、2、0四特征值，特定实暗相干态的四阶振幅也为0，五阶非零；经典混合p＝7t^6／36＋O(t^8)，暗态p＝t^10／3600＋O(t^12)。一般耦合展开显示暗态三阶系数范数平方为J^4(J−kappa)^2，不将指定比例的抵消当普遍规律。

40阶精确有理多项式与||H||≤9的余项证书证明：在模型时间t＝1，全部远图相干及任意参考都有p≥10^−4。浮点最小值约0.0001051107仅为复核；完整384维演化与四维纠缠参考直接读数一致。没有图后选择；没有宣称覆盖其他未知主体数据、全规模传播、内部自主编码／读出、物理秒、光速或三维。

本轮关闭的是当前有限自主分叉模型的实际数据信号接口。原始6个数据qubit、15个潜在边qubit及端点资源仍在；压缩记录不是免费物理内存。下一项回到可变伙伴的内部承载：复用438已有N log N记忆结论，核物理因子化及自主相互作用，避免继续六树或Taylor参数优化。6项检查、12公式，累计2089项、640份编号科学文件、675份保护证据；旧文件保持原字节，阶段未结项。
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
