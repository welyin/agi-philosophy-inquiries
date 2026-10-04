"""Guarded update for the exchange-participation source audit."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第432轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第432轮完成：** [交换参与关系的来源边界]({prefix}research_note_432.md)'
           '证明任意有限N、d的固定两体交换中，全部共同换轴不变关系冻结，当且仅当所有可能主体对等强。'
           '改初态不能解除该恒等式；相同独立态不动允许任意权重，仍未选参与表。'
           '直接用关系期望值作权重的特定反馈有精确1／4仿射性违例。'
           '6项检查，累计2020项；607份编号科学文件、642份保护证据。'
           '保留429—431正结果，下一项检验真实内部关系变量怎样产生参与结构；三维与GR仍开放。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第431轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary.format(prefix=prefix), 1)
    texts[path] = texts[path].replace('231—431轮', '231—432轮')
for path in (BASE/'research_direction.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：交换参与关系和初态资源的内部来源。**'))
    texts[path] = texts[path].replace(old,
        '**下一项：由真实内部关系变量产生交换参与条件。** 432已核均匀参与的精确边界：固定两体交换中，'
        '所有共同换轴不变读数冻结恰好对应全体等强，初态不能解除。'
        '不能直接把未知态期望值代为耦合而冒充固定量子过程。优先对接动态量子图和真实关系载体，'
        '核控制状态、生成与更新、429合同兼容性及未知输入参考。仅给边再贴一个控制标签不算来源已解决。'
        '复用350、356—357及398的内部源范围；不重做均匀扫描、有限链读数峰值或一般非线性反例。'
        '参与结构、初态资源、坐标及三维仍开放。\n\n**431后来源任务（已完成均匀版本边界）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace('最新科学轮次与检查数为431／2014', '最新科学轮次与检查数为432／2020')
texts[BASE/'RESEARCH_STATE.md'] = texts[BASE/'RESEARCH_STATE.md'].replace('最新431轮及累计2014项见本文开头', '最新432轮及累计2020项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：', '当前复算与冻结入口：[432科学核验](archive_231_/verify_exchange_participation_round.py)、[432整合核验](archive_231_/verify_round432_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2014项；旧科学证据保持原字节', '第三阶段累计2020项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成431轮', '当前完成432轮').replace('当前431不作为预定终点', '当前432不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [431：'))
texts[archive] = texts[archive].replace(row, row+'\n| [432：交换参与关系的来源边界](research_note_432.md) | 任意N、d的等强冻结充要条件；同态条件不选权重；特定均值反馈的仿射性违例 | [代码](exchange_participation_audit.py)、[结果](exchange_participation_audit_results.json)、[核验](research_round_432_checks.json)；6项检查；内部参与结构仍待产生 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：', '当前整合入口：[432整合核验](verify_round432_integration.py)、[整合记录](round432_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前2014项科学检查；639份保护证据，其中本阶段编号科学文件604份', '当前2020项科学检查；642份保护证据，其中本阶段编号科学文件607份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 60.' not in texts[audit]
texts[audit] += '''
## 60. 第432轮：参与结构不能由均匀交换的初态单独补出

2026-09-24。前一轮431已完成实质内部转移结果；本轮沿其来源问题，而非继续读数优化。首先回查429—431及350、356—357，区别指定源反馈、平均场和完全均匀封闭交换。对接Björnberg、Rosengren、Ryan 2023关于全体换位和的中心元事实，不宣称新群论。

在H=ΣJᵢⱼSᵢⱼ的有限同维主体模型中，证明全部共同换轴不变关系在所有t不动，当且仅当所有可能主体对J相同。必要性由单激发子空间的非对角矩阵元直接恢复J，适用于任意N≥3、d≥2；有限域证书支持具体实现，不承担外推。均匀条件下任意初态均不能改变这些关系效果；431两种初态若补齐全部十条等强边，则整体冻结、读者概率恒1。

一般关系流由相邻权重差乘三方交换子给出。相同独立乘积ρ^{⊗N}对所有权重都不动，故429条件不选择参与表。进一步审计Jᵢⱼ(ρ)=Tr(ρSᵢⱼ)的特定直接均值规则：两个静止分量的混合却变化，在明确效果上产生1／4的仿射性违例。此处只关闭这一候选，不否定保留真实记录的反馈、量子控制器或已证明条件下的平均场近似。

6项检查、9个编号公式通过；累计2020项、607份编号科学文件、642份保护证据。新JSON只读复现、旧证据原字节保持。没有图像检查、旧实验重跑或独立代理终审声明。本轮是具体范围内的新来源约束，不是完整认知世界或GR的反证。

下一项检验把参与条件放入真实内部关系变量的成熟模型；优先核动态量子图／图性研究。必须说明变量产生、更新、与429作用合同的关系及未知输入和参考，而非仅将旧连线表改名为量子标签。有限连续交换分支继续，三维与总目标保持活动。
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
