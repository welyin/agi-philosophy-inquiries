"""Guarded navigation update for round 449 protected processing and common readout window."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第449轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = '**第449轮完成：** [记录保护、伙伴重组与条件处理的兼容定理](%sresearch_note_449.md)在448共同载体同时保留完整记录势与较弱条件处理，给全未知匹配／逻辑／参考的残差证书及有限N全时间漏码界。四主体同一读数窗口的原伙伴消息差>1／16、两消息合法重配概率均>3／8；纯运输对该原伙伴的理想信号全时为零。7项检查，累计2127项；658份编号科学文件、693份保护证据。处理按g²／Δ变慢，作用、接触、准备及读数仍为输入，未生成三维。'
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第448轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—448轮', '231—449轮')
next_item = '**下一项：组合主体的交互合同与记录作用来源。** 449已证明同载体记录保护、伙伴重组与条件处理兼容，并核算强度／时间代价，不继续调整ε或扩增同类模型。先审计429“同态不动”从原始主体到复合主体应保留的条件，再核交换关系自由度能否承担互指比较；复用430、432、434—436已证的关系动力学与一般模拟。新增轮次必须给组合相容条件、来源障碍或减少独立输入，不能再次把另加对角势当来源。共同身份字典、参与规则、实际定位、三维与GR仍开放。'
for path in (BASE/'research_direction.md', BASE/'RESEARCH_STATE.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：记录保护与条件处理的兼容条件及作用来源。**'))
    texts[path] = texts[path].replace(old, next_item+
        '\n\n**448后兼容任务（449已完成同窗记录保护、重组与条件读数）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace(
    '最新科学轮次与检查数为448／2120', '最新科学轮次与检查数为449／2127')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新448轮及累计2120项见本文开头', '最新449轮及累计2127项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[449科学核验](archive_231_/verify_protected_partner_processing_round.py)、[449整合核验](archive_231_/verify_round449_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2120项；旧科学证据保持原字节', '第三阶段累计2127项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成448轮', '当前完成449轮').replace('当前448不作为预定终点', '当前449不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [448：'))
texts[archive] = texts[archive].replace(row, row+
    '\n| [449：记录保护、伙伴重组与条件处理的兼容定理](research_note_449.md) | 奖励分类；全未知误差；同窗条件信号与合法重配；能量时间代价 | [代码](protected_partner_processing_audit.py)、[结果](protected_partner_processing_audit_results.json)、[核验](research_round_449_checks.json)；7项检查；两种势及尺度仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[449整合核验](verify_round449_integration.py)、[整合记录](round449_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace(
    '当前2120项科学检查；690份保护证据，其中本阶段编号科学文件655份',
    '当前2127项科学检查；693份保护证据，其中本阶段编号科学文件658份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 77.' not in texts[audit]
texts[audit] += '\n## 77. 第449轮：受保护伙伴重组与条件处理的同一窗口\n\n2026-09-24。接续448相容任务，不重做434—436模拟存在性。对称pair对角奖励q=diag(a,b,b,c)分解为αI＋β(Z_i＋Z_j)＋γZ_iZ_j。共同内部翻转只要求a=c，独立翻转或每匹配全部逻辑同能则要求a=b=c。正奖励与非零ZZ可共存，但正奖励自身不保证全部逻辑／非匹配谱带分开。未知信息保存不被偷换为能量简并。\n\n在448同一完整原始载体取H=Δ[h_R＋εG＋ε²(jh_A＋ℓL)]，τ=ε²Δt。D保持匹配子空间；W=Z＋εW1＋ε²W2用445虚交换映射，精确残差为ε³(GW2＋DW1−W1B)＋ε⁴(DW2−W2B)，B=K＋D_P。均匀L与Wk交织，残差ℓ项完全相消；这不免除ℓ能量及读数代价。四主体裸码全未知／参考误差为(2＋|τ|(4＋4|j|))|ε|＋[13／5＋|τ|(4＋26|j|／5)]ε²。完整非因子化目标包含伙伴／逻辑纠缠。\n\n能量守恒另给全时间、任意有限偶数N漏码幅度≤m|ε|＋N(|j|／2＋|ℓ|)ε²，m=N(N−1)／2。不是无限规模一致、自动纠错或全时间有效生成元逼近。\n\nN4、j=ℓ=1、ε1／1024，初始M0=(01)(23)，源0的0／1消息、其余0，实际原伙伴1始终读A空值。48维两列40阶Gaussian有理Taylor、||B||≤7及exp7<2187给τ1有效概率差>9／100、两消息合法重配率>39／100。导数界6／2与完整误差合成τ∈[.999,1.001]上实际信号>1／16、两合法重配率>3／8，非互指部门不计成功且无后选。j0精确因子化时，源标签1只能到π(1)=1的非互指主体，故理想匹配内原伙伴信号全时为0，实际至多36ε²；本轮信号不能归为纯搬运。\n\n代价：J_A=jg²／Δ、λ=ℓg²／Δ、t=τΔ／g²；保护增强时所证家族变慢。g1、Δ1024的物理窗口[1022.976,1025.024]，原始8份5维寄存器、六个潜在主体对及两种作用均列明，代码／消息准备和实际读出仍为输入。429若重新施加于完整闭合二逻辑位，会排除非恒定对角势；因此不能将其直接称为429推论，也不强加所有复合层级均满足同态完全不动。\n\n7项检查、14公式，累计2127项；658份编号科学文件、693份保护证据。独立只读审核，旧科学文件保持原字节，无图像。后继只核组合交互合同及比较记录来源，复用已知一般交换模拟；不继续调ε或加任意势。共同身份、接触选择、实际定位、三维与GR继续开放，阶段未结项。\n'
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
