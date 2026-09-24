"""Guarded one-time navigation update for round 424; historical science stays read-only."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
ROOT = BASE.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第424轮完成' not in texts[BASE/'RESEARCH_STATE.md'], 'Already integrated'

summary = ('**第424轮完成：** [单一收缩替代全比例缩放]({prefix}research_note_424.md)在保组合重定向及完整qubit接口保留时，'
           '用一个收缩自同构替代385整套连续缩放，删除预算精确齐次性，仍得条件性三维及413方向坐标。'
           '非齐次度量给严格区别例，旧预算幂律不再沿用。9项检查，累计1969项；583份编号科学文件、618份保护证据。'
           '实际端点、连通部门、预算及单一收缩的来源仍开放，尚未无条件生成三维。')
for path in PATHS[:5]:
    prefix = 'research_cognition_physics/archive_231_/' if path == ROOT/'README.md' else ('' if path == HERE/'README.md' else 'archive_231_/')
    marker = next(line for line in texts[path].splitlines() if '**第423轮完成：**' in line)
    insertion = ('> ' if path == ROOT/'README.md' else '') + summary.format(prefix=prefix)
    texts[path] = texts[path].replace(marker, marker+'\n\n'+insertion, 1)
    texts[path] = texts[path].replace('231—423轮', '231—424轮')
    texts[path] = texts[path].replace('不能自动替代端点重放、保组合缩放和完整方向合同。',
        '不能自动替代端点重放、相容收缩和完整方向合同；424已削去全比例缩放的要求。')

direction = BASE/'research_direction.md'
texts[direction] = texts[direction].replace('最新科学轮次与检查数为423／1960', '最新科学轮次与检查数为424／1969')
next_step = ('**下一项：一个端点收缩的实际来源。** 424已明确削去全比例缩放及精确预算齐次性；接续检查实际可接续过程能否产生一个在端点等价类上良定义、保组合、可逆且反复收缩的操作。'
             '不把程序振幅缩小、概率混合或记录遗忘直接当成它。连通性、proper预算、保组合的全方向重定向及完整qubit读取仍须各自核实；'
             '已有385／415反例直接复用，不扫描收缩率、预算函数或坐标误差。像与透镜仍暂缓，三维及阶段结项未完成。')
texts[direction] = texts[direction].replace('**下一项：实际关系的位移来源。**',
    next_step+'\n\n**位置来源的总接口（423后的背景保留）：**', 1)
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新423轮及累计1960项见本文开头', '最新424轮及累计1969项见本文开头')

readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[424科学核验](archive_231_/verify_single_contraction_round.py)、[424整合核验](archive_231_/verify_round424_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计1960项；旧科学证据保持原字节', '第三阶段累计1969项；旧科学证据保持原字节')

archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成423轮', '当前完成424轮').replace('当前423不作为预定终点', '当前424不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [423：'))
texts[archive] = texts[archive].replace(row, row+'\n| [424：单一收缩替代全比例缩放](research_note_424.md) | 一个收缩自同构与极小重定向给椭球；保留条件性三维／方向坐标；非齐次预算严格区别例 | [代码](single_contraction_coordinate_audit.py)、[结果](single_contraction_coordinate_audit_results.json)、[核验](research_round_424_checks.json)；9项检查；全位移效果访问仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[424整合核验](verify_round424_integration.py)、[整合记录](round424_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace('当前1960项科学检查；615份保护证据，其中本阶段编号科学文件580份',
    '当前1969项科学检查；618份保护证据，其中本阶段编号科学文件583份')
texts[archive] = texts[archive].replace('## 下一步\n', '## 下一步\n\n'+next_step+'\n', 1)
texts[archive] += '\n424以423为科学基线，只读复算：`python -B -X utf8 research_cognition_physics/archive_231_/verify_single_contraction_round.py`。9项检查、14个编号公式；新的收缩路线不能沿用412的预算幂次校准。\n'

audit = HERE/'spatial_premise_closure_audit.md'
texts[audit] += '''
## 51. 第424轮：减少位移前提，而非再次定义位置

2026-09-24。上一目标回合完成423及冻结核验，记为实质进展。本回合先读导航、最新笔记和结果，两个独立只读候选分别回查共同读数、局部控制与已有位置商，主代理检查缩放合同。共同关系一致、配对Gram恢复、商下降及满秩端点等候选均找到221、287—288、371、376—379等既有覆盖，取消重复实验。

### 51.1 实际删去的条件

[424](research_note_424.md)在385第6节的保组合重定向路线内，删除整套连续正比例缩放和精确预算齐次性。改为一个拓扑群自同构α满足B(αg)≤qB(g)、q<1；连通的实际位置部门现在明确列入合同。proper预算给局部紧与收缩，已核的Siebert分类给幂零Lie结构；中心内外射线的预算穿壳代替连续缩放归一化，极小重定向迫使交换。紧矩阵群平均及完整射线穿壳再给椭球，由383—384取得条件性三维。

新证明需要保组合极小重定向才得到球壳，不能声称单一α替代385基础合同的所有输出。原385整体三维路线确实蕴含新合同；B(r)=r＋min(r,1)的例子满足新合同却没有任何非平凡精确预算齐次自同构，所以不是换名。

### 51.2 哪些旧结论可以接，哪些不能

412第3节的紧闭包、qubit表示延拓和稳定子论证仍成立；413的全差位移效果查询若另行给定，旧方向交会代码继续恢复坐标。实际新例使用非共形、负行列式的收缩，21次概率向量查询、零预算查询，误差小于2.39×10⁻¹⁵。

412第4节的预算幂律及倍增校准不能沿用；新例预算1、2、3对应一次、二次、四次位移，幂律会错报第三项为4。度量也未被等同最短传播长度；负行列式位移操作并不是未知qubit的普适反转，效果酉协变只用于H。重查实际效果与让旧摘要自主更新仍是不同合同。

### 51.3 交付与后继

9项新检查、14个编号公式通过，累计1969项；583份编号科学文件、618份保护证据。两个独立代理核查证明与文献，其中一个又终审实际稿件并只读复算。旧科学文件、结果和核验快照均保持原字节；旧坐标代码只复用函数，不重跑旧测试，不做图像检查。

这是一项真实的条件削减，不是认知原则已经无条件生成位移。后继集中于一个端点收缩的实际来源：操作能否下降到真实端点、保组合、可逆并反复收缩，或在明确有效尺度上是否存在足够的替代。385与415的旧失败例直接引用。端点身份、重放、连通性、预算及qubit方向的来源仍开放，不因本轮改善便触发阶段结项或完整GR完成；不再靠改变q或更换径向函数延长轮次。
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
