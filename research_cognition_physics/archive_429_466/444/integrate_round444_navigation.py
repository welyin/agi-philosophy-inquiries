"""Guarded navigation update for the reciprocal subject-port realization."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE, ROOT = HERE.parent, HERE.parent.parent
PATHS = [ROOT/'README.md', BASE/'README.md', BASE/'research_direction.md',
         BASE/'RESEARCH_STATE.md', HERE/'README.md', HERE/'spatial_premise_closure_audit.md']
original = {p: p.read_bytes() for p in PATHS}
texts = {p: raw.decode('utf-8').replace('\r\n', '\n') for p, raw in original.items()}
assert '**第444轮完成：**' not in texts[BASE/'RESEARCH_STATE.md']
summary = ('**第444轮完成：** [伙伴地址载体与四主体作用边界](%sresearch_note_444.md)'
           '给主体内部互指端口的显式张量实现，精确承接440自主H及主体读数，'
           '覆盖编码内未知图、数据和参考。固定容量时地址预算O(N log N)，'
           '不再要求每个潜在主体对一份独立边qubit；四主体作用及访问架构仍为输入。'
           '证明固定逐度代码中至多三主体的直接作用不能改变图。'
           '6项检查，累计2095项；643份编号科学文件、678份保护证据。'
           '没有生成三维或把一致性代码升级为认知公理。')
for path in PATHS[:5]:
    prefix = ('research_cognition_physics/archive_231_/' if path == ROOT/'README.md'
              else ('' if path == HERE/'README.md' else 'archive_231_/'))
    old = next(line for line in texts[path].splitlines() if '**第443轮完成：**' in line)
    texts[path] = texts[path].replace(old, old+'\n\n'+('> ' if path == ROOT/'README.md' else '')+summary % prefix, 1)
    texts[path] = texts[path].replace('231—443轮', '231—444轮')
next_item = (
    '**下一项：直接交换怎样维持伙伴一致性。** 444已给完整端口实现，'
    '但四主体换边规则仍未由429推出。不能把每一时刻严格留在固定度代码新增为认知公理；'
    '检验地址SWAP在短暂不一致中间态后形成的有效作用，是否只产生一般重配对，'
    '还是能选择440沿已有三路径的flip。须写出实际生成元、未知输入范围、约束能及访问代价。'
    '434—436的一般静态模拟存在性已足够，不再为重复证明可模拟性或提高精度另开轮。'
    '若一致性只允许跨分量重配对，准确保留几何局域选择的独立缺口；'
    '初始关系、测距、三维及新主体接入仍开放。')
for path in (BASE/'research_direction.md', BASE/'RESEARCH_STATE.md', HERE/'README.md'):
    old = next(line for line in texts[path].splitlines() if line.startswith('**下一项：可变伙伴的内部关系载体。**'))
    texts[path] = texts[path].replace(old, next_item+
        '\n\n**443后载体任务（444已补显式实现，作用来源仍开放）：** '+old.split('** ', 1)[1], 1)
texts[BASE/'research_direction.md'] = texts[BASE/'research_direction.md'].replace(
    '最新科学轮次与检查数为443／2089', '最新科学轮次与检查数为444／2095')
state = BASE/'RESEARCH_STATE.md'
texts[state] = texts[state].replace('最新443轮及累计2089项见本文开头', '最新444轮及累计2095项见本文开头')
readme = BASE/'README.md'
texts[readme] = texts[readme].replace('当前复算与冻结入口：',
    '当前复算与冻结入口：[444科学核验](archive_231_/verify_partner_port_carrier_round.py)、[444整合核验](archive_231_/verify_round444_integration.py)；历史入口：', 1)
texts[readme] = texts[readme].replace('第三阶段累计2089项；旧科学证据保持原字节', '第三阶段累计2095项；旧科学证据保持原字节')
archive = HERE/'README.md'
texts[archive] = texts[archive].replace('当前完成443轮', '当前完成444轮').replace('当前443不作为预定终点', '当前444不作为预定终点')
row = next(line for line in texts[archive].splitlines() if line.startswith('| [443：'))
texts[archive] = texts[archive].replace(row, row+
    '\n| [444：伙伴地址载体与四主体作用边界](research_note_444.md) | 显式主体因子；全未知代码交织；合法四主体扩展；固定度低体边界 | [代码](partner_port_carrier_audit.py)、[结果](partner_port_carrier_audit_results.json)、[核验](research_round_444_checks.json)；6项检查；地址空间和四主体规则仍为输入 |', 1)
texts[archive] = texts[archive].replace('当前整合入口：',
    '当前整合入口：[444整合核验](verify_round444_integration.py)、[整合记录](round444_integration_checks.json)；历史入口：', 1)
texts[archive] = texts[archive].replace(
    '当前2089项科学检查；675份保护证据，其中本阶段编号科学文件640份',
    '当前2095项科学检查；678份保护证据，其中本阶段编号科学文件643份')
audit = HERE/'spatial_premise_closure_audit.md'
assert '## 72.' not in texts[audit]
texts[audit] += '''
## 72. 第444轮：内部伙伴地址的精确承载与作用来源边界

2026-09-24。继续443明确的伙伴载体问题。438已给带标签匹配的N log N记忆量级，257／269已区分地址与通道，434—436已给一般交换模拟；本轮不重算这些旧结论。对接Arrighi等Addressable Quantum Gates（2023）的内部量子地址思想，但不把其运输原语当作429自动推论。

固定N、容量C，每主体保留原数据及C个(NC＋1)维地址寄存器。合法记录为无自主体、无重边的互指端口匹配。每图的全部端口赋值数为∏C!／(C−d_i)!，均匀相干赋值给等距W。主体数据因子不变，固定容量时地址预算NC ceil(log2(NC＋1))为O(N log N)，省去独立潜在边qubit；小模型可能反而较大，且地址维数随N增长。

对已有三路径a-b-c-d，把ab和cd的端口配对同时改成ac和bd，只读取四主体。源端外边必须各恰有一对互指，新边为空、bc存在；这一局部双射在源、像两端保持全局合法性等价，故加伴随后不漏出合法部门。若省略唯一性门，非法ab重边可被修成合法态，伴随就会泄漏；代码给明确反例。没有用全局代码投影伪称四主体支撑。

每条图路径在等大小赋值集合间双射，保留多路径权重，得到F_port W＝WF、完整H和逐主体数据读数交织；未知图／数据／参考均覆盖在W代码范围内。4主体41图对应441地址态、384计权转换；六树对应17496地址态、69984转换。512维完整数据模型整数交织，有限时误差约4.05×10^−15。任意非均匀地址态不自动等价于旧图，全部旧关系操作的低体访问也未证明。

固定逐顶点度数且精确保代码时，至多三主体作用不能改变图：外部互指记录锁定跨界边，三个内部顶点的固定内度唯一决定简单图。四主体flip达到下界。349920份外部记录检查通过；该限制不排除非一致中间态、虚过程、变度规则或新增辅助。代码约束不是新认知公理。

实际消去的是每潜在对独立关系载体的要求，未消去四主体潜在作用表、地址访问、代码准备及其来源。下一项检验直接地址SWAP与内部一致性是否选择现有路径flip，或仅产生一般跨分量重配对；不重述436通用模拟存在性。6项检查、12公式，累计2095项、643份编号科学文件、678份保护证据；三维与GR继续开放，阶段未结项。
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
