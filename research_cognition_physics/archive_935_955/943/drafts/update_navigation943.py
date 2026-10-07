"""Update eight living navigation documents; retain all frozen evidence."""
from pathlib import Path
import re
STAGE = Path(__file__).resolve().parents[2]
RESEARCH = STAGE.parent
paths = [RESEARCH/n for n in ('README.md', 'research_direction.md', 'RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md', '文件索引.md', '阶段成果总览.md', '跨阶段主题索引.md', '_shared/notes/unified_physics_condition_ledger_current.md')]
original = {p:p.read_bytes() for p in paths}
updates = {}
heading = '## 943：内部记录能否进入共同物理交互'
for p, raw in original.items():
    bom = raw.startswith(b'\xef\xbb\xbf')
    crlf = b'\r\n' in raw
    s = raw.decode('utf-8-sig').replace('\r\n', '\n')
    oldheading = '## 942：明确Dirac材料中的内部更新与几何来源'
    assert heading not in s and s.count(oldheading) == 1, p
    pre = 'archive_764_/' if p.parent == RESEARCH else '../../' if p == paths[-1] else ''
    block = (heading+'\n\n'+f'[943报告]({pre}research_note_943.md)复用935谱与939访问引理，核942实际材料：无额外χ味参考时，仅质量／几何耦合不能把929逻辑内容传给中性探针。内部事件差0.60，指定任务最坏误差至少0.30。'
        +f'[结果]({pre}943/material_content_access_results.json) · [核验]({pre}943/research_round_943_checks.json)。正式943／累计3728。'
        +'新增与质量对易的味矩阵可给有限内容信号，但属于额外物理输入；原场读出未签收。保留材料工具，停止编码细化，'
        +f'接[944]({pre}944/drafts/STATUS.md)选择共同交互；本候选访问边界不是纲领失败。\n\n')
    s = s.replace(oldheading, block+oldheading, 1)
    if p == paths[0]:
        assert s.count('001—942轮共942份') == 1
        s = s.replace('001—942轮共942份', '001—943轮共943份')
        old = '939整体恢复表及940—941接口保持；942用新增中性Dirac材料实现内部更新及同一应力，但物种成本和实际访问仍为独立条件。停止通用编码细化，接[943](archive_764_/943/drafts/STATUS.md)优先核共同模型的真实组织／访问划分，不能把内部角色等同空间主体。完整阶段尚未完成。'
        new = '939整体恢复表及940—942接口保持原范围；943判定仅质量／几何耦合不能自动输出929逻辑内容，明确新交互可解除该访问限制。保留材料工具、停止编码细化，接[944](archive_764_/944/drafts/STATUS.md)选择同一记录—传播—来源的有效交互，不把修好某候选作为整个纲领门槛。完整阶段尚未完成。'
        assert old in s
        s = s.replace(old, new, 1)
    if p in paths[1:3]:
        new = '[943](archive_764_/research_note_943.md)已判定942材料的内容访问边界：在无额外χ味参考、运动准备固定的合同中，仅质量／几何耦合给929逻辑输入恒定的中性探针输出；指定事件最坏误差至少0.30。新味矩阵可与质量对易并给有限指针区别，但它是新增交互输入，原场读出未签收。保留942工具，停止编码与仪器细化，接[944](archive_764_/944/drafts/STATUS.md)选择共同交互；不得将这个候选的访问限制升级成整个认知纲领失败。'
        s, count = re.subn(r'^\[942\]\(archive_764_/research_note_942\.md\)已给新增中性Dirac材料[^\n]*$', lambda m:new, s, count=1, flags=re.M)
        assert count == 1, p
    if p == paths[5]:
        assert s.count('# 231—942轮阶段成果总览') == 1 and s.count('231—942的712份') == 1
        s = s.replace('# 231—942轮阶段成果总览', '# 231—943轮阶段成果总览').replace('231—942的712份', '231—943的713份')
    if p == paths[-1]:
        old = '六条共同协议：全局缺口对应与检验优先级（截至942）'
        assert s.count(old) == 1
        s = s.replace(old, '六条共同协议：全局缺口对应与检验优先级（截至943）')
        rows = {
            'C02':'|C02 子系统、约束、边界拼接|1、3、4、5；360—365、705—706、942—943|942内部角色不等于空间主体；943给指定中性访问合同下的恒定逻辑输出|实际访问依赖所选交互／参考；不另加全部参与者空间分离的要求|',
            'C04':'|C04 同一内部动力学与控制|1、6、H2；935—943|942质量混合可承载h，943表明这不能自动提供跨材料逻辑通信|先选择共同内容敏感交互；新矩阵是物理输入，不继续编译器优化|',
            'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；853、929、939、942—943|仅质量／几何耦合在943合同下不输出逻辑内容；新保自由质量的味矩阵有有限正向见证|原p—R实际接收尚未运输，先判定共同交互收益；本候选边界不等于纲领失败|'
        }
        for key, row in rows.items():
            s, count = re.subn(r'^\|'+key+r' [^\n]*$', lambda m:row, s, count=1, flags=re.M)
            assert count == 1
        old = '4. 942新Dirac材料接通更新与同一物理来源，通用编码成本及空间访问边界明确。保留此存在性工具，停止物种／仪器优化；接[943](../../943/drafts/STATUS.md)优先核共同对象中实际参与者与访问划分。939恢复表及940—941接口保持原范围，不拼接成已完成模型。'
        new = '4. 943已完成材料内容访问的决定性检验：942的纯质量／几何接口不能自动输出929逻辑记录，新增味交互可给有限区别。保留存在性工具，停止编码／仪器修补；接[944](../../944/drafts/STATUS.md)选择共同记录—传播—来源的有效交互。939整体恢复表及940—942接口不拼接成已完成模型，也不以修复此候选作为纲领前提。'
        assert old in s
        s = s.replace(old, new, 1)
    if crlf:
        s = s.replace('\n', '\r\n')
    updates[p] = (b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes() == raw for p, raw in original.items())
for p, raw in updates.items():
    p.write_bytes(raw)
print('Updated eight living navigation documents; preserved frozen science.')
