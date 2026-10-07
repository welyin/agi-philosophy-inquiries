"""951 living navigation update; preserve all frozen historical reports."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1]
RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths}
updates={}
heading='## 951：有限物体内部状态的保源表示'
def line(s,pattern,replacement):
    s,n=re.subn(pattern,lambda m:replacement,s,count=1,flags=re.M)
    assert n==1,pattern
    return s
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    oldheading='## 950：记录材料的真空计数与有限匹配'
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[951报告]({pre}research_note_951.md)将947实际H、已声明来源与完整仪器精确运输到固定物体数有效表示，原误差保持；内部态无需逐个解释为基本Dirac物种。'
        +f'[结果]({pre}951/body_sector_bridge_results.json) · [核验]({pre}951/research_round_951_checks.json)。正式951／累计3736。'
        +f'空物体扇区不添该材料圈，高能匹配仍保留；完整父物理未认证。按[共同菜单]({pre}951/drafts/joint_minimum_menu.md)停止材料表示细化，接[952]({pre}952/drafts/STATUS.md)审同一保留阶的物理来源。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        s=line(s,r'^\[取舍准则\][^\n]*$',
            '[共同菜单v0.4](archive_764_/951/drafts/joint_minimum_menu.md)保留全部阶段目标。951给保947实际H与来源的有效物体表示，免于把每个内部态列作基本物种；原有限误差不变。正式951／累计3736。接[952](archive_764_/952/drafts/STATUS.md)审较广物理的保留阶来源，停止材料与编码细化。')
        assert s.count('001—950轮共950份')==1
        s=s.replace('001—950轮共950份','001—951轮共951份')
        s,n=re.subn(r'950区分实际新增味的单圈低能余项.*?完整阶段尚未完成。',lambda m:
            '951以固定物体数有效表示保留947实际操作与来源，不必把内部态维数当基本物种数；完整共同物理仍须按[菜单v0.4](archive_764_/951/drafts/joint_minimum_menu.md)验收。接[952](archive_764_/952/drafts/STATUS.md)审同一保留阶来源，不继续物体编码细化。完整阶段尚未完成。',s,count=1);assert n==1
    if p in paths[1:3]:
        s=line(s,r'^\[950\]\(archive_764_/research_note_950\.md\)用实际材料计数[^\n]*$',
            '[951](archive_764_/research_note_951.md)给固定物体数有效表示：原947实际H、源族和仪器等距运输，误差不变；空物体扇区不添材料圈，但微观匹配仍为输入。按[共同菜单v0.4](archive_764_/951/drafts/joint_minimum_menu.md)，量子、实际3+1、引力、规范物质、记录和跨尺度全部保留为验收部门。接[952](archive_764_/952/drafts/STATUS.md)审较广物理的保留阶来源，停止材料表示细化。正式951／累计3736；完整目标未完成，A1、应用目标和定时设置保持。')
    if p==paths[3]:
        s=line(s,r'^\[950投入决定\][^\n]*$',
            '[共同菜单v0.4](951/drafts/joint_minimum_menu.md)固定全部阶段验收部门。951保947实际H与来源到有效物体扇区，停止基本多味与编码细化；完整共同物理未认证。正式951／累计3736，951已结项。接[952](952/drafts/STATUS.md)审同一保留阶的物理来源。')
    if p==paths[4]:
        s=s.replace('当前正式950／累计3735，950已结项；最新取舍见[950投入决定](950/drafts/material_loop_decision.md)',
                    '当前正式951／累计3736，951已结项；最新取舍见[共同菜单v0.4](951/drafts/joint_minimum_menu.md)')
    if p==paths[5]:
        assert s.count('# 231—950轮阶段成果总览')==1 and s.count('231—950的720份')==1
        s=s.replace('# 231—950轮阶段成果总览','# 231—951轮阶段成果总览').replace('231—950的720份','231—951的721份')
        s=line(s,r'^\[950投入决定\][^\n]*$',
            '[共同菜单v0.4](951/drafts/joint_minimum_menu.md)保留全部目标部门。951提供保源的固定物体数有效表示，原947有限证据可继续使用，不必修完950基本多味候选。停止编码与材料表示细化，接952审共同保留阶物理来源。')
    if p==paths[6]:
        s=line(s,r'^\[948取舍复核\][^\n]*$',
            '[948取舍复核](948/drafts/priority_scope_reaudit.md)和[950投入决定](950/drafts/material_loop_decision.md)区分候选问题与纲领门槛；[951共同菜单](951/drafts/joint_minimum_menu.md)保留全部部门，并给固定物体数表示的实际H／来源运输。接[952](952/drafts/STATUS.md)审物理保留阶，不继续材料编码。')
    if p==paths[-1]:
        s=line(s,r'^按\[950投入决定\][^\n]*$',
            '按[951共同菜单](../../951/drafts/joint_minimum_menu.md)，947实际H、来源和仪器可精确运输到固定物体数有效表示，原误差保持；正式951／累计3736。内部态不必逐个作为基本味，微观匹配仍开放。接[952](../../952/drafts/STATUS.md)审较广物理的保留阶来源；不删除任何目标部门，不继续材料编码，A1与原目标保持。')
        s=s.replace('六条共同协议：全局缺口对应与检验优先级（截至950）','六条共同协议：全局缺口对应与检验优先级（截至951）',1)
        rows={
            'C01':'|C01 状态、概率、仪器、复合|1、3、4；001—230、923、947、951|951以精确扇区等距保947完整cq仪器、未知输入与旧参考，原界保持|表示身份不制造全部仪器，也不认证完整父物理的有限匹配|',
            'C15':'|C15 表示、粒子谱与统计|1、3、4；942、947、950—951|951固定物体数有效表示不把内部态数作为基本Dirac真空种数|守物体数工作域、稳定材料与高能匹配仍为输入；不是SM束缚形成证明|',
            'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；592、854、898、944—951|951在同一有效对象内同时保实际H与源族，947界原样运输|精确重表示不等于微观到有效匹配；共同菜单中各部门的有限误差仍须同域成立|',
            'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；929、939、942—951|951保原完整有限协议于有效物体内部态，避免强制基本多味解释|稳定物体数、准备与末读仍输入；停止编码／器件细化，不自动签收真实材料|',
            'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—951|951源无关扇区等距同时保已给源，空物体不抹去有物体时的交换|没有给出的完整应力不能由表示身份领取；下一项核较广物理的共同保留阶|'}
        for key,row in rows.items():s=line(s,r'^\|'+key+r' [^\n]*$',row)
        s=line(s,r'^4\. 950完成[^\n]*$',
            '4. 951固定共同验收菜单并给保源的有效物体表示，原947界保持，基本多味不是必需实现。停止编码与材料细化，接[952](../../952/drafts/STATUS.md)审较广物理的同阶来源；不把微观制造或全部UV变成先决条件，也不删去3+1、规范或引力部门。')
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,raw in updates.items():p.write_bytes(raw)
print('Updated eight living navigation documents for round 951.')

