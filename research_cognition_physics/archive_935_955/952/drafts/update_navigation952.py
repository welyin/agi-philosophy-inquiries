"""952: update only living navigation; frozen reports remain intact."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent;STAGE=HERE.parents[1];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 952：共同交换核的局部时间来源'
def replace(s,pattern,value):
    s,n=re.subn(pattern,lambda m:value,s,count=1,flags=re.M);assert n==1,pattern
    return s
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    old='## 951：有限物体内部状态的保源表示'
    assert s.count(old)==1 and heading not in s
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[952报告]({pre}research_note_952.md)将948原记录—W／费米交换核接到同一lapse来源；端点及传播场必须共同变分，有限变化有解析余项界。'
       +f'[结果]({pre}952/portal_lapse_source_results.json) · [核验]({pre}952/research_round_952_checks.json)。正式952／累计3737。'
       +f'只签收静态二次交换与指定几何来源，不是完整Einstein或947动态概率。按[范围决定]({pre}952/drafts/shared_source_decision.md)停止静态精度优化，接[953]({pre}953/drafts/STATUS.md)核静态／有限时间过程的共同字典。\n\n')
    s=s.replace(old,block+old,1)
    if p==paths[0]:
        s=replace(s,r'^\[共同菜单v0\.4\][^\n]*$',
            '[共同菜单v0.4](archive_764_/951/drafts/joint_minimum_menu.md)保留全部目标。952已接通原共同物质核的lapse来源，明确同阶端点／传播场贡献；正式952／累计3737。接[953](archive_764_/953/drafts/STATUS.md)核静态有源态与有限时间实际准备的共同字典，不继续静态精度优化。')
        assert '001—951轮共951份' in s;s=s.replace('001—951轮共951份','001—952轮共952份',1)
        s,n=re.subn(r'951以固定物体数有效表示保留947实际操作与来源.*?完整阶段尚未完成。',lambda m:
            '951保原有效过程的材料表示；952将948共同物质核与其lapse来源接通，并排除同阶漏场项的接法。接[953](archive_764_/953/drafts/STATUS.md)核静态来源与有限时间过程的同一准备／读口，全部共同菜单保持。完整阶段尚未完成。',s,count=1);assert n==1
    if p in paths[1:3]:
        s=replace(s,r'^\[951\]\(archive_764_/research_note_951\.md\)给固定物体数[^\n]*$',
            '[952](archive_764_/research_note_952.md)将同一记录—W／费米交换接到指定lapse来源，端点和传播场必须同阶变分；有限变化有受控余项。按[范围决定](archive_764_/952/drafts/shared_source_decision.md)，停止静态精度优化，接[953](archive_764_/953/drafts/STATUS.md)核静态有源态与有限时间正过程的共同字典。原947界不自动用于新准备；全部共同菜单、A1、应用目标及定时设置保持。正式952／累计3737，整体目标未完成。')
    if p==paths[3]:
        s=replace(s,r'^\[共同菜单v0\.4\][^\n]*$',
            '[共同菜单v0.4](951/drafts/joint_minimum_menu.md)保持；952核原共同物质交换的lapse来源及同阶传播场项，正式952／累计3737，952已结项。接[953](953/drafts/STATUS.md)核静态／有限时间的同一字典，停止静态核精度优化。')
    if p==paths[4]:
        s=s.replace('当前正式951／累计3736，951已结项；最新取舍见[共同菜单v0.4](951/drafts/joint_minimum_menu.md)',
            '当前正式952／累计3737，952已结项；最新取舍见[952范围决定](952/drafts/shared_source_decision.md)')
    if p==paths[5]:
        assert '# 231—951轮阶段成果总览' in s and '231—951的721份' in s
        s=s.replace('# 231—951轮阶段成果总览','# 231—952轮阶段成果总览',1).replace('231—951的721份','231—952的722份',1)
        s=replace(s,r'^\[共同菜单v0\.4\][^\n]*$',
            '[共同菜单v0.4](951/drafts/joint_minimum_menu.md)保全部目标；952连接同一记录—物质交换的lapse来源，排除端点单独变分。停止静态精度优化，接953核静态准备与有限时间过程，未认证完整共同物理。')
    if p==paths[6]:
        s=replace(s,r'^\[948取舍复核\][^\n]*$',
            '[948取舍复核](948/drafts/priority_scope_reaudit.md)、[950投入决定](950/drafts/material_loop_decision.md)和[951共同菜单](951/drafts/joint_minimum_menu.md)保持。952用原K把记录—物质交换与lapse来源接通，明确保留阶；接[953](953/drafts/STATUS.md)核静态／有限时间字典，不继续材料编码或静态精度。')
    if p==paths[-1]:
        s=replace(s,r'^按\[951共同菜单\][^\n]*$',
            '按[共同菜单](../../951/drafts/joint_minimum_menu.md)，952核948同核交换与lapse来源，端点及传播場同阶保留；正式952／累计3737。只签收指定静态来源，接[953](../../953/drafts/STATUS.md)核有限时间共同字典；全部目标部门、A1和应用目标保持。')
        s=s.replace('传播場','传播场')
        s=s.replace('六条共同协议：全局缺口对应与检验优先级（截至951）','六条共同协议：全局缺口对应与检验优先级（截至952）',1)
        rows={
          'C14':'|C14 内部规范群、全局形式、连接|1、4、5；375、930—934、946—952|952原W／费米质量交换共用同一lapse插入及互易核|群、质量机制及源合同仍输入；二次静态部门不是完整规范动力学|',
          'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；592、854、898、944—952|952给同一静态交换核对有限lapse变化的正余项与统一解析界|不将静态有源基态误作947无源有限时间准备；该界不覆盖整个物理菜单|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—952|952明确同一记录—物质交换的lapse来源，端点与传播场同阶；独立消元／变分相符|仅固定空间几何及源合同；完整应力、Einstein及动态共同误差未签收|'}
        for k,row in rows.items():s=replace(s,r'^\|'+k+r' [^\n]*$',row)
        s=replace(s,r'^4\. 951固定[^\n]*$',
            '4. 952完成原交换核—lapse三方连接及有限余项，停止静态核精度。接[953](../../953/drafts/STATUS.md)核同一K及来源在有限时间正过程中的字典，优先保实际准备；不直接合并静态结果与旧仪器界，也不将全部高阶引力作为门槛。')
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Updated eight living navigation documents for round 952.')

