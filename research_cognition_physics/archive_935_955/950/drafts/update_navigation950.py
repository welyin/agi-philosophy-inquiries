"""950 navigation update. Validate every edit before changing living documents."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
STAGE=HERE.parents[1]
RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths}
updates={}
heading='## 950：记录材料的真空计数与有限匹配'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf'); crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    oldheading='## 949：保持质量与协议的共同弱耦合族'
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[950报告]({pre}research_note_950.md)复用630，核实际947字面材料的4,063,232个真空味；固定局部匹配的一个圈低能余项很小，但局部参数仍可改变有限共同响应。'
           +f'[结果]({pre}950/material_vacuum_matching_results.json) · [核验]({pre}950/research_round_950_checks.json)。正式950／累计3735。'
           +f'按[投入决定]({pre}950/drafts/material_loop_decision.md)，停止字面多味材料的高阶修补，保留有效共存见证；未认证全阶父模型。接[951]({pre}951/drafts/STATUS.md)收敛最低共同验收菜单。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        s,n=re.subn(r'^\[取舍准则\][^\n]*$',lambda m:
            '[取舍准则](archive_764_/950/drafts/material_loop_decision.md)已用于实际候选选择：字面百万基本味降为备选，停止高阶修补；947有效过程和948／949接口证据保留。正式950／累计3735。接[951](archive_764_/951/drafts/STATUS.md)先核最低共同验收菜单，优先决定方向，不把候选内部问题升级为纲领门槛。',
            s,count=1,flags=re.M);assert n==1
        assert s.count('001—949轮共949份')==1
        s=s.replace('001—949轮共949份','001—950轮共950份')
        replacement='950区分实际新增味的单圈低能余项和有限局部匹配；字面百万味实现降为备选，保留947有效见证并停止其高阶修补。接[951](archive_764_/951/drafts/STATUS.md)收敛最低共同菜单，判断跨部门接合价值。完整阶段尚未完成。'
        s,n=re.subn(r'949的联动参数族保持质量谱和内部h[^。\n]*.*?完整阶段尚未完成。',lambda m:replacement,s,count=1);assert n==1
    if p in paths[1:3]:
        current='[950](archive_764_/research_note_950.md)用实际材料计数区分单圈低能余项与局部匹配：重材料不自动固定有限物理预测，基准大N也不支持直接认证全阶。按[投入决定](archive_764_/950/drafts/material_loop_decision.md)，字面百万基本味降为备选并停止其修补；947有效过程、948共同树级核、949质量保持族继续保留。接[951](archive_764_/951/drafts/STATUS.md)先收敛最低共同验收菜单，不将全部未算项串成前置任务。正式950／累计3735；整体目标未完成，A1、应用目标和定时设置保持。'
        s,n=re.subn(r'^\[949\]\(archive_764_/research_note_949\.md\)已给[^\n]*$',lambda m:current,s,count=1,flags=re.M);assert n==1
    if p==paths[3]:
        s,n=re.subn(r'^\[截至947的取舍审计\][^\n]*$',lambda m:
            '[950投入决定](950/drafts/material_loop_decision.md)将字面百万味材料降为备选，停止高阶修补；保留947有效过程与948／949连接。正式950／累计3735，950已结项。接[951](951/drafts/STATUS.md)先收敛最低共同验收菜单，原冻结入口与旧成果保持。',
            s,count=1,flags=re.M);assert n==1
    if p==paths[4]:
        s=s.replace('正式949／累计3734，949已结项','当前正式950／累计3735，950已结项；最新取舍见[950投入决定](950/drafts/material_loop_decision.md)')
    if p==paths[5]:
        assert s.count('# 231—949轮阶段成果总览')==1 and s.count('231—949的719份')==1
        s=s.replace('# 231—949轮阶段成果总览','# 231—950轮阶段成果总览').replace('231—949的719份','231—950的720份')
        s,n=re.subn(r'^\[截至947的取舍审计\][^\n]*$',lambda m:
            '[950投入决定](950/drafts/material_loop_decision.md)在948树级接合、949参数族后，只检验一次实际材料反馈。字面百万味降为备选；停止其高阶修补，保留有限共同见证。下一项先审最低共同菜单，未将候选全部自洽性升级为纲领门槛。',
            s,count=1,flags=re.M);assert n==1
    if p==paths[6]:
        s,n=re.subn(r'^\[截至947的取舍审计\][^\n]*$',lambda m:
            '[948取舍复核](948/drafts/priority_scope_reaudit.md)区分候选与纲领、有限预测与额外连续要求；[950投入决定](950/drafts/material_loop_decision.md)以630旧谱工具核实际重味，停止字面百万味高阶修补。接[951](951/drafts/STATUS.md)收敛最低共同菜单，不重做一般连续极限。',
            s,count=1,flags=re.M);assert n==1
    if p==paths[-1]:
        s,n=re.subn(r'^按\[取舍审计\][^\n]*$',lambda m:
            '按[950投入决定](../../950/drafts/material_loop_decision.md)，实际重味单圈低能余项可小，局部匹配仍影响有限预测；正式950／累计3735。字面百万味降为备选，停止其高阶修补；不从此推断纲领失败或父理论完成。接[951](../../951/drafts/STATUS.md)先收敛最低共同菜单，A1及原目标保持。',
            s,count=1,flags=re.M);assert n==1
        s=s.replace('六条共同协议：全局缺口对应与检验优先级（截至949）','六条共同协议：全局缺口对应与检验优先级（截至950）',1)
        rows={
          'C15':'|C15 表示、粒子谱与统计|1、3、4；942、947、950|950核字面多味实现的真空重数4,063,232；不能用任务初态概率替代|百万基本味降为备选；复合内部态不自动等于基本物种，未认证替代材料|',
          'C18':'|C18 耦合、质量、混合与自由参数|1、4、5、6；931、941、949—950|950同谱而不同局部K_pp可令有限共同物质响应改变约4.65%|允许自由匹配参数但须一次固定并跨部门一致；固定物理K后无此方案自由|',
          'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；630、854、895、898、944—950|950把实际重味单圈余项界接到948同核两种物质响应，固定匹配时低能余项小|单圈小余项不等于全阶或完整父匹配；停止该材料精度修补，返回共同菜单|',
          'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；929、939、942—950|947保有限协议与物理信号；950将字面百万基本味降为可替换接口|保有效共存见证，原终端权限及物理匹配仍有边界；不要求这一材料必须成功|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—950|950局部匹配仍须写入同一协变作用并共同变分，两个树级来源保互易|平直p二点核不认证完整应力；按实际采用的有限守恒／反馈结论核必要项，不先补全全部背景|'}
        for key,row in rows.items():
            s,n=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert n==1
        s,n=re.subn(r'^4\. 949排除了[^\n]*$',lambda m:
            '4. 950完成一次决定材料投入的检验：保有限有效见证，字面百万味降为备选，不继续环图、真空能和仪器修补。接[951](../../951/drafts/STATUS.md)先审最低共同菜单，只做能改变共同覆盖判断的检验。有效验收允许自由物理参数；尚未取得完整共同误差证据。',
            s,count=1,flags=re.M);assert n==1
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,raw in updates.items():p.write_bytes(raw)
print('Updated eight living navigation documents for round 950.')

