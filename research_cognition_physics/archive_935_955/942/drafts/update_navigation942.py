"""Update living navigation; keep all frozen rounds unchanged."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 942：明确Dirac材料中的内部更新与几何来源'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    oldheading='## 941：组织更新、复合钟与质量来源的共同扩展'
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[942报告]({pre}research_note_942.md)以明确新增的中性Dirac多重态承载925实际更新，同一作用给传播及应力；正能任务、标架来源和实际质量交换核通过。'
        +f'[结果]({pre}942/dirac_material_embedding_results.json) · [核验]({pre}942/research_round_942_checks.json)。正式942／累计3727。'
        +'物种与质量混合为输入，内部角色不自动成为空间主体；929通用编码需4,063,232个分量，局域读取未签收。停止本编码细化，'
        +f'接[943]({pre}943/drafts/STATUS.md)先判断共同模型的实际组织与访问划分。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        assert s.count('001—941轮共941份')==1
        s=s.replace('001—941轮共941份','001—942轮共942份')
        old='939整体恢复表与940弱场接口保持原范围；941用925实际更新接通复合组织的钟速、质量和来源。等效耦合与材料匹配为输入，不能直接签收原场论；按[共同方案](archive_764_/941/drafts/common_model_with_material_roles.md)，接[942](archive_764_/942/drafts/STATUS.md)选择同一物理对象及操作材料。停止钟与弱场示例细化，完整阶段尚未完成。'
        new='939整体恢复表及940—941接口保持；942用新增中性Dirac材料实现内部更新及同一应力，但物种成本和实际访问仍为独立条件。停止通用编码细化，接[943](archive_764_/943/drafts/STATUS.md)优先核共同模型的真实组织／访问划分，不能把内部角色等同空间主体。完整阶段尚未完成。'
        assert old in s;s=s.replace(old,new,1)
    if p in paths[1:3]:
        new='[942](archive_764_/research_note_942.md)已给新增中性Dirac材料的明确作用：原内部h进入质量混合，同一作用产生更新、传播和应力，来源运输复用720。它改变物种与父作用；直接运输929需4,063,232个分量，内部角色尚非空间主体，数学效果未自动成为局部仪器。因此保留为存在性工具，停止物种／钟／读口细化，接[943](archive_764_/943/drafts/STATUS.md)优先检验实际组织与访问划分；不以本候选全部修复为阶段门槛。'
        s,count=re.subn(r'^941的\[共同扩展\][^\n]*$',lambda m:new,s,count=1,flags=re.M);assert count==1
    if p==paths[5]:
        assert s.count('# 231—941轮阶段成果总览')==1 and s.count('231—941的711份')==1
        s=s.replace('# 231—941轮阶段成果总览','# 231—942轮阶段成果总览').replace('231—941的711份','231—942的712份')
    if p==paths[-1]:
        assert s.count('六条共同协议：全局缺口对应与检验优先级（截至941）')==1
        s=s.replace('六条共同协议：全局缺口对应与检验优先级（截至941）','六条共同协议：全局缺口对应与检验优先级（截至942）')
        rows={
          'C02':'|C02 子系统、约束、边界拼接|1、3、4、5；360—365、705—706、942|942内部任务可嵌入新物理场，逻辑张量分工不等于空间分离的参与者|优先核实际访问划分和共同任务；不以内部编码替代物理主体|',
          'C04':'|C04 同一内部动力学与控制|1、6、H2；935—942|942新增Dirac质量混合承载原有限h，同一作用给物理来源|工程化质量及新物种仍输入；局部作用不自动实现全部操作权限|',
          'C15':'|C15 表示、粒子谱与统计|1、3、4；既有旋量与CAR、934、942|942明确新增中性矢量样Dirac分量；未冒充原标准模型物种|实际物种成本、物理谱和匹配单列；百万分量通用编码只留存在性工具|',
          'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；853、929、939、941—942|942明确材料可承载内部协议，但正能效果一般随动量变，实际读口未实现|先判定共同组织／访问，而非继续设计本编码所有仪器|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；873、935—942|942同一新作用给应力及质量插入；720标架项不可略，实际质量通过940交换身份|新增量子材料改变来源和匹配，旧曲背景量子结论不能免费搬入；全共同反馈未签收|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert count==1
        old='4. 939恢复表、940弱场实例和941复合组织扩展保持各自范围；按[共同方案](../../941/drafts/common_model_with_material_roles.md)接[942](../../942/drafts/STATUS.md)把操作材料放入一份共同物理对象。停止钟、参考、占据及高阶修补；允许明确有效材料输入，但不得声称已由原物种产生。'
        new='4. 942新Dirac材料接通更新与同一物理来源，通用编码成本及空间访问边界明确。保留此存在性工具，停止物种／仪器优化；接[943](../../943/drafts/STATUS.md)优先核共同对象中实际参与者与访问划分。939恢复表及940—941接口保持原范围，不拼接成已完成模型。'
        assert old in s;s=s.replace(old,new,1)
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,raw in updates.items():p.write_bytes(raw)
print('Updated eight living navigation documents; preserved frozen science.')
