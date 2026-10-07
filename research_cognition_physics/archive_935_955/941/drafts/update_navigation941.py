"""Update only living navigation after the 941 scientific calculation."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 941：组织更新、复合钟与质量来源的共同扩展'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    oldheading='## 940：同一物质的弱场交换与正量子联合接口'
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[941报告]({pre}research_note_941.md)将925实际内部更新接入复合物体质量壳，同一h给钟速、质量与非对易来源；漏掉相互作用能使红移传递概率差约0.00235829。'
        +f'[结果]({pre}941/composite_operation_bridge_results.json) · [核验]({pre}941/research_round_941_checks.json)。正式941／累计3726。'
        +'等效耦合和复合材料匹配为物理输入，原场材料实现及动态几何未签收。停止钟精度和器件优化，'
        +f'接[942]({pre}942/drafts/STATUS.md)把明确操作材料放入一份共同物理对象。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        assert s.count('001—940轮共940份')==1
        s=s.replace('001—940轮共940份','001—941轮共941份')
        old='939保留整体恢复表和关系读口；940把实际Dirac来源接入弱场传播与约束的共同正量子接口，不能直接迁入原E_rec或全标准模型。停止该示例精度／高阶细化，接[941](archive_764_/941/drafts/STATUS.md)先比较整体共同对象的覆盖与净输入。完整阶段尚未完成。'
        new='939整体恢复表与940弱场接口保持原范围；941用925实际更新接通复合组织的钟速、质量和来源。等效耦合与材料匹配为输入，不能直接签收原场论；按[共同方案](archive_764_/941/drafts/common_model_with_material_roles.md)，接[942](archive_764_/942/drafts/STATUS.md)选择同一物理对象及操作材料。停止钟与弱场示例细化，完整阶段尚未完成。'
        assert old in s;s=s.replace(old,new,1)
    if p in paths[1:3]:
        new='941的[共同扩展](archive_764_/research_note_941.md)把925内部更新、操作钟速和质量来源连到同一h；复合物体及等效耦合是明示物理输入。漏相互作用源已有有限记录差，精确重组可共同保响应。停止钟精度和器件优化，按[整体方案](archive_764_/941/drafts/common_model_with_material_roles.md)接[942](archive_764_/942/drafts/STATUS.md)将操作材料放入一份共同物理对象；不把从原物种微观制造全部装置设为普遍门槛，也不把尚未匹配的材料称作原物种。'
        s,count=re.subn(r'^939的\[整体恢复表\][^\n]*$',lambda m:new,s,count=1,flags=re.M);assert count==1
    if p==paths[5]:
        assert s.count('# 231—940轮阶段成果总览')==1 and s.count('231—940的710份')==1
        s=s.replace('# 231—940轮阶段成果总览','# 231—941轮阶段成果总览').replace('231—940的710份','231—941的711份')
    if p==paths[-1]:
        assert s.count('六条共同协议：全局缺口对应与检验优先级（截至940）')==1
        s=s.replace('六条共同协议：全局缺口对应与检验优先级（截至940）','六条共同协议：全局缺口对应与检验优先级（截至941）')
        rows={
          'C09':'|C09 钟尺、参考态与可访问性|1、3、6；522—530、861、934—941|941复合材料合同中同一h给钟速与质量，925未知记录红移运输可行|复合匹配和等效耦合为输入；未构造原材料末读及boost仪器，停止钟优化|',
          'C10':'|C10 共同传播几何与普适耦合|1、4、5；334—340、924、930、940—941|940守恒来源交换保传播与约束；941将操作、惯性及引力内能统一为M=mI+h|共同锥、量子等效耦合仍为物理输入；不能只比能源均值|',
          'C18':'|C18 耦合、质量、混合与自由参数|1、4、5、6；930—941|941复合扩展保完整相互作用质量；非对易来源须用同一谱函数微分|原物理材料与参数字典尚需选定；自由参数不妨碍阶段存在性|',
          'C20':'|C20 粗化、有效极限与重整化|1、3、4、6、H3；854、898、936、939—941|941精确闭合重组保质量壳及来源，固定动量共同钟速有任意输入误差界|该界只比两份有效模型，不代替原场到复合材料匹配；940有限耦合外推仍未认证|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；873、935—941|941同一M给钟速及来源，漏相互作用产生有限记录差；940物理源结果保留|完整动态几何与原材料来源尚未在一份对象里共同验收；不能把局部扩展当统一完成|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert count==1
        old='4. 939[共同恢复表](../../939/drafts/common_model_recovery_map.md)与940弱场正量子接口分别保原身份；940给实际物理来源而非任意几何惯量。有限占据外推未认证，停止该示例精度／高阶工作。接[941](../../941/drafts/STATUS.md)先比较整体共同对象和认知操作的覆盖，原生H不是普遍前提；不重开尘埃、参考及旧反馈修补链。'
        new='4. 939恢复表、940弱场实例和941复合组织扩展保持各自范围；按[共同方案](../../941/drafts/common_model_with_material_roles.md)接[942](../../942/drafts/STATUS.md)把操作材料放入一份共同物理对象。停止钟、参考、占据及高阶修补；允许明确有效材料输入，但不得声称已由原物种产生。'
        assert old in s;s=s.replace(old,new,1)
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,raw in updates.items():p.write_bytes(raw)
print('Updated eight living navigation documents; preserved historical science.')
