"""One-time guarded update of the eight living navigation documents."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
heading='## 937：正能源内部参考与原全场约束共同接入'
oldheading='## 936：原生有效生成元共同保来源与几何力'
original={p:p.read_bytes() for p in paths}
updates={}
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    assert heading not in s and s.count(oldheading)==1,p
    prefix='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[937报告]({prefix}research_note_937.md)采用成熟正能源材料参考，把内部钟尺精确接入原753非Abelian—Einstein初值；原物质资料与Gauss保持，参考应力及质量来源由同一约化根产生。'
       +f'[结果]({prefix}937/dust_common_model_results.json) · [核验]({prefix}937/research_round_937_checks.json)。正式937／累计3722。'
       +'新增参考介质是物理输入；经典共同骨架未签收全量子Ward／有效匹配，873原h钟结果仍成立。停止尘埃和精确钟细化，'
       +f'接[938]({prefix}938/drafts/STATUS.md)直接审共同父作用的有限量子恢复合同；不把候选内部问题升级为纲领门槛。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        assert s.count('001—936轮共936份')==1
        s=s.replace('001—936轮共936份','001—937轮共937份')
        old='936给正规有限模式符号的原生正过程、同源力与可计算占据尾界；其校准不是完整物理模型。停止一般量子化和校准细化，接[937](archive_764_/937/drafts/STATUS.md)填实际全物理符号及匹配。'
        new='936给原生有效量子化接口；937以新增正能源参考接通原全场经典约束、应力和质量来源。两者尚未组成全物种正量子过程。停止钟／图册优化，接[938](archive_764_/938/drafts/STATUS.md)直接审有限量子恢复合同。'
        assert old in s;s=s.replace(old,new,1)
    if p in paths[1:3]:
        old='935已停止钟／编译器扩展；936完成原生量子化的来源与几何力构造接口，并明确完整物理符号仍缺。937直接填入共同物理内容及成熟低能匹配，不继续校准精度；不要求所有旧候选成功。应用目标和定时设置保持。'
        new='935—936已停止编译器和通用量子化校准；937采用新增正能源参考，实际接通原全场经典约束、参考应力与质量源。它尚非全量子共同模型。938直接审有限量子恢复合同，不继续尘埃图册／精确钟细化；不要求所有旧候选成功。应用目标和定时设置保持。'
        assert old in s;s=s.replace(old,new,1)
    if p==paths[5]:
        assert s.count('# 231—936轮阶段成果总览')==1
        assert s.count('231—936的706份')==1
        s=s.replace('# 231—936轮阶段成果总览','# 231—937轮阶段成果总览').replace('231—936的706份','231—937的707份')
    if p==paths[-1]:
        assert s.count('六条共同协议：全局缺口对应与检验优先级（截至936）')==1
        s=s.replace('六条共同协议：全局缺口对应与检验优先级（截至936）','六条共同协议：全局缺口对应与检验优先级（截至937）')
        rows={
        'C04':'|C04 同一内部动力学与控制|1、6、H2；935—937|936给有限物理符号量子化接口；937给新增尘埃参考的全场经典根，两者尚未共同量子实现|直接匹配实际父作用的有限预测，不继续通用编译器或任意符号校准|',
        'C09':'|C09 钟尺、参考态与可访问性|1、3、6；522—530、861、934—937|937正能源尘埃钟尺与原753初值共同成立；新增介质是物理输入，实际量子仪器仍未覆盖|选定参考及有限任务；不把原h/Y或低密度图册的全部问题当共同前提|',
        'C11':'|C11 引力自由度及约束闭合|1、4、5；344、351—352、366、764、861、936—937|937四约束约化及原全场接入，经典Gauss保持；全费米量子约化／Ward未交付|在同一父作用上核要采信的量子任务；经典根不能自动成为完整量子符号|',
        'C18':'|C18 耦合、质量、混合与自由参数|1、4、5、6；930—931、935—937|937新参考的质量源明确；静止尘埃四腿系数与873原h钟不同，不能混用父作用|允许参数及参考实现选择，保同一来源；量子排序和有限匹配仍须交代|',
        'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；585—590、873、935—937|937参考正应力、几何根变分与原约束同源；全物种量子反作用仍未签收|先核有限量子恢复合同，不以全部无截断量子引力完成为共同门槛|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert count==1
        old='4. 936给正规有限模式符号的原生正过程、来源、canonical力及占据尾界，实际量子化字典是实现选择，不另计新认知推导。完整物理符号及其匹配未完成；接[937](../../937/drafts/STATUS.md)直接交付物理填充，停止校准精度和编译器修补。'
        new='4. 936给原生有效量子化接口；937以新增正能源参考接通原全场经典约束及来源。参考物种、作用与维数仍是物理输入；未完成共同量子模型。接[938](../../938/drafts/STATUS.md)直接审有限量子恢复合同，停止尘埃／钟和通用量子化校准细化。'
        assert old in s;s=s.replace(old,new,1)
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
# Verify no other writer changed any source before mutating the set.
assert all(p.read_bytes()==data for p,data in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Updated 8 navigation documents; preserved original encoding and line endings.')
