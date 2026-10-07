"""Guarded one-time update: retain dust result, change quantum priority."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 938：新增参考的物理模式与量子复用边界'
oldheading='## 937：正能源内部参考与原全场约束共同接入'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[938报告]({pre}research_note_938.md)核新增参考增加四个物理canonical对；原非Abelian背景的低波数正密度初值族，已有非零参考流及生成元差，不能直接继承旧量子对象。'
        +f'[结果]({pre}938/reference_physical_mode_audit_results.json) · [核验]({pre}938/research_round_938_checks.json)。正式938／累计3723。'
        +'保留937经典结论；这不是尘埃量子化或认知纲领的反证。尘埃转为备选，停止继续修补。'
        +f'接[939]({pre}939/drafts/STATUS.md)整合已有同一E分支的共同模型及物理恢复范围，不再反复更换参考。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        assert s.count('001—937轮共937份')==1
        s=s.replace('001—937轮共937份','001—938轮共938份')
        old='936给原生有效量子化接口；937以新增正能源参考接通原全场经典约束、应力和质量来源。两者尚未组成全物种正量子过程。停止钟／图册优化，接[938](archive_764_/938/drafts/STATUS.md)直接审有限量子恢复合同。'
        new='936给原生有效量子化接口；937保留新增参考的经典接入，但938确认新增物理模式不能直接承接原量子对象，故尘埃降为备选。停止其技术修补，接[939](archive_764_/939/drafts/STATUS.md)整合已有同一E分支的整体恢复映射；共同模型尚未完成。'
        assert old in s;s=s.replace(old,new,1)
    if p in paths[1:3]:
        old='935—936已停止编译器和通用量子化校准；937采用新增正能源参考，实际接通原全场经典约束、参考应力与质量源。它尚非全量子共同模型。938直接审有限量子恢复合同，不继续尘埃图册／精确钟细化；不要求所有旧候选成功。应用目标和定时设置保持。'
        new='935—936已停止编译器和通用量子化校准；937经典参考成果保留。938确认新增参考增加物理模式，低波数即影响生成元，不能直接继承原量子对象。尘埃降为备选，停止其技术修补；939整合已有同一E分支的共同模型与恢复范围，不反复更换参考，不恢复全部旧精度支线为前提。应用目标和定时设置保持。'
        assert old in s;s=s.replace(old,new,1)
    if p==paths[5]:
        assert s.count('# 231—937轮阶段成果总览')==1 and s.count('231—937的707份')==1
        s=s.replace('# 231—937轮阶段成果总览','# 231—938轮阶段成果总览').replace('231—937的707份','231—938的708份')
    if p==paths[-1]:
        assert s.count('六条共同协议：全局缺口对应与检验优先级（截至937）')==1
        s=s.replace('六条共同协议：全局缺口对应与检验优先级（截至937）','六条共同协议：全局缺口对应与检验优先级（截至938）')
        rows={
        'C04':'|C04 同一内部动力学与控制|1、6、H2；935—938|936给原生量子化接口；937给经典根；938确认新参考低波数模式影响生成元，旧量子身份不能直接移植|优先整合同一已有E的模型与恢复范围，不再换参考／编译器来绕过共同身份|',
        'C09':'|C09 钟尺、参考态与可访问性|1、3、6；522—530、861、934—938|937经典内部参考成立；938增加四个物理canonical对及有限来源，不能当免费标签删除|尘埃降为备选；保留原材料参考及明确范围，不继续新参考的全部量子技术|',
        'C11':'|C11 引力自由度及约束闭合|1、4、5；764—803、861、936—938|原E量子结果各守其范围；938新总约束核与旧非尘埃物理空间不同|一份候选的约束／态／来源共同映射；不因经典约化简单就搬用另一量子对象|',
        'C18':'|C18 耦合、质量、混合与自由参数|1、4、5、6；930—931、935—938|937质量源属于新参考；938低波数初值族的生成元和变分保同一根|参数、初值族与作用源分列；允许选择参数，不把来源差异藏在坐标重命名中|',
        'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；585—590、873、935—938|938全场低波数参考流由总约束补偿，差异属有限尺度；全量子反作用仍未签收|继续原同一量子分支的有限共同恢复，不无限修补新增候选|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert count==1
        old='4. 936给原生有效量子化接口；937以新增正能源参考接通原全场经典约束及来源。参考物种、作用与维数仍是物理输入；未完成共同量子模型。接[938](../../938/drafts/STATUS.md)直接审有限量子恢复合同，停止尘埃／钟和通用量子化校准细化。'
        new='4. 937经典新增参考保留；938确认新增物理模式不能免费继承旧量子空间，且低波数生成元已有差异。尘埃降为备选，停止技术扩展。接[939](../../939/drafts/STATUS.md)整合已有同一E分支的整体模型与恢复证据，不反复更换参考或重开全部精度修补。'
        assert old in s;s=s.replace(old,new,1)
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==data for p,data in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Updated 8 living navigation documents; historical results retained.')
