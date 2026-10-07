"""Checked navigation update; historical reports and science stay frozen."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 947：完整有限协议与同一场—引力过程的共同运输'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    oldheading='## 946：规范不变门户中的记录、径向物质与弱场引力'
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[947报告]({pre}research_note_947.md)把929两次有限协议接入946同一有效H，实际事件的质量对易耦合保条件量子后态；全未知输入与16故障扇区的联合仪器误差≤0.000633467，三种物理信号仍有正下界。'
        +f'[结果]({pre}947/protocol_field_transport_results.json) · [核验]({pre}947/research_round_947_checks.json)。正式947／累计3732。'
        +'原终端逻辑权限与新增材料仍为输入，完整SM／Einstein匹配、全部仪器未签收。'
        +f'停止协议扩展，接[948]({pre}948/drafts/STATUS.md)收敛共同模型及关键物理恢复。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        assert s.count('001—946轮共946份')==1
        s=s.replace('001—946轮共946份','001—947轮共947份')
        new='947将929原两次有限六协议合同与946场／Newton反作用共同运输；联合仪器误差受控，实际事件有非零物理读数。原终端权限和材料成本仍输入，不能冒充全部仪器。[存在性范围](archive_764_/945/drafts/common_model_scope_v0_2.md)保持，停止协议细化，接[948](archive_764_/948/drafts/STATUS.md)收敛同一模型并选择关键SM／Einstein恢复证据。完整阶段尚未完成。'
        s,count=re.subn(r'946以明确规范不变门户，[^\n]*?完整阶段尚未完成。',lambda m:new,s,count=1);assert count==1,p
    if p in paths[1:3]:
        new='[947](archive_764_/research_note_947.md)已将929原有限六协议合同接到同一场—Newton有效H：实际事件、条件Q—参考、两次运行和16故障扇区共同运输，联合仪器误差≤0.000633467。原终端权限、工程化材料及物理输入保留；全部仪器和完整SM／Einstein匹配仍未签收。沿[存在性范围](archive_764_/945/drafts/common_model_scope_v0_2.md)，停止协议扩展，接[948](archive_764_/948/drafts/STATUS.md)收敛共同模型并选择关键物理恢复，不先修完当前候选。'
        s,count=re.subn(r'^\[946\]\(archive_764_/research_note_946\.md\)已在[^\n]*$',lambda m:new,s,count=1,flags=re.M);assert count==1,p
    if p==paths[5]:
        assert s.count('# 231—946轮阶段成果总览')==1 and s.count('231—946的716份')==1
        s=s.replace('# 231—946轮阶段成果总览','# 231—947轮阶段成果总览').replace('231—946的716份','231—947的717份')
    if p==paths[-1]:
        old='六条共同协议：全局缺口对应与检验优先级（截至946）';assert s.count(old)==1
        s=s.replace(old,'六条共同协议：全局缺口对应与检验优先级（截至947）')
        rows={
          'C03':'|C03 事件身份、关系记录、访问|1、2、3、4；391、923、929、934、946—947|947实际事件进入共同物理场，条件Q—旧参考保留；三物理效果对易|原终端逻辑权限仍输入，不能从有效PVM自动领取任意空间读取装置|',
          'C04':'|C04 同一内部动力学与控制|1、6、H2；929、935—943、947|947原h执行两次协议并进入实际质量；B仅修饰头两个钟位即可与h对易|工程化h、时标和事件敏感耦合仍输入；不重新制造通用处理器或优化编码|',
          'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；854、898、944—947|947联合仪器含未知输入、旧参考、16故障和有限时间窗，误差≤0.000633467|有界操作运输已核；完整父理论来源及物理恢复的误差仍需同一字典|',
          'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；853、929、939、942—947|929原有限六协议合同已在场／Newton有效过程运输，实际事件有外向信号|原终端权限及百万维中性材料保留；全部实际仪器、空间主体与完整SM恢复未签收|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—947|947同一h承担真实内部更新和质量，内部引力钟速未删；场及能源账沿946|完整场应力和非线性父理论匹配仍开放；概率界不自动成为无界来源界|',
          'C23':'|C23 统计、退相干与时间箭头|1、4、6；929、944—947|947同一过程两次运行、内部空白和故障资料保持；环境可记录已知事件|仅原有限故障／重复合同；不把条件后态保持说成测前未知态完全不受扰动|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert count==1
        new='4. 947已把原两次有限六协议合同接入同一场—Newton过程，联合仪器与三项物理信号有共同误差界；原访问和材料成本保留。停止协议／门户／仪器优化，接[948](../../948/drafts/STATUS.md)收敛共同模型及关键物理恢复。[存在性口径](../../945/drafts/common_model_scope_v0_2.md)保持；候选全部技术修复不是纲领前提。'
        s,count=re.subn(r'^4\. 946使记录[^\n]*$',lambda m:new,s,count=1,flags=re.M);assert count==1
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,raw in updates.items():p.write_bytes(raw)
print('Updated eight living navigation documents for round 947.')
