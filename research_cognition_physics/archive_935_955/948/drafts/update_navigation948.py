"""Update living navigation only, preserving all frozen history."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 948：中性记录与规范、费米物质的共同响应'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    oldheading='## 947：完整有限协议与同一场—引力过程的共同运输'
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[948报告]({pre}research_note_948.md)把947实际中性事件源经946同一两模核接到标准W／费米质量顶点；两种树级响应及反向来源共用一个消元作用，无需给全部记录味添加规范荷。'
        +f'[结果]({pre}948/neutral_matter_bridge_results.json) · [核验]({pre}948/research_round_948_checks.json)。正式948／累计3733。'
        +'质量机制与场表仍为物理输入；新增物质探测器的有限时间联合概率及完整SM／Einstein匹配未签收。'
        +f'停止门户细化，接[949]({pre}949/drafts/STATUS.md)核共同有效域与整体存在性证据。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    # Earlier pre-948 decision text stays in its saved audit; current headings
    # must not continue reporting that round 948 is unfinished.
    prefix_end=s.index(heading)
    prefix=s[:prefix_end]
    prefix=prefix.replace('正式947／累计3732保持，948未结项','正式948／累计3733，948已结项')
    prefix=prefix.replace('正式947／累计3732保持，948尚未结项','正式948／累计3733，948已结项')
    prefix=prefix.replace('正式947／累计3732保持，948未结项','正式948／累计3733，948已结项')
    prefix=prefix.replace('本次没有新增科学成果、轮次或条件关闭，948未结项','取舍本身不计轮次，随后948已完成树级物质接合检验')
    prefix=prefix.replace('本次不新增轮次；948先做方向接合判断','取舍本身不计轮次；948已完成接合检验，下一项为949')
    prefix=prefix.replace('当前执行顺序以该审计为准；正式948／累计3733，948已结项，本次没有新科学结论或条件关闭。',
        '该审计决定的最小检验现已完成；正式948／累计3733，树级物质接口已核，完整父理论有限联合匹配仍开放。下一步以949入口为准。')
    if p in paths[1:3]:
        prefix+=f'[948共同对象表]({pre}948/drafts/common_model_bridge_v0_3.md)现已区分E_rec、M947与P948：保留中性门户树级连接，停止逐个增加响应顶点；接[949]({pre}949/drafts/STATUS.md)审同一有限任务域的共同证据，不把全阶修补设为前提。\n\n'
    s=prefix+s[prefix_end:]
    if p==paths[0]:
        assert s.count('001—947轮共947份')==1
        s=s.replace('001—947轮共947份','001—948轮共948份')
        new='948将同一实际中性事件源接到W及费米质量顶点，正反响应由共同两模核约束；物理输入保留，树级系数不是新增探测概率。[共同对象表](archive_764_/948/drafts/common_model_bridge_v0_3.md)区分各分支，停止门户细化，接[949](archive_764_/949/drafts/STATUS.md)审共同有效域。完整阶段尚未完成。'
        s,n=re.subn(r'947将929[^。\n]*.*?完整阶段尚未完成。',lambda m:new,s,count=1);assert n==1
    if p==paths[5]:
        assert s.count('# 231—947轮阶段成果总览')==1 and s.count('231—947的717份')==1
        s=s.replace('# 231—947轮阶段成果总览','# 231—948轮阶段成果总览').replace('231—947的717份','231—948的718份')
    if p==paths[-1]:
        s=s.replace('六条共同协议：全局缺口对应与检验优先级（截至947）','六条共同协议：全局缺口对应与检验优先级（截至948）',1)
        rows={
          'C14':'|C14 内部规范群、全局形式、连接|1、4、5；375、930—934、946—948|948同一中性事件源经规范不变门户接到W质量顶点，反向来源由同一作用给出|规范群及质量机制仍输入；树级静态系数不等于完整规范量子动力学的有限匹配|',
          'C17':'|C17 质量、Higgs、Yukawa等|1、3、6；925—927、931、941—948|948同一两模核约束W及费米质量响应的比例，未新增记录规范荷|门户、真空、质量机制与参数为物理输入；完整物质后态与现实参数匹配未签收|',
          'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；854、898、944—948|947有限仪器误差已核；948给父作用树级源字典及单模有限动量匹配差|948未把947概率界延伸到新增探测器；同一父理论到有限过程的联合匹配仍开放|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—948|948三源静态Hessian固定正反物质响应；947内部h仍承担同一质量|全应力含物质及门户，不能仅据静态互易或领先Newton项签收完整Einstein反馈|'}
        for key,row in rows.items():
            s,n=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert n==1
        new='4. 948确认中性事件源可经同一门户共同作用于规范及费米质量顶点，正反响应不能分别选定；保留此方向，停止门户细化。接[949](../../949/drafts/STATUS.md)按[共同对象表](../../948/drafts/common_model_bridge_v0_3.md)审有限共同域；不能把树级响应冒充完整父理论联合概率，也不恢复候选全部技术修复为门槛。'
        s,n=re.subn(r'^4\. 947已[^\n]*$',lambda m:new,s,count=1,flags=re.M);assert n==1
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,raw in updates.items():p.write_bytes(raw)
print('Updated eight living navigation documents for round 948.')
