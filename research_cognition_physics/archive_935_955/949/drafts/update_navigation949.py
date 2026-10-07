"""949 navigation update; preserve all scientific history and entry files."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 949：保持质量与协议的共同弱耦合族'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    oldheading='## 948：中性记录与规范、费米物质的共同响应'
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[949报告]({pre}research_note_949.md)构造联动物理参数族：规范／Yukawa／记录交互及引力减弱，物质质量、946两标量谱和947实际h保持；跨部门领先响应同为ε²。'
        +f'[结果]({pre}949/joint_weak_domain_results.json) · [核验]({pre}949/research_round_949_checks.json)。正式949／累计3734。'
        +'已有条件性共同窗口充分条件，尚无父模型实际联合预算或现实参数验收。'
        +f'停止参数扫描，接[950]({pre}950/drafts/STATUS.md)核一份物理菜单的共同系数及来源。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        s,n=re.subn(r'^\[截至947的取舍审计\][^\n]*$',lambda m:
            '[取舍准则](archive_764_/948/drafts/priority_scope_reaudit.md)保持。948完成树级物质接合；949证明可在减弱交互时保持质量与原内部操作。正式949／累计3734。下一步接[950](archive_764_/950/drafts/STATUS.md)，审同一物理菜单的共同系数／来源与实际预算；不继续门户、参数或仪器细化。',
            s,count=1,flags=re.M);assert n==1
        assert s.count('001—948轮共948份')==1
        s=s.replace('001—948轮共948份','001—949轮共949份')
        replacement='949的联动参数族保持质量谱和内部h，同时保留非零物质／弱引力领先响应；没有认证完整父理论的共同误差窗口。停止参数扫描，接[950](archive_764_/950/drafts/STATUS.md)审同一菜单的系数、来源与预算。完整阶段尚未完成。'
        s,n=re.subn(r'948将同一实际中性事件源[^。\n]*.*?完整阶段尚未完成。',lambda m:replacement,s,count=1);assert n==1
    if p in paths[1:3]:
        current='[949](archive_764_/research_note_949.md)已给保持物质质量、两标量谱和947内部h的联动参数族；948物质响应及领先Newton相位同为ε²。它支持继续寻找共同有限域，尚未签收父模型实际预算、完整物理菜单或现实参数。沿[取舍审计](archive_764_/948/drafts/priority_scope_reaudit.md)，停止参数扫描和局部细化，接[950](archive_764_/950/drafts/STATUS.md)核共同系数与来源。正式949／累计3734；A1、应用目标和定时设置保持。'
        s,n=re.subn(r'^\[947\]\(archive_764_/research_note_947\.md\)已将929[^\n]*$',lambda m:current,s,count=1,flags=re.M);assert n==1
        s=re.sub(r'^\[948共同对象表\][^\n]*\n\n','',s,count=1,flags=re.M)
    if p in paths[3:7]:
        # Preserve the earlier decision as a dated record; show current status.
        s=s.replace('正式948／累计3733，948已结项','正式949／累计3734，949已结项')
        if p==paths[3]:
            s=s.replace('下一项先比较两条分支的对象、来源和观测字典','948树级接合与949参数族已核，下一项审共同菜单与来源')
        if p==paths[4]:
            s=s.replace('## 当前决策记录（不新增轮次）','## 948工作期决策记录（决策本身不计轮次）')
        if p==paths[5]:
            assert s.count('# 231—948轮阶段成果总览')==1 and s.count('231—948的718份')==1
            s=s.replace('# 231—948轮阶段成果总览','# 231—949轮阶段成果总览').replace('231—948的718份','231—949的719份')
        if p==paths[6]:
            s=s.replace('948已完成接合检验，下一项为949','948接合与949参数族已核，下一项为950')
    if p==paths[-1]:
        s,n=re.subn(r'^按\[截至947的取舍审计\][^\n]*$',lambda m:
            '按[取舍审计](../../948/drafts/priority_scope_reaudit.md)，949给保持质量及内部h的共同弱耦合参数族；正式949／累计3734。窗口条件仍依赖同一父对象的实际菜单、系数和来源预算，完整条件未关闭。接[950](../../950/drafts/STATUS.md)核这项共同连接；不以调小参数代替预算，不要求先完成全部UV。A1未修订。',
            s,count=1,flags=re.M);assert n==1
        s=s.replace('六条共同协议：全局缺口对应与检验优先级（截至948）','六条共同协议：全局缺口对应与检验优先级（截至949）',1)
        rows={
          'C04':'|C04 同一内部动力学与控制|1、6、H2；929、935—949|949的联动族保947实际h及材料质量，不把原内部协议一并缩为恒等|仅自由内部块及既有合同；全交互操作运输、终端权限与工程化h仍须明示|',
          'C17':'|C17 质量、Higgs、Yukawa等|1、3、6；931、941—949|949完整Higgs势、规范质量矩阵和Yukawa矩阵保持自由谱，交互可共同变弱|维数、群、质量机制及参数族为物理输入；未认证现实参数或全物质后态|',
          'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；854、898、944—949|949使原物质响应及Newton相位同为ε²，并给保留信号的条件性共同窗口|须取得同一物理菜单的实际系数、归一化次数损失与预算；无完整父理论窗口|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—949|949材料m+h保持，G按ε²缩放；来源仍对固定参数的同一作用变分|沿参数族求导不是应力；全源／接触与有限反作用误差没有由缩放自动认证|'}
        for key,row in rows.items():
            s,n=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert n==1
        s,n=re.subn(r'^4\. 948确认[^\n]*$',lambda m:
            '4. 949排除了“共同减弱交互必然消灭质量或内部协议”的疑虑；未取得父模型的实际共同预算。停止参数扫描，接[950](../../950/drafts/STATUS.md)核一份完整有限菜单的系数与来源，复用854／898而不搬用其梯子数字。不把候选全部自洽性设为纲领前提。',
            s,count=1,flags=re.M);assert n==1
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,raw in updates.items():p.write_bytes(raw)
print('Updated eight living navigation documents for round 949.')
