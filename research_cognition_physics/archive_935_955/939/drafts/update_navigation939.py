"""One-time update of living navigation; preserve historical round files."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 939：共同恢复表与接收记录的有限内部核验'
oldheading='## 938：新增参考的物理模式与量子复用边界'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[939报告]({pre}research_note_939.md)完成[整体恢复表]({pre}939/drafts/common_model_recovery_map.md)，核853末读的实际访问条件：中性无参考读口不能自动读取原Z；有限内部关系参考保原三方信号的1/8。'
        +f'[结果]({pre}939/receiver_relational_access_results.json) · [核验]({pre}939/research_round_939_checks.json)。正式939／累计3724。'
        +'新增准备与保荷读出操作为明确输入，未证明原作用自主产生装置；原853写入结论保持。停止参考器件细化，'
        +f'接[940]({pre}940/drafts/STATUS.md)回到共同有效过程与物理来源。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        assert s.count('001—938轮共938份')==1
        s=s.replace('001—938轮共938份','001—939轮共939份')
        old='936给原生有效量子化接口；937保留新增参考的经典接入，但938确认新增物理模式不能直接承接原量子对象，故尘埃降为备选。停止其技术修补，接[939](archive_764_/939/drafts/STATUS.md)整合已有同一E分支的整体恢复映射；共同模型尚未完成。'
        new='939完成同一E分支的整体恢复表及接收记录的有限关系核验；末读准备和作用仍需明示。尘埃保持备选，停止参考器件细化，接[940](archive_764_/940/drafts/STATUS.md)回到共同有效过程与来源。完整阶段尚未完成。'
        assert old in s;s=s.replace(old,new,1)
    if p in paths[1:3]:
        old='939的[具体取舍与下一次计算标准](archive_764_/939/drafts/common_model_identity_working.md#4-本次提醒后的实际取舍2026-10-07)已落实：E_rec暂作证据比较基线，先完成整体恢复表，再判断其共同连接价值；不先承诺补齐该候选的全场Hamiltonian，也不另开只复现已有二点响应的技术支线。有效描述须保声明任务的共同预测与来源，不要求与原全场生成元逐点相等。本次不增加科学轮次。'
        new='按939的[计算取舍](archive_764_/939/drafts/common_model_identity_working.md)已完成[整体恢复表](archive_764_/939/drafts/common_model_recovery_map.md)，并检验原853写入与实际末读的接口。有限关系参考保非零信号，但准备和读出操作仍为输入，原作用的自主实现未签收。停止参考器件细化，接[940](archive_764_/940/drafts/STATUS.md)共同有效过程与来源；不要求全场Hamiltonian逐点相等，也不重开仅复现已有二点响应的支线。'
        assert old in s;s=s.replace(old,new,1)
    if p==paths[5]:
        assert s.count('# 231—938轮阶段成果总览')==1 and s.count('231—938的708份')==1
        s=s.replace('# 231—938轮阶段成果总览','# 231—939轮阶段成果总览').replace('231—938的708份','231—939的709份')
    if p==paths[-1]:
        assert s.count('六条共同协议：全局缺口对应与检验优先级（截至938）')==1
        s=s.replace('六条共同协议：全局缺口对应与检验优先级（截至938）','六条共同协议：全局缺口对应与检验优先级（截至939）')
        rows={
          'C04':'|C04 同一内部动力学与控制|1、6、H2；935—939|939整体恢复表区分原作用的实际写入与另选末读；936仍为物理符号之后的工具|共同过程和来源优先；不把任一候选完整Hamiltonian修复当普遍前提|',
          'C09':'|C09 钟尺、参考态与可访问性|1、3、6；522—530、861、934—939|939原853味中性末读有限制；同种材料的关系参考保持指定信号|参考准备及读出作用单列输入；无需先做无限精度参考，尘埃仍备选|',
          'C20':'|C20 粗化、有效极限与重整化|1、3、4、6、H3；854、898、936、939|939关系CAR读口可进入同一有限菜单；共同恢复表保持各分支身份|原有限物理系数及来源预算仍需共同映射；不以单个读口签收全部有效物理|',
          'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；853、929、932、939|939同种接收材料的有限关系核验存在，原指定三方信号保留1/8；独立平均丢失关系|新增准备和保荷读口是输入，原场作用生成装置尚未证明；停止器件细化|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；873、935—939|939原味荷保持可核；参考准备改变绝对来源，不能免费移植旧源账|阶段仍需同一有限过程与几何反作用；味守恒并不等于完整能源或Einstein约束验收|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert count==1
        old='4. 937经典新增参考保留；938确认新增物理模式不能免费继承旧量子空间，且低波数生成元已有差异。尘埃降为备选，停止技术扩展。接[939](../../939/drafts/STATUS.md)整合已有同一E分支的整体模型与恢复证据，不反复更换参考或重开全部精度修补。'
        new='4. 939已完成[共同恢复表](../../939/drafts/common_model_recovery_map.md)并检验853记录的实际核验接口。有限关系参考保信号，准备／操作／来源仍需明示。停止本参考器件细化；接[940](../../940/drafts/STATUS.md)共同有效过程与来源，不重开尘埃及旧反馈修补链。'
        assert old in s;s=s.replace(old,new,1)
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==data for p,data in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Updated eight living navigation documents; historical science preserved.')
