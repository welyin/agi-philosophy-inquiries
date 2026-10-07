"""One-time update of living navigation; preserve frozen scientific files."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 940：同一物质的弱场交换与正量子联合接口'
oldheading='## 939：共同恢复表与接收记录的有限内部核验'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[940报告]({pre}research_note_940.md)用标准弱场物理传播模式及约束势共同恢复931实际Dirac来源的电磁／引力交换；同一有限正H给概率、能源和来源。'
        +f'[结果]({pre}940/native_exchange_bridge_results.json) · [核验]({pre}940/research_round_940_checks.json)。正式940／累计3725。'
        +'这是输入物理后的树级连接；有限占据边界项显著，未认证该数值过程的连续物理误差，未迁入原曲背景或完整记录。停止本示例精度和高阶修补，'
        +f'接[941]({pre}941/drafts/STATUS.md)先判断整体共同对象的覆盖；完整原生H不是纲领的普遍门槛。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        assert s.count('001—939轮共939份')==1
        s=s.replace('001—939轮共939份','001—940轮共940份')
        old='939完成同一E分支的整体恢复表及接收记录的有限关系核验；末读准备和作用仍需明示。尘埃保持备选，停止参考器件细化，接[940](archive_764_/940/drafts/STATUS.md)回到共同有效过程与来源。完整阶段尚未完成。'
        new='939保留整体恢复表和关系读口；940把实际Dirac来源接入弱场传播与约束的共同正量子接口，不能直接迁入原E_rec或全标准模型。停止该示例精度／高阶细化，接[941](archive_764_/941/drafts/STATUS.md)先比较整体共同对象的覆盖与净输入。完整阶段尚未完成。'
        assert old in s;s=s.replace(old,new,1)
    if p in paths[1:3]:
        old='按939的[计算取舍](archive_764_/939/drafts/common_model_identity_working.md)已完成[整体恢复表](archive_764_/939/drafts/common_model_recovery_map.md)，并检验原853写入与实际末读的接口。有限关系参考保非零信号，但准备和读出操作仍为输入，原作用的自主实现未签收。停止参考器件细化，接[940](archive_764_/940/drafts/STATUS.md)共同有效过程与来源；不要求全场Hamiltonian逐点相等，也不重开仅复现已有二点响应的支线。'
        new='939的[整体恢复表](archive_764_/939/drafts/common_model_recovery_map.md)保持原分支身份；[940](archive_764_/research_note_940.md)新增实际Dirac来源的弱场树级共同量子接口。有限数值的占据边界未忽略，但没有认证其有限耦合物理外推；这只限制该示例，不升级为整个纲领失败。停止精度／高阶修补，接[941](archive_764_/941/drafts/STATUS.md)先比较整体共同对象和阶段覆盖。完整原生H并非普遍门槛；成熟有效物理可复用，实际采信的共同预测及来源仍须同域验收。'
        assert old in s;s=s.replace(old,new,1)
    if p==paths[5]:
        assert s.count('# 231—939轮阶段成果总览')==1 and s.count('231—939的709份')==1
        s=s.replace('# 231—939轮阶段成果总览','# 231—940轮阶段成果总览').replace('231—939的709份','231—940的710份')
    if p==paths[-1]:
        assert s.count('六条共同协议：全局缺口对应与检验优先级（截至939）')==1
        s=s.replace('六条共同协议：全局缺口对应与检验优先级（截至939）','六条共同协议：全局缺口对应与检验优先级（截至940）')
        rows={
          'C04':'|C04 同一内部动力学与控制|1、6、H2；935—940|940用实际物理来源与自由模式构造正量子联合H；仅为明确弱场分支|比较整体共同覆盖；原生H是工具，不是其它有效表达必须跨过的门槛|',
          'C10':'|C10 共同传播几何与普适耦合|1、4、5；334—340、924、927、930、933、940|940同一Dirac来源的电磁及引力协变交换由物理传播模式和约束共同恢复；共锥仍为输入|保持真实物质响应；不把只保TT或共同锥当完整应力与反作用|',
          'C11':'|C11 引力自由度及约束闭合|1、4、5；764—803、861、936—940|940明示线性约束势与两个正TT偏振，无需负范数物理模式；未验非线性全约束|按实际声明阶与分支验收；有限模式正性不代替原曲背景约束闭合|',
          'C20':'|C20 粗化、有效极限与重整化|1、3、4、6、H3；854、898、936、939—940|940树级核匹配通过；有限占据力有明确边界项，有限耦合物理误差未认证|该示例数值暂不作物理外推，停止精度支线；有效阶段不要求UV完成|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；873、935—940|940完整Dirac应力和电流共同匹配，实际有限H保自身能源并含约束来源及占据力修正|同一实际物理菜单的来源与反馈仍须同域验收；不搬用到不同背景或完整标准模型|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert count==1
        old='4. 939已完成[共同恢复表](../../939/drafts/common_model_recovery_map.md)并检验853记录的实际核验接口。有限关系参考保信号，准备／操作／来源仍需明示。停止本参考器件细化；接[940](../../940/drafts/STATUS.md)共同有效过程与来源，不重开尘埃及旧反馈修补链。'
        new='4. 939[共同恢复表](../../939/drafts/common_model_recovery_map.md)与940弱场正量子接口分别保原身份；940给实际物理来源而非任意几何惯量。有限占据外推未认证，停止该示例精度／高阶工作。接[941](../../941/drafts/STATUS.md)先比较整体共同对象和认知操作的覆盖，原生H不是普遍前提；不重开尘埃、参考及旧反馈修补链。'
        assert old in s;s=s.replace(old,new,1)
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==data for p,data in original.items())
for p,data in updates.items():p.write_bytes(data)
print('Updated eight living navigation documents; frozen science preserved.')
