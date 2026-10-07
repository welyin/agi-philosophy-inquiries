"""Update living navigation for the finite-field connection; preserve frozen rounds."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 944：有限场通信与同一反作用的受控连接'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    oldheading='## 943：内部记录能否进入共同物理交互'
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[944报告]({pre}research_note_944.md)在明确重粒子—平滑KG场有效模型中，由同一H给记录、噪声及机械交换；有限质量概率差解析下界>0.0321150，定中心相位梯度与冲量一致。'
        +f'[结果]({pre}944/finite_field_communication_results.json) · [核验]({pre}944/research_round_944_checks.json)。正式944／累计3729。'
        +'三维、材料和交互仍为物理输入；原Dirac／几何匹配与无界应力误差未签收。停止本模型优化，'
        +f'接[945]({pre}945/drafts/STATUS.md)按整体覆盖选择共同匹配。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        assert s.count('001—943轮共943份')==1
        s=s.replace('001—943轮共943份','001—944轮共944份')
        old='939整体恢复表及940—942接口保持原范围；943判定仅质量／几何耦合不能自动输出929逻辑内容，明确新交互可解除该访问限制。保留材料工具、停止编码细化，接[944](archive_764_/944/drafts/STATUS.md)选择同一记录—传播—来源的有效交互，不把修好某候选作为整个纲领门槛。完整阶段尚未完成。'
        new='939整体恢复表及940—943接口保持原范围；944在明确有效材料中连接有限记录、场噪声和机械交换，并给有限质量误差界。停止材料优化，接[945](archive_764_/945/drafts/STATUS.md)审计同一物理对象的覆盖与匹配，不能把不同模型各自成功相加签收。完整阶段尚未完成。'
        assert old in s;s=s.replace(old,new,1)
    if p in paths[1:3]:
        new='[944](archive_764_/research_note_944.md)给重粒子—平滑KG场有效分支的有限共同连接：同一H产生记录、噪声、场能与机械力；有限质量全态误差≤0.000451823，实际模型关系读口概率差>0.0321150。该实现为新增物理输入，原Dirac／Einstein匹配及无界应力误差仍未签收。停止质量、形状与仪器细化，接[945](archive_764_/945/drafts/STATUS.md)从整体覆盖选择共同匹配，不把这个候选的全部修复列为纲领前提。'
        s,count=re.subn(r'^\[943\]\(archive_764_/research_note_943\.md\)已判定942材料[^\n]*$',lambda m:new,s,count=1,flags=re.M);assert count==1,p
    if p==paths[5]:
        assert s.count('# 231—943轮阶段成果总览')==1 and s.count('231—943的713份')==1
        s=s.replace('# 231—943轮阶段成果总览','# 231—944轮阶段成果总览').replace('231—943的713份','231—944的714份')
    if p==paths[-1]:
        old='六条共同协议：全局缺口对应与检验优先级（截至943）';assert s.count(old)==1
        s=s.replace(old,'六条共同协议：全局缺口对应与检验优先级（截至944）')
        rows={
          'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；854、898、944|944对明确重粒子有效H给保未知输入及旧参考的有限质量全态误差界|只运输有界概率；原Dirac匹配与无界来源误差不由该界自动取得|',
          'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；853、939、943—944|944同一标量场交互产生有限关系读数，运动模型概率差有正下界|准备与读口沿939为输入；全部六协议及原物种匹配未签收，停止材料优化|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—944|944同一H给场能及机械交换，定中心相位梯度等于冲量；运动模型总动量守恒|机械反作用不等于已完成几何反馈；共同作用匹配与有界任务所需来源预算分开验收|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert count==1
        old='4. 943已完成材料内容访问的决定性检验：942的纯质量／几何接口不能自动输出929逻辑记录，新增味交互可给有限区别。保留存在性工具，停止编码／仪器修补；接[944](../../944/drafts/STATUS.md)选择共同记录—传播—来源的有效交互。939整体恢复表及940—942接口不拼接成已完成模型，也不以修复此候选作为纲领前提。'
        new='4. 944已取得有限场通信与机械反作用的共同连接，有限质量概率差有严格正下界。停止本有效材料优化；接[945](../../945/drafts/STATUS.md)以整体覆盖选择物理匹配。939局部形式分支、942材料工具与944有效H各有范围，未构成已完成的统一模型；任何单一候选的全部修复均非纲领普遍门槛。'
        assert old in s;s=s.replace(old,new,1)
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,raw in updates.items():p.write_bytes(raw)
print('Updated eight living navigation documents; preserved frozen science.')
