"""Prepare719 evidence and next joint-spin/propagation audit."""
import hashlib
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


ledger=(HERE/'unified_physics_condition_ledger_718.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：719局部关系记录、内部参考与原质量反馈\n'+rest
ledger=ledger.replace(
    '接[717全账](unified_physics_condition_ledger_717.md)，回填[718报告](research_note_718.md)。[结果](joint_record_history_transport_results.json)、[核验](research_round_718_checks.json)。',
    '接[718全账](unified_physics_condition_ledger_718.md)，回填[719报告](research_note_719.md)。[结果](joint_local_relational_record_results.json)、[核验](research_round_719_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**719当前增量：** 原两种sterile自旋可分别承担系统和共享参考；局部偶关系读口重建原模式均值但不实现原Luders后态。读取删去参考的指定相干，原Dirac／跳跃给实际反馈；Majorana即时保持却在等待中破坏独立参考。节点操作输入和准备仍明示，连续因果实现及共同尺度未完成。\n\n'+marker)
updates={
    'C03 事件记录':'719原节点三值关系仪器保真实结果；统计重建与原分布模式Luders任务严格区分',
    'C15 物种手征':'719无需新物种的sterile双自旋参考；本仪器保Majorana但改变Dirac／边能量',
    'C19 参考态':'719共享参考的指定相干被读取消去；原Majorana等待产生系统参考关联，不能免费重置',
    'C22 来源反作用':'719实际局部投影的能源／几何差量由原质量和跳跃固定；准备及重建权仍同步'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+'|'):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 719原物质承担局部关系记录和参考

- 同一原sterile双自旋给节点局部偶三值投影，可在指定共享准备下重建原分布模式均值。
- 九个联合结果和真实后态保留；节点非信号不代表连续装置完成，数据须在共同未来汇集。
- 指定参考的r-dagger-s相干在每个结果后消失；关系记录仍保信息，不宣称全部资源消失或给最小擦除功。
- 完整原H的即时反馈仅含相关Dirac和sterile边项，Majorana配对保持；原束缚势控制固定图资源。
- Majorana本身在原真空核态的一阶演化生成系统参考关联。更新两个边缘不自动恢复独立准备或原标定。
- 三组校准不替代完整玻色／Gauss Gibbs或连续物理过程。709参考框架、706宇称及633限制直接复用。

## 本轮合并与下一项

C03／C15／C19／C22共用原两自旋、真实参考、质量及等待。消除了另添参考物种的需要，但不消除准备、轴、仪器、通信和完整历史要求。

接[720](round720_drafts/STATUS.md)：核原全自旋联合结构是否支持可保留参考，同时检验与空间传播及手征字典的一致性。工程设计后置，不把复量子的SU2代数当作任意Hamiltonian的空间旋转对称。
'''
write('unified_physics_condition_ledger_719.md',ledger)
write('round719_drafts/research_note_719_draft.md',(HERE/'research_note_719.md').read_text('utf8'))
write('round720_drafts/STATUS.md','''# 第720轮入口：联合自旋参考与实际空间传播

接[719](../research_note_719.md)、[全账](../unified_physics_condition_ledger_719.md)。原局部关系读数可恢复统计量，但参考不是不变背景；Majorana保局部自旋标量却把两自旋子系统关联起来。

1. 先回查604、614、633、666—667、706及所有相关自旋／方向结果，检查整体自旋共同变换的已有证明，禁止把旧SU2结论重算新轮次。
2. 若原自旋标量跳跃和质量留下联合参考结构，必须同时核该结构怎样接入633的实际Weyl传播和既有空间方向仪器；固定图内部对称不自动等于物理旋转。
3. 将参考、同一等待、可测关系和几何方向作为共同接口。不要只为保参考而关闭原Dirac、Majorana或跳跃，不另添免费外部轴或真空。
4. 先整合条件，装置设计后置。保留382—386、425、522—524、633／634及699准确范围；不重做参考精度／局部相位／一般群表示扫描。
''')
write('round719_drafts/literature_scope_audit.json',json.dumps(dict(
    sources=[
        dict(url='https://arxiv.org/html/2006.03087',authors='Szalay et al.',
             read='Sections7.1-7.4, explicitly local even maps and graded embeddings.',
             scope='Fermionic parity and locality framework, not a continuum detector theorem.'),
        dict(url='https://arxiv.org/html/quant-ph/0610030',authors='Bartlett, Rudolph, Spekkens',
             read='SectionsIII.3 andVI.4; shared reference as consumable relational resource.',
             scope='General reference framework; no imported particle-number SSR for the Majorana model.')],
    inherited='524 scalar probe; 598 coefficients; 633 causal witness; 706 parity; 709 twirling; 718 history.',
    own_mapping='Original sterile spins, actual three-valued instruments, Majorana invariance, Dirac/hop response and reference correlations.',
    no_image_checks=True),ensure_ascii=False,indent=2))
write('round719_drafts/scope_and_dedup_review.md','''# 719范围和推导审查

复用709的参考／群平均及706偶局部映射，未把普通LOCC或超选择写成新普遍发现。数据需共同未来汇集；给定节点Kraus不自动提供连续局域装置。

系统与参考为原sterile的两种自旋，不新增物种。CAR重排使XA XB的相干项为负；实际Kraus含零、正、负三结果。每个结果子空间有确定节点总宇称，故参考r-dagger-s被删，但记录仍可保关系信息。

统计重建的权重不是合法二值后处理概率，且不保原Luders后态。既有719入口的误差界继续限制原任务；本轮没有用另一任务冒充完成。

原完整H中局部旋转保Majorana singlet，独立节点夹断删相关Dirac和sterile跨边。原高势界控制资源，来源包括准备响应。读取能量可为负，未称全部为消耗。

等待中原Majorana一阶生成ar关联，真空切向符号用原256矩阵另核，误差0。一般联合态的四点关联不由两个边缘确定；数值条件二次Gaussian例不推广至全玻色态。

六组64模式系数、四CAR的16维Fock代数及原八CAR的256维演化均保真实矩阵。没有完整Gauss Gibbs数值模拟。旧空间和699结论范围不变，无独立代理复核。
''')
names=('research_note_719.md','joint_local_relational_record.py',
       'joint_local_relational_record_results.json','unified_physics_condition_ledger_719.md')
write('round719_drafts/final_review.txt','Main-agent review only. CAR signs, local projectors, '
      'reference loss, original Majorana/Dirac/hopping, vacuum correlation tangent and source scope checked.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names)+'\n')
print('Prepared719 evidence and720 interface.')
