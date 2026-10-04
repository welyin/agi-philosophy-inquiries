"""Generate721 verification and conflict-safe publication."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m.group()],text)


def write(name,text):
    with (HERE/name).open('x',encoding='utf8',newline='\n') as f:f.write(text)


v=replace((HERE/'verify_round720.py').read_text('utf8'),{
    '720':'721','719':'720','721':'722','3074':'3088','3088':'3102',
    '3365':'3368','1472':'1475',
    'joint_momentum_spin_reference':'joint_velocity_record_energy',
    'joint_local_relational_record':'joint_momentum_spin_reference',
    'spin_propagation_entry':'material_velocity_entry'})
write('verify_round721.py',v)
p=replace((HERE/'publish_round720.py').read_text('utf8'),{
    '720':'721','719':'720','721':'722','3365':'3368','3362':'3365',
    '1472':'1475','3088':'3102','366.':'367.','271.':'272.',
    'joint_momentum_spin_reference':'joint_velocity_record_energy'})
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**第721轮完成：** [原物质速度读口的完整能源与互补质量反馈]({p}research_note_721.md)原参考流保完整固定图能源形式域，实际两结果与原等待形成有限能源历史；s／T速度对Dirac／Majorana有互补反馈。固定instrument来源保留，变化读口的总导数另明域。三组、十六式通过，最新721／3368，1475份编号科学文件、3102份保护证据。[核验]({p}research_round_721_checks.json)、[全条件账]({p}unified_physics_condition_ledger_721.md)。自主装置、空间细化及统一目标开放。'")
old=next(x for x in p.splitlines() if x.startswith('order='))
p=p.replace(old,"order='**当前执行顺序（721后，优先于下方历史安排）：** 接[722导数平方参考与关系任务]({p}round722_drafts/STATUS.md)，复用554、652及704—708，核完整参考、实际记录与方向任务的同一性。停止速度函数及常数优化，工程设计后置；旧空间及699范围保持。'")
p=p.replace('传播参考、原质量与实际记录','原速度记录、完整能源与互补质量')
p=p.replace('旧空间合同复用，FW参考不自动是原局部轴','旧空间合同复用，有限能源不自动给联合坐标')
p=p.replace('[传播中的内部参考、原质量与实际记录](research_note_721.md)',
            '[原物质速度读口的完整能源与互补质量反馈](research_note_721.md)')
write('publish_round721.py',p)
p=replace((HERE/'postcheck_round720.py').read_text('utf8'),{
    '720':'721','721':'722','3088':'3102'})
write('postcheck_round721.py',p)
print('Generated721 verification/publication scripts.')

