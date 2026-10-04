"""Generate720 verification and conflict-safe publication."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m.group()],text)


def write(name,text):
    with (HERE/name).open('x',encoding='utf8',newline='\n') as f:f.write(text)


v=replace((HERE/'verify_round719.py').read_text('utf8'),{
    '719':'720','718':'719','720':'721','3060':'3074','3074':'3088',
    '3362':'3365','1469':'1472',
    'joint_local_relational_record':'joint_momentum_spin_reference',
    'joint_record_history_transport':'joint_local_relational_record',
    'causal_accuracy_entry':'spin_propagation_entry'})
v=v.replace("text['display_formulas']==20","text['display_formulas']==16")
write('verify_round720.py',v)
p=replace((HERE/'publish_round719.py').read_text('utf8'),{
    '719':'720','718':'719','720':'721','3362':'3365','3359':'3362',
    '1469':'1472','3074':'3088','365.':'366.','270.':'271.',
    'joint_local_relational_record':'joint_momentum_spin_reference'})
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**第720轮完成：** [传播中的内部参考、原质量与实际记录]({p}research_note_720.md)原两质量允许自由均匀分支的FW守恒自旋；局部字典、实际记录与几何来源须同步改变。全部原径向质量背景下，共同平移不变内部参考不能承载SU2。三组、十六式通过，最新720／3365，1472份编号科学文件、3088份保护证据。[核验]({p}research_round_720_checks.json)、[全条件账]({p}unified_physics_condition_ledger_720.md)。结论限于明确参考类别，统一目标开放。'")
old=next(x for x in p.splitlines() if x.startswith('order='))
p=p.replace(old,"order='**当前执行顺序（720后，优先于下方历史安排）：** 接[721原物质坐标与方向记录]({p}round721_drafts/STATUS.md)，复用548—554、573—574、647—651及523／524，核真正缺失的同一量子参考与instrument映射。停止保护spin分类，工程设计后置；旧空间及699范围保持。'")
p=p.replace('局部关系记录、原物质参考及质量反馈','传播参考、原质量与实际记录')
p=p.replace('旧空间合同复用，关系统计不替代原仪器','旧空间合同复用，FW参考不自动是原局部轴')
p=p.replace('[局部关系记录、内部参考与原质量反馈](research_note_720.md)',
            '[传播中的内部参考、原质量与实际记录](research_note_720.md)')
write('publish_round720.py',p)
p=replace((HERE/'postcheck_round719.py').read_text('utf8'),{
    '719':'720','720':'721','3074':'3088'})
write('postcheck_round720.py',p)
print('Generated720 verification/publication scripts.')

