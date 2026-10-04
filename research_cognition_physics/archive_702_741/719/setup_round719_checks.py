"""Generate719 verification and conflict-safe publication."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m.group()],text)


def write(name,text):
    with (HERE/name).open('x',encoding='utf8',newline='\n') as f:f.write(text)


v=replace((HERE/'verify_round718.py').read_text('utf8'),{
    '718':'719','717':'718','719':'720','3046':'3060','3060':'3074',
    '3359':'3362','1466':'1469','2026-10-03':'2026-10-04',
    'joint_record_history_transport':'joint_local_relational_record',
    'joint_record_mass_feedback':'joint_record_history_transport',
    'reflected_history_entry':'causal_accuracy_entry'})
write('verify_round719.py',v)
p=replace((HERE/'publish_round718.py').read_text('utf8'),{
    '718':'719','717':'718','719':'720','3359':'3362','3356':'3359',
    '1466':'1469','3060':'3074','364.':'365.','269.':'270.',
    '20261003':'20261004','2026-10-03':'2026-10-04',
    'joint_record_history_transport':'joint_local_relational_record'})
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**第719轮完成：** [局部关系记录、内部参考与原质量反馈]({p}research_note_719.md)原sterile双自旋可承担关系读取及共享参考；能重建均值却不保原Luders历史。读取消去指定参考相干，Dirac／跳跃给反馈，Majorana在等待中生成关联。三组、二十式通过，最新719／3362，1469份编号科学文件、3074份保护证据。[核验]({p}research_round_719_checks.json)、[全条件账]({p}unified_physics_condition_ledger_719.md)。连续因果装置及共同尺度仍开放。'")
old=next(x for x in p.splitlines() if x.startswith('order='))
p=p.replace(old,"order='**当前执行顺序（719后，优先于下方历史安排）：** 接[720联合自旋参考与空间传播]({p}round720_drafts/STATUS.md)，先回查604／614／633／666—667等，核同一参考结构与实际Weyl方向传播的兼容性。工程设计后置，停止一般局部相位和参考精度扫描；旧空间及699范围保持。'")
p=p.replace('原真实记录、有限时间及几何来源','局部关系记录、原物质参考及质量反馈')
p=p.replace('旧空间合同复用，有限时间界不代替光锥','旧空间合同复用，关系统计不替代原仪器')
p=p.replace('[保留真实记录的有限时间与几何输送](research_note_719.md)',
            '[局部关系记录、内部参考与原质量反馈](research_note_719.md)')
write('publish_round719.py',p)
p=replace((HERE/'postcheck_round718.py').read_text('utf8'),{
    '718':'719','719':'720','3060':'3074','20261003':'20261004','2026-10-03':'2026-10-04'})
write('postcheck_round719.py',p)
print('Generated719 verification/publication scripts.')
