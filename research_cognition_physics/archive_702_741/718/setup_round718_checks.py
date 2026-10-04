"""Generate718 verification and conflict-safe navigation publication."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent


def write(name,text):
    with (HERE/name).open('x',encoding='utf8',newline='\n') as f:f.write(text)


def replace(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m.group()],text)


v=replace((HERE/'verify_round717.py').read_text('utf8'),{
    '717':'718','716':'717','718':'719','3031':'3046','3046':'3060',
    '3356':'3359','1463':'1466','joint_record_mass_feedback':'joint_record_history_transport',
    'joint_smooth_mode_contract':'joint_record_mass_feedback',
    "text['display_formulas']==22":"text['display_formulas']==20",
    'postrecord_reference_entry':'reflected_history_entry',
    "==(4,0,0)":"==(3,0,0)","fresh_tests=dict(run=4":"fresh_tests=dict(run=3"})
v=v.replace(",'round718_drafts/force_results_before_bounded_record.json'","")
write('verify_round718.py',v)
p=replace((HERE/'publish_round717.py').read_text('utf8'),{
    '717':'718','716':'717','718':'719','3356':'3359','3352':'3356',
    '1463':'1466','3046':'3060','363.':'364.','268.':'269.',
    'joint_record_mass_feedback':'joint_record_history_transport'})
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**第718轮完成：** [保留真实记录的有限时间与几何输送]({p}research_note_718.md)同一局部资源控制首次／后续记录和末态的有限时间差异；移动仪器、交叉端点及热准备响应须同步保留。三组、二十式通过，最新718／3359，1466份编号科学文件、3060份保护证据。[核验]({p}research_round_718_checks.json)、[全条件账]({p}unified_physics_condition_ledger_718.md)。共同空间极限及因果实现仍开放。'")
old=next(x for x in p.splitlines() if x.startswith('order='))
p=p.replace(old,"order='**当前执行顺序（718后，优先于下方历史安排）：** 接[719真实记录与区域因果]({p}round719_drafts/STATUS.md)，复用524／633／634／667及区域Gauss合同，核同一过程的准确度与因果限制。停止一般距离、反射和高矩扫描；工程设计后置，旧空间及699范围保持。'")
p=p.replace('原记录、质量力及实际有界效果','原真实记录、有限时间及几何来源')
p=p.replace('旧空间合同复用，反馈系数不代替物理连续','旧空间合同复用，有限时间界不代替光锥')
p=p.replace('[原记录的质量反馈与有界读数](research_note_718.md)',
            '[保留真实记录的有限时间与几何输送](research_note_718.md)')
write('publish_round718.py',p)
p=replace((HERE/'postcheck_round717.py').read_text('utf8'),{'717':'718','718':'719','3046':'3060'})
write('postcheck_round718.py',p)
print('Generated718 verification and publication scripts.')
