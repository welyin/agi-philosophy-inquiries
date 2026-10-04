"""Generate717 verification and conflict-safe navigation publication."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent


def write(name,text):
    with (HERE/name).open('x',encoding='utf8',newline='\n') as f:f.write(text)


def replace(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m.group()],text)


v=replace((HERE/'verify_round716.py').read_text('utf8'),{
    '716':'717','715':'716','717':'718','3016':'3031','3031':'3046',
    '3352':'3356','1460':'1463','joint_smooth_mode_contract':'joint_record_mass_feedback',
    'joint_sharp_record_domain':'joint_smooth_mode_contract',
    "text['display_formulas']==20":"text['display_formulas']==22",
    'sterile_mode_entry':'postrecord_reference_entry',
    'first_action_diagnostic.json':'force_results_before_bounded_record.json',
    "==(3,0,0)":"==(4,0,0)","fresh_tests=dict(run=3":"fresh_tests=dict(run=4"})
write('verify_round717.py',v)

p=replace((HERE/'publish_round716.py').read_text('utf8'),{
    '716':'717','715':'716','717':'718','3352':'3356','3349':'3352',
    '1460':'1463','3031':'3046','362.':'363.','267.':'268.',
    'joint_smooth_mode_contract':'joint_record_mass_feedback'})
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**第717轮完成：** [原记录的质量反馈与有界读数]({p}research_note_717.md)原质量力决定读后玻色二阶反馈；有界能源不能控制全部无界势响应，但旧sin s效果有态无关的精确二阶界。四组、二十二式通过，最新717／3356，1463份编号科学文件、3046份保护证据。[核验]({p}research_round_717_checks.json)、[全条件账]({p}unified_physics_condition_ledger_717.md)。有限时间共同历史与物理连续仍开放。'")
old=next(x for x in p.splitlines() if x.startswith('order='))
p=p.replace(old,"order='**当前执行顺序（717后，优先于下方历史安排）：** 接[718原有界记录的有限时间历史]({p}round718_drafts/STATUS.md)，核同一读取的H与RHR比较，复用623—625、634及704。停止高矩反例／波包优化；旧空间、633因果及699范围保持。'")
p=p.replace('原物质模的体积、资源与几何共同条件','原记录、质量力及实际有界效果')
p=p.replace('旧空间合同复用，波包归一不选择物理维数','旧空间合同复用，反馈系数不代替物理连续')
p=p.replace('[光滑物质模的共同资源与来源](research_note_717.md)',
            '[原记录的质量反馈与有界读数](research_note_717.md)')
write('publish_round717.py',p)

p=replace((HERE/'postcheck_round716.py').read_text('utf8'),{'716':'717','717':'718','3031':'3046'})
write('postcheck_round717.py',p)
print('Generated717 verification and publication scripts.')
