"""Generate723 verification and conflict-safe publication."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m.group()],text)


def write(name,text):
    with (HERE/name).open('x',encoding='utf8',newline='\n') as f:f.write(text)


v=replace((HERE/'verify_round722.py').read_text('utf8'),{
    '722':'723','721':'722','723':'724','3102':'3116','3116':'3130',
    '3371':'3374','1478':'1481',
    'joint_squared_reference_readout':'joint_reference_process_transport',
    'joint_velocity_record_energy':'joint_squared_reference_readout',
    'squared_reference_entry':'joint_reference_entry',
    "text['display_formulas']==22":"text['display_formulas']==18"})
write('verify_round723.py',v)
p=replace((HERE/'publish_round722.py').read_text('utf8'),{
    '722':'723','721':'722','723':'724','3371':'3374','3368':'3371',
    '1478':'1481','3116':'3130','368.':'369.','273.':'274.',
    'joint_squared_reference_readout':'joint_reference_process_transport'})
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**第723轮完成：** [完整物质参考记录与同一有限量子过程]({p}research_note_723.md)原Gauss态可保四参考全部边缘与平均H而改变实际顺序记录；新参考菜单的完整有限过程保记录后态、能源和末端来源，立即联合几何一阶响应也共同收敛。三组、十八式通过，最新723／3374，1481份编号科学文件、3130份保护证据。[核验]({p}research_round_723_checks.json)、[全条件账]({p}unified_physics_condition_ledger_723.md)。变化等待高阶、关系区域和空间极限仍分范围。'")
old=next(x for x in p.splitlines() if x.startswith('order='))
p=p.replace(old,"order='**当前执行顺序（723后，优先于下方历史安排）：** 接[724实际参考与关系区域]({p}round724_drafts/STATUS.md)，回查647—652原选面、因果类型和完整边界来源，使用同一参考过程核实际对象映射。停止单读口、尾部与低秩优化；旧空间及699范围保持，统一目标开放。'")
p=p.replace('原平方参考、读取与能源域','完整物质参考、实际历史与来源')
p=p.replace('旧空间合同复用，单参考读取不自动给完整坐标','旧空间合同复用，完整内部过程仍须接实际区域')
p=p.replace('[原平方参考的能源障碍与可容许读取](research_note_723.md)',
            '[完整物质参考记录与同一有限量子过程](research_note_723.md)')
write('publish_round723.py',p)
p=replace((HERE/'postcheck_round722.py').read_text('utf8'),{
    '722':'723','723':'724','3116':'3130'})
write('postcheck_round723.py',p)
print('Generated723 verification/publication scripts.')
