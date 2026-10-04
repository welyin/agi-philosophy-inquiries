"""Generate722 verification and conflict-safe publication from721."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m.group()],text)


def write(name,text):
    with (HERE/name).open('x',encoding='utf8',newline='\n') as f:f.write(text)


v=replace((HERE/'verify_round721.py').read_text('utf8'),{
    '721':'722','720':'721','722':'723','3088':'3102','3102':'3116',
    '3368':'3371','1475':'1478',
    'joint_velocity_record_energy':'joint_squared_reference_readout',
    'joint_momentum_spin_reference':'joint_velocity_record_energy',
    'material_velocity_entry':'squared_reference_entry',
    "text['display_formulas']==16":"text['display_formulas']==22"})
write('verify_round722.py',v)
p=replace((HERE/'publish_round721.py').read_text('utf8'),{
    '721':'722','720':'721','722':'723','3368':'3371','3365':'3368',
    '1475':'1478','3102':'3116','367.':'368.','272.':'273.',
    'joint_velocity_record_energy':'joint_squared_reference_readout'})
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**第722轮完成：** [原平方参考的能源障碍与可容许读取]({p}research_note_722.md)指定平方相位读法使正常Gauss有限能源态出原能源域；完整Q含空间差分的有理instrument却保原能源、实际记录及几何一阶总响应。三组、二十二式通过，最新722／3371，1478份编号科学文件、3116份保护证据。[核验]({p}research_round_722_checks.json)、[全条件账]({p}unified_physics_condition_ledger_722.md)。固定图结论，联合坐标、跨尺度及统一目标开放。'")
old=next(x for x in p.splitlines() if x.startswith('order='))
p=p.replace(old,"order='**当前执行顺序（722后，优先于下方历史安排）：** 接[723同一完整参考与关系任务]({p}round723_drafts/STATUS.md)，复用554、647—652、704及707—708，核原四参考、共同态、区域和实际记录映射。停止谱函数、尾部及常数优化，工程设计后置；旧空间及699范围保持。'")
p=p.replace('原速度记录、完整能源与互补质量','原平方参考、读取与能源域')
p=p.replace('旧空间合同复用，有限能源不自动给联合坐标','旧空间合同复用，单参考读取不自动给完整坐标')
p=p.replace('[原物质速度读口的完整能源与互补质量反馈](research_note_722.md)',
            '[原平方参考的能源障碍与可容许读取](research_note_722.md)')
write('publish_round722.py',p)
p=replace((HERE/'postcheck_round721.py').read_text('utf8'),{
    '721':'722','722':'723','3102':'3116'})
write('postcheck_round722.py',p)
print('Generated722 verification/publication scripts.')
