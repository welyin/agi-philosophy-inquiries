"""Publish the unnumbered shared-event review with CAS navigation snapshots."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent
ROOT = RESEARCH.parent
checks = json.loads((HERE/'joint_event_conjecture_checks.json').read_text('utf8'))
assert checks['all_reported_checks_passed']
assert checks['total_protected_including_this_review'] == 1013
paths = [ROOT/'README.md', RESEARCH/'README.md', RESEARCH/'research_direction.md',
         RESEARCH/'RESEARCH_STATE.md', HERE/'README.md',
         HERE/'spatial_premise_closure_audit.md',
         HERE/'three_dimensional_four_conditions_review.md']
baseline = {p:p.read_bytes() for p in paths}
snapshot = HERE/'joint_event_conjecture_navigation_before'
snapshot.mkdir(exist_ok=False)
planned = {}
for p, raw in baseline.items():
    enc = 'utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl = '\r\n' if b'\r\n' in raw else '\n'
    body = raw.decode(enc).replace('\r\n', '\n')
    assert '2026-09-30共同指认猜想审查' not in body
    prefix = ('research_cognition_physics/archive_231_/' if p.parent == ROOT else
              'archive_231_/' if p.parent == RESEARCH else '')
    entry = (
        f'**2026-09-30共同指认猜想审查（不增轮次）：** '
        f'[完整评议]({prefix}joint_event_conjecture_review.md)、'
        f'[诊断结果]({prefix}joint_event_conjecture_probe_results.json)、'
        f'[核验入口]({prefix}verify_joint_event_conjecture_review.py)。'
        '采用多主体共同确认事件作为任务入口，允许记录、传播与几何联合检验；'
        '六条自动物理推论及“共同指认等价FUCP”尚不成立。'
        '纠正502固定采样、515—516同H不同读口、520完整仪器的引用范围；'
        '经典冗余、GHZ验证范围与私有量子信息的3项诊断通过。'
        '不撤回原FUCP，不把经典短协议当强化自主宇宙反例。'
        '编号仍520／2552，871份科学文件；本审查5份证据含两版稿，保护总数1013。')
    if p.name == 'spatial_premise_closure_audit.md':
        body += '\n\n## 159. 共同指认猜想的蕴涵与来源审查（未编号）\n\n'
    elif p.name == 'three_dimensional_four_conditions_review.md':
        body += '\n\n## 64. 共同事件协议仍须给出维数选择条件\n\n'
    else:
        body += '\n\n'
    body += entry + '\n'
    if p in (RESEARCH/'research_direction.md', RESEARCH/'RESEARCH_STATE.md',
             HERE/'spatial_premise_closure_audit.md'):
        body += (
            '\n**当前接续：** 先复用391已完成的最大共同尖锐记录与287的真实地图重叠工具，'
            '检验动态模型中的实际双端接触、记录来源、有限寿命及跨主体事件对应；'
            '只有关闭来源或空间选择缺口才新增编号。仅生成共同记录、重命名协变或一般滤波合并，'
            '不足以签收三维。三维仍为近期主线；允许动力学联合约束，不将标准模型完成设为前置。'
            '520后的局部读取候选尚未冻结为521，不把原本有条件的局部任务改成全态可知公理。'
            '后继核验须承接本次verify_joint_event_conjecture_review的1013份证据。\n')
    if p.name == 'three_dimensional_four_conditions_review.md':
        body += ('\n有限共同事件协议、被动重命名与经典冗余均未提供局部位置群、'
                 '一致半幅、保组合重定向和完整方向合同。不能以洛伦兹群名称、'
                 '共享标签或GHZ资源代替这些空间输入；尚无三维／GR结项。\n')
    planned[p] = body.replace('\n', nl).encode(enc)
    label = ('project_README.md' if p == ROOT/'README.md' else
             'research_'+p.name if p.parent == RESEARCH else 'archive_'+p.name)
    with (snapshot/label).open('xb') as f:
        f.write(raw)
for p, raw in baseline.items():
    assert p.read_bytes() == raw, str(p)
for p, data in planned.items():
    assert p.read_bytes() == baseline[p], str(p)
    p.write_bytes(data)
with (snapshot/'manifest.json').open('x', encoding='utf8') as f:
    json.dump({str(p.relative_to(ROOT)):{
        'before':hashlib.sha256(baseline[p]).hexdigest(),
        'after':hashlib.sha256(data).hexdigest()} for p,data in planned.items()},
        f, ensure_ascii=False, indent=2)
print(json.dumps(dict(navigation_files_updated=len(planned), snapshots_preserved=len(baseline))))
