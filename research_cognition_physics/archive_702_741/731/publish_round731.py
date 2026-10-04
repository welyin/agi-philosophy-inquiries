"""Publish731 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round731_navigation_checks.json'
checked=core.read(HERE/'research_round_731_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round731_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第731轮完成：** [原物质上的联合初始约束补偿与量子来源接口]({p}research_note_731.md)原非平坦场允许全电弱Gauss及三方向动量共同补偿，全部代价放回新的正共形解；颜色条件和同态一阶范围明确。三组、二十式通过，最新731／3398，1505份编号科学文件、3242份保护证据。[核验]({p}research_round_731_checks.json)、[全条件账]({p}unified_physics_condition_ledger_731.md)。实际绝对量子源、自洽反馈与全时间发展仍开放。'
order='**当前执行顺序（731后，优先于下方历史安排）：** 接[732实际同态来源与修正背景]({p}round732_drafts/STATUS.md)，复用原初始右逆，核完整源、参考处方、Ward及高阶资料的共同闭合。停止补偿菜单、ε和椭圆精度优化；旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第731轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 377.' not in text
        text+='\n\n## 377. 原物质上的联合初始约束补偿与量子来源接口\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 282.' not in text
        text+='\n\n## 282. 旧空间合同保持，初始约束补偿不替代完整自洽发展\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—730轮。' in text and '最新科学轮次与检查数为730／3395' in text
        text=text.replace('完成231—730轮。','完成231—731轮。').replace('最新科学轮次与检查数为730／3395','最新科学轮次与检查数为731／3398')
    if p==HERE/'README.md':
        text+='\n\n## 第731轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|731|[原物质上的联合初始约束补偿与量子来源接口](research_note_731.md)|[代码](joint_source_constraint_response.py)、[结果](joint_source_constraint_response_results.json)、[核验](research_round_731_checks.json)|\n'
    planned[p]=text.replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir(exist_ok=False);manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round731.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=731,cumulative_tests=3398,numbered_scientific_files=1505,
    unique_protected_evidence_files=3242,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=732,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
