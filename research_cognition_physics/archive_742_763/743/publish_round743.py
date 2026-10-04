"""Publish743 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round743_navigation_checks.json'
checked=core.read(HERE/'research_round_743_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round743_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第743轮完成：** [原量子标量与严格Gauss关联]({p}research_note_743.md)原完整H的正常Gauss准备产生非高斯四点缺陷；同一标量方差与Ys固定内部关联及纠缠，717读口／来源桥直接复用。无需为解除742限制增加物种，但理想仪器、末读自治及连续参考未签收。两组、十八式通过，最新743／3424，1541份编号科学文件、3401份保护证据。[核验]({p}research_round_743_checks.json)、[条件账]({p}unified_physics_condition_ledger_743.md)。'
order='**当前执行顺序（743后，优先于下方历史安排）：** 接[744原作用诱导的实际仪器]({p}round744_drafts/STATUS.md)，核未知偶编码输入、原sin s末读、完整后态与共同来源；不只比较均值。旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第743轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 389.' not in text
        text+='\n\n## 389. 原量子标量与严格Gauss记录关联\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 294.' not in text
        text+='\n\n## 294. 旧空间合同保持，原相互作用与实际记录共同核验\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—742轮。' in text and '最新科学轮次与检查数为742／3422' in text
        text=text.replace('完成231—742轮。','完成231—743轮。').replace('最新科学轮次与检查数为742／3422','最新科学轮次与检查数为743／3424')
    if p==HERE/'README.md':
        text+='\n\n## 第743轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|743|[原量子标量与严格Gauss记录关联](research_note_743.md)|[代码](joint_native_quantum_record.py)、[结果](joint_native_quantum_record_results.json)、[核验](research_round_743_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round743.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=743,cumulative_tests=3424,numbered_scientific_files=1541,
    unique_protected_evidence_files=3401,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=744,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
