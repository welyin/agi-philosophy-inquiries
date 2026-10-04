"""Publish751 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round751_navigation_checks.json'
checked=core.read(HERE/'research_round_751_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round751_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第751轮完成：** [共同报告与内部来源账]({p}research_note_751.md)原未知输入可有相同报告而不同来源，甚至同一完整cq后态却注能相反；限定精确端点可加实现受原Majorana的守恒障碍。保留有限精度／联合关联方案，不扩大到全GR。两组、十六式通过，最新751／3441，1565份编号科学文件、3521份保护证据。[核验]({p}research_round_751_checks.json)、[条件账]({p}unified_physics_condition_ledger_751.md)。'
order='**当前执行顺序（751后，优先于下方历史安排）：** 接[752联合假说的共同对象]({p}round752_drafts/STATUS.md)，沿障碍分类、共同原因、假说组合审计统一候选；复用旧源／约束边界，不再延伸读口精度反例。[工作报告]({p}round751_drafts/joint_hypothesis_obstacle_review.md)。旧空间、604、649／699与目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第751轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 397.' not in text
        text+='\n\n## 397. 共同报告与内部来源账\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 302.' not in text
        text+='\n\n## 302. 旧空间合同保持，联合假说与完整过程审计\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—750轮。' in text and '最新科学轮次与检查数为750／3439' in text
        text=text.replace('完成231—750轮。','完成231—751轮。').replace('最新科学轮次与检查数为750／3439','最新科学轮次与检查数为751／3441')
    if p==HERE/'README.md':
        text+='\n\n## 第751轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|751|[共同报告与内部来源账](research_note_751.md)|[代码](joint_record_source_retention.py)、[结果](joint_record_source_retention_results.json)、[核验](research_round_751_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round751.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=751,cumulative_tests=3441,numbered_scientific_files=1565,
    unique_protected_evidence_files=3521,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=752,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
