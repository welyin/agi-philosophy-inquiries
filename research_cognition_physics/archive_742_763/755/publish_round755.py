"""Publish755 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round755_navigation_checks.json'
checked=core.read(HERE/'research_round_755_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round755_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第755轮完成：** [共同物理态、来源与联合记录]({p}research_note_755.md)原闭Gauss全部夸克态的忠实准自由替换受严格支撑限制；两模式填充而其余态不变的接法，至少一份边缘误差≥1/2。原完整singlet态可保瞬时二次来源，却改变联合记录；未完成连续或量子引力桥。三组、十六式通过，最新755／3452，1577份编号科学文件、3562份保护证据。[核验]({p}research_round_755_checks.json)、[条件账]({p}unified_physics_condition_ledger_755.md)。'
order='**当前执行顺序（755后，优先于下方历史安排）：** 按[认知观察与共同候选]({p}round755_drafts/cognitive_joint_candidate_working_report.md)，接[756同一关联态与连续来源]({p}round756_drafts/STATUS.md)，保原Gauss态、实际记录和全部来源，不以Gaussian摘要替换整个过程；不继续中心标签或谱精度支线。[范围审计]({p}round755_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及应用目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第755轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 401.' not in text
        text+='\n\n## 401. 共同物理态、来源与联合记录\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 306.' not in text
        text+='\n\n## 306. 旧空间合同保持，Gauss物理态与来源任务\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—754轮。' in text and '最新科学轮次与检查数为754／3449' in text
        text=text.replace('完成231—754轮。','完成231—755轮。').replace('最新科学轮次与检查数为754／3449','最新科学轮次与检查数为755／3452')
    if p==HERE/'README.md':
        text+='\n\n## 第755轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|755|[共同物理态、来源与联合记录](research_note_755.md)|[代码](joint_gaussian_physical_state_bridge.py)、[结果](joint_gaussian_physical_state_bridge_results.json)、[核验](research_round_755_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round755.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=755,cumulative_tests=3452,numbered_scientific_files=1577,
    unique_protected_evidence_files=3562,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=756,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
