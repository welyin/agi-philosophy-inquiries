"""Publish739 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round739_navigation_checks.json'
checked=core.read(HERE/'research_round_739_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round739_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第739轮完成：** [同一中性—几何退迟方程与线性约束]({p}research_note_739.md)原局部jet与完整谱组成协变方程；互补非退化、零过去相容源和固定有限空间频段下，短时线性退迟解与引力约束共同成立。两组、十八式通过，最新739／3417，1529份编号科学文件、3345份保护证据。[核验]({p}research_round_739_checks.json)、[条件账]({p}unified_physics_condition_ledger_739.md)。内部规范、有效分支、任意尺度及非线性仍开放。'
order='**当前执行顺序（739后，优先于下方历史安排）：** 接[740总反馈与共同有效分支]({p}round740_drafts/STATUS.md)，核同一总符号的极点和601有效阶次；纯谱正性及短时可解不保证物理稳定。随后接实际动态参考。旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第739轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 385.' not in text
        text+='\n\n## 385. 同一中性—几何退迟方程与有限空间频段的约束闭合\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 290.' not in text
        text+='\n\n## 290. 旧空间合同保持，线性约束不替代非线性自洽\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—738轮。' in text and '最新科学轮次与检查数为738／3415' in text
        text=text.replace('完成231—738轮。','完成231—739轮。').replace('最新科学轮次与检查数为738／3415','最新科学轮次与检查数为739／3417')
    if p==HERE/'README.md':
        text+='\n\n## 第739轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|739|[同一中性—几何退迟方程与有限空间频段的约束闭合](research_note_739.md)|[代码](joint_covariant_response_closure.py)、[结果](joint_covariant_response_closure_results.json)、[核验](research_round_739_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round739.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=739,cumulative_tests=3417,numbered_scientific_files=1529,
    unique_protected_evidence_files=3345,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=740,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
