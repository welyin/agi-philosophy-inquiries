"""Publish721 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round721_navigation_checks.json'
checked=core.read(HERE/'research_round_721_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round721_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第721轮完成：** [原物质速度读口的完整能源与互补质量反馈]({p}research_note_721.md)原参考流保完整固定图能源形式域，实际两结果与原等待形成有限能源历史；s／T速度对Dirac／Majorana有互补反馈。固定instrument来源保留，变化读口的总导数另明域。三组、十六式通过，最新721／3368，1475份编号科学文件、3102份保护证据。[核验]({p}research_round_721_checks.json)、[全条件账]({p}unified_physics_condition_ledger_721.md)。自主装置、空间细化及统一目标开放。'
order='**当前执行顺序（721后，优先于下方历史安排）：** 接[722导数平方参考与关系任务]({p}round722_drafts/STATUS.md)，复用554、652及704—708，核完整参考、实际记录与方向任务的同一性。停止速度函数及常数优化，工程设计后置；旧空间及699范围保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第721轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 367.' not in text
        text+='\n\n## 367. 原速度记录、完整能源与互补质量\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 272.' not in text
        text+='\n\n## 272. 旧空间合同复用，有限能源不自动给联合坐标\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—720轮。' in text and '最新科学轮次与检查数为720／3365' in text
        text=text.replace('完成231—720轮。','完成231—721轮。').replace('最新科学轮次与检查数为720／3365','最新科学轮次与检查数为721／3368')
    if p==HERE/'README.md':
        text+='\n\n## 第721轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|721|[原物质速度读口的完整能源与互补质量反馈](research_note_721.md)|[代码](joint_velocity_record_energy.py)、[结果](joint_velocity_record_energy_results.json)、[核验](research_round_721_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round721.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=721,cumulative_tests=3368,numbered_scientific_files=1475,
    unique_protected_evidence_files=3102,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=722,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
