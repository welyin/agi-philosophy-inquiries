"""Publish717 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round717_navigation_checks.json'
checked=core.read(HERE/'research_round_717_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round717_20261003'
assert not folder.exists() and not TARGET.exists()
summary='**第717轮完成：** [原记录的质量反馈与有界读数]({p}research_note_717.md)原质量力决定读后玻色二阶反馈；有界能源不能控制全部无界势响应，但旧sin s效果有态无关的精确二阶界。四组、二十二式通过，最新717／3356，1463份编号科学文件、3046份保护证据。[核验]({p}research_round_717_checks.json)、[全条件账]({p}unified_physics_condition_ledger_717.md)。有限时间共同历史与物理连续仍开放。'
order='**当前执行顺序（717后，优先于下方历史安排）：** 接[718原有界记录的有限时间历史]({p}round718_drafts/STATUS.md)，核同一读取的H与RHR比较，复用623—625、634及704。停止高矩反例／波包优化；旧空间、633因果及699范围保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第717轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 363.' not in text
        text+='\n\n## 363. 原记录、质量力及实际有界效果\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 268.' not in text
        text+='\n\n## 268. 旧空间合同复用，反馈系数不代替物理连续\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—716轮。' in text and '最新科学轮次与检查数为716／3352' in text
        text=text.replace('完成231—716轮。','完成231—717轮。').replace('最新科学轮次与检查数为716／3352','最新科学轮次与检查数为717／3356')
    if p==HERE/'README.md':
        text+='\n\n## 第717轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|717|[原记录的质量反馈与有界读数](research_note_717.md)|[代码](joint_record_mass_feedback.py)、[结果](joint_record_mass_feedback_results.json)、[核验](research_round_717_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round717.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-03',latest_round=717,cumulative_tests=3356,numbered_scientific_files=1463,
    unique_protected_evidence_files=3046,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=718,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
