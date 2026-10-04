"""Publish735 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round735_navigation_checks.json'
checked=core.read(HERE/'research_round_735_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round735_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第735轮完成（限定一圈框架）：** [联合守恒处方与实际来源]({p}research_note_735.md)在声明的局部UV规范化框架内，原完整变化质量经一圈规范化和有限jet匹配接到共同守恒来源；保实际过去与记录，响应仍需完整接触。两组、十六式通过，最新735／3409，1517份编号科学文件、3304份保护证据。[核验]({p}research_round_735_checks.json)、[条件账]({p}unified_physics_condition_ledger_735.md)。有限物理常数及非线性自洽解仍开放。'
order='**当前执行顺序（735后，优先于下方历史安排）：** 接[736实际反馈正则性]({p}round736_drafts/STATUS.md)，先复用旧高频／来源结果，区分局部项、非局部态响应和降阶合同；不再重复荷表或有限矩阵优化。旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第735轮完成（限定一圈框架）：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 381.' not in text
        text+='\n\n## 381. 一圈联合守恒处方的局部提升与同一实际来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 286.' not in text
        text+='\n\n## 286. 旧空间合同保持，退迟来源与接触须共同变分\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—734轮。' in text and '最新科学轮次与检查数为734／3407' in text
        text=text.replace('完成231—734轮。','完成231—735轮。').replace('最新科学轮次与检查数为734／3407','最新科学轮次与检查数为735／3409')
    if p==HERE/'README.md':
        text+='\n\n## 第735轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|735|[一圈联合守恒处方的局部提升与同一实际来源](research_note_735.md)|[代码](joint_local_source_normalization.py)、[结果](joint_local_source_normalization_results.json)、[核验](research_round_735_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round735.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=735,cumulative_tests=3409,numbered_scientific_files=1517,
    unique_protected_evidence_files=3304,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=736,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
