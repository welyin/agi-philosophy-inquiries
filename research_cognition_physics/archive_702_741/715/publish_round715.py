"""Publish715 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round715_navigation_checks.json'
checked=core.read(HERE/'research_round_715_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round715_20261003'
assert not folder.exists() and not TARGET.exists()
summary='**第715轮完成：** [实际热读取与共同能源域]({p}research_note_715.md)指定环路锐化在原完整Gibbs上的平均注能为线性阶；记录／后态虽收敛，极限却有无限能源，共形来源同步发散。三组、十六式通过，最新715／3349，1457份编号科学文件、3016份保护证据。[核验]({p}research_round_715_checks.json)、[全条件账]({p}unified_physics_condition_ledger_715.md)。仅关闭固定图锐化接法，统一目标未完成。'
order='**当前执行顺序（715后，优先于下方历史安排）：** 接[716光滑物质模与同一原过程]({p}round716_drafts/STATUS.md)，复用633及666—667，核原有限图sterile读取、共同参考、质量与来源；停止锐化函数扫描，保留自由分支与因果实现限制。旧空间及699范围保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第715轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 361.' not in text
        text+='\n\n## 361. 原真实热态中的记录收敛与能源域\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 266.' not in text
        text+='\n\n## 266. 旧空间合同复用，锐化极限不代替空间连续\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—714轮。' in text and '最新科学轮次与检查数为714／3346' in text
        text=text.replace('完成231—714轮。','完成231—715轮。').replace('最新科学轮次与检查数为714／3346','最新科学轮次与检查数为715／3349')
    if p==HERE/'README.md':
        text+='\n\n## 第715轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|715|[实际热读取与共同能源域](research_note_715.md)|[代码](joint_sharp_record_domain.py)、[结果](joint_sharp_record_domain_results.json)、[核验](research_round_715_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round715.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-03',latest_round=715,cumulative_tests=3349,numbered_scientific_files=1457,
    unique_protected_evidence_files=3016,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=716,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
