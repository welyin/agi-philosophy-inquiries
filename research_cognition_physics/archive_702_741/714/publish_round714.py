"""Publish714 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round714_navigation_checks.json'
checked=core.read(HERE/'research_round_714_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round714_20261003'
assert not folder.exists() and not TARGET.exists()
summary='**第714轮完成：** [中心信息、带荷关联与共同热参考]({p}research_note_714.md)原中心变换的完整H差量仅含带荷边项；原Gibbs平均代价严格正，实际记录、热关联及几何来源受同一约束，来源导数须含参考响应。三组、十六式通过，最新714／3346，1454份编号科学文件、3002份保护证据。[核验]({p}research_round_714_checks.json)、[全条件账]({p}unified_physics_condition_ledger_714.md)。未计算原热环路均值或完成共同连续，目标未完成。'
order='**当前执行顺序（714后，优先于下方历史安排）：** 接[715原热环路信息与共同尺度]({p}round715_drafts/STATUS.md)，核实际观测、原参考与物理尺度映射；不把固定图严格正性当统一细化界，不重复被动性。旧空间及四分支保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第714轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 360.' not in text
        text+='\n\n## 360. 原热参考、带荷差量与几何响应\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 265.' not in text
        text+='\n\n## 265. 旧空间合同复用，热差量不反推空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—713轮。' in text and '最新科学轮次与检查数为713／3343' in text
        text=text.replace('完成231—713轮。','完成231—714轮。').replace('最新科学轮次与检查数为713／3343','最新科学轮次与检查数为714／3346')
    if p==HERE/'README.md':
        text+='\n\n## 第714轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|714|[中心信息与共同热参考](research_note_714.md)|[代码](joint_center_thermal_reference.py)、[结果](joint_center_thermal_reference_results.json)、[核验](research_round_714_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round714.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-03',latest_round=714,cumulative_tests=3346,numbered_scientific_files=1454,
    unique_protected_evidence_files=3002,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=715,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
