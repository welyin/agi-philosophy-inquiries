"""Publish verified 669, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round669_navigation_checks.json'
checked=core.read(HERE/'research_round_669_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round669_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第669轮完成：** [原动态标量、手征物质与共同规范过程的条件]({p}research_note_669.md)'
         '原热过程控制全部有限时刻标量多项式；动态质量候选可积并有未归一反射正性。'
         '保自由费米动能只积分原玻色变量，不能恢复局部Gauss；质量和链路须共同匹配。'
         '两组、十八式通过，最新669／3245，1319份编号科学文件、2364份保护证据。'
         '[核验]({p}research_round_669_checks.json)、[全条件账]({p}unified_physics_condition_ledger_669.md)。'
         '指定归一、完整规范过程、连续及量子GR仍开放。')
order=('**当前执行顺序（669后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[670真实规范链路与手征物质]({p}round670_drafts/STATUS.md)，'
       '核非平坦背景的测度、质量和观测共同字典，再接动态正性与记录；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第669轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 315.' not in text
        text+='\n\n## 315. 原动态标量与规范传播的联合条件\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 220.' not in text
        text+='\n\n## 220. 原热标量可积不替代共同规范动力学\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—668轮。' in text and '最新科学轮次与检查数为668／3243' in text
        text=text.replace('完成231—668轮。','完成231—669轮。').replace('最新科学轮次与检查数为668／3243','最新科学轮次与检查数为669／3245')
    if p==HERE/'README.md':
        text+='\n\n## 第669轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|669|[原动态标量、手征物质与共同规范过程的条件](research_note_669.md)|[代码](joint_dynamic_scalar_gauge_interface.py)、[结果](joint_dynamic_scalar_gauge_interface_results.json)、[核验](research_round_669_checks.json)|\n'
    planned[p]=text.replace('\n',newline).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir(exist_ok=False)
manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as stream:stream.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as stream:json.dump(manifest,stream,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round669.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=669,cumulative_tests=3245,numbered_scientific_files=1319,
            unique_protected_evidence_files=2364,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=670,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
