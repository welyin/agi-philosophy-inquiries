"""Publish verified 662, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round662_navigation_checks.json'
checked=core.read(HERE/'research_round_662_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round662_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第662轮完成：** [完整费米物质、局部时间来源与原记录的共同接口]({p}research_note_662.md)'
         '原完整固定图的正局部lapse族共用算符域、Gauss与原记录；完整能源流保留电微分及CAR矩阵补项。'
         '两组、十四式通过，最新662／3229，1298份编号科学文件、2289份保护证据。'
         '[核验]({p}research_round_662_checks.json)、[条件账]({p}unified_physics_condition_ledger_662.md)。'
         '局部分配、连续手征映射及量子GR仍开放。')
order=('**当前执行顺序（662后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[663局部来源与连续手征、几何的共同字典]({p}round663_drafts/STATUS.md)，'
       '核同一传播、态及几何来源的实际映射；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第662轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 308.' not in text
        text+='\n\n## 308. 完整局部时间来源与原记录的共同域\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 213.' not in text
        text+='\n\n## 213. 完整费米来源不能由纯玻色流代替\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—661轮。' in text and '最新科学轮次与检查数为661／3227' in text
        text=text.replace('完成231—661轮。','完成231—662轮。').replace('最新科学轮次与检查数为661／3227','最新科学轮次与检查数为662／3229')
    if p==HERE/'README.md':
        text+='\n\n## 第662轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|662|[完整费米物质、局部时间来源与原记录的共同接口](research_note_662.md)|[代码](joint_local_lapse_fermion_current.py)、[结果](joint_local_lapse_fermion_current_results.json)、[核验](research_round_662_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round662.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=662,cumulative_tests=3229,numbered_scientific_files=1298,
            unique_protected_evidence_files=2289,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=663,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
