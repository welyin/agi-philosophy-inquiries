"""Publish verified 664, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round664_navigation_checks.json'
checked=core.read(HERE/'research_round_664_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round664_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第664轮完成：** [原物质、联络选择与共同几何来源的匹配条件]({p}research_note_664.md)'
         '原Jordan最小独立联络同时改变标量目标和总轴流作用；指定匹配项恢复原经典公共作用。'
         '原CAR核共同接触、排序及几何来源。三组、十六式通过，最新664／3234，1304份编号科学文件、2306份保护证据。'
         '[核验]({p}research_round_664_checks.json)、[条件账]({p}unified_physics_condition_ledger_664.md)。'
         '保留无挠基线；量子测度、连续映射及量子GR仍开放。')
order=('**当前执行顺序（664后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[665同一CAR过程与联络辅助权重]({p}round665_drafts/STATUS.md)，'
       '核接触作用、有序权重、辅助测度及共同来源；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第664轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 310.' not in text
        text+='\n\n## 310. 原物质与联络分支的共同匹配\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 215.' not in text
        text+='\n\n## 215. 联络选择同时约束标量几何、物质与来源\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—663轮。' in text and '最新科学轮次与检查数为663／3231' in text
        text=text.replace('完成231—663轮。','完成231—664轮。').replace('最新科学轮次与检查数为663／3231','最新科学轮次与检查数为664／3234')
    if p==HERE/'README.md':
        text+='\n\n## 第664轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|664|[原物质、联络选择与共同几何来源的匹配条件](research_note_664.md)|[代码](joint_connection_matter_matching.py)、[结果](joint_connection_matter_matching_results.json)、[核验](research_round_664_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round664.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=664,cumulative_tests=3234,numbered_scientific_files=1304,
            unique_protected_evidence_files=2306,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=665,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
