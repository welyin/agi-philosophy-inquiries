"""Publish verified 663, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round663_navigation_checks.json'
checked=core.read(HERE/'research_round_663_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round663_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第663轮完成：** [完整旋量的局部时间组合、空间标架与共同几何来源]({p}research_note_663.md)'
         '原连续完整物质的局部时间组合接到spin切向输送；同步几何字典恢复有限组合及来源导数。'
         '两组、十四式通过，最新663／3231，1301份编号科学文件、2299份保护证据。'
         '[核验]({p}research_round_663_checks.json)、[条件账]({p}unified_physics_condition_ledger_663.md)。'
         '给定几何与spin结构；真实图态映射及量子GR仍开放。')
order=('**当前执行顺序（663后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[664原物质的联络来源与几何变分]({p}round664_drafts/STATUS.md)，'
       '核无挠／独立联络分支、完整来源和有效相互作用；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第663轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 309.' not in text
        text+='\n\n## 309. 完整旋量的局部来源与跨背景几何组合\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 214.' not in text
        text+='\n\n## 214. 空间组合须同时输送度规、标架和来源\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—662轮。' in text and '最新科学轮次与检查数为662／3229' in text
        text=text.replace('完成231—662轮。','完成231—663轮。').replace('最新科学轮次与检查数为662／3229','最新科学轮次与检查数为663／3231')
    if p==HERE/'README.md':
        text+='\n\n## 第663轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|663|[完整旋量的局部时间组合、空间标架与共同几何来源](research_note_663.md)|[代码](joint_spinorial_geometry_sources.py)、[结果](joint_spinorial_geometry_sources_results.json)、[核验](research_round_663_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round663.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=663,cumulative_tests=3231,numbered_scientific_files=1301,
            unique_protected_evidence_files=2299,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=664,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
