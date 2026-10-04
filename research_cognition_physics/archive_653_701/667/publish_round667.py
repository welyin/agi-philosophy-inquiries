"""Publish verified 667, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round667_navigation_checks.json'
checked=core.read(HERE/'research_round_667_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round667_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第667轮完成：** [局部参考、空间传播与共同量子过程的相容条件]({p}research_note_667.md)'
         '局部参考字典保区域和原过程，空间传播、质量及几何来源须共同输送；无限空Fock实现并非必要。'
         '两组、十六式通过，最新667／3241，1313份编号科学文件、2342份保护证据。'
         '[核验]({p}research_round_667_checks.json)、[全条件账]({p}unified_physics_condition_ledger_667.md)。'
         '明确复用382—386、425、522—523；旧条件性空间成果保留，实际共同尺度映射与量子GR仍开放。')
order=('**当前执行顺序（667后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[668原非零质量与手征辅助测度]({p}round668_drafts/STATUS.md)，'
       '核同一物理观测、权重、反射和来源，复用既有空间证明；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第667轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 313.' not in text
        text+='\n\n## 313. 局部参考与原空间传播的共同连接\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 218.' not in text
        text+='\n\n## 218. 复用旧维数与坐标证明，核实际物质和尺度映射\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—666轮。' in text and '最新科学轮次与检查数为666／3239' in text
        text=text.replace('完成231—666轮。','完成231—667轮。').replace('最新科学轮次与检查数为666／3239','最新科学轮次与检查数为667／3241')
    if p==HERE/'README.md':
        text+='\n\n## 第667轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|667|[局部参考、空间传播与共同量子过程的相容条件](research_note_667.md)|[代码](joint_reference_spatial_process.py)、[结果](joint_reference_spatial_process_results.json)、[核验](research_round_667_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round667.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=667,cumulative_tests=3241,numbered_scientific_files=1313,
            unique_protected_evidence_files=2342,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=668,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
