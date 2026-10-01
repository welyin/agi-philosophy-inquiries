"""Publish verified 618, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round618_navigation_checks.json'
checked=core.read(HERE/'research_round_618_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round618_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第618轮完成：** [原非最小耦合、几何边界与标量来源的共同匹配]({p}research_note_618.md)'
         '原F同时固定几何与五标量边界动量，混合Schur量恒M；'
         '无源类时拼接须同时匹配法向数据，框架变换须保边界标量来源。'
         '三组、十四式通过，最新618／3098，1166份编号科学文件、1922份保护证据。'
         '[核验]({p}research_round_618_checks.json)、[条件账]({p}unified_physics_condition_ledger_618.md)。'
         '限给定二导数经典作用；全物质边界、尺度及量子GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（618后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[619原物质边界电流与共同拼接]({p}round619_drafts/STATUS.md)，'
       '回查旋量与手征范围，再核原规范／费米边界及来源；目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第618轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 264.' not in text
        text+='\n\n## 264. 原F、共同边界动量与来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 169.' not in text
        text+='\n\n## 169. 无源接合的条件性匹配不选择维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—617轮。' in text and '最新科学轮次与检查数为617／3095' in text
        text=text.replace('完成231—617轮。','完成231—618轮。').replace('最新科学轮次与检查数为617／3095','最新科学轮次与检查数为618／3098')
    if p==HERE/'README.md':
        text+='\n\n## 第618轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|618|[原非最小耦合、几何边界与标量来源的共同匹配](research_note_618.md)|[代码](joint_geometric_boundary_matching.py)、[结果](joint_geometric_boundary_matching_results.json)、[核验](research_round_618_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round618.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=618,cumulative_tests=3098,numbered_scientific_files=1166,
            unique_protected_evidence_files=1922,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=619,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
