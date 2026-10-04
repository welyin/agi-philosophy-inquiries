"""Publish verified 651, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round651_navigation_checks.json'
checked=core.read(HERE/'research_round_651_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round651_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第651轮完成：** [原物质参考、因果区域与完整边界变分的共同连接]({p}research_note_651.md)'
         '同一原解的固定s面类空，原时钟组合给局部类时侧壁；'
         '原参考约束关系度规与F，同一嵌入输送完整边界及透射匹配。'
         '三组、十六式通过，最新651／3191，1265份编号科学文件、2175份保护证据。'
         '[核验]({p}research_round_651_checks.json)、[条件账]({p}unified_physics_condition_ledger_651.md)。'
         '限给定作用的局部经典接口；诊断族离壳，量子参考及全尺度统一仍开放。')
order=('**当前执行顺序（651后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[652原关系参考、完整量子过程与共同区域表示]({p}round652_drafts/STATUS.md)，'
       '回查原量子域、来源和参考的真实映射，不继续侧壁或坐标精度优化；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第651轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 297.' not in text
        text+='\n\n## 297. 原物质图册与共同边界输送\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 202.' not in text
        text+='\n\n## 202. 物质标签须与因果边界共同相容\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—650轮。' in text and '最新科学轮次与检查数为650／3188' in text
        text=text.replace('完成231—650轮。','完成231—651轮。').replace('最新科学轮次与检查数为650／3188','最新科学轮次与检查数为651／3191')
    if p==HERE/'README.md':
        text+='\n\n## 第651轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|651|[原物质参考、因果区域与完整边界变分的共同连接](research_note_651.md)|[代码](joint_material_boundary_transport.py)、[结果](joint_material_boundary_transport_results.json)、[核验](research_round_651_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round651.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=651,cumulative_tests=3191,numbered_scientific_files=1265,
            unique_protected_evidence_files=2175,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=652,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
