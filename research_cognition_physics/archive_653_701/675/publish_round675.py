"""Publish verified 675, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round675_navigation_checks.json'
checked=core.read(HERE/'research_round_675_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round675_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第675轮完成：** [物理标量边缘、完整辅助平均与共同正性判据]({p}research_note_675.md)'
         '657精确球面多项式接入原非平坦投影与质量；实际标量边缘保留完整双Gauss平均。'
         '591旧基态给此候选的受控大τ必要正性条件；固定λ尚非原物理时间。'
         '完整积分符号、RP及原过程身份仍开放。'
         '两组、十四式通过，最新675／3257，1337份编号科学文件、2444份保护证据。'
         '[核验]({p}research_round_675_checks.json)、[全条件账]({p}unified_physics_condition_ledger_675.md)。'
         '连续及量子GR未完成。')
order=('**当前执行顺序（675后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[676完整平均后的共同物理过程]({p}round676_drafts/STATUS.md)，'
       '回查原有序过程、实际观测与时间／来源字典，不继续扩张同类抽样；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第675轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 321.' not in text
        text+='\n\n## 321. 实际标量边缘与旧基态的共同正性限制\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 226.' not in text
        text+='\n\n## 226. 非平坦辅助收缩与实际共同时间的边界\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—674轮。' in text and '最新科学轮次与检查数为674／3255' in text
        text=text.replace('完成231—674轮。','完成231—675轮。').replace('最新科学轮次与检查数为674／3255','最新科学轮次与检查数为675／3257')
    if p==HERE/'README.md':
        text+='\n\n## 第675轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|675|[物理标量边缘、完整辅助平均与共同正性判据](research_note_675.md)|[代码](joint_physical_boundary_average.py)、[结果](joint_physical_boundary_average_results.json)、[核验](research_round_675_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round675.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=675,cumulative_tests=3257,numbered_scientific_files=1337,
            unique_protected_evidence_files=2444,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=676,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
