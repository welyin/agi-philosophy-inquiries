"""Publish verified 659, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round659_navigation_checks.json'
checked=core.read(HERE/'research_round_659_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round659_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第659轮完成：** [原手征权重、局部边界表示与共同来源的连接]({p}research_note_659.md)'
         '原带相位辅助权重接到指定domain-wall真实边界；受控双极限保留整个球面平均。'
         '边界变换与体减除来源不可另选或删去。'
         '四组、十六式通过，最新659／3221，1289份编号科学文件、2256份保护证据。'
         '[核验]({p}research_round_659_checks.json)、[全条件账]({p}unified_physics_condition_ledger_659.md)。'
         '共同物理时间、原CAR、一般规范场及量子GR仍开放。')
order=('**当前执行顺序（659后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[660原测度与共同物理时间]({p}round660_drafts/STATUS.md)，'
       '核成熟反射／Hamiltonian结果的实际映射；不优化第五方向层数或装置，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第659轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 305.' not in text
        text+='\n\n## 305. 原测度的真实边界连接\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 210.' not in text
        text+='\n\n## 210. 原测度、边界与来源须共用同一映射\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—658轮。' in text and '最新科学轮次与检查数为658／3217' in text
        text=text.replace('完成231—658轮。','完成231—659轮。').replace('最新科学轮次与检查数为658／3217','最新科学轮次与检查数为659／3221')
    if p==HERE/'README.md':
        text+='\n\n## 第659轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|659|[原手征权重、局部边界表示与共同来源的连接](research_note_659.md)|[代码](joint_domain_wall_boundary_mapping.py)、[结果](joint_domain_wall_boundary_mapping_results.json)、[核验](research_round_659_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round659.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=659,cumulative_tests=3221,numbered_scientific_files=1289,
            unique_protected_evidence_files=2256,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=660,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
