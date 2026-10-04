"""Publish712 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round712_navigation_checks.json'
checked=core.read(HERE/'research_round_712_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round712_20261003'
assert not folder.exists() and not TARGET.exists()
summary='**第712轮完成：** [原插入的共同演化与记录方向]({p}research_note_712.md)同一原H固定质量、跨边及Majorana伴随通道，并约束完整热谱；原singlet读口二阶消去，Higgs比较量非零且共用几何来源。三组、十八式通过，最新712／3339，1448份编号科学文件、2971份保护证据。[核验]({p}research_round_712_checks.json)、[全条件账]({p}unified_physics_condition_ledger_712.md)。没有生成耦合或拓扑来源，统一目标未完成。'
order='**当前执行顺序（712后，优先于下方历史安排）：** 接[713复合CAR与共同物理来源]({p}round713_drafts/STATUS.md)，回查688及已有Weyl来源字典，核实际插入、接触和时间拼接；不将代数映射当699正性恢复，不扫描高阶系数。旧空间及统一目标保留。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第712轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 358.' not in text
        text+='\n\n## 358. 原复合插入的共同演化、热谱及实际读口\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 263.' not in text
        text+='\n\n## 263. 旧空间合同复用，场空间记录方向不当空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—711轮。' in text and '最新科学轮次与检查数为711／3336' in text
        text=text.replace('完成231—711轮。','完成231—712轮。').replace('最新科学轮次与检查数为711／3336','最新科学轮次与检查数为712／3339')
    if p==HERE/'README.md':
        text+='\n\n## 第712轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|712|[原插入的共同演化与记录方向](research_note_712.md)|[代码](joint_vertex_shared_evolution.py)、[结果](joint_vertex_shared_evolution_results.json)、[核验](research_round_712_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round712.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-03',latest_round=712,cumulative_tests=3339,numbered_scientific_files=1448,
    unique_protected_evidence_files=2971,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=713,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
