"""Publish711 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round711_navigation_checks.json'
checked=core.read(HERE/'research_round_711_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round711_20262003'
assert not folder.exists() and not TARGET.exists()
summary='**第711轮完成：** [平坦历史相位与弱拓扑匹配边界]({p}research_note_711.md)原完整构形中的标量平坦相位可接共同H、新热参考、原记录和几何源；全部此类相位来自阿贝尔图环路，对纯弱历史平凡且保Nq，不能按原字典补成弱拓扑反常。三组、二十式通过，最新711／3336，1445份编号科学文件、2957份保护证据。[核验]({p}research_round_711_checks.json)、[全条件账]({p}unified_physics_condition_ledger_711.md)。限定路线结论，统一目标未完成。'
order='**当前执行顺序（711后，优先于下方历史安排）：** 接[712拓扑域、费米子与共同过程]({p}round712_drafts/STATUS.md)，回查原Wilson零集合、可容许域和物理费米子字典；改域或拼接时共同验热参考、时间正性和来源。停止平坦角扫描，认知约束继续协调条件，旧空间及699范围保留。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第711轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 357.' not in text
        text+='\n\n## 357. 原平坦历史相位与共同过程，弱拓扑匹配受限\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 262.' not in text
        text+='\n\n## 262. 旧空间合同复用，图环路不当物理拓扑数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—710轮。' in text and '最新科学轮次与检查数为710／3333' in text
        text=text.replace('完成231—710轮。','完成231—711轮。').replace('最新科学轮次与检查数为710／3333','最新科学轮次与检查数为711／3336')
    if p==HERE/'README.md':
        text+='\n\n## 第711轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|711|[平坦历史相位与弱拓扑边界](research_note_711.md)|[代码](joint_flat_history_phase.py)、[结果](joint_flat_history_phase_results.json)、[核验](research_round_711_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round711.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-03',latest_round=711,cumulative_tests=3336,numbered_scientific_files=1445,
    unique_protected_evidence_files=2957,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=712,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
