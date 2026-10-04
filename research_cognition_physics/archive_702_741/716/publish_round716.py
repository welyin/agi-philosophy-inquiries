"""Publish716 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round716_navigation_checks.json'
checked=core.read(HERE/'research_round_716_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round716_20261003'
assert not folder.exists() and not TARGET.exists()
summary='**第716轮完成：** [光滑物质模的共同资源与来源]({p}research_note_716.md)原模式的局部势矩与跳跃共同控制读取、时间变化及指定几何导数；非均匀体积改变仪器，来源需同步。三组、二十式通过，最新716／3352，1460份编号科学文件、3031份保护证据。[核验]({p}research_round_716_checks.json)、[全条件账]({p}unified_physics_condition_ledger_716.md)。实际共同连续与因果实现仍开放。'
order='**当前执行顺序（716后，优先于下方历史安排）：** 接[717原记录后参考与共同传播]({p}round717_drafts/STATUS.md)，先核实际非平稳态和有限时间过程，复用623—625、632—637及667。停止波包／范数优化，旧空间和699范围保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第716轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 362.' not in text
        text+='\n\n## 362. 原物质模的体积、资源与几何共同条件\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 267.' not in text
        text+='\n\n## 267. 旧空间合同复用，波包归一不选择物理维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—715轮。' in text and '最新科学轮次与检查数为715／3349' in text
        text=text.replace('完成231—715轮。','完成231—716轮。').replace('最新科学轮次与检查数为715／3349','最新科学轮次与检查数为716／3352')
    if p==HERE/'README.md':
        text+='\n\n## 第716轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|716|[光滑物质模的共同资源与来源](research_note_716.md)|[代码](joint_smooth_mode_contract.py)、[结果](joint_smooth_mode_contract_results.json)、[核验](research_round_716_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round716.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-03',latest_round=716,cumulative_tests=3352,numbered_scientific_files=1460,
    unique_protected_evidence_files=3031,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=717,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
