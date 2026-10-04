"""Publish709 only after scientific verification; preserve navigation snapshots."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core

HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round709_navigation_checks.json'
checked=core.read(HERE/'research_round_709_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round709_20261003'
assert not folder.exists() and not TARGET.exists()
summary='**第709轮完成：** [荷粗化、跨区相位与反常接口]({p}research_note_709.md)原完整H有保Gauss、热参考、中性真实历史和能源的正常量子条件期望；逐区独立使用会删原跨边流，629反常限制全U(1)外推。四组、十八式通过，最新709／3330，1439份编号科学文件、2929份保护证据。[核验]({p}research_round_709_checks.json)、[全条件账]({p}unified_physics_condition_ledger_709.md)。未删空间变量，旧空间与统一目标不变。'
order='**当前执行顺序（709后，优先于下方历史安排）：** 以少量认知假说共同约束缺口，具体工程设计后置。接[710拓扑相位与区域参考]({p}round710_drafts/STATUS.md)，核原规范中心、残余物理相位和跨区相干能否同时保留；四分支及699范围继续区分，不把有限图守恒当全尺度公理。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第709轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 355.' not in text
        text+='\n\n## 355. 原量子荷粗化、区域输运及反常共同验收\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 260.' not in text
        text+='\n\n## 260. 旧空间接口保留，有限荷对称不当连续生成\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—708轮。' in text and '最新科学轮次与检查数为708／3326' in text
        text=text.replace('完成231—708轮。','完成231—709轮。').replace('最新科学轮次与检查数为708／3326','最新科学轮次与检查数为709／3330')
    if p==HERE/'README.md':
        text+='\n\n## 第709轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|709|[荷粗化、跨区相位与反常接口](research_note_709.md)|[代码](joint_charge_quantum_coarse.py)、[结果](joint_charge_quantum_coarse_results.json)、[核验](research_round_709_checks.json)|\n'
    planned[p]=text.replace('\n',nl).encode(enc)
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
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round709.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-03',latest_round=709,cumulative_tests=3330,numbered_scientific_files=1439,
    unique_protected_evidence_files=2929,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=710,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
