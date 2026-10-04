"""Publish710 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round710_navigation_checks.json'
checked=core.read(HERE/'research_round_710_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round710_20261003'
assert not folder.exists() and not TARGET.exists()
summary='**第710轮完成：** [允许的荷变化、共同参考与拓扑相位边界]({p}research_note_710.md)原CAR／Gauss上的明确QQQL扩展可承接荷变化、原区域流、读口预算和新热参考；相位及体积源受共同约束。但单一系数相位仍可消去，不等于反常生成。三组、十八式通过，最新710／3333，1442份编号科学文件、2943份保护证据。[核验]({p}research_round_710_checks.json)、[全条件账]({p}unified_physics_condition_ledger_710.md)。新增相互作用及匹配明确记为输入；目标未完成。'
order='**当前执行顺序（710后，优先于下方历史安排）：** 以少量认知假说共同约束缺口，具体工程设计后置。接[711独立拓扑贡献与同一过程]({p}round711_drafts/STATUS.md)，先核原规范构形、区域拼接和实际来源中的拓扑载体；不以追加拟合相位替代生成。旧空间、四分支及699范围保留。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第710轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 356.' not in text
        text+='\n\n## 356. 荷变化与共同来源，有限相位不冒充拓扑生成\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 261.' not in text
        text+='\n\n## 261. 旧空间合同复用，新增体积匹配须同源\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—709轮。' in text and '最新科学轮次与检查数为709／3330' in text
        text=text.replace('完成231—709轮。','完成231—710轮。').replace('最新科学轮次与检查数为709／3330','最新科学轮次与检查数为710／3333')
    if p==HERE/'README.md':
        text+='\n\n## 第710轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|710|[荷变化与拓扑相位边界](research_note_710.md)|[代码](joint_charge_changing_vertex.py)、[结果](joint_charge_changing_vertex_results.json)、[核验](research_round_710_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round710.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-03',latest_round=710,cumulative_tests=3333,numbered_scientific_files=1442,
    unique_protected_evidence_files=2943,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=711,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
