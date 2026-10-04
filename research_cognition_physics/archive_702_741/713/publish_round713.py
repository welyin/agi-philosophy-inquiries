"""Publish713 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round713_navigation_checks.json'
checked=core.read(HERE/'research_round_713_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round713_20261003'
assert not folder.exists() and not TARGET.exists()
summary='**第713轮完成：** [共同尺度中的环路读取与状态电能]({p}research_note_713.md)原群平坦整环可在逐边趋单位时保对比；固定宽度平均instrument有一致注能界，但正常态的非零原字符信号另有裸电能下界，参考减除不能省略。四组、二十式通过，最新713／3343，1451份编号科学文件、2986份保护证据。[核验]({p}research_round_713_checks.json)、[全条件账]({p}unified_physics_condition_ledger_713.md)。未输送699反射型或完成连续态，目标未完成。'
order='**当前执行顺序（713后，优先于下方历史安排）：** 接[714共同参考、相对资源与来源]({p}round714_drafts/STATUS.md)，回查原完整参考与跨尺度运输，区分裸电能、参考背景和实际准备代价；几何导数与总H同步，不任意减常数。旧空间及四分支保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第713轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 359.' not in text
        text+='\n\n## 359. 原环路信息、读取预算与正常态电能\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 264.' not in text
        text+='\n\n## 264. 旧空间合同复用，读取面积不反推空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—712轮。' in text and '最新科学轮次与检查数为712／3339' in text
        text=text.replace('完成231—712轮。','完成231—713轮。').replace('最新科学轮次与检查数为712／3339','最新科学轮次与检查数为713／3343')
    if p==HERE/'README.md':
        text+='\n\n## 第713轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|713|[环路读取与状态电能](research_note_713.md)|[代码](joint_holonomy_readout_scale.py)、[结果](joint_holonomy_readout_scale_results.json)、[核验](research_round_713_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round713.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-03',latest_round=713,cumulative_tests=3343,numbered_scientific_files=1451,
    unique_protected_evidence_files=2986,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=714,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
