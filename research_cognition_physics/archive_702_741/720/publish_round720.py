"""Publish720 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round720_navigation_checks.json'
checked=core.read(HERE/'research_round_720_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round720_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第720轮完成：** [传播中的内部参考、原质量与实际记录]({p}research_note_720.md)原两质量允许自由均匀分支的FW守恒自旋；局部字典、实际记录与几何来源须同步改变。全部原径向质量背景下，共同平移不变内部参考不能承载SU2。三组、十六式通过，最新720／3365，1472份编号科学文件、3088份保护证据。[核验]({p}research_round_720_checks.json)、[全条件账]({p}unified_physics_condition_ledger_720.md)。结论限于明确参考类别，统一目标开放。'
order='**当前执行顺序（720后，优先于下方历史安排）：** 接[721原物质坐标与方向记录]({p}round721_drafts/STATUS.md)，复用548—554、573—574、647—651及523／524，核真正缺失的同一量子参考与instrument映射。停止保护spin分类，工程设计后置；旧空间及699范围保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第720轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 366.' not in text
        text+='\n\n## 366. 传播参考、原质量与实际记录\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 271.' not in text
        text+='\n\n## 271. 旧空间合同复用，FW参考不自动是原局部轴\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—719轮。' in text and '最新科学轮次与检查数为719／3362' in text
        text=text.replace('完成231—719轮。','完成231—720轮。').replace('最新科学轮次与检查数为719／3362','最新科学轮次与检查数为720／3365')
    if p==HERE/'README.md':
        text+='\n\n## 第720轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|720|[传播中的内部参考、原质量与实际记录](research_note_720.md)|[代码](joint_momentum_spin_reference.py)、[结果](joint_momentum_spin_reference_results.json)、[核验](research_round_720_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round720.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=720,cumulative_tests=3365,numbered_scientific_files=1472,
    unique_protected_evidence_files=3088,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=721,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
