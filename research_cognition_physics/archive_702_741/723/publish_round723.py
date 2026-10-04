"""Publish723 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round723_navigation_checks.json'
checked=core.read(HERE/'research_round_723_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round723_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第723轮完成：** [完整物质参考记录与同一有限量子过程]({p}research_note_723.md)原Gauss态可保四参考全部边缘与平均H而改变实际顺序记录；新参考菜单的完整有限过程保记录后态、能源和末端来源，立即联合几何一阶响应也共同收敛。三组、十八式通过，最新723／3374，1481份编号科学文件、3130份保护证据。[核验]({p}research_round_723_checks.json)、[全条件账]({p}unified_physics_condition_ledger_723.md)。变化等待高阶、关系区域和空间极限仍分范围。'
order='**当前执行顺序（723后，优先于下方历史安排）：** 接[724实际参考与关系区域]({p}round724_drafts/STATUS.md)，回查647—652原选面、因果类型和完整边界来源，使用同一参考过程核实际对象映射。停止单读口、尾部与低秩优化；旧空间及699范围保持，统一目标开放。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第723轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 369.' not in text
        text+='\n\n## 369. 完整物质参考、实际历史与来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 274.' not in text
        text+='\n\n## 274. 旧空间合同复用，完整内部过程仍须接实际区域\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—722轮。' in text and '最新科学轮次与检查数为722／3371' in text
        text=text.replace('完成231—722轮。','完成231—723轮。').replace('最新科学轮次与检查数为722／3371','最新科学轮次与检查数为723／3374')
    if p==HERE/'README.md':
        text+='\n\n## 第723轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|723|[完整物质参考记录与同一有限量子过程](research_note_723.md)|[代码](joint_reference_process_transport.py)、[结果](joint_reference_process_transport_results.json)、[核验](research_round_723_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round723.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=723,cumulative_tests=3374,numbered_scientific_files=1481,
    unique_protected_evidence_files=3130,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=724,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
