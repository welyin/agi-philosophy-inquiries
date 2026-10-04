"""Publish726 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round726_navigation_checks.json'
checked=core.read(HERE/'research_round_726_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round726_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第726轮完成：** [同源记录后的经典壁与完整物质涨落边界]({p}research_note_726.md)原真实记录在固定图半经典窗口共同恢复位置、法向谱及玻色来源；空CAR参考保均值却留正Majorana能源涨落，完整物质不能直接签收。三组、十八式通过，最新726／3383，1490份编号科学文件、3172份保护证据。[核验]({p}research_round_726_checks.json)、[全条件账]({p}unified_physics_condition_ledger_726.md)。原差分误差、量子连续及引力边界保留。'
order='**当前执行顺序（726后，优先于下方历史安排）：** 接[727原物质参考与空间传播]({p}round727_drafts/STATUS.md)，回查602—604、614、630—633及719—720，核费米参考、同一来源、质量与传播能否共同相容。停止波包与法向优化；旧空间、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第726轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 372.' not in text
        text+='\n\n## 372. 读后经典壁与完整物质参考\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 277.' not in text
        text+='\n\n## 277. 旧空间合同复用，玻色经典壁不等于完整来源\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—725轮。' in text and '最新科学轮次与检查数为725／3380' in text
        text=text.replace('完成231—725轮。','完成231—726轮。').replace('最新科学轮次与检查数为725／3380','最新科学轮次与检查数为726／3383')
    if p==HERE/'README.md':
        text+='\n\n## 第726轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|726|[同源记录后的经典壁与完整物质涨落边界](research_note_726.md)|[代码](joint_recorded_classical_wall.py)、[结果](joint_recorded_classical_wall_results.json)、[核验](research_round_726_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round726.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=726,cumulative_tests=3383,numbered_scientific_files=1490,
    unique_protected_evidence_files=3172,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=727,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
