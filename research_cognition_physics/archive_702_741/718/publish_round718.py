"""Publish718 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round718_navigation_checks.json'
checked=core.read(HERE/'research_round_718_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round718_20261003'
assert not folder.exists() and not TARGET.exists()
summary='**第718轮完成：** [保留真实记录的有限时间与几何输送]({p}research_note_718.md)同一局部资源控制首次／后续记录和末态的有限时间差异；移动仪器、交叉端点及热准备响应须同步保留。三组、二十式通过，最新718／3359，1466份编号科学文件、3060份保护证据。[核验]({p}research_round_718_checks.json)、[全条件账]({p}unified_physics_condition_ledger_718.md)。共同空间极限及因果实现仍开放。'
order='**当前执行顺序（718后，优先于下方历史安排）：** 接[719真实记录与区域因果]({p}round719_drafts/STATUS.md)，复用524／633／634／667及区域Gauss合同，核同一过程的准确度与因果限制。停止一般距离、反射和高矩扫描；工程设计后置，旧空间及699范围保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第718轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 364.' not in text
        text+='\n\n## 364. 原真实记录、有限时间及几何来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 269.' not in text
        text+='\n\n## 269. 旧空间合同复用，有限时间界不代替光锥\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—717轮。' in text and '最新科学轮次与检查数为717／3356' in text
        text=text.replace('完成231—717轮。','完成231—718轮。').replace('最新科学轮次与检查数为717／3356','最新科学轮次与检查数为718／3359')
    if p==HERE/'README.md':
        text+='\n\n## 第718轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|718|[保留真实记录的有限时间与几何输送](research_note_718.md)|[代码](joint_record_history_transport.py)、[结果](joint_record_history_transport_results.json)、[核验](research_round_718_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round718.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-03',latest_round=718,cumulative_tests=3359,numbered_scientific_files=1466,
    unique_protected_evidence_files=3060,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=719,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
