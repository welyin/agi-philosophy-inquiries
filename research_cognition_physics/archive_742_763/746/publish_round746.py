"""Publish746 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round746_navigation_checks.json'
checked=core.read(HERE/'research_round_746_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round746_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第746轮完成：** [原连通图中的占据读出]({p}research_note_746.md)旧Higgs幅度读口在原固定边图上有严格四阶占据差，图势和跳跃仅在此阶差值中抵消；全演化、Gauss及同态资源保留。无需弱化边或新增探针，尚非异地通信、无扰测量或自治终端。两组、十八式通过，最新746／3430，1550份编号科学文件、3449份保护证据。[核验]({p}research_round_746_checks.json)、[条件账]({p}unified_physics_condition_ledger_746.md)。'
order='**当前执行顺序（746后，优先于下方历史安排）：** 接[747原占据输入与异地记录]({p}round747_drafts/STATUS.md)，复用577通道并保未知态、实际后态及共同代价，避免新增免费控制。旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第746轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 392.' not in text
        text+='\n\n## 392. 原连通图中的占据读出\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 297.' not in text
        text+='\n\n## 297. 旧空间合同保持，原固定图读出与资源合并\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—745轮。' in text and '最新科学轮次与检查数为745／3428' in text
        text=text.replace('完成231—745轮。','完成231—746轮。').replace('最新科学轮次与检查数为745／3428','最新科学轮次与检查数为746／3430')
    if p==HERE/'README.md':
        text+='\n\n## 第746轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|746|[原连通图中的占据读出](research_note_746.md)|[代码](joint_connected_population_readout.py)、[结果](joint_connected_population_readout_results.json)、[核验](research_round_746_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round746.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=746,cumulative_tests=3430,numbered_scientific_files=1550,
    unique_protected_evidence_files=3449,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=747,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
