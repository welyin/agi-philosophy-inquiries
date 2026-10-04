"""Publish744 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round744_navigation_checks.json'
checked=core.read(HERE/'research_round_744_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round744_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第744轮完成：** [原读口与联合历史]({p}research_note_744.md)原全H的联合对称性严格限制单次读取；两次原读取的偶报告允许占据信息。原一节点完整量子H给六阶差约−0.55005的数值证据，严格符号、有限时及连通图结论仍待证。两组、十六式通过，最新744／3426，1544份编号科学文件、3421份保护证据。[核验]({p}research_round_744_checks.json)、[条件账]({p}unified_physics_condition_ledger_744.md)。'
order='**当前执行顺序（744后，优先于下方历史安排）：** 接[745联合历史的严格信号]({p}round745_drafts/STATUS.md)，核原数值信号的误差来源和有限时接口；随后接正等待及连通图，保完整后态与来源。旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第744轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 390.' not in text
        text+='\n\n## 390. 原读口与联合历史\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 295.' not in text
        text+='\n\n## 295. 旧空间合同保持，单次与联合读取范围分开\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—743轮。' in text and '最新科学轮次与检查数为743／3424' in text
        text=text.replace('完成231—743轮。','完成231—744轮。').replace('最新科学轮次与检查数为743／3424','最新科学轮次与检查数为744／3426')
    if p==HERE/'README.md':
        text+='\n\n## 第744轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|744|[原读口与联合历史](research_note_744.md)|[代码](joint_native_history_readout.py)、[结果](joint_native_history_readout_results.json)、[核验](research_round_744_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round744.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=744,cumulative_tests=3426,numbered_scientific_files=1544,
    unique_protected_evidence_files=3421,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=745,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
