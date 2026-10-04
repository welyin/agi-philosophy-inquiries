"""Publish750 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round750_navigation_checks.json'
checked=core.read(HERE/'research_round_750_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round750_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第750轮完成：** [实际标记来源与共同概率]({p}research_note_750.md)原条件后态固定完整来源和概率权重响应；分别合法的背景/准备分支仍须满足同一instrument归一。原历史试接出现非零概率缺陷，未宣称实际GR分支失败。三组、十六式通过，最新750／3439，1562份编号科学文件、3506份保护证据。[核验]({p}research_round_750_checks.json)、[条件账]({p}unified_physics_condition_ledger_750.md)。'
order='**当前执行顺序（750后，优先于下方历史安排）：** 接[751共同过去与未来响应]({p}round751_drafts/STATUS.md)，保原未知输入及共同记录权重，核来源/初值/未来的同一过程；不把条件解族当真实实施。统一目标和旧空间、604、649／699保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第750轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 396.' not in text
        text+='\n\n## 396. 实际标记来源与共同概率\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 301.' not in text
        text+='\n\n## 301. 旧空间合同保持，原条件源与共同概率相容\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—749轮。' in text and '最新科学轮次与检查数为749／3436' in text
        text=text.replace('完成231—749轮。','完成231—750轮。').replace('最新科学轮次与检查数为749／3436','最新科学轮次与检查数为750／3439')
    if p==HERE/'README.md':
        text+='\n\n## 第750轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|750|[实际标记来源与共同概率](research_note_750.md)|[代码](joint_marked_source_completion.py)、[结果](joint_marked_source_completion_results.json)、[核验](research_round_750_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round750.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=750,cumulative_tests=3439,numbered_scientific_files=1562,
    unique_protected_evidence_files=3506,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=751,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
