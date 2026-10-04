"""Publish724 with exclusive navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round724_navigation_checks.json'
checked=core.read(HERE/'research_round_724_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round724_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第724轮完成：** [原物质选区、规范拼接与移动记录来源]({p}research_note_724.md)原精确物质壁经同次有限读数选择区域；记录合并与全Gauss拼接交换，保原能源及来源；移动记录边界须计分支通量。三组、十八式通过，最新724／3377，1484份编号科学文件、3144份保护证据。[核验]({p}research_round_724_checks.json)、[全条件账]({p}unified_physics_condition_ledger_724.md)。未完成法向、空间细化或量子引力拼接。'
order='**当前执行顺序（724后，优先于下方历史安排）：** 接[725实际区域与原壁法向]({p}round725_drafts/STATUS.md)，回查563及651—652的导数域，核原h壁位置与法向资料能否共同实现。停止区间、阈值和读口优化；旧空间、649／699范围及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第724轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 370.' not in text
        text+='\n\n## 370. 原物质区域、规范组织与移动来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 275.' not in text
        text+='\n\n## 275. 旧空间合同复用，实际区域仍须接法向与因果类型\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—723轮。' in text and '最新科学轮次与检查数为723／3374' in text
        text=text.replace('完成231—723轮。','完成231—724轮。').replace('最新科学轮次与检查数为723／3374','最新科学轮次与检查数为724／3377')
    if p==HERE/'README.md':
        text+='\n\n## 第724轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|724|[原物质选区、规范拼接与移动记录来源](research_note_724.md)|[代码](joint_material_region_records.py)、[结果](joint_material_region_records_results.json)、[核验](research_round_724_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round724.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=724,cumulative_tests=3377,numbered_scientific_files=1484,
    unique_protected_evidence_files=3144,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=725,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
