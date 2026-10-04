"""Publish741 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round741_navigation_checks.json'
checked=core.read(HERE/'research_round_741_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round741_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第741轮完成：** [实际绝对来源与一阶共同反作用]({p}research_note_741.md)原全部线性场方程、初始约束及变背景记录态族与同一绝对来源接通，完整方程／约束残差O(λ²)。固定零阶准备时首阶不依赖辅助一阶扩展，二阶仍保历史差。两组、十八式通过，最新741／3421，1535份编号科学文件、3373份保护证据。[核验]({p}research_round_741_checks.json)、[条件账]({p}unified_physics_condition_ledger_741.md)。精确解、实际误差与内部准备仍开放。'
order='**当前执行顺序（741后，优先于下方历史安排）：** 接[742现有记录的局域过程与共同来源]({p}round742_drafts/STATUS.md)，先回查524、575、593、620、633—634及716，核同一操作和总来源。旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第741轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 387.' not in text
        text+='\n\n## 387. 实际绝对来源、联合初值与一阶共同反作用\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 292.' not in text
        text+='\n\n## 292. 旧空间合同保持，一阶绝对发展不替代内部准备\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—740轮。' in text and '最新科学轮次与检查数为740／3419' in text
        text=text.replace('完成231—740轮。','完成231—741轮。').replace('最新科学轮次与检查数为740／3419','最新科学轮次与检查数为741／3421')
    if p==HERE/'README.md':
        text+='\n\n## 第741轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|741|[实际绝对来源、联合初值与一阶共同反作用](research_note_741.md)|[代码](joint_absolute_source_development.py)、[结果](joint_absolute_source_development_results.json)、[核验](research_round_741_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round741.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=741,cumulative_tests=3421,numbered_scientific_files=1535,
    unique_protected_evidence_files=3373,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=742,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
