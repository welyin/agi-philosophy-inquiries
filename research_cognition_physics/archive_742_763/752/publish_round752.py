"""Publish752 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round752_navigation_checks.json'
checked=core.read(HERE/'research_round_752_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round752_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第752轮完成：** [关系定位、实际交互与共同光锥]({p}research_note_752.md)原导数参考窗口若进入指定局部作用，两径向精确共锥只允许A/B仿射；非零四标签紧支撑接法被排除，正目标度规的仿射类保留。数值限原初始jet，未求新全约束解。两组、十六式通过，最新752／3443，1568份编号科学文件、3530份保护证据。[核验]({p}research_round_752_checks.json)、[条件账]({p}unified_physics_condition_ledger_752.md)。'
order='**当前执行顺序（752后，优先于下方历史安排）：** 接[753原过程与关系报告]({p}round753_drafts/STATUS.md)，回到原动力学检验描述性窗口、实际交互和共同来源；旧双历史及关系读口直接复用，不继续窗口优化。[范围审计]({p}round752_drafts/scope_and_dedup_review.md)。旧空间、604、649／699与目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第752轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 398.' not in text
        text+='\n\n## 398. 关系定位、实际交互与共同光锥\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 303.' not in text
        text+='\n\n## 303. 旧空间合同保持，原参考交互与传播相容审计\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—751轮。' in text and '最新科学轮次与检查数为751／3441' in text
        text=text.replace('完成231—751轮。','完成231—752轮。').replace('最新科学轮次与检查数为751／3441','最新科学轮次与检查数为752／3443')
    if p==HERE/'README.md':
        text+='\n\n## 第752轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|752|[关系定位、实际交互与共同光锥](research_note_752.md)|[代码](joint_relational_interaction_cones.py)、[结果](joint_relational_interaction_cones_results.json)、[核验](research_round_752_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round752.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=752,cumulative_tests=3443,numbered_scientific_files=1568,
    unique_protected_evidence_files=3530,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=753,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
