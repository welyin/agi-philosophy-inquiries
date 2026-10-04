"""Publish verified 692, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round692_navigation_checks.json'
checked=core.read(HERE/'research_round_692_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round692_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第692轮完成：** [完整物理泛函的辅助轨道与九点精确求积]({p}research_note_692.md)'
         '保留原H_b、全配置及双Gauss平均，对全部E独立Gauss物理来源，每节点S9精确化为九点正权求和。'
         '两组、十四式通过，最新692／3291，1388份编号科学文件、2655份保护证据。'
         '[核验]({p}research_round_692_checks.json)、[全条件账]({p}unified_physics_condition_ledger_692.md)。'
         '不替代固定背景球面积分，不凭求积正权领取RP；原H_F、共同连续及量子GR仍开放。')
order=('**当前执行顺序（692后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[693辅助消元后的原规范来源]({p}round693_drafts/STATUS.md)，'
       '优先核完整平均的符号控制或精确分块，不重复纯时间链与无控制抽样；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第692轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 338.' not in text
        text+='\n\n## 338. 完整物理泛函的辅助轨道与精确求积\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 243.' not in text
        text+='\n\n## 243. 旧空间合同复用，辅助轨道代表不是物理参考\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—691轮。' in text and '最新科学轮次与检查数为691／3289' in text
        text=text.replace('完成231—691轮。','完成231—692轮。').replace('最新科学轮次与检查数为691／3289','最新科学轮次与检查数为692／3291')
    if p==HERE/'README.md':
        text+='\n\n## 第692轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|692|[完整物理泛函的辅助轨道与九点精确求积](research_note_692.md)|[代码](joint_auxiliary_orbit_quadrature.py)、[结果](joint_auxiliary_orbit_quadrature_results.json)、[核验](research_round_692_checks.json)|\n'
    planned[p]=text.replace('\n',newline).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir(exist_ok=False)
manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as stream:stream.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as stream:json.dump(manifest,stream,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round692.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=692,cumulative_tests=3291,numbered_scientific_files=1388,
            unique_protected_evidence_files=2655,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=693,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
