"""Publish728 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round728_navigation_checks.json'
checked=core.read(HERE/'research_round_728_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round728_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第728轮完成：** [完整物质基态的Gauss提升与原来源分支]({p}research_note_728.md)原非平坦源的基态字符障碍关闭；平坦内部链路分支另证全图正隙及完整一代电荷抵消，局部Gauss参考成立。三组、十四式通过，最新728／3389，1496份编号科学文件、3200份保护证据。[核验]({p}research_round_728_checks.json)、[全条件账]({p}unified_physics_condition_ledger_728.md)。两分支不混作同一Einstein来源。'
order='**当前执行顺序（728后，优先于下方历史安排）：** 接[729原非平坦来源与物质参考]({p}round729_drafts/STATUS.md)，恢复569—574实际规范链路、原字段及来源核费米谱；不以平坦证书替代原源。旧方差与响应结果复用，旧空间、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第728轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 374.' not in text
        text+='\n\n## 374. 完整物质基态的Gauss提升与来源分支\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 279.' not in text
        text+='\n\n## 279. 旧空间合同保持，平坦带隙不替代原非平坦来源\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—727轮。' in text and '最新科学轮次与检查数为727／3386' in text
        text=text.replace('完成231—727轮。','完成231—728轮。').replace('最新科学轮次与检查数为727／3386','最新科学轮次与检查数为728／3389')
    if p==HERE/'README.md':
        text+='\n\n## 第728轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|728|[完整物质基态的Gauss提升与原来源分支](research_note_728.md)|[代码](joint_ground_gauss_lift.py)、[结果](joint_ground_gauss_lift_results.json)、[核验](research_round_728_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round728.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=728,cumulative_tests=3389,numbered_scientific_files=1496,
    unique_protected_evidence_files=3200,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=729,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
