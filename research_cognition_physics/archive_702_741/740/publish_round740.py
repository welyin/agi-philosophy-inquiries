"""Publish740 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round740_navigation_checks.json'
checked=core.read(HERE/'research_round_740_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round740_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第740轮完成：** [原增长极点与一致有限阶因果响应]({p}research_note_740.md)原全频一圈重求和有增长极点；一致有限阶保完整记忆、因果与线性约束，并有二阶方程残差界。原校准脉冲修正约51%，实际物理窗口尚未认证。两组、十八式通过，最新740／3419，1532份编号科学文件、3359份保护证据。[核验]({p}research_round_740_checks.json)、[条件账]({p}unified_physics_condition_ledger_740.md)。'
order='**当前执行顺序（740后，优先于下方历史安排）：** 接[741实际历史与绝对来源的共同阶次]({p}round741_drafts/STATUS.md)，复用730—735，核实际初值、守恒及有效窗口。方程残差不等于不稳定精确解误差。旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第740轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 386.' not in text
        text+='\n\n## 386. 原总反馈的增长极点与一致有限阶因果响应\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 291.' not in text
        text+='\n\n## 291. 旧空间合同保持，有限阶残差不替代完整解误差\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—739轮。' in text and '最新科学轮次与检查数为739／3417' in text
        text=text.replace('完成231—739轮。','完成231—740轮。').replace('最新科学轮次与检查数为739／3417','最新科学轮次与检查数为740／3419')
    if p==HERE/'README.md':
        text+='\n\n## 第740轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|740|[原总反馈的增长极点与一致有限阶因果响应](research_note_740.md)|[代码](joint_causal_eft_branch.py)、[结果](joint_causal_eft_branch_results.json)、[核验](research_round_740_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round740.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=740,cumulative_tests=3419,numbered_scientific_files=1532,
    unique_protected_evidence_files=3359,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=741,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
