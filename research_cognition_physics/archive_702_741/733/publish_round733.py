"""Publish733 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round733_navigation_checks.json'
checked=core.read(HERE/'research_round_733_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round733_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第733轮完成：** [有限记录模式的非线性边界与完整参考响应]({p}research_note_733.md)有限记录可合法嵌入新参考，但明确非共形变化的完整参考差非紧，不能用有限光滑模式吸收；跨背景反馈须计参考及减除响应。三组、十八式通过，最新733／3404，1511份编号科学文件、3270份保护证据。[核验]({p}research_round_733_checks.json)、[全条件账]({p}unified_physics_condition_ledger_733.md)。732一阶结果保持，绝对反馈开放。'
order='**当前执行顺序（733后，优先于下方历史安排）：** 接[734同一过去准备的因果参考响应]({p}round734_drafts/STATUS.md)，核完整态、源算符、共同减除及初态项；区别保记录的新准备与真实退迟响应。停止有限补丁、秩及精度优化；旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第733轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 379.' not in text
        text+='\n\n## 379. 有限记录模式的非线性边界与完整参考响应\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 284.' not in text
        text+='\n\n## 284. 旧空间合同保持，完整参考响应不由有限记录替代\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—732轮。' in text and '最新科学轮次与检查数为732／3401' in text
        text=text.replace('完成231—732轮。','完成231—733轮。').replace('最新科学轮次与检查数为732／3401','最新科学轮次与检查数为733／3404')
    if p==HERE/'README.md':
        text+='\n\n## 第733轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|733|[有限记录模式的非线性边界与完整参考响应](research_note_733.md)|[代码](joint_reference_polarization_boundary.py)、[结果](joint_reference_polarization_boundary_results.json)、[核验](research_round_733_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round733.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=733,cumulative_tests=3404,numbered_scientific_files=1511,
    unique_protected_evidence_files=3270,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=734,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
