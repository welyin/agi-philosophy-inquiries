"""Publish732 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round732_navigation_checks.json'
checked=core.read(HERE/'research_round_732_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round732_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第732轮完成：** [同一记录的相对量子来源与短时联合发展]({p}research_note_732.md)原实际记录的有限光滑模式给共同应力、电流与质量力；联合Ward接原初始补偿及短时线性发展，约束沿时间保持。三组、十八式通过，最新732／3401，1508份编号科学文件、3256份保护证据。[核验]({p}research_round_732_checks.json)、[全条件账]({p}unified_physics_condition_ledger_732.md)。绝对源与完整非线性反馈仍开放。'
order='**当前执行顺序（732后，优先于下方历史安排）：** 接[733变化背景的实际参考与相对模式]({p}round733_drafts/STATUS.md)，核有限相对资料能否共同嵌入Hadamard态和同一记录，再研究反馈。停止重复Ward、线性波及小矩阵精度优化；旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第732轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 378.' not in text
        text+='\n\n## 378. 同一记录的相对量子来源与短时联合发展\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 283.' not in text
        text+='\n\n## 283. 旧空间合同保持，相对线性发展不替代完整非线性反馈\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—731轮。' in text and '最新科学轮次与检查数为731／3398' in text
        text=text.replace('完成231—731轮。','完成231—732轮。').replace('最新科学轮次与检查数为731／3398','最新科学轮次与检查数为732／3401')
    if p==HERE/'README.md':
        text+='\n\n## 第732轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|732|[同一记录的相对量子来源与短时联合发展](research_note_732.md)|[代码](joint_relative_source_development.py)、[结果](joint_relative_source_development_results.json)、[核验](research_round_732_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round732.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=732,cumulative_tests=3401,numbered_scientific_files=1508,
    unique_protected_evidence_files=3256,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=733,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
