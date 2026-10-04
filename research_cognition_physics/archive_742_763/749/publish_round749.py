"""Publish749 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round749_navigation_checks.json'
checked=core.read(HERE/'research_round_749_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round749_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第749轮完成（限定反例）：** [合并映射与非线性几何约束]({p}research_note_749.md)原731两份合法外给来源解，遗忘标签后直接平均psi及canonical资料，会有严格正二阶Hamiltonian约束缺陷；只排除该合并映射，不把外给源当实际记录或否定全部粗化。两组、十二式通过，最新749／3436，1559份编号科学文件、3493份保护证据。[核验]({p}research_round_749_checks.json)、[条件账]({p}unified_physics_condition_ledger_749.md)。'
order='**当前执行顺序（749后，优先于下方历史安排）：** 接[750实际标记来源与共同几何响应]({p}round750_drafts/STATUS.md)，复用634、731及741，核同一记录、来源、补偿与几何的对象身份；不重新证明已知记录协方差。统一目标与旧空间、604、649／699保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第749轮完成（限定反例）：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 395.' not in text
        text+='\n\n## 395. 合并映射与非线性几何约束\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 300.' not in text
        text+='\n\n## 300. 旧空间合同保持，合并资料与非线性约束共同检验\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—748轮。' in text and '最新科学轮次与检查数为748／3434' in text
        text=text.replace('完成231—748轮。','完成231—749轮。').replace('最新科学轮次与检查数为748／3434','最新科学轮次与检查数为749／3436')
    if p==HERE/'README.md':
        text+='\n\n## 第749轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|749|[合并映射与非线性几何约束](research_note_749.md)|[代码](joint_record_geometry_averaging.py)、[结果](joint_record_geometry_averaging_results.json)、[核验](research_round_749_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round749.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=749,cumulative_tests=3436,numbered_scientific_files=1559,
    unique_protected_evidence_files=3493,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=750,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
