"""Publish verified 652, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round652_navigation_checks.json'
checked=core.read(HERE/'research_round_652_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round652_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第652轮完成：** [原物质量子参考、共同能量与区域表示的联合连接]({p}research_note_652.md)'
         '原光滑物质参考及速度平方共享自伴表示、原总能量矩控制和规范区域映射；'
         '同一谱子空间以正形式保留平方漏出项及几何偏导。'
         '三组、十六式通过，最新652／3194，1268份编号科学文件、2186份保护证据。'
         '[核验]({p}research_round_652_checks.json)、[条件账]({p}unified_physics_condition_ledger_652.md)。'
         '限原固定图共同表示；联合量子坐标、连续手征和动态引力仍开放。')
order=('**当前执行顺序（652后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[653完整图物质与连续手征分支]({p}round653_drafts/STATUS.md)，'
       '核同一状态、规范测度和动力学的真实连接，不继续参考装置或精度优化；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第652轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 298.' not in text
        text+='\n\n## 298. 原量子参考与同一能量表示\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 203.' not in text
        text+='\n\n## 203. 参考量子表示须与原过程共同相容\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—651轮。' in text and '最新科学轮次与检查数为651／3191' in text
        text=text.replace('完成231—651轮。','完成231—652轮。').replace('最新科学轮次与检查数为651／3191','最新科学轮次与检查数为652／3194')
    if p==HERE/'README.md':
        text+='\n\n## 第652轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|652|[原物质量子参考、共同能量与区域表示的联合连接](research_note_652.md)|[代码](joint_quantum_reference_forms.py)、[结果](joint_quantum_reference_forms_results.json)、[核验](research_round_652_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round652.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=652,cumulative_tests=3194,numbered_scientific_files=1268,
            unique_protected_evidence_files=2186,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=653,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
