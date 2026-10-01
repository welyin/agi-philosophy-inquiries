"""Publish verified 623, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round623_navigation_checks.json'
checked=core.read(HERE/'research_round_623_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round623_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第623轮完成：** [原完整量子过程的共同算符域与真实热响应]({p}research_note_623.md)'
         '原束缚势共同控制动能、边耦合与CAR质量；指定参数族在同一Gauss热态中具有实际二阶过程展开，'
         '接通622的噪声、迟致及接触项。'
         '三组、十八式通过，最新623／3112，1181份编号科学文件、1958份保护证据。'
         '[核验]({p}research_round_623_checks.json)、[条件账]({p}unified_physics_condition_ledger_623.md)。'
         '限固定图与正外部几何；连续、动态引力与GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（623后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[624原实际记录与共同量子来源的接口]({p}round624_drafts/STATUS.md)，'
       '先回查原记录及尺度合同，再补真实共同过程的缺口；目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第623轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 269.' not in text
        text+='\n\n## 269. 共同算符域、原物质与真实热过程\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 174.' not in text
        text+='\n\n## 174. 原完整过程展开不构成维数选择\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—622轮。' in text and '最新科学轮次与检查数为622／3109' in text
        text=text.replace('完成231—622轮。','完成231—623轮。').replace('最新科学轮次与检查数为622／3109','最新科学轮次与检查数为623／3112')
    if p==HERE/'README.md':
        text+='\n\n## 第623轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|623|[原完整量子过程的共同算符域与真实热响应](research_note_623.md)|[代码](joint_operator_domain_completion.py)、[结果](joint_operator_domain_completion_results.json)、[核验](research_round_623_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round623.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=623,cumulative_tests=3112,numbered_scientific_files=1181,
            unique_protected_evidence_files=1958,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=624,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
