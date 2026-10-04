"""Publish verified 683, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round683_navigation_checks.json'
checked=core.read(HERE/'research_round_683_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round683_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第683轮完成：** [同时间片物理读口、完整来源极限与辅助反射边界]({p}research_note_683.md)'
         '原辅助观测可改为严格同时间片读口，规范协变且以层数一致O(a)误差恢复全部原未归一来源；'
         '指定补偿反射有精确自由障碍。两组、十六式通过，最新683／3273，1361份编号科学文件、2555份保护证据。'
         '[核验]({p}research_round_683_checks.json)、[全条件账]({p}unified_physics_condition_ledger_683.md)。'
         '原Q0正性、H_F身份、共同连续及量子GR仍开放，旧空间接口复用。')
order=('**当前执行顺序（683后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[684辅助反射的对合与共同正时间条件]({p}round684_drafts/STATUS.md)，'
       '核实际反射、全部来源及补偿，不把单独辅助障碍当原物理反例；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第683轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 329.' not in text
        text+='\n\n## 329. 同时间片读口接回原来源，辅助反射仍需核验\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 234.' not in text
        text+='\n\n## 234. 旧空间合同复用，不把辅助反射障碍扩为物理反例\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—682轮。' in text and '最新科学轮次与检查数为682／3271' in text
        text=text.replace('完成231—682轮。','完成231—683轮。').replace('最新科学轮次与检查数为682／3271','最新科学轮次与检查数为683／3273')
    if p==HERE/'README.md':
        text+='\n\n## 第683轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|683|[同时间片物理读口、完整来源极限与辅助反射边界](research_note_683.md)|[代码](joint_time_local_source_interface.py)、[结果](joint_time_local_source_interface_results.json)、[核验](research_round_683_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round683.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=683,cumulative_tests=3273,numbered_scientific_files=1361,
            unique_protected_evidence_files=2555,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=684,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
