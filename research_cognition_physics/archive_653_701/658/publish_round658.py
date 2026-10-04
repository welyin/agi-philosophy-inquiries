"""Publish verified 658, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round658_navigation_checks.json'
checked=core.read(HERE/'research_round_658_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round658_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第658轮完成：** [原球面平均的零模支撑与严格归一]({p}research_note_658.md)'
         '原半区零模只在空间零动量，原S⁹平均在强制占据部门有精确正开链系数。'
         '全部有限自由空间／偶数时间辅助积分严格正，657反射内积可实际归一。'
         '三组、十六式通过，最新658／3217，1286份编号科学文件、2246份保护证据。'
         '[核验]({p}research_round_658_checks.json)、[条件账]({p}unified_physics_condition_ledger_658.md)。'
         '共同单步过程、物理CAR、一般规范场与连续引力仍开放。')
order=('**当前执行顺序（658后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[659跨时间共同演化]({p}round659_drafts/STATUS.md)，'
       '核实际不同时间盒及原过程映射；不优化归一下界或装置，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第658轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 304.' not in text
        text+='\n\n## 304. 原球面平均的零模支撑\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 209.' not in text
        text+='\n\n## 209. 原平均、端点支撑与归一反射内积共同成立\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—657轮。' in text and '最新科学轮次与检查数为657／3214' in text
        text=text.replace('完成231—657轮。','完成231—658轮。').replace('最新科学轮次与检查数为657／3214','最新科学轮次与检查数为658／3217')
    if p==HERE/'README.md':
        text+='\n\n## 第658轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|658|[原球面平均的零模支撑与严格归一](research_note_658.md)|[代码](joint_auxiliary_normalization_support.py)、[结果](joint_auxiliary_normalization_support_results.json)、[核验](research_round_658_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round658.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=658,cumulative_tests=3217,numbered_scientific_files=1286,
            unique_protected_evidence_files=2246,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=659,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
