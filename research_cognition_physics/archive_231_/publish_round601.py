"""Publish verified 601, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round601_navigation_checks.json'
checked=core.read(HERE/'research_round_601_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round601_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第601轮完成：** [同一有效阶的物质演化与引力约束]({p}research_note_601.md)'
         '将原已知标量有效项的lapse/shift来源、物质演化和截断阶数共同核验：'
         '只改引力来源留下同阶约束漂移，同步修正物质后降至下一阶；一阶Hamiltonian代数余项也明确为二阶。'
         '三组、十四式核验通过，最新601／3048，1115份编号科学文件、1796份保护证据。'
         '[核验]({p}research_round_601_checks.json)、[条件账]({p}unified_physics_condition_ledger_601.md)。'
         '具体复算限原标量部门及给定Einstein领先作用，未完成全量子约束；未取得独立代理审查。')
order=('**当前执行顺序（601后，优先于下方历史安排）：** 继续先整合条件，后置认知实现设计。'
       '接[602共同量子态、来源与尺度]({p}round602_drafts/STATUS.md)，'
       '先去重既有图与连续模型接口，再选择真实匹配；不把局部有效约束相容当成统一模型完成，目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第601轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 247.' not in text
        text+='\n\n## 247. 共同有效阶的约束与物质演化\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 152.' not in text
        text+='\n\n## 152. 同阶有效约束相容仍依赖给定几何与维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—600轮。' in text and '最新科学轮次与检查数为600／3045' in text
        text=text.replace('完成231—600轮。','完成231—601轮。').replace('最新科学轮次与检查数为600／3045','最新科学轮次与检查数为601／3048')
    if p==HERE/'README.md':
        text+='\n\n## 第601轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|601|[同一有效阶的物质演化与约束](research_note_601.md)|[代码](joint_eft_constraint_order.py)、[结果](joint_eft_constraint_order_results.json)、[核验](research_round_601_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round601.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=601,cumulative_tests=3048,numbered_scientific_files=1115,
            unique_protected_evidence_files=1796,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=602,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
