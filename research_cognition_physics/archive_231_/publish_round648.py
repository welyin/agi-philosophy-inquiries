"""Publish verified 648, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round648_navigation_checks.json'
checked=core.read(HERE/'research_round_648_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round648_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第648轮完成：** [原物质参考、引力角点与区域量子拼接的共同条件]({p}research_note_648.md)'
         '原F同时固定角点面积荷与标量来源，原物质参考给同源局部角点；'
         '指定正则boost拼接无非零正常不变态，有限规范区间近似没有强极限。'
         '三组、十四式通过，最新648／3182，1256份编号科学文件、2144份保护证据。'
         '[核验]({p}research_round_648_checks.json)、[条件账]({p}unified_physics_condition_ledger_648.md)。'
         '结论限经典角点及声明的量子接法，未推出完整引力量子态或统计面积律。')
order=('**当前执行顺序（648后，优先于下方历史安排）：** 先整合共同条件，认知系统设计后置。'
       '接[649物理区域内积与同一过程的下降]({p}round649_drafts/STATUS.md)，'
       '核非紧约化的替代方案与原过程是否相容；不扫描规范区间精度，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第648轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 294.' not in text
        text+='\n\n## 294. 物质参考、角点荷与区域匹配\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 199.' not in text
        text+='\n\n## 199. 紧群区域拼接不能直接替代引力约化\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—647轮。' in text and '最新科学轮次与检查数为647／3179' in text
        text=text.replace('完成231—647轮。','完成231—648轮。').replace('最新科学轮次与检查数为647／3179','最新科学轮次与检查数为648／3182')
    if p==HERE/'README.md':
        text+='\n\n## 第648轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|648|[原物质参考、引力角点与区域量子拼接的共同条件](research_note_648.md)|[代码](joint_relational_corner_gluing.py)、[结果](joint_relational_corner_gluing_results.json)、[核验](research_round_648_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round648.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=648,cumulative_tests=3182,numbered_scientific_files=1256,
            unique_protected_evidence_files=2144,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=649,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
