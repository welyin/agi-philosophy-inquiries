"""Publish verified 671, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round671_navigation_checks.json'
checked=core.read(HERE/'research_round_671_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round671_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第671轮完成：** [静态规范曲率、原手征物质与共同反射泛函]({p}research_note_671.md)'
         '非交换矩阵谱将自由反射输入接到原静态空间规范曲率；物理Weyl、辅助及质量共用未归一正泛函。'
         '严格归一及动态规范过程仍单独验收；时间变化新增交换项已定位。'
         '两组、十八式通过，最新671／3249，1325份编号科学文件、2389份保护证据。'
         '[核验]({p}research_round_671_checks.json)、[全条件账]({p}unified_physics_condition_ledger_671.md)。'
         '原完整过程、动态规范正性、连续及量子GR仍开放。')
order=('**当前执行顺序（671后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[672动态规范链路与原共同反射过程]({p}round672_drafts/STATUS.md)，'
       '核时变链路的完整反射拼接及原过程映射；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第671轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 317.' not in text
        text+='\n\n## 317. 静态曲率与原物理辅助质量的共同反射条件\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 222.' not in text
        text+='\n\n## 222. 静态反射泛函与真实动态规范过程的边界\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—670轮。' in text and '最新科学轮次与检查数为670／3247' in text
        text=text.replace('完成231—670轮。','完成231—671轮。').replace('最新科学轮次与检查数为670／3247','最新科学轮次与检查数为671／3249')
    if p==HERE/'README.md':
        text+='\n\n## 第671轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|671|[静态规范曲率、原手征物质与共同反射泛函](research_note_671.md)|[代码](joint_static_gauge_reflection.py)、[结果](joint_static_gauge_reflection_results.json)、[核验](research_round_671_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round671.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=671,cumulative_tests=3249,numbered_scientific_files=1325,
            unique_protected_evidence_files=2389,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=672,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
