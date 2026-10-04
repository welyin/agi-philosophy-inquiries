"""Publish verified 626, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round626_navigation_checks.json'
checked=core.read(HERE/'research_round_626_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round626_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第626轮完成：** [共同参考态对局部有效作用的记忆限制]({p}research_note_626.md)'
         '原完整Gauss热态的均匀来源保跨时间关联；'
         '消去全部物质后只保源局部二次项，不能匹配两分离脉冲。'
         '同一初始能谱精确补全该来源部门。'
         '三组、十四式通过，最新626／3121，1190份编号科学文件、1980份保护证据。'
         '[核验]({p}research_round_626_checks.json)、[条件账]({p}unified_physics_condition_ledger_626.md)。'
         '数值限原32模式条件物质；不排除保留物质的局部EFT，连续与GR仍开放。')
order=('**当前执行顺序（626后，优先于下方历史安排）：** 继续整合共同条件，认知设计后置。'
       '接[627共同有效作用的保留与消去部门]({p}round627_drafts/STATUS.md)，'
       '核原规范、物质及几何来源的共同分工和态；不继续脉冲或记忆资源优化，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第626轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 272.' not in text
        text+='\n\n## 272. 共同参考态与局部作用的记忆限制\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 177.' not in text
        text+='\n\n## 177. 时间来源记忆不选择空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—625轮。' in text and '最新科学轮次与检查数为625／3118' in text
        text=text.replace('完成231—625轮。','完成231—626轮。').replace('最新科学轮次与检查数为625／3118','最新科学轮次与检查数为626／3121')
    if p==HERE/'README.md':
        text+='\n\n## 第626轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|626|[共同参考态对局部有效作用的记忆限制](research_note_626.md)|[代码](joint_equilibrium_memory_matching.py)、[结果](joint_equilibrium_memory_matching_results.json)、[核验](research_round_626_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round626.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=626,cumulative_tests=3121,numbered_scientific_files=1190,
            unique_protected_evidence_files=1980,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=627,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
