"""Publish verified 704, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round704_navigation_checks.json'
checked=core.read(HERE/'research_round_704_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round704_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第704轮完成：** [准备几何与真实记录历史的共同响应]({p}research_note_704.md)'
         '同一625有限CP族的自身准备Gibbs态、真实演化与记录来源共同收敛到总阶数二，'
         '包括准备／动力混合项；夹权热导数接入原弱响应。'
         '两组、十八式通过，最新704／3315，1424份编号科学文件、2855份保护证据。'
         '[核验]({p}research_round_704_checks.json)、[全条件账]({p}unified_physics_condition_ledger_704.md)。'
         '固定图有限过程接口接通，空间区域／连续及动态量子反馈仍开放；旧空间与目标不变。')
order=('**当前执行顺序（704后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[705局部区域与共同尺度映射]({p}round705_drafts/STATUS.md)，'
       '回查617／637—644等，区分精确切分、全谱逼近与真实局部截断；'
       '保留四分支与699范围，停止热导数和有限谱精度优化。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第704轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 350.' not in text
        text+='\n\n## 350. 准备、演化与记录来源共用一个有限过程\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 255.' not in text
        text+='\n\n## 255. 旧空间接口继承，固定图混合响应不替代空间映射\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—703轮。' in text and '最新科学轮次与检查数为703／3313' in text
        text=text.replace('完成231—703轮。','完成231—704轮。').replace('最新科学轮次与检查数为703／3313','最新科学轮次与检查数为704／3315')
    if p==HERE/'README.md':
        text+='\n\n## 第704轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|704|[准备几何与真实记录历史的共同响应](research_note_704.md)|[代码](joint_preparation_history_limit.py)、[结果](joint_preparation_history_limit_results.json)、[核验](research_round_704_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round704.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=704,cumulative_tests=3315,numbered_scientific_files=1424,
            unique_protected_evidence_files=2855,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=705,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
