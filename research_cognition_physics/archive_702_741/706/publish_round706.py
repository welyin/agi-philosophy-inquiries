"""Publish verified 706, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round706_navigation_checks.json'
checked=core.read(HERE/'research_round_706_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round706_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第706轮完成：** [保边界局部近似与真实记录来源]({p}research_note_706.md)'
         '同一保边界局部CP族接通原真实历史、准备／动力总阶数二来源及零阶cq态；'
         '仍使用原完整H及原Gibbs，整体无限维。'
         '两组、十八式通过，最新706／3320，1430份编号科学文件、2885份保护证据。'
         '[核验]({p}research_round_706_checks.json)、[全条件账]({p}unified_physics_condition_ledger_706.md)。'
         '705宇称补充已核，空间共同尺度仍开放，旧空间与目标不变。')
order=('**当前执行顺序（706后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[707实际共同空间尺度]({p}round707_drafts/STATUS.md)，'
       '回查637—646、663—667，区分区域合并、图细化与固定图谱精度；'
       '保留四分支及699范围，不重复固定图热迹、局部矩或容量反例。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第706轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 352.' not in text
        text+='\n\n## 352. 局部区域、真实记录与来源共同验收\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 257.' not in text
        text+='\n\n## 257. 旧空间接口继承，局部来源极限不当空间细化\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—705轮。' in text and '最新科学轮次与检查数为705／3318' in text
        text=text.replace('完成231—705轮。','完成231—706轮。').replace('最新科学轮次与检查数为705／3318','最新科学轮次与检查数为706／3320')
    if p==HERE/'README.md':
        text+='\n\n## 第706轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|706|[保边界局部近似与真实记录来源](research_note_706.md)|[代码](joint_local_history_source_limit.py)、[结果](joint_local_history_source_limit_results.json)、[核验](research_round_706_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round706.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=706,cumulative_tests=3320,numbered_scientific_files=1430,
            unique_protected_evidence_files=2885,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=707,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
