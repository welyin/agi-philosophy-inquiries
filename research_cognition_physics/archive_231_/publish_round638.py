"""Publish verified 638, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round638_navigation_checks.json'
checked=core.read(HERE/'research_round_638_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round638_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第638轮完成：** [原记录、共同参考与区域相对熵的联合近似]({p}research_note_638.md)'
         '原完整Gibbs与记录共用信息—注能预算；'
         '有限参考配原读取产生无限相对熵，共同保Gauss粗化恢复正确区域极限。'
         '三组、十六式通过，最新638／3154，1226份编号科学文件、2067份保护证据。'
         '[核验]({p}research_round_638_checks.json)、[条件账]({p}unified_physics_condition_ledger_638.md)。'
         '限固定图；统一空间尺度、连续几何熵及GR仍开放。')
order=('**当前执行顺序（638后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[639固定物理支撑、共同几何与信息预算]({p}round639_drafts/STATUS.md)，'
       '核原区域读取能否保尺度统一预算及几何来源；因果实现单列，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第638轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 284.' not in text
        text+='\n\n## 284. 原记录、参考和区域相对熵的共同近似\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 189.' not in text
        text+='\n\n## 189. 谱相对熵极限不替代空间连续匹配\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—637轮。' in text and '最新科学轮次与检查数为637／3151' in text
        text=text.replace('完成231—637轮。','完成231—638轮。').replace('最新科学轮次与检查数为637／3151','最新科学轮次与检查数为638／3154')
    if p==HERE/'README.md':
        text+='\n\n## 第638轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|638|[原记录、共同参考与区域相对熵的联合近似](research_note_638.md)|[代码](joint_record_relative_entropy.py)、[结果](joint_record_relative_entropy_results.json)、[核验](research_round_638_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round638.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=638,cumulative_tests=3154,numbered_scientific_files=1226,
            unique_protected_evidence_files=2067,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=639,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
