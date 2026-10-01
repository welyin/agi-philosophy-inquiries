"""Publish verified 653, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round653_navigation_checks.json'
checked=core.read(HERE/'research_round_653_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round653_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第653轮完成：** [原手征辅助测度、规范状态与时间拼接的共同表示]({p}research_note_653.md)'
         '原32CAR直接迹接法有电荷谱障碍；原整个商群的辅助积分有正表示，'
         '实际空间平凡时间链给同一正转移、规范投影与来源。'
         '四组、二十式通过，最新653／3198，1271份编号科学文件、2196份保护证据。'
         '[核验]({p}research_round_653_checks.json)、[条件账]({p}unified_physics_condition_ledger_653.md)。'
         '限零标量时间分支；原全图动力学、连续手征和引力仍开放。')
order=('**当前执行顺序（653后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[654原空间传播、质量与正转移]({p}round654_drafts/STATUS.md)，'
       '检验原全图对象与本次手征状态／演化的实际连接；不优化辅助维数，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第653轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 299.' not in text
        text+='\n\n## 299. 原辅助测度与规范时间拼接\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 204.' not in text
        text+='\n\n## 204. 手征状态须与测度及时间拼接共同相容\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—652轮。' in text and '最新科学轮次与检查数为652／3194' in text
        text=text.replace('完成231—652轮。','完成231—653轮。').replace('最新科学轮次与检查数为652／3194','最新科学轮次与检查数为653／3198')
    if p==HERE/'README.md':
        text+='\n\n## 第653轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|653|[原手征辅助测度、规范状态与时间拼接的共同表示](research_note_653.md)|[代码](joint_chiral_character_state.py)、[结果](joint_chiral_character_state_results.json)、[核验](research_round_653_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round653.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=653,cumulative_tests=3198,numbered_scientific_files=1271,
            unique_protected_evidence_files=2196,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=654,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
