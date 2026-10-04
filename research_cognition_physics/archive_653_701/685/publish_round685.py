"""Publish verified 685, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round685_navigation_checks.json'
checked=core.read(HERE/'research_round_685_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round685_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第685轮完成：** [原体权重、非零物理支撑与匹配边界]({p}research_note_685.md)'
         '复用615全S9积分，证明未减除体因子改变原非零物理支撑，曲率或体积常数匹配不足。'
         '两组、十六式通过，最新685／3277，1367份编号科学文件、2579份保护证据。'
         '[核验]({p}research_round_685_checks.json)、[全条件账]({p}unified_physics_condition_ledger_685.md)。'
         '原Q0正性、H_F身份、共同连续及量子GR仍开放，旧空间接口复用。')
order=('**当前执行顺序（685后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[686谱体匹配与原完整物理来源]({p}round686_drafts/STATUS.md)，'
       '核残余界、原完整平均与来源；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第685轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 331.' not in text
        text+='\n\n## 331. 原体权重接到实际非零物理支撑\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 236.' not in text
        text+='\n\n## 236. 旧空间合同复用，匹配原权重与来源\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—684轮。' in text and '最新科学轮次与检查数为684／3275' in text
        text=text.replace('完成231—684轮。','完成231—685轮。').replace('最新科学轮次与检查数为684／3275','最新科学轮次与检查数为685／3277')
    if p==HERE/'README.md':
        text+='\n\n## 第685轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|685|[原体权重、非零物理支撑与匹配边界](research_note_685.md)|[代码](joint_body_weight_source_matching.py)、[结果](joint_body_weight_source_matching_results.json)、[核验](research_round_685_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round685.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=685,cumulative_tests=3277,numbered_scientific_files=1367,
            unique_protected_evidence_files=2579,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=686,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
