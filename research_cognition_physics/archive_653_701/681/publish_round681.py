"""Publish verified 681, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round681_navigation_checks.json'
checked=core.read(HERE/'research_round_681_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round681_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第681轮完成：** [规范流、物理时间与共同正性接口的边界]({p}research_note_681.md)'
         '完整自由规范物理模式给任意正时空流的负反射范数；反射协变不足以把流后字段当作局部物理观测。'
         '仅空间流保已有正性，但不完成原电背景衰减。两组、十四式通过，最新681／3269，1355份编号科学文件、2525份保护证据。'
         '[核验]({p}research_round_681_checks.json)、[全条件账]({p}unified_physics_condition_ledger_681.md)。'
         '本轮不判定原全Gauss/S9候选；H_F身份、连续及量子GR仍开放，旧空间接口复用。')
order=('**当前执行顺序（681后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[682完整量子体、条件涨落与原边界来源]({p}round682_drafts/STATUS.md)，'
       '核真实体积分的涨落、权重及全部来源，和678约束体分别记账；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第681轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 327.' not in text
        text+='\n\n## 327. 规范流不自动给共同物理时间\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 232.' not in text
        text+='\n\n## 232. 旧空间合同复用，流后观测须核正性\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—680轮。' in text and '最新科学轮次与检查数为680／3267' in text
        text=text.replace('完成231—680轮。','完成231—681轮。').replace('最新科学轮次与检查数为680／3267','最新科学轮次与检查数为681／3269')
    if p==HERE/'README.md':
        text+='\n\n## 第681轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|681|[规范流、物理时间与共同正性接口的边界](research_note_681.md)|[代码](joint_gauge_flow_time_interface.py)、[结果](joint_gauge_flow_time_interface_results.json)、[核验](research_round_681_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round681.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=681,cumulative_tests=3269,numbered_scientific_files=1355,
            unique_protected_evidence_files=2525,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=682,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
