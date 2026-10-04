"""Publish verified 655, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round655_navigation_checks.json'
checked=core.read(HERE/'research_round_655_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round655_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第655轮完成：** [原完整图的质量转移、共同时间来源与几何归一边界]({p}research_note_655.md)'
         '正质量分割接回原全图、Gauss与固定原态记录，正时间常lapse来源共同收敛。'
         '辅助欧氏投影不等于真实遗忘；固定计数归一不可直接认作协变真空体积项。'
         '四组、二十二式通过，最新655／3206，1277份编号科学文件、2218份保护证据。'
         '[核验]({p}research_round_655_checks.json)、[条件账]({p}unified_physics_condition_ledger_655.md)、'
         '[654式15勘误]({p}round655_drafts/round654_formula15_erratum.md)。'
         '限原有限图；实际空间手征测度、连续与量子引力仍开放。')
order=('**当前执行顺序（655后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[656实际空间overlap测度与正图过程]({p}round656_drafts/STATUS.md)，'
       '恢复原辅助Pfaffian的真实空间依赖，核共同映射；不继续单独优化时间步，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第655轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 301.' not in text
        text+='\n\n## 301. 原全图过程与几何归一\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 206.' not in text
        text+='\n\n## 206. 全图演化、时间来源与体积归一须共同签收\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—654轮。' in text and '最新科学轮次与检查数为654／3202' in text
        text=text.replace('完成231—654轮。','完成231—655轮。').replace('最新科学轮次与检查数为654／3202','最新科学轮次与检查数为655／3206')
    if p==HERE/'README.md':
        text+='\n\n## 第655轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|655|[原完整图的质量转移、共同时间来源与几何归一边界](research_note_655.md)|[代码](joint_full_graph_transfer_sources.py)、[结果](joint_full_graph_transfer_sources_results.json)、[核验](research_round_655_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round655.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=655,cumulative_tests=3206,numbered_scientific_files=1277,
            unique_protected_evidence_files=2218,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=656,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
