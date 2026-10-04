"""Publish730 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round730_navigation_checks.json'
checked=core.read(HERE/'research_round_730_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round730_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第730轮完成：** [原动态背景上的连续参考、真实记录与同一来源]({p}research_note_730.md)保原非零几何／物质时间jet，完整Dirac／Majorana可接Hadamard参考族及同一光滑记录来源差；无需未来统一瞬时谱隙。三组、十八式通过，最新730／3395，1502份编号科学文件、3228份保护证据。[核验]({p}research_round_730_checks.json)、[全条件账]({p}unified_physics_condition_ledger_730.md)。原图映射、全Gauss及自洽引力仍开放。'
order='**当前执行顺序（730后，优先于下方历史安排）：** 接[731同一量子来源与原初始约束]({p}round731_drafts/STATUS.md)，回查325、572—574、578—583、630—635及649—651，核全部来源、核空间和反作用相容。停止谱隙及辅助过去参数优化；旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第730轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 376.' not in text
        text+='\n\n## 376. 原动态背景上的连续参考、真实记录与同一来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 281.' not in text
        text+='\n\n## 281. 旧空间合同保持，背景连续参考不替代自洽引力\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—729轮。' in text and '最新科学轮次与检查数为729／3392' in text
        text=text.replace('完成231—729轮。','完成231—730轮。').replace('最新科学轮次与检查数为729／3392','最新科学轮次与检查数为730／3395')
    if p==HERE/'README.md':
        text+='\n\n## 第730轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|730|[原动态背景上的连续参考、真实记录与同一来源](research_note_730.md)|[代码](joint_dynamic_continuum_reference.py)、[结果](joint_dynamic_continuum_reference_results.json)、[核验](research_round_730_checks.json)|\n'
    planned[p]=text.replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir(exist_ok=False);manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round730.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=730,cumulative_tests=3395,numbered_scientific_files=1502,
    unique_protected_evidence_files=3228,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=731,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
