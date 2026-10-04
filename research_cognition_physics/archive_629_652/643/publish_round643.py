"""Publish verified 643, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round643_navigation_checks.json'
checked=core.read(HERE/'research_round_643_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round643_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第643轮完成：** [完整Gauss物质、共同热参考与来源的有序历史表示]({p}research_note_643.md)'
         '原全模型满足矩阵热桥与Gauss闭合的共同可积条件；原32模式验证有序费米权重及来源插入。'
         '三组、十四式通过，最新643／3168，1241份编号科学文件、2104份保护证据。'
         '[核验]({p}research_round_643_checks.json)、[全条件账]({p}unified_physics_condition_ledger_643.md)。'
         '数值为条件因子，未积分全图；连续手征、动态几何及GR仍开放。')
order=('**当前执行顺序（643后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[644共同区域态、物质与几何熵]({p}round644_drafts/STATUS.md)，'
       '核同一热历史的区域拼接、规范边界及来源，不继续条件因子精度扫描；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第643轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 289.' not in text
        text+='\n\n## 289. 完整Gauss物质的共同热历史与来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 194.' not in text
        text+='\n\n## 194. 完整热桥表示不替代连续几何推导\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—642轮。' in text and '最新科学轮次与检查数为642／3165' in text
        text=text.replace('完成231—642轮。','完成231—643轮。').replace('最新科学轮次与检查数为642／3165','最新科学轮次与检查数为643／3168')
    if p==HERE/'README.md':
        text+='\n\n## 第643轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|643|[完整Gauss物质、共同热参考与来源的有序历史表示](research_note_643.md)|[代码](joint_gauss_fermion_influence.py)、[结果](joint_gauss_fermion_influence_results.json)、[核验](research_round_643_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round643.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=643,cumulative_tests=3168,numbered_scientific_files=1241,
            unique_protected_evidence_files=2104,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=644,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
