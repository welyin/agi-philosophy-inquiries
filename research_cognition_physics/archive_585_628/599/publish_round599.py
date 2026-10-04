"""Publish verified 599, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round599_navigation_checks.json'
checked=core.read(HERE/'research_round_599_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round599_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第599轮完成：** [同一量子处方的标量、费米与高阶几何来源]({p}research_note_599.md)'
         '在明示固定Einstein背景、canonical费米动能的一圈部门中，'
         '新增费米质量环不能抵消原标量四梯度及RS项，其共同几何来源须保留。'
         '两组、十式核验通过，最新599／3042，1109份编号科学文件、1782份保护证据。'
         '[核验]({p}research_round_599_checks.json)、[条件账]({p}unified_physics_condition_ledger_599.md)。'
         '连续处方仍输入；未算完整规范／引力环、手征相位或固定ℏ图匹配，未取得独立代理审查。')
order=('**当前执行顺序（599后，优先于下方历史安排）：** 继续先整合物理条件，后置认知实现设计。'
       '接[600共同有效作用与引力约束]({p}round600_drafts/STATUS.md)，'
       '联查协变来源、约束和原内部几何的接口；保持各量子部门与尺度边界，统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第599轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 245.' not in text
        text+='\n\n## 245. 相同量子处方中的高阶物质与几何来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 150.' not in text
        text+='\n\n## 150. 共同物质的量子修正仍不选择时空维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—598轮。' in text and '最新科学轮次与检查数为598／3040' in text
        text=text.replace('完成231—598轮。','完成231—599轮。').replace('最新科学轮次与检查数为598／3040','最新科学轮次与检查数为599／3042')
    if p==HERE/'README.md':
        text+='\n\n## 第599轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|599|[标量与费米环的共同高阶来源](research_note_599.md)|[代码](joint_fermion_scalar_loop_matching.py)、[结果](joint_fermion_scalar_loop_matching_results.json)、[核验](research_round_599_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round599.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=599,cumulative_tests=3042,numbered_scientific_files=1109,
            unique_protected_evidence_files=1782,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=600,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
