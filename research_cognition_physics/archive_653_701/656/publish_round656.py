"""Publish verified 656, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round656_navigation_checks.json'
checked=core.read(HERE/'research_round_656_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round656_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第656轮完成：** [原空间手征测度与共同传播投影子的连接]({p}research_note_656.md)'
         '原自由P给全S⁹辅助权重模长、平面精确正性及九方向共同Laplacian。'
         '恢复空间同时改写时间耦合；权重与其来源不能另拟。'
         '四组、二十式通过，最新656／3210，1280份编号科学文件、2229份保护证据。'
         '[核验]({p}research_round_656_checks.json)、[条件账]({p}unified_physics_condition_ledger_656.md)。'
         '限自由被积式；全S⁹积分、正物理过程、连续与量子引力仍开放。')
order=('**当前执行顺序（656后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[657原空间积分与正时间拼接]({p}round657_drafts/STATUS.md)，'
       '核实际测度能否共用正过程与原态；不把平面正性等同完整量子重建，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第656轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 302.' not in text
        text+='\n\n## 302. 原空间测度与共同传播\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 207.' not in text
        text+='\n\n## 207. 空间测度、时间耦合与来源须共同签收\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—655轮。' in text and '最新科学轮次与检查数为655／3206' in text
        text=text.replace('完成231—655轮。','完成231—656轮。').replace('最新科学轮次与检查数为655／3206','最新科学轮次与检查数为656／3210')
    if p==HERE/'README.md':
        text+='\n\n## 第656轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|656|[原空间手征测度与共同传播投影子的连接](research_note_656.md)|[代码](joint_spatial_auxiliary_geometry.py)、[结果](joint_spatial_auxiliary_geometry_results.json)、[核验](research_round_656_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round656.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=656,cumulative_tests=3210,numbered_scientific_files=1280,
            unique_protected_evidence_files=2229,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=657,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
