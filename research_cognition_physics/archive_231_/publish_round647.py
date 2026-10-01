"""Publish verified 647, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round647_navigation_checks.json'
checked=core.read(HERE/'research_round_647_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round647_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第647轮完成：** [原物质参考、关系读取与几何定位来源]({p}research_note_647.md)'
         '原读口与原物质参考共同定位后，完整来源须保读量、体积和定位三项；'
         '冻结定位有严格纯坐标反例，原初片复算验证三项抵消并连接关系体积。'
         '两组、十六式通过，最新647／3179，1253份编号科学文件、2135份保护证据。'
         '[核验]({p}research_round_647_checks.json)、[全条件账]({p}unified_physics_condition_ledger_647.md)。'
         '限经典关系接口；量子instrument、连续映射及动态量子几何仍开放。')
order=('**当前执行顺序（647后，优先于下方历史安排）：** 先尽量整合已有条件，认知系统设计后置。'
       '接[648关系区域、共同量子来源与实际操作]({p}round648_drafts/STATUS.md)，'
       '优先检查同一过程、同一区域和同一来源的连接；不另开装置或参考优化，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第647轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 293.' not in text
        text+='\n\n## 293. 原物质参考与关系定位的共同来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 198.' not in text
        text+='\n\n## 198. 读取、体积与定位必须共同输送\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—646轮。' in text and '最新科学轮次与检查数为646／3177' in text
        text=text.replace('完成231—646轮。','完成231—647轮。').replace('最新科学轮次与检查数为646／3177','最新科学轮次与检查数为647／3179')
    if p==HERE/'README.md':
        text+='\n\n## 第647轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|647|[原物质参考、关系读取与几何定位来源](research_note_647.md)|[代码](joint_relational_readout_source.py)、[结果](joint_relational_readout_source_results.json)、[核验](research_round_647_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round647.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=647,cumulative_tests=3179,numbered_scientific_files=1253,
            unique_protected_evidence_files=2135,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=648,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
