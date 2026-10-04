"""Publish verified 628, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round628_navigation_checks.json'
checked=core.read(HERE/'research_round_628_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round628_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第628轮完成：** [原商群、旋量几何与全局反常的共同条件]({p}research_note_628.md)'
         '原Z₆商群、普通spin与完整一代通过闭五维费米反常相位检查，包括非提升束；'
         '同一物种的辅助背景指数8同时区分无反常与无零模。'
         '两组、十四式通过，最新628／3125，1196份编号科学文件、1994份保护证据。'
         '[核验]({p}research_round_628_checks.json)、[条件账]({p}unified_physics_condition_ledger_628.md)。'
         '成熟定理的条件性应用；全局权重、连续区域测度、实时过程及GR仍开放。')
order=('**当前执行顺序（628后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[629拓扑权重与共同作用]({p}round629_drafts/STATUS.md)，'
       '回查原相位、质量与区域条件，核剩余共同资料；不重复反常或同类束检查，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第628轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 274.' not in text
        text+='\n\n## 274. 原全局群、物质与spin反常的合并\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 179.' not in text
        text+='\n\n## 179. 反常工具的辅助维数不是时空生成\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—627轮。' in text and '最新科学轮次与检查数为627／3123' in text
        text=text.replace('完成231—627轮。','完成231—628轮。').replace('最新科学轮次与检查数为627／3123','最新科学轮次与检查数为628／3125')
    if p==HERE/'README.md':
        text+='\n\n## 第628轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|628|[原商群、旋量几何与全局反常的共同条件](research_note_628.md)|[代码](joint_global_anomaly_bundle.py)、[结果](joint_global_anomaly_bundle_results.json)、[核验](research_round_628_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round628.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=628,cumulative_tests=3125,numbered_scientific_files=1196,
            unique_protected_evidence_files=1994,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=629,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
