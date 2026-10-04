"""Publish verified 609, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round609_navigation_checks.json'
checked=core.read(HERE/'research_round_609_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round609_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第609轮完成：** [一粒子协变补全、Fock组合与共同质量来源]({p}research_note_609.md)'
         '同一投影的dΓ联络及正项保原允许形式并兼容独立组合；'
         '质量须同步保正常pp／qq及仅pp配对，补模式成对产生有明确泄漏反例。'
         '三组、十四式通过，最新609／3072，1139份编号科学文件、1856份保护证据。'
         '[核验]({p}research_round_609_checks.json)、[条件账]({p}unified_physics_condition_ledger_609.md)。'
         '固定图Gauss／热态结论有条件；原物种嵌入、传播、测度及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（609后，优先于下方历史安排）：** 继续用共同条件连接各部门，认知实现后置。'
       '接[610真实谱投影、链路导数与共同配置域]({p}round610_drafts/STATUS.md)，'
       '先核一致谱隙的局部系数界及热态适用域，不把系数衰减当无界传播证明，目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第609轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 255.' not in text
        text+='\n\n## 255. 一粒子补全、Fock组合与质量共同约束\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 160.' not in text
        text+='\n\n## 160. 一粒子组合相容仍须实际空间传播\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—608轮。' in text and '最新科学轮次与检查数为608／3069' in text
        text=text.replace('完成231—608轮。','完成231—609轮。').replace('最新科学轮次与检查数为608／3069','最新科学轮次与检查数为609／3072')
    if p==HERE/'README.md':
        text+='\n\n## 第609轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|609|[一粒子协变补全、Fock组合与共同质量来源](research_note_609.md)|[代码](joint_fock_covariant_completion.py)、[结果](joint_fock_covariant_completion_results.json)、[核验](research_round_609_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round609.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=609,cumulative_tests=3072,numbered_scientific_files=1139,
            unique_protected_evidence_files=1856,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=610,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
