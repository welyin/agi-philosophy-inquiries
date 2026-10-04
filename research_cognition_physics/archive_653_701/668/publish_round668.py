"""Publish verified 668, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round668_navigation_checks.json'
checked=core.read(HERE/'research_round_668_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round668_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第668轮完成：** [原非零质量、物理观测与辅助测度的共同反射结构]({p}research_note_668.md)'
         '原常数质量接入同一物理Weyl观测，有限自由盒有全多项式反射正性和严格归一；来源保留观测映射导数。'
         '两组、十六式通过，最新668／3243，1316份编号科学文件、2354份保护证据。'
         '[核验]({p}research_round_668_checks.json)、[全条件账]({p}unified_physics_condition_ledger_668.md)。'
         '原动态标量、完整Gauss过程、连续及量子GR仍开放；旧空间接口直接继承。')
order=('**当前执行顺序（668后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[669原动态标量、物理质量与记录]({p}round669_drafts/STATUS.md)，'
       '核原标量热过程、反射和积分，再接原s读口；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第668轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 314.' not in text
        text+='\n\n## 314. 原质量与物理观测的共同正泛函\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 219.' not in text
        text+='\n\n## 219. 固定原质量的反射接口与动态标量缺口\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—667轮。' in text and '最新科学轮次与检查数为667／3241' in text
        text=text.replace('完成231—667轮。','完成231—668轮。').replace('最新科学轮次与检查数为667／3241','最新科学轮次与检查数为668／3243')
    if p==HERE/'README.md':
        text+='\n\n## 第668轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|668|[原非零质量、物理观测与辅助测度的共同反射结构](research_note_668.md)|[代码](joint_mass_auxiliary_reflection.py)、[结果](joint_mass_auxiliary_reflection_results.json)、[核验](research_round_668_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round668.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=668,cumulative_tests=3243,numbered_scientific_files=1316,
            unique_protected_evidence_files=2354,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=669,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
