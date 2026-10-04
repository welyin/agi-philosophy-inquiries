"""Publish verified 674, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round674_navigation_checks.json'
checked=core.read(HERE/'research_round_674_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round674_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第674轮完成：** [互补手征配对、完整规范边界与标量泛函的实性]({p}research_note_674.md)'
         '互补Pfaffian证明原标量／实质量来源的反射身份，含变秩和奇异质量；'
         '接673可积性，完整双Gauss标量核Hermitian、总权重为实。'
         '正性、非零归一、任意费米来源反射与原过程身份仍待证明。'
         '两组、十六式通过，最新674／3255，1334份编号科学文件、2433份保护证据。'
         '[核验]({p}research_round_674_checks.json)、[全条件账]({p}unified_physics_condition_ledger_674.md)。'
         '连续及量子GR仍开放。')
order=('**当前执行顺序（674后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[675完整标量边界核的正性与共同时间]({p}round675_drafts/STATUS.md)，'
       '核全群平均、原H_b与实际时间参数；旧空间接口直接复用，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第674轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 320.' not in text
        text+='\n\n## 320. 原配对、完整边界与标量反射的共同合同\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 225.' not in text
        text+='\n\n## 225. 完整标量泛函的实性与正性的区别\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—673轮。' in text and '最新科学轮次与检查数为673／3253' in text
        text=text.replace('完成231—673轮。','完成231—674轮。').replace('最新科学轮次与检查数为673／3253','最新科学轮次与检查数为674／3255')
    if p==HERE/'README.md':
        text+='\n\n## 第674轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|674|[互补手征配对、完整规范边界与标量泛函的实性](research_note_674.md)|[代码](joint_gauss_reflection_reality.py)、[结果](joint_gauss_reflection_reality_results.json)、[核验](research_round_674_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round674.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=674,cumulative_tests=3255,numbered_scientific_files=1334,
            unique_protected_evidence_files=2433,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=675,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
