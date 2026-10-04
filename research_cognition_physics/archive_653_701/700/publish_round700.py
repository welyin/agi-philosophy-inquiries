"""Publish verified 700, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round700_navigation_checks.json'
checked=core.read(HERE/'research_round_700_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round700_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第700轮完成：** [共同时间、物种与空间局域性]({p}research_note_700.md)'
         '原Wilson的声明各向异性路径中，有限速度固定时间／空间系数比；'
         '时间单独细化会重现轻角点，或失去实际投影的统一空间局域性界。'
         '完整AP时间及原16通道已核；共同细化和动态RP仍开放。'
         '两组、十六式通过，最新700／3307，1412份编号科学文件、2795份保护证据。'
         '[核验]({p}research_round_700_checks.json)、[全条件账]({p}unified_physics_condition_ledger_700.md)。'
         '旧空间接口及统一目标不变。')
order=('**当前执行顺序（700后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[701共同尺度下的物理观测与正性]({p}round701_drafts/STATUS.md)，'
       '核固定盒来源极限到实际时空共同极限的统一量词；'
       '不重复699积分或700角点，四分支继续分别记账。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第700轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 346.' not in text
        text+='\n\n## 346. 共同时间缩放约束物种与实际空间投影\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 251.' not in text
        text+='\n\n## 251. 旧空间合同复用，跨尺度局域性与维数定理分开\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—699轮。' in text and '最新科学轮次与检查数为699／3305' in text
        text=text.replace('完成231—699轮。','完成231—700轮。').replace('最新科学轮次与检查数为699／3305','最新科学轮次与检查数为700／3307')
    if p==HERE/'README.md':
        text+='\n\n## 第700轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|700|[共同时间、物种与空间局域性](research_note_700.md)|[代码](joint_anisotropic_time_contract.py)、[结果](joint_anisotropic_time_contract_results.json)、[核验](research_round_700_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round700.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=700,cumulative_tests=3307,numbered_scientific_files=1412,
            unique_protected_evidence_files=2795,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=701,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
