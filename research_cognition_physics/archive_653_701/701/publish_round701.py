"""Publish verified 701, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round701_navigation_checks.json'
checked=core.read(HERE/'research_round_701_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round701_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第701轮完成：** [观测分辨率与共同正性极限]({p}research_note_701.md)'
         '681自由磁场的全空间模式及有限来源多项式有受控正极限；'
         '同族缩放放大来源仍保留严格负范数，故物理观测映射与尺度量词必须共同验收。'
         '不修复699原手征候选。两组、十六式通过，最新701／3309，1415份编号科学文件、2809份保护证据。'
         '[核验]({p}research_round_701_checks.json)、[全条件账]({p}unified_physics_condition_ledger_701.md)。'
         '旧空间与524局域探针结果直接复用，目标不变。')
order=('**当前执行顺序（701后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[702原候选反射误差与来源映射]({p}round702_drafts/STATUS.md)，'
       '回到完整H_b／Gauss／S9和原H_F的受控连接；'
       '停止自由Gaussian系数优化，保留四分支及699反例的准确范围。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第701轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 347.' not in text
        text+='\n\n## 347. 物理来源分辨率与正性需要共同取极限\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 252.' not in text
        text+='\n\n## 252. 旧空间及局域探针继承，来源极限不替代维数前提\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—700轮。' in text and '最新科学轮次与检查数为700／3307' in text
        text=text.replace('完成231—700轮。','完成231—701轮。').replace('最新科学轮次与检查数为700／3307','最新科学轮次与检查数为701／3309')
    if p==HERE/'README.md':
        text+='\n\n## 第701轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|701|[观测分辨率与共同正性极限](research_note_701.md)|[代码](joint_flow_resolution_limit.py)、[结果](joint_flow_resolution_limit_results.json)、[核验](research_round_701_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round701.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=701,cumulative_tests=3309,numbered_scientific_files=1415,
            unique_protected_evidence_files=2809,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=702,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
