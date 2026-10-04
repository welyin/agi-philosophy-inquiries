"""Publish verified 707, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round707_navigation_checks.json'
checked=core.read(HERE/'research_round_707_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round707_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第707轮完成：** [曲目标空间分块与原共同过程]({p}research_note_707.md)'
         '原H⁵中点／相对向量给保Gauss的真实标量分块，原热来源、中点记录预算与全部质量共同输送；'
         '保内部S⁴及全部CAR，仍有记忆。'
         '三组、二十式通过，最新707／3323，1433份编号科学文件、2900份保护证据。'
         '[核验]({p}research_round_707_checks.json)、[全条件账]({p}unified_physics_condition_ledger_707.md)。'
         '不当裸单节点或连续极限，旧空间与目标不变。')
order=('**当前执行顺序（707后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[708嵌套合并与共同质量资料]({p}round708_drafts/STATUS.md)，'
       '核中点以外的幅度是否同时承担组合及原集体质量；'
       '不重复固定图精度、一般偏迹反例或Gaussian记忆，四分支及699范围保留。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第707轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 353.' not in text
        text+='\n\n## 353. 原非线性空间分块、Gauss热态与质量共同验收\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 258.' not in text
        text+='\n\n## 258. 旧空间接口继承，内部曲目标不当物理空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—706轮。' in text and '最新科学轮次与检查数为706／3320' in text
        text=text.replace('完成231—706轮。','完成231—707轮。').replace('最新科学轮次与检查数为706／3320','最新科学轮次与检查数为707／3323')
    if p==HERE/'README.md':
        text+='\n\n## 第707轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|707|[曲目标空间分块与原共同过程](research_note_707.md)|[代码](joint_geodesic_spatial_block.py)、[结果](joint_geodesic_spatial_block_results.json)、[核验](research_round_707_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round707.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=707,cumulative_tests=3323,numbered_scientific_files=1433,
            unique_protected_evidence_files=2900,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=708,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
