"""Publish verified 640, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round640_navigation_checks.json'
checked=core.read(HERE/'research_round_640_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round640_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第640轮完成：** [一次真实空间合并的共同热态、时间核与几何来源]({p}research_note_640.md)'
         '原中性二次模相邻合并后，精确热态Hamiltonian不能代替二频过程；'
         '同一空间消元须保记忆核及几何归一化来源。'
         '三组、十四式通过，最新640／3160，1232份编号科学文件、2081份保护证据。'
         '[核验]({p}research_round_640_checks.json)、[条件账]({p}unified_physics_condition_ledger_640.md)。'
         '限树级二次分支；完整Gauss、跨图连续及GR仍开放。')
order=('**当前执行顺序（640后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[641原空间消元后的共同区域过程]({p}round641_drafts/STATUS.md)，'
       '核原记录、共同参考及来源能否共用有记忆过程；保留尺度正则性条件，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第640轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 286.' not in text
        text+='\n\n## 286. 原空间合并的共同状态、时间核与来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 191.' not in text
        text+='\n\n## 191. 二次空间合并不等于完整连续物理\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—639轮。' in text and '最新科学轮次与检查数为639／3157' in text
        text=text.replace('完成231—639轮。','完成231—640轮。').replace('最新科学轮次与检查数为639／3157','最新科学轮次与检查数为640／3160')
    if p==HERE/'README.md':
        text+='\n\n## 第640轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|640|[一次真实空间合并的共同热态、时间核与几何来源](research_note_640.md)|[代码](joint_spatial_block_reference.py)、[结果](joint_spatial_block_reference_results.json)、[核验](research_round_640_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round640.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=640,cumulative_tests=3160,numbered_scientific_files=1232,
            unique_protected_evidence_files=2081,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=641,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
