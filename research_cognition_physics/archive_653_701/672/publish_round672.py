"""Publish verified 672, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round672_navigation_checks.json'
checked=core.read(HERE/'research_round_672_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round672_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第672轮完成：** [独立规范历史、原手征权重与电动能的共同拼接条件]({p}research_note_672.md)'
         '原带相位权重在独立半区历史间出现数值负方向；解析主子式条件限制共同电热核，原商群实际热核复现通过与失败分支。'
         '固定边界候选不等于完整Gauss过程。'
         '两组、十六式通过，最新672／3251，1328份编号科学文件、2403份保护证据。'
         '[核验]({p}research_round_672_checks.json)、[全条件账]({p}unified_physics_condition_ledger_672.md)。'
         '共享边界、原完整过程、连续及量子GR仍开放。')
order=('**当前执行顺序（672后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[673共享边界与原Gauss闭合]({p}round673_drafts/STATUS.md)，'
       '核负方向的实际物理代数与完整有序过程，停止热时扫描；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第672轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 318.' not in text
        text+='\n\n## 318. 独立规范历史与原电动能的共同限制\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 223.' not in text
        text+='\n\n## 223. 固定边界核与完整Gauss过程的边界\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—671轮。' in text and '最新科学轮次与检查数为671／3249' in text
        text=text.replace('完成231—671轮。','完成231—672轮。').replace('最新科学轮次与检查数为671／3249','最新科学轮次与检查数为672／3251')
    if p==HERE/'README.md':
        text+='\n\n## 第672轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|672|[独立规范历史、原手征权重与电动能的共同拼接条件](research_note_672.md)|[代码](joint_gauge_history_kernel.py)、[结果](joint_gauge_history_kernel_results.json)、[核验](research_round_672_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round672.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=672,cumulative_tests=3251,numbered_scientific_files=1328,
            unique_protected_evidence_files=2403,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=673,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
