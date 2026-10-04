"""Publish verified 691, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round691_navigation_checks.json'
checked=core.read(HERE/'research_round_691_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round691_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第691轮完成：** [中性相位、完整正性约化与奇Gauss测试]({p}research_note_691.md)'
         '原中性相位在已正反射商上有酉实现；完整无质量RP精确约化到保留全部辅助／规范平均的剩余泛函，偶全代数仍需奇Gauss测试。'
         '两组、十六式通过，最新691／3289，1385份编号科学文件、2641份保护证据。'
         '[核验]({p}research_round_691_checks.json)、[全条件账]({p}unified_physics_condition_ledger_691.md)。'
         '酉对称不冒充局部电荷或仪器；剩余动态RP、H_F、共同尺度及量子GR仍开放。')
order=('**当前执行顺序（691后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[692原奇Gauss复合与完整泛函]({p}round692_drafts/STATUS.md)，'
       '保留全部余子式与原平均，复用675—680，不继续修补已分离的中性部门；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第691轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 337.' not in text
        text+='\n\n## 337. 完整物理正性的中性分离与奇Gauss条件\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 242.' not in text
        text+='\n\n## 242. 旧空间合同复用，酉对称不冒充测量仪器\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—690轮。' in text and '最新科学轮次与检查数为690／3287' in text
        text=text.replace('完成231—690轮。','完成231—691轮。').replace('最新科学轮次与检查数为690／3287','最新科学轮次与检查数为691／3289')
    if p==HERE/'README.md':
        text+='\n\n## 第691轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|691|[中性相位、完整正性约化与奇Gauss测试](research_note_691.md)|[代码](joint_neutral_symmetry_reduction.py)、[结果](joint_neutral_symmetry_reduction_results.json)、[核验](research_round_691_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round691.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=691,cumulative_tests=3289,numbered_scientific_files=1385,
            unique_protected_evidence_files=2641,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=692,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
