"""Publish verified 646, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round646_navigation_checks.json'
checked=core.read(HERE/'research_round_646_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round646_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第646轮完成：** [原传播核的共同连续窗口与重合来源边界]({p}research_note_646.md)'
         '原自由overlap谱在非零物理时间窗口共同匹配关联及两阶能源插入；'
         '同一谱的重合二阶能源矩有发散下界，分离时间收敛不能代替接触匹配。'
         '三组、十六式通过，最新646／3177，1250份编号科学文件、2127份保护证据。'
         '[核验]({p}research_round_646_checks.json)、[条件账]({p}unified_physics_condition_ledger_646.md)。'
         '限自由核；完整Gauss区域态、手征连续和动态引力仍开放。')
order=('**当前执行顺序（646后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[647固定图过程与动态几何的共同条件]({p}round647_drafts/STATUS.md)，'
       '回填共同对象和来源，先查旧约束；不继续自由核精度优化，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第646轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 292.' not in text
        text+='\n\n## 292. 同一原谱的有效连续窗口与来源边界\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 197.' not in text
        text+='\n\n## 197. 分离时间连续不替代重合来源匹配\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—645轮。' in text and '最新科学轮次与检查数为645／3174' in text
        text=text.replace('完成231—645轮。','完成231—646轮。').replace('最新科学轮次与检查数为645／3174','最新科学轮次与检查数为646／3177')
    if p==HERE/'README.md':
        text+='\n\n## 第646轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|646|[原传播核的共同连续窗口与重合来源边界](research_note_646.md)|[代码](joint_overlap_continuum_sources.py)、[结果](joint_overlap_continuum_sources_results.json)、[核验](research_round_646_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round646.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=646,cumulative_tests=3177,numbered_scientific_files=1250,
            unique_protected_evidence_files=2127,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=647,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
