"""Publish verified 615, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round615_navigation_checks.json'
checked=core.read(HERE/'research_round_615_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round615_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第615轮完成：** [原子群的精确辅助测度与共同响应]({p}research_note_615.md)'
         '同一原子群平坦族的辅助Pfaffian积分已精确完成且严格为正；'
         '它与物理行列式共同产生来源，删除辅助项会漏掉非零响应。'
         '三组、十四式通过，最新615／3090，1157份编号科学文件、1900份保护证据。'
         '[核验]({p}research_round_615_checks.json)、[条件账]({p}unified_physics_condition_ledger_615.md)。'
         '限明确有限欧氏族；一般局域性、非零质量、重建及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（615后，优先于下方历史安排）：** 继续先合并共同条件，认知设计后置。'
       '接[616原物质、全局几何与量子态]({p}round616_drafts/STATUS.md)，'
       '核spin存在／选择、原群及状态来源，不继续孤立优化本轮积分；目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第615轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 261.' not in text
        text+='\n\n## 261. 同一测度、物理权重与共同来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 166.' not in text
        text+='\n\n## 166. 有限测度正性与全局物理前提\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—614轮。' in text and '最新科学轮次与检查数为614／3087' in text
        text=text.replace('完成231—614轮。','完成231—615轮。').replace('最新科学轮次与检查数为614／3087','最新科学轮次与检查数为615／3090')
    if p==HERE/'README.md':
        text+='\n\n## 第615轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|615|[原子群的精确辅助测度与共同响应](research_note_615.md)|[代码](joint_subgroup_measure_source.py)、[结果](joint_subgroup_measure_source_results.json)、[核验](research_round_615_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round615.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=615,cumulative_tests=3090,numbered_scientific_files=1157,
            unique_protected_evidence_files=1900,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=616,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
