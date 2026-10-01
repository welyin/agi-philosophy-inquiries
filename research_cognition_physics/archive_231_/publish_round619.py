"""Publish verified 619, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round619_navigation_checks.json'
checked=core.read(HERE/'research_round_619_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round619_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第619轮完成：** [原手征物种、边界电流与共同区域拼接]({p}research_note_619.md)'
         '原Nambu物种的指定对称局部反射需32维酉匹配而交织秩至多2；'
         '整体透射可保原质量与来源，同一F还固定旋量通量和质量权重。'
         '三组、十四式通过，最新619／3101，1169份编号科学文件、1930份保护证据。'
         '[核验]({p}research_round_619_checks.json)、[条件账]({p}unified_physics_condition_ledger_619.md)。'
         '限明示连续分支与边界规则；重建、尺度及量子GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（619后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[620区域对象、共同量子态与几何来源]({p}round620_drafts/STATUS.md)，'
       '回查旧态、时间谱和来源结果，补真实联合接口；目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第619轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 265.' not in text
        text+='\n\n## 265. 手征区域、法向流与共同来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 170.' not in text
        text+='\n\n## 170. 连续旋量边界要求不推出维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—618轮。' in text and '最新科学轮次与检查数为618／3098' in text
        text=text.replace('完成231—618轮。','完成231—619轮。').replace('最新科学轮次与检查数为618／3098','最新科学轮次与检查数为619／3101')
    if p==HERE/'README.md':
        text+='\n\n## 第619轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|619|[原手征物种、边界电流与共同区域拼接](research_note_619.md)|[代码](joint_chiral_boundary_gluing.py)、[结果](joint_chiral_boundary_gluing_results.json)、[核验](research_round_619_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round619.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=619,cumulative_tests=3101,numbered_scientific_files=1169,
            unique_protected_evidence_files=1930,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=620,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
