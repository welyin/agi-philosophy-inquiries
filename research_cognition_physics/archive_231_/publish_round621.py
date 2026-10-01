"""Publish verified 621, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round621_navigation_checks.json'
checked=core.read(HERE/'research_round_621_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round621_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第621轮完成：** [原全Gauss热态的有限时间窗联合来源]({p}research_note_621.md)'
         '原固定图形式来源在同一Gauss热态经有限时间窗涂抹后有联合二阶矩及谱截断尾界；'
         '一般形式前提仍不足以保证瞬时噪声。'
         '两组、十二式通过，最新621／3106，1175份编号科学文件、1944份保护证据。'
         '[核验]({p}research_round_621_checks.json)、[条件账]({p}unified_physics_condition_ledger_621.md)。'
         '限固定图与明示来源族；完整因果响应、连续及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（621后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[622共同有限窗量子来源与二阶因果过程]({p}round622_drafts/STATUS.md)，'
       '核噪声、迟致核及接触项能否出自同一保Gauss过程；目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第621轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 267.' not in text
        text+='\n\n## 267. 全Gauss热态、时间窗与共同来源域\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 172.' not in text
        text+='\n\n## 172. 有限窗量子来源不构成维数选择\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—620轮。' in text and '最新科学轮次与检查数为620／3104' in text
        text=text.replace('完成231—620轮。','完成231—621轮。').replace('最新科学轮次与检查数为620／3104','最新科学轮次与检查数为621／3106')
    if p==HERE/'README.md':
        text+='\n\n## 第621轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|621|[原全Gauss热态的有限时间窗联合来源](research_note_621.md)|[代码](joint_smeared_gauss_noise.py)、[结果](joint_smeared_gauss_noise_results.json)、[核验](research_round_621_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round621.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=621,cumulative_tests=3106,numbered_scientific_files=1175,
            unique_protected_evidence_files=1944,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=622,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
