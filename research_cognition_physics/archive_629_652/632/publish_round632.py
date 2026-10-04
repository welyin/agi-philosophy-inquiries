"""Publish verified 632, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round632_navigation_checks.json'
checked=core.read(HERE/'research_round_632_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round632_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第632轮完成：** [同一真空势、局部接触项与背景方程的共同匹配]({p}research_note_632.md)'
         '原完整质量的有限势给出零导数Weyl项与q—σ接触矩阵；'
         '声明的三个背景匹配条件同时约束物质驻定和几何零动量项。'
         '三组、十六式通过，最新632／3136，1208份编号科学文件、2022份保护证据。'
         '[核验]({p}research_round_632_checks.json)、[条件账]({p}unified_physics_condition_ledger_632.md)。'
         '背景条件为输入；全动量Ward、实际记录后态、图映射及GR仍开放。')
order=('**当前执行顺序（632后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[633连续记录更新、态与共同来源]({p}round633_drafts/STATUS.md)，'
       '复用旧记录结果，核同一连续物质的局域操作和来源；不继续真空参数扫描，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第632轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 278.' not in text
        text+='\n\n## 278. 原真空势与几何接触项的共同匹配\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 183.' not in text
        text+='\n\n## 183. 共同背景匹配不选择空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—631轮。' in text and '最新科学轮次与检查数为631／3133' in text
        text=text.replace('完成231—631轮。','完成231—632轮。').replace('最新科学轮次与检查数为631／3133','最新科学轮次与检查数为632／3136')
    if p==HERE/'README.md':
        text+='\n\n## 第632轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|632|[同一真空势、局部接触项与背景方程的共同匹配](research_note_632.md)|[代码](joint_background_contact_matching.py)、[结果](joint_background_contact_matching_results.json)、[核验](research_round_632_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round632.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=632,cumulative_tests=3136,numbered_scientific_files=1208,
            unique_protected_evidence_files=2022,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=633,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
