"""Publish verified 629, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round629_navigation_checks.json'
checked=core.read(HERE/'research_round_629_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round629_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第629轮完成：** [原质量相位、拓扑权重与共同量子测度的压缩]({p}research_note_629.md)'
         '原群的四个拓扑角与原单代五个非零质量相位共用指数Jacobian；'
         '对常量模块重定相取商后剩三个角不变量，单独擦掉质量相位会改变它们。'
         '两组、十四式通过，最新629／3127，1199份编号科学文件、2001份保护证据。'
         '[核验]({p}research_round_629_checks.json)、[条件账]({p}unified_physics_condition_ledger_629.md)。'
         '仅压缩明示参数族；角值、连续区域过程、尺度及GR仍开放。')
order=('**当前执行顺序（629后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[630完整过程与连续物理的共同接口]({p}round630_drafts/STATUS.md)，'
       '回填623—629并核同一尺度的过程、测度及几何来源；不继续相位优化，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第629轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 275.' not in text
        text+='\n\n## 275. 原质量、拓扑权重与参数冗余\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 180.' not in text
        text+='\n\n## 180. 相位参数商不是空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—628轮。' in text and '最新科学轮次与检查数为628／3125' in text
        text=text.replace('完成231—628轮。','完成231—629轮。').replace('最新科学轮次与检查数为628／3125','最新科学轮次与检查数为629／3127')
    if p==HERE/'README.md':
        text+='\n\n## 第629轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|629|[原质量相位、拓扑权重与共同量子测度的压缩](research_note_629.md)|[代码](joint_topological_mass_phases.py)、[结果](joint_topological_mass_phases_results.json)、[核验](research_round_629_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round629.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=629,cumulative_tests=3127,numbered_scientific_files=1199,
            unique_protected_evidence_files=2001,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=630,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
