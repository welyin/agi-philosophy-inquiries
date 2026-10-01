"""Publish verified 625, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round625_navigation_checks.json'
checked=core.read(HERE/'research_round_625_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round625_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第625轮完成：** [记录与二阶来源共用一个有限尺度过程]({p}research_note_625.md)'
         '同一显式有限CP族共同逼近原有限历史的记录及零、一、二阶来源系数；'
         '初态、仪器尾项与Hessian必须共同输送。'
         '三组、二十式通过，最新625／3118，1187份编号科学文件、1973份保护证据。'
         '[核验]({p}research_round_625_checks.json)、[条件账]({p}unified_physics_condition_ledger_625.md)。'
         '限固定图与有限矩、有限来源菜单；数值仅为64维诊断，局域连续及GR仍开放。')
order=('**当前执行顺序（625后，优先于下方历史安排）：** 先尽量整合剩余共同条件，认知设计后置。'
       '接[626共同量子过程与连续作用的条件]({p}round626_drafts/STATUS.md)，'
       '按条件账核共同对象、态、来源与尺度，再处理相容性卡点；不继续谱精度优化，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第625轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 271.' not in text
        text+='\n\n## 271. 同一有限过程的记录与来源匹配\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 176.' not in text
        text+='\n\n## 176. 有限过程共同匹配不选择空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—624轮。' in text and '最新科学轮次与检查数为624／3115' in text
        text=text.replace('完成231—624轮。','完成231—625轮。').replace('最新科学轮次与检查数为624／3115','最新科学轮次与检查数为625／3118')
    if p==HERE/'README.md':
        text+='\n\n## 第625轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|625|[记录与二阶来源共用一个有限尺度过程](research_note_625.md)|[代码](joint_source_preserving_compression.py)、[结果](joint_source_preserving_compression_results.json)、[核验](research_round_625_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round625.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=625,cumulative_tests=3118,numbered_scientific_files=1187,
            unique_protected_evidence_files=1973,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=626,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
