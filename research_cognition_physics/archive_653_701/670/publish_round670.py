"""Publish verified 670, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round670_navigation_checks.json'
checked=core.read(HERE/'research_round_670_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round670_20261002'
assert not folder.exists() and not TARGET.exists()
summary=('**第670轮完成：** [非平坦规范场中的共同观测，以及无质量零模处的正则延伸]({p}research_note_670.md)'
         '原配对消去观测／质量字典中的逆无质量传播块；非交换非平坦链路、原质量与全部物理来源共同相容。'
         '无质量零模不再是该字典的排除条件；Wilson谱隙和动态正性边界保留。'
         '两组、十八式通过，最新670／3247，1322份编号科学文件、2376份保护证据。'
         '[核验]({p}research_round_670_checks.json)、[全条件账]({p}unified_physics_condition_ledger_670.md)。'
         '原完整过程、动态规范正性、连续及量子GR仍开放。')
order=('**当前执行顺序（670后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[671共同规范背景与物理反射泛函]({p}round671_drafts/STATUS.md)，'
       '核反射兼容的真实空间曲率背景，再接动态规范正时间组合；目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第670轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 316.' not in text
        text+='\n\n## 316. 共同规范字典消去无质量零模限制\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 221.' not in text
        text+='\n\n## 221. 正则观测字典与动态规范正性仍须区分\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—669轮。' in text and '最新科学轮次与检查数为669／3245' in text
        text=text.replace('完成231—669轮。','完成231—670轮。').replace('最新科学轮次与检查数为669／3245','最新科学轮次与检查数为670／3247')
    if p==HERE/'README.md':
        text+='\n\n## 第670轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|670|[非平坦规范场中的共同观测，以及无质量零模处的正则延伸](research_note_670.md)|[代码](joint_nonflat_mass_measure.py)、[结果](joint_nonflat_mass_measure_results.json)、[核验](research_round_670_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round670.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-02',latest_round=670,cumulative_tests=3247,numbered_scientific_files=1322,
            unique_protected_evidence_files=2376,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=671,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
