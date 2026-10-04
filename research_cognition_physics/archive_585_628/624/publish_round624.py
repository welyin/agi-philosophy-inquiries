"""Publish verified 624, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round624_navigation_checks.json'
checked=core.read(HERE/'research_round_624_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round624_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第624轮完成：** [实际记录、条件响应与共同来源过程]({p}research_note_624.md)'
         '原仪器与完整CAR来源共享有限历史二阶过程；条件筛选须保记录概率导数，'
         '实际来源得分限制记录合并后的分辨性。'
         '三组、十四式通过，最新624／3115，1184份编号科学文件、1966份保护证据。'
         '[核验]({p}research_round_624_checks.json)、[条件账]({p}unified_physics_condition_ledger_624.md)。'
         '限固定图与有限矩历史；连续、动态引力与GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（624后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[625同一有限尺度映射的记录与来源合同]({p}round625_drafts/STATUS.md)，'
       '回查原谱近似和虚跃迁，核记录、噪声、迟致及接触项的共同匹配；目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第624轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 270.' not in text
        text+='\n\n## 270. 实际记录、条件选择与共同来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 175.' not in text
        text+='\n\n## 175. 来源分辨性不构成维数选择\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—623轮。' in text and '最新科学轮次与检查数为623／3112' in text
        text=text.replace('完成231—623轮。','完成231—624轮。').replace('最新科学轮次与检查数为623／3112','最新科学轮次与检查数为624／3115')
    if p==HERE/'README.md':
        text+='\n\n## 第624轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|624|[实际记录、条件响应与共同来源过程](research_note_624.md)|[代码](joint_record_source_response.py)、[结果](joint_record_source_response_results.json)、[核验](research_round_624_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round624.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=624,cumulative_tests=3115,numbered_scientific_files=1184,
            unique_protected_evidence_files=1966,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=625,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
