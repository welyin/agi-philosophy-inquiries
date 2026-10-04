"""Publish verified 633, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round633_navigation_checks.json'
checked=core.read(HERE/'research_round_633_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round633_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第633轮完成：** [连续记录、态正则性与同一物质的几何来源]({p}research_note_633.md)'
         '原混合中性场的光滑记录后态保短距结构与有限能源，并共同给质量和应力来源；'
         '同一场两分离模排除瞬时分布式读取。'
         '四组、十八式通过，最新633／3140，1211份编号科学文件、2029份保护证据。'
         '[核验]({p}research_round_633_checks.json)、[条件账]({p}unified_physics_condition_ledger_633.md)。'
         '限自由固定背景；真实仪器、图映射和GR仍开放。')
order=('**当前执行顺序（633后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[634记录条件化、来源与几何约束]({p}round634_drafts/STATUS.md)，'
       '核同一过程的状态边界与守恒；不继续波包或装置优化，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第633轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 279.' not in text
        text+='\n\n## 279. 原连续记录与物质几何来源的共同连接\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 184.' not in text
        text+='\n\n## 184. 记录条件的整合不选择空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—632轮。' in text and '最新科学轮次与检查数为632／3136' in text
        text=text.replace('完成231—632轮。','完成231—633轮。').replace('最新科学轮次与检查数为632／3136','最新科学轮次与检查数为633／3140')
    if p==HERE/'README.md':
        text+='\n\n## 第633轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|633|[连续记录、态正则性与同一物质的几何来源](research_note_633.md)|[代码](joint_continuum_record_sources.py)、[结果](joint_continuum_record_sources_results.json)、[核验](research_round_633_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round633.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=633,cumulative_tests=3140,numbered_scientific_files=1211,
            unique_protected_evidence_files=2029,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=634,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
