"""Publish verified 644, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round644_navigation_checks.json'
checked=core.read(HERE/'research_round_644_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round644_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第644轮完成：** [同一Gauss区域态的费米复制、二阶熵与几何来源]({p}research_note_644.md)'
         '原完整区域态、规范边界及CAR符号共用复制热历史，完整来源同时给二阶熵变分；'
         '两节点64模式的条件计算三路相符，原合法Gauss态另核符号。'
         '三组、十四式通过，最新644／3171，1244份编号科学文件、2111份保护证据。'
         '[核验]({p}research_round_644_checks.json)、[条件账]({p}unified_physics_condition_ledger_644.md)。'
         '数值非全Gauss热积分；面积律、连续匹配和GR未完成。')
order=('**当前执行顺序（644后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[645区域态与连续几何熵的共同尺度]({p}round645_drafts/STATUS.md)，'
       '核参考、区域、复制参数和尺度的真实映射；不继续小型热因子扫描，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第644轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 290.' not in text
        text+='\n\n## 290. 完整Gauss区域态的二阶熵与共同来源\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 195.' not in text
        text+='\n\n## 195. 区域复制一致不等于面积熵或GR成立\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—643轮。' in text and '最新科学轮次与检查数为643／3168' in text
        text=text.replace('完成231—643轮。','完成231—644轮。').replace('最新科学轮次与检查数为643／3168','最新科学轮次与检查数为644／3171')
    if p==HERE/'README.md':
        text+='\n\n## 第644轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|644|[同一Gauss区域态的费米复制、二阶熵与几何来源](research_note_644.md)|[代码](joint_region_replica_source.py)、[结果](joint_region_replica_source_results.json)、[核验](research_round_644_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round644.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=644,cumulative_tests=3171,numbered_scientific_files=1244,
            unique_protected_evidence_files=2111,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=645,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
