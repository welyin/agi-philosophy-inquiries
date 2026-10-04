"""Publish727 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round727_navigation_checks.json'
checked=core.read(HERE/'research_round_727_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round727_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第727轮完成：** [原物质基态、同一经典来源与保留的量子涨落]({p}research_note_727.md)局部有隙及Gauss提升条件下，同次记录的完整能源集中于Hb+e；原质量和三方向传播共同基态保非零真实来源噪声，不能用带内压缩替代。三组、十六式通过，最新727／3386，1493份编号科学文件、3186份保护证据。[核验]({p}research_round_727_checks.json)、[全条件账]({p}unified_physics_condition_ledger_727.md)。原经典来源、全图带隙及动态引力仍须核。'
order='**当前执行顺序（727后，优先于下方历史安排）：** 接[728宏观来源与真实时间尺度]({p}round728_drafts/STATUS.md)，回查591及630—633的亚隙谱／涂抹，核同一参考的实际宏观响应与变化过程。停止占据、单带函数及波包优化，不要求每个微观瞬时源零方差；旧空间、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第727轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 373.' not in text
        text+='\n\n## 373. 原共同基态、带间来源与宏观任务\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 278.' not in text
        text+='\n\n## 278. 旧空间合同保持，基态能源不等于零来源噪声\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—726轮。' in text and '最新科学轮次与检查数为726／3383' in text
        text=text.replace('完成231—726轮。','完成231—727轮。').replace('最新科学轮次与检查数为726／3383','最新科学轮次与检查数为727／3386')
    if p==HERE/'README.md':
        text+='\n\n## 第727轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|727|[原物质基态、同一经典来源与保留的量子涨落](research_note_727.md)|[代码](joint_matter_ground_source.py)、[结果](joint_matter_ground_source_results.json)、[核验](research_round_727_checks.json)|\n'
    planned[p]=text.replace('\n',nl).encode(enc)
links=0
for p,raw in planned.items():
    for link in core.link_parser()(raw.decode('utf-8-sig')):
        assert (p.parent/link).resolve().exists(),(p,link)
        links+=1
assert all(p.read_bytes()==raw for p,raw in before.items())
folder.mkdir(exist_ok=False);manifest={}
for i,(p,raw) in enumerate(before.items()):
    name=f'{i}_{p.name}'
    with (folder/name).open('xb') as f:f.write(raw)
    manifest[p.relative_to(ROOT).as_posix()]=dict(snapshot=name,sha256=hashlib.sha256(raw).hexdigest())
with (folder/'manifest.json').open('x',encoding='utf8') as f:json.dump(manifest,f,ensure_ascii=False,indent=2)
for p,raw in planned.items():
    assert p.read_bytes()==before[p]
    temp=p.with_name(p.name+'.round727.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=727,cumulative_tests=3386,numbered_scientific_files=1493,
    unique_protected_evidence_files=3186,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=728,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
