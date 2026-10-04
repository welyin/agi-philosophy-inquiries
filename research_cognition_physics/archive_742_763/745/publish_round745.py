"""Publish745 with preserved navigation snapshots and conflict checks."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round745_navigation_checks.json'
checked=core.read(HERE/'research_round_745_checks.json');assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths};folder=HERE/'navigation_before_round745_20261004'
assert not folder.exists() and not TARGET.exists()
summary='**第745轮完成：** [原联合历史的严格信号]({p}research_note_745.md)原一节点完整H的联合占据差经精确系数和带解析余项的区间积分认证；任意短邻域内存在可区分等待，足够小正读取间隔保差异。未给可用时间窗口或连通图结果。两组、十六式通过，最新745／3428，1547份编号科学文件、3436份保护证据。[核验]({p}research_round_745_checks.json)、[条件账]({p}unified_physics_condition_ledger_745.md)。'
order='**当前执行顺序（745后，优先于下方历史安排）：** 接[746原连通图与共同过程]({p}round746_drafts/STATUS.md)，保完整跳跃、Gauss、未知输入及资源，区分弱耦合存在性和原固定边系数。旧空间、604、649／699及统一目标保持。'
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8';nl='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n');assert '**第745轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 391.' not in text
        text+='\n\n## 391. 原联合历史的严格信号\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 296.' not in text
        text+='\n\n## 296. 旧空间合同保持，严格原历史信号进入图接口\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—744轮。' in text and '最新科学轮次与检查数为744／3426' in text
        text=text.replace('完成231—744轮。','完成231—745轮。').replace('最新科学轮次与检查数为744／3426','最新科学轮次与检查数为745／3428')
    if p==HERE/'README.md':
        text+='\n\n## 第745轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|745|[原联合历史的严格信号](research_note_745.md)|[代码](joint_native_history_certificate.py)、[结果](joint_native_history_certificate_results.json)、[核验](research_round_745_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round745.tmp')
    with temp.open('xb') as f:f.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-04',latest_round=745,cumulative_tests=3428,numbered_scientific_files=1547,
    unique_protected_evidence_files=3436,navigation_files=7,navigation_links=links,broken_links=0,
    navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
    next_round=746,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
