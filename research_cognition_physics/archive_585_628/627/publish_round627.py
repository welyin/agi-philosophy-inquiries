"""Publish verified 627, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round627_navigation_checks.json'
checked=core.read(HERE/'research_round_627_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round627_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第627轮完成：** [物质尺度划分、反常相位与共同参考域]({p}research_note_627.md)'
         '原一代完整模块的无补偿删减受反常联立限制，64子集仅4个通过所列必要检查；'
         '原全热态也不支持夸克质量的统一正点态隙。'
         '消去部门须保其测度相位、适用态域及来源。'
         '两组、十二式通过，最新627／3123，1193份编号科学文件、1987份保护证据。'
         '[核验]({p}research_round_627_checks.json)、[条件账]({p}unified_physics_condition_ledger_627.md)。'
         '未完成全局反常、连续手征过程或GR；不排除有补偿项的有效理论。')
order=('**当前执行顺序（627后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[628完整物质、手征测度与共同区域过程]({p}round628_drafts/STATUS.md)，'
       '核原一代、配置域、态、拼接及真实过程的共同接口；不重复阈值扫描，目标不改。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第627轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 273.' not in text
        text+='\n\n## 273. 原物质划分、测度相位与质量域\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 178.' not in text
        text+='\n\n## 178. 物质子集约束不选择空间维数\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—626轮。' in text and '最新科学轮次与检查数为626／3121' in text
        text=text.replace('完成231—626轮。','完成231—627轮。').replace('最新科学轮次与检查数为626／3121','最新科学轮次与检查数为627／3123')
    if p==HERE/'README.md':
        text+='\n\n## 第627轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|627|[物质尺度划分、反常相位与共同参考域](research_note_627.md)|[代码](joint_anomaly_scale_partition.py)、[结果](joint_anomaly_scale_partition_results.json)、[核验](research_round_627_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round627.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=627,cumulative_tests=3123,numbered_scientific_files=1193,
            unique_protected_evidence_files=1987,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=628,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
