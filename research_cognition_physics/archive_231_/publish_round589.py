"""Publish verified 589, preserve all navigation history."""
import hashlib
import json
import os
from pathlib import Path
import verify_interaction_rounds as core
HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parent;ROOT=RESEARCH.parent
TARGET=HERE/'round589_navigation_checks.json'
checked=core.read(HERE/'research_round_589_checks.json')
assert checked['all_reported_checks_passed']
for key in ('new_file_hashes','preserved_draft_hashes'):
    for name,digest in checked[key].items():assert core.digest(HERE/name)==digest,name
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',RESEARCH/'RESEARCH_STATE.md',
       HERE/'README.md',HERE/'spatial_premise_closure_audit.md',HERE/'three_dimensional_four_conditions_review.md']
before={p:p.read_bytes() for p in paths}
folder=HERE/'navigation_before_round589_20261001'
assert not folder.exists() and not TARGET.exists()
summary=('**第589轮完成：** [完整空间度规、同一量子物质与实际记录]({p}research_note_589.md)'
         '一般正空间度规下保原商群物质、正自伴性、Gauss与实际仪器，共形分支精确恢复；'
         '同一来源的剪切响应已计算。六组、十四式及主代理核验通过，'
         '最新589／3005，1079份编号科学文件、1678份保护证据。'
         '[核验]({p}research_round_589_checks.json)、[条件账]({p}unified_physics_condition_ledger_589.md)。'
         '限固定背景与经典一阶光滑主部；动态量子几何、固定ℏ连续和自主装置仍开放。未取得新的独立代理审查。')
order=('**当前执行顺序（589后，优先于下方历史安排）：** 接续'
       '[590变化几何、记录与能源平衡]({p}round590_drafts/STATUS.md)，'
       '核时变正度规下同一量子过程、实际仪器注能与几何功，再定位内部动力几何缺口；统一目标不变。')
planned={}
for p,raw in before.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    newline='\r\n' if b'\r\n' in raw else '\n'
    text=raw.decode(enc).replace('\r\n','\n')
    assert '**第589轮完成：**' not in text
    prefix='research_cognition_physics/archive_231_/' if p.parent==ROOT else 'archive_231_/' if p.parent==RESEARCH else ''
    if p in paths[:5]:
        head,rest=text.split('\n\n',1)
        text=head+'\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n\n'+rest
    elif p.name=='spatial_premise_closure_audit.md':
        assert '\n## 235.' not in text
        text+='\n\n## 235. 完整空间度规与同一有限量子过程\n\n'+summary.format(p=prefix)+'\n\n'+order.format(p=prefix)+'\n'
    else:
        assert '\n## 140.' not in text
        text+='\n\n## 140. 完整固定背景不等于三维来源与动态引力\n\n'+summary.format(p=prefix)+'\n'
    if p.name=='research_direction.md':
        assert '完成231—588轮。' in text and '最新科学轮次与检查数为588／2999' in text
        text=text.replace('完成231—588轮。','完成231—589轮。').replace('最新科学轮次与检查数为588／2999','最新科学轮次与检查数为589／3005')
    if p==HERE/'README.md':
        text+='\n\n## 第589轮研究索引\n\n|轮次|报告|复算与核验|\n|---|---|---|\n|589|[完整空间度规与同一量子物质](research_note_589.md)|[代码](joint_full_spatial_metric.py)、[结果](joint_full_spatial_metric_results.json)、[核验](research_round_589_checks.json)|\n'
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
    temp=p.with_name(p.name+'.round589.tmp')
    with temp.open('xb') as stream:stream.write(raw)
    os.replace(temp,p)
report=dict(date='2026-10-01',latest_round=589,cumulative_tests=3005,numbered_scientific_files=1079,
            unique_protected_evidence_files=1678,navigation_files=7,navigation_links=links,broken_links=0,
            navigation_snapshots_preserved=True,complete_report_published=True,existing_results_preserved=True,
            next_round=590,active_goal_unchanged=True,stage_complete=False,all_checks_passed=True)
with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(report))
